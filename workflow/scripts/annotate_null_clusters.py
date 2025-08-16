#!/usr/bin/env python3
"""
Null Cluster Generation and Annotation Script

This script generates null clusters for statistical comparison with real gene expression clusters.
It creates clusters of neighbroing genes that match the size distribution of real clusters,
then annotates them with the same genomic and functional features for fair comparison.
"""

import logging
import sys
import os
import numpy as np
import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sp
from tqdm import tqdm
from pathlib import Path

# Import local modules
try:
    from annotate_clusters import *
    from get_pcs import get_pc_bed
    from annotate_pcs import get_annotate_pcs
    from residualize import calculate_residual
except ImportError:
    # Fallback for when script is run directly
    import sys
    sys.path.append(str(Path(__file__).parent))
    from annotate_clusters import *
    from get_pcs import get_pc_bed
    from annotate_pcs import get_annotate_pcs
    from residualize import calculate_residual

# Configuration - use environment variables with defaults
PREFIX = os.getenv('PCQTL_PREFIX', '/home/klawren/oak/pcqtls')


def setup_logging(level=logging.INFO):
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def get_null_clusters(expressed_gencode, cluster_size, cluster_df=None):
    """
    Generate null clusters by sliding window across chromosomes.
    
    Creates clusters of consecutive genes that can be used as null controls
    for statistical comparison with real gene expression clusters.
    
    Args:
        expressed_gencode (pd.DataFrame): GENCODE annotations for expressed genes
        cluster_size (int): Number of genes per null cluster
        cluster_df (pd.DataFrame, optional): Real clusters to exclude overlapping genes
    
    Returns:
        pd.DataFrame: Null clusters with cluster_id, cluster_size, chr, num_genes columns
    """
    logger = logging.getLogger(__name__)
    logger.info(f'Generating null clusters of size {cluster_size}')
    
    # Sort genes by genomic position
    expressed_gencode = expressed_gencode.sort_values(['chr', 'start', 'end'])
    
    # Generate null clusters on a per-chromosome basis
    null_cluster_dfs = []
    for chr_id in range(1, 23):
        logger.debug(f'Processing chromosome {chr_id}')
        chr_subset_gencode = expressed_gencode[expressed_gencode['chr'] == f'chr{chr_id}']
        
        if len(chr_subset_gencode) < cluster_size:
            logger.debug(f'Chromosome {chr_id} has fewer than {cluster_size} genes, skipping')
            continue
        
        # Generate cluster IDs with underscore joining and alphanumeric sorting
        cluster_ids = []
        for i in range(len(chr_subset_gencode) - cluster_size + 1):
            cluster_transcripts = chr_subset_gencode['gene_id'].iloc[i:i+cluster_size].astype(str).tolist()
            cluster_ids.append('_'.join(sorted(cluster_transcripts)))
        
        # Calculate cluster sizes (genomic span)
        chr_sizes = chr_subset_gencode['end'].shift(-(cluster_size-1)) - chr_subset_gencode['start']
        chr_sizes = chr_sizes.iloc[:-(cluster_size-1)]

        chr_null_df = pd.DataFrame({'cluster_id': cluster_ids, 'cluster_size': chr_sizes, 'chr': chr_id})

        # Select those that are not already clusters (exclude overlapping genes)
        if cluster_df is not None:
            try:
                cluster_ids = np.concatenate(cluster_df['cluster_id'].str.split('_').values)
                chr_null_df['transcript_list'] = chr_null_df['cluster_id'].str.split('_')
                
                # Exclude any null clusters that contain two genes both in a real cluster
                null_clusters_exploded = chr_null_df.explode('transcript_list')
                null_clusters_exploded['in_cluster'] = null_clusters_exploded['transcript_list'].isin(cluster_ids).astype(int)
                null_clusters_cluster_df_count = null_clusters_exploded.groupby('cluster_id').agg({'in_cluster': sum})
                excluded_clusters = null_clusters_cluster_df_count[null_clusters_cluster_df_count['in_cluster'] > 1].index
                
                # Add to the nulls
                chr_null_df = chr_null_df[~chr_null_df['cluster_id'].isin(excluded_clusters)]
                logger.debug(f'Excluded {len(excluded_clusters)} clusters with overlapping genes on chromosome {chr_id}')
            except TypeError:
                # No subselection wanted, i.e. None passed for cluster_df
                pass
        
        null_cluster_dfs.append(chr_null_df)

    # Combine all chromosomes
    null_df = pd.concat(null_cluster_dfs)
    null_df.reset_index(drop=True, inplace=True)
    null_df['chr'] = null_df['chr'].apply(lambda x: f'chr{x}')
    null_df['num_genes'] = cluster_size
    
    logger.info(f'Generated {len(null_df)} null clusters')
    return null_df


def resample_dist(target, candidate_pool, n, seed=126124):   
    """
    Match a target distribution via weighted sampling from a candidate pool.
    
    Uses kernel density estimation to create sampling weights that match
    the target distribution from the candidate pool.
    
    Args:
        target (np.ndarray): Target distribution (1D array, values 0-1)
        candidate_pool (np.ndarray): Candidate pool (1D array, values 0-1)
        n (int): Number of indices to return
        seed (int): Random seed for reproducibility
    
    Returns:
        np.ndarray: n indices to elements in candidate_pool
    """
    rng = np.random.default_rng(seed)
    target_prob = sp.stats.gaussian_kde(target)
    candidate_prob = sp.stats.gaussian_kde(candidate_pool)

    bins = np.arange(0, 1, 0.0001)
    sampling_weight = target_prob(bins) / candidate_prob(bins)
    pool_bins = np.searchsorted(bins, candidate_pool) - 1
    pool_probability = sampling_weight[pool_bins]/sampling_weight[pool_bins].sum()

    return rng.choice(candidate_pool.size, size=n, replace=True, p=pool_probability)


def get_resamp_null_cluster(null_df, cluster_df, plot=False, number_null=5000):
    """
    Resample null clusters to match the size distribution of real clusters.
    
    Args:
        null_df (pd.DataFrame): Null clusters to resample
        cluster_df (pd.DataFrame): Real clusters to match distribution
        plot (bool): Whether to plot distributions before/after resampling
        number_null (int): Number of null clusters to return
    
    Returns:
        pd.DataFrame: Resampled null clusters matching real cluster distribution
    """
    logger = logging.getLogger(__name__)
    logger.info(f'Resampling {number_null} null clusters to match real cluster distribution')
    
    # Note: cluster_df should be restricted only to clusters with matching number of transcripts
    
    join_df = pd.concat([cluster_df, null_df], keys=['cluster', 'null'], names=['is_cluster', 'id'])
    
    if plot:
        # Size distribution before resampling
        sns.kdeplot(join_df.reset_index(), x='cluster_tss_size', hue='is_cluster', 
                   bw_adjust=.3, common_norm=False, log_scale=(10, None))
        plt.title('Distance distribution before resampling')
        plt.show()

    # Add a normalized cluster size column to resample on
    cluster_df['normed_log_size'] = np.log10(cluster_df['cluster_tss_size'])/np.log10(join_df['cluster_tss_size'].max())
    null_df['normed_log_size'] = np.log10(null_df['cluster_tss_size'])/np.log10(join_df['cluster_tss_size'].max())

    # Resample to match distance distribution
    resamp_idxs = resample_dist(cluster_df['normed_log_size'], null_df['normed_log_size'], n=number_null)
    resamp_null_df = null_df.reset_index().iloc[resamp_idxs]

    if plot:
        # Size distribution after resampling
        join_df = pd.concat([cluster_df, resamp_null_df], keys=['cluster', 'null'], names=['is_cluster', 'id'])
        sns.kdeplot(join_df, x='cluster_tss_size', hue='is_cluster', 
                   bw_adjust=.3, common_norm=False, log_scale=(10, None))
        plt.title('Distance distribution after resampling')
        plt.show()

    logger.info(f'Resampled {len(resamp_null_df)} null clusters')
    return resamp_null_df


def get_null_pcs_annotated(config, cluster_size, tissue_id):
    """
    Generate and annotate null principal components.
    
    This function creates null clusters, generates PCs from them, and annotates
    the PCs with genomic features for comparison with real cluster PCs.
    
    Args:
        config (dict): Configuration dictionary with file paths
        cluster_size (int): Number of genes per null cluster
        tissue_id (str): GTEx tissue identifier
    
    Returns:
        pd.DataFrame: Annotated null principal components
    """
    logger = logging.getLogger(__name__)
    logger.info(f'Generating annotated null PCs for tissue {tissue_id}, cluster size {cluster_size}')
    
    # Load data
    gid_gencode, full_gencode = load_gencode()
    expression_path = "{}/{}/{}.v8.normalized_expression.bed".format(PREFIX, config["expression_dir"], tissue_id)
    expression_df = pd.read_table(expression_path)
            expressed_gencode = full_gencode[full_gencode['gene_id'].isin(expression_df['gene_id'])]
    cluster_df = pd.read_csv('{}/{}/{}.clusters.txt'.format(PREFIX, config['clusters_dir'], tissue_id), index_col=0)
    
    # Generate null clusters
    null_clusters = get_null_clusters(expressed_gencode, cluster_size, cluster_df=cluster_df)
    
    # Load covariates and residualize expression
    covariates_path = "{}/{}/{}.v8.covariates.txt".format(PREFIX, config["covariates_dir"], tissue_id)
    covariates_df = pd.read_table(covariates_path, index_col=0).T
    residual_exp = calculate_residual(expression_df[covariates_df.index], covariates_df, center=True)
    residual_exp = pd.DataFrame(residual_exp, columns=covariates_df.index, index=expression_df['gene_id'])

    # Residualize the expression 
    expression_df_gid = expression_df.set_index('gene_id')
    expression_df_res = expression_df_gid.copy()
    expression_df_res[expression_df.columns[4:]] = residual_exp

    # Get cluster-based expression (from filter_expression_clusters)
    logger.info('Creating null cluster expression files...')
    cluster_sets = []
    for idx, row in tqdm(null_clusters.iterrows(), total=len(null_clusters)):
        cluster_set = expression_df_res.loc[row['cluster_id'].split('_')].copy()
        cluster_set['start'] = cluster_set['start'].min()
        cluster_set['end'] = cluster_set['end'].max()
        cluster_set.insert(3, 'gene_id', ['_'.join([*sorted(cluster_set.index.values), 'e', gid]) for gid in cluster_set.index])
        cluster_sets.append(cluster_set)

    cluster_expression = pd.concat(cluster_sets)
    
    # Get the null PCs
    logger.info('Generating null principal components...')
    null_pcs = get_pc_bed(null_clusters, cluster_expression, covariates_df)
    null_pcs['cluster_id'] = null_pcs['gene_id'].str.split('_pc').str[0]
    null_pcs['pc_id'] = null_pcs['gene_id'].str.split('_pc').str[1].astype('float')

    # Annotate null PCs
    logger.info('Annotating null principal components...')
    annotated_null_pcs = get_annotate_pcs(null_pcs.reset_index(drop=True), cluster_expression.drop(columns='gene_id'))
    
    logger.info(f'Generated {len(annotated_null_pcs)} annotated null PCs')
    return annotated_null_pcs


def main():
    """
    Main function to generate and annotate null clusters.
    
    Parses command-line arguments and runs the null cluster generation pipeline.
    """
    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description='Generate and annotate null clusters for statistical comparison',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('-t', '--tissue-id', required=True,
                       help='GTEx tissue identifier')
    parser.add_argument('-c', '--clusters', required=True,
                       help='Real cluster CSV file')
    parser.add_argument('-e', '--expression', required=True,
                       help='Normalized expression BED file')
    parser.add_argument('-co', '--covariates', required=True,
                       help='Covariates file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output file for annotated null clusters')
    
    # Optional arguments
    parser.add_argument('--cluster-size', type=int, default=2,
                       help='Number of genes per null cluster (default: 2)')
    parser.add_argument('--exclude-cluster-genes', type=int, default=1,
                       help='Exclude genes from real clusters? 0 for no, 1 for yes (default: 1)')
    parser.add_argument('--distance-matched', type=int, default=0,
                       help='Match to distances in real clusters? 0 for no, 1 for yes (default: 0)')
    parser.add_argument('--gencode', 
                       help='GENCODE annotation file')
    parser.add_argument('--full-abc', 
                       help='ABC predictions file')
    parser.add_argument('--abc-match', 
                       help='ABC-GTEx tissue matching file')
    parser.add_argument('--ctcf-match', 
                       help='CTCF-GTEx tissue matching file')
    parser.add_argument('--ctcf-dir', 
                       help='Directory containing CTCF binding site files')
    parser.add_argument('--paralog', 
                       help='Paralog data file')
    parser.add_argument('--go', 
                       help='GO term file')
    parser.add_argument('--cross-map', 
                       help='Cross-mappability file')
    parser.add_argument('--tad', 
                       help='TAD boundary file')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')

    # Parse arguments
    args = parser.parse_args()
    
    if args.verbose:
        setup_logging(logging.DEBUG)
    
    try:
        logger.info(f'Starting null cluster generation for tissue: {args.tissue}')
        logger.info(f'Cluster size: {args.cluster_size}, Exclude cluster genes: {args.exclude_cluster_genes}')
        
        # Load the real clusters
        logger.info(f'Loading real clusters from {args.cluster_path}')
        cluster_df = pd.read_csv(args.cluster_path, index_col=0)

        # Load expressed GENCODE (to make null clusters out of)
        # This assumes only protein coding genes
        logger.info('Loading GENCODE annotations...')
        gid_gencode, full_gencode = load_gencode(args.gencode_path)
        expression_df = pd.read_table(args.expression_path)
        expressed_gencode = full_gencode[full_gencode['gene_id'].isin(expression_df['gene_id'])]
        expressed_gencode = expressed_gencode.sort_values(['chr', 'start', 'end'])
        logger.info(f'Found {len(expressed_gencode)} expressed genes')

        # Generate the null clusters
        logger.info('Generating null clusters...')
        if args.exclude_cluster_genes == 1:
            # Excluding genes from real clusters
            logger.info('Excluding genes from real clusters')
            null_clusters = get_null_clusters(expressed_gencode, args.cluster_size, cluster_df=cluster_df).sample(5000)
        else:
            # Including all genes
            logger.info('Including all genes in null clusters')
            null_clusters = get_null_clusters(expressed_gencode, args.cluster_size, cluster_df=None).sample(5000)

        if args.distance_matched:
            # Distance matched null clusters
            logger.info('Matching null cluster distances to real clusters')
            null_clusters = get_resamp_null_cluster(null_clusters, cluster_df[cluster_df['num_genes']==args.cluster_size], plot=False)

        # Annotate the null clusters
        logger.info('Annotating null clusters...')
        load_and_annotate(null_clusters, args.tissue, args.covariates_path, args.expression_path, 
                          args.gencode_path, args.full_abc_path, args.abc_match_path, args.ctcf_match_path,
                          args.ctcf_dir, args.paralog_path, args.go_path, args.cross_map_path, args.tad_path)   
        
        # Write output
        logger.info(f'Writing {len(null_clusters)} annotated null clusters to {args.out_path}')
        null_clusters.to_csv(args.out_path, index=False)
        logger.info('Null cluster generation complete!')
        
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        sys.exit(1)
    except ValueError as e:
        logger.error(f"Invalid input: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
