#!/usr/bin/env python3
"""
Gene Expression Cluster Annotation Script

This script annotates gene expression clusters with various genomic and functional features
including ABC enhancer-gene connections, CTCF binding sites, TAD boundaries, paralogous genes,
GO terms. 
"""

import argparse
import ast
import logging
import os
import sys
from pathlib import Path
from ast import literal_eval
from typing import Any, Callable, Dict, List, Optional, Union

import numpy as np
import pandas as pd
import yaml
from tqdm import tqdm

# Import local modules - use relative imports if possible
try:
    from residualize import calculate_residual
except ImportError:
    # Fallback for when script is run directly
    sys.path.append(str(Path(__file__).parent))
    from residualize import calculate_residual

# Configuration - use environment variables with defaults
PREFIX = os.getenv('PCQTL_PREFIX', '/home/klawren/oak/pcqtls')

# Constants for genomic analysis
BIDIRECTIONAL_PROMOTER_DISTANCE = 1000  # Base pairs

def setup_logging(level: int = logging.INFO) -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def get_residual_expression(covariates_path: str, expression_path: str) -> pd.DataFrame:
    """
    Calculate residualized expression data by removing covariate effects.
    
    Args:
        covariates_path (str): Path to covariates file
        expression_path (str): Path to expression data file
    
    Returns:
        pd.DataFrame: Residualized expression data
    """
    logger = logging.getLogger(__name__)
    logger.debug(f'Calculating residualized expression from {expression_path}')
    
    # Load covariates and expression data
    covariates_df = pd.read_table(covariates_path, index_col=0).T
    expression_df = pd.read_table(expression_path)
    
    logger.debug(f'Loaded {len(covariates_df)} covariates and {len(expression_df)} genes')
    
    # Set up expression data for residualization
    gid_expression = expression_df.set_index('gene_id')[covariates_df.index]
    
    # Calculate residuals
    residual_exp = calculate_residual(gid_expression[covariates_df.index], covariates_df, center=True)
    residual_exp = pd.DataFrame(residual_exp, columns=covariates_df.index, index=gid_expression.index)
    
    logger.debug(f'Calculated residuals for {len(residual_exp)} genes')
    return residual_exp

def annotate_correlation(cluster_df: pd.DataFrame, residual_exp: pd.DataFrame) -> pd.DataFrame:
    """
    Annotate clusters with correlation statistics.
    
    Args:
        cluster_df (pd.DataFrame): Cluster DataFrame to annotate
        residual_exp (pd.DataFrame): Residualized expression data
    """
    for idx, row in tqdm(cluster_df.iterrows(), total=len(cluster_df)):
        transcript_list = row['cluster_id'].split('_')
        
        # Calculate Spearman correlations between genes in cluster
        cluster_expression = residual_exp.loc[transcript_list].T.corr('spearman').to_numpy()
        cluster_corr = cluster_expression[np.triu_indices(len(cluster_expression), k=1)]
        
        # Store correlation statistics
        cluster_df.loc[idx, 'mean_corr'] = cluster_corr.mean()
        cluster_df.loc[idx, 'corr_list'] = str(cluster_corr)
        cluster_df.loc[idx, 'mean_pos_corr'] = cluster_corr[cluster_corr > 0].mean()
        cluster_df.loc[idx, 'mean_neg_corr'] = cluster_corr[cluster_corr < 0].mean()


def clusters_from_membership(cluster_df, annotation_df):
    """
    Annotate clusters based on gene pair membership in an annotation.
    
    Args:
        cluster_df: DataFrame with cluster_id column containing gene_ids seperated by '_'
        annotation_df: DataFrame with gene_id_1, gene_id_2 columns for pairs that should have the annotation
    """
    # Explode clusters to get individual genes
    cluster_df['gene_id_1'] = cluster_df['cluster_id'].str.split('_')
    cluster_df['gene_id_2'] = cluster_df['cluster_id'].str.split('_')

    clusters_exploded = cluster_df.explode('gene_id_1').explode('gene_id_2')
    clusters_exploded['gene_id_1'] = clusters_exploded['gene_id_1'].str.split('.').str[0]
    clusters_exploded['gene_id_2'] = clusters_exploded['gene_id_2'].str.split('.').str[0]
    clusters_exploded =  clusters_exploded[~(clusters_exploded['gene_id_1'] == clusters_exploded['gene_id_2'])]
 
    # Merge with annotation (check both directions)
    matched_clusters = set()
    
    # Check direct matches
    matches1 = pd.merge(clusters_exploded, annotation_df, on=['gene_id_1', 'gene_id_2'])
    matched_clusters.update(matches1['cluster_id'].unique())
    
    # Check reverse matches
    matches2 = pd.merge(clusters_exploded, annotation_df, 
                       left_on=['gene_id_1', 'gene_id_2'], 
                       right_on=['gene_id_2', 'gene_id_1'])
    matched_clusters.update(matches2['cluster_id'].unique())
    
    return matched_clusters


def clusters_from_boundaries(clusters_df, boundaries_df):
    """
    Annotate clusters with whether any boundary positions fall within start/end coordinates.
    
    Parameters:
    clusters_df: DataFrame with columns ['chr', 'start', 'end']
    boundaries_df: DataFrame with columns ['chr', 'pos']
    
    Returns:
    cluster_ids that cross boundaries
    """   
    merged = pd.merge(clusters_df, boundaries_df, on='chr', how='left')
    merged['pos_in_range'] = (merged['pos'] >= merged['start']) & (merged['pos'] <= merged['end'])
    return merged[merged['pos_in_range']]['cluster_id'].unique()



def annotate_cluster_edges(clusters_df, gencode_path):
    """
    Add start and end to clusters.

    Args:
        clusters: DataFrame with cluster_id column containing gene_ids seperated by '_'
        gencode: DataFrame with gene_id, chr, start, end columns

    Returns:
    """
    if 'start' in clusters_df.columns and 'end' in clusters_df.columns:
        print("Start and end already exist, nothing to do")
    else:
        gencode = pd.read_table(gencode_path)
        clusters_df['gene_id'] = clusters_df['cluster_id'].str.split('_')
        clusters_exploded = pd.merge(clusters_df.explode('gene_id')[['gene_id', 'cluster_id']], gencode, left_on='gene_id', right_on='gene_id', how='left')
        clusters_df = pd.merge(clusters_df, clusters_exploded.groupby('cluster_id').agg({'start': 'min', 'end': 'max'}), on='cluster_id')

    return clusters_df



def get_pairs_from_names(names_df):
    pairs = []
    for chr_id, names_chr in tqdm(names_df.groupby('chr')):
        # group by enhancer
        names_chr_grouped = names_chr.groupby('element_name').agg({'gene_id':list})
        names_chr_grouped = names_chr_grouped[names_chr_grouped['gene_id'].apply(len) > 1]

        # explode list ot pairs
        names_chr_grouped['gene_id_1'] = names_chr_grouped['gene_id']
        names_chr_grouped['gene_id_2'] = names_chr_grouped['gene_id']
        pairs_chr = names_chr_grouped.explode('gene_id_1').explode('gene_id_2')[['gene_id_1', 'gene_id_2']]
        pairs_chr = pairs_chr[pairs_chr['gene_id_1'] != pairs_chr['gene_id_2']]
        
        pairs.append(pairs_chr)

    pairs_df = pd.concat(pairs).reset_index(drop=True)[['gene_id_1', 'gene_id_2']].drop_duplicates()
    pairs_df['gene_id_1'] = pairs_df['gene_id_1'].str.split('.').str[0]
    pairs_df['gene_id_2'] = pairs_df['gene_id_2'].str.split('.').str[0]
    return pairs_df


def get_paralog_pairs(paralog_path: str) -> pd.DataFrame:
    paralogs = pd.read_table(paralog_path)
    paralogs.dropna(subset=['Human paralogue gene stable ID'], inplace=True)
    paralogs['gene_id_1'] = paralogs['Gene stable ID']
    paralogs['gene_id_2'] = paralogs['Human paralogue gene stable ID']

    paralog_pairs = paralogs[['gene_id_1', 'gene_id_2']].drop_duplicates()
    return paralog_pairs

def load_abc(abc_path, gencode_path, score_threshold=0.1):
    abc = pd.read_table(abc_path)
    abc = abc[(abc['ABC.Score'] > score_threshold) & (abc['class']!='promoter')]
    gencode = pd.read_table(gencode_path)
    # get the columns we need
    abc = abc[['chr', 'name', 'TargetGene', 'CellType']]
    # add gene ids 
    abc = pd.merge(abc, gencode[['gene_name', 'gene_id']], left_on='TargetGene', right_on='gene_name', how='left')
    # drop all that aren't protien coding
    abc.dropna(subset=['gene_name'], inplace=True)
    return abc

def get_abc_pairs(abc_path, gencode_path, score_threshold=0.1):
    abc = load_abc(abc_path, gencode_path, score_threshold)
    # get enhancers that are multi-gene
    abc_multi = abc[abc['name'].isin(abc[abc.duplicated('name')]['name'].unique())]
    abc_multi.drop_duplicates(subset=['gene_id', 'name'], inplace=True)
    abc_pairs = get_pairs_from_names(abc_multi.rename(columns={'name': 'element_name'}))
    return abc_pairs

def get_go_pairs(go_path, gencode_path):
    # Load GO terms
    go_df = pd.read_table(go_path, header=None,
                        names = ['Gene stable ID', 'Gene stable ID version', 'Gene start (bp)', 'Gene end', 'Strand', 
                                 'tss', 'gencode_annotation', 'gene_name', 'transcript_type', 'go_accession', 'go_name', 'go_domain'])
    
    go_df = go_df[go_df['go_domain'] == 'biological_process'].rename(columns={'go_accession': 'element_name', 'Gene stable ID': 'gene_id'})
    gencode = pd.read_table(gencode_path)
    gencode['gene_id'] = gencode['gene_id'].str.split('.').str[0]
    go_df = pd.merge(go_df, gencode[['gene_id', 'chr', 'strand']], on='gene_id', how='left').reset_index(drop=True)
    go_pairs = get_pairs_from_names(go_df[['gene_id', 'element_name', 'chr']])
    return go_pairs

def get_tad_boundries(tad_path):
    tad_df = pd.read_table(tad_path, names=['chr', 'start', 'end'])
    tad_boundries_df = pd.concat([tad_df[['chr', 'start']].rename(columns={'start': 'pos'}),
        tad_df[['chr', 'end']].rename(columns={'end': 'pos'})]) 
    return tad_boundries_df

def get_ctcf_boundaries(my_tissue_id, ctcf_match_path, ctcf_dir):
    # Load CTCF data
    ctcf_gtex_match = pd.read_csv(ctcf_match_path)
    ctcf_file = ctcf_gtex_match[ctcf_gtex_match['GTEX'] == my_tissue_id].iloc[0]['ctcf']
    
    ctcf_df = pd.read_table(f'{ctcf_dir}/{ctcf_file}.bed.gz', 
                           names=['chr', 'start', 'end', 'name', 'score', 'strand', 
                                 'signal_value', 'p_value', 'q_value', 'peak'])

    ctcf_boundries_df = pd.concat([ctcf_df[['chr', 'start']].rename(columns={'start': 'pos'}),
        ctcf_df[['chr', 'end']].rename(columns={'end': 'pos'})])

    return ctcf_boundries_df

def get_promoter_pairs(gencode_path, same_strand, dist_threshold=1000):
    gencode = pd.read_table(gencode_path)
    gencode['tss_pos'] = gencode['tss_pos'].apply(literal_eval)
    gencode_exploded = gencode.explode('tss_pos')
    # must be on opposite strands
    if same_strand:
        positive_orr = pd.merge(gencode[gencode['strand'] == '+'], gencode[gencode['strand'] == '+'], on='chr', how='left', suffixes=('_1', '_2'))
        negative_orr = pd.merge(gencode[gencode['strand'] == '-'], gencode[gencode['strand'] == '-'], on='chr', how='left', suffixes=('_1', '_2'))
        correct_orr = pd.concat([positive_orr, negative_orr])
    else:
        correct_orr = pd.merge(gencode[gencode['strand'] == '+'], gencode[gencode['strand'] == '-'], on='chr', how='left', suffixes=('_1', '_2'))
    # only consider genes close to each other
    correct_orr_close = correct_orr[np.where(correct_orr['end_1'] < correct_orr['start_2'], (correct_orr['start_2']-correct_orr['end_1']<dist_threshold), 
            np.where(correct_orr['end_2'] < correct_orr['start_1'], (correct_orr['start_1']-correct_orr['end_2']<dist_threshold), True))]
    # check if tss are within 1kb
    correct_orr_close_explode = correct_orr_close.explode('tss_pos_1').explode('tss_pos_2')
    correct_orr_close_explode = correct_orr_close_explode[correct_orr_close_explode['gene_id_1'] != correct_orr_close_explode['gene_id_2']]
    correct_orr_close_explode['shared_promoter'] = abs(correct_orr_close_explode['tss_pos_1'] - correct_orr_close_explode['tss_pos_2']) < dist_threshold
    promotor_pairs = correct_orr_close_explode.groupby(['gene_id_1', 'gene_id_2']).agg({'shared_promoter':'any'}).reset_index()
    promotor_pairs =  promotor_pairs[promotor_pairs['shared_promoter']][['gene_id_1', 'gene_id_2']]
    promotor_pairs['gene_id_1'] = promotor_pairs['gene_id_1'].str.split('.').str[0]
    promotor_pairs['gene_id_2'] = promotor_pairs['gene_id_2'].str.split('.').str[0]
    return promotor_pairs


def get_overlapping_pairs(gencode_path, same_strand):
    gencode = pd.read_table(gencode_path)
    gencode['tss_pos'] = gencode['tss_pos'].apply(literal_eval)
    gencode_exploded = gencode.explode('tss_pos')
    # must be on opposite strands
    if same_strand:
        positive_orr = pd.merge(gencode[gencode['strand'] == '+'], gencode[gencode['strand'] == '+'], on='chr', how='left', suffixes=('_1', '_2'))
        negative_orr = pd.merge(gencode[gencode['strand'] == '-'], gencode[gencode['strand'] == '-'], on='chr', how='left', suffixes=('_1', '_2'))
        correct_orr = pd.concat([positive_orr, negative_orr])
    else:
        correct_orr = pd.merge(gencode[gencode['strand'] == '+'], gencode[gencode['strand'] == '-'], on='chr', how='left', suffixes=('_1', '_2'))

    overlapping = correct_orr[(correct_orr['start_1'] <= correct_orr['end_2']) & (correct_orr['start_2'] <= correct_orr['end_1'])]

    overlapping = overlapping[overlapping['gene_id_1'] != overlapping['gene_id_2']]
    overlapping = overlapping[['gene_id_1', 'gene_id_2']].drop_duplicates()
    overlapping['gene_id_1'] = overlapping['gene_id_1'].str.split('.').str[0]
    overlapping['gene_id_2'] = overlapping['gene_id_2'].str.split('.').str[0]
    return overlapping


def get_cross_map_pairs(cross_map_path):
    cross_map_df = pd.read_table(cross_map_path, names=['gene_id_1', 'gene_id_2', 'score'])
    cross_map_pairs = cross_map_df[cross_map_df['score'] > 100]
    cross_map_pairs['gene_id_1'] = cross_map_pairs['gene_id_1'].str.split('.').str[0]
    cross_map_pairs['gene_id_2'] = cross_map_pairs['gene_id_2'].str.split('.').str[0]
    return cross_map_pairs   


# function to add all annotations given data paths
def add_annotations(cluster_df, tissue_id, gencode_path, paralog_path, abc_path, go_path, tad_path, ctcf_match_path, ctcf_dir, cross_map_path, covariates_path, expression_path):
    """
    Add all annotations to the clusters.        
    Args:
        cluster_df (pd.DataFrame): Cluster DataFrame to annotate
        tissue_id (str): GTEx tissue identifier
        gencode_path (str): Path to GENCODE annotation file
        paralog_path (str): Path to paralog data file
        abc_path (str): Path to ABC predictions file
        go_path (str): Path to GO term file
        tad_path (str): Path to TAD boundary file
        ctcf_match_path (str): Path to CTCF-GTEx tissue matching file
        ctcf_dir (str): Path to directory containing CTCF binding site files
        cross_map_path (str): Path to cross-mappability file
        covariates_path (str): Path to covariates file
        expression_path (str): Path to expression file
    Returns:
        pd.DataFrame: Annotated cluster DataFrame
    """
    logger = logging.getLogger(__name__)
    logger.info(f'Starting general annotation of {len(cluster_df)} clusters')

    logger.info('Adding start and end to clusters...')
    cluster_df = annotate_cluster_edges(cluster_df, gencode_path)

    logger.info('Loading paralog pairs...')
    paralog_pairs = get_paralog_pairs(paralog_path)
    logger.info('Adding paralog annotations...')
    paralog_clusters = clusters_from_membership(cluster_df, paralog_pairs)
    cluster_df['has_paralog'] = cluster_df['cluster_id'].isin(paralog_clusters)

    logger.info('Loading abc pairs...')
    abc_pairs = get_abc_pairs(abc_path, gencode_path)
    logger.info('Adding abc annotations...')
    abc_clusters = clusters_from_membership(cluster_df, abc_pairs)
    cluster_df['has_abc'] = cluster_df['cluster_id'].isin(abc_clusters)

    logger.info('Loading GO pairs...')
    go_pairs = get_go_pairs(go_path, gencode_path)
    logger.info('Adding GO annotations...')
    go_clusters = clusters_from_membership(cluster_df, go_pairs)
    cluster_df['has_go'] = cluster_df['cluster_id'].isin(go_clusters)

    logger.info('Loading TAD boundaries...')
    tad_boundries = get_tad_boundries(tad_path)
    logger.info('Adding TAD annotations...')
    tad_clusters = clusters_from_boundaries(cluster_df, tad_boundries)
    cluster_df['has_tad'] = cluster_df['cluster_id'].isin(tad_clusters)

    logger.info('Loading CTCF boundaries...')
    ctcf_boundries = get_ctcf_boundaries(tissue_id, ctcf_match_path, ctcf_dir)
    logger.info('Adding CTCF annotations...')
    ctcf_clusters = clusters_from_boundaries(cluster_df, ctcf_boundries)
    cluster_df['has_ctcf'] = cluster_df['cluster_id'].isin(ctcf_clusters)

    logger.info('Loading opp-strand (bidirectional) promoter pairs...')
    opp_promoter_pairs = get_promoter_pairs(gencode_path, False)
    logger.info('Adding opp-strand (bidirectional) promoter annotations...')
    opp_promoter_clusters = clusters_from_membership(cluster_df, opp_promoter_pairs)
    cluster_df['has_shared_opp_promoter'] = cluster_df['cluster_id'].isin(opp_promoter_clusters)

    logger.info('Loading same-strand promoter pairs...')
    same_promoter_pairs = get_promoter_pairs(gencode_path, True)
    logger.info('Adding same-strand promoter annotations...')
    same_promoter_clusters = clusters_from_membership(cluster_df, same_promoter_pairs)
    cluster_df['has_shared_same_promoter'] = cluster_df['cluster_id'].isin(same_promoter_clusters)

    logger.info('Loading opp-strand overlapping pairs...')
    opp_overlapping_pairs = get_overlapping_pairs(gencode_path, False)
    logger.info('Adding opp-strand overlapping annotations...')
    opp_overlapping_clusters = clusters_from_membership(cluster_df, opp_overlapping_pairs)
    cluster_df['has_opp_overlapping'] = cluster_df['cluster_id'].isin(opp_overlapping_clusters)

    logger.info('Loading same-strand overlapping pairs...')
    same_overlapping_pairs = get_overlapping_pairs(gencode_path, True)
    logger.info('Adding same-strand overlapping annotations...')
    same_overlapping_clusters = clusters_from_membership(cluster_df, same_overlapping_pairs)
    cluster_df['has_same_overlapping'] = cluster_df['cluster_id'].isin(same_overlapping_clusters)

    logger.info('Loading cross-mappability pairs...')
    cross_map_pairs = get_cross_map_pairs(cross_map_path)
    logger.info('Adding cross-mappability annotations...')
    cross_map_clusters = clusters_from_membership(cluster_df, cross_map_pairs)
    cluster_df['has_cross_map'] = cluster_df['cluster_id'].isin(cross_map_clusters)

    # Correlation annotation
    if 'mean_corr' not in cluster_df.columns:
        logger.info('Adding correlation annotations...')
        residual_exp = get_residual_expression(covariates_path, expression_path)
        annotate_correlation(cluster_df, residual_exp)

    logger.info(f'Completed general annotation of {len(cluster_df)} clusters')
    return cluster_df


def annotate_clusters_from_config(config: dict, tissue_id: str) -> pd.DataFrame:
    """
    Annotate clusters using paths derived from a config dict.
    Builds all paths by prefixing config['working_dir'] to configured paths.
    Returns the annotated DataFrame (no writing).
    """
    logger = logging.getLogger(__name__)
    wd = config['working_dir']

    clusters_path = f"{wd}/{config['clusters_dir']}/{tissue_id}.clusters.txt"
    expression_path = f"{wd}/{config['filtered_expression_output_dir']}/{tissue_id}.v8.normalized_residualized_expression.cluster_genes.bed"
    covariates_path = f"{wd}/{config['covariates_dir']}/{tissue_id}.v8.covariates.txt"
    gencode_path = f"{wd}/{config['gencode_path']}"
    abc_path = f"{wd}/{config['abc_path']}"
    ctcf_match_path = f"{wd}/{config['ctcf_match_path']}"
    ctcf_dir = f"{wd}/{config['ctcf_dir']}"
    paralog_path = f"{wd}/{config['paralog_path']}"
    go_path = f"{wd}/{config['go_path']}"
    cross_map_path = f"{wd}/{config['cross_map_path']}"
    tad_path = f"{wd}/{config['tad_path']}"

    cluster_df = pd.read_table(clusters_path)

    annotated = add_annotations(
        cluster_df,
        tissue_id,
        gencode_path,
        paralog_path,
        abc_path,
        go_path,
        tad_path,
        ctcf_match_path,
        ctcf_dir,
        cross_map_path,
        covariates_path,
        expression_path,
    )

    return annotated


def main() -> None:
    """
    Main function to handle command-line interface and execute cluster annotation.
    
    Parses command-line arguments and runs the cluster annotation pipeline.
    """
    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description='Annotate gene expression clusters with genomic and functional features',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Arguments in the same order as add_annotations signature
    # add_annotations(cluster_df, tissue_id, gencode_path, paralog_path, abc_path, go_path, tad_path,
    #                 ctcf_match_path, ctcf_dir, cross_map_path, covariates_path, expression_path)
    parser.add_argument('-t', '--tissue-id', required=True,
                       help='GTEx tissue identifier')
    parser.add_argument('-c', '--clusters', required=True,
                       help='Cluster CSV file with a cluster_id column')
    parser.add_argument('--gencode', required=True,
                       help='GENCODE annotation file')
    parser.add_argument('--paralog', required=True,
                       help='Paralog data file')
    parser.add_argument('--abc', required=True,
                       help='ABC predictions file')
    parser.add_argument('--go', required=True,
                       help='GO term file')
    parser.add_argument('--tad', required=True,
                       help='TAD boundary file')
    parser.add_argument('--ctcf-match', required=True,
                       help='CTCF-GTEx tissue matching file')
    parser.add_argument('--ctcf-dir', required=True,
                       help='Directory containing CTCF binding site files')
    parser.add_argument('--cross-map', required=True,
                       help='Cross-mappability file')
    parser.add_argument('-co', '--covariates', required=True,
                       help='Covariates file')
    parser.add_argument('-e', '--expression', required=True,
                       help='Normalized expression BED/TSV file')

    # Output and misc
    parser.add_argument('-o', '--output', required=True,
                       help='Output file for annotated clusters')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')

    # Parse arguments
    args = parser.parse_args()
    
    if args.verbose:
        setup_logging(logging.DEBUG)
        logger = logging.getLogger(__name__)
    
    try:
        # Validate required input files/directories
        logger.info('Validating input files...')
        required_paths = {
            'clusters': args.clusters,
            'gencode': args.gencode,
            'paralog': args.paralog,
            'abc': args.abc,
            'go': args.go,
            'tad': args.tad,
            'ctcf_match': args.ctcf_match,
            'cross_map': args.cross_map,
            'covariates': args.covariates,
            'expression': args.expression,
        }
        missing = [name for name, path in required_paths.items() if not path or not os.path.exists(path)]
        if args.ctcf_dir and not os.path.isdir(args.ctcf_dir):
            missing.append('ctcf_dir')
        if missing:
            raise FileNotFoundError(f"Missing or non-existent required inputs: {', '.join(sorted(set(missing)))}")

        logger.info(f"Starting cluster annotation for tissue: {args.tissue_id}")

        # Load clusters
        cluster_df = pd.read_table(args.clusters)

        # Call annotation using the exact parameter order of add_annotations
        cluster_df_annotated = add_annotations(
            cluster_df,
            args.tissue_id,
            args.gencode,
            args.paralog,
            args.abc,
            args.go,
            args.tad,
            args.ctcf_match,
            args.ctcf_dir,
            args.cross_map,
            args.covariates,
            args.expression,
        )

        # Write output
        logger.info(f"Writing annotated clusters to {args.output}")
        cluster_df_annotated.to_csv(args.output, index=False, sep='\t')
        logger.info(f"Successfully wrote {len(cluster_df_annotated)} annotated clusters to {args.output}")
        
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


