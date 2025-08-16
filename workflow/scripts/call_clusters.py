#!/usr/bin/env python3
"""
Gene Expression Cluster Detection Script

This script identifies clusters of co-expressed neighboring genes within chromosomes using
correlation analysis on residualized expression data. It supports both value-based
and p-value-based correlation thresholds with Bonferroni correction.
"""

import logging
import sys
import pandas as pd
import numpy as np
from residualize import calculate_residual
from scipy.stats import spearmanr
import argparse


def setup_logging(level=logging.INFO):
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def validate_input_files(expression_path, covariates_path):
    """Validate that input files exist and are readable."""
    from pathlib import Path
    
    files_to_check = [
        (expression_path, "expression"),
        (covariates_path, "covariates")
    ]
    
    for file_path, file_type in files_to_check:
        if not Path(file_path).exists():
            raise FileNotFoundError(f"{file_type.capitalize()} file not found: {file_path}")
        if not Path(file_path).is_file():
            raise ValueError(f"{file_type.capitalize()} path is not a file: {file_path}")


def get_clusters_chr(chr_id, expression_df, residual_exp, total_pairs, tissue_id, 
                    min_cluster_size=2, max_cluster_size=50, min_corr_cutoff=.01, 
                    percent_corr_cutoff=.7, cutoff_type='pvalue', trim=False):
    """
    Identify gene expression clusters on a per-chromosome basis.
    
    Args:
        chr_id (int): Chromosome number
        expression_df (pd.DataFrame): Expression data with gene_id and #chr columns
        residual_exp (pd.DataFrame): Residualized expression data
        total_pairs (int): Total number of gene pairs for Bonferroni correction
        tissue_id (str): Tissue identifier
        min_cluster_size (int): Minimum number of genes per cluster
        max_cluster_size (int): Maximum number of genes per cluster
        min_corr_cutoff (float): Minimum correlation threshold
        percent_corr_cutoff (float): Minimum percentage of correlations above threshold
        cutoff_type (str): 'pvalue' for p-value based or 'value' for correlation value based
        trim (bool): Whether to trim cluster edges with no correlations
    
    Returns:
        pd.DataFrame: DataFrame with columns: chr, cluster_id, tissue_id, num_genes, 
                     percent_correlated, mean_corr, mean_pos_corr, mean_neg_corr
    """
    logger = logging.getLogger(__name__)
    logger.debug(f'Processing chromosome {chr_id}')
    
    # Get genes on this chromosome
    chr_gene_ids = expression_df[expression_df['#chr'] == f'chr{chr_id}']['gene_id']
    chr_residual_exp = residual_exp.loc[chr_gene_ids]
    
    logger.debug(f'Chromosome {chr_id} contains {len(chr_gene_ids)} genes')

    # Calculate correlation matrix and p-values
    try:
        chr_corr, chr_pvalue = spearmanr(chr_residual_exp, axis=1)
    except IndexError:
        # No genes on this chromosome
        logger.warning(f'No genes found on chromosome {chr_id}')
        return pd.DataFrame({})
    
    chr_corr = pd.DataFrame(chr_corr, index=chr_residual_exp.index, columns=chr_residual_exp.index)
    chr_pvalue = pd.DataFrame(chr_pvalue, index=chr_residual_exp.index, columns=chr_residual_exp.index)

    logger.info(f'Calculated correlations for chromosome {chr_id}')

    # Track transcripts already assigned to clusters
    transcript_blacklist = []
    cluster_output = []

    # Iterate through cluster sizes from largest to smallest
    for cluster_size in np.arange(max_cluster_size, min_cluster_size-1, -1):
        logger.debug(f'Checking clusters of size {cluster_size} on chromosome {chr_id}')
        
        # Total number of pairs considered for this cluster size
        number_pairs = sum(sum(np.triu(np.ones((cluster_size, cluster_size)), k=1)))

        # Check each possible cluster start position
        for cluster_start_idx in range(0, len(chr_corr)-cluster_size):
            # Extract cluster candidate matrix
            cluster_candidate = chr_corr.iloc[cluster_start_idx:cluster_start_idx+cluster_size, 
                                            cluster_start_idx:cluster_start_idx+cluster_size]
            cluster_values = cluster_candidate.values[np.triu_indices(cluster_size, k=1)]
            
            # Count correlations above threshold
            if cutoff_type == 'value':
                number_corrs_above_cutoff = sum(abs(cluster_values) > min_corr_cutoff)
            elif cutoff_type == 'pvalue':
                # Calculate p-values below Bonferroni-corrected threshold
                cluster_candidate_pvalues = chr_pvalue.iloc[cluster_start_idx:cluster_start_idx+cluster_size, 
                                                          cluster_start_idx:cluster_start_idx+cluster_size]
                cluster_pvalues = cluster_candidate_pvalues.values[np.triu_indices(cluster_size, k=1)]
                number_corrs_above_cutoff = sum(cluster_pvalues < (0.05/total_pairs))
            else:
                raise ValueError(f"Invalid cutoff_type: {cutoff_type}")

            # Check if transcripts are already in clusters
            cluster_transcripts = cluster_candidate.index.values
            in_blacklist = sum(cluster_candidate.index.isin(transcript_blacklist))
            percent_corr = number_corrs_above_cutoff/number_pairs

            # Record cluster if criteria are met
            if (percent_corr > percent_corr_cutoff) & (in_blacklist == 0):
                # Initialize trimmed cluster as full cluster
                if cutoff_type == 'value':
                    cluster_bool_df = abs(cluster_candidate) > min_corr_cutoff
                elif cutoff_type == 'pvalue':
                    cluster_bool_df = cluster_candidate_pvalues < (0.05/total_pairs)

                trimmed_cluster = cluster_candidate
                trimmed_cluster_bool = cluster_bool_df

                # Trim cluster edges if requested
                if trim:
                    # Trim genes from start with no correlations
                    first_gene_inclusion = sum(cluster_bool_df.iloc[1:,0])
                    while first_gene_inclusion == 0:
                        trimmed_cluster_bool = trimmed_cluster_bool.iloc[1:, 1:]
                        trimmed_cluster = trimmed_cluster.iloc[1:,1:]
                        first_gene_inclusion = sum(trimmed_cluster_bool.iloc[1:,0])

                    # Trim genes from end with no correlations
                    last_gene_inclusion = sum(trimmed_cluster_bool.iloc[-1, :-1])
                    while last_gene_inclusion == 0:
                        trimmed_cluster_bool = trimmed_cluster_bool.iloc[:-1, :-1]
                        trimmed_cluster = trimmed_cluster.iloc[:-1, :-1]
                        last_gene_inclusion = sum(trimmed_cluster_bool.iloc[-1, :-1])

                # Update cluster information after trimming
                trimmed_cluster_size = len(trimmed_cluster)
                cluster_transcripts = trimmed_cluster.index.values
                number_corrs_above_cutoff = sum(trimmed_cluster_bool.values[np.triu_indices(trimmed_cluster_size, k=1)])
                percent_corr = number_corrs_above_cutoff/sum(sum(np.triu(np.ones((trimmed_cluster_size, trimmed_cluster_size)), k=1)))
                trimmed_cluster_values = trimmed_cluster.values[np.triu_indices(trimmed_cluster_size, k=1)]
                
                # Create cluster record
                cluster_series = pd.Series({
                    'chr': f'chr{chr_id}',
                    'cluster_id': '_'.join(sorted(cluster_transcripts)),
                    'tissue_id': tissue_id,
                    'num_genes': trimmed_cluster_size,
                    'percent_correlated': percent_corr,
                    'mean_corr': np.mean(trimmed_cluster_values),
                    'mean_pos_corr': np.mean(trimmed_cluster_values[trimmed_cluster_values > 0]),
                    'mean_neg_corr': np.mean(trimmed_cluster_values[trimmed_cluster_values < 0])
                })
                
                cluster_output.append(cluster_series)
                # Add transcripts to blacklist to prevent overlap
                [transcript_blacklist.append(id) for id in cluster_transcripts]

    # Create output DataFrame
    out_df = pd.DataFrame(cluster_output)
    logger.info(f'Found {len(cluster_output)} clusters on chromosome {chr_id}')
    
    if len(cluster_output) == 0:
        # Return empty DataFrame with correct columns
        return pd.DataFrame(columns=['chr', 'cluster_id', 'tissue_id', 'num_genes', 
                                   'percent_correlated', 'mean_corr', 'mean_pos_corr', 'mean_neg_corr'])
    
    if trim:
        # Remove clusters smaller than minimum after trimming
        out_df = out_df[out_df['num_genes'] >= min_cluster_size]
    
    # Sort by cluster size (largest first)
    return out_df.sort_values('num_genes', ascending=False)


def get_clusters_from_paths(expression_path, covariates_path, tissue_id, 
                           min_cluster_size=2, max_cluster_size=50, min_corr_cutoff=0.1, 
                           percent_corr_cutoff=.7, cutoff_type='pvalue', trim=True):
    """
    Load data from files and identify gene expression clusters.
    
    Args:
        expression_path (str): Path to expression data file (.bed format)
        covariates_path (str): Path to covariates file
        tissue_id (str): Tissue identifier
        min_cluster_size (int): Minimum number of genes per cluster
        max_cluster_size (int): Maximum number of genes per cluster
        min_corr_cutoff (float): Minimum correlation threshold
        percent_corr_cutoff (float): Minimum percentage of correlations above threshold
        cutoff_type (str): 'pvalue' for p-value based or 'value' for correlation value based
        trim (bool): Whether to trim cluster edges with no correlations
    
    Returns:
        pd.DataFrame: DataFrame with columns: chr, cluster_id, tissue_id, num_genes, 
                     percent_correlated, mean_corr, mean_pos_corr, mean_neg_corr
    """
    # Load expression and covariates data
    expression_df = pd.read_table(expression_path)
    logger.info(f'Loaded expression data: {expression_df.shape}')
    logger.debug(expression_df.head())
    
    covariates_df = pd.read_table(covariates_path, index_col=0).T
    logger.info(f'Loaded covariates data: {covariates_df.shape}')
    logger.debug(covariates_df.head())
    logger.info('Data loading complete')

    # Residualize expression data
    residual_exp = calculate_residual(expression_df[covariates_df.index], covariates_df, center=True)
    residual_exp = pd.DataFrame(residual_exp, columns=covariates_df.index, index=expression_df['gene_id'])
    logger.info('Expression data residualized')

    return get_clusters(expression_df, residual_exp, tissue_id, 
                       min_cluster_size=min_cluster_size, max_cluster_size=max_cluster_size, 
                       min_corr_cutoff=min_corr_cutoff, percent_corr_cutoff=percent_corr_cutoff, 
                       cutoff_type=cutoff_type, trim=trim)


def get_clusters(expression_df, residual_exp, tissue_id, min_cluster_size=2, 
                max_cluster_size=50, min_corr_cutoff=0.1, percent_corr_cutoff=.7, 
                cutoff_type='pvalue', trim=True):
    """
    Identify gene expression clusters across all chromosomes.
    
    Args:
        expression_df (pd.DataFrame): Expression data with gene_id and #chr columns
        residual_exp (pd.DataFrame): Residualized expression data
        tissue_id (str): Tissue identifier
        min_cluster_size (int): Minimum number of genes per cluster
        max_cluster_size (int): Maximum number of genes per cluster
        min_corr_cutoff (float): Minimum correlation threshold
        percent_corr_cutoff (float): Minimum percentage of correlations above threshold
        cutoff_type (str): 'pvalue' for p-value based or 'value' for correlation value based
        trim (bool): Whether to trim cluster edges with no correlations
    
    Returns:
        pd.DataFrame: DataFrame with columns: chr, cluster_id, tissue_id, num_genes, 
                     percent_correlated, mean_corr, mean_pos_corr, mean_neg_corr
    """
    # Calculate total number of pairs for Bonferroni correction
    total_pairs = 0
    for i in range(1, 23):
        chr_gene_ids = expression_df[expression_df['#chr'] == f'chr{i}']['gene_id']
        upper_corner_idxs = np.triu(np.ones(len(chr_gene_ids)), k=1)
        excluded_cluster_size_idxs = np.triu(np.ones(len(chr_gene_ids)), k=max_cluster_size)
        total_pairs += upper_corner_idxs.sum() - excluded_cluster_size_idxs.sum()

    logger.info(f'Total pairs for Bonferroni correction: {total_pairs:,}')

    # Process each chromosome
    clusters_all_chr = []
    for i in range(1, 23):
        logger.info(f'Processing chromosome {i}...')
        clusters_all_chr.append(get_clusters_chr(i, expression_df, residual_exp, total_pairs, 
                                               tissue_id, min_cluster_size=min_cluster_size, 
                                               max_cluster_size=max_cluster_size, 
                                               min_corr_cutoff=min_corr_cutoff, 
                                               percent_corr_cutoff=percent_corr_cutoff, 
                                               cutoff_type=cutoff_type, trim=trim))

    # Combine results from all chromosomes
    return pd.concat(clusters_all_chr)


def main():
    """Main function to parse arguments and run cluster detection."""
    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    parser = argparse.ArgumentParser(
        description="Gene Expression Cluster Detection Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('-e', '--expression', required=True,
                       help='Normalized expression BED file (not residualized)')
    parser.add_argument('-co', '--covariates', required=True,
                       help='Covariates file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output file for cluster results')
    
    # Optional arguments
    parser.add_argument('--max-cluster-size', type=int, default=50,
                       help='Maximum number of genes per cluster (default: 50)')
    parser.add_argument('--min-cluster-size', type=int, default=2,
                       help='Minimum number of genes per cluster (default: 2)')
    parser.add_argument('--min-corr-cutoff', type=float, default=.01,
                       help='Minimum correlation value threshold (default: 0.01)')
    parser.add_argument('--percent-corr-cutoff', type=float, default=.7,
                       help='Minimum percentage of correlations above threshold (default: 0.7)')
    parser.add_argument('--cutoff-type', type=str, default='pvalue', 
                       choices=['pvalue', 'value'],
                       help='Correlation threshold type: pvalue (Bonferroni-corrected) or value (default: pvalue)')
    parser.add_argument('--tissue-id', type=str, default='Unknown',
                       help='Tissue identifier for output (default: Unknown)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')

    args = parser.parse_args()
    
    if args.verbose:
        setup_logging(logging.DEBUG)
    
    try:
        # Validate input files
        logger.info("Validating input files...")
        validate_input_files(args.expression_path, args.covariates_path)
        
        logger.info('Starting gene expression cluster detection...')
        logger.info(f'Tissue: {args.tissue_id}')
        logger.info(f'Expression file: {args.expression_path}')
        logger.info(f'Covariates file: {args.covariates_path}')
        logger.info(f'Output file: {args.out_path}')
        logger.info(f'Parameters: min_size={args.min_cluster_size}, max_size={args.max_cluster_size}, '
              f'cutoff_type={args.cutoff_type}, percent_cutoff={args.percent_corr_cutoff}')

        # Run cluster detection
        clusters_all_chr = get_clusters_from_paths(
            args.expression_path, args.covariates_path, args.tissue_id, 
            args.min_cluster_size, args.max_cluster_size, args.min_corr_cutoff, 
            args.percent_corr_cutoff, args.cutoff_type
        )
        
        clusters_all_chr.reset_index(drop=True, inplace=True)
        
        # Write results
        logger.info(f'Writing {len(clusters_all_chr)} clusters to {args.out_path}')
        clusters_all_chr.to_csv(args.out_path, index=False)
        logger.info(f'Cluster detection complete! Found {len(clusters_all_chr)} total clusters')
        logger.info(f'Results saved to: {args.out_path}')
        
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