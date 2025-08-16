#!/usr/bin/env python3
"""
Expression Cluster Filtering Script

This script filters expression data to include only genes that are part of identified
clusters. It residualizes the expression data to remove covariate effects, then
creates cluster-based expression files where each gene is renamed to include the
cluster identifier. This is used to prepare data for downstream QTL analysis on
cluster genes.
"""

import argparse
import logging
import sys
from typing import List, Tuple
import pandas as pd
import numpy as np
from residualize import calculate_residual


def setup_logging(level=logging.INFO):
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def load_data(cluster_path: str, expression_path: str, covariates_path: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load cluster, expression, and covariate data.
    
    Args:
        cluster_path: Path to cluster CSV file
        expression_path: Path to expression BED file
        covariates_path: Path to covariates file
        
    Returns:
        Tuple of (cluster_df, expression_df, covariates_df)
    """
    logger = logging.getLogger(__name__)
    logger.info('Loading cluster and expression data...')
    
    cluster_df = pd.read_csv(cluster_path)
    logger.info(f'Loaded {len(cluster_df)} clusters from {cluster_path}')
    
    expression_df = pd.read_table(expression_path)
    logger.info(f'Loaded {len(expression_df)} genes from {expression_path}')
    
    covariates_df = pd.read_table(covariates_path, index_col=0).T
    logger.info(f'Loaded {len(covariates_df)} covariates from {covariates_path}')
    
    return cluster_df, expression_df, covariates_df


def residualize_expression(expression_df: pd.DataFrame, covariates_df: pd.DataFrame) -> pd.DataFrame:
    """
    Residualize expression data to remove covariate effects.
    
    Args:
        expression_df: Expression data in BED format
        covariates_df: Covariates data
        
    Returns:
        Residualized expression data
    """
    logger = logging.getLogger(__name__)
    logger.info('Residualizing expression data...')
    
    # Set gene_id as index for easier manipulation
    expression_df_gid = expression_df.set_index('gene_id')

    # Residualize the expression data to remove covariate effects
    residual_exp = calculate_residual(expression_df[expression_df.columns[4:]], covariates_df, center=True)
    expression_df_res = expression_df_gid.copy()
    expression_df_res[expression_df.columns[4:]] = residual_exp
    
    logger.info(f'Residualized expression data for {len(expression_df_res)} genes')
    return expression_df_res


def create_cluster_expression_files(cluster_df: pd.DataFrame, expression_df_res: pd.DataFrame) -> pd.DataFrame:
    """
    Create cluster-based expression files.
    
    Args:
        cluster_df: Cluster information DataFrame
        expression_df_res: Residualized expression data
        
    Returns:
        Combined cluster expression data
    """
    logger = logging.getLogger(__name__)
    logger.info('Creating cluster-based expression files...')
    
    # Filter to just the genes that are in the clusters
    # Change gene_id to be the cluster_id_e_gene_id for each gene
    # Expand the search window to be the same size for all genes in the cluster
    cluster_sets: List[pd.DataFrame] = []
    
    logger.info(f'Processing {len(cluster_df)} clusters...')
    for idx, row in cluster_df.iterrows():
        # Get expression data for genes in this cluster
        cluster_set = expression_df_res.loc[row['cluster_id'].split('_')].copy()
        
        logger.debug(f'Cluster {row["cluster_id"]} contains {len(cluster_set)} genes')
        
        # Use the same genomic coordinates for all genes in the cluster
        cluster_set['start'] = cluster_set['start'].min()
        cluster_set['end'] = cluster_set['end'].max()
        
        # Rename genes to include cluster identifier
        cluster_set.insert(3, 'gene_id', [f"{row['cluster_id']}_e_{gid}" for gid in cluster_set.index])
        cluster_sets.append(cluster_set)

    # Combine all cluster expression data
    out_df = pd.concat(cluster_sets)
    logger.info(f'Combined expression data from {len(cluster_sets)} clusters')
    
    # Sort by genomic coordinates (required by tensorQTL)
    out_df = out_df.sort_values(by=['#chr', 'start'])
    
    return out_df


def main():
    """Main function to process cluster expression data."""
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Filter and residualize expression data for cluster-based eQTL analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--clusters', required=True,
                       help='Path to cluster CSV file')
    parser.add_argument('--expression', required=True,
                       help='Path to expression BED file')
    parser.add_argument('--covariates', required=True,
                       help='Path to covariates file')
    parser.add_argument('--output', required=True,
                       help='Output path for filtered expression data')
    
    # Optional arguments
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Set up logging
    if args.verbose:
        setup_logging(logging.DEBUG)
    else:
        setup_logging()
    
    logger = logging.getLogger(__name__)
    
    try:
        # Load input data
        cluster_df, expression_df, covariates_df = load_data(
            args.clusters, 
            args.expression, 
            args.covariates
        )
        
        # Residualize expression data
        expression_df_res = residualize_expression(expression_df, covariates_df)
        
        # Create cluster-based expression files
        out_df = create_cluster_expression_files(cluster_df, expression_df_res)
        
        # Write output
        logger.info('Writing filtered expression data...')
        out_df.to_csv(args.output, sep='\t', index=False)
        logger.info(f'Filtered expression data written to: {args.output}')
        logger.info(f'Total genes in clusters: {len(out_df)}')
        
    except Exception as e:
        logger.error(f"Error processing cluster expression data: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
