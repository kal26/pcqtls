#!/usr/bin/env python3
"""
Principal Component Analysis (PCA) Generation Script

This script generates principal components (PCs) from gene expression clusters for QTL analysis.
It performs PCA on the expression data within each cluster, then residualizes the PCs to remove
the effects of covariates. The resulting PCs can be used as phenotypes in QTL mapping.
"""

import logging
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import argparse
from sklearn.decomposition import PCA
import os
from residualize import calculate_residual
from typing import Dict, List, Optional, Union, Callable


# Define BED format column order as a constant
BED_COLUMNS = ['#chr', 'start', 'end', 'gene_id']


def setup_logging(level: int = logging.INFO) -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def validate_input_files(cluster_path: str, expression_path: str, covariates_path: str) -> None:
    """Validate that input files exist and are readable."""
    files_to_check = [
        (cluster_path, "cluster"),
        (expression_path, "expression"),
        (covariates_path, "covariates")
    ]
    
    for file_path, file_type in files_to_check:
        if not Path(file_path).exists():
            raise FileNotFoundError(f"{file_type.capitalize()} file not found: {file_path}")
        if not Path(file_path).is_file():
            raise ValueError(f"{file_type.capitalize()} path is not a file: {file_path}")


def make_bed_order(df: pd.DataFrame) -> pd.DataFrame:
    """
    Reorder DataFrame columns to match BED file format requirements.
    
    BED files require specific column order: #chr, start, end, gene_id, followed by data columns.
    This function ensures the DataFrame has the correct column ordering.
    
    Args:
        df (pd.DataFrame): DataFrame with BED columns that need reordering
    
    Returns:
        pd.DataFrame: DataFrame with columns in BED format order
    """
    # Get all columns that are not in the required BED format
    data_columns = [col for col in df.columns if col not in BED_COLUMNS]
    
    # Create the proper column order: BED columns first, then data columns
    proper_order = BED_COLUMNS + data_columns
    
    # Reorder the DataFrame
    return df[proper_order]


def pc_bed_from_paths(cluster_path: str, expression_path: str, covariates_path: str, pc_out_path: str) -> None:
    """
    Generate principal components from cluster data and write to BED file.
    
    This function serves as the main entry point for PC generation from file paths.
    It loads the cluster data, expression data, and covariates, then generates
    principal components for each cluster and writes the results to a BED file.
    
    Args:
        cluster_path (str): Path to cluster CSV file
        expression_path (str): Path to expression BED file (already residualized)
        covariates_path (str): Path to covariates file
        pc_out_path (str): Output path for PC BED file
    
    Returns:
        None: Writes PC data to specified output file
    """
    logger = logging.getLogger(__name__)
    
    # Load input data
    logger.info('Loading input data...')
    cluster_df = pd.read_csv(cluster_path)
    logger.info(f'Loaded {len(cluster_df)} clusters from {cluster_path}')
    
    # Expression data is already residualized from previous processing
    expression_df = pd.read_table(expression_path)
    logger.info(f'Loaded expression data with {len(expression_df)} genes from {expression_path}')
    
    covariates_df = pd.read_table(covariates_path, index_col=0).T
    logger.info(f'Loaded {len(covariates_df)} covariates from {covariates_path}')

    # Generate principal components for all clusters
    logger.info('Generating principal components...')
    cluster_pcs_df = get_pc_bed(cluster_df, expression_df, covariates_df)

    # Write output BED file
    logger.info(f'Writing PC data to {pc_out_path}')
    cluster_pcs_df.to_csv(pc_out_path, sep='\t', index=False)
    logger.info(f'Successfully wrote {len(cluster_pcs_df)} PC phenotypes to {pc_out_path}')


def get_pc_bed(cluster_df: pd.DataFrame, expression_df: pd.DataFrame, covariates_df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate principal components from gene expression clusters.
    
    This function performs PCA on the expression data within each cluster to create
    principal components that capture the main patterns of co-expression. The PCs
    are then residualized to remove the effects of technical covariates.
    
    Args:
        cluster_df (pd.DataFrame): DataFrame containing cluster information
        expression_df (pd.DataFrame): Expression data in BED format with gene_id column
        covariates_df (pd.DataFrame): Covariates data for residualization
    
    Returns:
        pd.DataFrame: BED-formatted DataFrame with principal components as phenotypes
    """
    logger = logging.getLogger(__name__)
    
    # Extract sample IDs from expression data (columns after BED format columns)
    sample_ids = expression_df.columns[4:]
    logger.info(f'Processing {len(sample_ids)} samples')

    # Parse gene IDs to extract individual gene names and cluster IDs
    # Expression data has format: cluster_id_e_gene_id
    expression_df['egene_id'] = expression_df['gene_id'].str.split('_e_').str[1]
    expression_df['cluster_id'] = expression_df['gene_id'].str.split('_e_').str[0]
    expression_df_gid = expression_df.set_index('egene_id')
    
    logger.info(f'Parsed {len(expression_df_gid)} genes from expression data')

    # Generate PCs for each cluster
    cluster_pcs_dfs = []
    logger.info(f'Generating PCs for {len(cluster_df)} clusters...')
    
    for idx, row in cluster_df.iterrows():
        # Get cluster identifier
        cluster_id = row['cluster_id']
        logger.debug(f'Processing cluster {cluster_id} ({idx + 1}/{len(cluster_df)})')

        # Extract expression data for genes in this cluster
        cluster_expression_df_gid = expression_df_gid[expression_df_gid['cluster_id'] == cluster_id]
        cluster = cluster_expression_df_gid.loc[row['cluster_id'].split('_')]
        
        logger.debug(f'Cluster {cluster_id} contains {len(cluster)} genes')

        # Perform PCA on cluster expression data
        # Transpose to get samples x genes format for PCA
        X = cluster[sample_ids].transpose()
        pca = PCA()
        pc_values = pca.fit_transform(X)
        
        logger.debug(f'Generated {pc_values.shape[1]} PCs for cluster {cluster_id}')

        # Note: Could filter PCs by explained variance threshold
        # This was considered but not implemented as it left too few PCs
        # pc_values = pc_values[:, pca.explained_variance_ > 0.1]

        # Generate unique IDs for each principal component
        gene_ids = []
        for pc_num in range(pc_values.shape[1]):
            gene_ids.append(f"{row['cluster_id']}_pc{pc_num+1}")

        # Residualize the principal components to remove covariate effects
        # This ensures PCs are not confounded by technical factors
        logger.debug(f'Residualizing PCs for cluster {cluster_id}')
        normed_residualized_pcs = calculate_residual(
            pd.DataFrame(pc_values.T, columns=sample_ids), 
            covariates_df, 
            center=True
        )

        # Create DataFrame for this cluster's PCs
        cluster_pcs_df = pd.DataFrame(
            normed_residualized_pcs, 
            columns=sample_ids, 
            index=gene_ids
        )

        # Add BED format information (genomic coordinates)
        cluster_pcs_df = cluster_pcs_df.reset_index().rename(columns={'index': 'gene_id'})
        cluster_pcs_df['start'] = cluster['start'].min()
        cluster_pcs_df['end'] = cluster['end'].max()
        cluster_pcs_df['#chr'] = cluster['#chr'].iloc[0]
        
        # Ensure proper BED column ordering
        cluster_pcs_dfs.append(make_bed_order(cluster_pcs_df))

    # Combine all cluster PCs into single DataFrame
    cluster_pcs_df = pd.concat(cluster_pcs_dfs)
    logger.info(f'Combined {len(cluster_pcs_df)} PC phenotypes from all clusters')
    
    # Sort by genomic coordinates (required for downstream analysis)
    cluster_pcs_df = cluster_pcs_df.sort_values(['#chr', 'start', 'end'])

    # Clean up infinite values that can cause downstream analysis errors
    # Replace inf values with NaN and drop rows with any NaN values
    cluster_pcs_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    dropped_count = sum(cluster_pcs_df.isna().sum(axis=1) > 0)
    if dropped_count > 0:
        logger.warning(f'Dropped {dropped_count} rows due to infinite values')
    cluster_pcs_df.dropna(inplace=True)
    
    logger.info(f'Final dataset contains {len(cluster_pcs_df)} PC phenotypes')

    return cluster_pcs_df


def main() -> None:
    """
    Main function to handle command-line interface and execute PC generation.
    
    Parses command-line arguments and calls the PC generation pipeline.
    """
    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description='Generate principal components from gene expression clusters for QTL analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python get_pcs.py -cl clusters.csv -e expression.bed -co covariates.txt -o pcs.bed
  python get_pcs.py -cl clusters.csv -e expression.bed -co covariates.txt -o pcs.bed --verbose
        """
    )
    
    # Define command-line arguments
    parser.add_argument('-cl', '--clusters', required=True,
                       help='Cluster CSV file containing cluster information')
    parser.add_argument('-e', '--expression', required=True,
                       help='Expression BED file (should be pre-residualized)')
    parser.add_argument('-co', '--covariates', required=True,
                       help='Covariates file for residualization')
    parser.add_argument('-o', '--output', required=True,
                       help='Output file for PC BED')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')

    # Parse arguments and execute PC generation
    args = parser.parse_args()
    
    if args.verbose:
        setup_logging(logging.DEBUG)
    
    try:
        # Validate input files
        logger.info("Validating input files...")
        validate_input_files(args.cluster_path, args.expression_path, args.covariates_path)
        
        logger.info('Starting principal component generation...')
        pc_bed_from_paths(args.cluster_path, args.expression_path, 
                         args.covariates_path, args.pc_out_path)
        logger.info('Principal component generation complete!')
        
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