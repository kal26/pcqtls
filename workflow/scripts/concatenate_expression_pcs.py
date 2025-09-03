#!/usr/bin/env python3
"""
Concatenate eQTL and pcQTL expression BED files for combined analysis.

This script combines individual gene expression data with principal component
data to create a unified phenotype file for combined QTL analysis.
"""

import argparse
import logging
import pandas as pd
import sys
from pathlib import Path


def setup_logging() -> logging.Logger:
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)


def load_bed_file(file_path: str, logger: logging.Logger) -> pd.DataFrame:
    """
    Load a BED file and return as DataFrame.
    
    Args:
        file_path (str): Path to BED file
        logger (logging.Logger): Logger instance
        
    Returns:
        pd.DataFrame: Loaded BED file data
    """
    try:
        logger.info(f'Loading BED file: {file_path}')
        df = pd.read_csv(file_path, sep='\t')
        logger.info(f'Loaded {len(df)} rows from {file_path}')
        return df
    except Exception as e:
        logger.error(f'Failed to load {file_path}: {e}')
        sys.exit(1)


def concatenate_expression_pcs(expression_df: pd.DataFrame, pcs_df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """
    Concatenate expression and PC BED files.
    
    Args:
        expression_df (pd.DataFrame): Expression BED data
        pcs_df (pd.DataFrame): Principal components BED data
        logger (logging.Logger): Logger instance
        
    Returns:
        pd.DataFrame: Combined BED data
    """
    logger.info('Concatenating expression and PC data...')
    
    # Ensure both DataFrames have the same structure
    required_columns = ['#chr', 'start', 'end', 'gene_id']
    
    # Check if required columns exist
    for col in required_columns:
        if col not in expression_df.columns:
            logger.error(f'Required column {col} not found in expression data')
            sys.exit(1)
        if col not in pcs_df.columns:
            logger.error(f'Required column {col} not found in PC data')
            sys.exit(1)
    
    # Get sample columns (all columns after the required BED columns)
    expression_samples = [col for col in expression_df.columns if col not in required_columns]
    pc_samples = [col for col in pcs_df.columns if col not in required_columns]
    
    logger.info(f'Expression data has {len(expression_samples)} samples')
    logger.info(f'PC data has {len(pc_samples)} samples')
    
    # Verify sample columns match between datasets
    if set(expression_samples) != set(pc_samples):
        logger.error('Sample columns do not match between expression and PC data')
        logger.error(f'Expression samples: {sorted(expression_samples[:5])}...')
        logger.error(f'PC samples: {sorted(pc_samples[:5])}...')
        sys.exit(1)
    
    # Concatenate the DataFrames
    combined_df = pd.concat([expression_df, pcs_df], ignore_index=True)
    
    # Sort by genomic coordinates (required for downstream analysis)
    combined_df = combined_df.sort_values(['#chr', 'start', 'end'])
    
    logger.info(f'Combined dataset contains {len(combined_df)} phenotypes')
    logger.info(f'Expression phenotypes: {len(expression_df)}')
    logger.info(f'PC phenotypes: {len(pcs_df)}')
    
    return combined_df


def main() -> None:
    """Main function to concatenate expression and PC BED files."""
    parser = argparse.ArgumentParser(
        description='Concatenate eQTL expression and pcQTL PC BED files for combined analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--expression', required=True,
                       help='Path to expression BED file')
    parser.add_argument('--pcs', required=True,
                       help='Path to PC BED file')
    parser.add_argument('--output', required=True,
                       help='Output path for combined BED file')
    
    args = parser.parse_args()
    logger = setup_logging()
    
    # Validate input files
    for file_path in [args.expression, args.pcs]:
        if not Path(file_path).exists():
            logger.error(f'Input file does not exist: {file_path}')
            sys.exit(1)
    
    # Load input data
    expression_df = load_bed_file(args.expression, logger)
    pcs_df = load_bed_file(args.pcs, logger)
    
    # Concatenate data
    combined_df = concatenate_expression_pcs(expression_df, pcs_df, logger)
    
    # Write output
    logger.info(f'Writing combined data to: {args.output}')
    combined_df.to_csv(args.output, sep='\t', index=False)
    logger.info('Successfully created combined BED file')


if __name__ == '__main__':
    main() 