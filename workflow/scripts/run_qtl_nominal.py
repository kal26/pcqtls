#!/usr/bin/env python3
"""
QTL Nominal Analysis Script

This script runs nominal QTL analysis using TensorQTL to identify associations
between genetic variants and gene expression phenotypes (either individual genes
or principal components). A wrapper is needed so that an empty file is written 
if there are no phenotypes on a chromosome.
"""

import argparse
import logging
import sys
import pandas as pd
import subprocess
import os
from pathlib import Path


def setup_logging(level=logging.INFO):
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def run_tensorqtl_nominal(genotype_stem: str, expression: str, output_dir: str, 
                         tissue: str, covariates: str, phenotype_type: str) -> None:
    """
    Run TensorQTL nominal analysis.
    
    Args:
        genotype_stem: PLINK genotype file stem
        expression: Path to expression BED file
        output_dir: Output directory for QTL results
        tissue: GTEx tissue identifier
        covariates: Path to covariates file
        phenotype_type: Type of phenotype ('.v8.cluster_genes' or '.v8.pcs')
    """
    logger = logging.getLogger(__name__)
    
    # Determine phenotype type description
    if phenotype_type == '.v8.cluster_genes':
        qtl_type = 'eQTL'
        phenotype_desc = 'cluster genes'
    elif phenotype_type == '.v8.pcs':
        qtl_type = 'pcQTL'
        phenotype_desc = 'principal components'
    else:
        raise ValueError(f"Invalid phenotype_type: {phenotype_type}. Must be '.v8.cluster_genes' or '.v8.pcs'")
    
    logger.info(f'Running TensorQTL nominal analysis for {qtl_type} ({phenotype_desc})...')
    
    # Construct the TensorQTL command
    output_prefix = f"{output_dir}{tissue}/{tissue}{phenotype_type}"
    nominal_command = f"python -m tensorqtl {genotype_stem} {expression} {output_prefix} --covariates {covariates} --mode cis_nominal"
    
    try:
        # Run the command
        result = subprocess.run(nominal_command, shell=True, check=True, text=True, capture_output=True)
        logger.info(f"TensorQTL {qtl_type} command completed successfully")
        logger.debug(f"Command output: {result.stdout}")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"TensorQTL {qtl_type} command failed: {e}")
        logger.error(f"Error output: {e.stderr}")
        # Continue to create missing files even if command fails


def create_missing_files(output_dir: str, tissue: str, phenotype_type: str) -> None:
    """
    Create missing parquet files for chromosomes with no phenotypes.
    
    Args:
        output_dir: Output directory for QTL results
        tissue: GTEx tissue identifier
        phenotype_type: Type of phenotype ('.v8.cluster_genes' or '.v8.pcs')
    """
    logger = logging.getLogger(__name__)
    logger.info('Checking for missing output files...')
    
    chr_range = range(1, 23)
    
    # Define the base filename
    base_filename = f'{output_dir}{tissue}/{tissue}{phenotype_type}.cis_qtl_pairs.chr'
    
    # Expected columns for TensorQTL output
    expected_columns = [
        "phenotype_id", "num_var", "beta_shape1", "beta_shape2", "true_df", 
        "pval_true_df", "variant_id", "start_distance", "end_distance", 
        "af", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se", 
        "pval_perm", "pval_beta", "beta", "beta_se"
    ]
    
    # Create empty DataFrame with expected columns
    empty_df = pd.DataFrame(columns=expected_columns)
    
    # Check each chromosome
    for chr_num in chr_range:
        expected_file = f"{base_filename}{chr_num}.parquet"
        
        if not os.path.exists(expected_file):
            logger.info(f"Creating empty file for chromosome {chr_num}: {expected_file}")
            
            # Ensure output directory exists
            os.makedirs(os.path.dirname(expected_file), exist_ok=True)
            
            # Write empty parquet file
            empty_df.to_parquet(expected_file, index=False)
            logger.info(f"Created empty parquet file: {expected_file}")
        else:
            logger.debug(f"File already exists for chromosome {chr_num}: {expected_file}")


def main():
    """
    Main function to run QTL nominal analysis.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Run nominal QTL analysis using TensorQTL for either eQTL or pcQTL',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--genotype-stem', required=True,
                       help='PLINK genotype file stem')
    parser.add_argument('--expression', required=True,
                       help='Path to expression BED file')
    parser.add_argument('--covariates', required=True,
                       help='Path to covariates file')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for QTL results')
    parser.add_argument('--tissue-id', required=True,
                       help='GTEx tissue identifier')
    parser.add_argument('--phenotype-type', required=True, 
                       choices=['.v8.cluster_genes', '.v8.pcs'],
                       help='Type of phenotype: .v8.cluster_genes for eQTL or .v8.pcs for pcQTL')
    
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
        # Determine QTL type for logging
        qtl_type = 'eQTL' if args.phenotype_type == '.v8.cluster_genes' else 'pcQTL'
        logger.info(f'Starting {qtl_type} nominal analysis...')
        
        logger.info(f'Processing tissue: {args.tissue}')
        logger.info(f'Genotype stem: {args.genotype_stem}')
        logger.info(f'Expression file: {args.expression}')
        logger.info(f'Output directory: {args.output_dir}')
        logger.info(f'Phenotype type: {args.phenotype_type}')
        
        # Run TensorQTL nominal analysis
        run_tensorqtl_nominal(
            args.genotype_stem, 
            args.expression, 
            args.output_dir, 
            args.tissue, 
            args.covariates, 
            args.phenotype_type
        )
        
        # Create missing files for chromosomes with no phenotypes
        create_missing_files(args.output_dir, args.tissue, args.phenotype_type)
        
        logger.info(f'{qtl_type} nominal analysis completed successfully!')
        
    except Exception as e:
        logger.error(f"Error during {qtl_type} nominal analysis: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 