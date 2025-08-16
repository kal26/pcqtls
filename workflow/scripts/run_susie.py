#!/usr/bin/env python3
"""
SuSiE Fine-mapping Script

This script performs SuSiE (Sum of Single Effects) fine-mapping on QTL data using tensorQTL.
It identifies credible sets of variants that are likely to contain the causal variant
for each QTL association. 
"""

import logging
import sys
from pathlib import Path
from tensorqtl import susie, genotypeio
import tensorqtl
import numpy as np
import pandas as pd
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


def validate_input_files(genotype_path, phenotype_path, covariates_path):
    """Validate that input files exist and are readable."""
    files_to_check = [
        (genotype_path, "genotype"),
        (phenotype_path, "phenotype"),
        (covariates_path, "covariates")
    ]
    
    for file_path, file_type in files_to_check:
        if not Path(file_path).exists():
            raise FileNotFoundError(f"{file_type.capitalize()} file not found: {file_path}")
        if not Path(file_path).is_file():
            raise ValueError(f"{file_type.capitalize()} path is not a file: {file_path}")


def main():
    """
    Main function to perform SuSiE fine-mapping on QTL data.
    
    Loads genotype, phenotype, and covariate data, then runs SuSiE fine-mapping
    to identify credible sets of variants for each QTL association.
    """
    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description='Perform SuSiE fine-mapping on QTL data',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Define required arguments
    parser.add_argument('--genotype-stem', required=True,
                       help='Path to PLINK genotype files (without extension)')
    parser.add_argument('--expression', required=True,
                       help='Path to phenotype BED file (normalized expression or PCs)')
    parser.add_argument('--output', required=True,
                       help='Output path for SuSiE results')
    parser.add_argument('--covariates', required=True,
                       help='Path to covariates file')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')

    # Parse arguments
    args = parser.parse_args()
    
    if args.verbose:
        setup_logging(logging.DEBUG)

    try:
        # Validate input files
        logger.info("Validating input files...")
        validate_input_files(args.genotype_stem, args.expression, args.covariates)

        logger.info('Loading phenotype and covariate data...')
        # Load phenotypes and covariates
        phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(args.expression)
        covariates_df = pd.read_table(args.covariates, index_col=0).T
        
        logger.info(f'Loaded {len(phenotype_df)} phenotypes and {len(covariates_df)} covariates')

        logger.info('Loading genotype data...')
        # Load genotypes using PLINK reader
        pgr = genotypeio.PlinkReader(args.genotype_stem)
        genotype_df = pgr.load_genotypes()
        variant_df = pgr.bim.set_index('snp')[['chrom', 'pos']]
        
        logger.info(f'Loaded {len(genotype_df)} variants and {len(genotype_df.columns)} samples')

        logger.info('Running SuSiE fine-mapping...')
        # Perform SuSiE fine-mapping
        susie_df = susie.map(genotype_df, variant_df, 
                            phenotype_df, 
                            phenotype_pos_df, 
                            covariates_df)
        
        logger.info('Writing results...')
        logger.info(f'SuSiE results shape: {susie_df.shape}')

        # Handle empty results (no credible sets found)
        if len(susie_df) == 0:
            logger.warning('No SuSiE credible sets found - writing empty results file')
            # Create empty DataFrame with correct columns
            empty_df = pd.DataFrame(columns=['idx', 'phenotype_id', 'variant_id', 'pip', 'af', 'cs_id'])
            empty_df.to_csv(args.output, sep='\t', index=False)
        else:
            # Write results to file
            susie_df.to_csv(args.output, sep='\t', index=False)
            logger.info(f'SuSiE results written to: {args.output}')
            logger.info(f'Found {len(susie_df)} credible set variants')
            
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