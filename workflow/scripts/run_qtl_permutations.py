#!/usr/bin/env python3
"""
QTL Permutation Analysis Script

This script runs permutation-based QTL analysis using TensorQTL to identify
independent associations and calculate empirical p-values for statistical significance.
"""

import argparse
import logging
import sys
import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis, genotypeio
from pathlib import Path
from typing import Dict, List, Optional, Union, Callable



def setup_logging(level: int = logging.INFO) -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def main() -> None:
    """
    Main function to run QTL permutation analysis.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Run permutation-based QTL analysis using TensorQTL',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--genotype-stem', required=True,
                       help='PLINK genotype file stem')
    parser.add_argument('--phenotype', required=True,
                       help='Path to phenotype BED file')
    parser.add_argument('--covariates', required=True,
                       help='Path to covariates file')
    parser.add_argument('--cis-results', required=True,
                       help='Path to cis QTL results file')
    parser.add_argument('--output', required=True,
                       help='Output file path')
    parser.add_argument('--tissue-id', required=True,
                       help='GTEx tissue identifier')

    
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
        logger.info('Starting QTL permutation analysis...')
        
        logger.info(f'Processing tissue: {args.tissue}')

        
        # Load phenotypes and covariates
        logger.info('Loading phenotype and covariate data...')
        phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(args.phenotype)
        covariates_df = pd.read_table(args.covariates, index_col=0).T
        
        logger.info(f'Loaded {phenotype_df.shape[1]} phenotypes and {covariates_df.shape[1]} covariates')
        
        # Load genotype data
        logger.info(f'Loading genotype data from: {args.genotype_stem}')
        pgr = genotypeio.PlinkReader(args.genotype_stem)
        genotype_df = pgr.load_genotypes()
        variant_df = pgr.bim.set_index('snp')[['chrom', 'pos']]
        
        logger.info(f'Loaded {genotype_df.shape[0]} variants and {genotype_df.shape[1]} samples')
        
        # Load cis results
        logger.info(f'Loading cis QTL results from: {args.cis_results}')
        cis_df = pd.read_table(args.cis_results, index_col=0)
        

        
        # Calculate q-values
        logger.info('Calculating q-values...')
        tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)
        
        # Run permutations
        logger.info('Running permutation analysis...')
        try:
            indep_df = cis.map_independent(
                genotype_df, variant_df, cis_df,
                phenotype_df, phenotype_pos_df, 
                covariates_df
            )
            
            # Write output
            logger.info(f'Writing {len(indep_df)} independent associations to: {args.output}')
            indep_df.to_csv(args.output, sep='\t', index=False, compression='gzip')
            logger.info('QTL permutation analysis complete!')
            
        except ValueError as e:
            logger.warning(f'Permutation analysis failed: {e}')
            logger.info('Writing empty results file...')
            
            # Create empty DataFrame with expected columns
            empty_df = pd.DataFrame(columns=[
                'phenotype_id', 'num_var', 'beta_shape1', 'beta_shape2', 'true_df',
                'pval_true_df', 'variant_id', 'start_distance', 'end_distance',
                'ma_samples', 'ma_count', 'af', 'pval_nominal', 'slope', 'slope_se',
                'pval_perm', 'pval_beta', 'rank'
            ])
            empty_df.to_csv(args.output, sep='\t', index=False, compression='gzip')
            logger.info('Empty results file written')
        
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