#!/usr/bin/env python3
"""
Subsample Generation Script

This script creates subsamples of expression and covariate data for testing and validation.
It can either use a random scramble order or take the first N samples from the dataset.
"""

import argparse
import logging
import sys
import numpy as np
import pandas as pd
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


def main():
    """
    Main function to create subsamples of expression and covariate data.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Create subsamples of expression and covariate data for testing and validation',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--expression', required=True,
                       help='Path to expression data file')
    parser.add_argument('--covariates', required=True,
                       help='Path to covariates data file')
    parser.add_argument('--output-expression', required=True,
                       help='Output path for subsampled expression data')
    parser.add_argument('--output-covariates', required=True,
                       help='Output path for subsampled covariates data')
    parser.add_argument('--use-scramble', action='store_true',
                       help='Use random scramble order instead of first N samples')
    parser.add_argument('--num-samples', type=int,
                       help='Number of samples to select (default: all available)')
    
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
        logger.info('Loading expression and covariate data...')
        
        # Load data
        expression_df = pd.read_table(args.expression)
        covariates_df = pd.read_table(args.covariates, index_col=0).T
        
        # Determine number of samples to select
        if args.num_samples is None:
            num_samples = len(covariates_df)
        else:
            num_samples = min(args.num_samples, len(covariates_df))
        
        logger.info(f'Dataset has {len(covariates_df)} samples, selecting {num_samples}')
        logger.info(f'Using scramble order: {args.use_scramble}')
        
        # Select samples based on parameter
        if args.use_scramble:
            logger.info('Selecting random sample subset...')
            selected_samples = covariates_df.sample(num_samples).index.values
        else:
            logger.info('Selecting first N samples...')
            selected_samples = covariates_df.index[:num_samples].values
        
        # Create subsampled expression data
        logger.info('Creating subsampled expression data...')
        sub_expression = expression_df[np.concatenate([expression_df.columns[:4].values, selected_samples])]
        sub_expression.to_csv(args.output_expression, sep='\t', index=None)
        
        # Create subsampled covariate data
        logger.info('Creating subsampled covariate data...')
        sub_covar = covariates_df.loc[selected_samples]
        sub_covar.T.to_csv(args.output_covariates, sep='\t')
        
        logger.info(f'Successfully created subsamples with {len(selected_samples)} samples')
        logger.info(f'Expression subsample written to: {args.output_expression}')
        logger.info(f'Covariate subsample written to: {args.output_covariates}')
        
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
