#!/usr/bin/env python3
"""
SuSiE to VCF Conversion Script

This script converts SuSiE fine-mapping results to VCF format for downstream analysis.
It extracts variant information from SuSiE output and formats it according to VCF specifications.
"""

import argparse
import logging
import sys
import pandas as pd
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
    Main function to convert SuSiE results to VCF format.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Convert SuSiE fine-mapping results to VCF format for downstream analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--e-susie', required=True,
                       help='Path to eQTL SuSiE results file')
    parser.add_argument('--pc-susie', required=True,
                       help='Path to pcQTL SuSiE results file')
    parser.add_argument('--output', required=True,
                       help='Output VCF file path')
    
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
        logger.info('Loading SuSiE fine-mapping results...')
        
        # Load SuSiE data
        pc_susie_df = pd.read_table(args.pc_susie, index_col=0)
        e_susie_df = pd.read_table(args.e_susie, index_col=0)
        
        logger.info(f'Loaded {len(pc_susie_df)} pcQTL and {len(e_susie_df)} eQTL SuSiE associations')
        
        # Get all variants in credible sets
        logger.info('Extracting unique variants from credible sets...')
        susie_vars = pd.DataFrame(pd.concat([e_susie_df, pc_susie_df])['variant_id'].drop_duplicates())
        
        # Parse variant information
        logger.info('Parsing variant information...')
        susie_vars['chr'] = susie_vars['variant_id'].str.split('_').str[0].str[3:]
        susie_vars['pos'] = susie_vars['variant_id'].str.split('_').str[1]
        susie_vars['ref'] = susie_vars['variant_id'].str.split('_').str[2]
        susie_vars['alt'] = susie_vars['variant_id'].str.split('_').str[3]
        susie_vars['blank'] = '.'
        
        # Create VCF format columns
        col_list = ['chr', 'pos', 'variant_id', 'ref', 'alt', 'blank', 'blank', 'blank']
        out_df = susie_vars[col_list]
        
        # Sort by chromosome and position
        out_df = out_df.sort_values(by=['chr', 'pos'])
        
        # Write VCF output
        logger.info(f'Writing {len(out_df)} variants to VCF format: {args.output}')
        out_df.to_csv(args.output, header=None, index=False, sep='\t')
        
        logger.info('SuSiE to VCF conversion complete!')
        
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