#!/usr/bin/env python3
"""
Principal Component Annotation Script

This script annotates principal components (PCs) derived from gene expression clusters
with individual gene expression data for downstream analysis.
"""

import argparse
import logging
import sys
import pandas as pd
import numpy as np
import ast
from tqdm import tqdm 
from scipy.stats import linregress
from pathlib import Path
from typing import List, Dict, Tuple


def setup_logging(level: int = logging.INFO) -> None:
    """Set up logging configuration.
    
    Args:
        level: Logging level (default: INFO)
    """
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def preprocess_pc(pc_path: str) -> pd.DataFrame:
    """Load and preprocess PC data.
    
    Args:
        pc_path: Path to PC data file
        
    Returns:
        DataFrame with PC data and derived columns
    """
    logger = logging.getLogger(__name__)
    logger.info(f'Loading PC data from: {pc_path}')
    
    pc_df = pd.read_table(pc_path)
    
    # Extract cluster and PC information
    pc_df['cluster_id'] = pc_df['gene_id'].str.split('_pc').str[0]
    pc_df['pc_id'] = pc_df['gene_id'].str.split('_pc').str[1].astype('float')
    pc_df['num_genes'] = pc_df['cluster_id'].str.split('_').apply(len)
    
    logger.info(f'Loaded {len(pc_df)} PC associations')
    logger.debug(f'PC data columns: {list(pc_df.columns)}')
    
    return pc_df





def preprocess_expression(pc_df: pd.DataFrame) -> pd.DataFrame:
    """Create expression dataframe for gene-level analysis.
    
    Args:
        pc_df: PC data DataFrame
        
    Returns:
        DataFrame with expression data structure
    """
    logger = logging.getLogger(__name__)
    logger.info('Creating expression dataframe for gene-level analysis...')
    
    expression_data = []
    for idx, row in pc_df.iterrows():
        cluster_id = row['cluster_id']
        for egene_id in cluster_id.split('_'):
            expression_data.append({
                'cluster_id': cluster_id,
                'egene_id': egene_id,
                'gene_id': f"{cluster_id}_e_{egene_id}"
            })
    
    expression_df = pd.DataFrame(expression_data)
    logger.info(f'Created expression dataframe with {len(expression_df)} entries')
    
    return expression_df


def calculate_gene_correlations(pc_values: np.ndarray, expression_values: np.ndarray) -> Tuple[float, float]:
    """Calculate correlation statistics between PC and gene expression.
    
    Args:
        pc_values: PC values for samples
        expression_values: Gene expression values for samples
        
    Returns:
        Tuple of (slope, r_squared)
    """
    slope, intercept, r_value, p_value, std_err = linregress(pc_values, expression_values)
    r_squared = r_value ** 2
    return slope, r_squared


def process_single_pc(pc_row: pd.Series, expression_df: pd.DataFrame, sample_ids: List[str]) -> Dict[str, List]:
    """Process a single PC to calculate gene-level correlations.
    
    Args:
        pc_row: Single row from PC DataFrame
        expression_df: Expression data DataFrame
        sample_ids: List of sample IDs
        
    Returns:
        Dictionary with gene slopes and variances
    """
    cluster_id = pc_row['cluster_id']
    expression_cluster = expression_df[expression_df['cluster_id'] == cluster_id].reset_index()
    
    gene_slopes = []
    gene_variances = []
    
    for egene_id in cluster_id.split('_'):
        # Get expression values for this gene
        gene_expression = expression_cluster[expression_cluster['egene_id'] == egene_id]
        if gene_expression.empty:
            # Handle case where gene is not found in expression data
            gene_slopes.append(np.nan)
            gene_variances.append(np.nan)
            continue
            
        expression_values = gene_expression[sample_ids].values.astype('float').flatten()
        pc_values = pc_row[sample_ids].astype('float').values
        
        # Calculate correlation
        slope, r_squared = calculate_gene_correlations(pc_values, expression_values)
        gene_slopes.append(slope)
        gene_variances.append(r_squared)
    
    return {
        'r2_list': str(gene_variances),
        'slope_list': str(gene_slopes)
    }


def annotate_pc_data(pc_df: pd.DataFrame, expression_df: pd.DataFrame) -> pd.DataFrame:
    """Annotate PC data with gene-level correlation statistics.
    
    Args:
        pc_df: PC data DataFrame
        expression_df: Expression data DataFrame
        
    Returns:
        Annotated PC DataFrame with gene-level statistics
    """
    logger = logging.getLogger(__name__)
    logger.info('Annotating PC data with gene-level correlations...')
    
    sample_ids = pc_df.columns[pc_df.columns.str.contains('GTEX')].tolist()
    
    # Process each PC
    for idx, row in tqdm(pc_df.iterrows(), total=pc_df.shape[0], desc="Processing PCs"):
        correlations = process_single_pc(row, expression_df, sample_ids)
        pc_df.loc[idx, 'r2_list'] = correlations['r2_list']
        pc_df.loc[idx, 'slope_list'] = correlations['slope_list']
    
    # Clean up and format results
    annotated_pc = pc_df.drop(columns=sample_ids)
    annotated_pc.rename(columns={'gene_id': 'pc_phenotype_id'}, inplace=True)
    
    # Create gene-level entries
    annotated_pc['egene_id'] = annotated_pc['cluster_id'].str.split('_')
    annotated_pc['egene_r2'] = annotated_pc['r2_list'].apply(ast.literal_eval)
    annotated_pc['egene_slope'] = annotated_pc['slope_list'].apply(ast.literal_eval)
    
    # Explode to create one row per gene
    result = annotated_pc.explode(['egene_id', 'egene_r2', 'egene_slope'])
    result = result.drop(columns=['r2_list', 'slope_list'])
    
    logger.info(f'Completed annotation with {len(result)} gene-level entries')
    return result








def main() -> None:
    """Main function to annotate principal components with gene expression data."""
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Annotate principal components with gene expression data',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--pcs', required=True,
                       help='Path to PC data file')
    parser.add_argument('--output', required=True,
                       help='Output path for annotated PCs')
    parser.add_argument('--tissue-id', required=True,
                       help='Tissue ID')
    
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
        # Load and process data
        pc_df = preprocess_pc(args.pcs)
        expression_df = preprocess_expression(pc_df)
        
        # Perform annotation
        annotated_pc = annotate_pc_data(pc_df, expression_df)
        annotated_pc['tissue_id'] = args.tissue
        
        # Save results
        logger.info(f'Saving {len(annotated_pc)} annotated PCs to: {args.output}')
        annotated_pc.to_csv(args.output, sep='\t', index=None)
        
        logger.info('PC annotation completed successfully!')
        
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

