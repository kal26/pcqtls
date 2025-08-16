#!/usr/bin/env python3
"""
SuSiE VEP Annotation and Merging Script

This script merges SuSiE fine-mapping results with VEP annotations and nominal QTL data.
It combines eQTL and pcQTL SuSiE results, adds variant effect predictions, and creates
a comprehensive dataset for downstream analysis.

"""

import argparse
import logging
import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path

# Import local modules for annotation functionality
try:
    from annotate_qtls import load_and_annotate
except ImportError:
    # Fallback for when script is run directly
    import sys
    sys.path.append(str(Path(__file__).parent))
    from annotate_qtls import load_and_annotate


def setup_logging(level=logging.INFO):
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def load_susie_data(e_susie_path, pc_susie_path):
    """
    Load and combine eQTL and pcQTL SuSiE results.
    
    Args:
        e_susie_path (str): Path to eQTL SuSiE results
        pc_susie_path (str): Path to pcQTL SuSiE results
    
    Returns:
        pd.DataFrame: Combined SuSiE results
    """
    logger = logging.getLogger(__name__)
    logger.info('Loading SuSiE fine-mapping results...')
    
    pc_susie_df = pd.read_table(pc_susie_path, index_col=0)
    pc_susie_df['cluster_id'] = pc_susie_df['phenotype_id'].str.split('_pc').str[0]
    
    e_susie_df = pd.read_table(e_susie_path, index_col=0)
    e_susie_df['cluster_id'] = e_susie_df['phenotype_id'].str.split('_e').str[0]
    
    combined_susie = pd.concat([e_susie_df, pc_susie_df], 
                              names=['type', 'idx'], 
                              keys=['eqtl', 'pcqtl']).reset_index(drop=0).drop(columns=['idx'])
    
    try:
        combined_susie.drop(columns=['af'], inplace=True)
        logger.debug('Removed AF column from SuSiE results')
    except KeyError:
        logger.debug('No AF column to drop (R SuSiE version)')
    
    return combined_susie


def load_nominal_data(e_nominal_path_list, pc_nominal_path_list):
    """
    Load nominal QTL data from multiple files.
    
    Args:
        e_nominal_path_list (list): List of eQTL nominal result paths
        pc_nominal_path_list (list): List of pcQTL nominal result paths
    
    Returns:
        tuple: (e_nominal, pc_nominal) DataFrames
    """
    logger = logging.getLogger(__name__)
    logger.info('Loading nominal QTL data...')
    
    e_nominal = pd.concat([pd.read_parquet(path) for path in e_nominal_path_list])
    pc_nominal = pd.concat([pd.read_parquet(path) for path in pc_nominal_path_list])
    
    return e_nominal, pc_nominal





def annotate_susie_results(e_susie_path, pc_susie_path, vep_path, annot_pc_path, 
                          e_nominal_path_list, pc_nominal_path_list, tissue_id,
                          gencode_path, full_abc_path, abc_match_path, ctcf_match_path,
                          ctcf_dir, tad_path):
    """
    Annotate SuSiE fine-mapping results with VEP predictions and functional data.
    
    This function processes SuSiE fine-mapping results from both eQTL and pcQTL analyses,
    merges them with VEP variant effect predictions, and combines with nominal QTL data
    to create a comprehensive dataset for downstream analysis.
    
    Args:
        e_susie_path (str): Path to eQTL SuSiE results
        pc_susie_path (str): Path to pcQTL SuSiE results
        vep_path (str): Path to VEP annotations
        annot_pc_path (str): Path to annotated PC data
        e_nominal_path_list (list): List of eQTL nominal result paths
        pc_nominal_path_list (list): List of pcQTL nominal result paths
        tissue_id (str): GTEx tissue identifier
        gencode_path (str): Path to GENCODE file
        full_abc_path (str): Path to ABC predictions file
        abc_match_path (str): Path to ABC-GTEx matching file
        ctcf_match_path (str): Path to CTCF-GTEx matching file
        ctcf_dir (str): Directory containing CTCF files
        tad_path (str): Path to TAD boundaries file

    
    Returns:
        pd.DataFrame: Annotated SuSiE results
    """
    logger = logging.getLogger(__name__)
    logger.info('Starting SuSiE annotation...')
    
    # Load data using direct file paths
    combined_susie = load_susie_data(e_susie_path, pc_susie_path)
    e_nominal, pc_nominal = load_nominal_data(e_nominal_path_list, pc_nominal_path_list)
    
    # Merge SuSiE with nominal data and add QTL variance
    logger.info('Merging SuSiE results with nominal QTL data...')
    qtls_nominal_merged = pd.merge(
        combined_susie, 
        pd.concat([pc_nominal, e_nominal]), 
        left_on=['phenotype_id', 'variant_id'],  
        right_on=['phenotype_id', 'variant_id'], 
        how='left'
    )
    qtls_nominal_merged['qtl_variance'] = qtls_nominal_merged['slope'].apply(np.square) * 100
    qtls_nominal_merged = qtls_nominal_merged.rename(columns={'slope': 'qtl_slope', 'slope_se': 'qtl_slope_se'})
    
    # Add eQTL nominal information aggregated by cluster and variant
    logger.info('Adding eQTL nominal information...')
    e_nominal_sub = e_nominal[e_nominal['variant_id'].isin(qtls_nominal_merged['variant_id'])]
    e_nominal_sub['cluster_id'] = e_nominal_sub['phenotype_id'].str.split('_e').str[0]
    e_nominal_sub['egene_id'] = e_nominal_sub['phenotype_id'].str.split('_e_').str[1]
    e_nominal_sub['variance'] = e_nominal_sub['slope'].apply(np.square) * 100
    
    egene_nominal = e_nominal_sub.groupby(['cluster_id', 'variant_id']).agg({
        'variance': list, 
        'egene_id': list, 
        'slope': list, 
        'slope_se': list
    })
    egene_nominal = egene_nominal.rename(columns={
        'variance': 'egene_variance_list', 
        'egene_id': 'egene_id_list', 
        'slope': 'egene_qtl_slope', 
        'slope_se': 'egene_qtl_slope_se'
    })
    
    qtls_nominal_merged = pd.merge(
        qtls_nominal_merged, 
        egene_nominal, 
        left_on=['cluster_id', 'variant_id'], 
        right_index=True, 
        how='left'
    )
    
    # Add PC annotations and sign flip QTL slopes
    logger.info('Adding PC annotations...')
    annotated_pcs = pd.read_table(annot_pc_path)
    annotated_pcs = annotated_pcs.rename(columns={
        'egene_r2': 'egene_pc_r2', 
        'egene_slope': 'egene_pc_slope'
    })
    
    qtl_annot_pc_merged = pd.merge(
        qtls_nominal_merged.explode(['egene_id_list', 'egene_qtl_slope', 'egene_variance_list', 'egene_qtl_slope_se']), 
        annotated_pcs[['pc_phenotype_id', 'egene_id', 'egene_pc_r2', 'egene_pc_slope']], 
        right_on=['pc_phenotype_id', 'egene_id'], 
        left_on=['phenotype_id', 'egene_id_list'], 
        how='left'
    )
    qtl_annot_pc_merged = qtl_annot_pc_merged.drop(columns=['pc_phenotype_id', 'egene_id'])
    qtl_annot_pc_merged['egene_qtl_slope_flipped'] = qtl_annot_pc_merged['egene_qtl_slope'] * qtl_annot_pc_merged['egene_pc_slope']
    
    # Regroup data by variant to consolidate eQTL information
    logger.info('Regrouping data by variant...')
    qtls_nominal_merged = qtl_annot_pc_merged.groupby(['phenotype_id', 'variant_id', 'cs_id']).agg({
        'type': 'first',
        'pip': 'first', 
        'cluster_id': 'first',
        'start_distance': 'first',
        'end_distance': 'first', 
        'af': 'first',
        'ma_samples': 'first',
        'ma_count': 'first',
        'pval_nominal': 'first',
        'qtl_slope': 'first',
        'qtl_slope_se': 'first',
        'qtl_variance': 'first',
        'egene_variance_list': list,
        'egene_id_list': list,
        'egene_qtl_slope': list,
        'egene_qtl_slope_se': list,
        'egene_pc_r2': list,
        'egene_pc_slope': list,
        'egene_qtl_slope_flipped': list
    }).reset_index()
    
    # Add VEP annotations
    logger.info('Adding VEP variant effect predictions...')
    vep = pd.read_table(vep_path, skiprows=4)
    qtls = pd.merge(
        qtls_nominal_merged, 
        vep[['ID', 'INFO', 'POS', '#CHROM']], 
        left_on='variant_id', 
        right_on='ID', 
        how='left'
    ).drop(columns=['ID'])
    
    qtls.rename(columns={
        'INFO': 'vep_info', 
        'POS': 'position', 
        '#CHROM': 'chr'
    }, inplace=True)
    qtls['tissue_id'] = tissue_id
    
    # Add additional QTL annotations using existing utilities
    logger.info('Adding genomic and functional annotations...')
    load_and_annotate(qtls, tissue_id, gencode_path, full_abc_path, abc_match_path,
                     ctcf_match_path, ctcf_dir, tad_path)
    
    logger.info(f'Completed SuSiE annotation with {len(qtls)} QTLs')
    return qtls


def main():
    """
    Main function to merge SuSiE results with VEP annotations.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Merge SuSiE fine-mapping results with VEP annotations and nominal QTL data',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--e-susie', required=True,
                       help='Path to eQTL SuSiE results file')
    parser.add_argument('--pc-susie', required=True,
                       help='Path to pcQTL SuSiE results file')
    parser.add_argument('--vep', required=True,
                       help='Path to VEP annotations file')
    parser.add_argument('--annot-pc', required=True,
                       help='Path to annotated PC data file')
    parser.add_argument('--e-nominal', nargs='+', required=True,
                       help='List of eQTL nominal result files')
    parser.add_argument('--pc-nominal', nargs='+', required=True,
                       help='List of pcQTL nominal result files')
    parser.add_argument('--tissue-id', required=True,
                       help='GTEx tissue identifier')
    parser.add_argument('--gencode', required=True,
                       help='Path to GENCODE annotation file')
    parser.add_argument('--full-abc', required=True,
                       help='Path to ABC predictions file')
    parser.add_argument('--abc-match', required=True,
                       help='Path to ABC-GTEx tissue matching file')
    parser.add_argument('--ctcf-match', required=True,
                       help='Path to CTCF-GTEx tissue matching file')
    parser.add_argument('--ctcf-dir', required=True,
                       help='Directory containing CTCF binding site files')
    parser.add_argument('--tad', required=True,
                       help='Path to TAD boundary file')

    parser.add_argument('--output', required=True,
                       help='Output file path')
    
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
        logger.info('Starting SuSiE VEP annotation and merging...')
        
        # Call the annotation function
        qtls = annotate_susie_results(
            e_susie_path=args.e_susie,
            pc_susie_path=args.pc_susie,
            vep_path=args.vep,
            annot_pc_path=args.annot_pc,
            e_nominal_path_list=args.e_nominal,
            pc_nominal_path_list=args.pc_nominal,
            tissue_id=args.tissue_id,
            gencode_path=args.gencode,
            full_abc_path=args.full_abc,
            abc_match_path=args.abc_match,
            ctcf_match_path=args.ctcf_match,
            ctcf_dir=args.ctcf_dir,
            tad_path=args.tad,

        )
        
        # Write output
        logger.info(f'Writing {len(qtls)} annotated QTLs to {args.output}')
        qtls.to_csv(args.output, sep='\t', index=False)
        logger.info('SuSiE VEP annotation complete!')
        
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


