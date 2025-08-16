#!/usr/bin/env python3
"""
QTL Annotation Script

This script annotates QTL associations with various genomic and functional features
including ABC enhancer-gene connections, CTCF binding sites, TAD boundaries,
bidirectional promoters, and expression statistics for downstream analysis.
"""

import logging
import sys
import os
import numpy as np
import pandas as pd
import ast
from tqdm import tqdm 
from scipy.stats import linregress
from pathlib import Path

# Import local modules
try:
    from utils import *
    from annotate_clusters import *
except ImportError:
    # Fallback for when script is run directly
    import sys
    sys.path.append(str(Path(__file__).parent))
    from utils import *
    from annotate_clusters import *

# Configuration - use environment variables with defaults
PREFIX = os.getenv('PCQTL_PREFIX', '/home/klawren/oak/pcqtls')

# Constants
BIDIRECTIONAL_PROMOTER_DISTANCE = 1000  # Base pairs


def setup_logging(level=logging.INFO):
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )





def annotate_enhancers_qtl(qtls, gene_enhancer_df):
    """
    Annotate QTLs with ABC enhancer information.
    
    Args:
        qtls (pd.DataFrame): QTL associations to annotate
        gene_enhancer_df (pd.DataFrame): ABC enhancer-gene connections
    """
    logger = logging.getLogger(__name__)
    logger.info('Annotating QTLs with ABC enhancer information...')
    
    for idx, qtl_row in tqdm(qtls.iterrows(), total=len(qtls)):
        # Handle both chr prefix and numeric chromosome formats
        chr_val = qtl_row['chr'] if isinstance(qtl_row['chr'], str) else f'chr{qtl_row["chr"]}'
        matched_enhancers = gene_enhancer_df[
            (gene_enhancer_df['chr'] == chr_val) & 
            (gene_enhancer_df['enhancer_start'] <= qtl_row['position']) & 
            (gene_enhancer_df['enhancer_end'] >= qtl_row['position'])
        ]
        
        # Is QTL in any ABC enhancer?
        qtls.loc[idx, 'qtl_num_abc_enhancers'] = matched_enhancers['enhancer'].nunique()
        # How many genes does the enhancer contact?
        qtls.loc[idx, 'qtl_num_abc_genes'] = matched_enhancers['TargetGene'].nunique()
        # How many cluster genes?
        qtls.loc[idx, 'qtl_matched_abc_genes'] = pd.Series(matched_enhancers.reset_index()['gene_id'].unique()).isin(qtl_row['cluster_id'].split('_')).sum()


def annotate_ctcf_tad_qtl(qtls, ctcf_df, tad_df, gid_gencode):
    """
    Annotate QTLs with CTCF binding sites and TAD boundary information.
    
    Args:
        qtls (pd.DataFrame): QTL associations to annotate
        ctcf_df (pd.DataFrame): CTCF binding site data
        tad_df (pd.DataFrame): TAD boundary data
        gid_gencode (pd.DataFrame): GENCODE annotations
    """
    logger = logging.getLogger(__name__)
    logger.info('Annotating QTLs with CTCF and TAD information...')
    
    # Create temporary TSS min/max for CTCF annotation
    qtls['tss_min'] = qtls['position'].astype(float)
    qtls['tss_max'] = qtls['position'].astype(float) + 1
    
    # Ensure chr column has consistent format
    qtls['chr'] = qtls['chr'].apply(lambda x: x if isinstance(x, str) else f'chr{x}')
    
    # Annotate CTCF binding sites
    annotate_ctcf(qtls, ctcf_df)
    qtls['qtl_in_ctcf'] = qtls['has_ctcf_peak']
    
    # Redo TSS min and max for proper position annotation
    annotate_positions(qtls, gid_gencode)
    
    # Add TAD annotations
    qtls['qtl_inter'] = pd.arrays.IntervalArray.from_arrays(
        qtls['position'].astype(float), 
        qtls['position'].astype(float) + 5000  # 10kb TAD resolution
    )
    qtls['num_tads_qtl'] = qtls.apply(count_tad_overlap, axis=1, args=(tad_df, 'qtl_inter'))
    qtls['qtl_in_tad'] = qtls['num_tads_qtl'] > 1
    qtls['between_tss'] = ((qtls['tss_min'] < qtls['position']) & (qtls['tss_max'] > qtls['position']))
    qtls['qtl_in_tss_ctcf'] = qtls['between_tss'] & qtls['qtl_in_ctcf']
    qtls['qtl_in_tad_ctcf'] = qtls['qtl_in_tad'] & qtls['qtl_in_ctcf']


def annotate_bidirectional_qtl(qtls, gid_gencode):
    """
    Annotate QTLs with bidirectional and shared promoter information.
    
    Args:
        qtls (pd.DataFrame): QTL associations to annotate
        gid_gencode (pd.DataFrame): GENCODE annotations
    """
    logger = logging.getLogger(__name__)
    logger.info('Annotating QTLs with bidirectional promoter information...')
    
    qtls['in_bidirectional_promoter'] = False
    qtls['in_shared_promoter'] = False
    
    for qtl_idx, row in tqdm(qtls.iterrows(), total=len(qtls)):
        gene_ids = row['cluster_id'].split('_')
        cluster_gencode = gid_gencode.loc[gene_ids]
        
        # Check if QTL is in the promoter (within 1kb of TSS)
        for idx, first_gene_row in cluster_gencode.iterrows():
            for idx, second_gene_row in cluster_gencode.iterrows():
                opp_strand = first_gene_row['strand'] != second_gene_row['strand']
                close = abs(first_gene_row['tss_start'] - second_gene_row['tss_start']) <= BIDIRECTIONAL_PROMOTER_DISTANCE
                
                if opp_strand & close:
                    # Found a bidirectional promoter
                    if ((row['position'] - first_gene_row['tss_start']) < BIDIRECTIONAL_PROMOTER_DISTANCE) & ((row['position'] - second_gene_row['tss_start']) < BIDIRECTIONAL_PROMOTER_DISTANCE):
                        qtls.loc[qtl_idx, 'in_bidirectional_promoter'] = True
                elif close:
                    # Found a shared promoter
                    if ((row['position'] - first_gene_row['tss_start']) < BIDIRECTIONAL_PROMOTER_DISTANCE) & ((row['position'] - second_gene_row['tss_start']) < BIDIRECTIONAL_PROMOTER_DISTANCE):
                        qtls.loc[qtl_idx, 'in_shared_promoter'] = True


def add_annotations_qtl(qtls, gid_gencode, gene_enhancer_df, ctcf_df, tad_df):
    """
    Add all QTL annotations to the DataFrame.
    
    Args:
        qtls (pd.DataFrame): QTL associations to annotate
        gid_gencode (pd.DataFrame): GENCODE annotations
        gene_enhancer_df (pd.DataFrame): ABC enhancer-gene connections
        ctcf_df (pd.DataFrame): CTCF binding site data
        tad_df (pd.DataFrame): TAD boundary data
    """
    logger = logging.getLogger(__name__)
    logger.info(f'Starting QTL annotation of {len(qtls)} associations')
    
    qtls.reset_index(drop=True, inplace=True)
    
    # Add enhancer annotations
    logger.info('Adding enhancer annotations...')
    annotate_enhancers_qtl(qtls, gene_enhancer_df)
    
    # Add CTCF and TAD annotations
    logger.info('Adding CTCF and TAD annotations...')
    annotate_ctcf_tad_qtl(qtls, ctcf_df, tad_df, gid_gencode)
    
    # Add bidirectional promoter annotations
    logger.info('Adding bidirectional promoter annotations...')
    annotate_bidirectional_qtl(qtls, gid_gencode)
    
    # Add distance annotations
    logger.info('Adding distance annotations...')
    annotate_distance(qtls, gid_gencode)
    
    logger.info(f'Completed QTL annotation of {len(qtls)} associations')


def load_and_annotate(qtls, my_tissue_id,
                      gencode_path='data/references/processed_gencode.v26.GRCh38.genes.txt', 
                      full_abc_path='data/references/functional_annotations/ABC_predictions/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz', 
                      abc_match_path='data/references/functional_annotations/ABC_predictions/ABC_matched_gtex.txt', 
                      ctcf_match_path='data/references/functional_annotations/ctcf_chip/ctcf_matched_gtex.txt', 
                      ctcf_dir='data/references/functional_annotations/ctcf_chip', 
                      tad_path='data/references/TAD_annotations/TADs_hg38/converted_HiC_IMR90_DI_10kb.txt'):
    """
    Load all annotation data and annotate QTLs.
    
    Args:
        qtls (pd.DataFrame): QTL associations to annotate
        my_tissue_id (str): GTEx tissue identifier
        gencode_path (str): Path to GENCODE file
        full_abc_path (str): Path to ABC predictions file
        abc_match_path (str): Path to ABC-GTEx matching file
        ctcf_match_path (str): Path to CTCF-GTEx matching file
        ctcf_dir (str): Directory containing CTCF files
        tad_path (str): Path to TAD boundaries file

    """
    logger = logging.getLogger(__name__)
    logger.info(f'Loading annotation data for tissue: {my_tissue_id}')
    
    # Load all annotation data
    
    logger.info('Loading GENCODE annotations...')
    gid_gencode, full_gencode = load_gencode(f'{PREFIX}/{gencode_path}')
    
    logger.info('Loading ABC enhancer-gene connections...')
    gene_enhancer_df = load_abc(my_tissue_id, full_gencode, f'{PREFIX}/{full_abc_path}', f'{PREFIX}/{abc_match_path}')
    
    logger.info('Loading CTCF binding sites...')
    ctcf_df = load_ctcf(my_tissue_id, f'{PREFIX}/{ctcf_match_path}', f'{PREFIX}/{ctcf_dir}')
    
    logger.info('Loading TAD boundaries...')
    tad_df = load_tad(f'{PREFIX}/{tad_path}')
    
    logger.info('All annotation data loaded successfully')
    
    # Annotate QTLs
    add_annotations_qtl(qtls, gid_gencode, gene_enhancer_df, ctcf_df, tad_df)


def annotate_distance(qtls, gid_gencode):
    qtls['position'] = qtls['variant_id'].str.split('_').str[1].astype(int)
    qtls['cluster_min_distance'] = qtls.apply(get_tss, axis=1, args=(gid_gencode,))


# distance to whichever gene in the cluster is closest
def get_tss(row, gid_gencode):
    cluster_gene_df = gid_gencode.loc[row['cluster_id'].split('_')]
    starts = cluster_gene_df['tss_start'].values
    distances = row['position'] - starts
    # return smallest absolute value distance
    idx = np.argmin(abs(distances))
    # make relative to gene orientation
    if cluster_gene_df.iloc[idx]['strand'] == '-':
        return -distances[idx]
    else:
        return distances[idx]