#!/usr/bin/env python3
"""
Gene Expression Cluster Annotation Script

This script annotates gene expression clusters with various genomic and functional features
including ABC enhancer-gene connections, CTCF binding sites, TAD boundaries, paralogous genes,
GO terms, cross-mappability, and correlation statistics. 
"""

import logging
import sys
import os
import numpy as np
import pandas as pd
import argparse
from tqdm import tqdm 
import yaml
import ast
from pathlib import Path

# Import local modules - use relative imports if possible
try:
    from residualize import calculate_residual
except ImportError:
    # Fallback for when script is run directly
    import sys
    sys.path.append(str(Path(__file__).parent))
    from residualize import calculate_residual

# Configuration - use environment variables with defaults
PREFIX = os.getenv('PCQTL_PREFIX', '/home/klawren/oak/pcqtls')

# Constants for genomic analysis
CHROMOSOMES = list(range(1, 23))  # Human autosomes
BIDIRECTIONAL_PROMOTER_DISTANCE = 1000  # Base pairs
HIGH_POSITIVE_CORRELATION_THRESHOLD = 0.5


def setup_logging(level=logging.INFO):
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def validate_input_files(cluster_path, expression_path, covariates_path):
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





def load_tad(tad_path=f'{PREFIX}/data/references/TAD_annotations/TADs_hg38/converted_HiC_IMR90_DI_10kb.txt'):
    """
    Load TAD (Topologically Associating Domain) boundary data.
    
    Args:
        tad_path (str): Path to TAD boundary file
    
    Returns:
        pd.DataFrame: TAD boundaries with chr, start, end columns
    """
    logger = logging.getLogger(__name__)
    logger.debug(f'Loading TAD boundaries from {tad_path}')
    tad_df = pd.read_table(tad_path, header=None, names=['chr', 'start', 'end'])
    tad_df['chr'] = tad_df['chr'].str.strip('chr')
    tad_df = tad_df[tad_df['chr'].isin([str(i) for i in CHROMOSOMES])]
    tad_df['chr'] = tad_df['chr'].astype(int)
    logger.debug(f'Loaded {len(tad_df)} TAD boundaries')
    return tad_df


def load_gencode(gencode_path=f'{PREFIX}/data/references/processed_gencode.v26.GRCh38.genes.txt', 
                protein_coding_only=True):
    """
    Load and process GENCODE gene annotations.
    
    Args:
        gencode_path (str): Path to GENCODE annotation file (TSV with columns: chr,start,end,strand,gene_id,tss_start)
        protein_coding_only (bool): Whether to filter to protein-coding genes only (not used with simplified format)
    
    Returns:
        tuple: (gid_gencode, full_gencode) where gid_gencode is indexed by gene_id
    """
    logger = logging.getLogger(__name__)
    logger.debug(f'Loading GENCODE annotations from {gencode_path}')
    
    # Load gene data (tab-separated)
    full_gencode = pd.read_table(gencode_path)
    logger.debug(f'Loaded {len(full_gencode)} total genes from GENCODE')
    
    # Add gene_name column if not present (you'll add this back)
    if 'gene_name' not in full_gencode.columns:
        full_gencode['gene_name'] = full_gencode['gene_id']
    
    gid_gencode = full_gencode.set_index('gene_id').drop_duplicates()
    logger.debug(f'Created gene-indexed GENCODE with {len(gid_gencode)} unique genes')
    return gid_gencode, full_gencode


def get_residual_expression(covariates_path, expression_path):
    """
    Calculate residualized expression data by removing covariate effects.
    
    Args:
        covariates_path (str): Path to covariates file
        expression_path (str): Path to expression data file
    
    Returns:
        pd.DataFrame: Residualized expression data
    """
    logger = logging.getLogger(__name__)
    logger.debug(f'Calculating residualized expression from {expression_path}')
    
    # Load covariates and expression data
    covariates_df = pd.read_table(covariates_path, index_col=0).T
    expression_df = pd.read_table(expression_path)
    
    logger.debug(f'Loaded {len(covariates_df)} covariates and {len(expression_df)} genes')
    
    # Set up expression data for residualization
    gid_expression = expression_df.set_index('gene_id')[covariates_df.index]
    
    # Calculate residuals
    residual_exp = calculate_residual(gid_expression[covariates_df.index], covariates_df, center=True)
    residual_exp = pd.DataFrame(residual_exp, columns=covariates_df.index, index=gid_expression.index)
    
    logger.debug(f'Calculated residuals for {len(residual_exp)} genes')
    return residual_exp


def load_abc(my_tissue_id, full_gencode=None, 
            full_abc_path=f'{PREFIX}/data/references/functional_annotations/ABC_predictions/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz', 
            abc_match_path=f'{PREFIX}/data/references/functional_annotations/ABC_predictions/ABC_matched_gtex.txt'):
    """
    Load ABC (Activity-By-Contact) enhancer-gene connection data for a specific tissue.
    
    Args:
        my_tissue_id (str): GTEx tissue identifier
        full_gencode (pd.DataFrame): GENCODE annotations (loaded if None)
        full_abc_path (str): Path to ABC predictions file
        abc_match_path (str): Path to ABC-GTEx tissue matching file
    
    Returns:
        pd.DataFrame: ABC enhancer-gene connections for the specified tissue
    """
    # Load GENCODE if not provided
    if full_gencode is None:
        gid_gencode, full_gencode = load_gencode()
    
    # Load ABC data for enhancer-gene connections
    full_abc_pred_df = pd.read_table(full_abc_path)
    
    # Load tissue matching for ABC-GTEx tissues
    abc_gtex_match = pd.read_csv(abc_match_path)
    
    # Get enhancer-gene connections for the matched tissue
    abc_df = full_abc_pred_df[full_abc_pred_df['CellType'] == 
                             abc_gtex_match[abc_gtex_match['GTEX_tissue'] == my_tissue_id]['ABC_biosample_id'].iloc[0]]
    
    # Add transcript IDs to ABC enhancer-gene connection columns
    gene_enhancer_df = pd.merge(full_gencode[['gene_id', 'gene_name']], 
                               abc_df[['chr', 'start', 'end', 'TargetGene', 'name', 'class', 'ABC.Score']], 
                               left_on='gene_name', right_on='TargetGene', how='left')
    gene_enhancer_df.rename(columns={'name': 'enhancer'}, inplace=True)
    gene_enhancer_df.set_index('gene_id', inplace=True)
    gene_enhancer_df.dropna(inplace=True)

    # Process chromosome and enhancer coordinates
    gene_enhancer_df['chr'] = gene_enhancer_df['chr'].str.split('chr').str[1]
    gene_enhancer_df = gene_enhancer_df[gene_enhancer_df['chr'].isin([str(i) for i in CHROMOSOMES])]
    gene_enhancer_df['chr'] = gene_enhancer_df['chr'].astype(int)
    gene_enhancer_df['enhancer_start'] = gene_enhancer_df['enhancer'].str.split(':').str[1].str.split('-').str[0].astype(int)
    gene_enhancer_df['enhancer_end'] = gene_enhancer_df['enhancer'].str.split(':').str[1].str.split('-').str[1].astype(int)
    
    return gene_enhancer_df


def load_ctcf(my_tissue_id, 
             ctcf_match_path=f'{PREFIX}/data/references/functional_annotations/ctcf_chip/ctcf_matched_gtex.txt', 
             ctcf_dir=f'{PREFIX}/data/references/functional_annotations/ctcf_chip'):
    """
    Load CTCF binding site data for a specific tissue.
    
    Args:
        my_tissue_id (str): GTEx tissue identifier
        ctcf_match_path (str): Path to CTCF-GTEx tissue matching file
        ctcf_dir (str): Directory containing CTCF binding site files
    
    Returns:
        pd.DataFrame: CTCF binding sites for the specified tissue
    """
    # Load CTCF data
    ctcf_gtex_match = pd.read_csv(ctcf_match_path)
    ctcf_file = ctcf_gtex_match[ctcf_gtex_match['GTEX'] == my_tissue_id].iloc[0]['ctcf']
    
    ctcf_df = pd.read_table(f'{ctcf_dir}/{ctcf_file}.bed.gz', 
                           names=['chr', 'start', 'end', 'name', 'score', 'strand', 
                                 'signal_value', 'p_value', 'q_value', 'peak'])
    
    # Create interval arrays for efficient overlap queries
    ctcf_df['peak_inter'] = pd.arrays.IntervalArray.from_arrays(ctcf_df['start'], ctcf_df['end'])
    ctcf_df['point_inter'] = pd.arrays.IntervalArray.from_arrays(ctcf_df['start'] + ctcf_df['peak'], 
                                                               ctcf_df['start'] + ctcf_df['peak'] + 1)
    ctcf_df['chr'] = ctcf_df['chr'].str.split('chr').str[1].astype(str)
    
    return ctcf_df


def load_paralogs(paralog_path=f'{PREFIX}/data/references/functional_annotations/paralogs_biomart_ensembl97.tsv.gz'):
    """
    Load paralogous gene data.
    
    Args:
        paralog_path (str): Path to paralog data file
    
    Returns:
        pd.Series: Grouped paralog data with gene IDs as index and sets of paralog IDs as values
    """
    # Load paralogs
    paralog_df = pd.read_table(paralog_path)
    
    # Drop genes that don't have paralogues
    paralog_df.dropna(subset=['Human paralogue gene stable ID'], inplace=True)
    
    # Group by the gene that has the paralogs (this is bidirectional)
    # Note: This uses gene IDs without versions as versions may differ between datasets
    paralog_df = paralog_df.groupby('Gene stable ID')['Human paralogue gene stable ID'].apply(set)
    
    return paralog_df


def load_go(go_path=f'{PREFIX}/data/references/functional_annotations/go_biomart_ensembl97.tsv.gz'):
    """
    Load GO (Gene Ontology) term data.
    
    Args:
        go_path (str): Path to GO term file
    
    Returns:
        pd.Series: Grouped GO terms with gene IDs as index and sets of GO accessions as values
    """
    # Load GO terms
    go_df = pd.read_table(go_path, header=None,
                        names = ['Gene stable ID', 'Gene stable ID version', 'Gene start (bp)', 'Gene end', 'Strand', 
                                 'tss', 'gencode_annotation', 'gene_name', 'transcript_type', 'go_accession', 'go_name', 'go_domain'])
    
    # Only consider matching biological process GO terms, as in Ribiero 2021
    go_df = go_df[go_df['go_domain'] == 'biological_process'].groupby('Gene stable ID')['go_accession'].apply(set)
    
    return go_df


def load_cross_map(cross_map_path=f'{PREFIX}/data/references/cross_mappability/cross_mappability_100_agg.csv'):
    """
    Load cross-mappability data.
    
    Args:
        cross_map_path (str): Path to cross-mappability file
    
    Returns:
        pd.DataFrame: Cross-mappability data with gene_1 as index and a dictionary of gene_2 and cross_mappability as columns
    """
    # Load cross mapablity (this is cleaned up in cross_mappability.ipynb)
    cross_mappability = pd.read_csv(cross_map_path)
    cross_mappability.set_index('gene_1', inplace=True)
    return cross_mappability



# Functions to annotate cluster properties

def get_cluster_size(row, gid_gencode):
    """
    Calculate the genomic span of a cluster.
    
    Args:
        row (pd.Series): Cluster row with cluster_id
        gid_gencode (pd.DataFrame): GENCODE annotations indexed by gene_id
    
    Returns:
        int: Genomic span (end - start) of the cluster
    """
    gene_ids = row['cluster_id'].split('_')
    cluster_gencode = gid_gencode.loc[gene_ids]
    return cluster_gencode['end'].max() - cluster_gencode['start'].min()


def get_cluster_tss_size(row, gid_gencode):
    """
    Calculate the TSS span of a cluster.
    
    Args:
        row (pd.Series): Cluster row with cluster_id
        gid_gencode (pd.DataFrame): GENCODE annotations indexed by gene_id
    
    Returns:
        int: TSS span (max TSS - min TSS) of the cluster
    """
    gene_ids = row['cluster_id'].split('_')
    cluster_gencode = gid_gencode.loc[gene_ids]
    return cluster_gencode['tss_start'].max() - cluster_gencode['tss_start'].min()


def get_num_overlapping(row, gid_gencode):
    """
    Count the number of overlapping TSS regions within a cluster.
    
    Args:
        row (pd.Series): Cluster row with cluster_id
        gid_gencode (pd.DataFrame): GENCODE annotations indexed by gene_id
    
    Returns:
        int: Number of overlapping TSS regions (within 1kb)
    """
    gene_ids = row['cluster_id'].split('_')
    cluster_gencode = gid_gencode.loc[gene_ids]
    return sum((cluster_gencode['tss_start'] - cluster_gencode['tss_start'].shift(1)) < 1000)

# def get_cluster_start_ids(cluster_df):
#     # the first for pairs, the first and second for threes, ect
#     cluster_start_ids = []
#     for i in range(cluster_df['N_genes'].max()):
#         out_ids = cluster_df[cluster_df['N_genes'] == i]['Transcripts'].str.split(',').str[:i-1].values
#         if len(out_ids)>0:
#             cluster_start_ids.append(np.concatenate(out_ids))
#         else:
#             cluster_start_ids.append([])
#     return cluster_start_ids

def annotate_sizes(cluster_df, gid_gencode):
    """
    Add cluster size annotations to the cluster DataFrame.
    
    Args:
        cluster_df (pd.DataFrame): Cluster DataFrame to annotate
        gid_gencode (pd.DataFrame): GENCODE annotations indexed by transcript_id
    """
    cluster_df['cluster_size'] = cluster_df.apply(get_cluster_size, axis=1, args=(gid_gencode,))
    cluster_df['cluster_tss_size'] = cluster_df.apply(get_cluster_tss_size, axis=1, args=(gid_gencode,))


def annotate_positions(cluster_df, gid_gencode):
    """
    Add genomic position annotations to the cluster DataFrame.
    
    Args:
        cluster_df (pd.DataFrame): Cluster DataFrame to annotate
        gid_gencode (pd.DataFrame): GENCODE annotations indexed by gene_id
    """
    for idx, row in cluster_df.iterrows():
        gene_ids = row['cluster_id'].split('_')
        cluster_gencode = gid_gencode.loc[gene_ids]
        cluster_df.loc[idx, 'start'] = cluster_gencode['start'].min()
        cluster_df.loc[idx, 'end'] = cluster_gencode['end'].max()
        cluster_df.loc[idx, 'tss_min'] = cluster_gencode['tss_start'].min()
        cluster_df.loc[idx, 'tss_max'] = cluster_gencode['tss_start'].max()


def get_bidirectional(row, gid_gencode):
    """
    Count bidirectional promoter pairs in a cluster.
    
    Args:
        row (pd.Series): Cluster row with cluster_id
        gid_gencode (pd.DataFrame): GENCODE annotations indexed by gene_id
    
    Returns:
        int: Number of bidirectional promoter pairs (divided by 2 to avoid double counting)
    """
    gene_ids = row['cluster_id'].split('_')
    cluster_gencode = gid_gencode.loc[gene_ids]
    num_bidirectional = 0
    
    # Check all pairwise combinations of genes
    for idx, first_gene_row in cluster_gencode.iterrows():
        for idx, second_gene_row in cluster_gencode.iterrows():
            opp_strand = first_gene_row['strand'] != second_gene_row['strand']
            close = abs(first_gene_row['tss_start'] - second_gene_row['tss_start']) <= BIDIRECTIONAL_PROMOTER_DISTANCE
            if opp_strand & close:
                # Found a bidirectional promoter
                num_bidirectional += 1

    # Return half to avoid double counting
    return num_bidirectional / 2


def annotate_bidirectional(cluster_df, gid_gencode):
    """
    Add bidirectional promoter annotations to the cluster DataFrame.
    
    Args:
        cluster_df (pd.DataFrame): Cluster DataFrame to annotate
        gid_gencode (pd.DataFrame): GENCODE annotations indexed by transcript_id
    """
    cluster_df['num_bidirectional_promoter'] = cluster_df.apply(get_bidirectional, axis=1, args=(gid_gencode,))
    cluster_df['has_bidirectional_promoter'] = cluster_df['num_bidirectional_promoter'] > 0



def annotate_ctcf(cluster_df, ctcf_df):
    cluster_df['interval'] = pd.arrays.IntervalArray.from_arrays(cluster_df['tss_min'], cluster_df['tss_max'])
    # ctcf intervals for each chromosome
    chr_ctcf_peaks={}
    chr_ctcf_points={}
    for chr_val in cluster_df['chr'].unique():
        # Handle both chr prefix and numeric chromosome formats
        chr_str = chr_val if isinstance(chr_val, str) else f'chr{chr_val}'
        ctcf_chr = ctcf_df[ctcf_df['chr'] == chr_str]
        chr_ctcf_peaks[chr_val] = pd.arrays.IntervalArray.from_arrays(ctcf_chr['start'], ctcf_chr['end'])
        chr_ctcf_points[chr_val] = pd.arrays.IntervalArray.from_arrays(ctcf_chr['start'] + ctcf_chr['peak'], ctcf_chr['start'] + ctcf_chr['peak'] + 1)
    # annotate each cluster
    for idx, row in tqdm(cluster_df.iterrows(), total=len(cluster_df)):
        try:
            num_ctcf_peak = sum(chr_ctcf_peaks[row['chr']].overlaps(row['interval']))
            num_ctcf_point = sum(chr_ctcf_points[row['chr']].overlaps(row['interval']))
        except KeyError:
            num_ctcf_peak = 0
            num_ctcf_point = 0
        cluster_df.loc[idx, 'num_ctcf_peak'] = num_ctcf_peak
        cluster_df.loc[idx, 'has_ctcf_peak'] = num_ctcf_peak > 0
        cluster_df.loc[idx, 'num_ctcf_point'] = num_ctcf_point
        cluster_df.loc[idx, 'has_ctcf_point'] = num_ctcf_point > 0


def count_tad_overlap(row, tad_df, inter_column):
    # Handle both chr prefix and numeric chromosome formats
    chr_val = row.chr if isinstance(row.chr, str) else f'chr{row.chr}'
    tad_chr = tad_df[tad_df['chr']==chr_val]
    # intervals for the TADs and TAD edges
    chr_tad_intervals = pd.arrays.IntervalArray.from_arrays(tad_chr['start'], tad_chr['end'])
    chr_tad_edges = pd.arrays.IntervalArray.from_arrays(pd.concat([tad_chr['start'], tad_chr['end']]), pd.concat([tad_chr['start']+1, tad_chr['end']+1]))

    # number of TADs and TAD edges overlapped
    try:
        tad_overlaps = sum(chr_tad_intervals.overlaps(row[inter_column]))
        tad_edge_overlaps = sum(chr_tad_edges.overlaps(row[inter_column]))
    except TypeError as e:
        tad_overlaps = 0
        tad_edge_overlaps = 0
    # 0: in no TADs
    if tad_edge_overlaps==0 and tad_overlaps==0:
        return 0
    # 2: crossing a TAD edge
    elif tad_edge_overlaps > 0:
        return 2
    # 1: containted within 1 TAD (can be nested inside multiple TADs)
    else:
        return 1

def annotate_tads(cluster_df, tad_df):
    """
    Annotate clusters with TAD boundary information.
    
    Args:
        cluster_df (pd.DataFrame): Cluster DataFrame to annotate
        tad_df (pd.DataFrame): TAD boundary data
    """
    # Create interval arrays for efficient overlap queries
    cluster_df['tss_inter'] = pd.arrays.IntervalArray.from_arrays(cluster_df['tss_min'], cluster_df['tss_max'])
    cluster_df['gene_inter'] = pd.arrays.IntervalArray.from_arrays(cluster_df['start'], cluster_df['end'])
    
    # Count TAD overlaps for gene and TSS regions
    cluster_df['num_tads_gene'] = cluster_df.apply(count_tad_overlap, axis=1, args=(tad_df, 'gene_inter'))
    cluster_df['num_tads_tss'] = cluster_df.apply(count_tad_overlap, axis=1, args=(tad_df, 'tss_inter'))
    
    # Create boolean flags for TAD spanning
    cluster_df['has_tads_gene'] = cluster_df['num_tads_gene'] > 1
    cluster_df['has_tads_tss'] = cluster_df['num_tads_tss'] > 1






def annotate_correlation(cluster_df, residual_exp):
    """
    Annotate clusters with correlation statistics.
    
    Args:
        cluster_df (pd.DataFrame): Cluster DataFrame to annotate
        residual_exp (pd.DataFrame): Residualized expression data
    """
    for idx, row in tqdm(cluster_df.iterrows(), total=len(cluster_df)):
        transcript_list = row['cluster_id'].split('_')
        
        # Calculate Spearman correlations between genes in cluster
        cluster_expression = residual_exp.loc[transcript_list].T.corr('spearman').to_numpy()
        cluster_corr = cluster_expression[np.triu_indices(len(cluster_expression), k=1)]
        
        # Store correlation statistics
        cluster_df.loc[idx, 'mean_corr'] = cluster_corr.mean()
        cluster_df.loc[idx, 'corr_list'] = str(cluster_corr)
        cluster_df.loc[idx, 'mean_pos_corr'] = cluster_corr[cluster_corr > 0].mean()
        cluster_df.loc[idx, 'mean_neg_corr'] = cluster_corr[cluster_corr < 0].mean()


def annotate_go(cluster_df, go_df):
    for idx, row in tqdm(cluster_df.iterrows(), total=len(cluster_df)):
        transcript_list_versions = row['cluster_id'].split('_')
        transcript_list_no_versions = [transcript.split('.')[0] for transcript in transcript_list_versions]

        go_list = go_df[go_df.index.isin(transcript_list_no_versions)]
        num_shared_go_all = sum(go_list.duplicated())
        # number genes that share all their go terms with another gene
        cluster_df.loc[idx, 'num_shared_go_all'] = num_shared_go_all
        cluster_df.loc[idx, 'has_shared_go_all'] = num_shared_go_all > 0

        # number go terms shared between any genes
        num_shared_go_any = sum(go_list.explode().duplicated())
        cluster_df.loc[idx, 'num_shared_go_any'] = num_shared_go_any
        cluster_df.loc[idx, 'has_shared_go_any'] = num_shared_go_any > 0


# annotate paralogs
def get_paralogs(row, paralog_df):
    transcript_list_versions = row['cluster_id'].split('_')
    transcript_list_no_versions = set([transcript.split('.')[0] for transcript in transcript_list_versions])

    paralogs = 0
    for transcript in transcript_list_no_versions:
        try:
            has_paralog = bool(paralog_df.loc[transcript] & transcript_list_no_versions)
            paralogs += has_paralog
        except KeyError:
            # if this isn't in the paralog df, it has no paralogs, so continue
            pass
    return paralogs

def annotate_paralogs(cluster_df, paralog_df):
    cluster_df['num_paralog'] = cluster_df.apply(get_paralogs, axis=1, args=(paralog_df,))
    cluster_df['has_paralog'] = cluster_df['num_paralog'] > 0


def get_cross_map(row, cross_mappability, cross_map_threshold=100):
    # number of transcripts that cross map to some other transcript in the cluster
    transcript_list = set(row['cluster_id'].split('_'))
    cross_maps = 0
    for transcript in transcript_list:
        try:
            cross_map_this_transcript = cross_mappability.loc[transcript]
            pass_threshold_mask = np.asarray(ast.literal_eval(cross_map_this_transcript['cross_mappability'])) > cross_map_threshold
            cross_map_this_transcript = np.asarray(ast.literal_eval(cross_map_this_transcript['gene_2_full']))[pass_threshold_mask]
            has_cross_map = bool(set(cross_map_this_transcript) & transcript_list)
            cross_maps += has_cross_map
        except KeyError:
            # if this isn't in the paralog df, it has no paralogs, so continue
            pass
    return cross_maps

def annotate_cross_maps(cluster_df, cross_mappability):
    cluster_df['num_cross_map'] = cluster_df.apply(get_cross_map, axis=1, args=(cross_mappability,))
    cluster_df['has_cross_map'] = cluster_df['num_cross_map'] > 0

def annotate_complexes(cluster_df, complex_df):
    for idx, row in tqdm(cluster_df.iterrows(), total=len(cluster_df)):
        complex_list = complex_df[complex_df.index.isin(row['cluster_id'].split('_'))]
        num_complexes = sum(complex_list.explode('ComplexID').duplicated())
        cluster_df.loc[idx, 'num_complexes'] = num_complexes
        cluster_df.loc[idx, 'has_complexes'] = num_complexes > 0




# function to add all annotations, give correctly loaded data
def add_annotations(cluster_df, gid_gencode, gene_enhancer_df, paralog_df, cross_mappability, go_df, ctcf_df, residual_exp, tad_df):
    """
    Add all annotations to the cluster DataFrame using the general annotation framework.
    
    This version is more efficient and easier to extend than the original add_annotations function.
    
    Args:
        cluster_df (pd.DataFrame): Cluster DataFrame to annotate
        gid_gencode (pd.DataFrame): GENCODE annotations
        gene_enhancer_df (pd.DataFrame): ABC enhancer-gene connections
        paralog_df (pd.Series): Paralog data
        cross_mappability (pd.DataFrame): Cross-mappability data
        go_df (pd.Series): GO term data
        ctcf_df (pd.DataFrame): CTCF binding site data
        residual_exp (pd.DataFrame): Residualized expression data
        tad_df (pd.DataFrame): TAD boundary data
        
    Returns:
        pd.DataFrame: Annotated cluster DataFrame
    """
    logger = logging.getLogger(__name__)
    logger.info(f'Starting general annotation of {len(cluster_df)} clusters')
    
    # Reset index for consistent processing
    cluster_df.reset_index(drop=True, inplace=True)
    
    # Add basic cluster properties first (these don't use the framework)
    logger.info('Adding cluster size annotations...')
    annotate_sizes(cluster_df, gid_gencode)
    
    # Add genomic positions (must be done before interval-based annotations)
    logger.info('Adding genomic position annotations...')
    annotate_positions(cluster_df, gid_gencode)
    
    # Add bidirectional promoter information (special case)
    logger.info('Adding bidirectional promoter annotations...')
    annotate_bidirectional(cluster_df, gid_gencode)
    
    # Create and add interval-based annotators
    logger.info('Setting up interval-based annotators...')
    
    # TAD annotator
    annotate_tads(cluster_df, tad_df)
    
    # CTCF annotator
    annotate_ctcf(cluster_df, ctcf_df)
    
    # Add gene-link-based annotations
    logger.info('Adding gene-link-based annotations...')
    
    # ABC enhancer annotation (detailed)
    annotate_abc_detailed(cluster_df, gene_enhancer_df)
    
    # Paralog annotation
    annotate_paralogs_general(cluster_df, paralog_df)
    
    # GO terms annotation
    annotate_go_general(cluster_df, go_df)
    
    # Cross-mappability annotation
    annotate_cross_maps_general(cluster_df, cross_mappability)
    
    # Correlation annotation
    logger.info('Adding correlation annotations...')
    annotate_correlation(cluster_df, residual_exp)
    cluster_df['has_neg_corr'] = ~cluster_df['mean_neg_corr'].isna()
    cluster_df['has_high_pos_corr'] = cluster_df['mean_pos_corr'] > HIGH_POSITIVE_CORRELATION_THRESHOLD
    
    logger.info(f'Completed general annotation of {len(cluster_df)} clusters')
    return cluster_df








def load_and_annotate(cluster_df, my_tissue_id, covariates_path, expression_path,
                      gencode_path=f'{PREFIX}/data/references/processed_gencode.v26.GRCh38.genes.txt', 
                      full_abc_path = f'{PREFIX}/data/references/functional_annotations/ABC_predictions/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz', 
                      abc_match_path=f'{PREFIX}/data/references/functional_annotations/ABC_predictions/ABC_matched_gtex.txt', 
                      ctcf_match_path=f'{PREFIX}/data/references/functional_annotations/ctcf_chip/ctcf_matched_gtex.txt', 
                      ctcf_dir=f'{PREFIX}/data/references/functional_annotations/ctcf_chip', 
                      paralog_path=f'{PREFIX}/data/references/functional_annotations/paralogs_biomart_ensembl97.tsv.gz', 
                      go_path=f'{PREFIX}/data/references/functional_annotations/go_biomart_ensembl97.tsv.gz', 
                      cross_map_path=f'{PREFIX}/data/references/cross_mappability/cross_mappability_100_agg.csv',
                      tad_path=f'{PREFIX}/data/references/TAD_annotations/TADs_hg38/converted_HiC_IMR90_DI_10kb.txt'):
    """
    Load all annotation data and annotate clusters using the general framework.
    
    Args:
        cluster_df: Cluster DataFrame to annotate
        my_tissue_id: GTEx tissue identifier
        covariates_path: Path to covariates file
        expression_path: Path to expression file
        annotation_config: Configuration for annotators (uses default if None)
        **kwargs: Additional file paths for annotation data
        
    Returns:
        pd.DataFrame: Annotated cluster DataFrame
    """
    logger = logging.getLogger(__name__)
    logger.info(f'Loading annotation data for tissue: {my_tissue_id}')
    
    # Load all annotation data
    logger.info('Loading GENCODE annotations...')
    gid_gencode, full_gencode = load_gencode(gencode_path)
    
    logger.info('Loading ABC enhancer-gene connections...')
    gene_enhancer_df = load_abc(my_tissue_id, full_gencode, full_abc_path, abc_match_path)
    
    logger.info('Loading CTCF binding sites...')
    ctcf_df = load_ctcf(my_tissue_id, ctcf_match_path, ctcf_dir)
    
    logger.info('Loading paralog data...')
    paralog_df = load_paralogs(paralog_path)
    
    logger.info('Loading GO terms...')
    go_df = load_go(go_path)
    
    logger.info('Loading cross-mappability data...')
    cross_mappability = load_cross_map(cross_map_path)
    
    logger.info('Calculating residualized expression...')
    residual_exp = get_residual_expression(covariates_path, expression_path)
    
    logger.info('Loading TAD boundaries...')
    tad_df = load_tad(tad_path)
    
    logger.info('All annotation data loaded successfully')
    
    # Annotate clusters using the general framework
    cluster_df.reset_index(drop=True, inplace=True)
    
    # Add basic cluster properties first
    logger.info('Adding basic cluster properties...')
    annotate_sizes(cluster_df, gid_gencode)
    annotate_positions(cluster_df, gid_gencode)
    annotate_bidirectional(cluster_df, gid_gencode)
    
    # Use general annotation functions
    cluster_df = add_annotations(cluster_df, gid_gencode, gene_enhancer_df, paralog_df, cross_mappability, go_df, ctcf_df, residual_exp, tad_df)
    
    return cluster_df


def parse_abc_data(cluster_df: pd.DataFrame) -> pd.DataFrame:
    """
    Parse stored ABC data back into a more usable format.
    
    This function converts the string-stored ABC data back into proper data structures
    for easier analysis.
    
    Args:
        cluster_df: Cluster DataFrame with ABC annotations
        
    Returns:
        pd.DataFrame: Cluster DataFrame with parsed ABC data
    """
    import ast
    
    logger = logging.getLogger(__name__)
    logger.info("Parsing ABC data for easier analysis...")
    
    # Create new columns for parsed data
    cluster_df['shared_enhancers_parsed'] = None
    cluster_df['shared_enhancer_scores_parsed'] = None
    cluster_df['shared_enhancer_cell_types_parsed'] = None
    
    for idx, row in cluster_df.iterrows():
        try:
            # Parse the string representations back to Python objects
            shared_enhancers = ast.literal_eval(row['shared_enhancers'])
            shared_scores = ast.literal_eval(row['shared_enhancer_scores'])
            shared_cell_types = ast.literal_eval(row['shared_enhancer_cell_types'])
            
            # Store parsed data
            cluster_df.loc[idx, 'shared_enhancers_parsed'] = shared_enhancers
            cluster_df.loc[idx, 'shared_enhancer_scores_parsed'] = shared_scores
            cluster_df.loc[idx, 'shared_enhancer_cell_types_parsed'] = shared_cell_types
            
        except (ValueError, SyntaxError) as e:
            logger.warning(f"Failed to parse ABC data for cluster {idx}: {e}")
            cluster_df.loc[idx, 'shared_enhancers_parsed'] = []
            cluster_df.loc[idx, 'shared_enhancer_scores_parsed'] = []
            cluster_df.loc[idx, 'shared_enhancer_cell_types_parsed'] = []
    
    return cluster_df


def get_abc_summary_stats(cluster_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate summary statistics for ABC data.
    
    Args:
        cluster_df: Cluster DataFrame with parsed ABC data
        
    Returns:
        pd.DataFrame: Summary statistics for ABC data
    """
    logger = logging.getLogger(__name__)
    logger.info("Calculating ABC summary statistics...")
    
    # Initialize summary columns
    cluster_df['abc_mean_score'] = 0.0
    cluster_df['abc_max_score'] = 0.0
    cluster_df['abc_min_score'] = 0.0
    cluster_df['abc_score_std'] = 0.0
    
    for idx, row in cluster_df.iterrows():
        if row['has_shared_enhancer'] and row['shared_enhancer_scores_parsed']:
            # Flatten all scores from all shared enhancers
            all_scores = []
            for score_list in row['shared_enhancer_scores_parsed']:
                all_scores.extend(score_list)
            
            if all_scores:
                cluster_df.loc[idx, 'abc_mean_score'] = np.mean(all_scores)
                cluster_df.loc[idx, 'abc_max_score'] = np.max(all_scores)
                cluster_df.loc[idx, 'abc_min_score'] = np.min(all_scores)
                cluster_df.loc[idx, 'abc_score_std'] = np.std(all_scores)
    
    return cluster_df


def main():
    """
    Main function to handle command-line interface and execute cluster annotation.
    
    Parses command-line arguments and runs the cluster annotation pipeline.
    """
    # Set up logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description='Annotate gene expression clusters with genomic and functional features',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('-t', '--tissue-id', required=True,
                       help='GTEx tissue identifier')
    parser.add_argument('-c', '--clusters', required=True,
                       help='Cluster CSV file')
    parser.add_argument('-e', '--expression', required=True,
                       help='Normalized expression BED file')
    parser.add_argument('-co', '--covariates', required=True,
                       help='Covariates file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output file for annotated clusters')
    
    # Optional annotation files
    parser.add_argument('--gencode', 
                       help='GENCODE annotation file')
    parser.add_argument('--full-abc', 
                       help='ABC predictions file')
    parser.add_argument('--abc-match', 
                       help='ABC-GTEx tissue matching file')
    parser.add_argument('--ctcf-match', 
                       help='CTCF-GTEx tissue matching file')
    parser.add_argument('--ctcf-dir', 
                       help='Directory containing CTCF binding site files')
    parser.add_argument('--paralog', 
                       help='Paralog data file')
    parser.add_argument('--go', 
                       help='GO term file')
    parser.add_argument('--cross-map', 
                       help='Cross-mappability file')
    parser.add_argument('--tad', 
                       help='TAD boundary file')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')


    # Parse arguments
    args = parser.parse_args()
    
    if args.verbose:
        setup_logging(logging.DEBUG)
    
    try:
        # Validate input files
        logger.info("Validating input files...")
        validate_input_files(args.cluster_path, args.expression_path, args.covariates_path)
        
        logger.info(f'Starting cluster annotation for tissue: {args.tissue}')
        
        # Run annotation
        cluster_df_annotated = load_and_annotate(cluster_df, my_tissue_id, covariates_path, expression_path)

        # Write output
        logger.info(f'Writing annotated clusters to {args.out_path}')
        cluster_df_annotated.to_csv(args.out_path, index=False)
        logger.info(f'Successfully wrote {len(cluster_df_annotated)} annotated clusters to {args.out_path}')
        
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

