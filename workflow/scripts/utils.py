# Data loading functions for PCQTL analysis
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt 
import seaborn as sns 
import networkx as nx
import os
from typing import Dict, List, Optional, Union, Callable
from pandas.errors import EmptyDataError 



def load_vep(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load VEP annotation results for a tissue."""
    sample_vep = pd.read_table('{}/{}/{}.v8.leadvars.vep.vcf'.format(config["working_dir"], config['annotations_output_dir'], tissue_id), skiprows=4)
    overlap_df = load_overlap(config, tissue_id)
    return pd.merge(sample_vep, overlap_df, left_on='ID', right_on='lead_variant_id', how='outer')

def load_susie_annotated(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load SuSiE annotated results for a tissue."""
    qtl_path = '{}/{}{}/{}.v8.susie_R_vars.annotated.csv'.format(config["working_dir"], config['annotations_output_dir'], tissue_id, tissue_id)
    susie_annotated = pd.read_table(qtl_path)
    susie_annotated['cs_num'] = susie_annotated['cs_id'] 
    susie_annotated['cs_id'] = susie_annotated['phenotype_id'] + '_cs_' + susie_annotated['cs_id'].astype(str)
    susie_annotated = add_lead_var(susie_annotated)
    return susie_annotated

def load_pc_annotated(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load PC annotation results for a tissue."""
    return pd.read_table("{}/{}{}/{}.v8.pcs_annotated.txt".format(config["working_dir"], config["annotations_output_dir"], tissue_id, tissue_id))

def load_clusters_annotated(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load annotated clusters for a tissue."""
    annot_cluster = pd.read_csv('{}/{}/{}/{}.clusters.annotated.txt'.format(config["working_dir"], config['annotations_output_dir'], tissue_id, tissue_id), index_col=0)
    annot_cluster['abs_cor'] = np.where(annot_cluster['mean_neg_corr'].isna(), annot_cluster['mean_pos_corr'], - annot_cluster['mean_neg_corr'])
    annot_cluster['abs_cor'] = np.where(annot_cluster['abs_cor'] < annot_cluster['mean_pos_corr'], annot_cluster['mean_pos_corr'], annot_cluster['abs_cor'])
    annot_cluster['max_pos_neg_cor'] = np.where(annot_cluster['mean_neg_corr'].isna(), annot_cluster['mean_pos_corr'], annot_cluster['mean_neg_corr'])
    annot_cluster['max_pos_neg_cor'] = np.where(abs(annot_cluster['max_pos_neg_cor']) < annot_cluster['mean_pos_corr'], annot_cluster['mean_pos_corr'], annot_cluster['max_pos_neg_cor'])
    return annot_cluster

def load_null_clusters_annotated(config: Dict[str, str], tissue_id: str, num_genes: int = 2) -> pd.DataFrame:
    """Load annotated null clusters for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.null_{}genes.annotated.txt'.format(config["working_dir"], config['annotations_output_dir'], tissue_id, tissue_id, num_genes), index_col=0)

def load_cluster(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load cluster data for a tissue."""
    cluster_df =  pd.read_csv('{}/{}/{}.clusters.txt'.format(config["working_dir"], config['clusters_dir'], tissue_id),index_col=0)
    return cluster_df

def load_pc_cis(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load pcQTL cis results for a tissue."""
    pc_cis_path = '{}/{}/{}/{}.v8.pcs.cis_qtl.txt.gz'.format(config["working_dir"], config['pcqtl_output_dir'], tissue_id, tissue_id)
    pc_cis_df = pd.read_table(pc_cis_path)
    pc_cis_df['cluster_id'] = pc_cis_df['phenotype_id'].str.split('_pc').str[0]
    annotate_pc_order(pc_cis_df)
    return pc_cis_df

def load_e_cis(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load eQTL cis results for a tissue."""
    e_cis_path = '{}/{}/{}/{}.v8.cluster_genes.cis_qtl.txt.gz'.format(config["working_dir"], config['eqtl_output_dir'], tissue_id, tissue_id)
    e_cis_df = pd.read_table(e_cis_path)
    e_cis_df['cluster_id'] = e_cis_df['phenotype_id'].str.split('_e').str[0]
    return e_cis_df

def load_e_nominal(config: Dict[str, str], tissue_id: str, chr_id: int = 22, get_var_position: bool = False) -> pd.DataFrame:
    """Load eQTL nominal results for a tissue and chromosome."""
    eqtl_output_dir = config['eqtl_output_dir']
    path = f'{config["working_dir"]}/{eqtl_output_dir}/{tissue_id}/{tissue_id}.v8.cluster_genes.cis_qtl_pairs.chr{chr_id}.parquet'
    e_nominal_df = pd.read_parquet(path)
    if get_var_position:
        e_nominal_df['variant_pos'] = var_pos(e_nominal_df)
    e_nominal_df['cluster_id'] = e_nominal_df['phenotype_id'].str.split('_e_').str[0]
    e_nominal_df['var_cluster'] = e_nominal_df['variant_id'] + '_' + e_nominal_df['cluster_id']
    e_nominal_df['egene_id'] = e_nominal_df['phenotype_id'].str.split('_e_').str[1]
    return e_nominal_df

def load_pc_nominal(config: Dict[str, str], tissue_id: str, chr_id: int = 22, get_var_position: bool = False) -> pd.DataFrame:
    """Load pcQTL nominal results for a tissue and chromosome."""
    pcqtl_output_dir = config['pcqtl_output_dir']
    path = f'{config["working_dir"]}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.cis_qtl_pairs.chr{chr_id}.parquet'
    pc_nominal_df = pd.read_parquet(path)
    if get_var_position:
        pc_nominal_df['variant_pos'] = var_pos(pc_nominal_df)
    pc_nominal_df['cluster_id'] = pc_nominal_df['phenotype_id'].str.split('_pc').str[0]
    pc_nominal_df['var_cluster'] = pc_nominal_df['variant_id'] + '_' + pc_nominal_df['cluster_id']
    pc_nominal_df['pc_num'] = pc_nominal_df['phenotype_id'].str.split('_pc').str[-1].astype(int)
    return pc_nominal_df

def load_pc_susie_r(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load pcQTL SuSiE R results for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.pcs.susie_R.txt'.format(config["working_dir"], config['pcqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_pc_susie(config: Dict[str, str], tissue_id: str, use_r: bool = False) -> pd.DataFrame:
    """Load pcQTL SuSiE results for a tissue."""
    if use_r:
        return load_pc_susie_r(config, tissue_id)
    else:
        return pd.read_csv('{}/{}/{}/{}.v8.pcs.susie.txt'.format(config["working_dir"], config['pcqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_e_susie_r(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load eQTL SuSiE R results for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.cluster_genes.susie_R.txt'.format(config["working_dir"], config['eqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_e_susie(config: Dict[str, str], tissue_id: str, use_r: bool = False) -> pd.DataFrame:
    """Load eQTL SuSiE results for a tissue."""
    if use_r:
        return load_e_susie_r(config, tissue_id)
    else:
        return pd.read_csv('{}/{}/{}/{}.v8.cluster_genes.susie.txt'.format(config["working_dir"], config['eqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_pairwise_coloc(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load pairwise co-localization results for a tissue."""
    return pd.read_csv('{}/{}/pairs/{}.v8.pairs_coloc.txt'.format(config["working_dir"], config['coloc_output_dir'], tissue_id), sep='\t')

def load_tissue_ids(config: Dict[str, str]) -> pd.DataFrame:
    """Load tissue IDs from configuration."""
    return pd.read_csv(f'{config["working_dir"]}/{config["tissue_id_path"]}')['Tissue'].tolist()

def load_tissue_df(config: Dict[str, str]) -> pd.DataFrame:
    """Load tissue information DataFrame."""
    return pd.read_csv(f'{config["working_dir"]}/{config["tissue_id_path"]}')

def load_pc_top_per_phenotype(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load top pcQTL per phenotype for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.pcs.cis_qtl.txt.gz'.format(config["working_dir"], config['pcqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_e_top_per_phenotype(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load top eQTL per phenotype for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.cluster_genes.cis_qtl.txt.gz'.format(config["working_dir"], config['eqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_pc_permutation(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load pcQTL permutation results for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.pcs.cis_independent_qtl.txt.gz'.format(config["working_dir"], config['pcqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_e_permutation(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load eQTL permutation results for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.cluster_genes.cis_independent_qtl.txt.gz'.format(config["working_dir"], config['eqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_e_nominal_all_chr(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load eQTL nominal results for all chromosomes for a tissue."""
    all_chr_data = []
    for chr_id in range(1, 23):
        print(f"Loading chromosome {chr_id} for tissue {tissue_id}...")
        try:
            chr_data = load_e_nominal(config, tissue_id, chr_id=chr_id)
            all_chr_data.append(chr_data)
            print(f"Successfully loaded chromosome {chr_id}")
        except Exception as e:
            print(f"Error loading chromosome {chr_id} for tissue {tissue_id}: {e}")
            continue
    
    if all_chr_data:
        return pd.concat(all_chr_data, ignore_index=True)
    else:
        return pd.DataFrame()

def load_pc_nominal_all_chr(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load pcQTL nominal results for all chromosomes for a tissue."""
    all_chr_data = []
    for chr_id in range(1, 23):
        print(f"Loading chromosome {chr_id} for tissue {tissue_id}...")
        try:
            chr_data = load_pc_nominal(config, tissue_id, chr_id=chr_id)
            all_chr_data.append(chr_data)
            print(f"Successfully loaded chromosome {chr_id}")
        except Exception as e:
            print(f"Error loading chromosome {chr_id} for tissue {tissue_id}: {e}")
            continue
    
    if all_chr_data:
        return pd.concat(all_chr_data, ignore_index=True)
    else:
        return pd.DataFrame()

def load_expression(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load expression data for a tissue."""
    return pd.read_csv('{}/{}/{}.v8.normalized_expression.bed'.format(config["working_dir"], config['expression_dir'], tissue_id), sep='\t')

def load_cluster_expression(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load cluster expression data for a tissue."""
    return pd.read_csv('{}/{}/{}.v8.normalized_residualized_expression.cluster_genes.bed'.format(config["working_dir"], config['filtered_expression_output_dir'], tissue_id), sep='\t')

def load_pc(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load principal components for a tissue."""
    return pd.read_csv('{}/{}/{}.pcs.bed'.format(config["working_dir"], config['pc_output_dir'], tissue_id), sep='\t')

def load_qtl_signal_groups(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load QTL signal groups for a tissue."""
    qtl_signal_groups_path = '{}/{}/qtl_signal_groups/{}.qtl_signal_groups.txt'.format(config["working_dir"], config['coloc_output_dir'], tissue_id)
    return pd.read_csv(qtl_signal_groups_path, sep='\t')

def load_gwas_signal_groups(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load GWAS signal groups for a tissue."""
    gwas_signal_groups_path = '{}/{}/gwas_signal_groups/{}.gwas_signal_groups.txt'.format(config["working_dir"], config['coloc_output_dir'], tissue_id)
    return pd.read_csv(gwas_signal_groups_path, sep='\t')

# Utility functions for data processing
def var_pos(df: pd.DataFrame, column: str = 'variant_id') -> pd.Series:
    """Extract variant positions from variant IDs."""
    return df[column].str.split('_').str[1].astype(int)

def add_lead_var(susie_df: pd.DataFrame) -> pd.DataFrame:
    """Add lead variant information to SuSiE results."""
    susie_df['lead_variant_id'] = susie_df['variant_id']
    return susie_df

def add_num_vars_cs(susie_df: pd.DataFrame) -> pd.DataFrame:
    """Add number of variants per credible set to SuSiE results."""
    susie_df['num_vars_cs'] = susie_df.groupby('cs_id')['variant_id'].transform('count')
    return susie_df

def load_across_tissues(config: Dict[str, str], load_func: Callable, tissue_ids: Optional[List[str]] = None) -> pd.DataFrame:
    """Load data across multiple tissues using a specified loading function."""
    if tissue_ids is None:
        tissue_ids = load_tissue_ids(config)
    
    results = []
    for tissue_id in tissue_ids:
        try:
            result = load_func(config, tissue_id)
            result['tissue_id'] = tissue_id
            results.append(result)
        except Exception as e:
            print(f"Error loading data for tissue {tissue_id}: {e}")
            continue
    
    if results:
        return pd.concat(results, ignore_index=True)
    else:
        return pd.DataFrame()

def annotate_pc_order(pc_df: pd.DataFrame) -> pd.DataFrame:
    """Annotate PC order in pcQTL results."""
    pc_df['pc_num'] = pc_df['phenotype_id'].str.split('_pc').str[-1].astype(int)
    pc_df['pc_order'] = 'middle'
    pc_df.loc[pc_df['pc_num'] == 1, 'pc_order'] = 'first'
    pc_df.loc[pc_df['pc_num'] == pc_df['num_genes'], 'pc_order'] = 'last'
    return pc_df

def remove_cross_map(df: pd.DataFrame, cross_map_ids: Optional[List[str]] = None, config: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """Remove cross-mappable genes from results."""
    if cross_map_ids is None and config is not None:
        # Load cross-mappable genes if not provided
        cross_map_df = pd.read_csv(config['cross_map_path'])
        cross_map_ids = cross_map_df['gene_id'].tolist()
    
    if cross_map_ids is not None:
        df = df[~df['gene_id'].isin(cross_map_ids)]
    
    return df