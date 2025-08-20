# Data loading functions for PCQTL analysis
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt 
import seaborn as sns 
import networkx as nx
import os
from pandas.errors import EmptyDataError 

# Working directory prefix
prefix = '/home/klawren/oak/pcqtls'

def load_vep(config, tissue_id):
    """Load VEP annotation results for a tissue."""
    sample_vep = pd.read_table('{}/{}/{}.v8.leadvars.vep.vcf'.format(prefix, config['annotations_output_dir'], tissue_id), skiprows=4)
    overlap_df = load_overlap(config, tissue_id)
    return pd.merge(sample_vep, overlap_df, left_on='ID', right_on='lead_variant_id', how='outer')

def load_susie_annotated(config, tissue_id):
    """Load SuSiE annotated results for a tissue."""
    qtl_path = '{}/{}{}/{}.v8.susie_R_vars.annotated.csv'.format(prefix, config['annotations_output_dir'], tissue_id, tissue_id)
    susie_annotated = pd.read_table(qtl_path)
    susie_annotated['cs_num'] = susie_annotated['cs_id'] 
    susie_annotated['cs_id'] = susie_annotated['phenotype_id'] + '_cs_' + susie_annotated['cs_id'].astype(str)
    susie_annotated = add_lead_var(susie_annotated)
    return susie_annotated

def load_pc_annotated(config, tissue_id):
    """Load PC annotation results for a tissue."""
    return pd.read_table("{}/{}{}/{}.v8.pcs_annotated.txt".format(prefix, config["annotations_output_dir"], tissue_id, tissue_id))

def load_clusters_annotated(config, tissue_id):
    """Load annotated clusters for a tissue."""
    annot_cluster = pd.read_csv('{}/{}/{}/{}.clusters.annotated.txt'.format(prefix, config['annotations_output_dir'], tissue_id, tissue_id), index_col=0)
    annot_cluster['abs_cor'] = np.where(annot_cluster['mean_neg_corr'].isna(), annot_cluster['mean_pos_corr'], - annot_cluster['mean_neg_corr'])
    annot_cluster['abs_cor'] = np.where(annot_cluster['abs_cor'] < annot_cluster['mean_pos_corr'], annot_cluster['mean_pos_corr'], annot_cluster['abs_cor'])
    annot_cluster['max_pos_neg_cor'] = np.where(annot_cluster['mean_neg_corr'].isna(), annot_cluster['mean_pos_corr'], annot_cluster['mean_neg_corr'])
    annot_cluster['max_pos_neg_cor'] = np.where(abs(annot_cluster['max_pos_neg_cor']) < annot_cluster['mean_pos_corr'], annot_cluster['mean_pos_corr'], annot_cluster['max_pos_neg_cor'])
    return annot_cluster

def load_null_clusters_annotated(config, tissue_id, num_genes=2):
    """Load annotated null clusters for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.null_{}genes.annotated.txt'.format(prefix, config['annotations_output_dir'], tissue_id, tissue_id, num_genes), index_col=0)

def load_cluster(config, tissue_id):
    """Load cluster data for a tissue."""
    cluster_df =  pd.read_csv('{}/{}/{}.clusters.txt'.format(prefix, config['clusters_dir'], tissue_id),index_col=0)
    return cluster_df

def load_pc_cis(config, tissue_id):
    """Load pcQTL cis results for a tissue."""
    pc_cis_path = '{}/{}/{}/{}.v8.pcs.cis_qtl.txt.gz'.format(prefix, config['pcqtl_output_dir'], tissue_id, tissue_id)
    pc_cis_df = pd.read_table(pc_cis_path)
    pc_cis_df['cluster_id'] = pc_cis_df['phenotype_id'].str.split('_pc').str[0]
    annotate_pc_order(pc_cis_df)
    return pc_cis_df

def load_e_cis(config, tissue_id):
    """Load eQTL cis results for a tissue."""
    e_cis_path = '{}/{}/{}/{}.v8.cluster_genes.cis_qtl.txt.gz'.format(prefix, config['eqtl_output_dir'], tissue_id, tissue_id)
    e_cis_df = pd.read_table(e_cis_path)
    e_cis_df['cluster_id'] = e_cis_df['phenotype_id'].str.split('_e').str[0]
    return e_cis_df

def load_e_nominal(config, tissue_id, chr_id=22, get_var_position=False):
    """Load eQTL nominal results for a tissue and chromosome."""
    eqtl_output_dir = config['eqtl_output_dir']
    path = f'{prefix}/{eqtl_output_dir}/{tissue_id}/{tissue_id}.v8.cluster_genes.cis_qtl_pairs.chr{chr_id}.parquet'
    e_nominal_df = pd.read_parquet(path)
    if get_var_position:
        e_nominal_df['variant_pos'] = var_pos(e_nominal_df)
    e_nominal_df['cluster_id'] = e_nominal_df['phenotype_id'].str.split('_e_').str[0]
    e_nominal_df['var_cluster'] = e_nominal_df['variant_id'] + '_' + e_nominal_df['cluster_id']
    e_nominal_df['egene_id'] = e_nominal_df['phenotype_id'].str.split('_e_').str[1]
    return e_nominal_df

def load_pc_nominal(config, tissue_id, chr_id=22, get_var_position=False):
    """Load pcQTL nominal results for a tissue and chromosome."""
    pcqtl_output_dir = config['pcqtl_output_dir']
    path = f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.cis_qtl_pairs.chr{chr_id}.parquet'
    pc_nominal_df = pd.read_parquet(path)
    if get_var_position:
        pc_nominal_df['variant_pos'] = var_pos(pc_nominal_df)
    pc_nominal_df['cluster_id'] = pc_nominal_df['phenotype_id'].str.split('_pc').str[0]
    pc_nominal_df['var_cluster'] = pc_nominal_df['variant_id'] + '_' + pc_nominal_df['cluster_id']
    pc_nominal_df['pc_num'] = pc_nominal_df['phenotype_id'].str.split('_pc').str[-1].astype(int)
    return pc_nominal_df

def load_pc_susie_r(config, tissue_id):
    """Load pcQTL SuSiE R results for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.pcs.susie_R.txt'.format(prefix, config['pcqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_pc_susie(config, tissue_id, use_r=False):
    """Load pcQTL SuSiE results for a tissue."""
    if use_r:
        return load_pc_susie_r(config, tissue_id)
    else:
        return pd.read_csv('{}/{}/{}/{}.v8.pcs.susie.txt'.format(prefix, config['pcqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_e_susie_r(config, tissue_id):
    """Load eQTL SuSiE R results for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.cluster_genes.susie_R.txt'.format(prefix, config['eqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_e_susie(config, tissue_id, use_r=False):
    """Load eQTL SuSiE results for a tissue."""
    if use_r:
        return load_e_susie_r(config, tissue_id)
    else:
        return pd.read_csv('{}/{}/{}/{}.v8.cluster_genes.susie.txt'.format(prefix, config['eqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_pairwise_coloc(config, tissue_id):
    """Load pairwise co-localization results for a tissue."""
    return pd.read_csv('{}/{}/pairs/{}.v8.pairs_coloc.txt'.format(prefix, config['coloc_output_dir'], tissue_id), sep='\t')

def load_tissue_ids(config):
    """Load tissue IDs from configuration."""
    return pd.read_csv(config['tissue_id_path'])['Tissue'].tolist()

def load_tissue_df(config):
    """Load tissue information DataFrame."""
    return pd.read_csv(config['tissue_id_path'])

def load_pc_top_per_phenotype(config, tissue_id):
    """Load top pcQTL per phenotype for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.pcs.cis_qtl.txt.gz'.format(prefix, config['pcqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_e_top_per_phenotype(config, tissue_id):
    """Load top eQTL per phenotype for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.cluster_genes.cis_qtl.txt.gz'.format(prefix, config['eqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_pc_permutation(config, tissue_id):
    """Load pcQTL permutation results for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.pcs.cis_independent_qtl.txt.gz'.format(prefix, config['pcqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_e_permutation(config, tissue_id):
    """Load eQTL permutation results for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.cluster_genes.cis_independent_qtl.txt.gz'.format(prefix, config['eqtl_output_dir'], tissue_id, tissue_id), sep='\t')

def load_e_nominal_all_chr(config, tissue_id):
    """Load eQTL nominal results for all chromosomes for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.cluster_genes.cis_qtl_pairs.all_chr.parquet'.format(prefix, config['eqtl_output_dir'], tissue_id, tissue_id))

def load_pc_nominal_all_chr(config, tissue_id):
    """Load pcQTL nominal results for all chromosomes for a tissue."""
    return pd.read_csv('{}/{}/{}/{}.v8.pcs.cis_qtl_pairs.all_chr.parquet'.format(prefix, config['pcqtl_output_dir'], tissue_id, tissue_id))

def load_expression(config, tissue_id):
    """Load expression data for a tissue."""
    return pd.read_csv('{}/{}/{}.v8.normalized_expression.bed'.format(prefix, config['expression_dir'], tissue_id), sep='\t')

def load_cluster_expression(config, tissue_id):
    """Load cluster expression data for a tissue."""
    return pd.read_csv('{}/{}/{}.v8.normalized_residualized_expression.cluster_genes.bed'.format(prefix, config['filtered_expression_output_dir'], tissue_id), sep='\t')

def load_pc(config, tissue_id):
    """Load principal components for a tissue."""
    return pd.read_csv('{}/{}/{}.pcs.bed'.format(prefix, config['pc_output_dir'], tissue_id), sep='\t')

def load_qtl_signal_groups(config, tissue_id):
    """Load QTL signal groups for a tissue."""
    qtl_signal_groups_path = '{}/{}/qtl_signal_groups/{}.qtl_signal_groups.txt'.format(prefix, config['coloc_output_dir'], tissue_id)
    return pd.read_csv(qtl_signal_groups_path, sep='\t')

def load_gwas_signal_groups(config, tissue_id):
    """Load GWAS signal groups for a tissue."""
    gwas_signal_groups_path = '{}/{}/gwas_signal_groups/{}.gwas_signal_groups.txt'.format(prefix, config['coloc_output_dir'], tissue_id)
    return pd.read_csv(gwas_signal_groups_path, sep='\t')

# Utility functions for data processing
def var_pos(df, column='variant_id'):
    """Extract variant positions from variant IDs."""
    return df[column].str.split('_').str[1].astype(int)

def add_lead_var(susie_df):
    """Add lead variant information to SuSiE results."""
    susie_df['lead_variant_id'] = susie_df['variant_id']
    return susie_df

def add_num_vars_cs(susie_df):
    """Add number of variants per credible set to SuSiE results."""
    susie_df['num_vars_cs'] = susie_df.groupby('cs_id')['variant_id'].transform('count')
    return susie_df

def load_across_tissues(config, load_func, tissue_ids=None):
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

def annotate_pc_order(pc_df):
    """Annotate PC order in pcQTL results."""
    pc_df['pc_num'] = pc_df['phenotype_id'].str.split('_pc').str[-1].astype(int)
    pc_df['pc_order'] = 'middle'
    pc_df.loc[pc_df['pc_num'] == 1, 'pc_order'] = 'first'
    pc_df.loc[pc_df['pc_num'] == pc_df['num_genes'], 'pc_order'] = 'last'
    return pc_df

def remove_cross_map(df, cross_map_ids=None, config=None):
    """Remove cross-mappable genes from results."""
    if cross_map_ids is None and config is not None:
        # Load cross-mappable genes if not provided
        cross_map_df = pd.read_csv(config['cross_map_path'])
        cross_map_ids = cross_map_df['gene_id'].tolist()
    
    if cross_map_ids is not None:
        df = df[~df['gene_id'].isin(cross_map_ids)]
    
    return df