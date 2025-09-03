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
    annot_cluster = pd.read_table('{}/{}/annotated_clusters/{}.clusters.annotated.txt'.format(config["working_dir"], config['annotations_output_dir'], tissue_id))
    return annot_cluster

def load_null_clusters_annotated(config: Dict[str, str], tissue_id: str, num_genes: int = 2) -> pd.DataFrame:
    """Load annotated null clusters for a tissue."""
    return pd.read_table('{}/{}/annotated_null_clusters/{}.null.{}genes.annotated.txt'.format(config["working_dir"], config['annotations_output_dir'], tissue_id, num_genes))

def load_cluster(config: Dict[str, str], tissue_id: str) -> pd.DataFrame:
    """Load cluster data for a tissue."""
    cluster_df =  pd.read_table('{}/{}/{}.clusters.txt'.format(config["working_dir"], config['clusters_dir'], tissue_id),index_col=0)
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

def load_gwas_signal_groups(config: Dict[str, str], tissue_id: str, use_susie: bool = True) -> pd.DataFrame:
    """Load GWAS signal groups for a tissue."""
    if use_susie:
        gwas_signal_groups_path = '{}/{}/gwas_signal_groups_susie_True/{}.gwas_signal_groups.txt'.format(config["working_dir"], config['coloc_output_dir'], tissue_id)
    else:
        gwas_signal_groups_path = '{}/{}/gwas_signal_groups_susie_False/{}.gwas_signal_groups.txt'.format(config["working_dir"], config['coloc_output_dir'], tissue_id)
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


def calculate_enrichment_with_covariates(clusters_df, null_df, annotation_cols, min_group_size=10, min_expected_freq=5):
    """
    Calculate enrichment odds ratios with covariates using logistic regression
    Skip only when sample sizes or expected frequencies are too small for reliable statistics
    """
    results = []
    
    # Combine datasets and add labels
    clusters_df = clusters_df.copy()
    null_df = null_df.copy()
    
    clusters_df['is_cluster'] = 1
    null_df['is_cluster'] = 0
    
    combined_df = pd.concat([clusters_df, null_df], ignore_index=True)
    
    for annotation in annotation_cols:
        try:
            # Prepare data for logistic regression
            # Remove rows with missing values for this annotation
            analysis_df = combined_df.dropna(subset=[annotation, 'cluster_length', 'num_genes'])
            
            if len(analysis_df) == 0:
                print(f"No data available for {annotation}")
                results.append({
                    'annotation': annotation,
                    'odds_ratio': np.nan,
                    'ci_lower': np.nan,
                    'ci_upper': np.nan,
                    'p_value': np.nan,
                    'reason': 'no_data'
                })
                continue
            
            # Check group sizes
            cluster_group = analysis_df[analysis_df['is_cluster'] == 1]
            null_group = analysis_df[analysis_df['is_cluster'] == 0]
            
            if len(cluster_group) < min_group_size or len(null_group) < min_group_size:
                print(f"Insufficient group sizes for {annotation}: clusters={len(cluster_group)}, null={len(null_group)}")
                results.append({
                    'annotation': annotation,
                    'odds_ratio': np.nan,
                    'ci_lower': np.nan,
                    'ci_upper': np.nan,
                    'p_value': np.nan,
                    'reason': f'small_groups_clusters_{len(cluster_group)}_null_{len(null_group)}'
                })
                continue
            
            # Check expected frequencies for 2x2 table (Fisher's exact test assumptions)
            cluster_pos = len(cluster_group[cluster_group[annotation] == 1])
            cluster_neg = len(cluster_group[cluster_group[annotation] == 0])
            null_pos = len(null_group[null_group[annotation] == 1])
            null_neg = len(null_group[null_group[annotation] == 0])
            
            # Calculate expected frequencies under null hypothesis
            total_pos = cluster_pos + null_pos
            total_neg = cluster_neg + null_neg
            total_cluster = len(cluster_group)
            total_null = len(null_group)
            total_all = len(analysis_df)
            
            exp_cluster_pos = (total_cluster * total_pos) / total_all
            exp_cluster_neg = (total_cluster * total_neg) / total_all
            exp_null_pos = (total_null * total_pos) / total_all
            exp_null_neg = (total_null * total_neg) / total_all
            
            min_expected = min(exp_cluster_pos, exp_cluster_neg, exp_null_pos, exp_null_neg)
            
            if min_expected < min_expected_freq:
                print(f"Low expected frequency for {annotation}: min_expected={min_expected:.2f}")
                results.append({
                    'annotation': annotation,
                    'odds_ratio': np.nan,
                    'ci_lower': np.nan,
                    'ci_upper': np.nan,
                    'p_value': np.nan,
                    'reason': f'low_expected_freq_{min_expected:.2f}'
                })
                continue
            
            # Check for complete separation (all positive in one group, all negative in another)
            if cluster_pos == 0 or cluster_neg == 0 or null_pos == 0 or null_neg == 0:
                print(f"Complete or quasi-complete separation for {annotation}")
                results.append({
                    'annotation': annotation,
                    'odds_ratio': np.nan,
                    'ci_lower': np.nan,
                    'ci_upper': np.nan,
                    'p_value': np.nan,
                    'reason': 'complete_separation'
                })
                continue
                
            # Define dependent and independent variables
            y = analysis_df[annotation]
            X = analysis_df[['is_cluster', 'cluster_length', 'num_genes']]
            X = sm.add_constant(X)  # Add intercept
            
            # Fit logistic regression
            model = sm.Logit(y, X)
            result = model.fit(disp=0)
            
            # Check if model converged
            if not result.mle_retvals['converged']:
                print(f"Model did not converge for {annotation}")
                results.append({
                    'annotation': annotation,
                    'odds_ratio': np.nan,
                    'ci_lower': np.nan,
                    'ci_upper': np.nan,
                    'p_value': np.nan,
                    'reason': 'no_convergence'
                })
                continue
            
            # Extract odds ratio and confidence interval for 'is_cluster'
            coef = result.params['is_cluster']
            odds_ratio = np.exp(coef)
            conf_int = np.exp(result.conf_int().loc['is_cluster'])
            
            results.append({
                'annotation': annotation,
                'odds_ratio': odds_ratio,
                'ci_lower': conf_int[0],
                'ci_upper': conf_int[1],
                'p_value': result.pvalues['is_cluster'],
                'reason': 'success',
                'contingency_table': f"cluster_pos:{cluster_pos}, cluster_neg:{cluster_neg}, null_pos:{null_pos}, null_neg:{null_neg}"
            })
            
        except Exception as e:
            print(f"Error calculating enrichment for {annotation}: {e}")
            results.append({
                'annotation': annotation,
                'odds_ratio': np.nan,
                'ci_lower': np.nan,
                'ci_upper': np.nan,
                'p_value': np.nan,
                'reason': f'error_{str(e)[:50]}'
            })
            continue
    return pd.DataFrame(results)

def categorize_correlation_type(df):
    """ Categorize clusters by correlation type """
    df = df.copy()
    conditions = [
        df['mean_neg_corr'].isna() & df['mean_pos_corr'].notna(),
        df['mean_pos_corr'].isna() & df['mean_neg_corr'].notna(),
        df['mean_pos_corr'].notna() & df['mean_neg_corr'].notna()
    ]
    choices = ['positive_only', 'negative_only', 'both_corr']
    df['corr_type'] = np.select(conditions, choices, default='unknown')
    return df