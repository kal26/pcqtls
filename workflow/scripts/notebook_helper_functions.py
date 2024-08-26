f# functions to load the requested file types given a config pointing the right directories
import pandas as pd

# working directory prefix
prefix = '/home/klawren/oak/pcqtls'

def load_overlap(config, tissue_id):
    overlap_df = pd.read_csv('{}/{}/{}.v8.overlap.txt'.format(prefix, config['overlap_output_dir'], tissue_id), sep='\t')
    overlap_df['var_cluster'] = overlap_df['lead_variant_id'] + '_' + overlap_df['cluster_id']
    return overlap_df

def load_clusters_annotated(config, tissue_id):
    return pd.read_csv('{}/{}/{}_clusters_annotated.csv'.format(prefix, config['annotations_output_dir'], tissue_id), index_col=0)

def load_null_clusters_annotated(config, tissue_id, num_genes=2):
    return pd.read_csv('{}/{}/{}_null_{}genes_annotated.csv'.format(prefix, config['annotations_output_dir'], tissue_id, num_genes), index_col=0)

def load_cluster(config, tissue_id):
    cluster_df =  pd.read_csv('{}/{}/{}_clusters_all_chr.csv'.format(prefix, config['clusters_dir'], tissue_id),index_col=0)
    for idx, row in cluster_df.iterrows():
        cluster_df.loc[idx, 'cluster_id_sort']  = '_'.join([*sorted(row['Transcripts'].split(','))])
    return cluster_df

# load in e nominal
def load_e_nominal(config, tissue_id, chr_id=22, get_var_position=False):
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
    pcqtl_output_dir = config['pcqtl_output_dir']
    path = f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.cis_qtl_pairs.chr{chr_id}.parquet'
    pc_nominal_df = pd.read_parquet(path)
    if get_var_position:
        pc_nominal_df['variant_pos'] = var_pos(pc_nominal_df)
    pc_nominal_df['cluster_id'] = pc_nominal_df['phenotype_id'].str[:-4]
    pc_nominal_df['var_cluster'] = pc_nominal_df['variant_id'] + '_' + pc_nominal_df['cluster_id']
    return pc_nominal_df

def load_pc_susie(config, tissue_id):
    pcqtl_output_dir = config['pcqtl_output_dir']
    pc_susie_df =  pd.read_csv(f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.susie.txt', sep='\t', index_col=0)
    pc_susie_df['var_cluster'] = pc_susie_df['variant_id'] + '_' + pc_susie_df['phenotype_id'].str.split('_pc1').str[0]
    pc_susie_df['cs_num'] = pc_susie_df['cs_id'] 
    pc_susie_df['cs_id'] = pc_susie_df['phenotype_id'] + '_' + pc_susie_df['cs_id'].astype(str)
    pc_susie_df['pc_num'] = pc_susie_df['phenotype_id'].str.split('_pc').str[-1].astype(int)
    pc_susie_df['cluster_id'] = pc_susie_df['phenotype_id'].str.split('_pc').str[0]
    pc_susie_df = add_lead_var(pc_susie_df)
    pc_susie_df = add_num_vars_cs(pc_susie_df)
    return pc_susie_df


def load_tissue_ids(config):
    tissue_id_path = config['tissue_id_path']
    tissue_df = pd.read_csv(f"{prefix}/{tissue_id_path}", header=0)
    return list(tissue_df['Tissue'])

def load_tissue_df(config):
    tissue_id_path = config['tissue_id_path']
    return pd.read_csv(f"{prefix}/{tissue_id_path}", header=0)


def load_pc_top_per_phenotype(config, tissue_id):
    pcqtl_output_dir = config['pcqtl_output_dir']
    return pd.read_csv(f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.cis_qtl.txt.gz', sep='\t')

def load_e_top_per_phenotype(config, tissue_id):
    eqtl_output_dir = config['eqtl_output_dir']
    return pd.read_csv(f'{prefix}/{eqtl_output_dir}/{tissue_id}/{tissue_id}.v8.cluster_genes.cis_qtl.txt.gz', sep='\t')

def load_pc_permutation(config, tissue_id, pc1_only=False):
    pcqtl_output_dir = config['pcqtl_output_dir']
    if pc1_only:
        pc_path = f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pc1_only.cis_independent_qtl.txt.gz'
    else:
        pc_path = f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.cis_independent_qtl.txt.gz'
    pc_permutation_df = pd.read_csv(pc_path, sep='\t')
    pc_permutation_df['cluster_id'] = pc_permutation_df['phenotype_id'].str.split('_pc').str[0]
    return pc_permutation_df

def load_e_permutation(config, tissue_id):
    eqtl_output_dir = config['eqtl_output_dir']
    e_permutation_df = pd.read_csv(f'{prefix}/{eqtl_output_dir}/{tissue_id}/{tissue_id}.v8.cluster_genes.cis_independent_qtl.txt.gz', sep='\t')
    e_permutation_df['cluster_id'] = e_permutation_df['phenotype_id'].str.split('_e').str[0]
    return e_permutation_df


def load_e_nominal_all_chr(config, tissue_id):
    e_nominal_dfs=[]
    for chr_id in range(1,23):
        e_nominal_dfs.append(load_e_nominal(config, tissue_id, chr_id=chr_id))
    return pd.concat(e_nominal_dfs)

def load_pc_nominal_all_chr(config, tissue_id):
    pc_nominal_dfs=[]
    for chr_id in range(1,23):
        pc_nominal_dfs.append(load_pc_nominal(config, tissue_id, chr_id=chr_id))
    return pd.concat(pc_nominal_dfs)



# functions to help with annotating loaded data
def var_pos(df, column='variant_id'):
    return df[column].str.split('_').str[1].astype(int)

def add_lead_var(susie_df):
    lead_vars = susie_df.loc[susie_df.groupby('cs_id')['pip'].idxmax(),['cs_id','variant_id']].set_index('cs_id')
    susie_df = susie_df.merge(lead_vars, how='left', left_on='cs_id', right_index=True)
    susie_df = susie_df.rename(columns={'variant_id_y':'lead_variant_id', 'variant_id_x':'variant_id'})
    return susie_df

def add_num_vars_cs(susie_df):
    num_vars = susie_df.groupby('cs_id').agg({'variant_id':'count'})
    num_vars = num_vars.rename(columns={'variant_id':'num_vars'})
    susie_df = susie_df.merge(num_vars, how='left', left_on='cs_id', right_index=True)
    return susie_df