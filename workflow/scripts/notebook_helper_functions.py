# functions to load the requested file types given a config pointing the right directories
import pandas as pd

# working directory prefix
prefix = '/home/klawren/oak/pcqtls'

def load_overlap(config, tissue_id):
    return pd.read_csv('{}/{}/{}.v8.overlap.txt'.format(prefix, config['overlap_output_dir'], tissue_id), sep='\t')

def load_clusters_annotated(config, tissue_id):
    return pd.read_csv('{}/{}/{}_clusters_annotated.csv'.format(prefix, config['annotations_output_dir'], tissue_id), index_col=0)

def load_null_clusters_annotated(config, tissue_id, num_genes=2):
    return pd.read_csv('{}/{}/{}_null_{}genes_annotated.csv'.format(prefix, config['annotations_output_dir'], tissue_id, num_genes), index_col=0)


# load in e nominal
def load_e_nominal(config, tissue_id, chr_id=22, get_var_position=False):
    eqtl_output_dir = config['eqtl_output_dir']
    path = f'{prefix}/{eqtl_output_dir}/{tissue_id}/{tissue_id}.v8.cluster_genes.cis_qtl_pairs.chr{chr_id}.parquet'
    e_nominal_df = pd.read_parquet(path)
    if get_var_position:
        e_nominal_df['variant_pos'] = var_pos(e_nominal_df)
    e_nominal_df['cluster_id'] = e_nominal_df['phenotype_id'].str.split('_e_').str[0]
    return e_nominal_df

def load_pc_nominal(config, tissue_id, chr_id=22, get_var_position=False):
    pcqtl_output_dir = config['pcqtl_output_dir']
    path = f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.cis_qtl_pairs.chr{chr_id}.parquet'
    pc_nominal_df = pd.read_parquet(path)
    if get_var_position:
        pc_nominal_df['variant_pos'] = var_pos(pc_nominal_df)
    pc_nominal_df['cluster_id'] = pc_nominal_df['phenotype_id'].str[:-4]
    return pc_nominal_df

def load_susie(config, tissue_id):
    pcqtl_output_dir = config['pcqtl_output_dir']
    return pd.read_csv(f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.susie.txt', sep='\t', index_col=0)

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

def load_pc_permutation(config, tissue_id):
    pcqtl_output_dir = config['pcqtl_output_dir']
    return pd.read_csv(f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.cis_independent_qtl.txt.gz', sep='\t')

def load_e_permutation(config, tissue_id):
    eqtl_output_dir = config['eqtl_output_dir']
    return pd.read_csv(f'{prefix}/{eqtl_output_dir}/{tissue_id}/{tissue_id}.v8.cluster_genes.cis_independent_qtl.txt.gz', sep='\t')


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
def var_pos(df):
    return df['variant_id'].str.split('_').str[1].astype(int)