import pandas as pd
from tqdm import tqdm 
import annotate_qtls

pc_path = snakemake.input[0]
expression_path = snakemake.input[1]
out_path = snakemake.output[0]


# read in data 
pc_df = pd.read_csv(pc_path, sep='\t')
pc_df['cluster_id'] = pc_df['gene_id'].str.split('_pc').str[0]
pc_df['pc_id'] = pc_df['gene_id'].str.split('_pc').str[1].astype('float')
pc_df['cluster_size'] = pc_df['cluster_id'].str.split('_').apply(len)

expression_df = pd.read_csv(expression_path, sep='\t')
expression_df['cluster_id'] = expression_df['gene_id'].str.split('_e_').str[0]
expression_df['egene_id'] = expression_df['gene_id'].str.split('_e_').str[1]

# annotate the pcs
print(pc_df.head())
print(expression_df.head())
annotated_pcs = annotate_qtls.get_annotate_pcs(pc_df, expression_df)

# write out results
annotated_pcs.to_csv(out_path, sep='\t', index=None)

