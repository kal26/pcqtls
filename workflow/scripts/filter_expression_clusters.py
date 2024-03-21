import pandas as pd
import numpy as np
from residualize import calculate_residual

# load in data
cluster_df = pd.read_csv(snakemake.input[0])
expression_df = pd.read_csv(snakemake.input[1], sep='\t')
covariates_df = pd.read_csv(snakemake.input[2], sep='\t', index_col=0).T


# set gene_id as index
expression_df_gid = expression_df.set_index('gene_id')

# residulize the expression 
residal_exp = calculate_residual(expression_df[expression_df.columns[4:]], covariates_df, center=True)
expression_df_res = expression_df_gid.copy()
expression_df_res[expression_df.columns[4:]] = residal_exp


# filter to just the genes that are in the clusters
# change gene_id to be the cluster id_e_gene_id for each gene
# expand the search window to be the same size for all genes in the cluster
cluster_sets = []
for idx, row in cluster_df.iterrows():
    cluster_set = expression_df_res.loc[row['Transcripts'].split(',')].copy()
    cluster_set['start'] = cluster_set['start'].min()
    cluster_set['end'] = cluster_set['end'].max()
    cluster_set.insert(3, 'gene_id', ['_'.join([*sorted(cluster_set.index.values), 'e', gid]) for gid in cluster_set.index])
    cluster_sets.append(cluster_set)


# combine the dfs
out_df = pd.concat(cluster_sets)
# sort (required by tensorqtl)
out_df = out_df.sort_values(by=['#chr', 'start'])
# write out
out_df.to_csv(snakemake.output[0], sep='\t', index=False)
