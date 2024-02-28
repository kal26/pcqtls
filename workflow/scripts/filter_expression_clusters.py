import pandas as pd
import numpy as np

# load in date
cluster_df = pd.read_csv(snakemake.input[0])
expression_df = pd.read_csv(snakemake.input[1], sep='\t')

# set gene_id as index
expression_df.set_index('gene_id', inplace=True)



# filter to just the genes that are in the clusters
# change gene_id to be the cluster id_e_gene_id for each gene
# expand the search window to be the same size for all genes in the cluster
cluster_sets = []
for idx, row in cluster_df.iterrows():
    cluster_set = expression_df.loc[row['Transcripts'].split(',')].copy()
    cluster_set['start'] = cluster_set['start'].min()
    cluster_set['end'] = cluster_set['end'].max()
    cluster_set.insert(3, 'gene_id', ['_'.join([*sorted(cluster_set.index.values), 'e', gid]) for gid in cluster_set.index])
    cluster_sets.append(cluster_set)

# write out
out_df = pd.concat(cluster_sets)
out_df = out_df.sort_values(by=['#chr', 'start'])


out_df.to_csv(snakemake.output[0], sep='\t', index=False)
