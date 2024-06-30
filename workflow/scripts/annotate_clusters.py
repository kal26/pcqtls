import numpy as np
import pandas as pd

# params
protien_coding_only = snakemake.params[0]

# paths to data for this cluster
cluster_path = snakemake.input[0]
expression_path = snakemake.input[1]

# paths to external data sources
gencode_path = snakemake.input[2]


# load in data
cluster_df = pd.read_csv(cluster_path, index_col=0)
expression_df = pd.read_csv(expression_path, sep='\t')
full_gencode = pd.read_csv(gencode_path)

# filter to protien coding
if protien_coding_only:
    non_protein_gencode = full_gencode.copy()
    full_gencode = full_gencode[full_gencode['gene_type'] == 'protein_coding']

# expressed genes in sample tissue
expressed_gencode = full_gencode[full_gencode['transcript_id'].isin(expression_df['gene_id'])]
expressed_gencode = expressed_gencode.sort_values(['chr', 'start', 'end'])


# load in abc
