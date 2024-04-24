# load data
import pandas as pd
import numpy as np

# load in data

expression_df = pd.read_csv(snakemake.input[0], sep='\t')
gtex_covariates_df = pd.read_csv(snakemake.input[1], sep='\t', index_col=0).T

# get covariates
# TODO
