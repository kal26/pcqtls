import pandas as pd
import numpy as np

# load data
expression_df = pd.read_csv(snakemake.input[0], sep='\t')
covariates_df = pd.read_csv(snakemake.input[1], sep='\t', index_col=0).T

use_scramble_order = snakemake.params[0]

if use_scramble_order:
    selected_samples = covariates_df.sample(num_samples).index.values
else:
    selected_samples = covariates_df.index[:num_samples].values

# write out first x as a subset
sub_expression = expression_df[np.concatenate([expression_df.columns[:4].values, selected_samples])]
sub_expression.to_csv(snakemake.output[0], sep='\t', index=None)

sub_covar = covariates_df.loc[selected_samples]
sub_covar.T.to_csv(snakemake.output[1], sep='\t')
