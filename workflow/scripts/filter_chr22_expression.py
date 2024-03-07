import pandas as pd
import numpy as np

# load in date
expression_df = pd.read_csv(snakemake.input[0], sep='\t')

out_df = expression_df[(expression_df['#chr'] == 'chr22') | (expression_df['#chr'] == 'chr2')| (expression_df['#chr'] == 'chr1')]

# write out
out_df = out_df.sort_values(by=['#chr', 'start'])
out_df.to_csv(snakemake.output[0], sep='\t', index=False)
