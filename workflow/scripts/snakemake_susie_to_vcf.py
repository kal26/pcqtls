import pandas as pd


e_susie_path = snakemake.input[0]
pc_susie_path = snakemake.input[1]
out_path = snakemake.output[0]

pc_susie_df = pd.read_csv(pc_susie_path, sep='\t', index_col=0)
e_susie_df =  pd.read_csv(e_susie_path, sep='\t', index_col=0)

# get all variants in a credible set
susie_vars = pd.DataFrame(pd.concat([e_susie_df, pc_susie_df])['variant_id'].drop_duplicates())
susie_vars['chr'] = susie_vars['variant_id'].str.split('_').str[0].str[3:]
susie_vars['pos'] = susie_vars['variant_id'].str.split('_').str[1]
susie_vars['ref'] = susie_vars['variant_id'].str.split('_').str[2]
susie_vars['alt'] = susie_vars['variant_id'].str.split('_').str[3]
susie_vars['blank'] = '.'

col_list = ['chr', 'pos', 'variant_id', 'ref', 'alt', 'blank', 'blank', 'blank']
out_df = susie_vars[col_list]
out_df = out_df.sort_values(by=['chr', 'pos'])

out_df.to_csv(out_path, header=None, index=False, sep='\t')