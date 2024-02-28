import pandas as pd



susie_path = snakemake.input[0]
out_path = snakemake.output[0]

susie_df = pd.read_csv(susie_path, sep='\t', index_col=0)

susie_df['chr'] = susie_df['variant_id'].str.split('_').str[0].str[3:]
susie_df['pos'] = susie_df['variant_id'].str.split('_').str[1]
susie_df['ref'] = susie_df['variant_id'].str.split('_').str[2]
susie_df['alt'] = susie_df['variant_id'].str.split('_').str[3]
susie_df['blank'] = '.'

col_list = ['chr', 'pos', 'variant_id', 'ref', 'alt', 'blank', 'blank', 'blank']
out_df = susie_df[col_list]
out_df = out_df.sort_values(by=['chr', 'pos'])

out_df.to_csv(out_path, header=None, index=False, sep='\t')