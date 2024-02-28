import pandas as pd

# reformat the lead variants from the credible sets 
# vep can then run annotation enrichment


df_path = snakemake.input[0]
out_path = snakemake.output[0]

df = pd.read_csv(df_path, sep='\t', index_col=0)

# drop duplicate lead variants
df = df.drop_duplicates(subset='lead_variant_id')

# pull vcf information
df['chr'] = df['lead_variant_id'].str.split('_').str[0].str[3:]
df['pos'] = df['lead_variant_id'].str.split('_').str[1]
df['ref'] = df['lead_variant_id'].str.split('_').str[2]
df['alt'] = df['lead_variant_id'].str.split('_').str[3]
df['blank'] = '.'

col_list = ['chr', 'pos', 'lead_variant_id', 'ref', 'alt', 'blank', 'blank', 'blank']
out_df = df[col_list]
out_df = out_df.sort_values(by=['chr', 'pos'])

out_df.to_csv(out_path, header=None, index=False, sep='\t')