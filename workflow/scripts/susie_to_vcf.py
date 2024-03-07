import pandas as pd


susie_path = snakemake.input[0]
out_path = snakemake.output[0]

susie_df = pd.read_csv(susie_path, sep='\t', index_col=0)

susie_df['cs_full_id'] = susie_df['phenotype_id'].astype(str) + '_e_cs' + susie_df['cs_id'].astype(str) 

def get_lead_var(susie_df):
    return susie_df.loc[susie_df.groupby('cs_full_id')['pip'].idxmax(),['cs_full_id','variant_id']].set_index('cs_full_id')


leadvar_df = pd.DataFrame(pd.Series(susie_df.groupby(['cs_full_id'])['variant_id'].apply(list), name='variant_list'))
leadvar_df['lead_variant_id'] = get_lead_var(susie_df)

leadvar_df['chr'] = leadvar_df['lead_variant_id'].str.split('_').str[0].str[3:]
leadvar_df['pos'] = leadvar_df['lead_variant_id'].str.split('_').str[1]
leadvar_df['ref'] = leadvar_df['lead_variant_id'].str.split('_').str[2]
leadvar_df['alt'] = leadvar_df['lead_variant_id'].str.split('_').str[3]
leadvar_df['blank'] = '.'

col_list = ['chr', 'pos', 'lead_variant_id', 'ref', 'alt', 'blank', 'blank', 'blank']
out_df = leadvar_df[col_list]
out_df = out_df.sort_values(by=['chr', 'pos'])

out_df.to_csv(out_path, header=None, index=False, sep='\t')