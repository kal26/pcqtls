import pandas as pd
eqtl = pd.read_parquet(snakemake.input[0])
cluster_list = pd.read_csv(snakemake.input[1], header=None, names=['cluster_id'])
eqtl['cluster_id'] = eqtl['phenotype_id'].str.split('_e').str[0]
tissue_id = snakemake.params[0]
chrom = snakemake.params[1]
output_dir = snakemake.params[2]
for idx, row in cluster_list.iterrows():
    snp_list = eqtl[eqtl['cluster_id'] == row['cluster_id']]['variant_id']
    snp_list.to_csv(f"{output_dir}/{tissue_id}.{chrom}.{row['cluster_id']}.snp_list.txt", index=False, sep='\t')
