import pandas as pd

# load in eqtl and get cluster id
eqtl = pd.read_parquet(snakemake.input[0])
eqtl['cluster_id'] = eqtl['phenotype_id'].str.split('_e').str[0]

tissue_id = snakemake.params[0]
chrom = snakemake.params[1]
output_dir = snakemake.params[2]

for cluster_id in eqtl['cluster_id'].unique():
    sub_eqtl = eqtl[eqtl['cluster_id'] == cluster_id]
    # just get vars for the first qtl of the cluster, as all have the same tested vars
    snp_list = sub_eqtl[sub_eqtl['phenotype_id'] == sub_eqtl['phenotype_id'].unique()[0]]['variant_id']
    outpath = f"{output_dir}/{tissue_id}.{chrom}.cluster_{row['cluster_id']}.snp_list.txt"
    print(outpath)
    snp_list.to_csv(outpath, index=False, sep='\t')
pd.DataFrame({}).to_csv(snakemake.output[0])
