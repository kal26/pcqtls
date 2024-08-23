import pandas as pd
import os

genotype_stem = snakemake.params[0]
cis_output_dir = snakemake.params[1]
tissue = snakemake.params[2]

expression_path = snakemake.input[1]
covariates_path = snakemake.input[2]
cis_results_path = snakemake.input[3]



os.system(f"python -m tensorqtl {genotype_stem} \
            {expression_path} \
            {cis_output_dir}{tissue}/{tissue}.v8.cluster_genes \
            --covariates {covariates_path} \
            --cis_output {cis_results_path} \
            --mode cis_independent \
            --maf_threshold .01")

outpath = snakemake.output[0]
# error raised when there are no signifigant results, and no file is written out
if not os.path.isfile(outpath):    
    print('writing out blank file')
    pd.DataFrame(columns=['phenotype_id', 'num_var', 'beta_shape1', 'beta_shape2', 'true_df',
       'pval_true_df', 'variant_id', 'start_distance', 'end_distance',
       'ma_samples', 'ma_count', 'af', 'pval_nominal', 'slope', 'slope_se',
       'pval_perm', 'pval_beta', 'rank']).to_csv(outpath, sep='\t', index=False)