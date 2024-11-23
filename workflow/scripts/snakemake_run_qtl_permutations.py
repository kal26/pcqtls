import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis, genotypeio

phenotype_path = snakemake.input[1]
covariates_path = snakemake.input[2]
cis_results_path = snakemake.input[3]

genotype_stem = snakemake.params[0]
tissue = snakemake.params[1]
pc1_only = snakemake.params[2]

outpath = snakemake.output[0]


# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_path)
covariates_df = pd.read_csv(covariates_path, sep='\t', index_col=0).T

# PLINK reader for genotypes
pgr = genotypeio.PlinkReader(genotype_stem)
genotype_df = pgr.load_genotypes()
variant_df = pgr.bim.set_index('snp')[['chrom', 'pos']]

# load in cis results
cis_df = pd.read_csv(cis_results_path, sep='\t', index_col=0)
# filter to only pc1 if requested and recalculate q values
if pc1_only:
    cis_df = cis_df[cis_df.index.str.contains('pc1')]

# calculate q values    
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)


try: 
    # run permutations
    indep_df = cis.map_independent(genotype_df, variant_df, cis_df,
                                phenotype_df, phenotype_pos_df, 
                                covariates_df)

    # write out 
    indep_df.to_csv(outpath, sep='\t', index=False, compression='gzip')
except ValueError as e:
    print(e)
    print('writing out blank file')
    pd.DataFrame(columns=['phenotype_id', 'num_var', 'beta_shape1', 'beta_shape2', 'true_df',
        'pval_true_df', 'variant_id', 'start_distance', 'end_distance',
        'ma_samples', 'ma_count', 'af', 'pval_nominal', 'slope', 'slope_se',
        'pval_perm', 'pval_beta', 'rank']).to_csv(outpath, sep='\t', index=False, compression='gzip')