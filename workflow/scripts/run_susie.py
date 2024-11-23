
from tensorqtl import susie, genotypeio
import tensorqtl
import numpy as np
import pandas as pd
import argparse



def main():
    # Parse arguments from cmd
    parser = argparse.ArgumentParser()
    parser.add_argument('genotype_path', help = 'path to genotypes')
    parser.add_argument('phenotype_path', help = 'path to .bed normalized expression or pcs')
    parser.add_argument('susie_out_path', help = 'path where susie results should be written out')
    parser.add_argument('covariates_path', help = 'path to covariates')
    parser.add_argument('--verbosity', type=int, default=0, help = 'output verbosity')

    args = parser.parse_args()

    # load in data
    # load phenotypes and covariates
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(args.phenotype_path)
    covariates_df = pd.read_csv(args.covariates_path, sep='\t', index_col=0).T

    # PLINK reader for genotypes
    pgr = genotypeio.PlinkReader(args.genotype_path)
    genotype_df = pgr.load_genotypes()
    variant_df = pgr.bim.set_index('snp')[['chrom', 'pos']]

    print('run_susie')
    susie_df = susie.map(genotype_df, variant_df, 
                    phenotype_df, 
                    phenotype_pos_df, 
                    covariates_df)
    # write out
    print('writing out')
    print(susie_df)

    if len(susie_df)==0:
        # if there are no credible sets, the output of susie.map is a list, 
        # and cannot be properly written out 
        print('No susie credible sets')
        pd.DataFrame(columns=['idx', 'phenotype_id', 'variant_id', 'pip', 'af', 'cs_id']).to_csv(args.susie_out_path, sep='\t')
    else:
        susie_df.to_csv(args.susie_out_path, sep='\t')
        

if __name__ == "__main__":
    main()