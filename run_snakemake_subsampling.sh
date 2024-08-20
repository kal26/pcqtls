#!/usr/bin/env bash
#conda activate snakemake_7.32

for i in 50 75 100 200 300
do
    echo "Running command for $i:"
    snakemake --snakefile workflow/Snakefile_subsampling --configfile config/subsampling/${i}_subsample.yaml --profile config/slurm_scg
done