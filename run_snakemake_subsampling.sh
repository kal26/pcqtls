#!/usr/bin/env bash
#conda activate snakemake_7.32

for i in 200 250 300 350 400 450
do
    echo "Running command for $i:"
    snakemake --snakefile workflow/Snakefile_subsampling --configfile config/subsampling/${i}_subsample.yaml --profile config/slurm_scg
done