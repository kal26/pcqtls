#!/usr/bin/env bash

#conda activate snakemake
snakemake --configfile config/pcqtl.yaml --cores 10 --use-conda --profile scg --resources -j 100 --cluster-cancel scancel