import pandas as pd
import subprocess
import os


# wrapper script is needed to make blank files for the chrs that don't end up with results

genotype_stem = snakemake.params[0]
expression = snakemake.input[1]
eqtl_output_dir = snakemake.params[1]
tissue_id = snakemake.params[2]
covariates = snakemake.input[2]



# run the shell command
nominal_comand = f"python -m tensorqtl {genotype_stem} {expression} {eqtl_output_dir}{tissue_id}/{tissue_id}.v8.cluster_genes --covariates {covariates} --mode cis_nominal --maf_threshold .01"
try:
    # Run the command
    result = subprocess.run(nominal_comand, shell=True, check=True, text=True, capture_output=True)
    # Print the output
    print("Command Output:\n", result.stdout)
except subprocess.CalledProcessError as e:
    print("Command failed with error:\n", e.stderr)


# check for missing parquet files
chr_range = range(1, 23)

# Define the base filename
base_filename = f'{eqtl_output_dir}{tissue_id}/{tissue_id}.v8.cluster_genes.cis_qtl_pairs.chr'

# Check for existing files and create missing ones
for chr_id in chr_range:
    parquet_filename = f"{base_filename}{i}.parquet"
    # Check if the output file exists
    if not os.path.exists(parquet_filename):
        # Create a blank DataFrame with the specified columns
        blank = pd.DataFrame(columns=["phenotype_id", "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df", "variant_id", "start_distance", "end_distance", "ma_samples", "ma_count", "af", "pval_nominal", "slope", "slope_se", "pval_perm", "pval_beta", "qval", "pval_nominal_threshold"])
        # Save the empty DataFrame as a Parquet file
        blank.to_parquet(parquet_filename)
        print(f"Created missing file: {parquet_filename}")
    else:
        print(f"File exists: {parquet_filename}")
