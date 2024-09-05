import os
import pandas as pd

# Directory path containing the CSV files
coloc_paths = snakemake.input[0]

# List to store individual DataFrames
cluster_colocs = []

# Iterate over files in the directory
for coloc_path in coloc_paths:
    this_cluster_coloc = pd.read_csv(coloc_paths)
    cluster_colocs.append(this_cluster_coloc)

# Concatenate all DataFrames into one
aggregated_df = pd.concat(cluster_colocs, ignore_index=True)

# Display the aggregated DataFrame
cluster_colocs.to_csv(snakemake.output[0])