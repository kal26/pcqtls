import os
import pandas as pd

# Directory path containing the CSV files
coloc_temp_directory = snakemake.params[0]

# List to store individual DataFrames
cluster_colocs = []

# Iterate over files in the directory
for file in os.listdir(coloc_temp_directory):
    file_path = os.path.join(directory, file)
    this_cluster_coloc = pd.read_csv(file_path)
    cluster_colocs.append(this_cluster_coloc)

# Concatenate all DataFrames into one
aggregated_df = pd.concat(cluster_colocs, ignore_index=True)

# Display the aggregated DataFrame
cluster_colocs.to_csv(snakemake.output[0])