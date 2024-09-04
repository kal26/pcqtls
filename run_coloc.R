
library(argparse)
library(coloc)
library(arrow)
library(tidyverse)

# Create Argument Parser
parser <- ArgumentParser()
parser$add_argument("--eqtl_path", help="Input file path for eQTL pairs")
parser$add_argument("--pcqtl_path", help="Input file path for PCQTL pairs")
parser$add_argument("--gwas_meta", help="Input file path for GWAS metadata")
parser$add_argument("--ld_path", help="Input file path for LD matrix")
parser$add_argument("--cluster", help="cluster id")
parser$add_argument("--gwas_folder", help="Folder path for GWAS data")
parser$add_argument("--output_path", help="Output file path for coloc results")

# Parse the Arguments
args <- parser$parse_args()

# Access the arguments
eqtl_path <- args$eqtl_path
pcqtl_path <- args$pcqtl_path
gwas_meta <- args$gwas_meta
ld_path <- args$ld_path
cluster_id <- args$cluster
gwas_folder <- args$gwas_folder
output_path <- args$output_path

# Print or use the arguments as needed in the script
cat("eQTL Path:", eqtl_path, "\n")
cat("PCQTL Path:", pcqtl_path, "\n")
cat("GWAS Metadata:", gwas_meta, "\n")
cat("LD Path:", ld_path, "\n")
cat("Cluster ID:", cluster_id, "\n")
cat("GWAS Folder:", gwas_folder, "\n")
cat("Output Path:", output_path, "\n")


# load in eqtl dataset
# load in pcqtl dataset
# load in gwas dataset 
# function to get qtl dataset 

# function to finemap qtl dataset

# function to get gwas dataset

# funciton to finemap gwas dataset

# funciton to loop through gwas and coloc with qtl

# function to coloc each eqtl and pcqtl pair

# function to aggregate results



# load in eqtl data

# load in pcqtl data
# load in ld data
# finemap eqtls for cluster
# finemap pcqtls for cluster
# coloc eqtls and pcqtls
# loop through gwas and coloc with e/pcqtl
# aggregate
# write out 