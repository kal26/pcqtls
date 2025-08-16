#!/usr/bin/env Rscript
# =============================================================================
# Script: coloc_run_pairs.R
# Description: Run colocalization analysis between eQTL and PCQTL pairs
# Author: Kate Lawrence
# Date: 2025-08-15
# =============================================================================

# Load required libraries
library(coloc)
library(nanoparquet)
library(tidyverse)
library(data.table)
library(Rfast)
library(argparse)
library(susieR)
library(dplyr)

# =============================================================================
# Argument parsing
# =============================================================================

# Create argument parser
parser <- ArgumentParser()
parser$add_argument("--code-dir", help="Directory containing R script functions", required=TRUE)
parser$add_argument("--eqtl-dir", help="Path for eQTL pairs")
parser$add_argument("--pcqtl-dir", help="Path for PCQTL pairs")
parser$add_argument("--gtex-meta", help="Input file path for GTEX metadata with sample size")
parser$add_argument("--tissue-id", help="Tissue ID")
parser$add_argument("--chr-id", help="Chromosome ID")
parser$add_argument("--ld-path-head", help="Directory file path for SNP list and LD")
parser$add_argument("--genotype-stem", help="Path to genotype stem")
parser$add_argument("--annotated-cluster", help="Path to position annotated clusters")
parser$add_argument("--output", help="Output file path for coloc results")
parser$add_argument("--coloc-temp-path-head", help="Tissue specific directory file path for susie and per chrom GWAS")

# Parse arguments
args <- parser$parse_args()

# Extract arguments
code_dir <- args$code_dir
eqtl_dir_path <- args$eqtl_dir
pcqtl_dir_path <- args$pcqtl_dir
gtex_meta_path <- args$gtex_meta
tissue_id <- args$tissue_id
chr_id <- args$chr_id
ld_path_head <- args$ld_path_head
genotype_stem <- args$genotype_stem
output_path <- args$output
annotated_cluster_path <- args$annotated_cluster
coloc_temp_path_head <- args$coloc_temp_path_head

# =============================================================================
# Source helper functions
# =============================================================================

# Source helper functions from the provided code directory
source(file.path(code_dir, "coloc_functions.R"))

# =============================================================================
# Print configuration
# =============================================================================

cat("Code directory:", code_dir, "\n")
cat("eQTL folder:", eqtl_dir_path, "\n")
cat("PCQTL folder:", pcqtl_dir_path, "\n")
cat("GTEX Metadata:", gtex_meta_path, "\n")
cat("LD and SNP list dir:", ld_path_head, "\n")
cat("Tissue ID:", tissue_id, "\n")
cat("Chromosome ID:", chr_id, "\n")
cat("Output Path:", output_path, "\n")
cat("Cluster Path:", annotated_cluster_path, "\n")

cat("Starting pair colocalizations\n")
start <- Sys.time()

# =============================================================================
# Load GTEX metadata for sample size
# =============================================================================

gtex_meta <- read.table(gtex_meta_path, sep='\t', header = TRUE)
num_gtex_samples <- gtex_meta[gtex_meta$tissue_id == tissue_id, 'sample_size']

# =============================================================================
# Initialize results
# =============================================================================

# Initialize empty results dataframe
all_cluster_coloc_results <- data.frame()
num_colocs <- 0

# =============================================================================
# Create output directories
# =============================================================================

# Create LD directory if it doesn't exist
if (!file.exists(ld_path_head)) {
  dir.create(ld_path_head, recursive = TRUE)
  cat("Directory created:", ld_path_head, "\n")
} else {
  cat("Directory already exists:", ld_path_head, "\n")
}

# Create coloc temp directory if it doesn't exist
if (!file.exists(coloc_temp_path_head)) {
  dir.create(coloc_temp_path_head, recursive = TRUE)
  cat("Directory created:", coloc_temp_path_head, "\n")
} else {
  cat("Directory already exists:", coloc_temp_path_head, "\n")
}

# =============================================================================
# Load clusters and process chromosome
# =============================================================================

# Load annotated clusters
cluster_df <- fread(annotated_cluster_path)

# Process specified chromosome
chr_id <- as.integer(sub("chr", "", chr_id))
chr_coloc_path <- paste(coloc_temp_path_head, tissue_id, ".v8.", 'chr_', chr_id, '.pairs_coloc.txt', sep="")

cat("Working on chromosome", chr_id, "\n")

# Subset clusters to current chromosome
cluster_df_chr <- cluster_df[cluster_df$Chromosome == chr_id]

# Load eQTL and PCQTL data for current chromosome
eqtl_chr <- get_eqtl_chr(eqtl_dir_path, chr_id, tissue_id)
pcqtl_chr <- get_pcqtl_chr(pcqtl_dir_path, chr_id, tissue_id)

# =============================================================================
# Run colocalization for each cluster
# =============================================================================

for (i in 1:nrow(cluster_df_chr)) {
  this_cluster <- cluster_df_chr[i]
  cluster_id <- this_cluster$cluster_id
  cat("Running colocalization on", cluster_id, "\n")
  
  num_colocs <- num_colocs + 1
  cluster_coloc <- coloc_pairs_cluster(eqtl_chr, pcqtl_chr, cluster_id, ld_path_head, 
                                      genotype_stem, num_gtex_samples, coloc_temp_path_head)
  
  # Add results to main dataframe
  cat("Adding to results list\n")
  cat(ncol(all_cluster_coloc_results), "vs", ncol(cluster_coloc), "\n")
  
  if (is.null(cluster_coloc)) {
    cat("\tResult is null\n")
  } else {
    all_cluster_coloc_results <- bind_rows(all_cluster_coloc_results, cluster_coloc) 
    cat(num_colocs, "colocalizations completed so far\n")
  }
}

# =============================================================================
# Add required columns for signal groups analysis
# =============================================================================

if (nrow(all_cluster_coloc_results) > 0) {
  all_cluster_coloc_results$cs_id_1 <- paste0(all_cluster_coloc_results$qtl1_id, "_cs_", all_cluster_coloc_results$idx1)
  all_cluster_coloc_results$cs_id_2 <- paste0(all_cluster_coloc_results$qtl2_id, "_cs_", all_cluster_coloc_results$idx2)
  all_cluster_coloc_results$cluster_id <- ifelse(grepl("_pc", all_cluster_coloc_results$qtl1_id), 
                                                 gsub("_pc.*", "", all_cluster_coloc_results$qtl1_id), 
                                                 gsub("_e_.*", "", all_cluster_coloc_results$qtl1_id))
}

# =============================================================================
# Write results
# =============================================================================

write.table(all_cluster_coloc_results, file = output_path, quote = FALSE, row.names = FALSE, sep = '\t')
cat("Colocalization analysis completed. Results written to:", output_path, "\n")