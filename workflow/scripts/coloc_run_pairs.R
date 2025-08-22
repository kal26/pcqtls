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
parser$add_argument("--eqtl-nominal-dir", help="Path for eQTL pairs")
parser$add_argument("--pcqtl-nominal-dir", help="Path for PCQTL pairs")
parser$add_argument("--gtex-meta", help="Input file path for GTEX metadata with sample size")
parser$add_argument("--tissue-id", help="Tissue ID")
parser$add_argument("--chr-id", help="Chromosome ID")
parser$add_argument("--ld-path-head", help="Directory file path for SNP list and LD")
parser$add_argument("--genotype-stem", help="Path to genotype stem")
parser$add_argument("--clusters", help="Path to clusters file")
parser$add_argument("--output", help="Output file path for coloc results")
parser$add_argument("--coloc-temp-path-head", help="Tissue specific directory file path for susie and per chrom GWAS")

# Parse arguments
args <- parser$parse_args()

# Extract arguments
code_dir <- args$code_dir
eqtl_nominal_dir <- args$eqtl_nominal_dir
pcqtl_nominal_dir <- args$pcqtl_nominal_dir
gtex_meta_path <- args$gtex_meta
tissue_id <- args$tissue_id
chr_id <- args$chr_id
ld_path_head <- args$ld_path_head
genotype_stem <- args$genotype_stem
output_path <- args$output
clusters_path <- args$clusters
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
cat("eQTL:", eqtl_nominal_dir, "\n")
cat("PCQTL:", pcqtl_nominal_dir, "\n")
cat("GTEX Metadata:", gtex_meta_path, "\n")
cat("LD and SNP list dir:", ld_path_head, "\n")
cat("Tissue ID:", tissue_id, "\n")
cat("Chromosome ID:", chr_id, "\n")
cat("Output Path:", output_path, "\n")
cat("Cluster Path:", clusters_path, "\n")

cat("Starting pair colocalizations for", tissue_id, "chromosome", chr_id, "\n")
start <- Sys.time()

# =============================================================================
# Load GTEX metadata for sample size
# =============================================================================

cat("Loading GTEX metadata for sample size...\n")
gtex_meta <- read.table(gtex_meta_path, sep='\t', header = TRUE)
num_gtex_samples <- gtex_meta[gtex_meta$tissue_id == tissue_id, 'sample_size']
cat("GTEX sample size for", tissue_id, ":", num_gtex_samples, "\n")

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
cat("Setting up output directories...\n")
if (!file.exists(ld_path_head)) {
  dir.create(ld_path_head, recursive = TRUE)
  cat("LD directory created:", ld_path_head, "\n")
} else {
  cat("LD directory already exists:", ld_path_head, "\n")
}

# Create coloc temp directory if it doesn't exist
if (!file.exists(coloc_temp_path_head)) {
  dir.create(coloc_temp_path_head, recursive = TRUE)
  cat("Coloc temp directory created:", coloc_temp_path_head, "\n")
} else {
  cat("Coloc temp directory already exists:", coloc_temp_path_head, "\n")
}

# =============================================================================
# Load clusters and process chromosome
# =============================================================================

# Load clusters
cat("Loading cluster information...\n")
cluster_df <- fread(clusters_path)
cat("Loaded", nrow(cluster_df), "total clusters\n")

# Process specified chromosome
chr_id <- as.integer(sub("chr", "", chr_id))
chr_coloc_path <- paste(coloc_temp_path_head, tissue_id, ".v8.", 'chr_', chr_id, '.pairs_coloc.txt', sep="")

cat("Working on chromosome", chr_id, "\n")

# Subset clusters to current chromosome
cluster_df_chr <- cluster_df[cluster_df$chr == paste0("chr", chr_id)]
cat("Found", nrow(cluster_df_chr), "clusters on chromosome", chr_id, "\n")

# Load eQTL and PCQTL data for current chromosome
cat("Loading QTL data for chromosome", chr_id, "...\n")
eqtl_chr <- get_eqtl_chr(eqtl_nominal_dir, chr_id, tissue_id)
pcqtl_chr <- get_pcqtl_chr(pcqtl_nominal_dir, chr_id, tissue_id)
cat("Loaded", nrow(eqtl_chr), "eQTLs and", nrow(pcqtl_chr), "PCQTLs for chromosome", chr_id, "\n")

# =============================================================================
# Run colocalization for each cluster
# =============================================================================

cat("Starting colocalization analysis for", nrow(cluster_df_chr), "clusters...\n")
for (i in seq_len(nrow(cluster_df_chr))) {
  this_cluster <- cluster_df_chr[i]
  cluster_id <- this_cluster$cluster_id
  cat("Processing cluster", i, "of", nrow(cluster_df_chr), ":", cluster_id, "\n")
  
  num_colocs <- num_colocs + 1
  cluster_coloc <- coloc_pairs_cluster(eqtl_chr, pcqtl_chr, cluster_id, ld_path_head, 
                                      genotype_stem, num_gtex_samples, coloc_temp_path_head)
  
  # Add results to main dataframe
  cat("Adding to results list\n")
  cat(ncol(all_cluster_coloc_results), "vs", ncol(cluster_coloc), "\n")
  
  if (is.null(cluster_coloc)) {
    cat("\tResult is null\n")
  } else {
    cat("\tResult is not null\n")
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
# Reorder columns and write results
# =============================================================================

# Reorder columns to put qtl1_id and qtl2_id first
if (nrow(all_cluster_coloc_results) > 0) {
  all_cluster_coloc_results <- all_cluster_coloc_results[, c("qtl1_id", "qtl2_id", 
                                                             setdiff(colnames(all_cluster_coloc_results), c("qtl1_id", "qtl2_id")))]
}

write.table(all_cluster_coloc_results, file = output_path, quote = FALSE, row.names = FALSE, sep = '\t')
cat("Colocalization analysis completed. Results written to:", output_path, "\n")