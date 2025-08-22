#!/usr/bin/env Rscript
# =============================================================================
# Script: coloc_run_gwas.R
# Description: Run colocalization analysis between GWAS and QTLs
# Author: Kate Lawrence
# Date: 2025-08-15
# 
# Based on tutorial at: https://hanruizhang.github.io/GWAS-eQTL-Colocalization/
# Requires R/4.2.2
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
parser$add_argument("--eqtl-nominal-dir", help="Path for eQTL nominal pairs")
parser$add_argument("--pcqtl-nominal-dir", help="Path for PCQTL nominal pairs")
parser$add_argument("--gwas", help="Path for GWAS summary statistics")
parser$add_argument("--gwas-id", help="GWAS ID")
parser$add_argument("--gwas-meta", help="Input file path for GWAS metadata")
parser$add_argument("--gtex-meta", help="Input file path for GTEX metadata with sample size")
parser$add_argument("--tissue-id", help="Tissue ID")
parser$add_argument("--ld-path-head", help="Directory file path for SNP list and LD")
parser$add_argument("--coloc-temp-path-head", help="Tissue specific directory file path for susie and per chrom GWAS")
parser$add_argument("--genotype-stem", help="Path to genotype stem")
parser$add_argument("--clusters", help="Path to clusters file")
parser$add_argument("--expression", help="Path to cluster expression file with coordinates")
parser$add_argument("--output", help="Output file path for coloc results")
parser$add_argument("--use-susie", help="Whether to use SuSiE for fine-mapping", default="True")

# Parse arguments
args <- parser$parse_args()

# Extract arguments
code_dir <- args$code_dir
eqtl_dir_path <- args$eqtl_nominal_dir
pcqtl_dir_path <- args$pcqtl_nominal_dir
gwas_path <- args$gwas
gwas_id <- args$gwas_id
gwas_meta_path <- args$gwas_meta
gtex_meta_path <- args$gtex_meta
tissue_id <- args$tissue_id
ld_path_head <- args$ld_path_head
coloc_temp_path_head <- args$coloc_temp_path_head
genotype_stem <- args$genotype_stem
output_path <- args$output
clusters_path <- args$clusters
expression_path <- args$expression
use_susie <- as.logical(args$use_susie)

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
cat("GWAS Metadata:", gwas_meta_path, "\n")
cat("GTEX Metadata:", gtex_meta_path, "\n")
cat("LD and SNP list dir:", ld_path_head, "\n")
cat("Tissue ID:", tissue_id, "\n")
cat("GWAS ID:", gwas_id, "\n")
cat("GWAS path:", gwas_path, "\n")
cat("Output Path:", output_path, "\n")
cat("Cluster Path:", clusters_path, "\n")
cat("Expression Path:", expression_path, "\n")
cat("Use SuSiE:", use_susie, "\n")

cat("Starting colocalizations for", tissue_id, "and GWAS ID:", gwas_id, "\n")
start <- Sys.time()

# =============================================================================
# Load metadata and GWAS data
# =============================================================================

# Load GTEX metadata for sample size
cat("Loading GTEX metadata for sample size...\n")
gtex_meta <- read.table(gtex_meta_path, sep = '\t', header = TRUE)
num_gtex_samples <- gtex_meta[gtex_meta$tissue_id == tissue_id, 'sample_size']
cat("GTEX sample size for", tissue_id, ":", num_gtex_samples, "\n")

# Load GWAS metadata and data
cat("Loading GWAS metadata and data...\n")
gwas_meta_df <- fread(gwas_meta_path)
gwas_with_meta <- load_gwas_from_path(gwas_path, gwas_meta_df, gwas_id)
gwas_data <- gwas_with_meta$gwas_data
cat("Total GWAS signals:", sum(gwas_data$pvalue < 1e-6), "\n")
cat("GWAS data loaded successfully\n")

# =============================================================================
# Initialize results
# =============================================================================

# Initialize empty results dataframe
gwas_all_cluster_coloc_results <- data.frame()
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
# Load clusters and process chromosomes
# =============================================================================

# Load clusters
cat("Loading cluster information...\n")
cluster_df <- fread(clusters_path)
cat("Loaded", nrow(cluster_df), "clusters\n")

# Load expression file to get coordinates
cat("Loading expression file for coordinates...\n")
expression_df <- fread(expression_path)
cat("Loaded expression data for", nrow(expression_df), "genes\n")

# Add cluster_id column to expression dataframe
cat("Processing expression data to extract cluster IDs...\n")
expression_df$cluster_id <- sapply(expression_df$gene_id, function(x) {
  parts <- strsplit(as.character(x), "_e_")[[1]]
  if (length(parts) > 0) return(parts[1]) else return(NA)
})
cat("Extracted cluster IDs from", sum(!is.na(expression_df$cluster_id)), "expression rows\n")

# Add start and end coordinates to clusters
cat("Adding genomic coordinates to clusters...\n")
cluster_df$start <- NA
cluster_df$end <- NA
cluster_df$chr_num <- as.integer(sub("chr", "", cluster_df$chr))


for (i in seq_len(nrow(cluster_df))) {
  this_cluster_id <- cluster_df$cluster_id[i]
  # Find matching rows in expression data
  cluster_expression <- expression_df %>% filter(cluster_id == this_cluster_id)
  if (nrow(cluster_expression) > 0) {
    # Use coordinates from the first matching row
    cluster_df$start[i] <- cluster_expression$start[1]
    cluster_df$end[i] <- cluster_expression$end[1]
  }
}

cat("Genomic coordinates added to clusters\n")
cat("Clusters with coordinates:", sum(!is.na(cluster_df$start) & !is.na(cluster_df$end)), "of", nrow(cluster_df), "\n")
missing_coords <- which(is.na(cluster_df$start) | is.na(cluster_df$end))
if (length(missing_coords) > 0) {
  cat("Example clusters missing coordinates (up to 5):", paste(head(cluster_df$cluster_id[missing_coords], 5), collapse=", "), "\n")
}

# Process each chromosome
cat("Starting chromosome-by-chromosome processing...\n")
for (chr_id in 1:22) {
  chr_coloc_path <- paste(coloc_temp_path_head, tissue_id, ".v8.", tissue_id, '.susie_', use_susie, 'chr_', chr_id, '.gwas_coloc.txt', sep = "")
  
  # Check if partial results exist
  if (file.exists(chr_coloc_path)) {
    # Check if file is empty or only contains whitespace
    file_size <- file.size(chr_coloc_path)
    if (file_size <= 1) {
      cat("Colocalization results file exists but is empty for chromosome", chr_id, "- regenerating\n")
      # Remove the empty file and regenerate
      unlink(chr_coloc_path)
    } else {
      cat("Colocalization results already exist for chromosome", chr_id, "- loading existing results\n")
      gwas_all_cluster_coloc_results <- read.table(chr_coloc_path, header = TRUE, sep = '\t')
    }
  }
  
  # If file doesn't exist or was empty, process the chromosome
  if (!file.exists(chr_coloc_path)) {
    cat("Time elapsed:", Sys.time() - start, "\n")
    cat("Processing chromosome", chr_id, "of 22\n")
    
    # Subset clusters to current chromosome
    cluster_df_chr <- cluster_df[cluster_df$chr_num == chr_id]
    
    # Subset GWAS data to current chromosome
    gwas_chr <- gwas_with_meta$gwas_data[gwas_with_meta$gwas_data$chromosome == paste('chr', chr_id, sep = "")] 
    
    # Initialize QTL data as NULL (will be loaded only if needed)
    pcqtl_chr <- NULL
    eqtl_chr <- NULL
    
    # =============================================================================
    # Process each cluster
    # =============================================================================
    
    for (i in 1:nrow(cluster_df_chr)) {
      this_cluster <- cluster_df_chr[i]
      cluster_id <- this_cluster$cluster_id
      cat("Running on", cluster_id, "\n")
      
      # Check if GWAS has a signal in this cluster
      if (check_gwas_cluster(gwas_chr, this_cluster)) {
        num_colocs <- num_colocs + 1
        cat("\tPossible colocalization for", this_cluster$cluster_id, "\n")
        cat(num_colocs, "colocalizations so far\n")
        
        # Load eQTL data if not already loaded
        if (is.null(eqtl_chr)) {
          eqtl_chr <- get_eqtl_chr(eqtl_dir_path, chr_id, tissue_id)
        }
        
        # Load PCQTL data if not already loaded
        if (is.null(pcqtl_chr)) {
          pcqtl_chr <- get_pcqtl_chr(pcqtl_dir_path, chr_id, tissue_id)
        }
        
        # Run colocalization for this cluster
        gwas_cluster_coloc <- coloc_gwas_cluster(gwas_with_meta, eqtl_chr, pcqtl_chr, cluster_id, 
                                                ld_path_head, coloc_temp_path_head, genotype_stem, 
                                                num_gtex_samples, use_susie = use_susie)
        
        cat("Finished cluster colocalization\n")
        
        # Ensure data types match for binding
        gwas_all_cluster_coloc_results$gwas_id <- as.character(gwas_all_cluster_coloc_results$gwas_id)
        gwas_cluster_coloc$gwas_id <- as.character(gwas_cluster_coloc$gwas_id)
        gwas_all_cluster_coloc_results$qtl_id <- as.character(gwas_all_cluster_coloc_results$qtl_id)
        gwas_cluster_coloc$qtl_id <- as.character(gwas_cluster_coloc$qtl_id)
        
        # Add results to main dataframe
        gwas_all_cluster_coloc_results <- bind_rows(gwas_all_cluster_coloc_results, gwas_cluster_coloc) 
      }
    }
    
    # Write partial results for this chromosome
    write.table(gwas_all_cluster_coloc_results, file = chr_coloc_path, quote = FALSE, row.names = FALSE, sep = '\t')
    cat("Time elapsed:", Sys.time() - start, "\n")
    cat("Wrote out partial results, up to chromosome", chr_id, "\n")
  }
}

# =============================================================================
# Write final results
# =============================================================================

write.table(gwas_all_cluster_coloc_results, file = output_path, quote = FALSE, row.names = FALSE, sep = '\t')
cat("Finished. Wrote out all results to:", output_path, "\n")
