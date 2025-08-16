#!/usr/bin/env Rscript
# =============================================================================
# Script: combine_RDS_susie.R
# Description: Combine SuSiE RDS files from colocalization analysis
# Author: Kate Lawrence
# Date: 2025-08-15
# =============================================================================

# Load required libraries
library(dplyr)
library(argparse)

# =============================================================================
# Argument parsing
# =============================================================================

# Create argument parser
parser <- ArgumentParser()
parser$add_argument("--code-dir", help="Directory containing R script functions", required=TRUE)
parser$add_argument("--qtl-dir", help="Folder for QTL pairs")
parser$add_argument("--tissue-id", help="Tissue ID")
parser$add_argument("--coloc-temp-path-head", help="Tissue specific directory file path for susie and per chrom GWAS")
parser$add_argument("--output", help="Output file path for coloc pair results")
parser$add_argument("--qtl-type", help="eQTL or PCQTL")

# Parse arguments
args <- parser$parse_args()

# Extract arguments
code_dir <- args$code_dir
coloc_temp_path_head <- args$coloc_temp_path_head
qtl_dir_path <- args$qtl_dir
tissue_id <- args$tissue_id
qtl_type <- args$qtl_type
output_path <- args$output

# =============================================================================
# Source helper functions
# =============================================================================

# Source helper functions from the provided code directory
source(file.path(code_dir, "coloc_functions.R"))

# =============================================================================
# Combine SuSiE results across chromosomes
# =============================================================================

# Initialize results dataframe
combined_results <- data.frame()

# Process each chromosome
for (chr_id in 1:22) {
  cat("Working on chromosome", chr_id, "\n")
  
  # Load nominal QTL results to get phenotype IDs
  if (qtl_type == 'eqtl') {
    qtl_chr <- get_eqtl_chr(qtl_dir_path, chr_id, tissue_id)
  } else {
    qtl_chr <- get_pcqtl_chr(qtl_dir_path, chr_id, tissue_id)
  }
  
  phenotype_ids <- unique(qtl_chr$phenotype_id)
  
  # Process each phenotype
  for (phenotype_id in phenotype_ids) {
    # Handle long QTL IDs
    if (nchar(phenotype_id) > 150) {
      out_qtl_id <- get_short_qtl_id(phenotype_id)
    } else {
      out_qtl_id <- phenotype_id
    }
    
    # Generate path to SuSiE RDS file
    susie_rds_path <- paste(coloc_temp_path_head, out_qtl_id, '.susie.rds', sep = "")
    
    # Check if file exists and extract data
    if (file.exists(susie_rds_path)) {
      # Extract data from RDS file and combine
      rds_data <- extract_data_from_rds(susie_rds_path)
      combined_results <- bind_rows(combined_results, rds_data)
    } else {
      # Expected that some files do not exist (no significant QTL signals)
      message(paste("File does not exist:", susie_rds_path))
    }
  }
}

# =============================================================================
# Write combined results
# =============================================================================

write.table(combined_results, file = output_path, sep = "\t", quote = FALSE, row.names = TRUE)
cat("Combined SuSiE results written to:", output_path, "\n")
