#!/usr/bin/env Rscript
# =============================================================================
# Script: coloc_functions.R
# Description: Helper functions for colocalization analysis between QTLs and GWAS
# Author: Kate Lawrence
# Date: 2025-08-15
# =============================================================================

# Set working directory
setwd('/home/klawren/oak/pcqtls/')

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
# QTL processing functions
# =============================================================================

#' Split QTL data for colocalization analysis
#' 
#' @param cluster_qtl QTL data for a cluster
#' @param ld_snp_set Set of SNPs in LD
#' @param cleaned_ld_matrix Cleaned LD matrix
#' @param num_gtex_samples Number of GTEX samples
#' @return List of QTL data formatted for colocalization
split_qtl_for_coloc <- function(cluster_qtl, ld_snp_set, cleaned_ld_matrix, num_gtex_samples) {
  qtl_filtered <- filter_qtl(cluster_qtl, ld_snp_set)
  
  # Split eQTL to each of the eGenes and clean up
  phenotype_ids <- unique(qtl_filtered$phenotype_id)
  qtls_for_coloc <- list()
  
  for (i in 1:length(phenotype_ids)) {
    phenotype_qtl <- qtl_filtered[qtl_filtered$phenotype_id == phenotype_ids[i], ]
    cat(paste("\tLooking for signals in", phenotype_ids[i]), "\n")
    cat(paste("\t\tQTL filtered to:", nrow(phenotype_qtl), "SNPs\n"))
    
    phenotype_qtl_for_coloc <- clean_qtl(phenotype_qtl, cleaned_ld_matrix, num_gtex_samples)
    qtls_for_coloc[[i]] <- phenotype_qtl_for_coloc
  }
  
  return(qtls_for_coloc)
}

#' Filter QTL data to include only SNPs in LD set
#' 
#' @param qtl QTL data
#' @param ld_snp_set Set of SNPs in LD
#' @return Filtered QTL data
filter_qtl <- function(qtl, ld_snp_set) {
  qtl_filtered <- qtl[qtl$variant_id %in% ld_snp_set, ]
  return(qtl_filtered)
}

#' Clean QTL data for colocalization analysis
#' 
#' @param qtl_filtered Filtered QTL data
#' @param cleaned_ld_matrix Cleaned LD matrix
#' @param num_gtex_samples Number of GTEX samples
#' @return Cleaned QTL data formatted for colocalization
clean_qtl <- function(qtl_filtered, cleaned_ld_matrix, num_gtex_samples) {
  # Remove SNPs with missing slope, slope_se, or zero allele frequency
  qtl_missing_snps <- qtl_filtered[is.na(qtl_filtered$slope) | 
                                   is.na(qtl_filtered$slope_se) | 
                                   (qtl_filtered$af == 0), 'variant_id']
  cat("\t\tNumber of SNPs with missing QTL data (removing):", length(qtl_missing_snps), "\n")
  
  cleaned_qtl <- qtl_filtered[!qtl_filtered$variant_id %in% qtl_missing_snps, ]
  
  # Check if there is a significant signal
  if (min(cleaned_qtl$pval_nominal, na.rm = TRUE) < 1e-6) {
    cat('\t\tSignal found\n')
    
    # Create position column for colocalization
    cleaned_qtl <- cleaned_qtl %>%
      mutate(position = as.integer(str_split(variant_id, "_") %>% sapply(pluck, 2)))
    
    # Subset LD matrix to include only SNPs in cleaned QTL data
    qtl_cleaned_ld_matrix <- cleaned_ld_matrix[rownames(cleaned_ld_matrix) %in% cleaned_qtl$variant_id, 
                                              colnames(cleaned_ld_matrix) %in% cleaned_qtl$variant_id]
    
    # Format data for colocalization
    cleaned_qtl_list <- list(beta = cleaned_qtl$slope, 
                            varbeta = cleaned_qtl$slope_se**2, 
                            snp = cleaned_qtl$variant_id, 
                            position = cleaned_qtl$position, 
                            pvalues = cleaned_qtl$pval_nominal,
                            type = "quant", 
                            N = num_gtex_samples, 
                            MAF = cleaned_qtl$af,
                            sdY = 1,
                            phenotype_id = cleaned_qtl$phenotype_id[1],
                            LD = as.matrix(qtl_cleaned_ld_matrix))
    return(cleaned_qtl_list)
  } else {
    cat('\t\tNo signal found\n')
    return(NULL)
  }    
}

# =============================================================================
# GWAS processing functions
# =============================================================================

#' Get GWAS data formatted for colocalization
#' 
#' @param gwas_with_meta GWAS data with metadata
#' @param ld_snp_set Set of SNPs in LD
#' @param snp_list List of SNPs
#' @param cleaned_ld_matrix Cleaned LD matrix
#' @return GWAS data formatted for colocalization
get_gwas_for_coloc <- function(gwas_with_meta, ld_snp_set, snp_list, cleaned_ld_matrix) {
  gwas_filtered <- filter_gwas(gwas_with_meta$gwas_data, snp_list, ld_snp_set)
  gwas_for_coloc <- clean_gwas(gwas_filtered, cleaned_ld_matrix, gwas_with_meta$gwas_type, gwas_with_meta$gwas_sample_size)
  return(gwas_for_coloc)
}

#' Filter GWAS data to include only SNPs in SNP list and LD set
#' 
#' @param gwas_data GWAS data
#' @param snp_list List of SNPs
#' @param ld_snp_set Set of SNPs in LD
#' @return Filtered GWAS data
filter_gwas <- function(gwas_data, snp_list, ld_snp_set) {
  gwas_filtered <- gwas_data[gwas_data$panel_variant_id %in% snp_list$variant_id, ]
  gwas_filtered <- gwas_filtered[gwas_filtered$panel_variant_id %in% ld_snp_set, ]
  cat("SNPs in GWAS (should match SNPs in QTL):", nrow(gwas_filtered), "\n")
  return(gwas_filtered)
}

#' Clean GWAS data for colocalization analysis
#' 
#' @param gwas_filtered Filtered GWAS data
#' @param cleaned_ld_matrix Cleaned LD matrix
#' @param gwas_type Type of GWAS (quantitative or case-control)
#' @param num_gwas_samples Number of GWAS samples
#' @return Cleaned GWAS data formatted for colocalization
clean_gwas <- function(gwas_filtered, cleaned_ld_matrix, gwas_type, num_gwas_samples) {
  # Remove SNPs with missing effect sizes, standard errors, or invalid frequencies
  gwas_missing_snps <- gwas_filtered[is.na(gwas_filtered$effect_size) | 
                                     is.na(gwas_filtered$standard_error) | 
                                     (gwas_filtered$frequency <= 0) | 
                                     (gwas_filtered$frequency >= 1), 'panel_variant_id']
  cat("\tNumber of SNPs with missing GWAS data (removing):", length(gwas_missing_snps$panel_variant_id), "\n")
  
  cleaned_gwas <- gwas_filtered[!gwas_filtered$panel_variant_id %in% gwas_missing_snps$panel_variant_id, ]
  
  if (min(cleaned_gwas$pvalue, na.rm = TRUE) < 1e-6) {
    gwas_cleaned_ld_matrix <- cleaned_ld_matrix[rownames(cleaned_ld_matrix) %in% cleaned_gwas$panel_variant_id, 
                                               colnames(cleaned_ld_matrix) %in% cleaned_gwas$panel_variant_id]
    
    # Format data for colocalization
    cleaned_gwas_list <- list(MAF = cleaned_gwas$frequency, 
                             snp = cleaned_gwas$panel_variant_id, 
                             position = cleaned_gwas$position, 
                             type = gwas_type, 
                             N = num_gwas_samples, 
                             pvalues = cleaned_gwas$pvalue,
                             LD = as.matrix(gwas_cleaned_ld_matrix), 
                             beta = cleaned_gwas$effect_size,
                             varbeta = cleaned_gwas$standard_error**2)
    cat('\tSignal found\n')
    return(cleaned_gwas_list)
  } else {
    cat('\tNo signal found, this is odd because I already checked and there should be a signal\n')
    return(NULL)
  }
}

# =============================================================================
# SuSiE functions
# =============================================================================

#' Run SuSiE with error handling
#' 
#' @param dataset Dataset for SuSiE analysis
#' @return SuSiE result or NULL if error occurs
runsusie_errorcatch <- function(dataset) {
  start <- Sys.time()
  susie <- tryCatch({
    runsusie(dataset, repeat_until_convergence = TRUE)
  }, error = function(e) {
    cat("An error occurred for QTL\n")
    cat("Error Message:", conditionMessage(e), "\n")
    NULL
  })
  cat("SuSiE runtime:", Sys.time() - start, "\n")
  return(susie)
}

# =============================================================================
# Colocalization functions
# =============================================================================

#' Run colocalization analysis for QTL pairs within a cluster
#' 
#' @param eqtl_chr eQTL data for chromosome
#' @param pcqtl_chr PCQTL data for chromosome
#' @param cluster_id Cluster ID
#' @param ld_path_head Path to LD files
#' @param genotype_stem Path to genotype files
#' @param num_gtex_samples Number of GTEX samples
#' @param coloc_temp_path_head Path to temporary colocalization files
#' @return Colocalization results for QTL pairs
coloc_pairs_cluster <- function(eqtl_chr, pcqtl_chr, cluster_id, ld_path_head, genotype_stem, 
                               num_gtex_samples, coloc_temp_path_head) {
  start <- Sys.time()
  
  # Subset eQTL and PCQTL to this cluster
  cluster_eqtl <- eqtl_chr[eqtl_chr$cluster_id == cluster_id, ]
  cluster_pcqtl <- pcqtl_chr[pcqtl_chr$cluster_id == cluster_id, ]
  
  # Get SNP list and LD matrix
  snp_list <- get_snp_list(cluster_eqtl, ld_path_head, cluster_id)
  cleaned_ld_matrix <- get_ld(ld_path_head, cluster_id, snp_list, genotype_stem)
  ld_snp_set <- rownames(cleaned_ld_matrix)

  # Clean the eQTL and PCQTL data for colocalization
  eqtls_for_coloc <- split_qtl_for_coloc(cluster_eqtl, ld_snp_set, cleaned_ld_matrix, num_gtex_samples)
  pcqtls_for_coloc <- split_qtl_for_coloc(cluster_pcqtl, ld_snp_set, cleaned_ld_matrix, num_gtex_samples)
  
  # Combine into one list and remove NULL entries
  qtls_for_coloc <- c(eqtls_for_coloc, pcqtls_for_coloc)
  qtls_for_coloc <- qtls_for_coloc[!sapply(qtls_for_coloc, is.null)]
  
  cat(length(qtls_for_coloc), "QTL phenotypes\n")
  cat("Processing time:", Sys.time() - start, "\n")
  
  if (length(qtls_for_coloc) == 0) {
    cat("No pairs of QTL signals in this cluster\n")
    return(NULL)
  }

  # Run SuSiE on each QTL phenotype
  qtl_susies <- list()
  for (i in 1:length(qtls_for_coloc)) {
    this_qtl_id <- qtls_for_coloc[[i]]$phenotype_id
    
    # Handle long QTL IDs
    if (nchar(this_qtl_id) > 150) {
      out_qtl_id <- get_short_qtl_id(this_qtl_id)
    } else {
      out_qtl_id <- this_qtl_id
    }
    
    susie_path <- paste(coloc_temp_path_head, out_qtl_id, '.susie.rds', sep = "")
    
    if (file.exists(susie_path)) {
      cat("SuSiE for", this_qtl_id, "already exists\n")
      qtl_susies[[i]] <- readRDS(file = susie_path)
    } else {
      cat("Running SuSiE for", this_qtl_id, "\n")
      cat("Time elapsed:", Sys.time() - start, "\n")
      
      this_qtl_susie <- runsusie_errorcatch(qtls_for_coloc[[i]])
      qtl_susies[[i]] <- this_qtl_susie
      
      cat("Fine-mapped", this_qtl_id, "\n")
      cat("Time elapsed:", Sys.time() - start, "\n")
      cat("Saving SuSiE to", susie_path, "\n")
      saveRDS(this_qtl_susie, file = susie_path)
    } 
  }
  
  # Remove NULL SuSiE results
  qtls_for_coloc <- qtls_for_coloc[!sapply(qtl_susies, is.null)]
  qtl_susies <- qtl_susies[!sapply(qtl_susies, is.null)]
  
  if (length(qtl_susies) < 2) {
    cat("No pairs of QTL signals fine-mapped in this cluster\n")
    return(NULL)
  } else {
    cat(length(qtl_susies), "SuSiE results found in this cluster\n")
  }
  
  # Run colocalization for each pair
  qtl_coloc_results <- get_qtl_pairwise_coloc(qtls_for_coloc, qtl_susies)
  return(qtl_coloc_results)
}

#' Get pairwise colocalization results for QTL pairs
#' 
#' @param qtls_for_coloc List of QTL data for colocalization
#' @param qtl_susies List of SuSiE results
#' @return Colocalization results for all QTL pairs
get_qtl_pairwise_coloc <- function(qtls_for_coloc, qtl_susies) {
  qtl_coloc_results <- data.frame()
  
  if (length(qtl_susies) > 1) {
    qtl_pair_idxs <- combn(seq_along(qtl_susies), 2, simplify = TRUE)
    
    for (i in 1:ncol(qtl_pair_idxs)) {
      qtl1 <- qtls_for_coloc[[qtl_pair_idxs[1, i]]]
      qtl2 <- qtls_for_coloc[[qtl_pair_idxs[2, i]]]
      susie1 <- qtl_susies[[qtl_pair_idxs[1, i]]]
      susie2 <- qtl_susies[[qtl_pair_idxs[2, i]]]
      
      cat("Time elapsed:", Sys.time() - start, "\n")
      cat(paste("Colocalization for", qtl1$phenotype_id, "-", qtl2$phenotype_id), "\n")
      
      this_coloc <- coloc.susie(susie1, susie2)$summary
      
      cat("Time elapsed:", Sys.time() - start, "\n")
      
      if (!is.null(this_coloc)) {
        this_coloc$qtl1_id <- qtl1$phenotype_id
        this_coloc$qtl2_id <- qtl2$phenotype_id
        qtl_coloc_results <- bind_rows(qtl_coloc_results, this_coloc)
      }
      
      if (is.null(this_coloc)) {
        cat("\tResult is null\n")
      }
    }
  }
  
  # Reorder columns to put qtl1_id and qtl2_id first
  if (nrow(qtl_coloc_results) > 0) {
    qtl_coloc_results <- qtl_coloc_results[, c("qtl1_id", "qtl2_id", 
                                               setdiff(colnames(qtl_coloc_results), c("qtl1_id", "qtl2_id")))]
  }
  
  return(qtl_coloc_results)
}

#' Run colocalization analysis between GWAS and QTLs in a cluster
#' 
#' @param gwas_with_meta GWAS data with metadata
#' @param eqtl_chr eQTL data for chromosome
#' @param pcqtl_chr PCQTL data for chromosome
#' @param cluster_id Cluster ID
#' @param ld_path_head Path to LD files
#' @param coloc_temp_path_head Path to temporary colocalization files
#' @param genotype_stem Path to genotype files
#' @param num_gtex_samples Number of GTEX samples
#' @param use_susie Whether to use SuSiE for fine-mapping
#' @return Colocalization results between GWAS and QTLs
coloc_gwas_cluster <- function(gwas_with_meta, eqtl_chr, pcqtl_chr, cluster_id, ld_path_head, 
                              coloc_temp_path_head, genotype_stem, num_gtex_samples, use_susie = FALSE) {
  start <- Sys.time()
  
  # Subset eQTL and PCQTL to this cluster
  cluster_eqtl <- eqtl_chr[eqtl_chr$cluster_id == cluster_id, ]
  cluster_pcqtl <- pcqtl_chr[pcqtl_chr$cluster_id == cluster_id, ]
  
  # Get SNP list and LD matrix
  snp_list <- get_snp_list(cluster_eqtl, ld_path_head, cluster_id)
  cleaned_ld_matrix <- get_ld(ld_path_head, cluster_id, snp_list, genotype_stem)
  ld_snp_set <- rownames(cleaned_ld_matrix)
  
  # Clean the eQTL and PCQTL data for colocalization
  eqtls_for_coloc <- split_qtl_for_coloc(cluster_eqtl, ld_snp_set, cleaned_ld_matrix, num_gtex_samples)
  pcqtls_for_coloc <- split_qtl_for_coloc(cluster_pcqtl, ld_snp_set, cleaned_ld_matrix, num_gtex_samples)
  
  # Combine into one list and remove NULL entries
  qtls_for_coloc <- c(eqtls_for_coloc, pcqtls_for_coloc)
  qtls_for_coloc <- qtls_for_coloc[!sapply(qtls_for_coloc, is.null)]
  
  cat(length(qtls_for_coloc), "QTL phenotypes\n")
  
  if (length(qtls_for_coloc) == 0) {
    cat("No QTL signals in this cluster\n")
    return(NULL)
  }
  
  # Clean the GWAS data for colocalization
  gwas_for_coloc <- get_gwas_for_coloc(gwas_with_meta, ld_snp_set, snp_list, cleaned_ld_matrix)  
  
  # Check if GWAS signal was filtered out
  if (is.null(gwas_for_coloc)) {
    cat("Significant GWAS signal filtered out\n")
    return(NULL)
  }
  
  if (use_susie) {
    # Run SuSiE on each QTL phenotype
    qtl_susies <- list()
    for (i in 1:length(qtls_for_coloc)) {
      this_qtl_id <- qtls_for_coloc[[i]]$phenotype_id
      
      # Handle long QTL IDs
      if (nchar(this_qtl_id) > 150) {
        out_qtl_id <- get_short_qtl_id(this_qtl_id)
      } else {
        out_qtl_id <- this_qtl_id
      }
      
      susie_path <- paste(coloc_temp_path_head, out_qtl_id, '.susie.rds', sep = "")
      
      if (file.exists(susie_path)) {
        cat("SuSiE for", this_qtl_id, "already exists\n")
        qtl_susies[[i]] <- readRDS(file = susie_path)
      } else {
        cat("Time elapsed:", Sys.time() - start, "\n")
        cat("Running SuSiE for", this_qtl_id, "\n")
        
        this_qtl_susie <- runsusie_errorcatch(qtls_for_coloc[[i]])
        qtl_susies[[i]] <- this_qtl_susie
        
        cat("Fine-mapped", this_qtl_id, "\n")
        cat("Time elapsed:", Sys.time() - start, "\n")
        saveRDS(this_qtl_susie, file = susie_path)
      } 
    }
    
    # Remove NULL SuSiE results
    qtls_for_coloc <- qtls_for_coloc[!sapply(qtl_susies, is.null)]
    qtl_susies <- qtl_susies[!sapply(qtl_susies, is.null)]
    
    if (length(qtl_susies) == 0) {
      cat("No QTL signals fine-mapped in this cluster\n")
      return(NULL)
    } else {
      cat(length(qtl_susies), "SuSiE results found in this cluster\n")
    }
    
    # Get the GWAS SuSiE result
    if (nchar(cluster_id) > 150) {
      out_cluster_id <- get_short_cluster_id(cluster_id)
    } else {
      out_cluster_id <- cluster_id
    }
    
    gwas_susie_path <- paste(coloc_temp_path_head, out_cluster_id, '.', gwas_with_meta$gwas_id, '.susie.rds', sep = "")
    
    if (file.exists(gwas_susie_path)) {
      gwas_susie <- readRDS(file = gwas_susie_path)
      cat("Loaded previous GWAS SuSiE result\n")
      
      if (is.null(gwas_susie)) {
        cat("No GWAS signals fine-mapped in this cluster\n")
        return(NULL)
      }
    } else {
      cat("Time elapsed:", Sys.time() - start, "\n")
      cat("Running SuSiE on GWAS\n")
      
      gwas_susie <- runsusie_errorcatch(gwas_for_coloc)
      
      if (is.null(gwas_susie)) {
        cat("No GWAS signals fine-mapped in this cluster\n")
        return(NULL)
      } else {
        cat("Fine-mapped GWAS\n")
        saveRDS(gwas_susie, file = gwas_susie_path)
      }
    }
  }
  
  cat(length(qtls_for_coloc), "colocalizations found in this cluster\n")
  gwas_coloc_results <- data.frame()

  if (use_susie) {
    for (i in 1:length(qtl_susies)) {
      this_qtl_id <- qtls_for_coloc[[i]]$phenotype_id
      cat("Time elapsed:", Sys.time() - start, "\n")
      cat(paste("\t\t\tColocalization for", gwas_with_meta$gwas_id, "and", this_qtl_id, "\n"))
      cat("\t\t\t\tUsing SuSiE to colocalize", i, "out of", length(qtl_susies), "total\n")
      
      this_qtl_susie <- qtl_susies[[i]]
      this_coloc <- coloc.susie(gwas_susie, this_qtl_susie)$summary
      
      if (!is.null(this_coloc)) {
        this_coloc$gwas_id <- gwas_with_meta$gwas_id
        this_coloc$qtl_id <- this_qtl_id
        gwas_coloc_results <- bind_rows(gwas_coloc_results, this_coloc)
      }
    }
  } else {
    for (i in 1:length(qtls_for_coloc)) {
      this_qtl_id <- qtls_for_coloc[[i]]$phenotype_id
      cat(paste("\t\t\tColocalization for", gwas_with_meta$gwas_id, "and", this_qtl_id, "\n"))
      cat("\t\t\t\tUsing ABF to colocalize", i, "out of", length(qtls_for_coloc), "total\n")
      
      this_coloc <- coloc.abf(gwas_for_coloc, qtls_for_coloc[[i]])$summary
      
      if (!is.null(this_coloc)) {
        this_coloc$gwas_id <- gwas_with_meta$gwas_id
        this_coloc$qtl_id <- this_qtl_id
      }
      gwas_coloc_results <- bind_rows(gwas_coloc_results, this_coloc)
    }
  }
  

  
  return(gwas_coloc_results)
}

#' Check if GWAS has a signal in a cluster
#' 
#' @param gwas_chr GWAS data for chromosome
#' @param this_cluster Cluster information
#' @return TRUE if GWAS has signal in cluster, FALSE otherwise
check_gwas_cluster <- function(gwas_chr, this_cluster) {
  # Check for GWAS signals within 1Mb of cluster boundaries
  gwas_cluster <- gwas_chr[gwas_chr$position > this_cluster$start - 1e6, ] 
  gwas_cluster <- gwas_cluster[gwas_cluster$position < this_cluster$end + 1e6, ]
  return(sum(gwas_cluster$pvalue < 1e-6) > 0)
}

# =============================================================================
# LD and SNP list functions
# =============================================================================

#' Get LD matrix for a cluster
#' 
#' @param ld_path_head Path to LD files
#' @param cluster_id Cluster ID
#' @param snp_list List of SNPs
#' @param genotype_stem Path to genotype files
#' @return Cleaned LD matrix
get_ld <- function(ld_path_head, cluster_id, snp_list, genotype_stem) {
  if (nchar(cluster_id) > 150) {
    cluster_id <- get_short_cluster_id(cluster_id)
  }
  
  # Check if LD matrix already exists
  ld_matrix_path <- paste(ld_path_head, cluster_id, '.ld', sep = "")
  
  if (file.exists(ld_matrix_path)) {
    cat("LD matrix already exists\n")
  } else {    
    cat("LD matrix does not exist, generating...\n")
    
    # Generate LD matrix using PLINK
    ld_plink_path <- paste(ld_path_head, cluster_id, sep = "")
    snp_path <- paste(ld_path_head, cluster_id, '.snp_list.txt', sep = "")
    plink_command <- sprintf("plink --bfile %s --extract %s --r square --out %s --keep-allele-order", 
                           genotype_stem, snp_path, ld_plink_path)
    cat("PLINK command:", plink_command, "\n")
    system(plink_command, intern = TRUE)
    cat("Generated LD matrix\n")
  }
  
  # Load LD matrix with error handling
  tryCatch({
    ld_matrix <- fread(ld_matrix_path)
  }, error = function(e) {
    cat("Error reading LD matrix:", e$message, "\n")
    cat("Regenerating the file...\n")
    
    ld_plink_path <- paste(ld_path_head, cluster_id, sep = "")
    snp_path <- paste(ld_path_head, cluster_id, '.snp_list.txt', sep = "")
    plink_command <- sprintf("plink --bfile %s --extract %s --r square --out %s --keep-allele-order", 
                           genotype_stem, snp_path, ld_plink_path)
    cat("PLINK command:", plink_command, "\n")
    system(plink_command, intern = TRUE)
    cat("Generated LD matrix\n")
    ld_matrix <- fread(ld_matrix_path)
  })

  # Convert to data frame and set row/column names
  ld_matrix <- data.frame(ld_matrix)
  rownames(ld_matrix) <- snp_list$variant_id
  colnames(ld_matrix) <- snp_list$variant_id
  
  # Remove SNPs with missing values from LD matrix
  ld_missing_snps <- get_ld_missing_snps(ld_matrix)
  cleaned_ld_matrix <- ld_matrix[!rownames(ld_matrix) %in% ld_missing_snps, 
                                !colnames(ld_matrix) %in% ld_missing_snps]
  
  cat("Total working SNPs:", nrow(cleaned_ld_matrix), "\n")
  return(cleaned_ld_matrix)
}

#' Get SNPs with missing values in LD matrix
#' 
#' @param ld_matrix LD matrix
#' @return Vector of SNP IDs with missing values
get_ld_missing_snps <- function(ld_matrix) {
  rows_with_missing <- rownames(ld_matrix)[which(rowSums(is.na(ld_matrix)) > 0)]
  cols_with_missing <- colnames(ld_matrix)[which(colSums(is.na(ld_matrix)) > 0)]
  ld_missing_snps <- union(rows_with_missing, cols_with_missing)
  return(ld_missing_snps)
}

#' Get SNP list for a cluster
#' 
#' @param cluster_eqtl eQTL data for cluster
#' @param ld_path_head Path to LD files
#' @param cluster_id Cluster ID
#' @return Data frame with SNP list
get_snp_list <- function(cluster_eqtl, ld_path_head, cluster_id) {
  if (nchar(cluster_id) > 150) {
    cluster_id <- get_short_cluster_id(cluster_id)
  }
  
  snp_path <- paste(ld_path_head, cluster_id, '.snp_list.txt', sep = "")
  
  # Generate SNP list if it doesn't exist
  if (file.exists(snp_path)) {
    cat("SNP list already exists\n")
  } else {
    snp_list <- unique(cluster_eqtl$variant_id)
    cat("Generated SNP list with", length(snp_list), "SNPs\n")
    write.table(snp_list, file = snp_path, row.names = FALSE, sep = "\t", 
                col.names = FALSE, quote = FALSE)
  }
  
  # Read SNP list
  snp_list <- read.table(snp_path, header = FALSE, sep = "\t")
  colnames(snp_list) <- "variant_id"
  cat(nrow(snp_list), "SNPs found for this cluster\n")
  return(snp_list)
}

# =============================================================================
# Utility functions
# =============================================================================

#' Get shortened cluster ID for file naming
#' 
#' @param long_cluster_id Long cluster ID
#' @return Shortened cluster ID
get_short_cluster_id <- function(long_cluster_id) {
  # Filenames can only be 260 characters, so shorten long cluster IDs
  # Take first 5 parts of the cluster ID
  long_cluster_id_split <- strsplit(long_cluster_id, "_")[[1]]
  short_cluster_id <- paste(long_cluster_id_split[1:5], collapse = "_")
  return(short_cluster_id)
}

#' Get shortened QTL ID for file naming
#' 
#' @param long_qtl_id Long QTL ID
#' @return Shortened QTL ID
get_short_qtl_id <- function(long_qtl_id) {
  # Filenames can only be 260 characters, so shorten long QTL IDs
  # Take last 5 parts of the QTL ID to avoid conflicts within clusters
  long_qtl_id_split <- strsplit(long_qtl_id, "_")[[1]]
  short_qtl_id <- paste(tail(long_qtl_id_split, 5), collapse = "_")
  return(short_qtl_id)
}

#' Load GWAS data from file path with metadata
#' 
#' @param gwas_path Path to GWAS data file
#' @param gwas_meta_df GWAS metadata dataframe
#' @param gwas_id GWAS ID
#' @return List containing GWAS data and metadata
load_gwas_from_path <- function(gwas_path, gwas_meta_df, gwas_id) {
  # Get GWAS metadata
  gwas_meta <- gwas_meta_df[gwas_meta_df$Tag == gwas_id, ]
  num_gwas_samples <- gwas_meta$Sample_Size
  
  # Determine GWAS type
  if (gwas_meta$Binary == 0) {
    gwas_type <- 'quant'
  } else {
    gwas_type <- 'cc'
  }
  
  cat(paste("\t\tReading GWAS:", Sys.time() - start), "\n")
  gwas <- fread(gwas_path)
  cat(paste("\t\tFinished reading GWAS:", Sys.time() - start), "\n")
  
  return(list("gwas_data" = gwas, 
              "gwas_id" = gwas_id, 
              "gwas_type" = gwas_type, 
              "gwas_sample_size" = num_gwas_samples))
}

#' Load eQTL data for a chromosome
#' 
#' @param eqtl_dir_path Path to eQTL directory
#' @param chr_id Chromosome ID
#' @param tissue_id Tissue ID
#' @return eQTL data with cluster IDs
get_eqtl_chr <- function(eqtl_dir_path, chr_id, tissue_id) {
  eqtl_path <- paste(eqtl_dir_path, "/", tissue_id, '.v8.cluster_genes.cis_qtl_pairs.chr', chr_id, '.parquet', sep = "")
  eqtl <- nanoparquet::read_parquet(eqtl_path)
  eqtl$cluster_id <- sapply(eqtl$phenotype_id, function(x) unlist(strsplit(as.character(x), '_e_'))[1])
  return(eqtl)
}

#' Load PCQTL data for a chromosome
#' 
#' @param pcqtl_dir_path Path to PCQTL directory
#' @param chr_id Chromosome ID
#' @param tissue_id Tissue ID
#' @return PCQTL data with cluster IDs
get_pcqtl_chr <- function(pcqtl_dir_path, chr_id, tissue_id) {
  pcqtl_path <- paste(pcqtl_dir_path, "/", tissue_id, '.v8.pcs.cis_qtl_pairs.chr', chr_id, '.parquet', sep = "")
  pcqtl <- nanoparquet::read_parquet(pcqtl_path)
  pcqtl$cluster_id <- sapply(pcqtl$phenotype_id, function(x) unlist(strsplit(as.character(x), '_pc'))[1])
  return(pcqtl)
}

#' Extract data from SuSiE RDS object
#' 
#' @param rds_path Path to SuSiE RDS file
#' @return Data frame with extracted SuSiE results
extract_data_from_rds <- function(rds_path) {
  # Load SuSiE object
  susie_object <- readRDS(rds_path)
  
  # Extract phenotype ID from filename
  phenotype_id <- gsub("\\.susie\\.rds$", "", basename(rds_path))
  
  # Initialize results list
  results <- list()
  
  # Loop through each credible set
  for (cs_name in names(susie_object$sets$cs)) {
    variant_ids <- names(susie_object$sets$cs[[cs_name]])
    pip_values <- susie_object$pip[variant_ids]
    
    # Assign CS IDs numerically
    cs_id <- as.numeric(gsub("L", "", cs_name))
    
    # Create data frame for current credible set
    cs_results <- data.frame(
      phenotype_id = phenotype_id,
      variant_id = variant_ids,
      pip = pip_values,
      cs_id = cs_id
    )
    
    # Add to results list
    results[[cs_name]] <- cs_results
  }
  
  # Combine all results
  combined_results <- bind_rows(results)
  return(combined_results)
}