setwd('/home/klawren/oak/pcqtls/')
library(coloc)
library(nanoparquet)
library(tidyverse)
library(data.table)
library(Rfast)
library(argparse)
library(tidyverse)
library(susieR)
library(dplyr)
#library(pecotmr)

####### funtions #########

split_qtl_for_coloc <- function(cluster_qtl, ld_snp_set, cleaned_ld_matrix, num_gtex_samples){
  qtl_filtered <- filter_qtl(cluster_qtl, ld_snp_set)
  # split eqtl to each of the egenes and clean up
  phenotype_ids <-  unique(qtl_filtered$phenotype_id)
  qtls_for_coloc <- list()
  for (i in 1:length(phenotype_ids)){
    phenotype_qtl <- qtl_filtered[qtl_filtered$phenotype_id == phenotype_ids[i], ]
    cat(paste("\tlooking for signals in",phenotype_ids[i]), "\n")
    cat(paste("\t\tqtl filtered to :", nrow(phenotype_qtl), "snps\n"))
    phenotype_qtl_for_coloc <- clean_qtl(phenotype_qtl, cleaned_ld_matrix, num_gtex_samples)
    qtls_for_coloc[[i]] <- phenotype_qtl_for_coloc
  }
  return(qtls_for_coloc)
}

# filter eqtl data
filter_qtl <- function(qtl, ld_snp_set){
  qtl_filtered <- qtl[qtl$variant_id %in% ld_snp_set, ]
  return(qtl_filtered)
}

clean_qtl <- function(qtl_filtered, cleaned_ld_matrix, num_gtex_samples){
  # some slope_se are na
  qtl_missing_snps <- qtl_filtered[is.na(qtl_filtered$slope) | is.na(qtl_filtered$slope_se) | (qtl_filtered$af==0), 'variant_id']
  cat("\t\tNumber snps with  qtl missing: ", length(qtl_missing_snps), "\n")
  cleaned_qtl <- qtl_filtered[!qtl_filtered %in% qtl_missing_snps, ]
  # if there isn't a signal, further analysis should not be run
  if (min(cleaned_qtl$pval_nominal) < 1e-6){
    cat('\t\tsignal found \n')
    # make some columns to play well with coloc
    cleaned_qtl <- cleaned_qtl %>%
      mutate(position = as.integer(str_split(variant_id, "_") %>% sapply(pluck, 2)))
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
                             LD = as.matrix(cleaned_ld_matrix))
    return(cleaned_qtl_list)
  } else{
    cat('\t no signal found \n')
    return(NULL)
  }    
}

get_gwas_for_coloc <- function(gwas_with_meta, ld_snp_set, snp_list, cleaned_ld_matrix){
  gwas_filtered <- filter_gwas(gwas_with_meta$gwas_data, snp_list, ld_snp_set)
  gwas_for_coloc <- clean_gwas(gwas_filtered, cleaned_ld_matrix, gwas_with_meta$gwas_type, gwas_with_meta$gwas_sample_size)
  return(gwas_for_coloc)
}

filter_gwas <- function(gwas_data, snp_list, ld_snp_set){
  gwas_filtered <- gwas_data[gwas_data$panel_variant_id %in% snp_list$variant_id, ]
  gwas_filtered <- gwas_filtered[gwas_filtered$panel_variant_id %in% ld_snp_set, ]
  cat("snps in gwas: ")
  cat(nrow(gwas_filtered))
  cat('\n')
  return(gwas_filtered)
}

clean_gwas <- function(gwas_filtered, cleaned_ld_matrix, gwas_type, num_gwas_samples){
  # susie needs effect sizes, so we must also drop the snps with na for gwas effect
  gwas_missing_snps <- gwas_filtered[is.na(gwas_filtered$effect_size) | is.na(gwas_filtered$standard_error) | (gwas_filtered$frequency<= 0) | (gwas_filtered$frequency >= 1), 'panel_variant_id']
  cat("\t Number snps with ld or gwas missing: ", length(gwas_missing_snps$panel_variant_id), "\n")
  cleaned_gwas <- gwas_filtered[!gwas_filtered$panel_variant_id %in% gwas_missing_snps$panel_variant_id, ]
  if(min(cleaned_gwas$pvalue) < 1e-6){
    gwas_cleaned_ld_matrix <- cleaned_ld_matrix[rownames(cleaned_ld_matrix) %in% cleaned_gwas$panel_variant_id, colnames(cleaned_ld_matrix) %in% cleaned_gwas$panel_variant_id]
    cleaned_gwas_list <- list(MAF = cleaned_gwas$frequency, 
                              snp = cleaned_gwas$panel_variant_id, 
                              position = cleaned_gwas$position, 
                              type = gwas_type, 
                              N = num_gwas_samples, 
                              pvalues=cleaned_gwas$pvalue,
                              LD = as.matrix(gwas_cleaned_ld_matrix), 
                              beta = cleaned_gwas$effect_size,
                              varbeta = cleaned_gwas$standard_error**2)
    cat('\t signal found \n')
    return(cleaned_gwas_list)
  } else {
    cat('\t no signal found, this is odd becuase I already checked and there should be a signal\n')
    return(NULL)
  }
}


runsusie_errorcatch <- function(dataset){
  start <- Sys.time()
  susie <- tryCatch({
    runsusie(dataset, repeat_until_convergence=TRUE) #DEBUG (should I run until convergence or not?) # Call runsusie function for each QTL
  }, error = function(e) {
    cat("An error occurred for QTL \n")
    cat("Error Message:", conditionMessage(e), "\n")  # Print the error message
    NULL  # Return NULL if an error occurs
  })
  cat(Sys.time()-start)
  return(susie)  # Store the result or NULL in the qtl_susies list
}


coloc_pairs_cluster <- function(eqtl_chr, pcqtl_chr, cluster_id, ld_path_head, genotype_stem, num_gtex_samples, coloc_temp_path_head){
  # subset eqtl and pcqtl to this cluster
  start <- Sys.time()
  cluster_eqtl <-eqtl_chr[eqtl_chr$cluster_id == cluster_id, ]
  cluster_pcqtl <- pcqtl_chr[pcqtl_chr$cluster_id == cluster_id, ]
  # get snp list and ld matrix
  snp_list <- get_snp_list(cluster_eqtl, ld_path_head, cluster_id)
  cleaned_ld_matrix <- get_ld(ld_path_head, cluster_id, snp_list, genotype_stem)
  ld_snp_set <- rownames(cleaned_ld_matrix)

  # clean the eqtl and pcqtl data into the right format for susie
  eqtls_for_coloc <- split_qtl_for_coloc(cluster_eqtl, ld_snp_set, cleaned_ld_matrix, num_gtex_samples)
  pcqtls_for_coloc <- split_qtl_for_coloc(cluster_pcqtl, ld_snp_set, cleaned_ld_matrix, num_gtex_samples)
  # combine them into one list
  qtls_for_coloc <- c(eqtls_for_coloc, pcqtls_for_coloc)
  qtls_for_coloc <- qtls_for_coloc[!sapply(qtls_for_coloc, is.null)]
  cat(length(qtls_for_coloc))
  cat(" qtl phenotypes\n")
  cat(Sys.time() - start, "\n")
  if (length(qtls_for_coloc)==0){
    cat("no pairs of qtl signals in this cluster\n")
    return(NULL)
  }

  # run susie on each qtl phenotype
  # check if the susie dataset alread exists 
  qtl_susies <- list()
  for (i in 1:length(qtls_for_coloc)){
    this_qtl_id <- qtls_for_coloc[[i]]$phenotype_id
    if (nchar(this_qtl_id) > 150){
      out_qtl_id <- get_short_qtl_id(this_qtl_id)
    } else {
      out_qtl_id <- this_qtl_id
    }
    susie_path <- paste(coloc_temp_path_head, out_qtl_id, '.susie.rds', sep="")
    if (file.exists(susie_path)) {
      cat("susie for ", this_qtl_id, "already exists\n")
      qtl_susies[[i]] <-  readRDS(file = susie_path)
    } else {
      cat("running susie for ", this_qtl_id,"\n")
      cat(Sys.time() - start, "\n")
      this_qtl_susie <- runsusie_errorcatch(qtls_for_coloc[[i]])
      qtl_susies[[i]] <- this_qtl_susie
      cat("finemapped ", this_qtl_id, "\n")
      cat(Sys.time() - start, "\n")
      saveRDS(this_qtl_susie, file = susie_path)
    } 
  }
  # some susies are null, remove those from consideration
  qtls_for_coloc <- qtls_for_coloc[!sapply(qtl_susies, is.null)]
  qtl_susies <- qtl_susies[!sapply(qtl_susies, is.null)]
  if (length(qtl_susies)<2){
    cat("no pairs of qtl signals finemapped in this cluster\n")
    return(NULL)
  } else {
    cat(length(qtl_susies), "susies found in this cluster\n")
  }
  # colocalize each pair
  qtl_coloc_results <- get_qtl_pairwise_coloc(qtls_for_coloc, qtl_susies)
  return(qtl_coloc_results)
}

get_qtl_pairwise_coloc <- function(qtls_for_coloc, qtl_susies){
  qtl_coloc_results <- data.frame()
  if(length(qtl_susies)>1){
    qtl_pair_idxs <- combn(seq_along(qtl_susies), 2, simplty=TRUE)
    for (i in 1:ncol(qtl_pair_idxs)) {
      qtl1 <- qtls_for_coloc[[qtl_pair_idxs[1, i]]]
      qtl2 <- qtls_for_coloc[[qtl_pair_idxs[2, i]]]
      susie1 <- qtl_susies[[qtl_pair_idxs[1, i]]]
      susie2 <- qtl_susies[[qtl_pair_idxs[2, i]]]
      cat(Sys.time() - start, "\n")
      cat(paste("coloc for", qtl1$phenotype_id , "-", qtl2$phenotype_id), "\n")
      this_coloc <- coloc.susie(susie1,susie2)$summary
      cat(Sys.time() - start, "\n")
      if (!is.null(this_coloc)){
        this_coloc$qtl1_id <- qtl1$phenotype_id
        this_coloc$qtl2_id <- qtl2$phenotype_id
        qtl_coloc_results <- bind_rows(qtl_coloc_results, this_coloc)
      }
      if(is.null(this_coloc)){
        cat("\tresult is null\n")
        cat(ncol(qtl_coloc_results), " vs ", ncol(this_coloc), "\n")
        print(summary(susie1))
        print(summary(susie2))
      } else {
      }
    }
  }
  return(qtl_coloc_results)
}


coloc_gwas_cluster <- function(gwas_with_meta, eqtl_chr, pcqtl_chr, cluster_id, ld_path_head, coloc_temp_path_head, genotype_stem, num_gtex_samples, use_susie = FALSE){
  # subset eqtl and pcqtl to this cluster
  start <- Sys.time()
  cluster_eqtl <-eqtl_chr[eqtl_chr$cluster_id == cluster_id, ]
  cluster_pcqtl <- pcqtl_chr[pcqtl_chr$cluster_id == cluster_id, ]
  # get snp list and ld matrix
  snp_list <- get_snp_list(cluster_eqtl, ld_path_head, cluster_id)
  cleaned_ld_matrix <- get_ld(ld_path_head, cluster_id, snp_list, genotype_stem)
  ld_snp_set <- rownames(cleaned_ld_matrix)
  
  # clean the eqtl and pcqtl data into the right format for susie
  eqtls_for_coloc <- split_qtl_for_coloc(cluster_eqtl, ld_snp_set, cleaned_ld_matrix, num_gtex_samples)
  pcqtls_for_coloc <- split_qtl_for_coloc(cluster_pcqtl, ld_snp_set, cleaned_ld_matrix, num_gtex_samples)
  # combine them into one list
  qtls_for_coloc <- c(eqtls_for_coloc, pcqtls_for_coloc)
  qtls_for_coloc <- qtls_for_coloc[!sapply(qtls_for_coloc, is.null)]
  cat(length(qtls_for_coloc))
  if (length(qtls_for_coloc)==0){
    cat("no qtl signals in this cluster\n")
    return(NULL)
  }
  cat(" qtl phenotypes\n")
  
  # clean the gwas data into the right format for susie
  gwas_for_coloc <- get_gwas_for_coloc(gwas_with_meta, ld_snp_set, snp_list, cleaned_ld_matrix)  
  # in rare cases, the snps with <1e-6 singal are not in the qtl snp_list or are nans in the ld matrix
  if(is.null(gwas_for_coloc)){
    cat("signifigant gwas signal filtered out\n")
    return(NULL)
  }
  if(use_susie){
    # run susie on each qtl phenotype
    # check if the susie dataset alread exists 
    qtl_susies <- list()
    for (i in 1:length(qtls_for_coloc)){
      this_qtl_id <- qtls_for_coloc[[i]]$phenotype_id
      if (nchar(this_qtl_id) > 150){
        out_qtl_id <- get_short_qtl_id(this_qtl_id)
      } else {
        out_qtl_id <- this_qtl_id
      }
      susie_path <- paste(coloc_temp_path_head, out_qtl_id, '.susie.rds', sep="")
      if (file.exists(susie_path)) {
        cat("susie for ", this_qtl_id, "already exists\n")
        qtl_susies[[i]] <-  readRDS(file = susie_path)
      } else {
        cat(Sys.time() - start, "\n")
        cat("susie for ", this_qtl_id, "\n")
        this_qtl_susie <- runsusie_errorcatch(qtls_for_coloc[[i]])
        qtl_susies[[i]] <- this_qtl_susie
        cat("finemapped ", this_qtl_id, "\n")
        cat(Sys.time() - start, "\n")
        saveRDS(this_qtl_susie, file = susie_path)
      } 
    }
    # some susies are null, remove those from consideration
    qtls_for_coloc <- qtls_for_coloc[!sapply(qtl_susies, is.null)]
    qtl_susies <- qtl_susies[!sapply(qtl_susies, is.null)]
    if (length(qtl_susies)==0){
      cat("no qtl signals finemapped in this cluster\n")
      return(NULL)
    } else {
      cat(length(qtl_susies), "susies found in this cluster\n")
    }
    # get the gwas susie (this is specific to this cluster and gwas, but maybe still write it out?)
    if (nchar(cluster_id) > 150){
      out_cluster_id <- get_short_cluster_id(cluster_id)
    }
    else{
      out_cluster_id <- cluster_id
    }
    gwas_susie_path <- paste(coloc_temp_path_head, out_cluster_id, '.', gwas_with_meta$gwas_id, '.susie.rds', sep="")
    if (file.exists(gwas_susie_path)) {
      gwas_susie <- readRDS(file = gwas_susie_path)
      cat("loaded previous gwas susie\n")
      if (is.null(gwas_susie)){
        cat("no gwas signals finemapped in this cluster\n")
        return(NULL)
      }
    } else {
      cat(Sys.time() - start, "\n")
      cat("running susie gwas\n")
      gwas_susie <- runsusie_errorcatch(gwas_for_coloc)
      if (is.null(gwas_susie)){
        cat("no gwas signals finemapped in this cluster\n")
        return(NULL)
      } else {
      cat("finemapped gwas\n")
      saveRDS(gwas_susie, file = gwas_susie_path)
      }
    }
  }
  cat(length(qtls_for_coloc), "colocs found in this cluster\n")
  gwas_coloc_results <- data.frame()

  if(use_susie){
    for (i in 1:length(qtl_susies)){
      this_qtl_id <- qtls_for_coloc[[i]]$phenotype_id
      cat(Sys.time() - start, "\n")
      cat(paste("\t\t\tcoloc for", gwas_with_meta$gwas_id, "and", this_qtl_id, "\n"))
      cat("\t\t\t\t using susie to coloc ", i, " out of ", length(qtl_susies), " total\n")
      this_qtl_susie <- qtl_susies[[i]]
      this_coloc <- coloc.susie(gwas_susie,this_qtl_susie)$summary
      if(!is.null(this_coloc)){
        this_coloc$gwas_id <- gwas_with_meta$gwas_id
        this_coloc$qtl_id <- this_qtl_id
        gwas_coloc_results <- bind_rows(gwas_coloc_results, this_coloc) # should this be outisde the if?
      }
    }
  } else {
    for (i in 1:length(qtls_for_coloc)){
      this_qtl_id <- qtls_for_coloc[[i]]$phenotype_id
      cat(paste("\t\t\tcoloc for", gwas_with_meta$gwas_id, "and", this_qtl_id, "\n"))
      cat("\t\t\t\t abf to coloc ", i, " out of ", length(qtls_for_coloc), " total\n")
      this_coloc <- coloc.abf(gwas_for_coloc, qtls_for_coloc[[i]])$summary
      if(!is.null(this_coloc)){
        this_coloc$gwas_id <- gwas_with_meta$gwas_id
        this_coloc$qtl_id <- this_qtl_id
      }
      gwas_coloc_results <- bind_rows(gwas_coloc_results, this_coloc)
    }
  }
  # return the results
  gwas_coloc_results$gwas_cs_is <- gwas_coloc_results$idx1 
  gwas_coloc_results$qtl_cs_is <- gwas_coloc_results$idx2
  return(gwas_coloc_results)
}


check_gwas_cluster <- function(gwas_chr, this_cluster){
  # return true if the gwas has a signal in this cluster, false otherwise
  gwas_cluster <- gwas_chr[gwas_chr$position > this_cluster$start - 1e6, ] 
  gwas_cluster <- gwas_cluster[gwas_cluster$position < this_cluster$end + 1e6, ]
  return(sum(gwas_cluster$pvalue < 1e-6)>0)
}

get_ld <- function(ld_path_head, cluster_id, snp_list, genotype_stem){
  if (nchar(cluster_id) > 150){
    cluster_id <- get_short_cluster_id(cluster_id)
  }
  # check if ld already exists
  ld_matrix_path <- paste(ld_path_head, cluster_id, '.ld', sep="")
  if (file.exists(ld_matrix_path)) {
    cat("ld matrix already exists\n")
  } else {    
    cat("ld matrix does not already exist\n")
    # get ld if not
    ld_plink_path <- paste(ld_path_head, cluster_id, sep="")
    snp_path <- paste(ld_path_head, cluster_id, '.snp_list.txt', sep="")
    plink_command <- sprintf("plink --bfile %s --extract %s --r square --out %s --keep-allele-order", genotype_stem, snp_path, ld_plink_path)
    cat(plink_command) 
    system(plink_command, intern=TRUE)
    cat("generated ld matrix\n")
  }
  # load in ld and return 
  #ld_matrix <- read.table(ld_matrix_path)
  ld_matrix <- fread(ld_matrix_path)
  ld_matrix <- data.frame(ld_matrix)
  rownames(ld_matrix) <- snp_list$variant_id
  colnames(ld_matrix) <- snp_list$variant_id
  
  # drop snps with missing values from LD
  ld_missing_snps <- get_ld_missing_snps(ld_matrix)
  cleaned_ld_matrix <- ld_matrix[!rownames(ld_matrix) %in% ld_missing_snps, !colnames(ld_matrix) %in% ld_missing_snps]
  cat("total working snps: ", nrow(cleaned_ld_matrix), "\n")
  return(cleaned_ld_matrix)
}

get_ld_missing_snps <- function(ld_matrix){
  # drop snps with missing values from LD for eqtl
  rows_with_missing <- rownames(ld_matrix)[which(rowSums(is.na(ld_matrix)) > 0)]
  cols_with_missing <- colnames(ld_matrix)[which(colSums(is.na(ld_matrix)) > 0)]
  ld_missing_snps = union(rows_with_missing, cols_with_missing)
  return(ld_missing_snps)
}

get_short_cluster_id <- function(long_cluster_id){
  # filenames can only be 260 characters
  # this leads to erros for large clusters
  # I just take the first 5 of the cluster and call that good. 
  # There shouldn't be any overlap in transcripts anyway
  long_cluster_id_split <- strsplit(long_cluster_id, "_")[[1]]
  # Take all text before the 5th '_'
  short_cluster_id <- paste(long_cluster_id_split[1:5], collapse = "_")
  return(short_cluster_id)
}


get_short_qtl_id <- function(long_qtl_id){
  # filenames can only be 260 characters
  # this leads to erros for large clusters
  # I just take the last 5 of the qtl and call that good.
  # needs to be last so multiple qtls in a cluster aren't messed up 
  # There shouldn't be any overlap in transcripts anyway
  long_qtl_id_split <- strsplit(long_qtl_id, "_")[[1]]
  # Take all text before the 5th '_'
  short_qtl_id <- paste(tail(long_qtl_id_split, 5), collapse = "_")
  return(short_qtl_id)
}

get_snp_list <- function(cluster_eqtl, ld_path_head, cluster_id){
  if (nchar(cluster_id) > 150){
    cluster_id <- get_short_cluster_id(cluster_id)
  }
  snp_path <- paste(ld_path_head, cluster_id, '.snp_list.txt', sep="")
  # make snp list if not
  if (file.exists(snp_path)) {
    cat("snp list already exists\n")
  } else {
    snp_list <- unique(cluster_eqtl$variant_id)
    length(snp_list)
    cat("generated snp list\n")
    write.table(snp_list, file=snp_path, row.names=FALSE, sep="\t", col.names = FALSE, quote=FALSE)
  }
  snp_list <- read.table(snp_path, header = FALSE, sep="\t")
  colnames(snp_list) <- "variant_id"
  cat(nrow(snp_list), " snps found for this cluster\n")
  return(snp_list)
}

load_gwas_from_path <- function(gwas_path, gwas_meta_df, gwas_id){
  # return gwas object with relevant metadata
  gwas_meta <- gwas_meta_df[gwas_meta_df$Tag == gwas_id, ]
  num_gwas_samples <- gwas_meta$Sample_Size
  if (gwas_meta$Binary == 0) {
    gwas_type <- 'quant'
  } else {
    gwas_type <- 'cc'
  }
  cat(paste("\t\treading gwas:", Sys.time()-start), "\n")
  gwas <- fread(gwas_path)
  cat(paste("\t\tfinished reading gwas:", Sys.time()-start), "\n")
  return(list("gwas_data" = gwas, "gwas_id" = gwas_id, "gwas_type" = gwas_type, "gwas_sample_size" = num_gwas_samples))
}

get_eqtl_chr <- function(eqtl_dir_path, chr_id, tissue_id){
  eqtl_path <- paste(eqtl_dir_path, "/", tissue_id, '.v8.cluster_genes.cis_qtl_pairs.chr', chr_id, '.parquet', sep="")
  eqtl <- nanoparquet::read_parquet(eqtl_path)
  eqtl$cluster_id <- sapply(eqtl$phenotype_id, function(x) unlist(strsplit(as.character(x), '_e_'))[1])
  return(eqtl)
}

get_pcqtl_chr <- function(pcqtl_dir_path, chr_id, tissue_id){
  pcqtl_path <- paste(pcqtl_dir_path, "/", tissue_id, '.v8.pcs.cis_qtl_pairs.chr', chr_id, '.parquet', sep="")
  pcqtl <- nanoparquet::read_parquet(pcqtl_path)
  pcqtl$cluster_id <- sapply(pcqtl$phenotype_id, function(x) unlist(strsplit(as.character(x), '_pc'))[1])
  return(pcqtl)
}


# Function to extract and format data from an RDS object
extract_data_from_rds <- function(rds_path) {
  # Load RDS object
  susie_object <- readRDS(rds_path)
  # Extract necessary components
  phenotype_id <- gsub("\\.susie\\.rds$", "", basename(rds_path))
  # Initialize a list to store results
  results <- list()
  # Loop through each credible set
  for (cs_name in names(susie_object$sets$cs)) {
    variant_ids <- names(susie_object$sets$cs[[cs_name]])  # Extract variant IDs for current credible set
    pip_values <- susie_object$pip[variant_ids]            # Extract corresponding PIP values
    # Assign CS IDs numerically
    cs_id <- as.numeric(gsub("L", "", cs_name))
    # Create a data frame for the current credible set
    cs_results <- data.frame(
      phenotype_id = phenotype_id,
      variant_id = variant_ids,
      pip = pip_values,
      cs_id = cs_id
    )
    # Combine current results into the results list
    results[[cs_name]] <- cs_results
  }
  # Combine all results from credible sets into one data frame
  combined_results <- bind_rows(results)
  return(combined_results)  # Return the combined data frame
}