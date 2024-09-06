
library(argparse)
library(coloc)
library(arrow)
library(tidyverse)
library(data.table)
library(Rfast)


########### functions #########

get_ld_missing_snps <- function(ld_matrix){
  # drop snps with missing values from LD for eqtl
  rows_with_missing <- rownames(ld_matrix)[which(rowSums(is.na(ld_matrix)) > 0)]
  cols_with_missing <- colnames(ld_matrix)[which(colSums(is.na(ld_matrix)) > 0)]
  ld_missing_snps = union(rows_with_missing, cols_with_missing)
  return(ld_missing_snps)
}

clean_eqtl <- function(eqtl_filtered, cleaned_ld_matrix, num_gtex_samples){
  # make some columns to play well with coloc
  eqtl_filtered <- eqtl_filtered %>%
    mutate(position = as.integer(str_split(variant_id, "_") %>% sapply(pluck, 2)))
  cleaned_eqtl_list <- list(beta = eqtl_filtered$slope, 
                            varbeta = eqtl_filtered$slope_se**2, 
                            snp = eqtl_filtered$variant_id, 
                            position = eqtl_filtered$position, 
                            pvalues = eqtl_filtered$pval_nominal,
                            type = "quant", 
                            N = num_gtex_samples, 
                            MAF = eqtl_filtered$af,
                            sdY = 1,
                            phenotype_id = eqtl_filtered$phenotype_id[1],
                            LD = as.matrix(cleaned_ld_matrix))
  #print(check_dataset(cleaned_eqtl_list,req="LD"))
  # if there isn't a signal, further analysis should not be run
  if (min(cleaned_eqtl_list$pvalues) < 1e-6){
    cat('\t signal found \n')
    return(cleaned_eqtl_list)
  } else{
    cat('\t no signal found \n')
  }    
}

clean_gwas <- function(gwas_filtered, cleaned_ld_matrix, gwas_type, num_gwas_samples){
  # susie needs effect sizes, so we must also drop the snps with na for gwas effect
  gwas_missing_snps <- gwas_filtered[is.na(gwas_filtered$effect_size), 'panel_variant_id']
  cat("\t Number snps with ld or gwas missing: ")
  cat(length(gwas_missing_snps$panel_variant_id))
  cat("\n")
  cleaned_gwas <- gwas_filtered[!gwas_filtered$panel_variant_id %in% gwas_missing_snps$panel_variant_id, ]
  gwas_cleaned_ld_matrix <- cleaned_ld_matrix[!rownames(cleaned_ld_matrix) %in% gwas_missing_snps$panel_variant_id, !colnames(cleaned_ld_matrix) %in% gwas_missing_snps$panel_variant_id]
  cleaned_gwas_list <- list(MAF = cleaned_gwas$frequency, 
                            snp = cleaned_gwas$panel_variant_id, 
                            position = cleaned_gwas$position, 
                            type = gwas_type, 
                            N = num_gwas_samples, 
                            pvalues=cleaned_gwas$pvalue,
                            LD = as.matrix(gwas_cleaned_ld_matrix), 
                            beta = cleaned_gwas$effect_size,
                            varbeta = cleaned_gwas$standard_error**2)
  #print(check_dataset(cleaned_gwas_list,req="LD"))
  # if there isn't a signal, further analysis should not be run
  if (min(cleaned_gwas_list$pvalues) < 1e-6){
    cat('\t signal found \n')
    return(cleaned_gwas_list)
  } else {
    cat('\t no signal found \n')
  }
}

run_coloc <- function(susie1, susie2) {
  # run susie finemapping on both datasets
  coloc <- coloc.susie(susie1,susie2)
  return(coloc$summary)
}

filter_gwas <- function(gwas, snp_list, ld_missing_snps){
  gwas_filtered <- subset(gwas, panel_variant_id %in% snp_list$variant_id & !(panel_variant_id %in% ld_missing_snps))
  return(gwas_filtered)
}

# filter eqtl data
filter_qtl <- function(qtl, snp_list, ld_missing_snps){
  qtl_filtered <- subset(qtl, variant_id %in% snp_list$variant_id & !(variant_id %in% ld_missing_snps))
  return(qtl_filtered)
}

split_qtl <- function(qtl, snp_list, ld_missing_snps, cleaned_ld_matrix, num_gtex_samples){
  qtl_filtered <- filter_qtl(qtl, snp_list, ld_missing_snps)
  # split eqtl to each of the egenes and clean up
  split_qtls <- split(qtl_filtered, qtl_filtered$phenotype_id)
  qtls_for_coloc <- list()
  for(i in seq_along(split_qtls)) {
    this_qtl_for_coloc <- clean_eqtl(split_qtls[[i]], cleaned_ld_matrix, num_gtex_samples)
    qtls_for_coloc[[i]] <- this_qtl_for_coloc
  }
  return(qtls_for_coloc)
}

gwas_from_path <- function(gwas_folder, gwas_meta, snp_list, ld_missing_snps, cleaned_ld_matrix){
  gwas_id <- gwas_meta$Tag
  cat(paste("working on gwas ", gwas_id, "\n"))
  num_gwas_samples <- gwas_meta$Sample_Size
  if (gwas_meta$Binary == 0) {
    gwas_type <- 'quant'
  } else {
    gwas_type <- 'cc'
  }
  gwas_path <- paste(gwas_folder, '/imputed_', gwas_id, '.txt.gz', sep = '')
  gwas <- fread(gwas_path)
  
  # filter gwas and clean gwas data
  gwas_filtered <- filter_gwas(gwas, snp_list, ld_missing_snps)
  gwas_for_coloc <- clean_gwas(gwas_filtered, cleaned_ld_matrix, gwas_type, num_gwas_samples)
  return(gwas_for_coloc)
}

get_empty_gwas_coloc <- function(){
  gwas_coloc_results <- data.frame(gwas_id = character(), 
                                   qtl_id = character(), 
                                   nsnps = numeric(),
                                   hit_gwas = character(),
                                   hit_qtl = character(),
                                   PP.H0.abf = numeric(),
                                   PP.H1.abf = numeric(),
                                   PP.H2.abf = numeric(),
                                   PP.H3.abf = numeric(),
                                   PP.H4.abf = numeric(),
                                   idx1 = numeric(),
                                   idx1 = numeric(),
                                   stringsAsFactors = FALSE)
  return(gwas_coloc_results)
} 


get_gwas_coloc_cluster <- function(eqtl, pcqtl, cluster_id, tissue_id, chr_id, snp_path_head, gwas_folder, gwas_meta, num_gtex_samples){
  
  # get snp list path from cluster id
  snp_list_path <-paste(snp_path_head, cluster_id, '.snp_list.txt', sep="")
  # load in snp list
  snp_list <- read_table(snp_list_path)
  cat("Total snps: ")
  cat(length(snp_list$variant_id))
  cat("\n")
  
  ld_matrix_path <- paste(snp_path_head, cluster_id, '.ld', sep="")
  # load in ld matrix
  # can't use fread here,  not sure why
  ld_matrix <- read.table(ld_matrix_path)
  rownames(ld_matrix) <- snp_list$variant_id
  colnames(ld_matrix) <- snp_list$variant_id
  
  # drop snps with missing values from LD 
  ld_missing_snps = get_ld_missing_snps(ld_matrix)
  cleaned_ld_matrix <- ld_matrix[!rownames(ld_matrix) %in% ld_missing_snps, !colnames(ld_matrix) %in% ld_missing_snps]
  cat("\t Number snps with ld missing: ")
  cat(length(ld_missing_snps))
  cat("\n")
  cleaned_ld_matrix <- ld_matrix[!rownames(ld_matrix) %in% ld_missing_snps, !colnames(ld_matrix) %in% ld_missing_snps]
  
  
  # filter eqtls and pcqtls
  eqtls_for_coloc <- split_qtl(eqtl, snp_list, ld_missing_snps, cleaned_ld_matrix, num_gtex_samples)
  pcqtls_for_coloc <- split_qtl(pcqtl, snp_list, ld_missing_snps, cleaned_ld_matrix, num_gtex_samples)
  cat(length(eqtls_for_coloc))
  cat(" qtl phenotypes \n")
  # conmbine them into one list
  qtls_for_coloc <- c(eqtls_for_coloc, pcqtls_for_coloc)
  # remove the nulls (from when there isn't a strong enough signal to finemap)
  qtls_for_coloc <- qtls_for_coloc[!sapply(qtls_for_coloc, is.null)]
  

  # initialize a null
  gwas_coloc_results <- get_empty_gwas_coloc()
  
  qtl_susies <- vector("list", length(qtls_for_coloc))
  if (length(qtls_for_coloc) != 0){
    cat(paste('running susie on', length(qtls_for_coloc), 'qtls \n'))
    for(i in seq_along(qtls_for_coloc)) {
      qtl_susie <- tryCatch({
        runsusie(qtls_for_coloc[[i]])  # Call runsusie function for each QTL
      }, error = function(e) {
        cat("An error occurred for QTL \n")
        cat("Error Message:", conditionMessage(e), "\n")  # Print the error message
        NULL  # Return NULL if an error occurs
      })
      qtl_susies[[i]] <- qtl_susie  # Store the result or NULL in the qtl_susies list
    }
    qtls_for_coloc <- qtls_for_coloc[!sapply(qtl_susies, is.null)]
    qtl_susies <- qtl_susies[!sapply(qtl_susies, is.null)]
    
    cat(paste('\t\t sucessful susie on', length(qtl_susies), 'qtls \n'))
    # write out the pairwise qtl colocs (we get these for free basically after running susie)
    qtl_coloc_results <- get_qtl_pairwise_coloc(qtls_for_coloc, qtl_susies)
    qtl_coloc_path <- paste(snp_path_head, cluster_id, '.qtl_coloc.txt', sep="")
    write.table(qtl_coloc_results, file=qtl_coloc_path, quote=FALSE, row.names=FALSE, sep='\t')
    
    if (length(qtl_susies) != 0){
      for (i in 1:nrow(gwas_meta)){
        # check if I've already done this gwas
        gwas_id <- gwas_meta[i]$Tag
          gwas_for_coloc <- gwas_from_path(gwas_folder, gwas_meta[i], snp_list, ld_missing_snps, cleaned_ld_matrix)
          if (!is.null(gwas_for_coloc)){
            gwas_susie <- tryCatch({
            runsusie(gwas_for_coloc)  # Call runsusie function for each QTL
          }, error = function(e) {
            cat("An error occurred for gwas:", gwas_id, "\n")
            cat("Error Message:", conditionMessage(e), "\n")  # Print the error message
            NULL  # Return NULL if an error occurs
          })
          if (!is.null(gwas_susie)){
          # coloc this gwas with each qtl
            for (i in seq_along(qtl_susies)) {
              this_qtl_susie <- qtl_susies[[i]]
              this_qtl_id <- qtls_for_coloc[[i]]$phenotype_id
              cat(paste("coloc for", gwas_id, "and", this_qtl_id, "\n"))
              this_coloc <- run_coloc(gwas_susie,this_qtl_susie)
              this_coloc$gwas_id <- gwas_id
              this_coloc$qtl_id <- this_qtl_id
              gwas_coloc_results <- rbind(gwas_coloc_results, this_coloc)
              cat(paste(dim(gwas_coloc_results), " total colocalizations run so far\n"))
            }
          }
        }
      }
    }
  }
  # rename and return results   
  gwas_coloc_results$gwas_cs_is <- gwas_coloc_results$idx1 
  gwas_coloc_results$qtl_cs_is <- gwas_coloc_results$idx2
  # write out this cluster coloc
  # this allows me to save my work a bit
  return(gwas_coloc_results)
}


make_snp_list <- function(eqtl, cluster_id, snp_path){
  # get a list of snps
  cluster_eqtl <- eqtl[eqtl$cluster_id == cluster_id, ]
  # write out snp list for use by plink to get ld
  snp_list <- cluster_eqtl[cluster_eqtl$phenotype_id == unique(cluster_eqtl$phenotype_id)[1], 'variant_id']
  write.csv(snp_list, file=snp_path, quote=FALSE, row.names=FALSE)
  
}


get_qtl_pairwise_coloc <- function(qtls_for_coloc, qtl_susies){
  qtl_coloc_results=list()
  qtl_coloc_results <- data.frame(qtl1_id = character(), 
                                  qtl2_id = character(), 
                                  nsnps = numeric(),
                                  hit1 = character(),
                                  hit2 = character(),
                                  PP.H0.abf = numeric(),
                                  PP.H1.abf = numeric(),
                                  PP.H2.abf = numeric(),
                                  PP.H3.abf = numeric(),
                                  PP.H4.abf = numeric(),
                                  idx1 = numeric(),
                                  idx1 = numeric(),
                                  stringsAsFactors = FALSE)
  
  qtl_pair_idxs <- combn(seq_along(qtl_susies), 2, simplty=TRUE)
  
  for (i in 1:ncol(qtl_pair_idxs)) {
    qtl1 <- qtls_for_coloc[[qtl_pair_idxs[1, i]]]
    qtl2 <- qtls_for_coloc[[qtl_pair_idxs[2, i]]]
    susie1 <- qtl_susies[[qtl_pair_idxs[1, i]]]
    susie2 <- qtl_susies[[qtl_pair_idxs[2, i]]]
    cat(paste("coloc for", qtl1$phenotype_id , "-", qtl2$phenotype_id), "\n")
    this_coloc <- run_coloc(susie1,susie2)
    if (!is.null(this_coloc)){
      this_coloc$qtl1_id <- qtl1$phenotype_id
      this_coloc$qtl2_id <- qtl2$phenotype_id
      qtl_coloc_results <- tryCatch({
          rbind(qtl_coloc_results, this_coloc)  
        }, error = function(e) {
          cat("An error occurred for QTL \n")
          cat("Error Message:", conditionMessage(e), "\n")  # Print the error message
          print(qtl_coloc_results)
          print(this_coloc)
          NULL  # Return NULL if an error occurs
      })
    }
  }
  return(qtl_coloc_results)
}


#############

# Create Argument Parser
parser <- ArgumentParser()
parser$add_argument("--eqtl_path", help="Input file path for eQTL pairs")
parser$add_argument("--pcqtl_path", help="Input file path for PCQTL pairs")
parser$add_argument("--gwas_meta", help="Input file path for GWAS metadata, with sample size and quant vs cc")
parser$add_argument("--gtex_meta", help="Input file path for GTEX metadata with sample size")
parser$add_argument("--chr_id", help="cluster id")
parser$add_argument("--tissue", help="tissue id")
parser$add_argument("--snp_path_head", help="directory file path for snp list and ld")
parser$add_argument("--genotype_stem", help="path to genotype_stem")
parser$add_argument("--gwas_folder", help="Folder path for GWAS data")
parser$add_argument("--output_path", help="Output file path for coloc results")




# Parse the Arguments
args <- parser$parse_args()

# Access the arguments
eqtl_path <- args$eqtl_path
pcqtl_path <- args$pcqtl_path
gwas_meta_path <- args$gwas_meta
gtex_meta_path <- args$gtex_meta
chr_id <- args$chr_id
tissue_id <- args$tissue
snp_path_head <- args$snp_path_head
genotype_stem <- args$genotype_stem
gwas_folder <- args$gwas_folder
output_path <- args$output_path



# Print or use the arguments as needed in the script
cat("eQTL Path:", eqtl_path, "\n")
cat("PCQTL Path:", pcqtl_path, "\n")
cat("GWAS Metadata:", gwas_meta_path, "\n")
cat("GTEX Metadata:", gtex_meta_path, "\n")
cat("LD and snplist dir:", snp_path_head, "\n")
cat("Chr ID:", chr_id, "\n")
cat("Tissue ID:", tissue_id, "\n")
cat("GWAS Folder:", gwas_folder, "\n")
cat("Output Path:", output_path, "\n")

# read in gtex meta
gtex_meta <- read.table(gtex_meta_path, sep='\t', header = T)
num_gtex_samples <- gtex_meta[gtex_meta$tissue_id == tissue_id, 'sample_size']

# read om gwas meta
gwas_meta <- fread(gwas_meta_path)
# subset for debugging
#gwas_meta <- gwas_meta[37:38]

#load in eqtl data
eqtl <- read_parquet(eqtl_path)
pcqtl <- read_parquet(pcqtl_path)

# get a list of clusters
eqtl$cluster_id <- sapply(eqtl$phenotype_id, function(x) unlist(strsplit(as.character(x), '_e_'))[1])
unique_cluster_ids <- unique(eqtl$cluster_id)

# subset for debug
#unique_cluster_ids <- unique_cluster_ids[1:2] 

# make the folder for output

# Check if the directory exists
if (!file.exists(snp_path_head)) {
  # If the directory does not exist, create it
  dir.create(snp_path_head, recursive = TRUE)
  cat("Directory created: ", snp_path_head, "\n")
} else {
  cat("Directory already exists: ", snp_path_head, "\n")
}

# initialize the output
gwas_all_cluster_coloc_results <- get_empty_gwas_coloc()


# run coloc for each cluster
for (cluster_id in unique_cluster_ids){
  cat("running on cluster")
  cat(cluster_id)
  cat("\n")
  snp_path <- paste(snp_path_head, cluster_id, '.snp_list.txt', sep="")
  # Check if the file exists
  if (file.exists(snp_path)) {
    cat("snp list already exists\n")
  } else {
    make_snp_list(eqtl, cluster_id, snp_path)
    cat("generated snp list\n")
  }
  
  ld_path <- paste(snp_path_head, cluster_id, sep="")
  ld_path_full <- paste(snp_path_head, cluster_id,'.ld', sep="")
  # Check if the file exists
  if (file.exists(ld_path_full)) {
    cat("ld matrix already exists\n")
  } else {
    plink_command <- sprintf("plink --bfile %s --extract %s --r2 square --out %s", genotype_stem, snp_path, ld_path)
    system(plink_command, intern=TRUE)
    cat("generated ld matrix\n")
  }
  
  cluster_coloc_path <-  paste(snp_path_head, cluster_id, '.gwas_coloc.txt', sep="")
  # Check if the file exists
  if (file.exists(cluster_coloc_path)) {
    cat("this cluster already colocalized\n")
    this_cluster_coloc <- read.table(cluster_coloc_path, sep='\t', header=T)
  } else {
    this_cluster_coloc <- get_gwas_coloc_cluster(eqtl, pcqtl, cluster_id, tissue_id, chr_id, snp_path_head, gwas_folder, gwas_meta, num_gtex_samples)
  }
  
  this_cluster_coloc <- get_gwas_coloc_cluster(eqtl, pcqtl, cluster_id, tissue_id, chr_id, snp_path_head, gwas_folder, gwas_meta, num_gtex_samples)
  cat(paste(dim(this_cluster_coloc), " total colocalizations run for", cluster_id, "\n"))
  gwas_all_cluster_coloc_results <- rbind(gwas_all_cluster_coloc_results, this_cluster_coloc) 
}

write.table(gwas_all_cluster_coloc_results, file=output_path, quote=FALSE, row.names=FALSE, sep='\t')


