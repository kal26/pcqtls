## following tutorial at https://hanruizhang.github.io/GWAS-eQTL-Colocalization/
## in r/4.2.2

setwd('/home/klawren/oak/pcqtls/')
library(coloc)
library(arrow)
library(tidyverse)
library(data.table)
library(Rfast)
library(argparse)
library(tidyverse)
library(susieR)


# test args for debugging
eqtl_dir_path <- 'output/proteincoding_main/control_eqtl/Thyroid'
pcqtl_dir_path <- 'output/proteincoding_main/pcqtl/Thyroid'
gwas_meta_path <- '/home/klawren/oak/pcqtls/data/references/gwas_metadata.txt'
gtex_meta_path <- '/home/klawren/oak/pcqtls/data/references/gtex_sample_sizes.csv'
tissue_id <- 'Thyroid'
snp_path_head <- 'output/proteincoding_main/gwas_coloc/Thyroid/temp/'
genotype_stem <- "data/processed/genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01"
gwas_path <- '/oak/stanford/groups/smontgom/shared/gwas_summary_stats/barbeira_gtex_imputed/imputed_gwas_hg38_1.1/imputed_UKB_20016_Fluid_intelligence_score.txt.gz'
output_path <- 'output/temp/test_coloc.txt'
gwas_id <- 'UKB_20016_Fluid_intelligence_score'
annotated_cluster_path <- 'output/proteincoding_main/annotations/Thyroid_clusters_annotated.csv'
cluster_id <-'ENSG00000161016.17_ENSG00000196378.11'

####### funtions #########

split_qtl_for_coloc <- function(cluster_qtl, ld_snp_set, cleaned_ld_matrix, num_gtex_samples){
  qtl_filtered <- filter_qtl(cluster_qtl, ld_snp_set)
  # split eqtl to each of the egenes and clean up
  phenotype_ids <-  unique(qtl_filtered$phenotype_id)
  qtls_for_coloc <- list()
  for (i in 1:length(phenotype_ids)){
    phenotype_qtl <- qtl_filtered[qtl_filtered$phenotype_id == phenotype_ids[i], ]
    cat(paste("\t\tlooking for signals in",phenotype_ids[i]), "\n")
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
  # if there isn't a signal, further analysis should not be run
  if (min(qtl_filtered$pval_nominal) < 1e-6){
    cat('\t signal found \n')
    # make some columns to play well with coloc
    qtl_filtered <- qtl_filtered %>%
      mutate(position = as.integer(str_split(variant_id, "_") %>% sapply(pluck, 2)))
    cleaned_qtl_list <- list(beta = qtl_filtered$slope, 
                             varbeta = qtl_filtered$slope_se**2, 
                             snp = qtl_filtered$variant_id, 
                             position = qtl_filtered$position, 
                             pvalues = qtl_filtered$pval_nominal,
                             type = "quant", 
                             N = num_gtex_samples, 
                             MAF = qtl_filtered$af,
                             sdY = 1,
                             phenotype_id = qtl_filtered$phenotype_id[1],
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
  gwas_missing_snps <- gwas_filtered[is.na(gwas_filtered$effect_size) | is.na(gwas_filtered$effect_size), 'panel_variant_id']
  cat("\t Number snps with ld or gwas missing: ", length(gwas_missing_snps$panel_variant_id), "\n")
  cleaned_gwas <- gwas_filtered[!gwas_filtered$panel_variant_id %in% gwas_missing_snps$panel_variant_id, ]
  if(min(cleaned_gwas$pvalue) < 1e-6){
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
    cat('\t signal found \n')
    return(cleaned_gwas_list)
  } else {
    cat('\t no signal found, this is odd becuase I already checked and there should be a signal\n')
    return(NULL)
  }
}

get_empty_gwas_coloc <- function(use_susie){
  if(use_susie){
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
  } else {
    gwas_coloc_results <- data.frame(gwas_id = character(), 
                                     qtl_id = character(), 
                                     nsnps = numeric(),
                                     PP.H0.abf = numeric(),
                                     PP.H1.abf = numeric(),
                                     PP.H2.abf = numeric(),
                                     PP.H3.abf = numeric(),
                                     PP.H4.abf = numeric(),
                                     stringsAsFactors = FALSE)
  }

  return(gwas_coloc_results)
} 

get_empty_qtl_coloc <- function(){
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
  return(qtl_coloc_results)                                
}

runsusie_errorcatch <- function(dataset){
  cat("running susie\n")
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

coloc_gwas_cluster <- function(gwas_with_meta, eqtl_chr, pcqtl_chr, cluster_id, snp_path_head, genotype_stem, num_gtex_samples, use_susie = FALSE){
  # subset eqtl and pcqtl to this cluster
  cluster_eqtl <-eqtl_chr[eqtl_chr$cluster_id == cluster_id, ]
  cluster_pcqtl <- pcqtl_chr[pcqtl_chr$cluster_id == cluster_id, ]
  # get snp list and ld matrix
  snp_list <- get_snp_list(cluster_eqtl, snp_path_head, cluster_id)
  cleaned_ld_matrix <- get_ld(snp_path_head, cluster_id, snp_list, genotype_stem)
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
  
  #DEBUG (trying to work out saving susie)
  #gwas_susie[[1]]
  #saveRDS(gwas_susie, file = "output/temp/gwas_susie_data.rds")
  # end debug
  
  # clean the gwas data into the right format for susie
  gwas_for_coloc <- get_gwas_for_coloc(gwas_with_meta, ld_snp_set, snp_list, cleaned_ld_matrix)
  # make the gwas susie
  
  if(use_susie){
    # run susie on each qtl phenotype
    qtl_susies <- lapply(qtls_for_coloc, runsusie_errorcatch)
    qtls_for_coloc <- qtls_for_coloc[!sapply(qtl_susies, is.null)]
    qtl_susies <- qtl_susies[!sapply(qtl_susies, is.null)]
    if (length(qtl_susies)==0){
      cat("no qtl signals finemapped in this cluster\n")
      return(NULL)
    }
    gwas_susie <- runsusie_errorcatch(gwas_for_coloc)
  }
  # colocalize each qtl susie to the gwas susie
  gwas_coloc_results <- get_empty_gwas_coloc(use_susie)
  for (i in 1:length(qtls_for_coloc)){
    this_qtl_id <- qtls_for_coloc[[i]]$phenotype_id
    cat(paste("\t\t\tcoloc for", gwas_with_meta$gwas_id, "and", this_qtl_id, "\n"))
    if(use_susie){
      cat("\t\t using susie to coloc \n")
      this_qtl_susie <- qtl_susies[[i]]
      this_coloc <- run_coloc(gwas_susie,this_qtl_susie)
      if(!is.null(this_coloc)){
        this_coloc$gwas_id <- gwas_id
        this_coloc$qtl_id <- this_qtl_id
      }
    } else{
      cat("\t\t not using susie to coloc \n")
      this_coloc <- coloc.abf(gwas_for_coloc, qtls_for_coloc[[i]])$summary
      this_coloc$gwas_id <- gwas_id
      this_coloc$qtl_id <- this_qtl_id
    }
    gwas_coloc_results <- rbind(gwas_coloc_results, this_coloc)
  }
  # return the results
  gwas_coloc_results$gwas_cs_is <- gwas_coloc_results$idx1 
  gwas_coloc_results$qtl_cs_is <- gwas_coloc_results$idx2
  return(gwas_coloc_results)
}


run_coloc <- function(susie1, susie2) {
  # run susie finemapping on both datasets
  coloc <- coloc.susie(susie1,susie2)
  return(coloc$summary)
}


check_gwas_cluster <- function(gwas_chr, this_cluster){
  # return true if the gwas has a signal in this cluster, false otherwise
  gwas_cluster <- gwas_chr[gwas_chr$position > this_cluster$start - 1e6, ] 
  gwas_cluster <- gwas_cluster[gwas_cluster$position < this_cluster$end + 1e6, ]
  return(sum(gwas_cluster$pvalue < 1e-6)>0)
}

get_ld <- function(snp_path_head, cluster_id, snp_list, genotype_stem){
  # check if ld already exists
  ld_matrix_path <- paste(snp_path_head, cluster_id, '.ld', sep="")
  if (file.exists(ld_matrix_path)) {
    cat("ld matrix already exists\n")
  } else {    
    cat("ld matrix does not already exist\n")
    # get ld if not
    ld_plink_path <- paste(snp_path_head, cluster_id, sep="")
    snp_path <- paste(snp_path_head, cluster_id, '.snp_list.txt', sep="")
    plink_command <- sprintf("plink --bfile %s --extract %s --r square --out %s", genotype_stem, snp_path, ld_plink_path)
    cat(plink_command) 
    system(plink_command, intern=TRUE)
    cat("generated ld matrix\n")
  }
  # load in ld and return 
  ld_matrix <- read.table(ld_matrix_path)
  length(snp_list$variant_id)
  nrow(ld_matrix)
  rownames(ld_matrix) <- snp_list$variant_id
  colnames(ld_matrix) <- snp_list$variant_id
  
  # drop snps with missing values from LD
  ld_missing_snps = get_ld_missing_snps(ld_matrix)
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

get_snp_list <- function(cluster_eqtl, snp_path_head, cluster_id){
  snp_path <- paste(snp_path_head, cluster_id, '.snp_list.txt', sep="")
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
  eqtl <- read_parquet(eqtl_path)
  eqtl$cluster_id <- sapply(eqtl$phenotype_id, function(x) unlist(strsplit(as.character(x), '_e_'))[1])
  return(eqtl)
}

get_pcqtl_chr <- function(pcqtl_dir_path, chr_id, tissue_id){
  pcqtl_path <- paste(pcqtl_dir_path, "/", tissue_id, '.v8.pcs.cis_qtl_pairs.chr', chr_id, '.parquet', sep="")
  pcqtl <- read_parquet(pcqtl_path)
  pcqtl$cluster_id <- sapply(pcqtl$phenotype_id, function(x) unlist(strsplit(as.character(x), '_pc'))[1])
  return(pcqtl)
}

###### main ######


# parse args
# Create Argument Parser
parser <- ArgumentParser()
parser$add_argument("--eqtl_dir_path", help="folder for eQTL pairs")
parser$add_argument("--pcqtl_dir_path", help="folder for PCQTL pairs")
parser$add_argument("--gwas_meta", help="Input file path for GWAS metadata, with sample size and quant vs cc")
parser$add_argument("--gtex_meta", help="Input file path for GTEX metadata with sample size")
parser$add_argument("--tissue_id", help="tissue id")
parser$add_argument("--snp_path_head", help="directory file path for snp list and ld")
parser$add_argument("--genotype_stem", help="path to genotype_stem")
parser$add_argument("--gwas_path", help="path to GWAS data")
parser$add_argument("--gwas_id", help="GWAS id")
parser$add_argument("--annotated_cluster_path", help="path to position annotated clusters")
parser$add_argument("--output_path", help="Output file path for coloc results")
parser$add_argument("--use_susie", help="use susie for finemapping?", default=FALSE)



# Parse the Arguments
args <- parser$parse_args()

# Access the arguments
eqtl_dir_path <- args$eqtl_dir_path
pcqtl_dir_path <- args$pcqtl_dir_path
gwas_meta_path <- args$gwas_meta
gtex_meta_path <- args$gtex_meta
tissue_id <- args$tissue_id
snp_path_head <- args$snp_path_head
genotype_stem <- args$genotype_stem
gwas_path <- args$gwas_path
output_path <- args$output_path
gwas_id <- args$gwas_id
annotated_cluster_path <- args$annotated_cluster_path
use_susie <- as.logical(args$use_susie)

# Print or use the arguments as needed in the script
cat("eQTL folder:", eqtl_dir_path, "\n")
cat("PCQTL folder:", pcqtl_dir_path, "\n")
cat("GWAS Metadata:", gwas_meta_path, "\n")
cat("GTEX Metadata:", gtex_meta_path, "\n")
cat("LD and snplist dir:", snp_path_head, "\n")
cat("Tissue ID:", tissue_id, "\n")
cat("GWAS path:", gwas_path, "\n")
cat("Output Path:", output_path, "\n")
cat("Cluster Path:", annotated_cluster_path, "\n")
cat("use susie", use_susie, "\n")

cat("starting colocs for ", gwas_id, "\n")
start <- Sys.time()

# read in gtex meta to get sample size
gtex_meta <- read.table(gtex_meta_path, sep='\t', header = T)
num_gtex_samples <- gtex_meta[gtex_meta$tissue_id == tissue_id, 'sample_size']

# load in one gwas (specified)
gwas_meta_df <- fread(gwas_meta_path)
gwas_with_meta <- load_gwas_from_path(gwas_path, gwas_meta_df, gwas_id)
gwas_data <- gwas_with_meta$gwas_data
cat("total gwas signals", sum(gwas_data$pvalue < 1e-6), "\n")

# initalize empty results to add on to
gwas_all_cluster_coloc_results <- get_empty_gwas_coloc(use_susie)
num_colocs <- 0


# Check if the directory exists
if (!file.exists(snp_path_head)) {
  # If the directory does not exist, create it
  dir.create(snp_path_head, recursive = TRUE)
  cat("Directory created: ", snp_path_head, "\n")
} else {
  cat("Directory already exists: ", snp_path_head, "\n")
}

# load in clusters
cluster_df <- fread(annotated_cluster_path)

# for each chromosome
for (chr_id in 1:22){
  #chr_id <- 8#debug
  cat("working on chr ", chr_id, "\n")
  cluster_df_chr <- cluster_df[cluster_df$Chromosome == chr_id]
  gwas_chr <- gwas_with_meta$gwas_data[gwas_with_meta$gwas_data$chromosome == paste('chr', chr_id, sep="")] 
  pcqtl_chr <- NULL
  eqtl_chr <- NULL
  # for each cluster, check if the gwas has a signal
  for (i in nrow(cluster_df_chr)){
    #i <- 34 #debug
    this_cluster <- cluster_df_chr[i]
    cluster_id <- this_cluster$cluster_id
    # if it does, load in eqtl and pcqtl and colocalize 
    if (check_gwas_cluster(gwas_chr, this_cluster)){
      num_colocs <- num_colocs + 1
      cat("possible coloc for", this_cluster$cluster_id, "\n")
      cat(num_colocs, " colocs so far \n")
      # load in eqtl if I haven't already
      if(is.null(eqtl_chr)){
        eqtl_chr <- get_eqtl_chr(eqtl_dir_path, chr_id, tissue_id)
      }
      # load in pcqtl if I haven't already
      if(is.null(pcqtl_chr)){
        pcqtl_chr <- get_pcqtl_chr(pcqtl_dir_path, chr_id, tissue_id)
      }
      gwas_cluster_coloc <- coloc_gwas_cluster(gwas_with_meta, eqtl_chr, pcqtl_chr, cluster_id, snp_path_head, genotype_stem, num_gtex_samples, use_susie=use_susie)
      # add to results list 
      gwas_all_cluster_coloc_results <- rbind(gwas_all_cluster_coloc_results, gwas_cluster_coloc) 
    }
  }
}

# write out (tissue_id.gwas_id)
write.table(gwas_all_cluster_coloc_results, file=output_path, quote=FALSE, row.names=FALSE, sep='\t')
