
library(argparse)
library(coloc)
library(arrow)
library(tidyverse)
library(data.table)


########### functions #########

get_ld_missing_snps <- function(ld_matrix){
  # drop snps with missing values from LD for eqtl
  rows_with_missing <- rownames(ld_matrix)[which(rowSums(is.na(ld_matrix)) > 0)]
  cols_with_missing <- colnames(ld_matrix)[which(colSums(is.na(ld_matrix)) > 0)]
  ld_missing_snps = union(rows_with_missing, cols_with_missing)
  return(ld_missing_snps)
}

clean_eqtl <- function(eqtl_filtered, cleaned_ld_matrix){
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
    cat('signal found \n')
    return(cleaned_eqtl_list)
  } else{
    cat('no signal found \n')
  }    
}

clean_gwas <- function(gwas_filtered, cleaned_ld_matrix){
  # susie needs effect sizes, so we must also drop the snps with na for gwas effect
  gwas_missing_snps <- gwas_filtered[is.na(gwas_filtered$effect_size), 'panel_variant_id']
  print("Number snps with ld or gwas missing")
  print(length(gwas_missing_snps$panel_variant_id))
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
    cat('signal found \n')
    return(cleaned_gwas_list)
  } else {
    cat('no signal found \n')
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
filter_qtl <- function(qtl, snplist, ld_missing_snps){
  qtl_filtered <- subset(qtl, variant_id %in% snp_list$variant_id & !(variant_id %in% ld_missing_snps))
  return(qtl_filtered)
}

split_qtl <- function(qtl, snplist, ld_missing_snps, cleaned_ld_matrix){
  qtl_filtered <- filter_qtl(qtl, snplist, ld_missing_snps)
  # split eqtl to each of the egenes and clean up
  split_qtls <- split(qtl_filtered, qtl_filtered$phenotype_id)
  qtls_for_coloc <- list()
  for(i in seq_along(split_qtls)) {
    this_qtl_for_coloc <- clean_eqtl(split_qtls[[i]], cleaned_ld_matrix)
    qtls_for_coloc[[i]] <- this_qtl_for_coloc
  }
  return(qtls_for_coloc)
}



#############

# Create Argument Parser
parser <- ArgumentParser()
parser$add_argument("--eqtl_path", help="Input file path for eQTL pairs")
parser$add_argument("--pcqtl_path", help="Input file path for PCQTL pairs")
parser$add_argument("--gwas_meta", help="Input file path for GWAS metadata, with sample size and quant vs cc")
parser$add_argument("--gtex_meta", help="Input file path for GTEX metadata with sample size")
parser$add_argument("--cluster", help="cluster id")
parser$add_argument("--tissue", help="tissue id")
parser$add_argument("--snp_list_path", help="Input file path for snp list")
parser$add_argument("--ld_path", help="Input file path for LD matrix")
parser$add_argument("--gwas_folder", help="Folder path for GWAS data")
parser$add_argument("--output_path", help="Output file path for coloc results")


# Parse the Arguments
args <- parser$parse_args()

# Access the arguments
eqtl_path <- args$eqtl_path
pcqtl_path <- args$pcqtl_path
gwas_meta <- args$gwas_meta
gtex_meta <- args$gtex_meta
ld_path <- args$ld_path
cluster_id <- args$cluster
tissue_id <- args$tissue
snp_list_path <- args$snp_list_path
gwas_folder <- args$gwas_folder
output_path <- args$output_path



# Print or use the arguments as needed in the script
cat("eQTL Path:", eqtl_path, "\n")
cat("PCQTL Path:", pcqtl_path, "\n")
cat("GWAS Metadata:", gwas_meta, "\n")
cat("GTEX Metadata:", gwas_meta, "\n")
cat("LD Path:", ld_path, "\n")
cat("Cluster ID:", cluster_id, "\n")
cat("Tissue ID:", cluster_id, "\n")
cat("GWAS Folder:", gwas_folder, "\n")
cat("Output Path:", output_path, "\n")

# read in gtex meta
gtex_meta <- read.table(gtex_meta_path, sep='\t', header = T)
num_gtex_samples <- gtex_meta[gtex_meta$tissue_id == tissue_id, 'sample_size']

# read om gwas meta
gwas_meta <- fread(gwas_meta_path)

#load in eqtl data
eqtl <- read_parquet(eqtl_path)
pcqtl <- read_parquet(pcqtl_path)

# load in snp list
snp_list <- read_table(snp_list_path)
print("Total snps")
print(length(snp_list$variant_id))

# load in ld matrix
# can't use fread here,  not sure why
ld_matrix <- read.table(ld_path)
rownames(ld_matrix) <- snp_list$variant_id
colnames(ld_matrix) <- snp_list$variant_id


# drop snps with missing values from LD 
ld_missing_snps = get_ld_missing_snps(ld_matrix)
cleaned_ld_matrix <- ld_matrix[!rownames(ld_matrix) %in% ld_missing_snps, !colnames(ld_matrix) %in% ld_missing_snps]
print("Number snps with ld missing")
print(length(ld_missing_snps))
cleaned_ld_matrix <- ld_matrix[!rownames(ld_matrix) %in% ld_missing_snps, !colnames(ld_matrix) %in% ld_missing_snps]



# combine pc and eqtls 
eqtls_for_coloc <- split_qtl(eqtl, snplist, ld_missing_snps, cleaned_ld_matrix)
pcqtls_for_coloc <- split_qtl(pcqtl, snplist, ld_missing_snps, cleaned_ld_matrix)
qtls_for_coloc <- c(eqtls_for_coloc, pcqtls_for_coloc)
qtls_for_coloc <- qtls_for_coloc[!sapply(qtls_for_coloc, is.null)]


# intialize a null result
gwas_coloc_results <-  list()
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



if (length(qtls_for_coloc) != 0){
  print(paste('running susie on', length(qtls_for_coloc), 'qtls'))
  qtl_susies <- lapply(qtls_for_coloc, runsusie)
  
  for (i in 1:nrow(gwas_meta)){
    this_gwas <- gwas_meta[i]
    gwas_id <- this_gwas$Tag
    print(paste("working on gwas ", gwas_id))
    num_gwas_samples <- this_gwas$Sample_Size
    if (this_gwas$Binary == 0) {
      gwas_type <- 'quant'
    } else {
      gwas_type <- 'cc'
    }
    gwas_path <- paste(gwas_folder, '/imputed_', gwas_id, '.txt.gz', sep = '')
    gwas <- fread(gwas_path)
    
    # filter gwas and clean gwas data
    gwas_filtered <- filter_gwas(gwas, snp_list, ld_missing_snps)
    gwas_for_coloc <- clean_gwas(gwas_filtered, cleaned_ld_matrix)
    if (!is.null(gwas_for_coloc)){
      gwas_susie <- runsusie(gwas_for_coloc)
      print(summary(gwas_susie))
      # coloc this gwas with each qtl
      for (i in seq_along(qtl_susies)) {
        this_qtl_susie <- qtl_susies[[i]]
        this_qtl_id <- short_qtls_for_coloc[[i]]$phenotype_id
        print(paste("coloc for", gwas_id, "and", this_qtl_id))
        this_coloc <- run_coloc(gwas_susie,this_qtl_susie)
        this_coloc$gwas_id <- gwas_id
        this_coloc$qtl_id <- short_qtls_for_coloc[[i]]$phenotype_id
        gwas_coloc_results <- rbind(gwas_coloc_results, this_coloc)
        print(paste(dim(gwas_coloc_results), " total colocalizations run"))
      }
    }
  }
}


gwas_coloc_results$gwas_cs_is <- gwas_coloc_results$idx1 
gwas_coloc_results$qtl_cs_is <- gwas_coloc_results$idx2
write.table(gwas_coloc_results, file=output_path, row.names=FALSE, sep='\t')


