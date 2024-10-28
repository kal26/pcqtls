## following tutorial at https://hanruizhang.github.io/GWAS-eQTL-Colocalization/
## in r/4.2.2

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



# test args for debugging
#eqtl_dir_path <- 'output/proteincoding_main/control_eqtl/Thyroid'
#pcqtl_dir_path <- 'output/proteincoding_main/pcqtl/Thyroid'
#gwas_meta_path <- '/home/klawren/oak/pcqtls/data/references/gwas_metadata.txt'
#gtex_meta_path <- '/home/klawren/oak/pcqtls/data/references/gtex_sample_sizes.csv'
#tissue_id <- 'Thyroid'
#snp_path_head <- 'output/proteincoding_main/gwas_coloc/Thyroid/temp/'
#genotype_stem <- "data/processed/genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01"
#gwas_path <- '/oak/stanford/groups/smontgom/shared/gwas_summary_stats/barbeira_gtex_imputed/imputed_gwas_hg38_1.1/imputed_GIANT_WHR_Combined_EUR.txt.gz'
#output_path <- 'output/temp/test_coloc.txt'
#gwas_id <- 'GIANT_WHR_Combined_EUR'
#annotated_cluster_path <- 'output/proteincoding_main/annotations/Thyroid_clusters_annotated.csv'
#cluster_id <-'ENSG00000129151.8_ENSG00000176971.3'
#use_susie <- FALSE
#i <- 52
#chr_id <- 11


source("workflow/scripts/coloc_functions.R")

# parse args
# Create Argument Parser
parser <- ArgumentParser()
parser$add_argument("--eqtl_dir_path", help="folder for eQTL pairs")
parser$add_argument("--pcqtl_dir_path", help="folder for PCQTL pairs")
parser$add_argument("--gwas_meta", help="Input file path for GWAS metadata, with sample size and quant vs cc")
parser$add_argument("--gtex_meta", help="Input file path for GTEX metadata with sample size")
parser$add_argument("--tissue_id", help="tissue id")
parser$add_argument("--ld_path_head", help="directory file path for snp list and ld")
parser$add_argument("--coloc_temp_path_head", help="tissue specific directory file path for susie and per chrom gwas")
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
ld_path_head <- args$ld_path_head
coloc_temp_path_head <- args$coloc_temp_path_head
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
cat("LD and snplist dir:", ld_path_head, "\n")
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
gwas_all_cluster_coloc_results <- data.frame()
num_colocs <- 0


# Check if the directory exists
if (!file.exists(ld_path_head)) {
  # If the directory does not exist, create it
  dir.create(ld_path_head, recursive = TRUE)
  cat("Directory created: ", ld_path_head, "\n")
} else {
  cat("Directory already exists: ", ld_path_head, "\n")
}

# Check if the directory exists
if (!file.exists(coloc_temp_path_head)) {
  # If the directory does not exist, create it
  dir.create(coloc_temp_path_head, recursive = TRUE)
  cat("Directory created: ", coloc_temp_path_head, "\n")
} else {
  cat("Directory already exists: ", coloc_temp_path_head, "\n")
}

# load in clusters
cluster_df <- fread(annotated_cluster_path)

# for each chromosome
for (chr_id in 1:22){
  chr_coloc_path <- paste(coloc_temp_path_head, tissue_id, ".v8.", gwas_id, '.susie_', use_susie,'chr_', chr_id, '.gwas_coloc.txt', sep="")
  # check if partial results exist
  if (file.exists(chr_coloc_path)) {
    cat("coloc results already exist up to chr", chr_id, "\n")
    gwas_all_cluster_coloc_results <- read.table(chr_coloc_path, header=TRUE, sep='\t')
  } else {
    cat("working on chr ", chr_id, "\n")
    cluster_df_chr <- cluster_df[cluster_df$Chromosome == chr_id]
    gwas_chr <- gwas_with_meta$gwas_data[gwas_with_meta$gwas_data$chromosome == paste('chr', chr_id, sep="")] 
    pcqtl_chr <- NULL
    eqtl_chr <- NULL
    # for each cluster, check if the gwas has a signal
    for (i in 1:nrow(cluster_df_chr)){
      this_cluster <- cluster_df_chr[i]
      cluster_id <- this_cluster$cluster_id
      cat("running on ", cluster_id, "\n")
      # if it does, load in eqtl and pcqtl and colocalize 
      if (check_gwas_cluster(gwas_chr, this_cluster)){
        num_colocs <- num_colocs + 1
        cat("\tpossible coloc for", this_cluster$cluster_id, "\n")
        cat(num_colocs, " colocs so far \n")
        # load in eqtl if I haven't already
        if(is.null(eqtl_chr)){
          eqtl_chr <- get_eqtl_chr(eqtl_dir_path, chr_id, tissue_id)
        }
        # load in pcqtl if I haven't already
        if(is.null(pcqtl_chr)){
          pcqtl_chr <- get_pcqtl_chr(pcqtl_dir_path, chr_id, tissue_id)
        }
        gwas_cluster_coloc <- coloc_gwas_cluster(gwas_with_meta, eqtl_chr, pcqtl_chr, cluster_id, ld_path_head, coloc_temp_path_head, genotype_stem, num_gtex_samples, use_susie=use_susie)
        cat("finished cluster coloc")
        # make sure types match (not sure why this needs to be here?)
        gwas_all_cluster_coloc_results$gwas_id <- as.character(gwas_all_cluster_coloc_results$gwas_id)
        gwas_cluster_coloc$gwas_id <- as.character(gwas_cluster_coloc$gwas_id)
        gwas_all_cluster_coloc_results$qtl_id <- as.character(gwas_all_cluster_coloc_results$qtl_id)
        gwas_cluster_coloc$qtl_id <- as.character(gwas_cluster_coloc$qtl_id)
        # add to results list 
        gwas_all_cluster_coloc_results <- bind_rows(gwas_all_cluster_coloc_results, gwas_cluster_coloc) 
      }
    }
    # finished with a chr, write out temp results
    write.table(gwas_all_cluster_coloc_results, file=chr_coloc_path, quote=FALSE, row.names=FALSE, sep='\t')
    cat("wrote out partial results, up to chr", chr_id, "\n")
  }
}

# write out (tissue_id.gwas_id)
write.table(gwas_all_cluster_coloc_results, file=output_path, quote=FALSE, row.names=FALSE, sep='\t')
cat("Finished. Wrote out all results.")
