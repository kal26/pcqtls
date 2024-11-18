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


source("workflow/scripts/coloc_functions.R")


# Create Argument Parser
parser <- ArgumentParser()
parser$add_argument("--eqtl_dir_path", help="path for eQTL pairs")
parser$add_argument("--pcqtl_dir_path", help="path for PCQTL pairs")
parser$add_argument("--gtex_meta", help="Input file path for GTEX metadata with sample size")
parser$add_argument("--tissue_id", help="tissue id")
parser$add_argument("--chr_id", help="chromosome id")
parser$add_argument("--ld_path_head", help="directory file path for snp list and ld")
parser$add_argument("--genotype_stem", help="path to genotype_stem")
parser$add_argument("--annotated_cluster_path", help="path to position annotated clusters")
parser$add_argument("--output_path", help="Output file path for coloc results")
parser$add_argument("--coloc_temp_path_head", help="tissue specific directory file path for susie and per chrom gwas")



# Parse the Arguments
args <- parser$parse_args()

# Access the arguments
eqtl_dir_path <- args$eqtl_dir_path
pcqtl_dir_path <- args$pcqtl_dir_path
gtex_meta_path <- args$gtex_meta
tissue_id <- args$tissue_id
chr_id <- args$chr_id
ld_path_head <- args$ld_path_head
genotype_stem <- args$genotype_stem
output_path <- args$output_path
annotated_cluster_path <- args$annotated_cluster_path
coloc_temp_path_head <- args$coloc_temp_path_head



# Print or use the arguments as needed in the script
cat("eQTL folder:", eqtl_dir_path, "\n")
cat("PCQTL folder:", pcqtl_dir_path, "\n")
cat("GTEX Metadata:", gtex_meta_path, "\n")
cat("LD and snplist dir:", ld_path_head, "\n")
cat("Tissue ID:", tissue_id, "\n")
cat("Chromosome ID:", chr_id, "\n")
cat("Output Path:", output_path, "\n")
cat("Cluster Path:", annotated_cluster_path, "\n")

cat("starting pair colocs \n")
start <- Sys.time()

# read in gtex meta to get sample size
gtex_meta <- read.table(gtex_meta_path, sep='\t', header = T)
num_gtex_samples <- gtex_meta[gtex_meta$tissue_id == tissue_id, 'sample_size']


# initalize empty results to add on to
all_cluster_coloc_results <- data.frame()
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
chr_id <- as.integer(sub("chr", "", chr_id))
chr_coloc_path <- paste(coloc_temp_path_head, tissue_id, ".v8.", 'chr_', chr_id, '.pairs_coloc.txt', sep="")
# check if partial results exist
# if (file.exists(chr_coloc_path)) {
#   cat("coloc results already exist up to chr", chr_id, "\n")
#   gwas_all_cluster_coloc_results <- read.table(chr_coloc_path, header=TRUE, sep='\t')
# } else {
cat("working on chr ", chr_id, "\n")
cluster_df_chr <- cluster_df[cluster_df$Chromosome == chr_id]
eqtl_chr <- get_eqtl_chr(eqtl_dir_path, chr_id, tissue_id)
pcqtl_chr <- get_pcqtl_chr(pcqtl_dir_path, chr_id, tissue_id)
# for each cluster, coloc the pc and eqtls
for (i in 1:nrow(cluster_df_chr)){
  this_cluster <- cluster_df_chr[i]
  cluster_id <- this_cluster$cluster_id
  cat("running on ", cluster_id, "\n")
  num_colocs <- num_colocs + 1
  cluster_coloc <- coloc_pairs_cluster(eqtl_chr, pcqtl_chr, cluster_id, ld_path_head, genotype_stem, num_gtex_samples, coloc_temp_path_head)
  # add to results list 
  cat("adding to results list\n")
  cat(ncol(all_cluster_coloc_results), " vs ", ncol(cluster_coloc), "\n")
  if(is.null(cluster_coloc)){
    cat("\tresult is null\n")
  } else {
    all_cluster_coloc_results <- bind_rows(all_cluster_coloc_results, cluster_coloc) 
    cat(num_colocs, " colocs so far \n")
  }
}
  # # finished with a chr, write out temp results
  # write.table(all_cluster_coloc_results, file=chr_coloc_path, quote=FALSE, row.names=FALSE, sep='\t')
  # cat("wrote out partial results, up to chr", chr_id, "\n")
#}

# write out (tissue_id.gwas_id)
write.table(all_cluster_coloc_results, file=output_path, quote=FALSE, row.names=FALSE, sep='\t')