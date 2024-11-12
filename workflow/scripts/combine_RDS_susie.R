setwd('/home/klawren/oak/pcqtls/')

library(dplyr)
library(argparse)

source("workflow/scripts/coloc_functions.R")

# Create Argument Parser
parser <- ArgumentParser()
parser$add_argument("--qtl_dir_path", help="folder for QTL pairs")
parser$add_argument("--tissue_id", help="tissue id")
parser$add_argument("--coloc_temp_path_head", help="tissue specific directory file path for susie and per chrom gwas")
parser$add_argument("--output_path", help="Output file path for coloc pair results")
parser$add_argument("--qtl_type", help="eqtl or pcqtl")

# Parse the Arguments
args <- parser$parse_args()

# assign values
coloc_temp_path_head <- args$coloc_temp_path_head
qtl_dir_path <- args$qtl_dir_path
tissue_id <- args$tissue_id
qtl_type <- args$qtl_type
output_path <- args$output_path


combined_results <- data.frame()  # To store combined results

for (chr_id in 1:22){
    cat("working on chr ", chr_id, "\n")
    # load in nominal pc and e resutls to get phenotype ids
    if (qtl_type == 'eqtl'){
        qtl_chr <- get_eqtl_chr(qtl_dir_path, chr_id, tissue_id)
    } else {
        qtl_chr <- get_pcqtl_chr(qtl_dir_path, chr_id, tissue_id)
    }
    phenotype_ids <- unique(qtl_chr$phenotype_id)
    # Generate file paths for the corresponding .susie.rds files
    for (phenotype_id in phenotype_ids) {
        if (nchar(phenotype_id) > 150){
            out_qtl_id <- get_short_qtl_id(phenotype_id)
        } else {
            out_qtl_id <- phenotype_id
        }
        susie_rds_path <- paste(coloc_temp_path_head, out_qtl_id, '.susie.rds', sep="")
        # Check if the file exists
        if (file.exists(susie_rds_path)) {
            # Extract data from the RDS file and combine it
            rds_data <- extract_data_from_rds(susie_rds_path)
            combined_results <- bind_rows(combined_results, rds_data)
        } else {
            # expected that some files do not exist 
            # there were no signifigant QTL signals for this phenotype
            message(paste("File does not exist:", susie_rds_path))
        }
    }
}

write.table(combined_results, file = output_path, sep = "\t", quote = FALSE, row.names = TRUE)
