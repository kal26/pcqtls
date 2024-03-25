#!/usr/bin/env Rscript

#-----------------------------------------------------------#
# CALCULATE CROSS-CORRELATIONS FOR USER-SPECIFIED:          #
#  - tissue      [character]                                #
#  - covars_yn   [character]  (yes_covars, no_covars)       #
#  - n_peers     [integer]    (0:60)                        #
#  - n_pcs       [integer]    (0:5)                         #
#  - chr_txt     [character]  (comma-separated,             #
#                              colon-separated, or int;     #
#                              e.g. 1:3,4,5 OR 1,2 OR 21    #
#                                                           #
#                                                           #
# FOR TESTING:                                              #
# tissue='Muscle_Skeletal'                                  #
# covars_yn = 'yes_covars'                                  #
# n_peers = 60                                              #
# n_pcs = 0                                                 #
# chr_txt = '22'                                            #
#-----------------------------------------------------------#

#---------------------------------------------------------------------------------------#
# Rscript calculate_coexpression_within_chr.R Muscle_Skeletal yes_covars 60 0 20:22     #
#---------------------------------------------------------------------------------------#

rm(list=ls())

startTime = Sys.time()

library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(corrplot)

#-------------------------------#
# ARGUMENTS FROM COMMAND LINE   #
#-------------------------------#
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
  
} else if (length(args) > 1) {
  tissue = args[1]
  covars_yn = args[2]
  n_peers = as.numeric(args[3])
  n_pcs = as.numeric(args[4])
  chr_txt = args[5]
  gene_annotations_path = args[6]
  normalized_expression_path = args[7]
  covariates_path = args[8]
  gencode_path = args[9]
  correlation_outpath = args[10]
  covariates_outpath = args[11]
  residualized_expression_outpath = args[12]
}



#-------------#
# FUNCTIONS   #
#-------------#
clean_tissue_rnaseq = function(rnaseq_norm_df){
  
  # Clean up gene expression matrix
  rownames(rnaseq_norm_df) = rnaseq_norm_df$gene_id
  
  # Order by transcript distance
  rnaseq_norm_df = rnaseq_norm_df[order(rnaseq_norm_df$start),]
  
  # Keep IDs
  ids = names(rnaseq_norm_df %>% select(starts_with("GTEX")))
  rnaseq_norm_df = rnaseq_norm_df[,ids]
  
  rnaseq_norm_df = data.frame(t(rnaseq_norm_df))
  
  return(rnaseq_norm_df)
  
}

make_ann_df = function(rnaseq_norm_df){
  
  ann_df = rnaseq_norm_df %>% select(X.chr:gene_id)
  
  names(ann_df) = c('Chromosome', 'Start', 'End', 'Gene_id')
  ann_df$Mean_position = (ann_df$Start + ann_df$End)/2
  
  ann_df$Distance_rank = seq(1, nrow(ann_df))
  
  return(ann_df)
}

clean_covars_df = function(covars_df){
  
  # Make ID the rownames 
  rownames(covars_df) = covars_df$ID
  covars_df = covars_df[,-1]
  
  # Mix row and col names 
  covars_df = data.frame(t(covars_df))
  
  return(covars_df)
}

subset_rnaseq_chr_n = function(gene_ann, rnaseq_df, chr_n){
  
  # Identify genes in chromosome n 
  genes_chr = (gene_ann %>%
                 filter(Chromosome == str_interp("chr${chr_n}")))$Gene_id
  
  # Replace - with . and keep genes for which RNA-seq is available. 
  genes_chr = genes_chr[genes_chr %in% names(rnaseq_df)]
  
  n_genes_chr = length(genes_chr)
  
  # Subset genes from RNA-seq data 
  rnaseq_chr_df = rnaseq_df %>%
    select(all_of(genes_chr))
  
  print(str_interp("There are ${n_genes_chr} gene IDs in chromosome ${chr_n}."))
  
  return(rnaseq_chr_df)
  
}

cor_mtest = function(mat, ...) {
  
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  
  return(p.mat)
  
}

calculate_corrs_and_p = function(tissue_chr_df, gene_ann, chr_n){
  
  # Subset genes in chromosome n
  gene_ann_chr = gene_ann %>% 
    filter(Chromosome == str_interp("chr${chr_n}")) %>%
    filter(Gene_id %in% names(tissue_chr_df))
  
  # Order by distance
  gene_ann_chr = gene_ann_chr[order(gene_ann_chr$Start),]
  tissue_chr_df = tissue_chr_df[,gene_ann_chr$Gene_id]
  
  n_transcripts = ncol(tissue_chr_df)
  n_individuals = nrow(tissue_chr_df)
  
  print(str_interp("Calculating correlations across ${n_transcripts} transcripts and ${n_individuals} individuals."))
  
  # Make numeric 
  tissue_chr_df = sapply(tissue_chr_df, as.numeric)
  
  # Calculate Spearman correlations
  cor_r_df = cor(tissue_chr_df, method = 'spearman', use='complete.obs')
  cor_p_df = cor_mtest(tissue_chr_df)
  
  results_list = list(cor_r_df, cor_p_df)
  
  return(results_list)
  
}

transform_mat_to_df = function(cor_r_m, cor_p_m){
  
  cor_r_df = data.frame(row=rownames(cor_r_m)[row(cor_r_m)[upper.tri(cor_r_m)]], 
                        col=colnames(cor_r_m)[col(cor_r_m)[upper.tri(cor_r_m)]], 
                        corr=cor_r_m[upper.tri(cor_r_m)])
  
  cor_p_df = data.frame(row=rownames(cor_p_m)[row(cor_p_m)[upper.tri(cor_p_m)]], 
                        col=colnames(cor_p_m)[col(cor_p_m)[upper.tri(cor_p_m)]], 
                        p=cor_p_m[upper.tri(cor_p_m)])
  
  cor_df = merge(cor_r_df, cor_p_df, by=c('row', 'col'))
  
  names(cor_df) = c('Gene1', 'Gene2', 'Corr_coefficient', 'P')
  
  return(cor_df)
  
}

keep_sig_and_add_gene_distance = function(cor_df, gene_ann, ref_df, p_thresh){
  
  cor_df = cor_df %>%
    filter(P < p_thresh)
  
  for (i in 1:nrow(cor_df)){
    
    # Annotate gene 1
    gene1_info = ref_df[which(ref_df$Transcript_name == cor_df[i,'Gene1']), ]
    
    gene1_name = gene1_info$Gene_name
    gene1_start = gene1_info$Start
    gene1_end = gene1_info$End
    gene1_distance_rank = gene1_info$Distance_rank
    gene1_distance_expressed_rank = gene_ann[which(gene_ann$Gene_id == cor_df[i,'Gene1']), 'Distance_rank']
    
    # Annotate gene 2
    gene2_info = ref_df[which(ref_df$Transcript_name == cor_df[i,'Gene2']), ]
    
    gene2_name = gene2_info$Gene_name
    gene2_start = gene2_info$Start
    gene2_end = gene2_info$End
    gene2_distance_rank = gene2_info$Distance_rank
    gene2_distance_expressed_rank = gene_ann[which(gene_ann$Gene_id == cor_df[i,'Gene2']), 'Distance_rank']
    
    cor_df[i, c('Gene1_name', 'G1_start', 'G1_end', 'G1_distance_rank',
                'G1_distance_exp_rank', 
                'Gene2_name', 'G2_start', 'G2_end', 'G2_distance_rank',
                'G2_distance_exp_rank')] = 
      c(gene1_name, gene1_start, gene1_end, gene1_distance_rank, gene1_distance_expressed_rank,
        gene2_name, gene2_start, gene2_end, gene2_distance_rank, gene2_distance_expressed_rank)
    
  }
  
  return(cor_df)
  
}

make_covariates_txt = function(covariates, n_pcs, n_peers){
  
  # Covariates 
  covars_list = c()
  
  # Number of PCs 
  if (n_pcs >= 1){
    for (i in 1:n_pcs){
      covars_list = c(covars_list, str_interp("PC${i}"))
    }
  }
  
  # Number of PEER factors 
  if (n_peers >= 1){
    for (i in 1:n_peers){
      covars_list = c(covars_list, str_interp("InferredCov${i}"))
    }
  }
  
  # Add standard covariates 
  covars_list = c(covars_list, covariates)
  
  # Convert covariates to text 
  covars_txt = paste(covars_list, collapse=' + ')
  
  return(covars_txt)
}

calculate_residuals = function(rnaseq_tissue_chr_df, covars_df,
                               covariates, n_pcs, n_peers){
  
  # Transcripts 
  transcripts = names(rnaseq_tissue_chr_df)
  
  # Merge RNA-seq with covariates 
  df_with_covar = merge(rnaseq_tissue_chr_df, covars_df, by=0)
  ids = df_with_covar$Row.names
  df_with_covar = df_with_covar[,-1]
  
  df_with_covar = as.data.frame(sapply(df_with_covar, as.numeric))
  rownames(df_with_covar) = ids
  
  # Make covariates
  covars_txt = make_covariates_txt(covariates, n_pcs, n_peers)

  # write out covariates
  write.csv(covars_txt, covariates_outpath,
            row.names=FALSE) 

  # Residualize transcripts 
  if (covars_txt != ''){
    
    for (t in transcripts){
      
      formula = str_interp("${t} ~ ${covars_txt}")
      
      results = lm(formula, data = df_with_covar)
      rnaseq_tissue_chr_df[,t] = results$residuals
      
    }
    
  }
  
  return(rnaseq_tissue_chr_df)
  
}

process_chr_txt = function(chr_txt){
  
  # Only comas 
  if (grepl(",", chr_txt) & !(grepl(":", chr_txt))){
    chr_list_clean = as.numeric(unlist(strsplit(chr_txt, ",")))
    
    # Only colon
  } else if (grepl(":", chr_txt) & !(grepl(",", chr_txt))){
    
    chr_list_clean = seq(unlist(strsplit(chr_txt, ":"))[1], 
                         unlist(strsplit(chr_txt, ":"))[2])
    
    # Comas + colon
  } else if (grepl(",", chr_txt) & (grepl(":", chr_txt))){
    chr_list = unlist(strsplit(chr_txt, ","))
    
    chr_list_clean = c()
    
    for (i in chr_list){
      if (grepl(":", i)){
        chr_list = seq(unlist(strsplit(i, ":"))[1], 
                       unlist(strsplit(i, ":"))[2])
        chr_list_clean = c(chr_list_clean, chr_list)
      } else {
        chr_list_clean = c(chr_list_clean, i)
      }
    }
    
    chr_list_clean = as.numeric(chr_list_clean)
    
  } else {
    chr_list_clean = as.numeric(chr_txt)
  }
  
  chr_list_clean
}



#-------------#
# LOAD DATA   #
#-------------#

# Protein-coding reference 
ref_proteincoding_df = read.csv(gene_annotations_path)

# Tissue RNA-seq data
rnaseq_norm_df = read.delim2(normalized_expression)

# Keep only protein-coding genes
rnaseq_norm_df = rnaseq_norm_df %>%
  filter(gene_id %in% ref_proteincoding_df$Transcript_name)


# Extract gene annotatons
ann_df = make_ann_df(rnaseq_norm_df)

# Clean tissue RNA-seq data 
rnaseq_norm_df = clean_tissue_rnaseq(rnaseq_norm_df)

# Covariates 
covars_df = read.delim2(covariates_path)

covars_df = clean_covars_df(covars_df)

# Reference (with extra info)
ref_df = read.csv(gencode_path)

ref_proteincoding_df = merge(ref_proteincoding_df, ref_df,
                             by = c('Transcript_name', 'Gene_name', 
                                    'Chromosome', 'Mean_position'),
                             all.x=TRUE, all.y=FALSE)


#--------------------------------------------------------------------#
# PIPELINE:                                                          #
# 1. Correct data for X covariates                                   #
# 2. For every chromosome:                                           #
#     a) Calculate correlations                                      #
#     b) Plot correlations                                           #
#     c) Calculate correlation (cor, distance)                       #
#--------------------------------------------------------------------#

## 2. For every chromosome:
results_df = data.frame(matrix(nrow=0, ncol=8))
names(results_df) = c('Chromosome', 'Covars', 'N_PC', 'N_PEERs', 
                      'Corrs_pearson', 'P_pearson', 
                      'Corrs_spearman', 'P_spearman')


chr_list = process_chr_txt(chr_txt)

for (chr_n in chr_list){
  
  print(str_interp("CHROMOSOME ${chr_n}:"))
  
  # a) Subset data for chromosome
  rnaseq_tissue_chr_df = subset_rnaseq_chr_n(ann_df, rnaseq_norm_df, chr_n)

  if (covars_yn == 'yes_covars'){ covariates = c('pcr', 'sex', 'platform')}
  if (covars_yn == 'no_covars'){ covariates = c()}
    
  covars_filename = str_interp("${covars_yn}_${n_pcs}pc_${n_peers}peers")
    
  # b) Residualize on covariates (covars, PCs, PEER factors)
  rnaseq_tissue_chr_resid_df = calculate_residuals(rnaseq_tissue_chr_df, covars_df,
                                                   covariates, n_pcs, n_peers)

  # write out residualized expresion
  write.csv(rnaseq_tissue_chr_resid_df, residualized_expression_outpath,
            row.names=FALSE)
                                                  
    
  # c) Calculate correlations and P-values
  list_cors = calculate_corrs_and_p(rnaseq_tissue_chr_resid_df, ann_df, chr_n)
    
  cor_r_df = list_cors[[1]]
  cor_p_df = list_cors[[2]]
    
  # d) Combine corrs and p values
  cor_df = transform_mat_to_df(cor_r_df, cor_p_df)
    
  # e) Plot correlogram 
  p_bonf = 0.05/nrow(cor_df)
  print(str_interp("There are ${nrow(cor_df)} independent tests."))
  print(str_interp("The p-value threshold I'm using to plot is ${p_bonf}."))
    
  # LARGE FILE, annotated with genes
  cor_r_ann_df = cor_r_df
  cor_p_ann_df = cor_p_df
    
  transcripts_df = data.frame(rownames(cor_r_df))
  names(transcripts_df) = 'transcript_id'
    
  transcripts_gene_df = merge(transcripts_df, 
                                ref_df %>% select(Transcript_name, Gene_name),
                                by.x='transcript_id', by.y='Transcript_name',
                                sort=FALSE)
    
  rownames(cor_r_ann_df) = transcripts_gene_df$Gene_name
  colnames(cor_r_ann_df) = transcripts_gene_df$Gene_name
    
  rownames(cor_p_ann_df) = transcripts_gene_df$Gene_name
  colnames(cor_p_ann_df) = transcripts_gene_df$Gene_name

    
  # f) Keep significant (P < P-bonferoni) & add distances between genes 
  cor_df = keep_sig_and_add_gene_distance(cor_df, ann_df, ref_proteincoding_df, p_bonf)
    
  print(str_interp("There are ${nrow(cor_df)} significant correlations."))
  
  # Save file with distances
  write.csv(cor_df, correlation_outpath,
            row.names=FALSE)
  
  print(str_interp("I saved the output at ${correlation_outpath}"))
  
}

endTime = Sys.time()

print(endTime - startTime)
