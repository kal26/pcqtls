#!/usr/bin/env Rscript

# from tami, I think this goes from expression to correlations?


#-----------------------------------------------------------------#
# EFFECT OF CORRECTING FOR N PEERs AND PCs ON CROSS-CORRELATION   #
#                                                                 #
# USER INPUTS:                                                    #
#   -tissue                  [character]                          #
#   -conditions_filename     [character]                          #
#                            (produced by make_conditions_df.R)   #
#                                                                 #
#  - chr_txt     [character]                                      #
#                (comma-separated,colon-separated, or int;        #
#                              e.g. 1:3,4,5 OR 1,2 OR 21)         #
#-----------------------------------------------------------------#

#-------------------------------------------------------------------#
# FOR TESTING:                                                      #
# tissue='Muscle_Skeletal'                                          #
# conditions_filename = 'conditions_both_covars_60peers_by3.csv'    #
# chr_txt = "22"                                                    #
#-------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------#
# Rscript effect_of_covars_on_coexpression.R Muscle_Skeletal conditions_both_covars_60peers_by3 22   #
#----------------------------------------------------------------------------------------------------#

rm(list=ls())

library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(corrplot)

startTime = Sys.time()

#-------------------------------#
# ARGUMENTS FROM COMMAND LINE   #
#-------------------------------#
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
  
} else if (length(args) > 1) {
  tissue = as.character(args[1])
  conditions_filename = as.character(args[2])
  chr_txt = as.character(args[3])
}


#-------------#
# CONSTANTS   #
#-------------#
dir_eqtl_data = '/oak/stanford/groups/smontgom/tami/eqtl_project/data/'
dir_eqtl_tmp = '/oak/stanford/groups/smontgom/tami/eqtl_project/tmp/'
dir_eqtl_output = '/oak/stanford/groups/smontgom/tami/eqtl_project/output/1_correlations/'


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
  
  cor_df[cor_df$P > p_thresh, 'Corr_coefficient'] = 0
  cor_df[cor_df$P > p_thresh, 'P'] = 1
  
  ref_df = ref_df %>% 
    select(Transcript_name, Start)
  
  cor_df1 = merge(cor_df, ref_df, by.x='Gene1', by.y='Transcript_name')
  cor_df2 = merge(cor_df, ref_df, by.x='Gene2', by.y='Transcript_name')
  
  cor_df = merge(cor_df1, cor_df2, by=c('Gene1', 'Gene2', 'Corr_coefficient', 'P'))
  cor_df$Distance = abs(cor_df$Start.x - cor_df$Start.y)

  cor_df = cor_df %>%
    select(Gene1, Gene2, Corr_coefficient, P, Distance)
  
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
ref_proteincoding_df = read.csv(
  str_interp("${dir_eqtl_data}/gene_annotations_protein_coding.csv"))

# Tissue RNA-seq data
rnaseq_norm_df = read.delim2(
  str_interp('${dir_eqtl_data}/${tissue}.v8.normalized_expression.txt'))

# Keep only protein-coding genes
rnaseq_norm_df = rnaseq_norm_df %>%
  filter(gene_id %in% ref_proteincoding_df$Transcript_name)


# Extract gene annotatons
ann_df = make_ann_df(rnaseq_norm_df)

# Clean tissue RNA-seq data 
rnaseq_norm_df = clean_tissue_rnaseq(rnaseq_norm_df)

# Covariates 
covars_df = read.delim2(
  str_interp('${dir_eqtl_data}/${tissue}.v8.covariates.txt'))

covars_df = clean_covars_df(covars_df)

# Reference file
ref_df = read.csv(
  str_interp("${dir_eqtl_data}/processed_gencode.v26.GRCh38.genes.gtf.csv"))

# Conditions with N_peers, N_pcs 
conditions_df = read.csv(str_interp("${dir_eqtl_tmp}/${conditions_filename}"))


#--------------------------------------------------------------------#
# PIPELINE:                                                          #
# 1. Correct data for X covariates                                   #
# 2. For every chromosome:                                           #
#     a) Calculate correlations                                      #
#     b) Plot correlations                                           #
#     c) Calculate correlation (cor, distance)                       #
#--------------------------------------------------------------------#

## For every chromosome:
results_df = data.frame(matrix(nrow=0, ncol=6))
names(results_df) = c('Chromosome', 'Covars', 'N_PC', 'N_PEERs', 
                      'Corrs_pearson', 'P_pearson')


# Set output dir and make tissue dir 
setwd(dir_eqtl_output)

if (!dir.exists(tissue)){
  dir.create(tissue)
}

setwd(tissue)

if (!dir.exists("effect_of_peers")){
  dir.create("effect_of_peers")
}

setwd("effect_of_peers")

chr_list = process_chr_txt(chr_txt)


for (chr_n in chr_list){
  
  print(str_interp("CHROMOSOME ${chr_n}:"))
  
  # a) Subset data for chromosome
  rnaseq_tissue_chr_df = subset_rnaseq_chr_n(ann_df, rnaseq_norm_df, chr_n)
  
  # Calculate stats for each condition 
  for (i in 1:nrow(conditions_df)){
    
    covars = conditions_df[i,'Covariates']
    if (covars == 'yes_covars'){ covariates = c('pcr', 'sex', 'platform')}
    if (covars == 'no_covars'){ covariates = c()}
    
    n_pcs = conditions_df[i,'N_PCs']
    n_peers = conditions_df[i, 'N_PEERs']
    
    print(str_interp("This is for ${n_peers} peers."))
    
    covars_filename = str_interp("${covars}_${n_pcs}pc_${n_peers}peers")
    
    # b) Residualize on covariates (covars, PCs, PEER factors)
    rnaseq_tissue_chr_resid_df = calculate_residuals(rnaseq_tissue_chr_df, covars_df,
                                                     covariates, n_pcs, n_peers)
    
    # c) Calculate correlations and P-values
    list_cors = calculate_corrs_and_p(rnaseq_tissue_chr_resid_df, ann_df, chr_n)
    
    cor_r_df = list_cors[[1]]
    cor_p_df = list_cors[[2]]
    
    # d) Combine corrs and p values
    cor_df = transform_mat_to_df(cor_r_df, cor_p_df)
    
    # e) Keep significant (P < P-bonferoni) & add distances between genes 
    p_bonf = 0.05/(nrow(cor_df))
    cor_df = keep_sig_and_add_gene_distance(cor_df, ann_df, ref_df, p_bonf)
    
    # Save file with distances
    write.csv(cor_df, str_interp('corr_df_chr${chr_n}_${covars_filename}_bonf_w_distances_proteincoding.csv'))
    
    # f) Plot 
    p_bonf = 0.05/nrow(cor_df)
    
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
    
    png(str_interp("corr_chr${chr_n}_${covars_filename}_LARGE_ANN_proteincoding.png"), width=1000, height=1000, type="cairo")
    corrplot(cor_r_ann_df, method="color", tl.col="black", tl.srt=45, tl.cex = 0.5,
             p.mat = cor_p_ann_df, sig.level = p_bonf, insig = "blank",
             type = 'upper', diag = FALSE, order = 'original',
             title = str_interp("${tissue} - Chromosome ${chr_n} - ${n_peers} PEERS - ${n_pcs} PCs"))
    dev.off()
    
  }
  
}

endTime = Sys.time()

print(endTime - startTime)
