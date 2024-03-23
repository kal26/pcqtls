#!/usr/bin/env Rscript

# from tami, I think this goes from correlations to clusters

#-----------------------------------------------------------#
# CALCULATE CROSS-CORRELATIONS FOR USER-SPECIFIED:          #
#  - tissue           [character]                           #
#  - perc_threshold   [integer]     (0:100)                 #
#  - c_size_max      [integer]     (0:60)                   #
#  - c_size_min       [integer]     (0:5)                   #
#                                                           #
# FOR TESTING:                                              #
# tissue='Muscle_Skeletal'                                  #
# perc_threshold = 70                                       #
# c_size_max = 50                                           #
# c_size_min = 2                                            #
# chr_txt = 22                                              #
#-----------------------------------------------------------#

#-------------------------------------------------------#
# Rscript call_clusters.R Muscle_Skeletal 70 50 2 22    #
#-------------------------------------------------------#

rm(list=ls())

library(dplyr)
library(stringr)

startTime = Sys.time()

#-------------------------------#
# ARGUMENTS FROM COMMAND LINE   #
#-------------------------------#
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
  
} else if (length(args) > 1) {
  
  tissue = args[1]
  perc_threshold = as.numeric(args[2])
  c_size_max = as.numeric(args[3])
  c_size_min = as.numeric(args[4])
  chr_txt = args[5]
  corr_path = args[6]
  clusters_outpath = args[7]
  clusters_meta_outpath = args[8]
  gene_annotations_path = args[9]
  normalized_expression_path = args[10]
}


#------------#
# FUNCTIONS  #
#------------#
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

make_ann_df = function(rnaseq_norm_df){
  
  ann_df = rnaseq_norm_df %>% select(X.chr:gene_id)
  
  names(ann_df) = c('Chromosome', 'Start', 'End', 'Gene_id')
  ann_df$Mean_position = (ann_df$Start + ann_df$End)/2
  
  ann_df$Distance_rank = seq(1, nrow(ann_df))
  
  return(ann_df)
}

prepare_transcript_names = function(rnaseq_norm_df, chr_n){
  
  # Extract gene annotatons
  ann_df = make_ann_df(rnaseq_norm_df)
  
  # Clean rnaseq data
  rnaseq_norm_df = clean_tissue_rnaseq(rnaseq_norm_df)
  
  # Extract chromosome data 
  rnaseq_tissue_chr_df = subset_rnaseq_chr_n(ann_df, rnaseq_norm_df, chr_n)
  
  # Extract genes and order them by rank of expression
  ann_df = ann_df %>%
    filter(Gene_id %in% colnames(rnaseq_tissue_chr_df))
  
  genes = ann_df[order(ann_df$Distance_rank), 'Gene_id']
  
  return(genes)
  
}

turn_corr_df_to_matrix = function(corr_df, transcript_names){
  
  df = corr_df[,1:3]
  
  mat = data.frame(matrix(nrow=length(transcript_names), 
                          ncol=length(transcript_names)))
  rownames(mat) = transcript_names
  colnames(mat) = transcript_names
  
  for (i in 1:nrow(mat)){
    for (j in 1:ncol(mat)){
      
      cor1 = df[df$Gene1 == rownames(mat)[i] &
                  df$Gene2 == colnames(mat)[j], 'Corr_coefficient']
      
      cor2 = df[df$Gene2 == rownames(mat)[i] &
                  df$Gene1 == colnames(mat)[j], 'Corr_coefficient']
      
      cor = unique(c(cor1, cor2))
      
      if (length(cor) == 0){
        mat[i,j] = NA
      } else {
        mat[i,j] = cor
      }
      
    }
  }
  
  mat[upper.tri(mat)] = NA
  
  return(mat)
}

convert_transcripts_to_genes = function(tmp, ref_df){
  
  transcripts_df = data.frame(names(tmp))
  names(transcripts_df) = 'Transcript_name'
  
  transcripts_genes_df = merge(transcripts_df, 
                               ref_df %>% select(Transcript_name, Gene_name),
                               all.x = TRUE, all.y = FALSE)
  
  genes_txt = paste(transcripts_genes_df$Gene_name, collapse=', ')
  
  return(genes_txt)
}

call_clusters = function(corr_mat, perc_threshold, c_sizes){
  
  clusters_df = data.frame(matrix(nrow=0, ncol=7))
  names(clusters_df) = c('N_genes', 'Transcripts', 'Genes', 'Perc_cor',
                         'Mean_cor', 'Mean_pos_cor', 'Mean_neg_cor')
  
  i=1
 
  for (c_size in c_sizes){
    
    while (i < nrow(corr_mat)){
      
      j = i+c_size-1
      
      # Subset tmp matrix
      if (j <= ncol(corr_mat)){
        tmp = corr_mat[i:j, i:j]
      } else {
        tmp = corr_mat[i:(nrow(corr_mat)), i:(ncol(corr_mat))]
      }
      
      cluster_size = nrow(tmp)
      
      # Count number of correlations 
      n_max = (cluster_size*(cluster_size-1))/2
      n_cors = (cluster_size)^2 - sum(is.na(tmp))
      
      # Percent correlated 
      perc_cor = round((n_cors/n_max)*100, 3)
      
      # If % correlated is above threshold 
      if (perc_cor >= perc_threshold){
        
        cors = tmp[is.na(tmp)==FALSE]
        
        cor = round(mean(cors), 3)
        cor_pos = round(mean(cors[cors > 0]), 3)
        cor_neg = round(mean(cors[cors < 0]), 3)
        
        # Has the cluster been recorded before as a part of a larger one?
        transcripts_txt = paste(names(tmp), collapse=',')
        genes_txt = convert_transcripts_to_genes(tmp, ref_df)
        
        # Only record cluster if not contained within a larger one 
        if (nrow(clusters_df) != 0){
          if (!any(grepl(transcripts_txt, clusters_df$Transcripts))){
            
            # Record cluster
            clusters_df[nrow(clusters_df)+1, ] =
              c(cluster_size, transcripts_txt, genes_txt,
                perc_cor, cor, cor_pos, cor_neg)
          }
        } else {

          # Record cluster
          clusters_df[nrow(clusters_df)+1, ] =
            c(cluster_size, transcripts_txt, genes_txt,
              perc_cor, cor, cor_pos, cor_neg)
        }

        i = i+cluster_size
        
      } else {
        i = i+1 
      }
      
    }
    
    i=1
    
  }
  
  return(clusters_df)
  
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


#------------#
# VARIABLES  #
#------------#
#dir_eqtl_output = '/oak/stanford/groups/smontgom/tami/eqtl_project/output/1_correlations/'
#dir_eqtl_data = '/oak/stanford/groups/smontgom/tami/eqtl_project/data/'
#dir_gtex_refs = '/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/references'

c_sizes = c(c_size_max:c_size_min)
chr_list = process_chr_txt(chr_txt)


#------------#
# LOAD DATA  #
#------------#

# Reference 
ref_df = read.csv(gene_annotations_path)

# Expression data for tissue 
rnaseq_norm_df = read.delim2(normalized_expression_path)

# Keep only protein-coding genes
rnaseq_norm_df = rnaseq_norm_df %>%
  filter(gene_id %in% ref_df$Transcript_name)


#------------#
# ANALYSIS   #
#------------#

# Clusters data 
main_clusters_df = data.frame(matrix(nrow=0, ncol=8))
names(main_clusters_df) = c('N_genes', 'Transcripts', 'Genes', 'Perc_cor',
                            'Mean_cor', 'Mean_pos_cor', 'Mean_neg_cor',
                            'Chromosome')

# Meta data
meta_df = data.frame(matrix(nrow=0, ncol=3))
names(meta_df) = c('Tissue', 'Chromosome', 'N_Genes')

for (chr_n in chr_list){
    
  # Correlation results 
  corr_df = read.csv(corr_path)
    
  # Get all transcript and gene names for chromosome
  transcript_names = prepare_transcript_names(rnaseq_norm_df, chr_n)
  
  # Record meta data 
  n_genes = length(transcript_names)
  meta_df[nrow(meta_df)+1, ] = c(tissue, chr_n, n_genes)
    
  # Turn correlation df to matrix
  corr_mat = turn_corr_df_to_matrix(corr_df, transcript_names)
    
  # Identify clusters
  clusters_df = call_clusters(corr_mat, perc_threshold, c_sizes)
  clusters_df$Chromosome = chr_n
    
  # Add to main clusters df 
  main_clusters_df = rbind(main_clusters_df, clusters_df)
    
}

main_clusters_df$Tissue = tissue

# Save clusters and meta data 
write.csv(main_clusters_df,
          clusters_outpath,
          row.names=FALSE)

write.csv(meta_df,
          clusters_meta_outpath, 
          row.names=FALSE)

# Print execution time
endTime = Sys.time()
print(endTime - startTime)

