## following tutorial at https://hanruizhang.github.io/GWAS-eQTL-Colocalization/
## in r/4.2.2

setwd('/home/klawren/oak/pcqtls/')
library(coloc)
library(arrow)
library(tidyverse)
library(data.table)
library(Rfast)

####### funtions #########

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



###### main #######

# paths to data
ld_path <- 'output/proteincoding_main/gwas_coloc/Lung/temp/ENSG00000275464.4_ENSG00000280071.3_ENSG00000280433.1.ld'
snp_list_path <- 'output/proteincoding_main/gwas_coloc/Lung/temp/ENSG00000275464.4_ENSG00000280071.3_ENSG00000280433.1.snp_list.txt'

cluster_id <- 'ENSG00000275464.4_ENSG00000280071.3_ENSG00000280433.1'

eqtl_path <- 'output/proteincoding_main/control_eqtl/Lung/Lung.v8.cluster_genes.cis_qtl_pairs.chr21.parquet'
pcqtl_path <- '/home/klawren/oak/pcqtls/output/proteincoding_main/pcqtl/Lung/Lung.v8.pcs.cis_qtl_pairs.chr21.parquet'



gwas_meta_path <- '/home/klawren/oak/pcqtls/data/references/gwas_metadata.txt'
gtex_meta_path <- '/home/klawren/oak/pcqtls/data/references/gtex_sample_sizes.csv'
tissue_id <- 'Lung'
gwas_folder <- '/oak/stanford/groups/smontgom/shared/gwas_summary_stats/barbeira_gtex_imputed/imputed_gwas_hg38_1.1'
snp_path_head <- 'output/temp/'


# get number gtex samples form gtex_meta and tisuse
gtex_meta <- read.table(gtex_meta_path, sep='\t', header = T)
num_gtex_samples <- gtex_meta[gtex_meta$tissue_id == tissue_id, 'sample_size']

# gwas metadata
gwas_meta <- fread(gwas_meta_path)
gwas_meta <- gwas_meta[37:38]


#load in eqtl data
eqtl <- read_parquet(eqtl_path)
pcqtl <- read_parquet(pcqtl_path)


gwas_coloc_results <- '/home/klawren/oak/pcqtls/output/proteincoding_main/gwas_coloc/Lung/temp/ENSG00000275464.4_ENSG00000280071.3_ENSG00000280433.1.coloc_qtl.txt'
this_cluster_coloc <- read.table(gwas_coloc_results, sep='\t', header=T)



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


###new

get_qtl_pairwise_coloc <- function(qtls_for_coloc, qtl_susies){
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
  
  qtl_pair_idxs <- combn(seq_along(qtls_for_coloc), 2, simplty=TRUE)

  for (i in 1:ncol(qtl_pair_idxs)) {
    qtl1 <- qtls_for_coloc[[qtl_pair_idxs[1, i]]]
    qtl2 <- qtls_for_coloc[[qtl_pair_idxs[2, i]]]
    susie1 <- qtl_susies[[qtl_pair_idxs[1, i]]]
    susie2 <- qtl_susies[[qtl_pair_idxs[2, i]]]
    cat(paste("coloc for", qtl1$phenotype_id , "-", qtl2$phenotype_id), "\n")
    this_coloc <- run_coloc(susie1,susie2)
    this_coloc$qtl1_id <- qtl1$phenotype_id
    this_coloc$qtl2_id <- qtl2$phenotype_id
    qtl_coloc_results <- rbind(qtl_coloc_results, this_coloc)
  }
  return(qtl_coloc_results)
}


qtl_susies <- lapply(qtls_for_coloc, runsusie)

qtl_coloc_results <- get_qtl_pairwise_coloc(qtls_for_coloc, qtl_susies)

write.table(qtl_coloc_results, file='output/temp/coloc_qtl.txt', quote=FALSE, row.names=FALSE, sep='\t')





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

qtls_for_coloc <- eqtls_for_coloc[1:4][!sapply(eqtls_for_coloc[1:4], is.null)]

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
write.table(gwas_coloc_results, file='output/temp/coloc.txt', row.names=FALSE, sep='\t')




dim(gwas_meta)[[0]]
for (i in seq_along(qtl_susies)){
  print(i)
}












susie_gwas = runsusie(gwas_for_coloc)
summary(susie_gwas)
plot_dataset(gwas_for_coloc, susie_obj = susie_gwas)


susie_eqtl = runsusie(eqtl_for_coloc)
summary(susie_eqtl)
plot_dataset(eqtl_for_coloc, susie_obj = susie_eqtl)


coloc=coloc.susie(susie_gwas,susie_eqtl)
print("test")
print(coloc$summary)


merged <- merge(test_eqtl, gwas_filtered, by="rsid", all=FALSE, suffixes=c("1","2"))
plot()


eqtl_for_coloc

#### debug ####


##### ploting ###

plot_merged <- function(gwas_filtered, test_eqtl){
  
  gwas_filtered$rsid <- gwas_filtered$panel_variant_id
  gwas_filtered$pval <- gwas_filtered$pvalue
  gwas_filtered$logp <- -log(gwas_filtered$pval)
  
  test_eqtl$rsid <- test_eqtl$variant_id
  test_eqtl$pval <- test_eqtl$pval_nominal
  test_eqtl$logp <- -log(test_eqtl$pval)
  
  merged <- merge(test_eqtl, gwas_filtered, by="rsid", all=FALSE, suffixes=c("1","2"))
  plot(x=merged$logp1, y=merged$logp2)
  
}



# susie version


print(summary(gwas_susie))

susie1 = runsusie(clean_eqtl_list)
print(summary(susie1))
plot_dataset(cleaned_eqtl_list, susie_obj = eqtl_susie)





# colocalize the susie resutls
# note that snps do not need to perfectly match for this
susie.res=coloc.susie(gwas_susie,eqlt_susie)
print("test")
print(susie.res$summary)
plot_dataset(cleaned_eqtl_list, susie_obj = eqtl_susie)
plot_dataset(susie.res)

library(locuscomparer)
cleaned_gwas$pval <- cleaned_gwas$pvalue
cleaned_gwas$rsid <- cleaned_gwas$panel_variant_id
cleaned_eqtl$pval <- cleaned_eqtl$pval_nominal
cleaned_eqtl$rsid <- cleaned_eqtl$variant_id

locuscompare(in_fn1 = c(cleaned_eqtl$rsid, cleaned_eqtl$pvall) , in_fn2 =c(cleaned_gwas$rsid, cleaned_gwas$pvall))






# filter to variants in both 
gwas_filtered <- gwas[gwas$panel_variant_id %in% eqtl$variant_id, ]
eqtl_filtered <- eqtl[eqtl$variant_id %in% gwas_filtered$panel_variant_id, ]

# make some columns to play well with coloc
eqtl_filtered <- eqtl_filtered %>%
  mutate(position = as.integer(str_split(variant_id, "_") %>% sapply(pluck, 2)))

# sdY is 1 as expression is normalized
# TODO file for eqtl sample size
eqtl_list <- list(beta = eqtl_filtered$slope, 
                  varbeta = eqtl_filtered$slope_se**2, 
                  snp = eqtl_filtered$variant_id, 
                  position = eqtl_filtered$position, 
                  type = "quant", 
                  N = 581, 
                  MAF = eqtl_filtered$af,
                  sdY = 1)
check_dataset(eqtl_list)

# plot the eqtl
plot_dataset(eqtl_list)

# check is case-control or quantitative for gwas 
# TODO make a file for this
# is sdY = 1?????????
gwas_list <- list(MAF = gwas_filtered$frequency, 
                  snp = gwas_filtered$panel_variant_id, 
                  position = gwas_filtered$position, 
                  type = 'quant', 
                  N = 337119, 
                  #sdY = 1, 
                  pvalues=gwas_filtered$pvalue)
check_dataset(gwas_list)


# plot the gwas
gwas_filtered$logpval <- -log10(gwas_filtered$pvalue)
plot_dataset(gwas_filtered, alty=gwas_filtered$logpval)

result <- coloc.abf(dataset1=gwas_list, dataset2=eqtl_list)
result

sensitivity(result,rule="H4 > 0.5") 



