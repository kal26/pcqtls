## following tutorial at https://hanruizhang.github.io/GWAS-eQTL-Colocalization/
## in r/4.2.2

setwd('/home/klawren/oak/pcqtls/')
library(coloc)
library(arrow)
library(tidyverse)

# test data
data(coloc_test_data)
coloc_test_data$D1

# paths to data
eqtl_path <- 'output/proteincoding_main/control_eqtl/Adipose_Subcutaneous/Adipose_Subcutaneous.v8.cluster_genes.cis_qtl_pairs.chr22.parquet'
gwas_path <- '/oak/stanford/groups/smontgom/shared/gwas_summary_stats/barbeira_gtex_imputed/imputed_gwas_hg38_1.1/imputed_UKB_50_Standing_height.txt.gz'



eqtl_path <- 'output/temp/test_eqtl.csv'
gwas_path <- 'output/temp/test_gwas.csv'

#load in eqtl data
#eqtl <- read_parquet(eqtl_path)
eqtl <- read.table(eqtl_path, header=T, as.is=T)
# load in gwas data
gwas <- read.table(gwas_path, header=T, as.is=T)

eqtl[1,'phenotype_id']
cluster_id = strsplit(as.character(eqtl[1,'phenotype_id']), '_e_')[[1]][1]
cluster_id

eqtl <- eqtl %>%
  mutate(cluster_id = sapply(strsplit(as.character(phenotype_id), '_e_'), function(x) x[[1]]))



# single variant version





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





# susie version



# get LD matrix
ld <- read.table('output/temp/test.ld')
snp_list <- eqtl$variant_id
print("total snps")
print(length(snp_list))
rownames(ld) <- snp_list
colnames(ld) <- snp_list

# drop snps with missing values from LD for eqtl
rows_with_missing <- rownames(ld)[which(rowSums(is.na(ld)) > 0)]
cols_with_missing <- colnames(ld)[which(colSums(is.na(ld)) > 0)]

ld_missing_snps = union(rows_with_missing, cols_with_missing)
print("Number snps with ld missing")
print(length(ld_missing_snps))

cleaned_ld_matrix <- ld[!rownames(ld) %in% ld_missing_snps, !colnames(ld) %in% ld_missing_snps]
cleaned_eqtl <- eqtl_filtered[!eqtl_filtered$variant_id %in% ld_missing_snps, ]

cleaned_eqtl_list <- list(beta = cleaned_eqtl$slope, 
                          varbeta = cleaned_eqtl$slope_se**2, 
                          snp = cleaned_eqtl$variant_id, 
                          position = cleaned_eqtl$position, 
                          type = "quant", 
                          N = 581, 
                          MAF = cleaned_eqtl$af,
                          sdY = 1,
                          LD = as.matrix(cleaned_ld_matrix))


check_dataset(cleaned_eqtl_list,req="LD")


# susie needs effect sizes, so we must also drop the snps with na for gwas effect
gwas_missing_snps <- gwas[is.na(gwas$effect_size), 'panel_variant_id']
missing_snps = union(ld_missing_snps, gwas_missing_snps)
print("Number snps with ld or gwas missing")
print(length(missing_snps))
cleaned_gwas <- gwas[!gwas$panel_variant_id %in% missing_snps, ]
gwas_cleaned_ld_matrix <- ld[!rownames(ld) %in% missing_snps, !colnames(ld) %in% missing_snps]

cleaned_gwas_list <- list(MAF = cleaned_gwas$frequency, 
                  snp = cleaned_gwas$panel_variant_id, 
                  position = cleaned_gwas$position, 
                  type = 'quant', 
                  N = 337119, 
                  pvalues=cleaned_gwas$pvalue,
                  LD = as.matrix(gwas_cleaned_ld_matrix), 
                  beta = cleaned_gwas$effect_size,
                  varbeta = cleaned_gwas$standard_error**2)

check_dataset(cleaned_gwas_list,req="LD")


# run susie finemapping on both datasets
eqlt_susie = runsusie(cleaned_eqtl_list)
summary(eqlt_susie)
plot_dataset(cleaned_eqtl_list, susie_obj = eqtl_susie)


gwas_susie = runsusie(cleaned_gwas_list)
summary(gwas_susie)


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

cleaned_gwas[,c("rsid","pval")]
cleaned_eqtl[,c("rsid","pval")]
