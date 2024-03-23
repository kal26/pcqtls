# from dan

# "Here is the code snippet for gene expression normalization and PCA"

library(corral)
library(scran)
library(BiocSingular)
library(PCAtools)
library(RNOmni)

library(vroom)
library(tidyverse)

# Read in counts with vroom.
counts_raw <- vroom("PATH_TO_COUNTS")

# Use two methods from scran to identify sufficiently variable genes for downstream analysis and subset the data.
set.seed(12345L)
counts_gene_var <- scran::modelGeneVarByPoisson(counts_raw)
counts_gene_var_top <- scran::getTopHVGs(counts_gene_var)
counts_gene_cv2 <- scran::modelGeneCV2(counts_filtered)
counts_gene_cv2_top <- scran::getTopHVGs(counts_gene_cv2, var.field = "ratio", var.threshold = 1L)
counts_gene_shared <- intersect(counts_gene_var_top, counts_gene_cv2_top)

counts_variable <- counts_raw[counts_gene_shared, ]

# Use corral to normalize counts for library depth and gene length, and use RNOmni for the inverse normal transformation.
# Use the standardized_data for QTL mapping directly.
standardized_data <- corral::corral_preproc(counts_variable) %>% apply(1, RNOmni::RankNorm) %>% t()

# Run PCA for hidden factor estimation.
# Use the Gavish-Donoho method to choose the number of PCs.
# Alternatively, can use chooseMarcenkoPastur() for the less conservative Marcenko-Pastur method.
pca_standardized <- PCAtools::pca(standardized_data, BSPARAM = BiocSingular::ExactParam()(), rank = init_npcs)
n_pcs <- PCAtools::chooseGavishDonoho(standardized_data, var.explained = pca_standardized$sdev^2, noise = 1)

# You will probably use this object for your hidden factors.
pca_eigenvectors <- pca_standardized$rotated[,1:n_pcs]

# Optionally, run rank normal transformation on PCs - this is not normally done but is provided here as an example.
pca_eigenvectors_standardized <- apply(pca_standardized$rotated[,1:n_pcs], 2, RNOmni::RankNorm)