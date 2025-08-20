
# Output File Formats

## Gene Clusters (`clusters_dir/{TISSUE}.clusters.txt`)
Tab-separated file containing gene expression clusters:
- `chr`: Chromosome (e.g., 'chr1')
- `cluster_id`: Underscore-separated list of transcripts in each cluster
- `tissue_id`: Tissue identifier
- `num_genes`: Number of genes in the cluster
- `percent_correlated`: Percentage of correlations above threshold
- `mean_corr`: Mean correlation value
- `mean_pos_corr`: Mean positive correlation value
- `mean_neg_corr`: Mean negative correlation value

## Principal Components (`pc_output_dir/{TISSUE}.pcs.bed`)
BED format file containing principal components for each cluster:
- `chr`: Chromosome (e.g., 'chr1')
- `start`: First start position of any gene in the cluster (1-based)
- `end`: Last end position of any gene in the cluster (1-based)
- `gene_id`: Principal component ID (format: `{cluster_id}_pc{pc_number}`)
- Sample columns: One column per sample with PC values

## eQTL Summary Stats (`eqtl_output_dir/{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.parquet`)
Parquet file with eQTL associations:
- `phenotype_id`: Gene ID (format: `{cluster_id}_e_{gene_id}`)
- `variant_id`: Variant ID (format: `{chr}_{pos}_{ref}_{alt}`)
- `start_distance`: Distance from variant to phenotype start
- `end_distance`: Distance from variant to phenotype end
- `af`: In-sample ALT allele frequency of the variant
- `ma_samples`: Number of samples with at least one minor allele
- `ma_count`: Minor allele count
- `pval_nominal`: Nominal p-value of the association between the phenotype and variant
- `slope`: ERegression slope (beta)
- `slope_se`: Standard error of Regression slope

## pcQTL Summary Stats  (`pcqtl_output_dir/{TISSUE}/{TISSUE}.v8.pcs.cis_qtl_pairs.parquet`)
Parquet file with pcQTL associations:
- `phenotype_id`: Principal component ID (format: `{cluster_id}_pc{pc_number}`)
- `variant_id`: Variant ID (format: `{chr}_{pos}_{ref}_{alt}`)
- `start_distance`: Distance from variant to phenotype start
- `end_distance`: Distance from variant to phenotype end
- `af`: In-sample ALT allele frequency of the variant
- `ma_samples`: Number of samples with at least one minor allele
- `ma_count`: Minor allele count
- `pval_nominal`: Nominal p-value of the association between the phenotype and variant
- `slope`: ERegression slope (beta)
- `slope_se`: Standard error of Regression slope

## SuSiE Fine-mapping Results (`{eqtl/pcqtl}_output_dir/{TISSUE}/{TISSUE}.v8.{cluster_genes/pcs}.susie_R.txt`)
Tab-separated file with SuSiE fine-mapping results:
- `phenotype_id`: Gene ID or PC ID
- `variant_id`: Variant ID (format: `{chr}_{pos}_{ref}_{alt}`)
- `pip`: Posterior inclusion probability
- `cs_id`: Credible set ID 

## Annotated Clusters (`annotations_output_dir/{TISSUE}/{TISSUE}.clusters.annotated.txt`)
Tab-separated file with functional annotations for clusters:
- `cluster_id`: Cluster identifier
- `gene_id`: Ensembl gene ID
- `gene_name`: Gene symbol
- `chromosome`: Chromosome
- `start`: Gene start position
- `end`: Gene end position
- `strand`: Gene strand
- `cluster_size`: Number of genes in cluster
- `cluster_center_chr`: Chromosome of cluster center
- `cluster_center_pos`: Position of cluster center
- `abc_enhancer_count`: Number of ABC enhancers connected to gene
- `ctcf_binding_sites`: Number of CTCF binding sites in cluster
- `tad_boundary_distance`: Distance to nearest TAD boundary
- `paralog_count`: Number of paralogs in cluster
- `go_terms`: Gene ontology terms (comma-separated)
- `expression_variance`: Variance in expression across samples


## QTL Colocalization (`coloc_output_dir/pairs/{TISSUE}.v8.pairs_coloc.{CHROM}.txt`):
Tab-separated file with colocalizations between QTLs:
- `qtl1_id`: Gene ID or PC ID for first phenotype
- `qtl2_id`: Gene ID or PC ID for second phenotype
- `nsnps`: Number of SNPs in the region
- `hit1`: Lead variant ID of first phenotype (format: `{chr}_{pos}_{ref}_{alt}`)
- `hit2`: Lead variant ID of second phenotype (format: `{chr}_{pos}_{ref}_{alt}`)
- `PP.H0.abf`: Posterior probability of no association
- `PP.H1.abf`: Posterior probability of association with first phenotype only
- `PP.H2.abf`: Posterior probability of association with second phenotype only
- `PP.H3.abf`: Posterior probability of association with both, different causal variants
- `PP.H4.abf`: Posterior probability of association with both, same causal variant (colocalization)
- `idx1`: Credible set ID of first phenotype
- `idx2`: Credible set ID of second phenotype


## GWAS Colocalization (`coloc_output_dir/gwas/{TISSUE}/{TISSUE}.v8.{GWAS}.susie_True.gwas_coloc.txt`):
Tab-separated file with colocalizations between GWAS and QTLs:
- `gwas_id`: GWAS identifier
- `qtl_id`: Gene ID or PC ID
- `nsnps`: Number of SNPs in the region
- `hit1`: Lead variant ID of GWAS (format: `{chr}_{pos}_{ref}_{alt}`)
- `hit2`: Lead variant ID of QTL (format: `{chr}_{pos}_{ref}_{alt}`)
- `PP.H0.abf`: Posterior probability of no association
- `PP.H1.abf`: Posterior probability of association with QTL only
- `PP.H2.abf`: Posterior probability of association with GWAS only
- `PP.H3.abf`: Posterior probability of association with both, different causal variants
- `PP.H4.abf`: Posterior probability of association with both, same causal variant (colocalization)
- `idx1`: Credible set ID of first phenotype
- `idx2`: Credible set ID of second phenotype


## QTL Signal Groups (`coloc_output_dir/qtl_signal_groups/{TISSUE}.qtl_signal_groups.txt`):
Tab-separated file with credible set groups for QTLs:
- `signal_id`: Unique identifier for signal group (dash-separated list of credible set IDs)
- `cluster_id`: Cluster identifier
- `num_e_coloc`: Number of eQTL signals in the group
- `num_pc_coloc`: Number of pcQTL signals in the group
- `cluster_id`: Cluster identifier
- `lead_var_set`: List of lead variants in the signal group

## GWAS Signal Groups (`coloc_output_dir/gwas_signal_groups/{TISSUE}.{GWAS}.gwas_signal_groups.txt`):
Tab-separated file with credible set groups for GWAS and QTLs:
- `signal_id`: Unique identifier for signal group (dash-separated list of credible set IDs)
- `cluster_id`: Cluster identifier
- `num_qtl_coloc`: Number of QTL signals in the group
- `num_gwas_coloc`: Number of GWAS signals in the group
- `num_e_coloc`: Number of eQTL signals in the group
- `num_pc_coloc`: Number of pcQTL signals in the group
