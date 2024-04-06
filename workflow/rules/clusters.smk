rule call_clusters: 
    input:
        normalized_expression
        full_covariates
    output:
        clusters

rule get_covariates:
    input: 
        gtex_provided_covariates
        normalized_expression
    output:
        full_covariates



rule get_covariate_conditions:
## if you want to compare different number of pcs/peers
    input:
    params:
        covars_yn = covars_yn,
        pc_or_peer = pc_or_peer,
        n_factors_max = n_pcs,
        factors_break = factors_break
    output:
        covariate_ conditions = corr_dir + 'conditions_covars_peers_60_by_15'
    shell:"""
        module load r/4.0
        Rscript get_covariate_conditions.R {params.covars_yn} \
        {params.pc_or_peer} \
        {params.n_factors_max} \
        {params.factors_break} \
        {output.conditions}
        """

rule get_correlations:
    input: 
        gene_annotations = 'data/references/gene_annotations_protein_coding.csv',
        normalized_expression = expression_dir + '{TISSUE}.v8.normalized_expression.bed',
        covariates = covariates_dir + '{TISSUE}.v8.covariates.txt'
        gencode = 'data/references/processed_gencode.v26.GRCh38.genes.gtf.csv'
    params:
        covars_yn = covars_yn,
        n_peers = n_peers,
        n_pcs = n_pcs
    output:
        correlations = corr_dir + '{TISSUE}/sig_corr_df_chr{CHROM}_{covars_yn}_{n_pcs}pcs_{n_peers}peers_w_distances.csv'
        #?? full_covariates = covariates_dir + '{TISSUE}.v8.full_covariates.txt',
        residualized_expression = expression_dir + '{TISSUE}.v8.residualized_expression.bed'

    shell:"""
        module load r/4.0
        Rscript ../scripts/calculate_expression_correlation.R {wildcars.TISSUE} \
        {params.covars_yn} {params.npeers} \
        {params.n_pcs} {wildcards.CHROM} \
        {input.gene_annotations} {input.normalized_expression} \
        {input.covariates} {input.gencode} \
        {output.correlations} {output.full_covariates} \
        {output.residualized_expression}
        """


rule call_clusters:
    input:
        correlations = corr_dir + '{TISSUE}/sig_corr_df_chr{CHROM}_{covars_yn}_{n_pcs}pcs_{n_peers}peers_w_distances.csv'
        gene_annotations = 'data/references/gene_annotations_protein_coding.csv'
        normalized_expression = expression_dir + '{TISSUE}.v8.normalized_expression.bed'
    params:
        perc_threshold = perc_threshold,
        c_size_max = c_size_max,
        c_size_min = c_size_min 
    output:
        clusters = clusters_dir + '{TISSUE}/clusters_chr{CHROM}.csv',
        clusters_metadata =  clusters_dir + '{TISSUE}/clusters_chr{CHROM}_metadata.csv'

    shell:"""
        module load r/4.0
        Rscript ../scripts/call_clusters.R {wildcards.TISSUE} \
        {params.perc_threshold} {params.c_size_max} {params.c_size_min} {wildcards.CHROM} \
        dir_eqtl_output \
        dir_eqtl_data}/gene_annotations_protein_coding.csv \
        dir_eqtl_data}/${tissue}.v8.normalized_expression.txt \ $### add args in right order
        """
   
