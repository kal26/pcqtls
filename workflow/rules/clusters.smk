
rule get_covariates:
    input: 
        normalized_expression =  expression_dir + '{TISSUE}.v8.normalized_expression.bed',
        gtex_provided_covariates =  covariates_dir + '{TISSUE}.v8.gtex_covariates.bed'
    output:
        full_covariates =  covariates_dir + '{TISSUE}.v8.covariates.txt'
    script:
        '../scripts/get_covariates.py'


rule call_clusters: 
    input:
        normalized_expression = expression_dir + '{TISSUE}.v8.normalized_expression.bed',
        full_covariates = covariates_dir + '{TISSUE}.v8.covariates.txt'
    params:
        max_cluster_size = 50,
        min_cluster_size = 2,
        min_corr_cutoff = 0.1,
        percent_corr_cutoff = .7,
        cutoff_type = 'pvalue' # or can be 'value'
    conda:
        'tensorqtl_r'
    output:
        clusters = clusters_dir + '{TISSUE}_clusters_all_chr.csv'
    script:
        '../scripts/call_clusters.py'
