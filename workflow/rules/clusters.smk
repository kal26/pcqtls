
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
        max_cluster_size = max_cluster_size,
        min_cluster_size = min_cluster_size,
        min_corr_cutoff = min_corr_cutoff,
        percent_corr_cutoff = percent_corr_cutoff,
        cutoff_type = 'pvalue', # or can be 'value'
    resources:
        mem = "80G", 
        time = "2:00:00",
    conda:
        'tensorqtl_r'
    output:
        clusters = clusters_dir + '{TISSUE}_clusters_all_chr.csv'
    shell:"""
        python workflow/scripts/call_clusters.py \
            -e {input.normalized_expression} \
            -co {input.full_covariates} \
            -o {output.clusters} \
            --max_cluster_size {params.max_cluster_size} \
            --min_cluster_size {params.min_cluster_size} \
            --min_corr_cutoff {params.min_corr_cutoff} \
            --percent_corr_cutoff {params.percent_corr_cutoff} \
            --cutoff_type {params.cutoff_type} \
            --tissue_id {wildcards.TISSUE}
        """
