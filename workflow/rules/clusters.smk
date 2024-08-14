
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


rule annotate_clusters:
    input:
        normalized_expression = expression_dir + '{TISSUE}.v8.normalized_expression.bed',
        full_covariates = covariates_dir + '{TISSUE}.v8.covariates.txt',
        clusters = clusters_dir + '{TISSUE}_clusters_all_chr.csv'
    resources:
        mem = "80G", 
        time = "2:00:00"
    conda:
        'tensorqtl_r'
    output:
        annotated_clusters = clusters_dir + '{TISSUE}_clusters_annotated.csv'
    shell:"""
        python workflow/scripts/annotate_clusters.py \
            -t {wildcards.TISSUE} \
            -c {input.clusters} \
            -e {input.normalized_expression} \
            -co {input.full_covariates} \
            -o {output.annotated_clusters} \
            --verbosity 1
        """

rule annotate_nulls:
    input:
        normalized_expression = expression_dir + '{TISSUE}.v8.normalized_expression.bed',
        full_covariates = covariates_dir + '{TISSUE}.v8.covariates.txt',
        clusters = clusters_dir + '{TISSUE}_clusters_all_chr.csv'
    resources:
        mem = "80G", 
        time = "2:00:00"
    params:
        gencode = 'data/references/processed_gencode.v26.GRCh38.genes.csv',
        cluster_size = 2,
        exclude_cluster_genes = 1,
        distance_matched = 0
    conda:
        'tensorqtl_r'
    output:
        annotated_nulls = clusters_dir + '{TISSUE}_null_annotated.csv'
    shell:"""
        python workflow/scripts/annotate_null_clusters.py \
            -t {wildcards.TISSUE} \
            -c {input.clusters} \
            -e {input.normalized_expression} \
            -co {input.full_covariates} \
            -o {output.annotated_nulls} \
            -g {params.gencode} \
            --cluster_size {params.cluster_size} \
            --exclude_cluster_genes {params.exclude_cluster_genes} \
            --distance_matched {params.distance_matched} \
            --verbosity 1
        """