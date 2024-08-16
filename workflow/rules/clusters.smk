
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
    params:
        gencode_path = gencode_path, 
        full_abc_path = full_abc_path, 
        abc_match_path = abc_match_path, 
        ctcf_match_path = ctcf_match_path, 
        ctcf_dir = ctcf_dir, 
        paralog_path = paralog_path,
        go_path = go_path, 
        cross_map_path = cross_map_path 
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
            --gencode_path {params.gencode_path} \
            --full_abc_path {params.full_abc_path} \
            --abc_match_path {params.abc_match_path} \
            --ctcf_match_path {params.ctcf_match_path} \
            --ctcf_dir {params.ctcf_dir} \
            --paralog_path {params.paralog_path} \
            --go_path {params.go_path} \
            --cross_map_path {params.cross_map_path} \
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
        exclude_cluster_genes = 1,
        distance_matched = 0,
        gencode_path = gencode_path, 
        full_abc_path = full_abc_path, 
        abc_match_path = abc_match_path, 
        ctcf_match_path = ctcf_match_path, 
        ctcf_dir = ctcf_dir, 
        paralog_path = paralog_path,
        go_path = go_path, 
        cross_map_path = cross_map_path 
    conda:
        'tensorqtl_r'
    output:
        annotated_nulls = clusters_dir + '{TISSUE}_null_{CLUSTER_SIZE}genes_annotated.csv'
    shell:"""
        python workflow/scripts/annotate_null_clusters.py \
            -t {wildcards.TISSUE} \
            -c {input.clusters} \
            -e {input.normalized_expression} \
            -co {input.full_covariates} \
            -o {output.annotated_nulls} \
            --gencode_path {params.gencode_path} \
            --full_abc_path {params.full_abc_path} \
            --abc_match_path {params.abc_match_path} \
            --ctcf_match_path {params.ctcf_match_path} \
            --ctcf_dir {params.ctcf_dir} \
            --paralog_path {params.paralog_path} \
            --go_path {params.go_path} \
            --cross_map_path {params.cross_map_path} \
            --cluster_size {wildcards.CLUSTER_SIZE} \
            --exclude_cluster_genes {params.exclude_cluster_genes} \
            --distance_matched {params.distance_matched} \
            --verbosity 1
        """