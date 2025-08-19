"""
Gene Expression Clustering Rules

This module contains rules for identifying co-expressed gene clusters
from normalized expression data.
"""

rule call_clusters:
    """
    Identify co-expressed gene clusters from normalized expression data.
    
    This rule uses correlation-based clustering to group genes with similar
    expression patterns across samples.
    """
    input:
        expression = config['expression_dir'] + '{TISSUE}.v8.normalized_expression.bed',
        covariates = config['covariates_dir'] + '{TISSUE}.v8.covariates.txt'
    
    output:
        clusters = config['clusters_dir'] + '{TISSUE}.clusters.txt'
    
    params:
        max_cluster_size = 50,
        min_cluster_size = 2,
        min_corr_cutoff = config['min_corr_cutoff'],
        percent_corr_cutoff = float(config['percent_corr_cutoff']),
        cutoff_type = config['cutoff_type']
    
    resources:
        mem = "80G",
        time = "2:00:00"
    
    threads: 20
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/call_clusters.py \
            -e {input.expression} \
            -co {input.covariates} \
            -o {output.clusters} \
            --max-cluster-size {params.max_cluster_size} \
            --min-cluster-size {params.min_cluster_size} \
            --min-corr-cutoff {params.min_corr_cutoff} \
            --percent-corr-cutoff {params.percent_corr_cutoff} \
            --cutoff-type {params.cutoff_type} \
            --tissue-id {wildcards.TISSUE}
        """