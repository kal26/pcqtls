"""
eQTL Analysis Rules

This module contains rules for eQTL analysis including expression filtering,
nominal QTL mapping, permutation analysis, and SuSiE fine-mapping.
"""

rule filter_expression_clusters:
    """
    Filter and residualize expression data for cluster-based eQTL analysis.
    
    This rule processes expression data to create cluster-specific expression
    files for downstream QTL analysis.
    """
    input:
        clusters = f"{config['clusters_dir']}{{TISSUE}}.clusters.txt",
        expression = f"{config['expression_dir']}{{TISSUE}}.v8.normalized_expression.bed",
        covariates = f"{config['covariates_dir']}{{TISSUE}}.v8.covariates.txt"
    
    output:
        filtered_expression = f"{config['filtered_expression_output_dir']}{{TISSUE}}.v8.normalized_residualized_expression.cluster_genes.bed"
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/filter_expression_clusters.py \
            --clusters {input.clusters} \
            --expression {input.expression} \
            --covariates {input.covariates} \
            --output {output.filtered_expression}
        """


rule run_eqtl_cis_nominal:
    """
    Run nominal eQTL analysis for all variant-phenotype pairs.
    
    This rule performs cis-QTL mapping using TensorQTL to identify associations
    between genetic variants and individual gene expression levels within clusters.
    """
    input:
        genotypes = f"{config['genotype_stem']}.fam",
        expression = f"{config['filtered_expression_output_dir']}{{TISSUE}}.v8.normalized_residualized_expression.cluster_genes.bed",
        covariates = f"{config['covariates_dir']}{{TISSUE}}.v8.covariates.txt"
    
    output:
        expand(f"{config['eqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.cluster_genes.cis_qtl_pairs.{{CHROM}}.parquet", 
               CHROM=chr_list, allow_missing=True)
    
    params:
        genotype_stem = config['genotype_stem'],
        eqtl_output_dir = config['eqtl_output_dir'],
        tissue = "{TISSUE}"
    
    resources:
        mem = "30G",
        time = "4:00:00"
    
    threads: 10
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/run_qtl_nominal.py \
            --genotype-stem {params.genotype_stem} \
            --expression {input.expression} \
            --covariates {input.covariates} \
            --output-dir {params.eqtl_output_dir} \
            --tissue-id {params.tissue} \
            --phenotype-type .v8.cluster_genes
        """


rule run_eqtl_cis:
    """
    Run eQTL permutation analysis to identify top variants per phenotype.
    
    This rule performs permutation-based QTL analysis to identify the most
    significant variant for each gene expression phenotype.
    """
    input:
        genotypes = f"{config['genotype_stem']}.fam",
        expression = f"{config['filtered_expression_output_dir']}{{TISSUE}}.v8.normalized_residualized_expression.cluster_genes.bed",
        covariates = f"{config['covariates_dir']}{{TISSUE}}.v8.covariates.txt"
    
    output:
        eqtl_results = f"{config['eqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.cluster_genes.cis_qtl.txt.gz"
    
    params:
        genotype_stem = config['genotype_stem'],
        eqtl_output_dir = config['eqtl_output_dir']
    
    resources:
        mem = "30G",
        time = "4:00:00"
    
    threads: 10
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python -m tensorqtl {params.genotype_stem} \
            {input.expression} \
            {params.eqtl_output_dir}{wildcards.TISSUE}/{wildcards.TISSUE}.v8.cluster_genes \
            --covariates {input.covariates} \
            --mode cis
        """


rule run_eqtl_cis_independent:
    """
    Run conditionally independent eQTL analysis.
    
    This rule identifies conditionally independent cis-QTLs using stepwise
    regression, building on the permutation analysis results.
    """
    input:
        genotypes = f"{config['genotype_stem']}.fam",
        expression = f"{config['filtered_expression_output_dir']}{{TISSUE}}.v8.normalized_residualized_expression.cluster_genes.bed",
        covariates = f"{config['covariates_dir']}{{TISSUE}}.v8.covariates.txt",
        cis_results = f"{config['eqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.cluster_genes.cis_qtl.txt.gz"
    
    output:
        independent_qtls = f"{config['eqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.cluster_genes.cis_independent_qtl.txt.gz"
    
    params:
        genotype_stem = config['genotype_stem'],
        tissue = "{TISSUE}"
    
    resources:
        mem = "60G",
        time = "6:00:00"
    
    threads: 20
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/run_qtl_permutations.py \
            --genotype-stem {params.genotype_stem} \
            --phenotype {input.expression} \
            --covariates {input.covariates} \
            --cis-results {input.cis_results} \
            --output {output.independent_qtls} \
            --tissue-id {params.tissue}
        """


rule run_eqtl_susie:
    """
    Run SuSiE fine-mapping for eQTL credible sets.
    
    This rule performs SuSiE fine-mapping to identify credible sets of variants
    that are likely to contain the causal variant for each eQTL association.
    """
    input:
        genotypes = f"{config['genotype_stem']}.fam",
        expression = f"{config['filtered_expression_output_dir']}{{TISSUE}}.v8.normalized_residualized_expression.cluster_genes.bed",
        covariates = f"{config['covariates_dir']}{{TISSUE}}.v8.covariates.txt"
    
    output:
        susie_results = f"{config['eqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.cluster_genes.susie.txt"
    
    params:
        genotype_stem = config['genotype_stem'],
        eqtl_output_dir = config['eqtl_output_dir']
    
    resources:
        mem = "30G",
        time = "4:00:00"
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/run_susie.py \
            --genotype-stem {params.genotype_stem} \
            --expression {input.expression} \
            --covariates {input.covariates} \
            --output {output.susie_results}
        """