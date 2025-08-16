"""
pcQTL Analysis Rules

This module contains rules for principal component QTL (pcQTL) analysis including
PC calculation, nominal QTL mapping, permutation analysis, and SuSiE fine-mapping.
"""

rule calculate_pcs:
    """
    Calculate principal components from gene expression clusters.
    
    This rule generates principal components from co-expressed gene clusters
    for downstream pcQTL analysis.
    """
    input:
        clusters = f"{config['clusters_dir']}{{TISSUE}}.clusters.txt",
        filtered_expression = f"{config['filtered_expression_output_dir']}{{TISSUE}}.v8.normalized_residualized_expression.cluster_genes.bed",
        covariates = f"{config['covariates_dir']}{{TISSUE}}.v8.covariates.txt"
    
    output:
        pcs = f"{config['pc_output_dir']}{{TISSUE}}.pcs.bed"
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/get_pcs.py \
            -cl {input.clusters} \
            -e {input.filtered_expression} \
            -co {input.covariates} \
            -o {output.pcs}
        """


rule run_pcqtl_cis_nominal:
    """
    Run nominal pcQTL analysis for all variant-phenotype pairs.
    
    This rule performs cis-QTL mapping using TensorQTL to identify associations
    between genetic variants and principal components derived from gene clusters.
    """
    input:
        genotypes = f"{config['genotype_stem']}.fam",
        pcs = f"{config['pc_output_dir']}{{TISSUE}}.pcs.bed",
        covariates = f"{config['covariates_dir']}{{TISSUE}}.v8.covariates.txt"
    
    output:
        expand(f"{config['pcqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.pcs.cis_qtl_pairs.{{CHROM}}.parquet", 
               CHROM=chr_list, allow_missing=True)
    
    params:
        genotype_stem = config['genotype_stem'],
        pcqtl_output_dir = config['pcqtl_output_dir'],
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
            --expression {input.pcs} \
            --covariates {input.covariates} \
            --output-dir {params.pcqtl_output_dir} \
            --tissue-id {params.tissue} \
            --phenotype-type .v8.pcs
        """


rule run_pcqtl_cis:
    """
    Run pcQTL permutation analysis to identify top variants per phenotype.
    
    This rule performs permutation-based QTL analysis to identify the most
    significant variant for each principal component phenotype.
    """
    input:
        genotypes = f"{config['genotype_stem']}.fam",
        pcs = f"{config['pc_output_dir']}{{TISSUE}}.pcs.bed",
        covariates = f"{config['covariates_dir']}{{TISSUE}}.v8.covariates.txt"
    
    output:
        pcqtl_results = f"{config['pcqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.pcs.cis_qtl.txt.gz"
    
    params:
        genotype_stem = config['genotype_stem'],
        pcqtl_output_dir = config['pcqtl_output_dir']
    
    resources:
        mem = "30G",
        time = "4:00:00"
    
    threads: 10
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python -m tensorqtl {params.genotype_stem} \
            {input.pcs} \
            {params.pcqtl_output_dir}{wildcards.TISSUE}/{wildcards.TISSUE}.v8.pcs \
            --covariates {input.covariates} \
            --mode cis
        """


rule run_pcqtl_cis_independent:
    """
    Run conditionally independent pcQTL analysis.
    
    This rule identifies conditionally independent cis-QTLs using stepwise
    regression, building on the permutation analysis results.
    """
    input:
        genotypes = f"{config['genotype_stem']}.fam",
        pcs = f"{config['pc_output_dir']}{{TISSUE}}.pcs.bed",
        covariates = f"{config['covariates_dir']}{{TISSUE}}.v8.covariates.txt",
        cis_results = f"{config['pcqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.pcs.cis_qtl.txt.gz"
    
    output:
        independent_qtls = f"{config['pcqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.pcs.cis_independent_qtl.txt.gz"
    
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
            --phenotype {input.pcs} \
            --covariates {input.covariates} \
            --cis-results {input.cis_results} \
            --output {output.independent_qtls} \
            --tissue-id {params.tissue}
        """





rule run_pcqtl_susie:
    """
    Run SuSiE fine-mapping for pcQTL credible sets.
    
    This rule performs SuSiE fine-mapping to identify credible sets of variants
    that are likely to contain the causal variant for each pcQTL association.
    """
    input:
        genotypes = f"{config['genotype_stem']}.fam",
        pcs = f"{config['pc_output_dir']}{{TISSUE}}.pcs.bed",
        covariates = f"{config['covariates_dir']}{{TISSUE}}.v8.covariates.txt"
    
    output:
        susie_results = f"{config['pcqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.pcs.susie.txt"
    
    params:
        genotype_stem = config['genotype_stem'],
        pcqtl_output_dir = config['pcqtl_output_dir']
    
    resources:
        mem = "30G",
        time = "4:00:00"
    
    threads: 10
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/run_susie.py \
            --genotype-stem {params.genotype_stem} \
            --expression {input.pcs} \
            --covariates {input.covariates} \
            --output {output.susie_results}
        """