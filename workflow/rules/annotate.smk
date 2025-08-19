"""
Annotation Rules

This module contains rules for functional annotation of QTL variants and clusters,
including VEP annotation, ABC enhancer analysis, CTCF binding, TAD boundaries,
and bidirectional promoter analysis.
"""

rule annotate_clusters:
    """
    Annotate gene expression clusters with functional information.
    
    This rule adds functional annotations to gene clusters including ABC enhancer
    connections, CTCF binding sites, TAD boundaries, and expression statistics.
    """
    input:
        clusters = config['clusters_dir'] + '{TISSUE}.clusters.txt',
        expression = config['filtered_expression_output_dir'] + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed',
        covariates = config['covariates_dir'] + '{TISSUE}.v8.covariates.txt',
        gencode = config['gencode_path'],
        abc = config['full_abc_path'],
        abc_match = config['abc_match_path'],
        ctcf_match = config['ctcf_match_path'],
        ctcf_dir = config['ctcf_dir'],
        paralog = config['paralog_path'],
        go = config['go_path'],
        cross_map = config['cross_map_path'],
        tad = config['tad_path']
    
    output:
        annotated_clusters = config['annotations_output_dir'] + '{TISSUE}/{TISSUE}.clusters.annotated.txt'
    
    resources:
        mem = "20G",
        time = "2:00:00"
    
    threads: 10
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/annotate_clusters.py \
            --tissue-id {wildcards.TISSUE} \
            --clusters {input.clusters} \
            --expression {input.expression} \
            --covariates {input.covariates} \
            --output {output.annotated_clusters} \
            --gencode {input.gencode} \
            --full-abc {input.abc} \
            --abc-match {input.abc_match} \
            --ctcf-match {input.ctcf_match} \
            --ctcf-dir {input.ctcf_dir} \
            --paralog {input.paralog} \
            --go {input.go} \
            --cross-map {input.cross_map} \
            --tad {input.tad}
        """


rule annotate_null_clusters:
    """
    Annotate null gene clusters for statistical comparison.
    
    This rule creates and annotates null clusters by randomly selecting genes
    from the same genomic regions as real clusters for statistical validation.
    """
    input:
        clusters = config['clusters_dir'] + '{TISSUE}.clusters.txt',
        expression = config['filtered_expression_output_dir'] + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed',
        covariates = config['covariates_dir'] + '{TISSUE}.v8.covariates.txt',
        gencode = config['gencode_path'],
        abc = config['full_abc_path'],
        abc_match = config['abc_match_path'],
        ctcf_match = config['ctcf_match_path'],
        ctcf_dir = config['ctcf_dir'],
        paralog = config['paralog_path'],
        go = config['go_path'],
        cross_map = config['cross_map_path'],
        tad = config['tad_path']
    
    output:
        annotated_nulls = config['annotations_output_dir'] + '{TISSUE}/{TISSUE}.null_{CLUSTER_SIZE}genes.annotated.txt'
    
    params:
        cluster_size = "{CLUSTER_SIZE}"
    
    resources:
        mem = "20G",
        time = "2:00:00"
    
    threads: 10
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/annotate_null_clusters.py \
            --tissue-id {wildcards.TISSUE} \
            --clusters {input.clusters} \
            --expression {input.expression} \
            --covariates {input.covariates} \
            --output {output.annotated_nulls} \
            --cluster-size {params.cluster_size} \
            --gencode {input.gencode} \
            --full-abc {input.abc} \
            --abc-match {input.abc_match} \
            --ctcf-match {input.ctcf_match} \
            --ctcf-dir {input.ctcf_dir} \
            --paralog {input.paralog} \
            --go {input.go} \
            --cross-map {input.cross_map} \
            --tad {input.tad}
        """


rule annotate_pcs:
    """
    Annotate principal components with gene expression correlations.
    
    This rule annotates principal components with gene-level correlation statistics
    including slopes and R-squared values for each PC-gene pair within clusters.
    """
    input:
        pcs = config['pc_output_dir'] + '{TISSUE}.pcs.bed',
        filtered_normed_expression = config['filtered_expression_output_dir'] + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed',
    conda:
        'tensorqtl_r'
    resources:
        mem = "20G"
    output:
        pc_annotated = config['annotations_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs_annotated.txt'
    shell:
        """
        python scripts/annotate_pcs.py \
            --pcs {input.pcs} \
            --output {output.pc_annotated} \
            --tissue-id {wildcards.TISSUE}
        """


rule convert_susie_to_vcf:
    """
    Convert SuSiE results to VCF format for VEP annotation.
    
    This rule converts SuSiE fine-mapping results to VCF format to enable
    variant effect prediction using VEP.
    """
    input:
        e_susie_results = config['eqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie.txt',
        pc_susie_results = config['pcqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs.susie.txt'
    
    output:
        vcf = config['annotations_output_dir'] + '{TISSUE}/{TISSUE}.v8.susie_R_vars.vcf'
    
    resources:
        mem = "5G",
        time = "0:30:00"
    
    threads: 1
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/susie_to_vcf.py \
            --e-susie {input.e_susie_results} \
            --pc-susie {input.pc_susie_results} \
            --output {output.vcf}
        """


rule run_vep:
    """
    Run Variant Effect Predictor (VEP) on SuSiE variants.
    
    This rule annotates variants with functional predictions including
    regulatory effects, protein consequences, and population frequencies.
    """
    input:
        vcf = config['annotations_output_dir'] + '{TISSUE}/{TISSUE}.v8.susie_R_vars.vcf'
    
    output:
        vep_results = config['annotations_output_dir'] + '{TISSUE}/{TISSUE}.v8.leadvars.vep.vcf'
    
    resources:
        mem = "30G",
        time = "4:00:00"
    
    threads: 10
    
    conda:
        "vep"
    
    shell:
        """
        vep --input_file {input.vcf} \
            --output_file {output.vep_results} \
            --format vcf \
            --species homo_sapiens \
            --cache \
            --offline \
            --force_overwrite \
            --vcf
        """


rule merge_susie_vep_annotations:
    """
    Merge SuSiE results with VEP annotations and functional data.
    
    This rule combines SuSiE fine-mapping results with VEP annotations and
    additional functional data including ABC enhancers, CTCF binding, and
    TAD boundaries.
    """
    input:
        e_susie_results = config['eqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie_R.txt',
        pc_susie_results = config['pcqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs.susie_R.txt',
        vep_results = config['annotations_output_dir'] + '{TISSUE}/{TISSUE}.v8.leadvars.vep.vcf',
        annotated_pcs = config['annotations_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs_annotated.txt',
        e_nominal = expand(config['eqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.chr{chr}.parquet', chr=range(1, 23), allow_missing=True),
        pc_nominal = expand(config['pcqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs.cis_qtl_pairs.chr{chr}.parquet', chr=range(1, 23), allow_missing=True),
        gencode = config['gencode_path'],
        abc = config['full_abc_path'],
        abc_match = config['abc_match_path'],
        ctcf_match = config['ctcf_match_path'],
        ctcf_dir = config['ctcf_dir'],
        tad = config['tad_path']
    
    output:
        annotated_susie = config['annotations_output_dir'] + '{TISSUE}/{TISSUE}.v8.susie_R_vars.annotated.csv'
    
    resources:
        mem = "20G",
        time = "1:00:00"
    
    threads: 10
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/annotate_susie.py \
            --e-susie {input.e_susie_results} \
            --pc-susie {input.pc_susie_results} \
            --vep {input.vep_results} \
            --annot-pc {input.annotated_pcs} \
            --e-nominal {input.e_nominal} \
            --pc-nominal {input.pc_nominal} \
            --tissue-id {wildcards.TISSUE} \
            --gencode {input.gencode} \
            --full-abc {input.abc} \
            --abc-match {input.abc_match} \
            --ctcf-match {input.ctcf_match} \
            --ctcf-dir {input.ctcf_dir} \
            --tad {input.tad} \
            --output {output.annotated_susie}
        """




