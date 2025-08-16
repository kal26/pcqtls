"""
Co-localization Analysis Rules

This module contains rules for co-localization analysis between QTL associations
and GWAS summary statistics to identify shared genetic signals.
"""




rule run_coloc_pairs:
    """
    Run co-localization analysis between eQTL and pcQTL pairs.
    
    This rule performs co-localization analysis between eQTL and pcQTL associations
    to identify shared genetic signals within the same tissue.
    """
    input:
        eqtl_nominal = f"{config['eqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.cluster_genes.cis_qtl_pairs.{{CHROM}}.parquet",
        pcqtl_nominal = f"{config['pcqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.pcs.cis_qtl_pairs.{{CHROM}}.parquet",
        gtex_meta = config['gtex_meta_path'],
        annotated_cluster = f"{config['annotations_output_dir']}{{TISSUE}}/{{TISSUE}}.clusters.annotated.txt"
    
    output:
        coloc_pairs = f"{config['coloc_output_dir']}pairs/{{TISSUE}}.v8.pairs_coloc.{{CHROM}}.txt"
    
    params:
        tissue = "{TISSUE}",
        chrom = "{CHROM}",
        ld_path_head = f"{config['coloc_output_dir']}{{TISSUE}}/temp/",
        genotype_stem = config['genotype_stem'],
        coloc_temp_path_head = f"{config['coloc_output_dir']}{{TISSUE}}/temp/",
        code_dir = config['code_dir']
    
    resources:
        mem = "50G",
        time = "10:00:00"
    
    threads: 20
    
    conda:
        "coloc"
    
    shell:
        """
        Rscript scripts/coloc_run_pairs.R \
            --code-dir {params.code_dir} \
            --eqtl-dir {input.eqtl_nominal} \
            --pcqtl-dir {input.pcqtl_nominal} \
            --gtex-meta {input.gtex_meta} \
            --tissue-id {params.tissue} \
            --chr-id {params.chrom} \
            --ld-path-head {params.ld_path_head} \
            --genotype-stem {params.genotype_stem} \
            --annotated-cluster {input.annotated_cluster} \
            --output {output.coloc_pairs} \
            --coloc-temp-path-head {params.coloc_temp_path_head}
        """


rule run_gwas_coloc:
    """
    Run co-localization analysis with GWAS summary statistics.
    
    This rule performs co-localization analysis between QTL associations and
    GWAS summary statistics to identify shared genetic signals.
    """
    input:
        eqtl_susie = f"{config['eqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.cluster_genes.susie.txt",
        pcqtl_susie = f"{config['pcqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.pcs.susie.txt",
        gwas_summary = f"{config['gwas_folder']}{{GWAS}}.txt",
        gwas_meta = config['gwas_meta_path'],
        gtex_meta = config['gtex_meta_path'],
        annotated_cluster = f"{config['annotations_output_dir']}{{TISSUE}}/{{TISSUE}}.clusters.annotated.txt"
    
    output:
        gwas_coloc = f"{config['coloc_output_dir']}gwas/{{TISSUE}}/{{TISSUE}}.v8.{{GWAS}}.susie_{USE_SUSIE}.gwas_coloc.txt"
    
    params:
        tissue = "{TISSUE}",
        gwas = "{GWAS}",
        use_susie = USE_SUSIE,
        ld_path_head = f"{config['coloc_output_dir']}{{TISSUE}}/temp/",
        genotype_stem = config['genotype_stem'],
        coloc_temp_path_head = f"{config['coloc_output_dir']}{{TISSUE}}/temp/",
        code_dir = config['code_dir']
    
    resources:
        mem = "100G",
        time = "48:00:00"
    
    threads: 40
    
    conda:
        "coloc"
    
    shell:
        """
        Rscript scripts/coloc_run_gwas.R \
            --code-dir {params.code_dir} \
            --eqtl-dir {input.eqtl_susie} \
            --pcqtl-dir {input.pcqtl_susie} \
            --gwas-meta {input.gwas_meta} \
            --gtex-meta {input.gtex_meta} \
            --tissue-id {params.tissue} \
            --ld-path-head {params.ld_path_head} \
            --coloc-temp-path-head {params.coloc_temp_path_head} \
            --genotype-stem {params.genotype_stem} \
            --gwas {input.gwas_summary} \
            --gwas-id {params.gwas} \
            --annotated-cluster {input.annotated_cluster} \
            --output {output.gwas_coloc} \
            --use-susie {params.use_susie}
        """


rule group_qtl_signals:
    """
    Group QTL signals from co-localization results.
    
    This rule identifies signal groups by connecting co-localized credible sets
    using network analysis to identify connected components.
    """
    input:
        pair_coloc = f"{config['coloc_output_dir']}pairs/{{TISSUE}}.v8.pairs_coloc.txt",
        pc_susie = f"{config['pcqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.pcs.susie_R.txt",
        e_susie = f"{config['eqtl_output_dir']}{{TISSUE}}/{{TISSUE}}.v8.cluster_genes.susie_R.txt"
    
    output:
        qtl_signal_groups = f"{config['coloc_output_dir']}qtl_signal_groups/{{TISSUE}}.qtl_signal_groups.txt"
    
    params:
        tissue = "{TISSUE}",
        coloc_cutoff = 0.75
    
    resources:
        mem = "20G",
        time = "2:00:00"
    
    threads: 10
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/group_signals.py \
            --mode qtl \
            --pair-coloc {input.pair_coloc} \
            --pc-susie {input.pc_susie} \
            --e-susie {input.e_susie} \
            --output {output.qtl_signal_groups} \
            --coloc-cutoff {params.coloc_cutoff} \
            --verbose
        """


rule group_gwas_signals:
    """
    Group GWAS signals from co-localization results.
    
    This rule identifies signal groups that include both QTL-QTL and QTL-GWAS
    co-localizations using network analysis.
    """
    input:
        gwas_coloc = f"{config['coloc_output_dir']}gwas/{{TISSUE}}/{{TISSUE}}.v8.{{GWAS}}.susie_{USE_SUSIE}.gwas_coloc.txt",
        pair_coloc = f"{config['coloc_output_dir']}pairs/{{TISSUE}}.v8.pairs_coloc.txt"
    
    output:
        gwas_signal_groups = f"{config['coloc_output_dir']}gwas_signal_groups/{{TISSUE}}.{{GWAS}}.gwas_signal_groups.txt"
    
    params:
        tissue = "{TISSUE}",
        gwas = "{GWAS}",
        coloc_cutoff = 0.75
    
    resources:
        mem = "30G",
        time = "4:00:00"
    
    threads: 15
    
    conda:
        "tensorqtl_r"
    
    shell:
        """
        python scripts/group_signals.py \
            --mode gwas \
            --gwas-coloc {input.gwas_coloc} \
            --pair-coloc {input.pair_coloc} \
            --output {output.gwas_signal_groups} \
            --coloc-cutoff {params.coloc_cutoff} \
            --verbose
        """


