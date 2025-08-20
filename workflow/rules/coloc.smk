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
        eqtl_nominal = config['eqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.{CHROM}.parquet',
        pcqtl_nominal = config['pcqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs.cis_qtl_pairs.{CHROM}.parquet',
        gtex_meta = config['gtex_meta'],
        clusters = config['clusters_dir'] + '{TISSUE}.clusters.txt'
    
    output:
        coloc_pairs = config['coloc_output_dir'] + 'pairs/{TISSUE}.v8.pairs_coloc.{CHROM}.txt'
    
    params:
        tissue = "{TISSUE}",
        chrom = "{CHROM}",
        ld_path_head = config['coloc_output_dir'] + '{TISSUE}/temp/',
        genotype_stem = config['genotype_stem'],
        coloc_temp_path_head = config['coloc_output_dir'] + '{TISSUE}/temp/',
        code_dir = config['code_dir']
    
    resources:
        mem = "50G",
        time = "10:00:00"
    
    threads: 20
    
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
            --clusters {input.clusters} \
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
        eqtl_susie = config['eqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie.txt',
        pcqtl_susie = config['pcqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs.susie.txt',
        gwas_summary = config['gwas_folder'] + 'imputed_{GWAS}.txt.gz',
        gwas_meta = config['gwas_meta'],
        gtex_meta = config['gtex_meta'],
        clusters = config['clusters_dir'] + '{TISSUE}.clusters.txt',
        expression = config['filtered_expression_output_dir'] + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed'
    
    output:
        gwas_coloc = config['coloc_output_dir'] + 'gwas/{TISSUE}/{TISSUE}.v8.{GWAS}.susie_' + USE_SUSIE + '.gwas_coloc.txt'
    
    params:
        tissue = "{TISSUE}",
        gwas = "{GWAS}",
        use_susie = USE_SUSIE,
        ld_path_head = config['coloc_output_dir'] + '{TISSUE}/temp/',
        genotype_stem = config['genotype_stem'],
        coloc_temp_path_head = config['coloc_output_dir'] + '{TISSUE}/temp/',
        code_dir = config['code_dir']
    
    resources:
        mem = "100G",
        time = "48:00:00"
    
    threads: 40
    
    shell:
        """
        Rscript scripts/coloc_run_gwas.R \
            --code-dir {params.code_dir} \
            --eqtl-dir {input.eqtl_susie} \
            --pcqtl-dir {input.pcqtl_susie} \
            --gwas {input.gwas_summary} \
            --gwas-meta {input.gwas_meta} \
            --gtex-meta {input.gtex_meta} \
            --tissue-id {params.tissue} \
            --ld-path-head {params.ld_path_head} \
            --coloc-temp-path-head {params.coloc_temp_path_head} \
            --genotype-stem {params.genotype_stem} \
            --clusters {input.clusters} \
            --expression {input.expression} \
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
        pair_coloc = expand(config['coloc_output_dir'] + 'pairs/{TISSUE}.v8.pairs_coloc.{CHROM}.txt', CHROM=chr_list, allow_missing=True),
        pc_susie = config['pcqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs.susie_R.txt',
        e_susie = config['eqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie_R.txt'
    
    output:
        qtl_signal_groups = config['coloc_output_dir'] + 'qtl_signal_groups/{TISSUE}.qtl_signal_groups.txt'
    
    params:
        tissue = "{TISSUE}",
        coloc_cutoff = 0.75,
        code_dir = config['code_dir']
    
    resources:
        mem = "20G",
        time = "2:00:00"
    
    threads: 10
    
    shell:
        """
        python {params.code_dir}/group_signals.py \
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
    
    This rule identifies signal groups that include a QTL-GWAS co-localization.
    """
    input:
        gwas_coloc = expand(config['coloc_output_dir'] + 'gwas/{TISSUE}/{TISSUE}.v8.{GWAS}.susie_' + USE_SUSIE + '.gwas_coloc.txt', GWAS=gwas_ids, allow_missing=True),
        pair_coloc = expand(config['coloc_output_dir'] + 'pairs/{TISSUE}.v8.pairs_coloc.{CHROM}.txt', CHROM=chr_list, allow_missing=True)
    
    output:
        gwas_signal_groups = config['coloc_output_dir'] + 'gwas_signal_groups/{TISSUE}.gwas_signal_groups.txt'
    
    params:
        tissue = "{TISSUE}",
        coloc_cutoff = 0.75,
        code_dir = config['code_dir']
    
    resources:
        mem = "30G",
        time = "4:00:00"
    
    threads: 15
    
    shell:
        """
        python {params.code_dir}/group_signals.py \
            --mode gwas \
            --gwas-coloc {input.gwas_coloc} \
            --pair-coloc {input.pair_coloc} \
            --output {output.gwas_signal_groups} \
            --coloc-cutoff {params.coloc_cutoff} \
            --verbose
        """


