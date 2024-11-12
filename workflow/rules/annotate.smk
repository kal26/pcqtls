# rule annotate_overlap:
#     input:
#         e_susie = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie.txt',
#         pc_susie = pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.susie.txt',
#     resources:
#         mem = "20G"
#     conda:
#         'tensorqtl_r'
#     output:
#         overlap_df = overlap_output_dir + '{TISSUE}.v8.overlap.txt'
#     script:
#         '../scripts/get_overlap.py'

# rule overlap_to_vcf:
#     input:
#         overlap_df = overlap_output_dir + '{TISSUE}.v8.overlap.txt'
#     resources:
#         mem = "5G"
#     conda:
#         'coloc'
#     output:
#         overlap_vcf = annotations_output_dir + '{TISSUE}.v8.leadvars.vcf'
#     script:
#         '../scripts/overlap_to_vcf.py'

rule annotate_pcs:
    input:
        pcs = pc_output_dir + '{TISSUE}.pcs.bed',
        filtered_normed_expression = filtered_expression_output_dir + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed',
    conda:
        'tensorqtl_r'
    resources:
        mem = "20G"
    output:
        pc_annotated = annotations_output_dir + '{TISSUE}/{TISSUE}.v8.pcs_annotated.txt'
    script: 
        '../scripts/snakemake_annotate_pcs.py'

rule susie_to_vcf:
    input:
        e_susie = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie_R.txt',
        pc_susie = pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.susie_R.txt'
    conda:
        'tensorqtl_r'
    resources:
        mem = "20G"
    output:
        susie_vcf = annotations_output_dir + '{TISSUE}/vep/{TISSUE}.v8.susie_R.vcf'
    script:
        '../scripts/snakemake_susie_to_vcf.py'


rule vep_annotate_susie:
    input:
        susie_vcf = annotations_output_dir + '{TISSUE}/vep/{TISSUE}.v8.susie_R.vcf'
    conda:
        'vep'
    resources:
        mem = "30G",
        time = "4:00:00"
    output:
        vep = annotations_output_dir + '{TISSUE}/vep/{TISSUE}.v8.susie_R.vep.vcf',
        vep_summary = annotations_output_dir + '{TISSUE}/vep/{TISSUE}.v8.susie_R.vep.vcf_summary.html'
    shell:"""
    vep -i {input.susie_vcf} \
        -o {output.vep} \
        --cache --vcf -v --force_overwrite \
        --dir_cache data/references/vep/ \
        --no_check_variants_order \
        --nearest transcript \
        --symbol --biotype --regulatory --af --pubmed --verbose
    """

rule qtl_annotate_susie:
    input:
        e_susie = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie_R.txt',
        pc_susie = pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.susie_R.txt',
        vep = annotations_output_dir + '{TISSUE}/vep/{TISSUE}.v8.susie_R.vep.vcf',
        pc_annotated = annotations_output_dir + '{TISSUE}/{TISSUE}.v8.pcs_annotated.txt',
        e_nominal_path_list = expand(eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.{CHROM}.parquet', CHROM=chr_list, allow_missing=True),
        pc_nominal_path_list = expand(pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.cis_qtl_pairs.{CHROM}.parquet', CHROM=chr_list,  allow_missing=True),
    params:
        tissue_id = '{TISSUE}',
        gencode_path=gencode_path,
        full_abc_path = full_abc_path,
        abc_match_path=abc_match_path,
        ctcf_match_path=ctcf_match_path,
        ctcf_dir=ctcf_dir,
        tad_path=tad_path,
        avg_expression_path=avg_expression_path
    conda:
        'tensorqtl_r'
    resources:
        mem = "20G",
        time = "1:00:00"
    output:
        var_annot = annotations_output_dir + '{TISSUE}/{TISSUE}.v8.susie_R_vars.annotated.csv'
    script:
        '../scripts/snakemake_susie_vep_and_annotate.py'





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
        cross_map_path = cross_map_path, 
        tad_path = tad_path
    conda:
        'tensorqtl_r'
    output:
        annotated_clusters = annotations_output_dir + '{TISSUE}/{TISSUE}_clusters_annotated.csv'
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
            --tad_path {params.tad_path} \
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
        cross_map_path = cross_map_path,
        tad_path = tad_path
    conda:
        'tensorqtl_r'
    output:
        annotated_nulls = annotations_output_dir + '{TISSUE}/{TISSUE}_null_{CLUSTER_SIZE}genes_annotated.csv'
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
            --tad_path {params.tad_path} \
            --cluster_size {wildcards.CLUSTER_SIZE} \
            --exclude_cluster_genes {params.exclude_cluster_genes} \
            --distance_matched {params.distance_matched} \
            --verbosity 1
        """

