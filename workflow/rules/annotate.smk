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
        e_susie = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie.txt',
        pc_susie = pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.susie.txt'
    conda:
        'tensorqtl_r'
    resources:
        mem = "20G"
    output:
        susie_vcf = annotations_output_dir + '{TISSUE}/vep/{TISSUE}.v8.susie.vcf'
    script:
        '../scripts/snakemake_susie_to_vcf.py'


rule vep_annotate_susie:
    input:
        susie_vcf = annotations_output_dir + '{TISSUE}/vep/{TISSUE}.v8.susie.vcf'
    conda:
        'vep'
    resources:
        mem = "30G",
        time = "4:00:00"
    output:
        vep = annotations_output_dir + '{TISSUE}/vep/{TISSUE}.v8.susie.vep.vcf',
        vep_summary = annotations_output_dir + '{TISSUE}/vep/{TISSUE}.v8.susie.vep.vcf_summary.html'
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
        e_susie = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie.txt',
        pc_susie = pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.susie.txt',
        vep = annotations_output_dir + '{TISSUE}/vep/{TISSUE}.v8.susie.vep.vcf',
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
        var_annot = annotations_output_dir + '{TISSUE}/{TISSUE}.v8.susie_vars.annotated.csv'
    script:
        '../scripts/snakemake_susie_vep_and_annotate.py'




