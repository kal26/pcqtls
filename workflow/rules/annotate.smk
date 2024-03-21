rule annotate_overlap:
    input:
        e_susie = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie.txt',
        pc_susie = pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.susie.txt',
    resources:
        mem = "20G"
    conda:
        'tensorqtl_r'
    output:
        overlap_df = overlap_output_dir + '{TISSUE}.v8.overlap.txt'
    script:
        '../scripts/get_overlap.py'

rule overlap_to_vcf:
    input:
        overlap_df = overlap_output_dir + '{TISSUE}.v8.overlap.txt'
    resources:
        mem = "5G"
    conda:
        'coloc'
    output:
        overlap_vcf = annotations_output_dir + '{TISSUE}.v8.leadvars.vcf'
    script:
        '../scripts/overlap_to_vcf.py'

rule annotate_top_vars:
    input:
        overlap_vcf = annotations_output_dir + '{TISSUE}.v8.leadvars.vcf'
    conda:
        'coloc'
    resources:
        mem = "20G"
    output:
        vep = annotations_output_dir + '{TISSUE}.v8.leadvars.vep.vcf',
        vep_summary = annotations_output_dir + '{TISSUE}.v8.leadvars.vep.vcf_summary.html'
    shell:"""
    vep -i {input.overlap_vcf} \
        -o {output.vep} \
        --cache --vcf -v --force_overwrite \
        --dir_cache data/references/vep/ \
        --no_check_variants_order \
        --nearest transcript \
        --symbol --biotype --regulatory --af --pubmed
    """


rule susie_to_vcf:
    input:
        susie_df = 'output/chr22_eqtl/{TISSUE}/{TISSUE}.v8.chr22_genes.susie.txt'
    resources:
        mem = "10G"
    conda:
        'coloc'
    output:
        overlap_vcf = 'output/chr22_eqtl_annotations/{TISSUE}.v8.chr22_genes.leadvars.vcf'
    script:
        '../scripts/susie_to_vcf.py'

rule annotate_top_vars_control:
    input:
        leadvar_vcf = 'output/chr22_eqtl_annotations/{TISSUE}.v8.chr22_genes.leadvars.vcf'
    conda:
        'coloc'
    resources:
        mem = "20G"
    output:
        vep = 'output/chr22_eqtl_annotations/{TISSUE}.v8.chr22_genes.leadvars.vep.vcf',
        vep_summary = 'output/chr22_eqtl_annotations/{TISSUE}.v8.chr22_genes.leadvars.vep.vcf_summary.html'
    shell:"""
    vep -i {input.leadvar_vcf} \
        -o {output.vep} \
        --cache --vcf -v --force_overwrite \
        --dir_cache data/references/vep/ \
        --no_check_variants_order \
        --nearest transcript \
        --symbol --biotype --regulatory --af --pubmed
    """