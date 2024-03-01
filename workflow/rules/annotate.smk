rule annotate_overlap:
    input:
        e_susie = 'output/control_eqtl/{TISSUE}/{TISSUE}.v8.cluster_genes.susie.txt',
        pc_1_susie = 'output/pcqtl/{TISSUE}/{TISSUE}.v8.pc_1.susie.txt',
        pc_2_susie = 'output/pcqtl/{TISSUE}/{TISSUE}.v8.pc_2.susie.txt'
    resources:
        mem = 30
    threads: 10
    conda:
        'tensorqtl_r'
    output:
        overlap_df = 'output/overlap/{TISSUE}.v8.overlap.txt'
    script:
        '../scripts/get_overlap.py'

rule overlap_to_vcf:
    input:
        overlap_df = 'output/overlap/{TISSUE}.v8.overlap.txt'
    resources:
        mem = 10
    conda:
        'coloc'
    output:
        overlap_vcf = 'output/annotations/{TISSUE}.v8.leadvars.vcf'
    script:
        '../scripts/overlap_to_vcf.py'

rule annotate_top_vars:
    input:
        overlap_vcf = 'output/annotations/{TISSUE}.v8.leadvars.vcf'
    conda:
        'coloc'
    resources:
        mem = 10
    output:
        vep = 'output/annotations/{TISSUE}.v8.leadvars.vep.vcf',
        vep_summary = 'output/annotations/{TISSUE}.v8.leadvars.vep.vcf_summary.html'
    shell:"""
    vep -i {input.overlap_vcf} \
        -o {output.vep} \
        --cache --vcf -v --force_overwrite \
        --dir_cache data/references/vep/ \
        --no_check_variants_order \
        --nearest transcript \
        --symbol --biotype --regulatory --af --pubmed
    """