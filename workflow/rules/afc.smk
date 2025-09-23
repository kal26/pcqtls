"""
Allelic Fold-Change (aFC) Rules

This module computes aFC for top QTLs using `scripts/calculate_afc.py` for both
cluster-based eQTLs and pcQTLs.

Requirements:
- Genotype VCFs must be bgzipped and tabix-indexed (.vcf.gz + .tbi)
- Phenotype BEDs must be bgzipped and tabix-indexed (.bed.gz + .tbi)
- QTL input tables must contain columns: 'pid' (phenotype id), 'sid' (SNP id);
  optional 'sid_chr' and 'sid_pos' speed up VCF fetching.
"""





rule convert_plink_to_vcf:
    """
    Convert PLINK genotype files to VCF format for aFC calculation.
    
    This rule preserves the original sample IDs from the .fam file without
    adding the '0_' prefix that PLINK normally adds during VCF conversion.
    """
    input:
        bed = config['genotype_stem'] + '.bed',
        bim = config['genotype_stem'] + '.bim',
        fam = config['genotype_stem'] + '.fam'
    params:
        genotype_stem = config['genotype_stem'],
    output:
        vcf = config['genotype_stem'] + '.vcf.gz',
        tbi = config['genotype_stem'] + '.vcf.gz.tbi'
    resources:
        mem = "30G",
        time = "4:00:00"
    shell:
        """
        # Convert PLINK to VCF format
        plink --bfile {params.genotype_stem} \
            --recode vcf \
            --out {params.genotype_stem}
        
        # Fix sample IDs by removing '0_' prefix to match phenotype file format
        # This preserves the original sample IDs from the .fam file
        sed 's/\\t0_/\\t/g' {params.genotype_stem}.vcf > {params.genotype_stem}_fixed.vcf
        
        # Compress the fixed VCF file
        bgzip -c {params.genotype_stem}_fixed.vcf > {output.vcf}
        
        # Clean up intermediate files
        rm {params.genotype_stem}.vcf {params.genotype_stem}_fixed.vcf
        
        # Create tabix index
        tabix -p vcf {output.vcf}
        """


rule bgzip_pheno:
    """
    Bgzip and tabix index the phenotype BED for a tissue.
    """
    input:
        bed = config['expression_dir'] + '{TISSUE}.v8.normalized_expression.bed'
    output:
        bed_gz = config['expression_dir'] + '{TISSUE}.v8.normalized_expression.bed.gz',
        tbi = config['expression_dir'] + '{TISSUE}.v8.normalized_expression.bed.gz.tbi'
    resources:
        mem = "30G",
        time = "4:00:00"
    shell:
        """
        bgzip -c {input.bed} > {output.bed_gz}
        tabix -p bed {output.bed_gz}
        """

rule calculate_eqtl_afc:
    """
    Calculate aFC for cluster-based eQTLs for one tissue.

    Expects a QTL table with 'pid' and 'sid' columns, phenotype as cluster
    expression BED, and genotype VCF.
    """
    input:
        vcf = config['genotype_stem'] + '.vcf.gz',
        vcf_tbi = config['genotype_stem'] + '.vcf.gz.tbi',
        pheno = config['expression_dir'] + '{TISSUE}.v8.normalized_expression.bed.gz',
        pheno_tbi = config['expression_dir'] + '{TISSUE}.v8.normalized_expression.bed.gz.tbi',
        qtl = config['eqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie_for_afc.tsv'
    output:
        afc = config['eqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.cluster_genes.afc.txt'
    params:
        code_dir = config['code_dir'],
    resources:
        mem = "30G",
        time = "4:00:00"
    shell:
        """
        python {params.code_dir}/calculate_afc.py \
            --vcf {input.vcf} \
            --pheno {input.pheno} \
            --qtl {input.qtl} \
            --log_xform 1 \
            --log_base 2 \
            --boot 0 \
            --output {output.afc}
        """


rule calculate_pcqtl_afc:
    """
    Calculate aFC for pcQTLs for one tissue.

    Expects a QTL table with 'pid' and 'sid' columns, phenotype as PCs BED,
    and genotype VCF.
    """
    input:
        vcf = config['genotype_stem'] + '.vcf.gz',
        vcf_tbi = config['genotype_stem'] + '.vcf.gz.tbi',
        pheno = config['expression_dir'] + '{TISSUE}.v8.normalized_expression.bed.gz',
        pheno_tbi = config['expression_dir'] + '{TISSUE}.v8.normalized_expression.bed.gz.tbi',
        qtl = config['pcqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs.susie_for_afc.tsv'
    output:
        afc = config['pcqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs.afc.txt'
    params:
        code_dir = config['code_dir'],
    resources:
        mem = "30G",
        time = "4:00:00"
    shell:
        """
        python {params.code_dir}/calculate_afc.py \
            --vcf {input.vcf} \
            --pheno {input.pheno} \
            --qtl {input.qtl} \
            --log_xform 1 \
            --log_base 2 \
            --boot 0 \
            --output {output.afc}
        """


rule make_eqtl_for_afc_from_susie:
    """
    Build eQTL qtl_for_afc table from SuSiE credible set variants.
    Produces columns: pid, sid, sid_chr, sid_pos.
    """
    input:
        susie = config['eqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie.txt'
    output:
        tsv = config['eqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie_for_afc.tsv'
    params:
        code_dir = config['code_dir'],
    resources:
        mem = "30G",
        time = "4:00:00"
    shell:
        """
        python {params.code_dir}/make_qtl_for_afc_from_susie.py \
            --susie {input.susie} \
            --output {output.tsv} \
            --qtl_type eqtl
        """


rule make_pcqtl_for_afc_from_susie:
    """
    Build pcQTL qtl_for_afc table from SuSiE credible set variants.
    Produces columns: pid, sid, sid_chr, sid_pos.
    """
    input:
        susie = config['pcqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs.susie.txt'
    output:
        tsv = config['pcqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs.susie_for_afc.tsv'
    params:
        code_dir = config['code_dir'],
    resources:
        mem = "30G",
        time = "4:00:00"
    shell:
        """
        python {params.code_dir}/make_qtl_for_afc_from_susie.py \
            --susie {input.susie} \
            --output {output.tsv} \
            --qtl_type pcqtl
        """ 