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
    """
    input:
        bed = config['genotype_stem'] + '.bed',
        bim = config['genotype_stem'] + '.bim',
        fam = config['genotype_stem'] + '.fam'
    output:
        vcf = config['genotype_stem'] + '.vcf.gz',
        tbi = config['genotype_stem'] + '.vcf.gz.tbi'
    shell:
        """
        plink2 --bfile {config['genotype_stem']} \
            --recode vcf \
            --out {config['genotype_stem']}
        
        bgzip {config['genotype_stem']}.vcf
        tabix -p vcf {output.vcf}
        """


rule bgzip_eqtl_pheno:
    """
    Bgzip and tabix index the eQTL phenotype BED for a tissue.
    """
    input:
        bed = config['filtered_expression_output_dir'] + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed'
    output:
        bed_gz = config['filtered_expression_output_dir'] + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed.gz',
        tbi = config['filtered_expression_output_dir'] + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed.gz.tbi'
    shell:
        """
        bgzip -c {input.bed} > {output.bed_gz}
        tabix -p bed {output.bed_gz}
        """


rule bgzip_pc_pheno:
    """
    Bgzip and tabix index the pcQTL phenotype BED (PCs) for a tissue.
    """
    input:
        bed = config['pc_output_dir'] + '{TISSUE}.pcs.bed'
    output:
        bed_gz = config['pc_output_dir'] + '{TISSUE}.pcs.bed.gz',
        tbi = config['pc_output_dir'] + '{TISSUE}.pcs.bed.gz.tbi'
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
        pheno = config['filtered_expression_output_dir'] + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed.gz',
        pheno_tbi = config['filtered_expression_output_dir'] + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed.gz.tbi',
        covariates = config['covariates_dir'] + '{TISSUE}.v8.covariates.txt',
        qtl = config['eqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie_for_afc.tsv'
    output:
        afc = config['eqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.cluster_genes.afc.tsv'
    params:
        code_dir = config['code_dir'],
        log_xform = config.get('afc_log_xform', 1),
        chr = config.get('afc_chr_limit', '')
    shell:
        """
        python {params.code_dir}/calculate_afc.py \
            --vcf {input.vcf} \
            --pheno {input.pheno} \
            --qtl {input.qtl} \
            --log_xform {params.log_xform} \
            --boot 0 \
            --cov {input.covariates} \
            {('--chr ' + params.chr) if len(params.chr) > 0 else ''} \
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
        pheno = config['pc_output_dir'] + '{TISSUE}.pcs.bed.gz',
        pheno_tbi = config['pc_output_dir'] + '{TISSUE}.pcs.bed.gz.tbi',
        qtl = config['pcqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs.susie_for_afc.tsv'
    output:
        afc = config['pcqtl_output_dir'] + '{TISSUE}/{TISSUE}.v8.pcs.afc.tsv'
    params:
        code_dir = config['code_dir'],
        log_xform = config.get('afc_log_xform', 1),
        chr = config.get('afc_chr_limit', '')
    shell:
        """
        python {params.code_dir}/calculate_afc.py \
            --vcf {input.vcf} \
            --pheno {input.pheno} \
            --qtl {input.qtl} \
            --log_xform {params.log_xform} \
            --boot 0 \
            --cov {input.covariates} \
            {('--chr ' + params.chr) if len(params.chr) > 0 else ''} \
            --output {output.afc}
        """


rule make_eqtl_qtl_for_afc_from_susie:
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
        min_pip = config.get('susie_min_pip', 0.0)
    shell:
        """
        python {params.code_dir}/make_qtl_for_afc_from_susie.py \
            --susie {input.susie} \
            --output {output.tsv} \
            --min-pip {params.min_pip}
        """


rule make_pcqtl_qtl_for_afc_from_susie:
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
        min_pip = config.get('susie_min_pip', 0.0)
    shell:
        """
        python {params.code_dir}/make_qtl_for_afc_from_susie.py \
            --susie {input.susie} \
            --output {output.tsv} \
            --min-pip {params.min_pip}
        """ 