
rule calculate_pcs:
    input:
        clusters = clusters_dir + '{TISSUE}_clusters_all_chr.csv',
        expression = expression_dir + '{TISSUE}.v8.normalized_expression.bed',
        covariates = covariates_dir + '{TISSUE}.v8.covariates.txt'
    conda:
        'tensorqtl_r'
    output:
        pc_output_dir + '{TISSUE}.pc_{PC_ID}.bed'
    shell:"""
        python workflow/scripts/get_pcs.py \
            -cl {input.clusters} \
            -e {input.expression} \
            -co {input.covariates} \
            -o {output} \
            --pc_num {wildcards.PC_ID}
        """

# PCQTLS

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
rule run_pcqtl_cis_nominal:
    input:
        genotypes = genotype_stem + '.fam',
        pcs = pc_output_dir + '{TISSUE}.pc_{PC_ID}.bed',
        covariates = covariates_dir + '{TISSUE}.v8.covariates.txt'
    params:
        genotype_stem = genotype_stem,
        pcqtl_output_dir = pcqtl_output_dir
    resources:
        mem = "30G",
        time = "2:00:00"
    threads: 10
    conda:
        'tensorqtl_r'
    output:
        expand(pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pc_{PC_ID}.cis_qtl_pairs.{CHROM}.parquet', CHROM=chr_list, allow_missing=True)
    shell:"""
        python -m tensorqtl {params.genotype_stem} \
            {input.pcs} \
            {params.pcqtl_output_dir}{wildcards.TISSUE}/{wildcards.TISSUE}.v8.pc_{wildcards.PC_ID} \
            --covariates {input.covariates} \
            --mode cis_nominal
        """

# cis eQTL mapping: permutations (i.e. top variant per phenotype group)
rule run_pcqtl_cis:
    input:
        genotypes = genotype_stem + '.fam',
        pcs = pc_output_dir + '{TISSUE}.pc_{PC_ID}.bed',
        covariates = covariates_dir + '{TISSUE}.v8.covariates.txt'
    params:
        genotype_stem = genotype_stem,
        pcqtl_output_dir = pcqtl_output_dir
    resources:
        mem = "30G",
        time = "2:00:00"
    threads: 10
    conda:
        'tensorqtl_r'
    output:
        pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pc_{PC_ID}.cis_qtl.txt.gz'
    shell:"""
        python -m tensorqtl {params.genotype_stem} \
            {input.pcs} \
            {params.pcqtl_output_dir}{wildcards.TISSUE}/{wildcards.TISSUE}.v8.pc_{wildcards.PC_ID} \
            --covariates {input.covariates} \
            --mode cis
        """

# cis-QTL mapping: conditionally independent QTLs
# This mode maps conditionally independent cis-QTLs using the stepwise regression procedure described in GTEx Consortium, 2017. 
# The output from the permutation step (see map_cis above) is required. 
rule run_pcqtl_cis_independent:
    input:
        genotypes = genotype_stem + '.fam',
        pcs = pc_output_dir + '{TISSUE}.pc_{PC_ID}.bed',
        covariates = covariates_dir + '{TISSUE}.v8.covariates.txt',
        cis_result = pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pc_{PC_ID}.cis_qtl.txt.gz'
    params:
        genotype_stem = genotype_stem, 
        pcqtl_output_dir = pcqtl_output_dir
    resources:
        mem = "30G",
        time = "2:00:00"
    threads: 10
    conda:
        'tensorqtl_r'
    output:
        pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pc_{PC_ID}.cis_independent_qtl.txt.gz'
    shell:"""
        python -m tensorqtl {params.genotype_stem} \
            {input.pcs} \
            {params.pcqtl_output_dir}{wildcards.TISSUE}/{wildcards.TISSUE}.v8.pc_{PC_ID} \
            --covariates {input.covariates} \
            --cis_output {input.cis_results} \
            --mode cis_independent
        """

# cis-QTL mapping: susie credible set summary stats
rule run_pcqtl_susie:
    input:
        genotypes = genotype_stem + '.fam',
        pcs = pc_output_dir + '{TISSUE}.pc_{PC_ID}.bed',
        covariates = covariates_dir + '{TISSUE}.v8.covariates.txt',
    params:
        genotype_stem = genotype_stem,
        pcqtl_output_dir = pcqtl_output_dir
    resources:
        mem = "30G",
        time = "2:00:00"
    threads: 10
    conda:
        'tensorqtl_r'
    output:
        pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pc_{PC_ID}.susie.txt'
    shell:"""
        python workflow/scripts/run_susie.py {params.genotype_stem} \
            {input.pcs} \
            {params.pcqtl_output_dir}{wildcards.TISSUE}/{wildcards.TISSUE}.v8.pc_{wildcards.PC_ID}.susie.txt \
            {input.covariates}
        """