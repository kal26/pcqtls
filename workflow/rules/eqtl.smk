rule filter_expression_clusters:
    input:
        clusters = clusters_dir + '{TISSUE}_clusters_all_chr.csv',
        expression = expression_dir + '{TISSUE}.v8.normalized_expression.bed',
        covariates = covariates_dir + '{TISSUE}.v8.covariates.txt'
    conda:
        'tensorqtl_r'
    output:
        filtered_normed_expression = filtered_expression_output_dir + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed'
    script:
        '../scripts/filter_expression_clusters.py'


# control expression qtls

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
rule run_eqtl_cis_nominal:
    input:
        genotypes = genotype_stem + '.fam',
        expression = filtered_expression_output_dir + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed',
        covariates = covariates_dir + '{TISSUE}.v8.covariates.txt'
    params:
        genotype_stem = genotype_stem,
        eqtl_output_dir = eqtl_output_dir
    resources:
        mem = "30G", 
        time = "2:00:00"
    threads: 10
    conda:
        'tensorqtl_r'
    output:
        expand(eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.{CHROM}.parquet', CHROM=chr_list, allow_missing=True)
    shell:"""
        python -m tensorqtl {params.genotype_stem} \
            {input.expression} \
            {params.eqtl_output_dir}{wildcards.TISSUE}/{wildcards.TISSUE}.v8.cluster_genes \
            --covariates {input.covariates} \
            --mode cis_nominal \
            --maf_threshold .01
        """

# cis eQTL mapping: permutations (i.e. top variant per phenotype group)
rule run_eqtl_cis:
    input:
        genotypes = genotype_stem + '.fam',
        expression = filtered_expression_output_dir + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed',
        covariates = covariates_dir + '{TISSUE}.v8.covariates.txt'
    params:
        genotype_stem = genotype_stem,
        eqtl_output_dir = eqtl_output_dir 
    resources:
        mem = "30G",
        time = "2:00:00"

    threads: 10
    conda:
        'tensorqtl_r'
    output:
        eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl.txt.gz'
    shell:"""
        python -m tensorqtl {params.genotype_stem} \
            {input.expression} \
            {params.eqtl_output_dir}{wildcards.TISSUE}/{wildcards.TISSUE}.v8.cluster_genes \
            --covariates {input.covariates} \
            --mode cis \
            --maf_threshold .01
        """

# cis-QTL mapping: conditionally independent QTLs
# This mode maps conditionally independent cis-QTLs using the stepwise regression procedure described in GTEx Consortium, 2017. 
# The output from the permutation step (see map_cis above) is required. 
rule run_eqtl_cis_independent:
    input:
        genotypes = genotype_stem + '.fam',
        expression = filtered_expression_output_dir + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed',
        covariates = covariates_dir + '{TISSUE}.v8.covariates.txt',
        cis_results = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl.txt.gz'
    params:
        genotype_stem = genotype_stem,
        eqtl_output_dir = eqtl_output_dir
    resources:
        mem = "30G",
        time = "2:00:00"
    threads: 10
    conda:
        'tensorqtl_r'
    output:
       eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_independent_qtl.txt.gz'
    shell:"""
        python -m tensorqtl {params.genotype_stem} \
            {input.expression} \
            {params.eqtl_output_dir}{wildcards.TISSUE}/{wildcards.TISSUE}.v8.cluster_genes \
            --covariates {input.covariates} \
            --cis_output {input.cis_results} \
            --mode cis_independent \
            --maf_threshold .01
        """

# cis-QTL mapping: susie credible set summary stats
rule run_eqtl_susie:
    input:
        genotypes = genotype_stem + '.fam',
        expression = filtered_expression_output_dir + '{TISSUE}.v8.normalized_residualized_expression.cluster_genes.bed',
        covariates = covariates_dir + '{TISSUE}.v8.covariates.txt',
    params:
        genotype_stem = genotype_stem,
        eqtl_output_dir = eqtl_output_dir 
    resources:
        mem = "30G",
        time = "2:00:00"
    threads: 10
    conda:
        'tensorqtl_r'
    output:
       eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie.txt'
    shell:"""
        python workflow/scripts/run_susie.py {params.genotype_stem} \
            {input.expression} \
            {params.eqtl_output_dir}{wildcards.TISSUE}/{wildcards.TISSUE}.v8.cluster_genes.susie.txt \
            {input.covariates}
        """


# # eqtl on all genes instead of cluster genes
# rule get_chr22_expression:
#     input:
#         expression = expression_dir + '{TISSUE}.v8.normalized_expression.bed'
#     conda:
#         'tensorqtl_r'
#     output:
#         'data/processed/chr22_expression/{TISSUE}.v8.normalized_expression.chr22_genes.bed',    
#     script:
#         '../scripts/filter_chr22_expression.py'


# # cis-QTL mapping: susie credible set summary stats
# rule run_eqtl_susie_all:
#     input:
#         genotypes = genotype_stem + '.fam',
#         expression = 'data/processed/chr22_expression/{TISSUE}.v8.normalized_expression.chr22_genes.bed',
#         covariates = covariates_dir + '{TISSUE}.v8.covariates.txt',
#     params:
#         genotype_stem = genotype_stem
#     resources:
#         mem = "30G",
#         time = "2:00:00"
#     threads: 10
#     conda:
#         'tensorqtl_r'
#     output:
#        'output/chr22_eqtl/{TISSUE}/{TISSUE}.v8.chr22_genes.susie.txt'
#     shell:"""
#         python workflow/scripts/run_susie.py {params.genotype_stem} \
#             {input.expression} \
#             output/chr22_eqtl/{wildcards.TISSUE}/{wildcards.TISSUE}.v8.chr22_genes.susie.txt \
#             {input.covariates}
#         """