# snp list for all eqtl snps for each cluster (so on a per tissue basis). 
# they will then be filtered down to overlapping gwas variants
# this prevents having to calculate a new ld matrix for each gwas trait
# output id df with cluster_id, snp_list

rule get_snp_lists:
    input: 
        eqtl_pairs = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.{CHROM}.parquet',
        cluster_list = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_list.txt'
    resources:
        mem = "3G", 
        time = "1:00:00"
    params:
        tissue_id = '{TISSUE}',
        chrom = '{CHROM}',
        output_dir = coloc_output_dir + '{TISSUE}/temp'
    threads: 10
    conda:
        'tensorqtl_r'
    output: 
        snp_lists = dynamic(coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_{CLUSTER}.snp_list.txt'),
    script:
        '../scripts/get_snp_list.py'
# output is ld matrixes for each cluster
# maybe store as multiindex?? 
rule get_ld_cluster:
    input: 
        snp_list = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_{CLUSTER}.snp_list.txt',
        genotypes = genotype_stem + '.fam'
    resources:
        mem = "30G", 
        time = "4:00:00"
    threads: 10
    params:
        genotype_stem = genotype_stem
    conda:
        'coloc'
    output: 
        ld_matrix = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_{CLUSTER}.ld'
    shell:
        """
        module load plink

        plink --bfile {params.genotype_stem} \
            --extract {input.snp_list} --r2 square \
            --out {output.ld_matrix}
        """

# per cluster colocalization with gwas
rule run_coloc_cluster:
    input:
        eqtl_pairs = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.{CHROM}.parquet',
        pcqtl_pairs = pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.cis_qtl_pairs.{CHROM}.parquet',
        gwas_meta = gwas_meta,
        gtex_meta = gtex_meta, 
        ld_matrix = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_{CLUSTER}.ld', 
        snp_list = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_{CLUSTER}.snp_list.txt',
    resources:
        mem = "128G", 
        time = "4:00:00"
    threads: 10
    conda:
        'coloc'
    params:
        gwas_folder = gwas_folder
    output:
        coloc_cluster = coloc_output_dir + '{TISSUE}/temp/colocs/{TISSUE}.v8.{CHROM}.cluster_{CLUSTER}.gwas_coloc.txt',
    shell:
        """
        module load r/4.2.2

        Rscript ../scripts/run_coloc.R \
            -eqtl_path {input.eqtl_pairs} \
            --pcqtl_path {input.pcqtl_pairs} \
            --gwas_meta {input.gwas_meta} \
            --gtex_meta {input.gtex_meta} \
            --cluster {wildcards.CLUSTER} \
            --tissue {wildcards.TISSUE} \
            --ld_path {input.ld_matrix} \
            --snp_list_path {input.snp_list} \
            --gwas_folder {params.gwas_folder} \
            --output_path {output.coloc_cluster} \
        """


# # Define the rule to aggregate coloc results into one file
# rule aggregate_coloc_results_chr:
#     input:
#         cluster_list = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_list.txt',
#         coloc_cluster = expand(coloc_output_dir + '{TISSUE}/temp/colocs/{TISSUE}.v8.{CHROM}.cluster_{CLUSTER}.gwas.coloc.txt', CLUSTER=clusters, allow_missing=True)
#     params:
#         coloc_temp_dir = coloc_output_dir + '{TISSUE}/temp/colocs/'
#     output:
#         aggregate_results = coloc_output_dir + '{TISSUE}/{TISSUE}.{CHROM}.coloc_results.txt'
#     script:
#         "../script/aggregate_coloc_results.py"

