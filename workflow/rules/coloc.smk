# snp list for all eqtl snps for each cluster (so on a per tissue basis). 
# they will then be filtered down to overlapping gwas variants
# this prevents having to calculate a new ld matrix for each gwas trait
# output id df with cluster_id, snp_list

# rule get_snp_lists:
#     input: 
#         eqtl_pairs = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.{CHROM}.parquet'
#     resources:
#         mem = "3G", 
#         time = "1:00:00"
#     params:
#         tissue_id = '{TISSUE}',
#         chrom = '{CHROM}',
#         output_dir = coloc_output_dir + '{TISSUE}/temp'
#     threads: 10
#     output: 
#         #made_snp_lists = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.made_snp_list.txt'
#         snp_lists = dynamic(coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_{CLUSTER}.snp_list.txt'),
#     script:
#         '../scripts/get_snp_list.py'


# # output is ld matrixes for each cluster
# # maybe store as multiindex?? 
# rule get_ld_cluster:
#     input: 
#         #made_snp_lists = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.made_snp_list.txt',
#         snp_list = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_{CLUSTER}.snp_list.txt',
#         genotypes = genotype_stem + '.fam'
#     resources:
#         mem = "30G", 
#         time = "4:00:00"
#     threads: 10
#     params:
#         snp_list = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_{CLUSTER}.snp_list.txt',
#         genotype_stem = genotype_stem,
#         ld_matrix_stem = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_{CLUSTER}'
#     conda:
#         'coloc'
#     output: 
#         ld_matrix = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_{CLUSTER}.ld'
#     shell:
#         """
#         module load plink

#         plink --bfile {params.genotype_stem} \
#             --extract {params.snp_list} --r2 square \
#             --out {params.ld_matrix_stem}
#         """

# def get_ld_paths(wildcards):
#     clusters = get_clusters(wildcards.TISSUE,wildcards.CHROM)
#     return [f'{coloc_output_dir}{wildcards.TISSUE}/temp/{wildcards.TISSUE}.{wildcards.CHROM}.cluster_{cluster}.ld'.strip() for cluster in clusters]

# def get_snp_list_paths(wildcards):
#     clusters = get_clusters(wildcards.TISSUE,wildcards.CHROM)
#     return [f'{coloc_output_dir}{wildcards.TISSUE}/temp/{wildcards.TISSUE}.{wildcards.CHROM}.cluster_{cluster}.snp_list.txt'.strip() for cluster in clusters]



# per cluster colocalization with gwas
rule run_coloc_chr:
    input:
        eqtl_pairs = expand(eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.{CHROM}.parquet', CHROM=chr_list, allow_missing=True),
        pcqtl_pairs = expand(pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.cis_qtl_pairs.{CHROM}.parquet', CHROM=chr_list, allow_missing=True),
        gwas_meta = gwas_meta,
        gtex_meta = gtex_meta, 
        genotypes = genotype_stem + '.fam',
        gwas_nominal = gwas_folder + '/imputed_{GWAS}.txt.gz',
        annotated_clusters = annotations_output_dir + '{TISSUE}_clusters_annotated.csv'

    resources:
        mem = "100G", 
        time = "20:00:00" ,
    params:
        eqtl_dir_path = eqtl_output_dir + '{TISSUE}',
        pcqtl_dir_path =  pcqtl_output_dir + '{TISSUE}',
        ld_path_head = coloc_output_dir + 'temp/',
        gwas_temp_path_head = coloc_output_dir + '{TISSUE}/temp/',
        genotype_stem = genotype_stem,
        use_susie = '{USE_SUSIE}'

    output:
        coloc_gwas = coloc_output_dir + '{TISSUE}/{TISSUE}.v8.{GWAS}.susie_{USE_SUSIE}.gwas_coloc.txt'
    shell:
        """
        module load r/4.2.2
        module load plink

        Rscript workflow/scripts/draft_coloc.R \
            --eqtl_dir_path {params.eqtl_dir_path} \
            --pcqtl_dir_path {params.pcqtl_dir_path} \
            --gwas_meta {input.gwas_meta} \
            --gtex_meta {input.gtex_meta} \
            --tissue_id {wildcards.TISSUE} \
            --ld_path_head {params.ld_path_head} \
            --gwas_temp_path_head {params.gwas_temp_path_head} \
            --genotype_stem {params.genotype_stem} \
            --gwas_path {input.gwas_nominal} \
            --gwas_id {wildcards.GWAS} \
            --annotated_cluster_path {input.annotated_clusters} \
            --output_path {output.coloc_gwas} \
            --use_susie {params.use_susie} \
        """


# # Define the rule to aggregate coloc results into one file
# rule aggregate_coloc_results_chr:
#     input:
#         coloc_clusters = expand(coloc_output_dir + '{TISSUE}/temp/colocs/{TISSUE}.v8.{CHROM}.cluster_{CLUSTER}.gwas.coloc.txt', zip, CLUSTER=clusters, CHROM=chr_zipped, allow_missing=True)
#     params:
#         coloc_temp_dir = coloc_output_dir + '{TISSUE}/temp/colocs/'
#     output:
#         aggregate_results = coloc_output_dir + '{TISSUE}/{TISSUE}.coloc_results.txt'
#     script:
#         "../script/aggregate_coloc_results.py"

