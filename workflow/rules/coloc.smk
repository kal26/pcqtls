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
rule run_coloc_gwas:
    input:
        eqtl_pairs = expand(eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.{CHROM}.parquet', CHROM=chr_list, allow_missing=True),
        pcqtl_pairs = expand(pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.cis_qtl_pairs.{CHROM}.parquet', CHROM=chr_list, allow_missing=True),
        gwas_meta = gwas_meta,
        gtex_meta = gtex_meta, 
        genotypes = genotype_stem + '.fam',
        gwas_nominal = gwas_folder + '/imputed_{GWAS}.txt.gz',
        annotated_clusters = annotations_output_dir + '{TISSUE}/{TISSUE}_clusters_annotated.csv'

    resources:
        mem = "200G", 
        time = "48:00:00" ,
    params:
        eqtl_dir_path = eqtl_output_dir + '{TISSUE}',
        pcqtl_dir_path =  pcqtl_output_dir + '{TISSUE}',
        ld_path_head = coloc_output_dir + 'temp/all_tissues/',
        coloc_temp_path_head = coloc_output_dir + 'temp/{TISSUE}/',
        genotype_stem = genotype_stem,
        use_susie = '{USE_SUSIE}'
    conda:
        "/oak/stanford/groups/smontgom/klawren/micromamba/envs/susie_r"
    output:
        coloc_gwas = coloc_output_dir + 'gwas/{TISSUE}/{TISSUE}.v8.{GWAS}.susie_{USE_SUSIE}.gwas_coloc.txt'
    shell:
        """
        Rscript workflow/scripts/coloc_run_gwas.R \
            --eqtl_dir_path {params.eqtl_dir_path} \
            --pcqtl_dir_path {params.pcqtl_dir_path} \
            --gwas_meta {input.gwas_meta} \
            --gtex_meta {input.gtex_meta} \
            --tissue_id {wildcards.TISSUE} \
            --ld_path_head {params.ld_path_head} \
            --coloc_temp_path_head {params.coloc_temp_path_head} \
            --genotype_stem {params.genotype_stem} \
            --gwas_path {input.gwas_nominal} \
            --gwas_id {wildcards.GWAS} \
            --annotated_cluster_path {input.annotated_clusters} \
            --output_path {output.coloc_gwas} \
            --use_susie {params.use_susie} \
        """

        # module load r/4.2.2
        # module load plink


# per cluster colocalization with gwas
rule run_coloc_pairs:
    input:
        eqtl_pairs = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.{CHROM}.parquet', 
        pcqtl_pairs = pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.cis_qtl_pairs.{CHROM}.parquet',
        gtex_meta = gtex_meta, 
        genotypes = genotype_stem + '.fam',
        annotated_clusters = annotations_output_dir + '{TISSUE}/{TISSUE}_clusters_annotated.csv'
    resources:
        mem = "200G", 
        time = "10:00:00" ,
    params:
        eqtl_dir_path = eqtl_output_dir + '{TISSUE}',
        pcqtl_dir_path =  pcqtl_output_dir + '{TISSUE}',
        ld_path_head = coloc_output_dir + 'temp/all_tissues/',
        genotype_stem = genotype_stem,
        coloc_temp_path_head = coloc_output_dir + 'temp/{TISSUE}/',
    conda:
        "/oak/stanford/groups/smontgom/klawren/micromamba/envs/susie_r"
    output:
        coloc_pairs = coloc_output_dir + 'pairs/{TISSUE}.v8.pairs_coloc.{CHROM}.txt'
    shell:
        """
        Rscript workflow/scripts/coloc_run_pairs.R \
            --eqtl_dir_path {params.eqtl_dir_path} \
            --pcqtl_dir_path {params.pcqtl_dir_path} \
            --gtex_meta {input.gtex_meta} \
            --tissue_id {wildcards.TISSUE} \
            --ld_path_head {params.ld_path_head} \
            --genotype_stem {params.genotype_stem} \
            --annotated_cluster_path {input.annotated_clusters} \
            --output_path {output.coloc_pairs} \
            --chr_id {wildcards.CHROM}
            --coloc_temp_path_head {params.coloc_temp_path_head}
        """


# combine eqtl results into a similar format to tensorqtl 
rule gather_eqtl_susie:
    input:
        eqtl_pairs = expand(eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.{CHROM}.parquet', CHROM=chr_list, allow_missing=True),
        coloc_pairs = expand(coloc_output_dir + 'pairs/{TISSUE}.v8.pairs_coloc.{CHROM}.txt', CHROM=chr_list, allow_missing=True) # indicates all the susies have been run
    resources:
        mem = "50G", 
        time = "4:00:00" ,
    params:
        eqtl_dir_path = eqtl_output_dir + '{TISSUE}',
        coloc_temp_path_head = coloc_output_dir + 'temp/{TISSUE}/'
    conda:
        "/oak/stanford/groups/smontgom/klawren/micromamba/envs/susie_r"
    output:
        eqtl_susie_pairs = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.susie_R.txt'
    shell:
        """
        Rscript workflow/scripts/combine_RDS_susie.R \
            --qtl_dir_path {params.eqtl_dir_path} \
            --output_path {output.eqtl_susie_pairs} \
            --coloc_temp_path_head {params.coloc_temp_path_head} \
            --qtl_type eqtl \
            --tissue_id {wildcards.TISSUE} \
        """


# combine pcqtl results into a similar format to tensorqtl 
rule gather_pcqtl_susie:
    input:
        pcqtl_pairs = expand(pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.cis_qtl_pairs.{CHROM}.parquet', CHROM=chr_list, allow_missing=True),
        coloc_pairs = expand(coloc_output_dir + 'pairs/{TISSUE}.v8.pairs_coloc.{CHROM}.txt', CHROM=chr_list, allow_missing=True) # indicates all the susies have been run
    resources:
        mem = "50G", 
        time = "4:00:00" ,
    params:
        pcqtl_dir_path =  pcqtl_output_dir + '{TISSUE}',
        coloc_temp_path_head = coloc_output_dir + 'temp/{TISSUE}/'
    conda:
        "/oak/stanford/groups/smontgom/klawren/micromamba/envs/susie_r"
    output:
        pcqtl_susie_pairs = pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.susie_R.txt'
    shell:
        """
        Rscript workflow/scripts/combine_RDS_susie.R \
            --qtl_dir_path {params.pcqtl_dir_path} \
            --output_path {output.pcqtl_susie_pairs} \
            --coloc_temp_path_head {params.coloc_temp_path_head} \
            --qtl_type pcqtl \
            --tissue_id {wildcards.TISSUE} \
        """


