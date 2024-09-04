# snp list for all eqtl snps for each cluster (so on a per tissue basis). 
# they will then be filtered down to overlapping gwas variants
# this prevents having to calculate a new ld matrix for each gwas trait
# output id df with cluster_id, snp_list


rule get_cluster_list:
    input:
        clusters = clusters_dir + '{TISSUE}_clusters_all_chr.csv'
    output:
        cluster_list = expand(coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_list.txt', CHROM=chr_list, allow_missing=True)
    params:
        tissue_id = TISSUE
    conda:
        'coloc'
    run:
        import pandas as pd
        cluster_df =  pd.read_csv(snakemake.input[0],index_col=0)
        for idx, row in cluster_df.iterrows():
            cluster_df.loc[idx, 'cluster_id']  = '_'.join([*sorted(row['Transcripts'].split(','))])
        output_head = snakemake.output[0].split('/temp/')[0]
        tissue_id = snakemake.params[0]
        for chrom in cluster_df['Chromosome'].unique():
            cluster_df[cluster_df['Chromosome']==chrom]['cluster_id'].to_csv('f{output_head}/temp/{tissue_id}.{chrom}.cluster_list.txt', index=False, header=False)


rule get_snp_lists:
    input: 
        clusters = clusters_dir + '{TISSUE}_clusters_all_chr.csv'
        eqtl_pairs = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.{CHROM}.parquet',
        cluster_list = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_list.txt'
    resources:
        mem = "3G", 
        time = "1:00:00"
    params:
        cluster = CLUSTER
    threads: 10
    conda:
        'coloc'
    output: 
        snp_lists = expand(coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.{CLUSTER}.snp_list.txt', CLUSTER=open(coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_list.txt').read().strip().split("\n"), allow_missing=True)
    script:
        '../scripts/get_snp_list.py'


# output is ld matrixes for each cluster
# maybe store as multiindex?? 
rule get_ld_cluster:
    input: 
        snp_list = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.{CLUSTER}.snp_list.txt'
        genotype_stem = genotype_stem
    resources:
        mem = "30G", 
        time = "4:00:00"
    threads: 10
    conda:
        'coloc'
    output: 
        ld_matrix = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.{CLUSTER}.ld'
    shell:
        """
        module load plink

        plink --bfile {input.genotype_stem} \
            --extract {input.snp_list} --r2 square \
            --out {output.ld_matrix}
        """

# per cluster colocalization with gwas
rule run_coloc_cluster:
    input:
        eqtl_pairs = eqtl_output_dir + '{TISSUE}/{TISSUE}.v8.cluster_genes.cis_qtl_pairs.{CHROM}.parquet',
        pcqtl_pairs = pcqtl_output_dir + '{TISSUE}/{TISSUE}.v8.pcs.cis_qtl_pairs.{CHROM}.parquet',
        gwas_meta = gwas_meta,
        ld_matrix = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CLUSTER}.ld'
    resources:
        mem = "128G", 
        time = "4:00:00"
    threads: 10
    conda:
        'coloc'
    params:
        gwas_folder = gwas_folder
    output:
        coloc_cluster = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.v8.{CLUSTER}.gwas.coloc.txt'
    shell:
        """
        module load r/4.2.2

        Rscript ../scripts/run_coloc.R --eqtl_path {intput.eqtl_pairs} \
            --pcqtl_path {input.pcqtl_pairs} \
            --gwas_meta {input.gwas_meta} \
            --cluster {wildcards.cluster} \
            --ld_path {input.ld_matrix} \
            --gwas_folder {params.gwas_folder} \
            --output_path {output.coloc_eqtl_gwas_cluster}
        """


# Define the rule to aggregate coloc results into one file
rule aggregate_coloc_results_chr:
    input:
        cluster_list = coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_list.txt'
        coloc_cluster = expand('{CLUSTER}_coloc.txt', CLUSTER=open(coloc_output_dir + '{TISSUE}/temp/{TISSUE}.{CHROM}.cluster_list.txt').read().strip().split("\n"), allow_missing=True)
    output:
        aggregate_results = '{TISSUE}/{TISSUE}.{CHROM}.coloc_results.txt'
    script:
        "../script/aggregate_coloc_results.py"

# maybe this can just be a run, not a full script??