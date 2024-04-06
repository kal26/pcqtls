import pandas as pd
import numpy as np
from residualize import calculate_residual

# set params
max_cluster_size = snakemake.params[0],
min_cluster_size = snakemake.params[1],
min_corr_cutoff = snakemake.params[2],
percent_corr_cutoff = snakemake.params[3],

# load data
expression_df = pd.read_csv(snakemake.input[0], sep='\t')
covariates_df = pd.read_csv(snakemake.input[1], sep='\t', index_col=0).T

# residulize the expression 
residal_exp = calculate_residual(expression_df[covariates_df.index], covariates_df, center=True)
residal_exp = pd.DataFrame(residal_exp, columns=covariates_df.index, index=expression_df['gene_id'])

# function for per chr clusters
def get_clusters_chr(chr_id, expression_df, residal_exp):
    # genes on this chr
    chr_gene_ids = expression_df[expression_df['#chr'] == f'chr{chr_id}']['gene_id']
    chr_residual_exp = residal_exp.loc[chr_gene_ids]

    # get correlation
    chr_corr = chr_residual_exp.T.corr()
    # zero the diagonal and everything in the lower half
    chr_corr_zeroed = pd.DataFrame(np.triu(chr_corr.values, k=1), index=chr_gene_ids, columns=chr_gene_ids)


    # iterate through and call the cluster


    # list of transcripts that are already in a cluster
    transcript_blacklist = []

    # list of clusters
    cluster_output = []

    for cluster_size in np.arange(max_cluster_size, min_cluster_size-1, -1):
        # total number of pairs we consizer for this cluster size 
        number_pairs = sum(sum(np.triu(np.ones((cluster_size, cluster_size)), k=1)))

        # each possible start along the genome
        for cluster_start_idx in range(0, len(chr_corr_zeroed)-cluster_size):
            # pull the cluster of this size stating at this index
            cluster_candidate = chr_corr_zeroed.iloc[cluster_start_idx:cluster_start_idx+cluster_size, cluster_start_idx:cluster_start_idx+cluster_size]

            # corr values in upper triangle
            cluster_values = cluster_candidate.values[np.triu_indices(cluster_size, k=1)]
            
            # number of corrs with abs value above cuttoff
            number_corrs_above_cutoff = sum(sum(abs(cluster_candidate.values)>min_corr_cuntoff))

            # transcripts
            cluster_transcripts = cluster_candidate.index.values
            # in blacklist?
            in_blacklist = sum(cluster_candidate.index.isin(transcript_blacklist))

            # record cluster if enough pairs have higher enough abs correlation, and not recoreded before
            if (number_corrs_above_cutoff/number_pairs > percent_corr_cutoff) & (in_blacklist==0):
                # this is the data on each cluster we need
                # N_genes,Transcripts,Genes,Perc_cor,Mean_cor,Mean_pos_cor,Mean_neg_cor,Chromosome,Tissue

                cluster_series = pd.Series({'N_genes':cluster_size,
                            'Transcripts':','.join(cluster_transcripts),
                            'Perc_cor':number_corrs_above_cutoff/number_pairs,
                            'Mean_cor':np.mean(cluster_values),
                            'Mean_pos_cor':np.mean(cluster_values[cluster_values>0]),
                            'Mean_neg_cor':np.mean(cluster_values[cluster_values<0]),
                            'Chromosome':chr_id,
                            'Tissue':tissue_id})
                # record cluster
                cluster_output.append(cluster_series)
                # recored transcript ids so they aren't included in another cluster
                [transcript_blacklist.append(id) for id in cluster_transcripts]

    # make one dataframe and write out
    return pd.DataFrame(cluster_output)

# cycle thorugh all chrs
clusters_all_chr = []
for i in np.arange(1,23,1):
    print(f'Working on chr{i}')
    clusters_all_chr.append(get_clusters_chr(i, expression_df, residal_exp))

# write out
pd.concat(clusters_all_chr).to_csv(snakemake.output[0])