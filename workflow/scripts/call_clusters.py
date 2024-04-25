print('running call clusters')
import pandas as pd
import numpy as np
from residualize import calculate_residual
from scipy.stats import spearmanr

print('starting')

# set params
max_cluster_size = snakemake.params[0]
min_cluster_size = snakemake.params[1]
min_corr_cutoff = snakemake.params[2]
percent_corr_cutoff = snakemake.params[3]
cutoff_type = snakemake.params[4]
tissue_id = snakemake.params[5]


# load data
expression_df = pd.read_csv(snakemake.input[0], sep='\t')
covariates_df = pd.read_csv(snakemake.input[1], sep='\t', index_col=0).T
print('loaded data')


# residulize the expression 
residal_exp = calculate_residual(expression_df[covariates_df.index], covariates_df, center=True)
residal_exp = pd.DataFrame(residal_exp, columns=covariates_df.index, index=expression_df['gene_id'])
print('residualized expression')

# calculate total number of pairs considered for bonferroni correction
total_pairs = 0
for i in np.arange(1,23,1):
    chr_gene_ids = expression_df[expression_df['#chr'] == f'chr{i}']['gene_id']
    upper_corner_idxs = np.triu(np.ones(len(chr_gene_ids)), k=1)
    excluded_cluster_size_idxs = np.triu(np.ones(len(chr_gene_ids)), k=max_cluster_size)
    total_pairs += upper_corner_idxs.sum()  - excluded_cluster_size_idxs.sum()


# function to avoid calculating far-off diagonal corrs that won't be needed for any clusters
def get_corr_diag(df, max_cluster_size):
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    corrs = dfcols.transpose().join(dfcols, how='outer')
    for i in range(len(df.columns)):
            gene_id_1 = df.columns[i]
            for j in range(len(df.columns)):
                gene_id_2 = df.columns[j]
                if (i > j) & (i < j + max_cluster_size):
                    result = spearmanr(df[gene_id_1], df[gene_id_2])
                else: 
                    result = (np.NAN, np.NAN)
                corrs[gene_id_1][gene_id_2] = result[0]
                pvalues[gene_id_1][gene_id_2] = result[1]
    return corrs, pvalues


# get expression clusters on a per-chrom basis
def get_clusters_chr(chr_id, expression_df, residal_exp, total_pairs):
    # genes on this chr
    chr_gene_ids = expression_df[expression_df['#chr'] == f'chr{chr_id}']['gene_id']
    chr_residual_exp = residal_exp.loc[chr_gene_ids]

    # get correlation and pvalues
    chr_corr, chr_pvalue = get_corr_diag(chr_residual_exp.T, max_cluster_size)
    print('calcuated correlations')

    # iterate through and call the cluster
    # list of transcripts that are already in a cluster
    transcript_blacklist = []

    # list of clusters
    cluster_output = []

    for cluster_size in np.arange(max_cluster_size, min_cluster_size-1, -1):
        # total number of pairs we consizer for this cluster size 
        number_pairs = sum(sum(np.triu(np.ones((cluster_size, cluster_size)), k=1)))

        # each possible start along the genome
        for cluster_start_idx in range(0, len(chr_corr)-cluster_size):
            # pull the cluster of this size stating at this index
            
            cluster_candidate = chr_corr.iloc[cluster_start_idx:cluster_start_idx+cluster_size, cluster_start_idx:cluster_start_idx+cluster_size]
            # corr values in upper triangle
            cluster_values = cluster_candidate.values[np.triu_indices(cluster_size, k=1)]
            
            # number of corrs with abs value above cuttoff or p value below cutoff
            if cutoff_type == 'value':
                number_corrs_above_cutoff = sum(sum(abs(cluster_candidate.values)>min_corr_cuntoff))
            elif cutoff_type == 'pvalue':
                # calculate p values below a bonferonni 0.05
                cluster_candidate_pvalues = chr_pvalue.iloc[cluster_start_idx:cluster_start_idx+cluster_size, cluster_start_idx:cluster_start_idx+cluster_size]
                # only take top corner 
                cluster_pvalues = cluster_candidate_pvalues.values[np.triu_indices(cluster_size, k=1)]
                # pvalues values in upper triangle less than 0.05/total number tests
                number_corrs_above_cutoff = sum(cluster_pvalues<(0.05/total_pairs))
            else:
                raise ValueError


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
    clusters_all_chr.append(get_clusters_chr(i, expression_df, residal_exp, total_pairs))

# write out
pd.concat(clusters_all_chr).to_csv(snakemake.output[0])