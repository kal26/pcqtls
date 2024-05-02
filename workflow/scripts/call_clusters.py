print('running call clusters')
import pandas as pd
import numpy as np
from residualize import calculate_residual
from scipy.stats import spearmanr
import argparse



# get expression clusters on a per-chrom basis
def get_clusters_chr(chr_id, expression_df, residal_exp, total_pairs, tissue_id, min_cluster_size, max_cluster_size, min_corr_cutoff, percent_corr_cutoff, cutoff_type):
    # genes on this chr
    chr_gene_ids = expression_df[expression_df['#chr'] == f'chr{chr_id}']['gene_id']
    chr_residual_exp = residal_exp.loc[chr_gene_ids]

    # get correlation and pvalues
    try:
        chr_corr, chr_pvalue = spearmanr(chr_residual_exp, axis=1)
    except IndexError:
        # no genes on this chr
        return pd.DataFrame({})
    chr_corr = pd.DataFrame(chr_corr, index=chr_residual_exp.index, columns=chr_residual_exp.index)
    chr_pvalue = pd.DataFrame(chr_pvalue, index=chr_residual_exp.index, columns=chr_residual_exp.index)

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
                number_corrs_above_cutoff = sum(sum(abs(cluster_candidate.values)>min_corr_cutoff))
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

def get_clusters_from_paths(expression_path, covariates_path, tissue_id, min_cluster_size=2, max_cluster_size=50, min_corr_cutoff=0.1, percent_corr_cutoff=.7, cutoff_type='pvalue'):
    # load data
    expression_df = pd.read_csv(expression_path, sep='\t')
    covariates_df = pd.read_csv(covariates_path, sep='\t', index_col=0).T
    print('loaded data')


    # residulize the expression 
    residal_exp = calculate_residual(expression_df[covariates_df.index], covariates_df, center=True)
    residal_exp = pd.DataFrame(residal_exp, columns=covariates_df.index, index=expression_df['gene_id'])
    print('residualized expression')

    return get_clusters(expression_df, residal_exp, tissue_id, min_cluster_size=min_cluster_size, max_cluster_size=max_cluster_size, min_corr_cutoff=min_corr_cutoff, percent_corr_cutoff=percent_corr_cutoff, cutoff_type=cutoff_type)


def get_clusters(expression_df, residal_exp, tissue_id, min_cluster_size=2, max_cluster_size=50, min_corr_cutoff=0.1, percent_corr_cutoff=.7, cutoff_type='pvalue'):
    # calculate total number of pairs considered for bonferroni correction
    total_pairs = 0
    for i in np.arange(1,23,1):
        chr_gene_ids = expression_df[expression_df['#chr'] == f'chr{i}']['gene_id']
        upper_corner_idxs = np.triu(np.ones(len(chr_gene_ids)), k=1)
        excluded_cluster_size_idxs = np.triu(np.ones(len(chr_gene_ids)), k=max_cluster_size)
        total_pairs += upper_corner_idxs.sum()  - excluded_cluster_size_idxs.sum()

    # cycle thorugh all chrs
    clusters_all_chr = []
    for i in np.arange(1,23,1):
        print(f'Working on chr{i}')
        clusters_all_chr.append(get_clusters_chr(i, expression_df, residal_exp, total_pairs, tissue_id, min_cluster_size, max_cluster_size, min_corr_cutoff, percent_corr_cutoff, cutoff_type))

    # concat and return 
    return pd.concat(clusters_all_chr)


def main():
    # Parse arguments from cmd
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--expression_path', help = 'path to .bed normalized expression (not residulaized)')
    parser.add_argument('-co', '--covariates_path', help = 'path to covariates')
    parser.add_argument('-o', '--out_path', help = 'path where cluster results should be written out')
    parser.add_argument('--verbosity', type=int, default=0, help = 'output verbosity')
        # set params

    parser.add_argument('--max_cluster_size', type=int, default=50, help = 'maximum number of genes per cluster')
    parser.add_argument('--min_cluster_size', type=int, default=2, help = 'minimum number of genes per cluster')
    parser.add_argument('--min_corr_cutoff', type=float, default=.01, help = 'minimum corr value to consider (if using value based cutoff)')
    parser.add_argument('--percent_corr_cutoff', type=float, default=.7, help = 'number of corrs above threshold required')
    parser.add_argument('--cutoff_type', type=str, default='pvalue', help = 'value for value based, pvalue for p value based (bonferroni corrected)')
    parser.add_argument('--tissue_id', type=str, default='None', help = 'id for the tissue')

    args = parser.parse_args()

    # call the clusters funciton
    clusters_all_chr = get_clusters(args.expression_path, args.covariates_path, args.tissue_id, args.min_cluster_size, args.max_cluster_size, args.min_corr_cutoff, args.percent_corr_cutoff, args.cutoff_type)

    # write out 
    clusters_all_chr.to_csv(args.out_path)


if __name__ == "__main__":
    main()