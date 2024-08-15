import sys
sys.path.append('/home/klawren/oak/pcqtls/workflow/scripts')
from annotate_clusters import *

import numpy as np
import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sp



def get_null_clusters(expressed_gencode, cluster_size, cluster_df=None):
    # sort
    expressed_gencode = expressed_gencode.sort_values(['chr', 'start', 'end'])
    # on a per chrom basis
    null_cluster_dfs = []
    for chr_id in range (1,23,1):
        chr_subset_gencode = expressed_gencode[expressed_gencode['chr'] == f'chr{chr_id}']
        transcripts = chr_subset_gencode['transcript_id'].astype(str) + ',' + chr_subset_gencode['transcript_id'].shift(-(cluster_size-1)).astype(str)
        chr_sizes = chr_subset_gencode['end'].shift(-(cluster_size-1)) - chr_subset_gencode['start']

        # trim off the blanks created from shifting
        transcripts = transcripts.iloc[:-(cluster_size-1)]
        chr_sizes = chr_sizes.iloc[:-(cluster_size-1)]

        # select those that are not already clusters
        try:
            cluster_start_ids = get_cluster_start_ids(cluster_df)
            this_cluster_size_start_ids = np.concatenate(cluster_start_ids[cluster_size:])
            in_cluster_bool = pd.Series(chr_subset_gencode['transcript_id'].iloc[:-(cluster_size-1)]).isin(this_cluster_size_start_ids).values
            null_cluster_dfs.append(pd.DataFrame({'Transcripts':transcripts[~in_cluster_bool], 'cluster_size':chr_sizes[~in_cluster_bool], 'chr':chr_id}))
        except TypeError:
            # no subselection wanted, i.e.None passed for cluster_df
            null_cluster_dfs.append(pd.DataFrame({'Transcripts':transcripts, 'cluster_size':chr_sizes, 'chr':chr_id}))

    null_df = pd.concat(null_cluster_dfs)
    null_df.reset_index(drop=True, inplace=True)
    null_df['Chromosome'] = null_df['chr']
    return null_df


# code from ben
# target has to be smaller or the two, n is length of target

def resample_dist(target, candidate_pool, n, seed=126124):   
    """Match a target distribution via weighted sampling from a candidate pool
    Args:
        target, candidate_pool: 1D numpy arrays with values ranging from 0 to 1
        n: integer number of indices to return
    Return:
        n indices to elements candidate_pool to use for the sample
    """
    rng = np.random.default_rng(seed)
    target_prob = sp.stats.gaussian_kde(target)
    candidate_prob = sp.stats.gaussian_kde(candidate_pool)

    bins = np.arange(0, 1, 0.0001)
    sampling_weight = target_prob(bins) / candidate_prob(bins)
    pool_bins = np.searchsorted(bins, candidate_pool) - 1
    pool_probability = sampling_weight[pool_bins]/sampling_weight[pool_bins].sum()

    return rng.choice(candidate_pool.size, size=n, replace=True, p=pool_probability)


def get_resamp_null_cluster(null_df, cluster_df, plot=False, number_null = 5000):
    # note that cluster_df should be resticted only to clusters with a matching number of transcripts already

    join_df = pd.concat([cluster_df, null_df], keys=['cluster', 'null'], names=['is_cluster', 'id'])
    if plot:
        # size distribution before
        sns.kdeplot(join_df.reset_index(), x='cluster_tss_size', hue='is_cluster', bw_adjust=.3, common_norm=False, log_scale=(10,None))
        plt.title('Distance distribution before resampling')
        plt.show()

    # add a normalized cluster size column to resample on
    cluster_df['normed_log_size'] = np.log10(cluster_df['cluster_tss_size'])/np.log10(join_df['cluster_tss_size'].max())
    null_df['normed_log_size'] = np.log10(null_df['cluster_tss_size'])/np.log10(join_df['cluster_tss_size'].max())

    # resample to match distance
    resamp_idxs = resample_dist(cluster_df['normed_log_size'], null_df['normed_log_size'], n=number_null)
    resamp_null_df = null_df.reset_index().iloc[resamp_idxs]

    if plot:
        # size distribution after resampling
        join_df = pd.concat([cluster_df, resamp_null_df], keys=['cluster', 'null'], names=['is_cluster', 'id'])
        sns.kdeplot(join_df, x='cluster_tss_size', hue='is_cluster', bw_adjust=.3, common_norm=False, log_scale=(10,None))
        plt.title('Distance distribution after resampling')
        plt.show()

    return resamp_null_df


def main():
    # Parse arguments from cmd
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tissue', help = 'which tissue we are in')
    parser.add_argument('-c', '--cluster_path', help = 'path to .csv clusters')
    parser.add_argument('-e', '--expression_path', help = 'path to .bed normalized expression')
    parser.add_argument('-co', '--covariates_path', help = 'path to covariates')
    parser.add_argument('-o', '--out_path', help='path to write out annotated null clusters')
    parser.add_argument('-g', '--gencode_path', help='path to gencode')
    parser.add_argument('--cluster_size', type=int, default=2, help = 'genes per cluster')
    parser.add_argument('--exclude_cluster_genes', type=int, default=1, help = 'exlude cluster genes? 0 for no, 1 for yes')
    parser.add_argument('--distance_matched', type=int, default=0, help = 'match to distances in cluster? 0 for no, 1 for yes')
    parser.add_argument('--verbosity', type=int, default=1, help = 'output verbosity')

    args = parser.parse_args()

    # load the clusters
    cluster_df = pd.read_csv(args.cluster_path, index_col=0)

    # load expressed gencode (to make null clusters out of)
    # this assumes only protien coding
    gid_gencode, full_gencode = load_gencode(args.gencode_path)
    expression_df = pd.read_csv(args.expression_path, sep='\t')
    expressed_gencode = full_gencode[full_gencode['transcript_id'].isin(expression_df['gene_id'])]
    expressed_gencode = expressed_gencode.sort_values(['chr', 'start', 'end'])

    # generate the null
    if args.exclude_cluster_genes ==1:
        # excluding clusters
        null_clusters = get_null_clusters(expressed_gencode, args.cluster_size, cluster_df=cluster_df).sample(5000)
    else:
        # including clusters
        null_clusters = get_null_clusters(expressed_gencode, args.cluster_size, cluster_df=None).sample(5000)

    if args.distance_matched:
        # distance matched null
        null_clusters = get_resamp_null_cluster(null_clusters, cluster_df[cluster_df['N_genes']==args.cluster_size], plot=False)

    # annotate the null
    load_and_annotate(null_clusters, args.tissue, args.covariates_path, args.expression_path, verbosity=args.verbosity)   
    # write out
    null_clusters.to_csv(args.out_path)

if __name__ == "__main__":
    main()
