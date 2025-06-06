import sys
sys.path.append('/home/klawren/oak/pcqtls/workflow/scripts')
from annotate_clusters import *

import numpy as np
import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sp
from get_pcs import get_pc_bed
from annotate_qtls import get_annotate_pcs


def get_null_clusters(expressed_gencode, cluster_size, cluster_df=None):
    # sort
    expressed_gencode = expressed_gencode.sort_values(['chr', 'start', 'end'])
    # on a per chrom basis
    null_cluster_dfs = []
    for chr_id in range (1,23,1):
        chr_subset_gencode = expressed_gencode[expressed_gencode['chr'] == f'chr{chr_id}']
        transcripts = chr_subset_gencode['transcript_id'].astype(str)
        for i in range(1,cluster_size):
            transcripts = transcripts + ',' + chr_subset_gencode['transcript_id'].shift(-i).astype(str)
        chr_sizes = chr_subset_gencode['end'].shift(-(cluster_size-1)) - chr_subset_gencode['start']

        # trim off the blanks created from shifting
        transcripts = transcripts.iloc[:-(cluster_size-1)]
        chr_sizes = chr_sizes.iloc[:-(cluster_size-1)]

        chr_null_df = pd.DataFrame({'Transcripts':transcripts, 'cluster_size':chr_sizes, 'chr':chr_id})

        # select those that are not already clusters
        try:
            cluster_ids = np.concatenate(cluster_df['Transcripts'].str.split(',').values)
            chr_null_df['transcript_list'] = chr_null_df['Transcripts'].str.split(',')
            # exclude any null clusters that contain two genes both in a cluster
            null_clusters_exploded = chr_null_df.explode('transcript_list')
            null_clusters_exploded['in_cluster'] = null_clusters_exploded['transcript_list'].isin(cluster_ids).astype(int)
            null_clusters_cluster_df_count = null_clusters_exploded.groupby('Transcripts').agg({'in_cluster':sum})
            exlcuded_clusters = null_clusters_cluster_df_count[null_clusters_cluster_df_count['in_cluster'] > 1].index
            # add to the nulls
            null_cluster_dfs.append(chr_null_df[~chr_null_df['Transcripts'].isin(exlcuded_clusters)])
        except TypeError:
            # no subselection wanted, i.e.None passed for cluster_df
            null_cluster_dfs.append(chr_null_df)

    null_df = pd.concat(null_cluster_dfs)
    null_df.reset_index(drop=True, inplace=True)
    null_df['Chromosome'] = null_df['chr']
    null_df['N_genes'] = cluster_size
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


# takes about 5 minutes

def get_null_pcs_annotated(config, cluster_size, tissue_id):
    gid_gencode, full_gencode = load_gencode()
    expression_path = "{}/{}/{}.v8.normalized_expression.bed".format(prefix, config["expression_dir"], tissue_id)
    expression_df = pd.read_csv(expression_path, sep='\t')
    expressed_gencode = full_gencode[full_gencode['transcript_id'].isin(expression_df['gene_id'])]
    cluster_df = pd.read_csv('{}/{}/{}_clusters_all_chr.csv'.format(prefix, config['clusters_dir'], tissue_id),index_col=0)
    null_clusters = get_null_clusters(expressed_gencode, cluster_size, cluster_df=cluster_df)
    covariates_path = "{}/{}/{}.v8.covariates.txt".format(prefix, config["covariates_dir"], tissue_id)
    covariates_df = pd.read_csv(covariates_path, sep='\t', index_col=0).T
    residal_exp = calculate_residual(expression_df[covariates_df.index], covariates_df, center=True)
    residal_exp = pd.DataFrame(residal_exp, columns=covariates_df.index, index=expression_df['gene_id'])

    # residulize the expression 
    expression_df_gid = expression_df.set_index('gene_id')
    expression_df_res = expression_df_gid.copy()
    expression_df_res[expression_df.columns[4:]] = residal_exp

    # get cluster based expression (from snakemake_filter_expression_clusters)
    cluster_sets = []
    for idx, row in tqdm(null_clusters.iterrows(), total=len(null_clusters)):
        cluster_set = expression_df_res.loc[row['Transcripts'].split(',')].copy()
        cluster_set['start'] = cluster_set['start'].min()
        cluster_set['end'] = cluster_set['end'].max()
        cluster_set.insert(3, 'gene_id', ['_'.join([*sorted(cluster_set.index.values), 'e', gid]) for gid in cluster_set.index])
        cluster_sets.append(cluster_set)

    cluster_expresison = pd.concat(cluster_sets)
    # get the null pcs
    null_pcs = get_pc_bed(null_clusters, cluster_expresison, covariates_df)
    null_pcs['cluster_id'] = null_pcs['gene_id'].str.split('_pc').str[0]
    null_pcs['pc_id'] = null_pcs['gene_id'].str.split('_pc').str[1].astype('float')
    null_pcs['cluster_size'] = null_pcs['cluster_id'].str.split('_').apply(len)
    annotated_null_pcs = get_annotate_pcs(null_pcs.reset_index(drop=True), cluster_expresison.drop(columns='gene_id'))
    return annotated_null_pcs


def main():
    # Parse arguments from cmd
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tissue', help = 'which tissue we are in')
    parser.add_argument('-c', '--cluster_path', help = 'path to .csv clusters')
    parser.add_argument('-e', '--expression_path', help = 'path to .bed normalized expression')
    parser.add_argument('-co', '--covariates_path', help = 'path to covariates')
    parser.add_argument('-o', '--out_path', help='path to write out annotated null clusters')
    parser.add_argument('--cluster_size', type=int, default=2, help = 'genes per cluster')
    parser.add_argument('--exclude_cluster_genes', type=int, default=1, help = 'exlude cluster genes? 0 for no, 1 for yes')
    parser.add_argument('--distance_matched', type=int, default=0, help = 'match to distances in cluster? 0 for no, 1 for yes')
    parser.add_argument('--gencode_path', help = 'path to gencode')
    parser.add_argument('--full_abc_path', help = 'path to abc')
    parser.add_argument('--abc_match_path', help = 'path to abc gtex matching')
    parser.add_argument('--ctcf_match_path', help = 'path to ctcf matching')
    parser.add_argument('--ctcf_dir', help = 'path to ctcf dir')
    parser.add_argument('--paralog_path', help = 'path to paralogs')
    parser.add_argument('--go_path', help = 'path to go')
    parser.add_argument('--cross_map_path', help = 'path to cross map')
    parser.add_argument('--tad_path', help = 'path to tad data')
    parser.add_argument('--verbosity', type=int, default=0, help = 'output verbosity')

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
    load_and_annotate(null_clusters, args.tissue, args.covariates_path, args.expression_path, 
                      args.gencode_path, args.full_abc_path, args.abc_match_path, args.ctcf_match_path,
                      args.ctcf_dir, args.paralog_path, args.go_path, args.cross_map_path, args.tad_path, args.verbosity)   
    # write out
    null_clusters.to_csv(args.out_path)

if __name__ == "__main__":
    main()
