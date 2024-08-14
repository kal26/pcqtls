import numpy as np
import pandas as pd
import argparse
from tqdm import tqdm 
import yaml


import sys
sys.path.append('/home/klawren/oak/pcqtls/workflow/scripts')
from residualize import calculate_residual

#working directory prefix
prefix = '/home/klawren/oak/pcqtls'


# functions to load in external data

def load_gencode(gencode_path, protien_coding_only = True):
    # load in gene data
    full_gencode=pd.read_csv(gencode_path)
    # filter to protien coding
    if protien_coding_only:
        non_protein_gencode = full_gencode.copy()
        full_gencode = full_gencode[full_gencode['gene_type'] == 'protein_coding']
    gid_gencode = full_gencode.set_index('transcript_id').drop_duplicates()
    return gid_gencode, full_gencode

def get_redidual_expression(covariates_path, expression_path):
    # residulized expression for correlations
    covariates_df = pd.read_csv(covariates_path, sep='\t', index_col=0).T
    expression_df = pd.read_csv(expression_path, sep='\t')
    gid_expression = expression_df.set_index('gene_id')[covariates_df.index]
    residal_exp = calculate_residual(gid_expression[covariates_df.index], covariates_df, center=True)
    residal_exp = pd.DataFrame(residal_exp, columns=covariates_df.index, index=gid_expression.index)
    return residal_exp


def load_abc(full_gencode, my_tissue_id, full_abc_path, abc_match_path):
    # load in ABC data for enhancer gene connections
    full_abc_pred_df = pd.read_csv(full_abc_path, sep='\t')
    # load in tissue matching for ABC-gtex tissues
    # this has to be done by hand
    # some were clear, for those that aren't I've asked for help and put None for now
    abc_gtex_match = pd.read_csv(abc_match_path)
    # get just the enhancer-gene connections for the matched tissue
    abc_df = full_abc_pred_df[full_abc_pred_df['CellType'] == abc_gtex_match[abc_gtex_match['GTEX_tissue'] == my_tissue_id]['ABC_biosample_id'].iloc[0]]
    # add transcript ids to relevant abc enhancer-gene connection columns these and set as index
    gene_enhancer_df = pd.merge(full_gencode[['transcript_id', 'gene_name']], abc_df[['TargetGene','name','class', 'ABC.Score']], left_on='gene_name', right_on='TargetGene', how='left')
    gene_enhancer_df.rename(columns={'name':'enhancer'}, inplace=True)
    gene_enhancer_df.set_index('transcript_id', inplace=True)
    gene_enhancer_df.dropna(inplace=True)
    return gene_enhancer_df


def load_ctcf(my_tissue_id, ctcf_match_path, ctcf_dir):
    # load in ctcf data
    ctcf_gtex_match = pd.read_csv(ctcf_match_path)
    ctcf_file = ctcf_gtex_match[ctcf_gtex_match['GTEX'] == my_tissue_id].iloc[0]['ctcf']
    ctcf_df = pd.read_csv(f'{ctcf_dir}/{ctcf_file}.bed.gz', sep='\t', 
                        names=['chr', 'start', 'end', 'name', 'score', 'strand', 'signal_value', 'p_value', 'q_value', 'peak'])
    ctcf_df['peak_inter'] = pd.arrays.IntervalArray.from_arrays(ctcf_df['start'], ctcf_df['end'])
    ctcf_df['point_inter'] = pd.arrays.IntervalArray.from_arrays(ctcf_df['start'] + ctcf_df['peak'], ctcf_df['start'] + ctcf_df['peak'] + 1)
    ctcf_df['Chromosome'] = ctcf_df['chr'].str.split('chr').str[1].astype(str)
    return ctcf_df


def load_paralogs(paralog_path):
    # load in paralogs
    paralog_df = pd.read_csv(paralog_path, sep='\t')
    # drop genes that don't have paralogues
    paralog_df.dropna(subset=['Human paralogue gene stable ID'], inplace=True)
    # group by the gene that has the paralogs (this is bidirectional)
    # this needs to be done on gene ids without versions, as this isn't the same version number (sadly the correct version isn't availble on biomart)
    paralog_df = paralog_df.groupby('Gene stable ID')['Human paralogue gene stable ID'].apply(set)
    return paralog_df


def load_go(go_path):
    # load in go terms
    go_df = pd.read_csv(go_path, sep='\t', header=None,
                        names = ['Gene stable ID', 'Gene stable ID version', 'Gene start (bp)', 'Gene end', 'Strand', 
                                 'tss', 'gencode_annotation', 'gene_name', 'transcript_type', 'go_accession', 'go_name', 'go_domain'])
    # only consider matching biological process go terms, as in Ribiero 2021
    go_df = go_df[go_df['go_domain'] == 'biological_process'].groupby('Gene stable ID')['go_accession'].apply(set)
    return go_df


def load_cross_map(cross_map_path):
    # load in cross mapablity (this is cleaned up in cross_mappability.ipynb)
    cross_mappability = pd.read_csv(cross_map_path)
    cross_mappability.set_index('gene_1', inplace=True)
    return cross_mappability



# functions to annotate a cluster df

def get_cluster_size(row, gid_gencode):
    transcript_ids = row['Transcripts'].split(',')
    cluster_gencode = gid_gencode.loc[transcript_ids]
    return  cluster_gencode['end'].max() - cluster_gencode['start'].min()

def get_cluster_tss_size(row, gid_gencode):
    transcript_ids = row['Transcripts'].split(',')
    cluster_gencode = gid_gencode.loc[transcript_ids]
    return  cluster_gencode['tss_start'].max() - cluster_gencode['tss_start'].min()

def get_cluster_start_ids(cluster_df):
    # the first for pairs, the first and second for threes, ect
    cluster_start_ids = []
    for i in range(cluster_df['N_genes'].max()):
        out_ids = cluster_df[cluster_df['N_genes'] == i]['Transcripts'].str.split(',').str[:i-1].values
        if len(out_ids)>0:
            cluster_start_ids.append(np.concatenate(out_ids))
        else:
            cluster_start_ids.append([])
    return cluster_start_ids

def annotate_sizes(cluster_df, gid_gencode):
    cluster_df['cluster_size'] = cluster_df.apply(get_cluster_size, axis=1, args=(gid_gencode,))
    cluster_df['cluster_tss_size'] = cluster_df.apply(get_cluster_tss_size, axis=1, args=(gid_gencode,))

def annotate_positions(cluster_df, gid_gencode):
    for idx, row in cluster_df.iterrows():
        transcript_ids = row['Transcripts'].split(',')
        cluster_gencode = gid_gencode.loc[transcript_ids]
        cluster_df.loc[idx, 'start'] = cluster_gencode['start'].min()
        cluster_df.loc[idx, 'end'] = cluster_gencode['end'].max()

def get_bidirectional(row, gid_gencode):
    transcript_ids = row['Transcripts'].split(',')
    cluster_gencode = gid_gencode.loc[transcript_ids]
    num_bidirectional = 0
    # check all pairwise combos of genes
    for idx, first_gene_row in cluster_gencode.iterrows():
        for idx, second_gene_row in cluster_gencode.iterrows():
            opp_strand = first_gene_row['strand'] != second_gene_row['strand']
            close = abs(first_gene_row['tss_start'] - second_gene_row['tss_start']) <= 1000
            if opp_strand & close:
                # found a bidirectional promotor
                num_bidirectional +=1

    # didn't find a bidirectional promotor
    return num_bidirectional/2

def annotate_bidirectional(cluster_df, gid_gencode):
    cluster_df['num_bidirectional_promoter'] = cluster_df.apply(get_bidirectional, axis=1, args=(gid_gencode,))
    cluster_df['has_bidirectional_promoter'] = cluster_df['num_bidirectional_promoter'] > 0

# annotate clusters with the number of shared enhancers

def annotate_enhancers(cluster_df, gene_enhancer_df):
    for idx, row in tqdm(cluster_df.iterrows(), total=len(cluster_df)):
        enhancer_list = gene_enhancer_df[gene_enhancer_df.index.isin(row['Transcripts'].split(','))]
        full_enhancer_list = enhancer_list['enhancer']
        strong_enhancer_list = enhancer_list[enhancer_list['ABC.Score']>=0.1]['enhancer']
        very_strong_enhancer_list = enhancer_list[enhancer_list['ABC.Score']>=0.25]['enhancer']
        num_shared_enhancers = sum(full_enhancer_list.duplicated())
        num_shared_strong_enhancers = sum(strong_enhancer_list.duplicated())
        num_shared_very_strong_enhancers = sum(very_strong_enhancer_list.duplicated())
        cluster_df.loc[idx, 'num_shared_enhancers'] = num_shared_enhancers
        cluster_df.loc[idx, 'num_shared_strong_enhancers'] = num_shared_strong_enhancers
        cluster_df.loc[idx, 'num_enhancers'] = len(full_enhancer_list)
        cluster_df.loc[idx, 'num_strong_enhancers'] = len(strong_enhancer_list)
        cluster_df.loc[idx, 'has_shared_enhancer'] = num_shared_enhancers > 0
        cluster_df.loc[idx, 'has_shared_strong_enhancer'] = num_shared_strong_enhancers > 0
        cluster_df.loc[idx, 'has_shared_very_strong_enhancer'] = num_shared_very_strong_enhancers > 0

def annotate_ctcf(cluster_df, ctcf_df):
    cluster_df['interval'] = pd.arrays.IntervalArray.from_arrays(cluster_df['start'], cluster_df['end'])
    # ctcf intervals for each chromosome
    chr_ctcf_peaks={}
    chr_ctcf_points={}
    for chr in cluster_df['Chromosome'].unique():
        ctcf_chr = ctcf_df[ctcf_df['Chromosome'] == chr.astype(str)]
        chr_ctcf_peaks[chr] = pd.arrays.IntervalArray.from_arrays(ctcf_chr['start'], ctcf_chr['end'])
        chr_ctcf_points[chr] = pd.arrays.IntervalArray.from_arrays(ctcf_chr['start'] + ctcf_chr['peak'], ctcf_chr['start'] + ctcf_chr['peak'] + 1)
    # annotate each cluster
    for idx, row in tqdm(cluster_df.iterrows(), total=len(cluster_df)):
        num_ctcf_peak = sum(chr_ctcf_peaks[row['Chromosome']].overlaps(row['interval']))
        num_ctcf_point = sum(chr_ctcf_points[row['Chromosome']].overlaps(row['interval']))

        cluster_df.loc[idx, 'num_ctcf_peak'] = num_ctcf_peak
        cluster_df.loc[idx, 'has_ctcf_peak'] = num_ctcf_peak > 0
        cluster_df.loc[idx, 'num_ctcf_point'] = num_ctcf_point
        cluster_df.loc[idx, 'has_ctcf_point'] = num_ctcf_point > 0



def annotate_enhancers_jaccard(cluster_df, gene_enhancer_df):
    for idx, row in tqdm(cluster_df.iterrows(), total=len(cluster_df)):
        transcript_list = row['Transcripts'].split(',')
        jaccards_unweighted=[]
        jaccards_weighted=[]
        for i in range(len(transcript_list)):
            for j in range(i):
                enhancer_list = gene_enhancer_df[gene_enhancer_df.index.isin([transcript_list[i], transcript_list[j]])]
                enhancer_list['ABC.Score_min'] = enhancer_list['ABC.Score']
                enhancer_min_max = enhancer_list.groupby('enhancer').agg({'ABC.Score':'max', 'ABC.Score_min':'min'})

                # zero out the mins for those elements that exist for only 1 gene
                enhancer_gene_counts = enhancer_list.groupby('enhancer').agg({'gene_name':'nunique'}) 
                single_enhancers = enhancer_gene_counts.index.values[enhancer_gene_counts['gene_name'] == 1]
                enhancer_min_max.loc[single_enhancers, 'ABC.Score_min'] = 0
                # jaccard without reweighting
                jaccards_unweighted.append(enhancer_min_max['ABC.Score_min'].sum()/enhancer_min_max['ABC.Score'].sum())

                # assuming these don't sum to 1 for a given gene becuase the promoter-self connections aren't listed
                # add in an element for each genes promotor to get the final weighting right
                reweightings = enhancer_list.groupby('gene_name').agg({'ABC.Score':sum})
                reweightings['ABC.Score'] = 1- reweightings['ABC.Score']
                reweightings['ABC.Score_min'] = 0

                enhancer_min_max = pd.concat([enhancer_min_max, reweightings])
                # jaccard with reweighting
                jaccards_weighted.append(enhancer_min_max['ABC.Score_min'].sum()/enhancer_min_max['ABC.Score'].sum())

        #jaccards_unweighted = np.nan_to_num(jaccards_unweighted)
        #jaccards_weighted = np.nan_to_num(jaccards_weighted)

        cluster_df.loc[idx, 'max_jaccard_unweighted'] = max(jaccards_unweighted)
        cluster_df.loc[idx, 'max_jaccard_weighted'] = max(jaccards_weighted)
        cluster_df.loc[idx, 'has_high_jaccard_unweighted'] = max(jaccards_unweighted) > 0.5
        cluster_df.loc[idx, 'has_high_jaccard_weighted'] = max(jaccards_weighted) > 0.1
        cluster_df.loc[idx, 'mean_jaccard_unweighted'] = np.average(jaccards_unweighted)
        cluster_df.loc[idx, 'mean_jaccard_weighted'] = np.average(jaccards_weighted)


def annotate_correlation(cluster_df, residal_exp):
    for idx, row in tqdm(cluster_df.iterrows(), total=len(cluster_df)):
        transcript_list = row['Transcripts'].split(',')
        cluster_expression = residal_exp.loc[transcript_list].T.corr('spearman').to_numpy()
        cluster_corr = cluster_expression[np.triu_indices(len(cluster_expression), k=1)]
        cluster_df.loc[idx, 'Mean_cor'] = cluster_corr.mean()
        cluster_df.loc[idx, 'Mean_pos_cor'] = cluster_corr[cluster_corr>0].mean()
        cluster_df.loc[idx, 'Mean_neg_cor'] = cluster_corr[cluster_corr<0].mean()


def annotate_go(cluster_df, go_df):
    for idx, row in tqdm(cluster_df.iterrows(), total=len(cluster_df)):
        transcript_list_versions = row['Transcripts'].split(',')
        transcript_list_no_versions = [transcript.split('.')[0] for transcript in transcript_list_versions]

        go_list = go_df[go_df.index.isin(transcript_list_no_versions)]
        num_shared_go_all = sum(go_list.duplicated())
        # number genes that share all their go terms with another gene
        cluster_df.loc[idx, 'num_shared_go_all'] = num_shared_go_all
        cluster_df.loc[idx, 'has_shared_go_all'] = num_shared_go_all > 0

        # number go terms shared between any genes
        num_shared_go_any = sum(go_list.explode().duplicated())
        cluster_df.loc[idx, 'num_shared_go_any'] = num_shared_go_any
        cluster_df.loc[idx, 'has_shared_go_any'] = num_shared_go_any > 0


# annotate paralogs
def get_paralogs(row, paralog_df):
    transcript_list_versions = row['Transcripts'].split(',')
    transcript_list_no_versions = set([transcript.split('.')[0] for transcript in transcript_list_versions])

    paralogs = 0
    for transcript in transcript_list_no_versions:
        try:
            has_paralog = bool(paralog_df.loc[transcript] & transcript_list_no_versions)
            paralogs += has_paralog
        except KeyError:
            # if this isn't in the paralog df, it has no paralogs, so continue
            pass
    return paralogs

def annotate_paralogs(cluster_df, paralog_df):
    cluster_df['num_paralog'] = cluster_df.apply(get_paralogs, axis=1, args=(paralog_df,))
    cluster_df['has_paralog'] = cluster_df['num_paralog'] > 0


def get_cross_map(row, cross_mappability, cross_map_threshold=100):
    # number of transcripts that cross map to some other transcript in the cluster
    transcript_list = set(row['Transcripts'].split(','))
    cross_maps = 0
    for transcript in transcript_list:
        try:
            cross_map_this_transcript = cross_mappability.loc[transcript]
            pass_threshold_mask = np.asarray(ast.literal_eval(cross_map_this_transcript['cross_mappability'])) > cross_map_threshold
            cross_map_this_transcript = np.asarray(ast.literal_eval(cross_map_this_transcript['gene_2_full']))[pass_threshold_mask]
            has_cross_map = bool(set(cross_map_this_transcript) & transcript_list)
            cross_maps += has_cross_map
        except KeyError:
            # if this isn't in the paralog df, it has no paralogs, so continue
            pass
    return cross_maps

def annotate_cross_maps(cluster_df, cross_mappability):
    cluster_df['num_cross_map'] = cluster_df.apply(get_cross_map, axis=1, args=(cross_mappability,))
    cluster_df['has_cross_map'] = cluster_df['num_cross_map'] > 0



# function to add all annotations, give correctly loaded data
def add_annotations(cluster_df, gid_gencode, gene_enhancer_df, paralog_df, cross_mappability, go_df, ctcf_df, residal_exp):
    annotate_sizes(cluster_df, gid_gencode)
    print('annotated sizes')
    annotate_positions(cluster_df, gid_gencode) # this must go before annotate ctcf
    print('annotated positions')
    annotate_bidirectional(cluster_df, gid_gencode)
    print('annotated bidirectional promoters')
    annotate_enhancers(cluster_df, gene_enhancer_df)
    print('annotated enhancers')
    annotate_paralogs(cluster_df, paralog_df)
    print('annotated paralogs')
    annotate_cross_maps(cluster_df,cross_mappability)
    print('annotated cross mappability')
    annotate_go(cluster_df, go_df)
    print('annotated go')
    annotate_enhancers_jaccard(cluster_df, gene_enhancer_df)
    print('annotated enhancers with jaccard index')
    annotate_ctcf(cluster_df, ctcf_df)
    print('annotated ctcf')
    try:
        cluster_df['has_neg_corr'] = ~cluster_df['Mean_neg_cor'].isna()
        cluster_df['has_high_pos_corr'] = cluster_df['Mean_pos_cor'] > .5
    except KeyError:
        annotate_correlation(cluster_df, residal_exp)
        cluster_df['has_neg_corr'] = ~cluster_df['Mean_neg_cor'].isna()
        cluster_df['has_high_pos_corr'] = cluster_df['Mean_pos_cor'] > .5
        print('annotated correlations')

def load_and_annotate(cluster_df, my_tissue_id, covariates_path, expression_path,
                      gencode_path='data/references/processed_gencode.v26.GRCh38.genes.csv', 
                      full_abc_path = 'data/references/functional_annotations/ABC_predictions/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz', 
                      abc_match_path='data/references/functional_annotations/ABC_predictions/ABC_matched_gtex.csv', 
                      ctcf_match_path='data/references/functional_annotations/ctcf_chip/ctcf_matched_gtex.csv', 
                      ctcf_dir='data/references/functional_annotations/ctcf_chip', 
                      paralog_path='/data/references/functional_annotations/paralogs_biomart_ensembl97.tsv.gz', 
                      go_path='data/references/functional_annotations/go_biomart_ensembl97.tsv.gz', 
                      cross_map_path='data/references/cross_mappability/cross_mappability_100_agg.csv'):
    gid_gencode, full_gencode = load_gencode(f'{prefix}/{gencode_path}')
    gene_enhancer_df = load_abc(full_gencode, my_tissue_id, f'{prefix}/{full_abc_path}', f'{prefix}/{abc_match_path}')
    ctcf_df = load_ctcf(my_tissue_id, f'{prefix}/{ctcf_match_path}', f'{prefix}/{ctcf_dir}')
    paralog_df = load_paralogs(f'{prefix}/{paralog_path}')
    go_df = load_go(f'{prefix}/{go_path}')
    cross_mappability = load_cross_map(f'{prefix}/{cross_map_path}')
    residal_exp = get_redidual_expression(covariates_path, expression_path)
    add_annotations(cluster_df, gid_gencode, gene_enhancer_df, paralog_df, cross_mappability, go_df, ctcf_df, residal_exp)


def run_annotate_from_config(config_path, my_tissue_id):
    # general paths from config
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    clusters_dir = config['clusters_dir']
    expression_dir = config['expression_dir']
    expression_path = f'{prefix}/{expression_dir}/{my_tissue_id}.v8.normalized_expression.bed'
    covariates_dir = config['covariates_dir']
    covariates_path = f'{prefix}/{covariates_dir}/{my_tissue_id}.v8.covariates.txt'

    cluster_df = pd.read_csv(f'{prefix}/{clusters_dir}/{my_tissue_id}_clusters_all_chr.csv', index_col=0)

    load_and_annotate(cluster_df, my_tissue_id, covariates_path, expression_path)
    return cluster_df

def run_annotate_from_paths(my_tissue_id, clusters_path, expression_path, covariates_path):
    # this version for use in snakemake
    cluster_df = pd.read_csv(clusters_path, index_col=0)
    load_and_annotate(cluster_df, my_tissue_id, covariates_path, expression_path)
    return cluster_df


def main():
    # Parse arguments from cmd
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tissue', help = 'which tissue we are in')
    parser.add_argument('-c', '--cluster_path', help = 'path to .csv clusters')
    parser.add_argument('-e', '--expression_path', help = 'path to .bed normalized expression')
    parser.add_argument('-co', '--covariates_path', help = 'path to covariates')
    parser.add_argument('-o', '--out_path', help='path to write out annotated clusters')
    parser.add_argument('--verbosity', type=int, default=0, help = 'output verbosity')

    args = parser.parse_args()
    # call the pc funciton
    cluster_df_annotated = run_annotate_from_paths(args.tissue, args.cluster_path, args.expression_path, args.covariates_path)
    cluster_df_annotated.to_csv(args.out_path)

if __name__ == "__main__":
    main()

