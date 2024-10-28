import numpy as np
import pandas as pd
import ast
from tqdm import tqdm 
from scipy.stats import linregress



from notebook_helper_functions import *
from annotate_clusters import *

def annotate_avg_expression_qtl(qtls, tissue_avg_expression):
    for idx, row in tqdm(qtls.iterrows(), total=len(qtls)):
        transcript_list = row['cluster_id'].split('_')
        cluster_avg_expression = tissue_avg_expression.loc[row['tissue_id'],transcript_list]
        qtls.loc[idx, 'avg_expression'] = np.mean(cluster_avg_expression)
        qtls.loc[idx, 'avg_log_expression'] = np.mean(np.log10(cluster_avg_expression))
        qtls.loc[idx, 'median_expression'] = np.median(cluster_avg_expression)
        qtls.loc[idx, 'min_expression'] = np.min(cluster_avg_expression)


def annotate_enhancers_qtl(qtls, gene_enhancer_df):
    for idx, qtl_row in tqdm(qtls.iterrows(), total=len(qtls)):
        matched_enhancers = gene_enhancer_df[(gene_enhancer_df['chr'] == qtl_row['chr'])&(gene_enhancer_df['enhancer_start'] <= qtl_row['position'])&(gene_enhancer_df['enhancer_end'] >= qtl_row['position'])]
        # is qtl in any abc enhancer?
        qtls.loc[idx, 'qtl_num_abc_enhancers'] = matched_enhancers['enhancer'].nunique()
        # how many genes does the enhancer contact?
        qtls.loc[idx, 'qtl_num_abc_genes'] = matched_enhancers['TargetGene'].nunique()
        # how many cluster genes?
        qtls.loc[idx, 'qtl_matched_abc_genes'] = pd.Series(matched_enhancers.reset_index()['transcript_id'].unique()).isin(qtl_row['cluster_id'].split('_')).sum()


def annotate_ctcf_tad_qtl(qtls, ctcf_df, tad_df, gid_gencode):
    # temp tss_min/tss_max for the ctcf annotation
    qtls['tss_min'] = qtls['position'].astype(float)
    qtls['tss_max'] = qtls['position'].astype(float) + 1
    qtls['Chromosome'] = qtls['chr']
    annotate_ctcf(qtls, ctcf_df)
    qtls['qtl_in_ctcf'] = qtls['has_ctcf_peak']
    # redo tss min and tss max
    annotate_positions(qtls, gid_gencode)
    # add tad annotations
    qtls['qtl_inter'] = pd.arrays.IntervalArray.from_arrays(qtls['position'].astype(float), qtls['position'].astype(float)+5000) # 10kb tad resolution
    qtls['num_tads_qtl'] = qtls.apply(count_tad_overlap, axis=1, args=(tad_df, 'qtl_inter'))
    qtls['qtl_in_tad'] = qtls['num_tads_qtl'] > 1
    qtls['between_tss'] = ((qtls['tss_min'] < qtls['position']) & (qtls['tss_max'] > qtls['position']))
    qtls['qtl_in_tss_ctcf'] = qtls['between_tss'] & qtls['qtl_in_ctcf']
    qtls['qtl_in_tad_ctcf'] = qtls['qtl_in_tad'] & qtls['qtl_in_ctcf']


def annotate_bidirectional_qtl(qtls, gid_gencode):
    qtls['in_bidirectional_promoter'] = False
    qtls['in_shared_promoter'] = False
    for qtl_idx, row in tqdm(qtls.iterrows(), total=len(qtls)):
        transcript_ids = row['cluster_id'].split('_')
        cluster_gencode = gid_gencode.loc[transcript_ids]
        # in the promoter if within 1kb of tss_start
        for idx, first_gene_row in cluster_gencode.iterrows():
            for idx, second_gene_row in cluster_gencode.iterrows():
                opp_strand = first_gene_row['strand'] != second_gene_row['strand']
                close = abs(first_gene_row['tss_start'] - second_gene_row['tss_start']) <= 1000
                if opp_strand & close:
                    # found a bidirectional promotor
                    if ((row['position'] - first_gene_row['tss_start']) < 1000) & ((row['position'] - second_gene_row['tss_start']) < 1000):
                        qtls.loc[qtl_idx, 'in_bidirectional_promoter'] = True
                elif close:
                    # found a shared promotor
                    if ((row['position'] - first_gene_row['tss_start']) < 1000) & ((row['position'] - second_gene_row['tss_start']) < 1000):
                        qtls.loc[qtl_idx, 'in_shared_promoter'] = True

# TODO add effect size annotation
# function to add all annotations, give correctly loaded data
def add_annotations_qtl(qtls, gid_gencode, gene_enhancer_df, ctcf_df, tad_df, tissue_avg_expression):
    qtls.reset_index(drop=True, inplace=True)
    print('started annotating qtls')
    annotate_avg_expression_qtl(qtls, tissue_avg_expression)
    print('annotated avg expression')
    annotate_enhancers_qtl(qtls, gene_enhancer_df)
    print('annotated enhancers')
    annotate_ctcf_tad_qtl(qtls, ctcf_df, tad_df, gid_gencode)
    print('annotated tads and ctcf')
    annotate_bidirectional_qtl(qtls, gid_gencode)
    print('annotated bidirectional promoters')
    annotate_distance(qtls, gid_gencode)


def load_and_annotate(qtls, my_tissue_id,
                      gencode_path='data/references/processed_gencode.v26.GRCh38.genes.csv', 
                      full_abc_path = 'data/references/functional_annotations/ABC_predictions/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz', 
                      abc_match_path='data/references/functional_annotations/ABC_predictions/ABC_matched_gtex.csv', 
                      ctcf_match_path='data/references/functional_annotations/ctcf_chip/ctcf_matched_gtex.csv', 
                      ctcf_dir='data/references/functional_annotations/ctcf_chip', 
                      tad_path='data/references/TAD_annotations/TADs_hg38/converted_HiC_IMR90_DI_10kb.txt',
                      avg_expression_path='data/processed/GTEx_Analysis_RSEMv1.gene_tpm.tissue_avg.csv',
                      verbosity=0):
    if verbosity:
        print('loading data')
    avg_expression = load_avg_exression(avg_expression_path)
    gid_gencode, full_gencode = load_gencode(f'{prefix}/{gencode_path}')
    gene_enhancer_df = load_abc(my_tissue_id, full_gencode, f'{prefix}/{full_abc_path}', f'{prefix}/{abc_match_path}')
    ctcf_df = load_ctcf(my_tissue_id, f'{prefix}/{ctcf_match_path}', f'{prefix}/{ctcf_dir}')
    tad_df = load_tad(f'{prefix}/{tad_path}')
    if verbosity:
        print('data loaded')
    qtls.reset_index(drop=True, inplace=True)
    add_annotations_qtl(qtls, gid_gencode, gene_enhancer_df, ctcf_df, tad_df, avg_expression)
    return qtls


# generate a pc annotation with the slope, r_value, p_value for each pc-egene pair
def get_annotate_pcs(pc_df, expression_df):
    sample_ids = pc_df.columns[pc_df.columns.str.contains('GTEX')]
    # add the slopes to the pc df
    for idx, row in tqdm(pc_df.iterrows(), total=pc_df.shape[0]):
        expression_cluster = expression_df[expression_df['cluster_id'] == row['cluster_id']].reset_index()
        gene_slopes = []
        gene_variances = []
        for egene_id in row['cluster_id'].split('_'):
            expression_values = expression_cluster[expression_cluster['egene_id']==egene_id][sample_ids].values.astype('float')            # pull the right data
            pc_values = row[sample_ids].astype('float')
            # get the r squared value
            slope, intercept, r_value, p_value, std_err = linregress(pc_values, expression_values)
            gene_slopes.append(slope)
            gene_variances.append(r_value**2)
        pc_df.loc[idx, 'r2_list'] = str(gene_variances)
        pc_df.loc[idx, 'slope_list'] = str(gene_slopes)
    annot_pc = pc_df.drop(columns=sample_ids)
    annot_pc.rename(columns={'gene_id':'pc_phenotype_id'}, inplace=True)
    annot_pc['egene_id'] = annot_pc['cluster_id'].str.split('_')
    annot_pc['egene_r2'] = annot_pc['r2_list'].apply(ast.literal_eval)
    annot_pc['egene_slope'] = annot_pc['slope_list'].apply(ast.literal_eval)
    return annot_pc.explode(['egene_id', 'egene_r2', 'egene_slope']).drop(columns=['r2_list', 'slope_list'])


def annotate_distance(qtls, gid_gencode):
    qtls['position'] = qtls['variant_id'].str.split('_').str[1].astype(int)
    qtls['cluster_min_distance'] = qtls.progress_apply(get_tss, axis=1, args=(gid_gencode,))


# distance to whichever gene in the cluster is closest
def get_tss(row, gid_gencode):
    cluster_gene_df = gid_gencode.loc[row['cluster_id'].split('_')]
    starts = cluster_gene_df['tss_start'].values
    distances = row['position'] - starts
    # return smallest absolute value distance
    idx = np.argmin(abs(distances))
    # make relative to gene orientation
    if cluster_gene_df.iloc[idx]['strand'] == '-':
        return -distances[idx]
    else:
        return distances[idx]