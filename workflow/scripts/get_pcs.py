import numpy as np
import pandas as pd
import argparse
from sklearn.decomposition import PCA
import os
import tensorqtl
import torch
from tensorqtl.core import Residualizer
from residualize import calculate_residual


def make_bed_order(df):
    # the column order matters, so rearrange columns
    cols = list(df)
    cols.insert(0, cols.pop(cols.index('#chr')))
    cols.insert(1, cols.pop(cols.index('start')))
    cols.insert(2, cols.pop(cols.index('end')))
    cols.insert(3, cols.pop(cols.index('gene_id')))
    df = df.loc[:, cols]
    return df

def pc_bed(cluster_path, expression_path, covariates_path, pc_out_path, verb=0):
   
    # load in data
    if verb:
        print('loading data')
    cluster_df = pd.read_csv(cluster_path)
    # expression is already residualized
    expression_df = pd.read_csv(expression_path, sep='\t')
    covariates_df = pd.read_csv(covariates_path, sep='\t', index_col=0).T

    # pull out sample ids
    sample_ids = expression_df.columns[4:]


    # add .bed info to cluster
    expression_df['egene_id'] = expression_df['gene_id'].str.split('_e_').str[1]
    expression_df['cluster_id'] = expression_df['gene_id'].str.split('_e_').str[0]
    expression_df_gid = expression_df.set_index('egene_id')

    # iterate through clusters and gather PCs
    cluster_pcs_dfs = []
    for idx, row in cluster_df.iterrows():

        # cluster id
        cluster_id = '_'.join([*sorted(row['Transcripts'].split(','))])

        # get the gene expresion for genes in this cluster
        cluster_expression_df_gid = expression_df_gid[expression_df_gid['cluster_id']==cluster_id]
        cluster = cluster_expression_df_gid.loc[row['Transcripts'].split(',')]

        # get all pcs 
        X = cluster[sample_ids].transpose()
        pca = PCA()
        pc_values = pca.fit_transform(X)

        # only take those with eigenvalue > .1? 
        # this leaves nothing - lower threshold? or do this post hoc?
        # as in wang et al 2022
        #pc_values = pc_values[:,pca.explained_variance_ > .1]

        # get an id for each pc
        gene_ids = []
        for pc_num in range(pc_values.shape[1]):
            gene_ids.append('_'.join([*sorted(row['Transcripts'].split(',')), f'pc{pc_num+1}']))

        # center normalized and residualize the pcs
        normed_residualized_pcs = calculate_residual(pd.DataFrame(pc_values.T, columns = sample_ids), covariates_df, center=True)

        # make a dataframe of the pcs for this cluster
        cluster_pcs_df = pd.DataFrame(normed_residualized_pcs, 
                    columns = sample_ids, 
                    index = gene_ids)


        # make the gene id index a column
        cluster_pcs_df = cluster_pcs_df.reset_index().rename(columns={'index': 'gene_id'})
        cluster_pcs_df['start'] = cluster['start'].min()
        cluster_pcs_df['end'] = cluster['end'].max()
        cluster_pcs_df['#chr'] = cluster['#chr'].iloc[0]
        # make the right order for bed
        cluster_pcs_dfs.append(make_bed_order(cluster_pcs_df))

    # sorting required by tensorqtl
    cluster_pcs_df = pd.concat(cluster_pcs_dfs)
    cluster_pcs_df = cluster_pcs_df.sort_values(['#chr', 'start', 'end'])

    # occasionally we get inf for all the values in a row.
    # Drop these as they cause susie to error
    cluster_pcs_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    print('Dropped {} rows due to inf'.format(sum(cluster_pcs_df.isna().sum(axis=1) > 0)))
    cluster_pcs_df.dropna(inplace=True)

    # write out bed pc file
    if verb:
        print('Writing out to {}'.format(pc_out_path))
    cluster_pcs_df.to_csv(pc_out_path, sep='\t', index=False)

def main():
    # Parse arguments from cmd
    parser = argparse.ArgumentParser()
    parser.add_argument('-cl', '--cluster_path', help = 'path to .csv clusters')
    parser.add_argument('-e', '--expression_path', help = 'path to .bed normalized expression')
    parser.add_argument('-co', '--covariates_path', help = 'path to covariates')
    parser.add_argument('-o', '--pc_out_path', help = 'path where pc results should be written out')
    parser.add_argument('--verbosity', type=int, default=0, help = 'output verbosity')

    args = parser.parse_args()
    # call the pc funciton
    pc_bed(args.cluster_path, args.expression_path, args.covariates_path, args.pc_out_path, verb=args.verbosity)

if __name__ == "__main__":
    main()