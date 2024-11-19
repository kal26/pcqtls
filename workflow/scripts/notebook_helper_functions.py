# functions to load the requested file types given a config pointing the right directories
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt 
import seaborn as sns 
import networkx as nx


# working directory prefix
prefix = '/home/klawren/oak/pcqtls'

def load_vep(config, tissue_id):
    sample_vep = pd.read_csv('{}/{}/{}.v8.leadvars.vep.vcf'.format(prefix, config['annotations_output_dir'], tissue_id), skiprows=4, sep='\t')
    overlap_df = load_overlap(config, tissue_id)
    return pd.merge(sample_vep, overlap_df, left_on='ID', right_on='lead_variant_id', how='outer')

def load_susie_annotated(config, tissue_id):
    qtl_path = '{}/{}{}/{}.v8.susie_R_vars.annotated.csv'.format(prefix, config['annotations_output_dir'], tissue_id, tissue_id)
    susie_annotated = pd.read_csv(qtl_path, sep='\t')
    susie_annotated['cs_num'] = susie_annotated['cs_id'] 
    susie_annotated['cs_id'] = susie_annotated['phenotype_id'] + '_cs_' + susie_annotated['cs_id'].astype(str)
    susie_annotated = add_lead_var(susie_annotated)
    return susie_annotated

def load_overlap(config, tissue_id):
    overlap_df = pd.read_csv('{}/{}/{}.v8.overlap.txt'.format(prefix, config['overlap_output_dir'], tissue_id), sep='\t')
    overlap_df['var_cluster'] = overlap_df['lead_variant_id'] + '_' + overlap_df['cluster_id']
    # add pc num column
    overlap_df['pc_num'] = np.where(overlap_df['orig_cs_dataset']=='pc_qtl',  overlap_df['cs_full_id'].str.split('_').str[-2].str.strip('pc'), 0)
    overlap_df['pc_num'] = overlap_df['pc_num'].astype(int)
    # add in cluster size
    overlap_df['cluster_size'] = overlap_df['cluster_id'].str.split('_').apply(len)
    # add in the variant position as a column
    overlap_df['lead_variant_pos'] = overlap_df['lead_variant_id'].str.split('_').str[1].astype(int)
    # split first, last, and middle pcs
    overlap_df['pc_order'] = 'middle'
    overlap_df.loc[overlap_df['pc_num'] == overlap_df['cluster_size'],'pc_order'] = 'last'
    overlap_df.loc[overlap_df['pc_num'] == 1,'pc_order'] = 'first'
    overlap_df.loc[overlap_df['orig_cs_dataset'] == 'control_eqtl','pc_order'] = 'eqtl'
    
    return overlap_df

def load_clusters_annotated(config, tissue_id):
    annot_cluster = pd.read_csv('{}/{}/{}/{}_clusters_annotated.csv'.format(prefix, config['annotations_output_dir'], tissue_id, tissue_id), index_col=0)
    for idx, row in annot_cluster.iterrows():
        annot_cluster.loc[idx, 'cluster_id']  = '_'.join([*sorted(row['Transcripts'].split(','))])
    annot_cluster['abs_cor'] = np.where(annot_cluster['Mean_neg_cor'].isna(), annot_cluster['Mean_pos_cor'], - annot_cluster['Mean_neg_cor'])
    annot_cluster['abs_cor'] = np.where(annot_cluster['abs_cor'] < annot_cluster['Mean_pos_cor'], annot_cluster['Mean_pos_cor'], annot_cluster['abs_cor'])
    annot_cluster['max_pos_neg_cor'] = np.where(annot_cluster['Mean_neg_cor'].isna(), annot_cluster['Mean_pos_cor'], annot_cluster['Mean_neg_cor'])
    annot_cluster['max_pos_neg_cor'] = np.where(abs(annot_cluster['max_pos_neg_cor']) < annot_cluster['Mean_pos_cor'], annot_cluster['Mean_pos_cor'], annot_cluster['max_pos_neg_cor'])
    return annot_cluster

def load_null_clusters_annotated(config, tissue_id, num_genes=2):
    return pd.read_csv('{}/{}/{}/{}_null_{}genes_annotated.csv'.format(prefix, config['annotations_output_dir'], tissue_id, tissue_id, num_genes), index_col=0)

def load_cluster(config, tissue_id):
    cluster_df =  pd.read_csv('{}/{}/{}_clusters_all_chr.csv'.format(prefix, config['clusters_dir'], tissue_id),index_col=0)
    for idx, row in cluster_df.iterrows():
        cluster_df.loc[idx, 'cluster_id']  = '_'.join([*sorted(row['Transcripts'].split(','))])
    return cluster_df

def load_pc_cis(config, tissue_id):
    pc_cis_path = '{}/{}/{}/{}.v8.pcs.cis_qtl.txt.gz'.format(prefix, config['pcqtl_output_dir'], tissue_id, tissue_id)
    pc_cis_df = pd.read_csv(pc_cis_path, sep='\t')
    pc_cis_df['cluster_id'] = pc_cis_df['phenotype_id'].str.split('_pc').str[0]
    annotate_pc_order(pc_cis_df)
    return pc_cis_df

def load_e_cis(config, tissue_id):
    e_cis_path = '{}/{}/{}/{}.v8.cluster_genes.cis_qtl.txt.gz'.format(prefix, config['eqtl_output_dir'], tissue_id, tissue_id)
    e_cis_df = pd.read_csv(e_cis_path, sep='\t')
    e_cis_df['cluster_id'] = e_cis_df['phenotype_id'].str.split('_e').str[0]
    return e_cis_df

# load in e nominal
def load_e_nominal(config, tissue_id, chr_id=22, get_var_position=False):
    eqtl_output_dir = config['eqtl_output_dir']
    path = f'{prefix}/{eqtl_output_dir}/{tissue_id}/{tissue_id}.v8.cluster_genes.cis_qtl_pairs.chr{chr_id}.parquet'
    e_nominal_df = pd.read_parquet(path)
    if get_var_position:
        e_nominal_df['variant_pos'] = var_pos(e_nominal_df)
    e_nominal_df['cluster_id'] = e_nominal_df['phenotype_id'].str.split('_e_').str[0]
    e_nominal_df['var_cluster'] = e_nominal_df['variant_id'] + '_' + e_nominal_df['cluster_id']
    e_nominal_df['egene_id'] = e_nominal_df['phenotype_id'].str.split('_e_').str[1]
    return e_nominal_df

def load_pc_nominal(config, tissue_id, chr_id=22, get_var_position=False):
    pcqtl_output_dir = config['pcqtl_output_dir']
    path = f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.cis_qtl_pairs.chr{chr_id}.parquet'
    pc_nominal_df = pd.read_parquet(path)
    if get_var_position:
        pc_nominal_df['variant_pos'] = var_pos(pc_nominal_df)
    pc_nominal_df['cluster_id'] = pc_nominal_df['phenotype_id'].str.split('_pc').str[0]
    pc_nominal_df['var_cluster'] = pc_nominal_df['variant_id'] + '_' + pc_nominal_df['cluster_id']
    pc_nominal_df['pc_num'] = pc_nominal_df['phenotype_id'].str.split('_pc').str[-1].astype(int)
    return pc_nominal_df

def load_pc_susie_r(config, tissue_id):
    return load_pc_susie(config, tissue_id, use_r=True)

def load_pc_susie(config, tissue_id, use_r=False):
    pcqtl_output_dir = config['pcqtl_output_dir']
    if use_r:
        pc_susie_path = f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.susie_R.txt'
    else:
        pc_susie_path = f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.susie.txt'
    pc_susie_df =  pd.read_csv(pc_susie_path, sep='\t', index_col=0).reset_index()
    pc_susie_df['var_cluster'] = pc_susie_df['variant_id'] + '_' + pc_susie_df['phenotype_id'].str.split('_pc').str[0]
    pc_susie_df['cs_num'] = pc_susie_df['cs_id'] 
    pc_susie_df['cs_id'] = pc_susie_df['phenotype_id'] + '_cs_' + pc_susie_df['cs_id'].astype(str)
    pc_susie_df['pc_num'] = pc_susie_df['phenotype_id'].str.split('_pc').str[-1].astype(int)
    pc_susie_df['cluster_id'] = pc_susie_df['phenotype_id'].str.split('_pc').str[0]
    pc_susie_df = add_lead_var(pc_susie_df)
    pc_susie_df = add_num_vars_cs(pc_susie_df)
    return pc_susie_df

def load_e_susie_r(config, tissue_id):
    return load_e_susie(config, tissue_id, use_r=True)

def load_e_susie(config, tissue_id, use_r=False):
    eqtl_output_dir = config['eqtl_output_dir']
    if use_r:
        e_susie_path = f'{prefix}/{eqtl_output_dir}/{tissue_id}/{tissue_id}.v8.cluster_genes.susie_R.txt'
    else:
        e_susie_path = f'{prefix}/{eqtl_output_dir}/{tissue_id}/{tissue_id}.v8.cluster_genes.susie.txt'
    e_susie_df =  pd.read_csv(e_susie_path, sep='\t', index_col=0).reset_index()
    e_susie_df['var_cluster'] = e_susie_df['variant_id'] + '_' + e_susie_df['phenotype_id'].str.split('_e').str[0]
    e_susie_df['cs_num'] = e_susie_df['cs_id'] 
    e_susie_df['cs_id'] = e_susie_df['phenotype_id'] + '_cs_' + e_susie_df['cs_id'].astype(str)
    e_susie_df['cluster_id'] = e_susie_df['phenotype_id'].str.split('_e').str[0]
    e_susie_df = add_lead_var(e_susie_df)
    e_susie_df = add_num_vars_cs(e_susie_df)
    return e_susie_df


def load_pairwise_coloc(config, tissue_id):
    pair_coloc = []
    for chr_id in range(1,23):
        pair_coloc_path = "{}/{}/pairs/{}.v8.pairs_coloc.chr{}.txt".format(prefix, config["coloc_output_dir"], tissue_id, chr_id)
        pair_coloc.append(pd.read_csv(pair_coloc_path, sep='\t'))
    pair_coloc = pd.concat(pair_coloc)
    pair_coloc['cs_id_1'] = pair_coloc['qtl1_id'] + '_cs_' + pair_coloc['idx1'].astype(str)
    pair_coloc['cs_id_2'] = pair_coloc['qtl2_id'] + '_cs_' + pair_coloc['idx2'].astype(str)
    return pair_coloc

def load_tissue_ids(config):
    tissue_id_path = config['tissue_id_path']
    tissue_df = pd.read_csv(f"{prefix}/{tissue_id_path}", header=0)
    return list(tissue_df['Tissue'])

def load_tissue_df(config):
    tissue_id_path = config['tissue_id_path']
    return pd.read_csv(f"{prefix}/{tissue_id_path}", header=0)


def load_pc_top_per_phenotype(config, tissue_id):
    pcqtl_output_dir = config['pcqtl_output_dir']
    return pd.read_csv(f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.cis_qtl.txt.gz', sep='\t')

def load_e_top_per_phenotype(config, tissue_id):
    eqtl_output_dir = config['eqtl_output_dir']
    return pd.read_csv(f'{prefix}/{eqtl_output_dir}/{tissue_id}/{tissue_id}.v8.cluster_genes.cis_qtl.txt.gz', sep='\t')

def load_pc_permutation(config, tissue_id, pc1_only=False):
    pcqtl_output_dir = config['pcqtl_output_dir']
    if pc1_only:
        pc_path = f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pc1_only.cis_independent_qtl.txt.gz'
    else:
        pc_path = f'{prefix}/{pcqtl_output_dir}/{tissue_id}/{tissue_id}.v8.pcs.cis_independent_qtl.txt.gz'
    pc_permutation_df = pd.read_csv(pc_path, sep='\t')
    pc_permutation_df['cluster_id'] = pc_permutation_df['phenotype_id'].str.split('_pc').str[0]
    return pc_permutation_df

def load_e_permutation(config, tissue_id):
    eqtl_output_dir = config['eqtl_output_dir']
    e_permutation_df = pd.read_csv(f'{prefix}/{eqtl_output_dir}/{tissue_id}/{tissue_id}.v8.cluster_genes.cis_independent_qtl.txt.gz', sep='\t')
    e_permutation_df['cluster_id'] = e_permutation_df['phenotype_id'].str.split('_e').str[0]
    return e_permutation_df


def load_e_nominal_all_chr(config, tissue_id):
    e_nominal_dfs=[]
    for chr_id in range(1,23):
        e_nominal_dfs.append(load_e_nominal(config, tissue_id, chr_id=chr_id))
    return pd.concat(e_nominal_dfs)

def load_pc_nominal_all_chr(config, tissue_id):
    pc_nominal_dfs=[]
    for chr_id in range(1,23):
        pc_nominal_dfs.append(load_pc_nominal(config, tissue_id, chr_id=chr_id))
    return pd.concat(pc_nominal_dfs)

def load_expression(config, tissue_id):
    expression_path = '/{}/{}/{}.v8.normalized_expression.bed'.format(prefix, config['expression_dir'], tissue_id)
    expression_df = pd.read_csv(expression_path, sep='\t')
    return expression_df

def load_cluster_expression(config, tissue_id):
    expression_path = '/{}/{}/{}.v8.normalized_residualized_expression.cluster_genes.bed'.format(prefix, config['filtered_expression_output_dir'], tissue_id)
    expression_df = pd.read_csv(expression_path, sep='\t')
    expression_df['cluster_id'] = expression_df['gene_id'].str.split('_e_').str[0]
    expression_df['egene_id'] = expression_df['gene_id'].str.split('_e_').str[1]
    return expression_df

def load_pc(config, tissue_id):
    pc_df = pd.read_csv('{}/{}/{}.pcs.bed'.format(prefix, config['pc_output_dir'], tissue_id), sep='\t')
    pc_df['cluster_id'] = pc_df['gene_id'].str.split('_pc').str[0]
    pc_df['pc_id'] = pc_df['gene_id'].str.split('_pc').str[1].astype('float')
    pc_df['cluster_size'] = pc_df['cluster_id'].str.split('_').apply(len)
    return pc_df


def load_signal_groups(config, tissue_id):
    pair_coloc = load_pairwise_coloc(config, tissue_id)
    pc_susie_r = load_pc_susie_r(config, tissue_id)
    e_susie_r = load_e_susie_r(config, tissue_id)
    signal_groups = get_signal_groups_tissue(pair_coloc, pc_susie_r, e_susie_r)
    signal_groups['tissue_id'] = tissue_id
    return signal_groups


# functions to help with annotating loaded data
def var_pos(df, column='variant_id'):
    return df[column].str.split('_').str[1].astype(int)

def add_lead_var(susie_df):
    lead_vars = susie_df.loc[susie_df.groupby('cs_id')['pip'].idxmax(),['cs_id','variant_id']].set_index('cs_id')
    susie_df = susie_df.merge(lead_vars, how='left', left_on='cs_id', right_index=True)
    susie_df = susie_df.rename(columns={'variant_id_y':'lead_variant_id', 'variant_id_x':'variant_id'})
    return susie_df

def add_num_vars_cs(susie_df):
    num_vars = susie_df.groupby('cs_id').agg({'variant_id':'count'})
    num_vars = num_vars.rename(columns={'variant_id':'num_vars'})
    susie_df = susie_df.merge(num_vars, how='left', left_on='cs_id', right_index=True)
    return susie_df

def load_across_tissues(config, load_func, tissue_ids = None):
    if tissue_ids == None:
        tissue_ids = load_tissue_ids(config)
    combined = [load_func(config, tissue_id) for tissue_id in tissue_ids]
    combined = pd.concat([df.assign(tissue_id=n) for df, n in zip(combined, tissue_ids)])
    combined.reset_index(inplace=True, drop=True)
    return combined

def annotate_pc_order(pc_df):
    pc_df['pc_num'] = pc_df['phenotype_id'].str.split('_pc').str[-1].astype(int)
    pc_df['cluster_id'] = pc_df['phenotype_id'].str.split('_pc').str[0]
    pc_df['cluster_size'] = pc_df['cluster_id'].str.split('_').apply(len)
    # split first, last, and middle pcs
    pc_df['pc_order'] = 'middle'
    pc_df.loc[pc_df['pc_num'] == pc_df['cluster_size'],'pc_order'] = 'last'
    pc_df.loc[pc_df['pc_num'] == 1,'pc_order'] = 'first'





# log odds stuff
# fit a logit model

def fit_logistic_model(joined_df, column, label_col, correct_on=False, filter_on_column=None, verb=True, correct_on_column='log_size'):
    if verb:
        print(column)
    try:
        final_joined_df = joined_df[joined_df[filter_on_column]]
    except KeyError:
        final_joined_df = joined_df

    if correct_on:
        x = final_joined_df[np.append(np.asarray(correct_on_column),column)].astype(float)
    else:
        x = final_joined_df[column].astype(float)

    y = final_joined_df.reset_index()[label_col].values.astype(bool)
    x_with_constant = sm.add_constant(x) # Add intercept term

    logit_model = sm.Logit(y, x_with_constant)
    try:
        result = logit_model.fit()
        
    except np.linalg.LinAlgError:
        # this can happen if all elements are the same
        return pd.Series({'lower_cb':np.nan, 'upper_cb':np.nan, 'odds_ratio':np.nan, 'p_value':np.nan, 'col':column,'lower_cb_diff':np.nan , 'upper_cb_diff':np.nan}, name=column)
    if verb:
        print(result.summary())

    coefficients = result.params
    conf_int = result.conf_int()
    odds_ratios = np.exp(coefficients)
    odds_ratios_ci = np.exp(conf_int)

    odds_ratios_ci.rename(columns={0:'lower_cb', 1:'upper_cb'}, inplace=True)
    odds_ratios_ci['odds_ratio'] = odds_ratios
    odds_ratios_ci['p_value'] = result.pvalues
    odds_ratios_ci['col'] = column
    odds_ratios_ci['lower_cb_diff'] = odds_ratios_ci['odds_ratio'] - odds_ratios_ci['lower_cb']
    odds_ratios_ci['upper_cb_diff'] = odds_ratios_ci['upper_cb'] - odds_ratios_ci['odds_ratio']
    
    return odds_ratios_ci.loc[column]


def get_odds_df(joined_df, label_col, verb=True, correct_on=False, correct_on_column='log_size', 
                column_list = [], filter_list=[], filter_on_column = 'has_multiple_abc_genes'):
    column_list = pd.Series(column_list)
    odds_ratios_no_filter = pd.DataFrame([fit_logistic_model(joined_df, c, label_col, verb=verb, correct_on=correct_on, correct_on_column=correct_on_column) for c in column_list[~column_list.isin(filter_list)]]) 
    if len(filter_list)>0:
        odds_ratios_filtered = pd.DataFrame([fit_logistic_model(joined_df, c, label_col, filter_on_column=filter_on_column, verb=verb, correct_on=correct_on, correct_on_column=correct_on_column) for c in column_list[column_list.isin(filter_list)]]) 
        return pd.concat([odds_ratios_no_filter, odds_ratios_filtered])
    else: 
        return odds_ratios_no_filter

def make_log_odds_plot(log_odds_df, ax=None, add_annotations=True):
    log_odds_df = log_odds_df.reset_index()
    if ax==None:
        fig, ax = plt.subplots(1, figsize=(9,9))

    # log odds plot
    ax.errorbar(y=log_odds_df['col'], x=log_odds_df['odds_ratio'], xerr=log_odds_df[['lower_cb_diff', 'upper_cb_diff']].values.transpose(), fmt="o", color='k')
    ax.axvline(1, color='k', linestyle='--')
    ax.set_xlabel('Log odds')
    if add_annotations:
        for idx,row in log_odds_df.iterrows():
            ax.annotate('OR: {:.2f},\np: {:.1E}'.format(row['odds_ratio'], row['p_value']), (row['odds_ratio'], idx+.2))
    ax.set_xscale(u'log')
    return ax

# log odds plot with multiple odds per category 
def make_log_odds_plot_multiple(odds_ratios_list, ax=None, labels=None, add_annotations=True, offset = 0.2, colors = sns.color_palette()):
    if ax==None:
        fig, ax = plt.subplots(1, figsize=(9,9))
    
    for idx, odds_ratio_df in enumerate(odds_ratios_list):
        odds_ratio_df = odds_ratio_df.reindex(odds_ratios_list[0].index).reset_index()
        color = colors[idx % len(colors)]  # cycle through colors if more than available
        ax.errorbar(y=odds_ratio_df.reset_index().index.values + idx*offset, x=odds_ratio_df['odds_ratio'], 
                    xerr=odds_ratio_df[['lower_cb_diff', 'upper_cb_diff']].values.transpose(), fmt="o", 
                    color=color, label=labels[idx] if labels else None, markersize=3)
        ax.axvline(1, color='k', linestyle='--')

        if add_annotations:
            for row_idx, row in odds_ratio_df.iterrows():
                ax.annotate('OR = {:.2f}, p={:.1E}'.format(row['odds_ratio'], row['p_value']), 
                            (row['odds_ratio'], row_idx + idx*offset + 0.05), fontsize=6)
            
    if labels:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(reversed(handles), reversed(labels))

    ax.set_xscale(u'log')
    ax.set_yticks(ticks=odds_ratio_df.index.values + (len(odds_ratios_list)-1)*offset/2, labels=(odds_ratio_df['col']))
    return ax


def get_signal_groups_tissue(pair_coloc, pc_susie_r, e_susie_r, coloc_cutoff=.75):
    # filter to coloc'd signals
    pair_coloc = pair_coloc[pair_coloc['PP.H4.abf'] > coloc_cutoff]

    # Create an undirected graph
    G = nx.Graph()

    # Add edges to the graph from the DataFrame
    for index, row in pair_coloc.iterrows():
        G.add_edge(row['cs_id_1'], row['cs_id_2'])

    # Get the connected components of the graph
    connected_components = list(nx.connected_components(G))

    # Generate underlying signal ids
    underlying_signals = ['-'.join(sorted(component)) for component in connected_components]

    # add in the unique ones
    all_credible_set_ids = set(pc_susie_r['cs_id']).union(set(e_susie_r['cs_id']))

    # Add standalone sets (those not in any connected component)
    for credible_set_id in all_credible_set_ids:
        if not any(credible_set_id in component for component in connected_components):
            underlying_signals.append(credible_set_id)

    underlying_signals = pd.DataFrame({'signal_id':underlying_signals})
    underlying_signals['num_e_coloc'] = underlying_signals['signal_id'].astype(str).str.count('_e_')
    underlying_signals['num_pc_coloc'] = underlying_signals['signal_id'].astype(str).str.count('_pc')
    underlying_signals['multiple_e'] = underlying_signals['num_e_coloc'] > 1
    underlying_signals['multiple_pc'] = underlying_signals['num_pc_coloc'] > 1
    return underlying_signals

def run_get_signal_groups(pair_coloc, pc_susie_r, e_susie_r, coloc_cutoff=.75):
    # Ensure all inputs have 'tissue_id' column for grouping
    if 'tissue_id' not in pair_coloc.columns or \
       'tissue_id' not in pc_susie_r.columns or \
       'tissue_id' not in e_susie_r.columns:
        raise ValueError("All input DataFrames must have a 'tissue_id' column.")

    # Group by tissue_id
    results = []
    tissue_groups = pair_coloc['tissue_id'].unique()

    for tissue in tissue_groups:
        # Filter the data for the current tissue_id
        pair_coloc_group = pair_coloc[pair_coloc['tissue_id'] == tissue]
        pc_susie_r_group = pc_susie_r[pc_susie_r['tissue_id'] == tissue]
        e_susie_r_group = e_susie_r[e_susie_r['tissue_id'] == tissue]
        
        # Call get_signal_groups_tissue for the current group
        tissue_results = get_signal_groups_tissue(pair_coloc_group, pc_susie_r_group, e_susie_r_group, coloc_cutoff=coloc_cutoff)
        
        # Add tissue_id to the results
        tissue_results['tissue_id'] = tissue
        results.append(tissue_results)

    # Concatenate all results into a single DataFrame
    final_results = pd.concat(results, ignore_index=True)
    return final_results