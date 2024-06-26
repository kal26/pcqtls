import numpy as np
import pandas as pd
import re

# load data
e_susie_df = pd.read_csv(snakemake.input[0], sep='\t', index_col=0)
pc_susie_df = pd.read_csv(snakemake.input[1], sep='\t', index_col=0)

print('Overlapping {} eqtls and {} pcqtls'.format(len(e_susie_df), len(pc_susie_df)))


# add ids for each credible set
pc_susie_df['cs_full_id'] = pc_susie_df['phenotype_id'].astype(str) + '_cs' + pc_susie_df['cs_id'].astype(str)
e_susie_df['cs_full_id'] = e_susie_df['phenotype_id'].astype(str) + '_e_cs' + e_susie_df['cs_id'].astype(str) 

def get_lead_var(susie_df):
    return susie_df.loc[susie_df.groupby('cs_full_id')['pip'].idxmax(),['cs_full_id','variant_id']].set_index('cs_full_id')

def generate_overlap_df(susie_df):
    overlap_df = pd.DataFrame(pd.Series(susie_df.groupby(['cs_full_id'])['variant_id'].apply(list), name='variant_list'))
    # add in the top variant id
    overlap_df['lead_variant_id'] = get_lead_var(susie_df)
    overlap_df['cluster_id'] = overlap_df.index.str.split('_e|_pc').str[0]
    return overlap_df

# make dfs with one row per credible set to annotate with overlap informaiton
e_overlap_df = generate_overlap_df(e_susie_df)
e_overlap_df['orig_cs_dataset'] = 'control_eqtl'
pc_overlap_df = generate_overlap_df(pc_susie_df)
pc_overlap_df['orig_cs_dataset'] = 'pc_qtl'


def shared_lead_var(row, overlap_df, match=True):
    if match:
        matched_cluster_df = overlap_df[overlap_df['cluster_id']==row.cluster_id]
        return matched_cluster_df[matched_cluster_df['lead_variant_id'] == row.lead_variant_id].index.values
    else:
        return overlap_df[overlap_df['lead_variant_id'] == row.lead_variant_id].index.values

def num_shared_lead_var(row, overlap_df, match=True):
    return len(shared_lead_var(row, overlap_df, match=match))

def generate_csoverlap_dict(susie_df):
    cs_overlap_dict = susie_df.groupby('variant_id')['cs_full_id'].apply(set).to_dict()
    return cs_overlap_dict

def get_overlap_ids_optimized(row, overlap_dict, match=True):
    # check if there is anything to overlap with
    if len(overlap_dict)==0:
        return []
    
    overlap_ids = set()
    for variant_id in row.variant_list:
        if variant_id in overlap_dict:
            overlap_ids.update(overlap_dict[variant_id])
    # check all are in the right cluster
    overlap_ids = list(overlap_ids)
    if match:
        in_same_cluster = [row['cluster_id'] == re.split('_e|_pc', phenotype_id)[0] for phenotype_id in overlap_ids]
        return list(pd.Series(overlap_ids)[in_same_cluster])
    else:
        return list(overlap_ids)

def num_csoverlap(row, overlap_dict, match=True):
    return len(get_overlap_ids_optimized(row, overlap_dict, match=match))

# dict for faster indexing
e_cs_overlap_dict = generate_csoverlap_dict(e_susie_df)
pc_cs_overlap_dict = generate_csoverlap_dict(pc_susie_df)

def annotate_overlap_df(overlap_df):
    # get the overlap with other css with the same lead variant
    overlap_df['e_samelead'] = overlap_df.apply(shared_lead_var, axis=1, args=(e_overlap_df,))
    overlap_df['pc_samelead'] = overlap_df.apply(shared_lead_var, axis=1, args=(pc_overlap_df,))
    overlap_df['num_e_samelead'] = overlap_df.apply(num_shared_lead_var, axis=1, args=(e_overlap_df,))
    overlap_df['num_pc_samelead'] = overlap_df.apply(num_shared_lead_var, axis=1, args=(pc_overlap_df,))
    # get the other css with the any credible set variant overlap 
    print('overlap_dict')
    print(e_cs_overlap_dict)
    print(pc_cs_overlap_dict)
    print('test overlap run')
    print(overlap_df.apply(get_overlap_ids_optimized, axis=1, args=(e_cs_overlap_dict,)))
    overlap_df['e_overlap'] = overlap_df.apply(get_overlap_ids_optimized, axis=1, args=(e_cs_overlap_dict,))
    overlap_df['pc_overlap'] = overlap_df.apply(get_overlap_ids_optimized, axis=1, args=(pc_cs_overlap_dict,))
    overlap_df['num_e_overlap'] = overlap_df.apply(num_csoverlap, axis=1, args=(e_cs_overlap_dict,))
    overlap_df['num_pc_overlap'] = overlap_df.apply(num_csoverlap, axis=1, args=(pc_cs_overlap_dict,))

    # overlap with same lead, not necessarily matching (for debugging)
    overlap_df['e_samelead_all'] = overlap_df.apply(shared_lead_var, axis=1, args=(e_overlap_df,), match=False)
    overlap_df['pc_samelead_all'] = overlap_df.apply(shared_lead_var, axis=1, args=(pc_overlap_df,), match=False)
    overlap_df['num_e_samelead_all'] = overlap_df.apply(num_shared_lead_var, axis=1, args=(e_overlap_df,), match=False)
    overlap_df['num_pc_samelead_all'] = overlap_df.apply(num_shared_lead_var, axis=1, args=(pc_overlap_df,), match=False)
    # overlap any cs variant, not necessarily matching (for debugging)
    overlap_df['e_overlap_all'] = overlap_df.apply(get_overlap_ids_optimized, axis=1, args=(e_cs_overlap_dict,), match=False)
    overlap_df['pc_overlap_all'] = overlap_df.apply(get_overlap_ids_optimized, axis=1, args=(pc_cs_overlap_dict,), match=False)
    overlap_df['num_e_overlap_all'] = overlap_df.apply(num_csoverlap, axis=1, args=(e_cs_overlap_dict,), match=False)
    overlap_df['num_pc_overlap_all'] = overlap_df.apply(num_csoverlap, axis=1, args=(pc_cs_overlap_dict,), match=False)



# add overlap information
annotate_overlap_df(e_overlap_df)
annotate_overlap_df(pc_overlap_df)


# combine and write out
full_overlap = pd.concat([e_overlap_df, pc_overlap_df])
full_overlap.to_csv(snakemake.output[0], sep='\t')