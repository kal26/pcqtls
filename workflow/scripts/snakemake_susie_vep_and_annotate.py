import pandas as pd
from annotate_qtls import *


# parse snakemake
e_susie_path = snakemake.input[0]
pc_susie_path = snakemake.input[1]
vep_path = snakemake.input[2]
annot_pc_path = snakemake.input[3]
e_nominal_path_list = snakemake.input[4:26]
pc_nominal_path_list = snakemake.input[26:]

print(e_nominal_path_list)
print(pc_nominal_path_list)

tissue_id = snakemake.params[0]
gencode_path=snakemake.params[1]
full_abc_path = snakemake.params[2]
abc_match_path=snakemake.params[3]
ctcf_match_path=snakemake.params[4]
ctcf_dir=snakemake.params[5]
tad_path=snakemake.params[6]
avg_expression_path=snakemake.params[7]

out_path = snakemake.output[0]

# load in susie data
pc_susie_df = pd.read_csv(pc_susie_path, sep='\t', index_col=0)
pc_susie_df['cluster_id'] = pc_susie_df['phenotype_id'].str.split('_pc').str[0]
e_susie_df =  pd.read_csv(e_susie_path, sep='\t', index_col=0)
e_susie_df['cluster_id'] = e_susie_df['phenotype_id'].str.split('_e').str[0]
combined_susie = pd.concat([e_susie_df, pc_susie_df], names=['type', 'idx'], keys=['eqtl', 'pcqtl']).reset_index(drop=0).drop(columns=['idx'])
# we will get af again from the nominals
combined_susie.drop(columns=['af'], inplace=True)

# load in the nominal data
e_nominal = pd.concat([pd.read_parquet(path) for path in e_nominal_path_list])
pc_nominal = pd.concat([pd.read_parquet(path) for path in pc_nominal_path_list])
qtls_nominal_merged = pd.merge(combined_susie, pd.concat([pc_nominal, e_nominal]), left_on=['phenotype_id', 'variant_id'],  right_on=['phenotype_id', 'variant_id'], how='left')
qtls_nominal_merged['qtl_variance'] = qtls_nominal_merged['slope'].apply(np.square) * 100
qtls_nominal_merged = qtls_nominal_merged.rename(columns={'slope':'qtl_slope', 'slope_se':'qtl_slope_se'})

# add e nominal info
e_nominal_sub = e_nominal[e_nominal['variant_id'].isin(qtls_nominal_merged['variant_id'])]
e_nominal_sub['cluster_id'] = e_nominal_sub['phenotype_id'].str.split('_e').str[0]
e_nominal_sub['egene_id'] = e_nominal_sub['phenotype_id'].str.split('_e_').str[1]
e_nominal_sub['variance'] = e_nominal_sub['slope'].apply(np.square) * 100
egene_nominal = e_nominal_sub.groupby(['cluster_id', 'variant_id']).agg({'variance':list, 'egene_id':list, 'slope':list, 'slope_se':list})
egene_nominal = egene_nominal.rename(columns={'variance':'egene_variance_list', 'egene_id':'egene_id_list', 'slope':'egene_qtl_slope', 'slope_se':'egene_qtl_slope_se'})
qtls_nominal_merged = pd.merge(qtls_nominal_merged, egene_nominal, left_on=['cluster_id', 'variant_id'], right_index=True, how='left')

# add in annoated pc info
annotated_pcs = pd.read_csv(annot_pc_path, sep='\t')
annotated_pcs = annotated_pcs.rename(columns={'egene_r2':'egene_pc_r2', 
                         'egene_slope':'egene_pc_slope'})
qtl_annot_pc_merged = pd.merge(qtls_nominal_merged.explode(['egene_id_list', 'egene_qtl_slope', 'egene_variance_list', 'egene_qtl_slope_se']), 
                               annotated_pcs[['pc_phenotype_id', 'egene_id', 'egene_pc_r2', 'egene_pc_slope']], 
                               right_on=['pc_phenotype_id', 'egene_id'], 
                               left_on=['phenotype_id', 'egene_id_list'], 
                               how='left')
qtl_annot_pc_merged = qtl_annot_pc_merged.drop(columns=['pc_phenotype_id', 'egene_id'])
# sign flipped qtl slopes

qtl_annot_pc_merged['egene_qtl_slope_flipped'] = qtl_annot_pc_merged['egene_qtl_slope'] * qtl_annot_pc_merged['egene_pc_slope']

# regroup for each variant
qtls_nominal_merged = qtl_annot_pc_merged.groupby(['phenotype_id', 'variant_id','cs_id']).agg({'type':'first',
                                                                         'pip':'first', 
                                                                         'cluster_id':'first',
                                                                         'start_distance':'first',
                                                                         'end_distance':'first', 
                                                                         'af':'first',
                                                                         'ma_samples':'first',
                                                                         'ma_count':'first',
                                                                         'pval_nominal':'first',
                                                                         'qtl_slope':'first',
                                                                         'qtl_slope_se':'first',
                                                                         'qtl_variance':'first',
                                                                         'egene_variance_list':list,
                                                                         'egene_id_list':list,
                                                                         'egene_qtl_slope':list,
                                                                         'egene_qtl_slope_se': list,
                                                                         'egene_pc_r2':list,
                                                                         'egene_pc_slope':list,
                                                                         'egene_qtl_slope_flipped':list}).reset_index()

# add in vep info to each variant
vep = pd.read_csv(vep_path, skiprows=4, sep='\t')
qtls = pd.merge(qtls_nominal_merged, vep[['ID', 'INFO', 'POS', '#CHROM']], left_on='variant_id', right_on='ID', how='left').drop(columns=['ID'])
qtls.rename(columns={'INFO':'vep_info', 
                     'POS':'position', 
                     '#CHROM':'chr'}, inplace=True)

# this is onyl meant to be for one tissue, 
# I only need a tissue ID for the ABC and CTCF tissue matching
qtls['tissue_id'] = tissue_id

# add in the other qtl annotaitons
qtls['Transcripts'] = qtls['cluster_id'].str.replace('_', ',')
qtls = load_and_annotate(qtls, tissue_id,
                      gencode_path=gencode_path,
                      full_abc_path = full_abc_path,
                      abc_match_path=abc_match_path,
                      ctcf_match_path=ctcf_match_path,
                      ctcf_dir=ctcf_dir,
                      tad_path=tad_path,
                      avg_expression_path=avg_expression_path, verbosity=1)

# write out
qtls.to_csv(out_path, index=False, sep='\t')


