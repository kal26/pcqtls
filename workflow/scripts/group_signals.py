#!/usr/bin/env python3
"""
Signal Groups Analysis

This script contains functions for analyzing signal groups from co-localization
results, including QTL-QTL and QTL-GWAS signal group identification.
"""
import sys; print(f"Using Python: {sys.executable}")

import logging
import sys
import os
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
from pandas.errors import EmptyDataError

# Import local modules
try:
    from utils import *
except ImportError:
    # Fallback for when script is run directly
    import sys
    sys.path.append(str(Path(__file__).parent))
    from utils import *

# Configuration - use environment variables with defaults
PREFIX = os.getenv('PCQTL_PREFIX', '/home/klawren/oak/pcqtls')


def setup_logging(level=logging.INFO):
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def group_qtl_signals(pair_coloc, pc_susie_r, e_susie_r, coloc_cutoff=.75, get_variants=False):
    """
    Group QTL signals for a single tissue.
    
    This function creates signal groups by connecting co-localized credible sets
    using network analysis to identify connected components.
    
    Args:
        pair_coloc: DataFrame with pairwise co-localization results
        pc_susie_r: DataFrame with pcQTL SuSiE results
        e_susie_r: DataFrame with eQTL SuSiE results
        coloc_cutoff: PP.H4 threshold for co-localization (default: 0.75)
        get_variants: Whether to include variant sets in output (default: False)
    
    Returns:
        DataFrame with QTL signal groups
    """
    # Filter to coloc'd signals
    pair_coloc = pair_coloc[pair_coloc['PP.H4.abf'] > coloc_cutoff]

    # Create an undirected graph
    G = nx.Graph()

    # Create CS IDs for pairwise coloc data
    pair_coloc['cs_id_1'] = pair_coloc['qtl1_id'] + '_cs_' + pair_coloc['idx1'].astype(int).astype(str)
    pair_coloc['cs_id_2'] = pair_coloc['qtl2_id'] + '_cs_' + pair_coloc['idx2'].astype(int).astype(str)
    
    # Add edges to the graph from the DataFrame
    for index, row in pair_coloc.iterrows():
        G.add_edge(row['cs_id_1'], row['cs_id_2'])

    # Get the connected components of the graph
    connected_components = list(nx.connected_components(G))

    # Generate underlying signal ids
    underlying_signals = ['-'.join(sorted(component)) for component in connected_components]

    # Add in the unique ones
    all_credible_set_ids = set(pc_susie_r['cs_id']).union(set(e_susie_r['cs_id']))

    # Add standalone sets (those not in any connected component)
    for credible_set_id in all_credible_set_ids:
        if not any(credible_set_id in component for component in connected_components):
            underlying_signals.append(credible_set_id)

    underlying_signals = pd.DataFrame({'signal_id': underlying_signals})
    # Extract cluster_id from signal_id - extract from the first part before _cs_
    def extract_cluster_id(signal_id):
        # Split by dash and look for elements containing '_cs_'
        parts = signal_id.split('-')
        for part in parts:
            if '_cs_' in part:
                # Extract the part before _cs_ and then get cluster_id
                base_id = part.split('_cs_')[0]
                if '_e_' in base_id:
                    return base_id.split('_e_')[0]
                elif '_pc' in base_id:
                    return base_id.split('_pc')[0]
        return None
    
    underlying_signals['cluster_id'] = underlying_signals['signal_id'].apply(extract_cluster_id)
    underlying_signals['num_e_coloc'] = underlying_signals['signal_id'].astype(str).str.count('_e_')
    underlying_signals['num_pc_coloc'] = underlying_signals['signal_id'].astype(str).str.count('_pc')

    # Add the set of lead variants for all signals in the group
    def get_var_set(row, lead=True):
        var_ids = []
        if lead:
            [var_ids.append(var) for var in pc_susie_r[pc_susie_r['cs_id'].isin(row['signal_id'].split('-'))]['lead_variant_id'].values]
            [var_ids.append(var) for var in e_susie_r[e_susie_r['cs_id'].isin(row['signal_id'].split('-'))]['lead_variant_id'].values]
        else:
            [var_ids.append(var) for var in pc_susie_r[pc_susie_r['cs_id'].isin(row['signal_id'].split('-'))]['variant_id'].values]
            [var_ids.append(var) for var in e_susie_r[e_susie_r['cs_id'].isin(row['signal_id'].split('-'))]['variant_id'].values]
        return list(set(var_ids))
    underlying_signals['lead_var_set'] = underlying_signals.apply(get_var_set, axis=1, args=(True,))
    if get_variants:
        underlying_signals['var_set'] = underlying_signals.apply(get_var_set, axis=1, args=(False,))

    # Reorder columns to match specification
    column_order = ['signal_id', 'cluster_id', 'num_e_coloc', 'num_pc_coloc', 'lead_var_set', 'var_set']
    if get_variants:
        underlying_signals = underlying_signals[column_order]
    else:
        underlying_signals = underlying_signals[column_order[:-1]]  # Remove variant columns if not requested

    return underlying_signals





def load_gwas_coloc(config):
    """
    Load GWAS co-localization results from all tissues.
    
    This function recursively searches for GWAS co-localization files and
    loads them into a single DataFrame with proper annotations.
    
    Args:
        config: Configuration dictionary with paths
    
    Returns:
        DataFrame with all GWAS co-localization results
    """
    # Recursively get list of files ending with .qtl_gwas.txt
    def get_files(directory):
        file_list = []
        for root, directories, files in os.walk(directory):
            if not 'temp' in root:
                for file in files:
                    if "susie_True" in file:
                        file_list.append(os.path.join(root, file))
        return file_list
    
    coloc_file_list = get_files(f'{PREFIX}/{config["coloc_output_dir"]}gwas')

    # Load each file into a DataFrame and concatenate them
    cluster_colocs = []
    for cluster_file in coloc_file_list:
        try:
            cluster_coloc = pd.read_csv(cluster_file, sep='\t')
            cluster_coloc['coloc_file'] = cluster_file
            cluster_colocs.append(cluster_coloc)
        except EmptyDataError as e:
            print(f"File is empty: {cluster_file}") 

    # Concatenate all DataFrames into a single DataFrame
    gwas_coloc = pd.concat(cluster_colocs, ignore_index=True)
    # Drop duplicate rows from intermediate write outs
    gwas_coloc = gwas_coloc.drop_duplicates()

    # Add information to dataframe
    gwas_coloc['cluster_id'] = np.where(gwas_coloc['qtl_id'].str.contains('_e'), 
                                       gwas_coloc['qtl_id'].str.split('_e').str[0], 
                                       gwas_coloc['qtl_id'].str.split('_pc').str[0])
    
    # Make ids for each credible set in the qtl and gwas
    gwas_coloc['gwas_cs_id'] = gwas_coloc['gwas_id'] + '_cs_' + gwas_coloc['idx1'].astype(int).astype(str) + '_cluster_' + gwas_coloc['cluster_id']
    gwas_coloc['qtl_cs_id'] = gwas_coloc['qtl_id'] + '_cs_' + gwas_coloc['idx2'].astype(int).astype(str) + '_cluster_' + gwas_coloc['cluster_id']
    
    # Set type as pcqtl or eqtl
    gwas_coloc['type'] = np.where(gwas_coloc['qtl_cs_id'].str.contains('_pc'), 'pcqtl', 'eqtl')
    
    # Fill in nas with 0
    gwas_coloc[['PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf']] = gwas_coloc[['PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf']].fillna(0)

    # Add cs ids
    gwas_coloc['gwas_cs_id'] = 'gwas_' + gwas_coloc['gwas_cs_id']
    gwas_coloc['qtl_cs_id'] = 'qtl_' + gwas_coloc['qtl_cs_id']
    
    return gwas_coloc


def get_gwas_signals(gwas_coloc_hits, pair_coloc_hits):
    """
    Identify GWAS signal groups from co-localization results.
    
    This function creates signal groups that include both QTL-QTL and QTL-GWAS
    co-localizations using network analysis.
    
    Args:
        gwas_coloc_hits: DataFrame with GWAS co-localization results
        pair_coloc_hits: DataFrame with pairwise QTL co-localization results
    
    Returns:
        DataFrame with GWAS signal groups and their properties
    """
    # Create CS IDs for pairwise coloc data
    pair_coloc_hits['cs_id_1'] = pair_coloc_hits['qtl1_id'] + '_cs_' + pair_coloc_hits['idx1'].astype(int).astype(str)
    pair_coloc_hits['cs_id_2'] = pair_coloc_hits['qtl2_id'] + '_cs_' + pair_coloc_hits['idx2'].astype(int).astype(str)

    # Create an undirected graph
    G = nx.Graph()
    
    # Add an edge for each gwas-qtl coloc
    for index, row in gwas_coloc_hits.iterrows():
        G.add_edge(row['gwas_cs_id'], row['qtl_cs_id'])

    # Add an edge for each qtl-qtl coloc
    for index, row in pair_coloc_hits.iterrows():
        G.add_edge(row['cs_id_1'], row['cs_id_2'])

    # Get the connected components of the graph
    connected_components = list(nx.connected_components(G))
    
    # Generate underlying signal ids
    underlying_signals = ['-'.join(sorted(component)) for component in connected_components]

    # Make a df
    underlying_signals = pd.DataFrame({'signal_id': underlying_signals})
    # Extract cluster_id from signal_id - extract from the first part before _cs_
    def extract_cluster_id(signal_id):
        # Split by dash and look for elements containing '_cs_'
        parts = signal_id.split('-')
        for part in parts:
            if '_cs_' in part:
                # Extract the part before _cs_ and then get cluster_id
                base_id = part.split('_cs_')[0]
                if '_e_' in base_id:
                    return base_id.split('_e_')[0]
                elif '_pc' in base_id:
                    return base_id.split('_pc')[0]
        return None
    
    underlying_signals['cluster_id'] = underlying_signals['signal_id'].apply(extract_cluster_id)
    underlying_signals['num_qtl_coloc'] = underlying_signals['signal_id'].astype(str).str.count('qtl_')
    underlying_signals['num_gwas_coloc'] = underlying_signals['signal_id'].astype(str).str.count('gwas_')
    underlying_signals['num_e_coloc'] = underlying_signals['signal_id'].astype(str).str.count('_e_')
    underlying_signals['num_pc_coloc'] = underlying_signals['signal_id'].astype(str).str.count('_pc')
    
    # Reorder columns to match specification
    column_order = ['signal_id', 'cluster_id', 'num_qtl_coloc', 'num_gwas_coloc', 'num_e_coloc', 'num_pc_coloc']
    underlying_signals = underlying_signals[column_order]
    
    return underlying_signals


def group_gwas_signals(gwas_coloc_files, pair_coloc_files, output_file, coloc_cutoff=0.75):
    """
    Group GWAS signals from co-localization files.
    
    This function loads GWAS and pairwise co-localization results and identifies
    signal groups that include both QTL-QTL and QTL-GWAS co-localizations.
    
    Args:
        gwas_coloc_files: Path(s) to GWAS co-localization results file(s)
        pair_coloc_files: Path(s) to pairwise co-localization results file(s)
        output_file: Path to save GWAS signal groups results
        coloc_cutoff: PP.H4 threshold for co-localization (default: 0.75)
    """
    logger = logging.getLogger(__name__)
    
    try:
        logger.info('Loading GWAS and pairwise co-localization results...')
        
        # Load data
        # Handle multiple gwas_coloc files (per GWAS)
        if isinstance(gwas_coloc_files, str):
            gwas_coloc_files = [gwas_coloc_files]
            
        if len(gwas_coloc_files) == 1:
            gwas_coloc = pd.read_csv(gwas_coloc_files[0], sep='\t')
        else:
            logger.info(f'Loading {len(gwas_coloc_files)} gwas_coloc files...')
            gwas_coloc_dfs = []
            for file_path in gwas_coloc_files:
                if os.path.exists(file_path):
                    df = pd.read_csv(file_path, sep='\t')
                    gwas_coloc_dfs.append(df)
                else:
                    logger.warning(f'File not found: {file_path}')
            
            if gwas_coloc_dfs:
                gwas_coloc = pd.concat(gwas_coloc_dfs, ignore_index=True)
                logger.info(f'Combined {len(gwas_coloc_dfs)} GWAS files into {len(gwas_coloc)} rows')
            else:
                raise ValueError("No valid gwas_coloc files found")
        
        # Handle multiple pair_coloc files (per chromosome)
        if isinstance(pair_coloc_files, str):
            pair_coloc_files = [pair_coloc_files]
            
        if len(pair_coloc_files) == 1:
            pair_coloc = pd.read_csv(pair_coloc_files[0], sep='\t')
        else:
            logger.info(f'Loading {len(pair_coloc_files)} pair_coloc files...')
            pair_coloc_dfs = []
            for file_path in pair_coloc_files:
                if os.path.exists(file_path):
                    df = pd.read_csv(file_path, sep='\t')
                    pair_coloc_dfs.append(df)
                else:
                    logger.warning(f'File not found: {file_path}')
            
            if pair_coloc_dfs:
                pair_coloc = pd.concat(pair_coloc_dfs, ignore_index=True)
                logger.info(f'Combined {len(pair_coloc_dfs)} files into {len(pair_coloc)} rows')
            else:
                raise ValueError("No valid pair_coloc files found")
        
        logger.info('Filtering to significant co-localizations...')
        
        # Filter to significant co-localizations
        gwas_coloc_hits = gwas_coloc[gwas_coloc['PP.H4.abf'] > coloc_cutoff]
        pair_coloc_hits = pair_coloc[pair_coloc['PP.H4.abf'] > coloc_cutoff]
        
        logger.info(f'Found {len(gwas_coloc_hits)} significant GWAS co-localizations')
        logger.info(f'Found {len(pair_coloc_hits)} significant pairwise co-localizations')
        
        logger.info('Running GWAS signal groups analysis...')
        
        # Run GWAS signal groups analysis
        gwas_signal_groups = get_gwas_signals(gwas_coloc_hits, pair_coloc_hits)
        
        logger.info(f'Saving results to {output_file}')
        gwas_signal_groups.to_csv(output_file, sep='\t', index=False)
        
        logger.info('GWAS signal groups analysis completed successfully')
        
    except Exception as e:
        logger.error(f'Error in GWAS signal groups analysis: {e}')
        raise


def main():
    """Main function for command-line usage."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze signal groups from co-localization results')
    parser.add_argument('--mode', choices=['qtl', 'gwas'], required=True, 
                       help='Analysis mode: qtl for QTL signal groups, gwas for GWAS signal groups')
    parser.add_argument('--pair-coloc', nargs='+', help='Path(s) to pairwise co-localization file(s)')
    parser.add_argument('--pc-susie', help='Path to pcQTL SuSiE results')
    parser.add_argument('--e-susie', help='Path to eQTL SuSiE results')
    parser.add_argument('--gwas-coloc', nargs='+', help='Path(s) to GWAS co-localization file(s)')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--coloc-cutoff', type=float, default=0.75, help='Co-localization threshold')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Set up logging
    if args.verbose:
        setup_logging(logging.DEBUG)
    else:
        setup_logging()
    
    logger = logging.getLogger(__name__)
    
    try:
        if args.mode == 'qtl':
            logger.info('Running QTL signal groups analysis...')
            
            if not all([args.pair_coloc, args.pc_susie, args.e_susie]):
                raise ValueError("QTL mode requires --pair-coloc, --pc-susie, and --e-susie arguments")
            
            # Load data
            # Handle multiple pair_coloc files (per chromosome)
            if len(args.pair_coloc) == 1:
                pair_coloc = pd.read_csv(args.pair_coloc[0], sep='\t')
            else:
                logger.info(f'Loading {len(args.pair_coloc)} pair_coloc files...')
                pair_coloc_dfs = []
                for file_path in args.pair_coloc:
                    if os.path.exists(file_path):
                        df = pd.read_csv(file_path, sep='\t')
                        pair_coloc_dfs.append(df)
                    else:
                        logger.warning(f'File not found: {file_path}')
                
                if pair_coloc_dfs:
                    pair_coloc = pd.concat(pair_coloc_dfs, ignore_index=True)
                    logger.info(f'Combined {len(pair_coloc_dfs)} files into {len(pair_coloc)} rows')
                else:
                    raise ValueError("No valid pair_coloc files found")
            
            pc_susie = pd.read_csv(args.pc_susie, sep='\t')
            e_susie = pd.read_csv(args.e_susie, sep='\t')
            
            # Run QTL signal group analysis
            signal_groups = group_qtl_signals(pair_coloc, pc_susie, e_susie, args.coloc_cutoff)
            
            logger.info(f'Saving results to {args.output}')
            signal_groups.to_csv(args.output, sep='\t', index=False)
            
            logger.info('QTL signal group analysis completed successfully')
            
        elif args.mode == 'gwas':
            logger.info('Running GWAS signal groups analysis...')
            
            if not all([args.gwas_coloc, args.pair_coloc]):
                raise ValueError("GWAS mode requires --gwas-coloc and --pair-coloc arguments")
            
            # Run GWAS signal groups analysis
            group_gwas_signals(args.gwas_coloc, args.pair_coloc, args.output, args.coloc_cutoff)
        
    except Exception as e:
        logger.error(f'Error in signal group analysis: {e}')
        sys.exit(1)


if __name__ == "__main__":
    main() 