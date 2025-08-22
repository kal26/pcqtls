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
from typing import List, Union, Optional
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

# Constants
DEFAULT_COLOC_CUTOFF = 0.75
DEFAULT_COLUMNS = ['PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf']


def setup_logging(level: int = logging.INFO) -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )





def group_qtl_signals(pair_coloc: pd.DataFrame, pc_susie_r: pd.DataFrame, e_susie_r: pd.DataFrame, 
                     coloc_cutoff: float = DEFAULT_COLOC_CUTOFF, get_variants: bool = False) -> pd.DataFrame:
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
    # Check if we have any data to work with
    if pair_coloc.empty and pc_susie_r.empty and e_susie_r.empty:
        # Return empty DataFrame with expected columns
        columns = ['signal_id', 'cluster_id', 'num_e_coloc', 'num_pc_coloc', 'lead_var_set']
        if get_variants:
            columns.append('var_set')
        return pd.DataFrame(columns=columns)
    
    # Add lead_variant_id and cs_id_full to SuSiE results (only if not empty)
    def add_lead_variant_id(susie_df: pd.DataFrame) -> pd.DataFrame:
        if susie_df.empty:
            return susie_df
        # Group by phenotype_id and cs_id, then find the variant with highest PIP
        lead_variants = susie_df.loc[susie_df.groupby(['phenotype_id', 'cs_id'])['pip'].idxmax()]
        lead_variants = lead_variants[['phenotype_id', 'cs_id', 'variant_id']].rename(columns={'variant_id': 'lead_variant_id'})
        # Merge back to original dataframe
        susie_df = susie_df.merge(lead_variants, on=['phenotype_id', 'cs_id'], how='left')
        # Add qtl_cs_id column to match pairwise coloc format
        susie_df['qtl_cs_id'] = susie_df['phenotype_id'] + '_cs_' + susie_df['cs_id'].astype(str)
        return susie_df
    
    pc_susie_r = add_lead_variant_id(pc_susie_r)
    e_susie_r = add_lead_variant_id(e_susie_r)
    
    # Filter to coloc'd signals (only if not empty)
    if not pair_coloc.empty:
        pair_coloc = pair_coloc[pair_coloc['PP.H4.abf'] > coloc_cutoff].copy()
    else:
        pair_coloc = pd.DataFrame()

    # Create an undirected graph
    G = nx.Graph()

    # Create qtl_cs_id for pairwise coloc data to match SuSiE files (only if not empty)
    if not pair_coloc.empty:
        pair_coloc = pair_coloc.copy()
        pair_coloc['qtl_cs_id_1'] = pair_coloc['qtl1_id'] + '_cs_' + pair_coloc['idx1'].astype(int).astype(str)
        pair_coloc['qtl_cs_id_2'] = pair_coloc['qtl2_id'] + '_cs_' + pair_coloc['idx2'].astype(int).astype(str)
    
    if not pc_susie_r.empty:
        pc_susie_r = pc_susie_r.copy()
        pc_susie_r['qtl_cs_id'] = pc_susie_r['phenotype_id'] + '_cs_' + pc_susie_r['cs_id'].astype(int).astype(str)
    
    if not e_susie_r.empty:
        e_susie_r = e_susie_r.copy()
        e_susie_r['qtl_cs_id'] = e_susie_r['phenotype_id'] + '_cs_' + e_susie_r['cs_id'].astype(int).astype(str)
    
    # Add edges to the graph from the DataFrame (only if not empty)
    if not pair_coloc.empty:
        for index, row in pair_coloc.iterrows():
            G.add_edge(row['qtl_cs_id_1'], row['qtl_cs_id_2'])

    # Get the connected components of the graph
    connected_components = list(nx.connected_components(G))

    # Generate underlying signal ids
    underlying_signals = ['-'.join(sorted(str(item) for item in component)) for component in connected_components]

    # Add in the unique ones (only if we have SuSiE data)
    if not pc_susie_r.empty or not e_susie_r.empty:
        all_credible_set_ids = set()
        if not pc_susie_r.empty:
            all_credible_set_ids.update(pc_susie_r['qtl_cs_id'])
        if not e_susie_r.empty:
            all_credible_set_ids.update(e_susie_r['qtl_cs_id'])

        # Add standalone sets (those not in any connected component)
        for credible_set_id in all_credible_set_ids:
            if not any(credible_set_id in component for component in connected_components):
                underlying_signals.append(str(credible_set_id))

    underlying_signals = pd.DataFrame({'signal_id': underlying_signals})
    # Extract cluster_id from signal_id - extract from the first part before _cs_
    def extract_cluster_id(signal_id: str) -> Optional[str]:
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
    def get_var_set(row: pd.Series, lead: bool = True) -> List[str]:
        var_ids = []
        signal_parts = row['signal_id'].split('-')
        
        for part in signal_parts:
            # Find matching variants in pc_susie_r using qtl_cs_id (only if not empty)
            if not pc_susie_r.empty:
                pc_matches = pc_susie_r[pc_susie_r['qtl_cs_id'] == part]
                if lead:
                    var_ids.extend(pc_matches['lead_variant_id'].dropna().values)
                else:
                    var_ids.extend(pc_matches['variant_id'].values)
            
            # Find matching variants in e_susie_r using qtl_cs_id (only if not empty)
            if not e_susie_r.empty:
                e_matches = e_susie_r[e_susie_r['qtl_cs_id'] == part]
                if lead:
                    var_ids.extend(e_matches['lead_variant_id'].dropna().values)
                else:
                    var_ids.extend(e_matches['variant_id'].values)
        
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





def _prepare_gwas_coloc_data(gwas_coloc: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare GWAS co-localization data by adding necessary columns and annotations.
    
    This helper function adds cluster_id, gwas_cs_id, qtl_cs_id, and other
    required columns to GWAS co-localization data.
    
    Args:
        gwas_coloc: DataFrame with GWAS co-localization results
    
    Returns:
        DataFrame with added columns and annotations
    
    Raises:
        ValueError: If required columns are missing from the input DataFrame
    """
    if gwas_coloc.empty:
        return gwas_coloc
    
    # Validate required columns
    required_cols = ['gwas_id', 'qtl_id', 'idx1', 'idx2']
    missing_cols = [col for col in required_cols if col not in gwas_coloc.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}. Available columns: {list(gwas_coloc.columns)}")
    
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
    gwas_coloc[DEFAULT_COLUMNS] = gwas_coloc[DEFAULT_COLUMNS].fillna(0)

    # Add cs ids
    gwas_coloc['gwas_cs_id'] = 'gwas_' + gwas_coloc['gwas_cs_id']
    gwas_coloc['qtl_cs_id'] = 'qtl_' + gwas_coloc['qtl_cs_id']
    
    return gwas_coloc


def _load_data_files(file_paths: Union[str, List[str]], logger: logging.Logger, file_type: str = "data") -> pd.DataFrame:
    """
    Load data from file(s) with proper error handling.
    
    This helper function handles loading single or multiple files and concatenating
    them into a single DataFrame.
    
    Args:
        file_paths: String or list of file paths
        logger: Logger instance for logging messages
        file_type: Type of files being loaded (for logging purposes)
    
    Returns:
        DataFrame with loaded data
    """
    # Handle single or multiple files
    if isinstance(file_paths, str):
        file_paths = [file_paths]
        
    if len(file_paths) == 1:
        try:
            data = pd.read_csv(file_paths[0], sep='\t')
            if data.empty:
                logger.warning(f'File is empty: {file_paths[0]}')
                data = pd.DataFrame()
        except pd.errors.EmptyDataError:
            logger.warning(f'File is empty or has no valid data: {file_paths[0]}')
            data = pd.DataFrame()
    else:
        logger.info(f'Loading {len(file_paths)} {file_type} files...')
        data_dfs = []
        empty_files = []
        for file_path in file_paths:
            if os.path.exists(file_path):
                try:
                    df = pd.read_csv(file_path, sep='\t')
                    if not df.empty:
                        data_dfs.append(df)
                    else:
                        empty_files.append(file_path)
                except pd.errors.EmptyDataError:
                    empty_files.append(file_path)
                    logger.debug(f'File is empty or has no valid data: {file_path}')
            else:
                logger.warning(f'File not found: {file_path}')
        
        if empty_files:
            logger.info(f'Skipped {len(empty_files)} empty files')
        
        if data_dfs:
            data = pd.concat(data_dfs, ignore_index=True)
            logger.info(f'Combined {len(data_dfs)} {file_type} files into {len(data)} rows')
        else:
            logger.warning(f"No valid {file_type} files found with data")
            data = pd.DataFrame()
    
    return data





def load_gwas_coloc(config: Dict[str, str]) -> pd.DataFrame:
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
    def get_files(directory: str) -> List[str]:
        file_list = []
        for root, directories, files in os.walk(directory):
            if not 'temp' in root:
                for file in files:
                    if "susie_True" in file:
                        file_list.append(os.path.join(root, file))
        return file_list
    
    coloc_file_list = get_files(f'{config["working_dir"]}/{config["coloc_output_dir"]}gwas')

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
    
    # Use the helper function to prepare the data
    return _prepare_gwas_coloc_data(gwas_coloc)


def group_gwas_signals(gwas_coloc: pd.DataFrame, pair_coloc: pd.DataFrame, coloc_cutoff: float = .75) -> pd.DataFrame:
    """
    Group GWAS signals from co-localization results.
    
    This function creates signal groups that include both QTL-QTL and QTL-GWAS
    co-localizations using network analysis to identify connected components.
    
    Args:
        gwas_coloc: DataFrame with GWAS co-localization results
        pair_coloc: DataFrame with pairwise QTL co-localization results
        coloc_cutoff: PP.H4 threshold for co-localization (default: 0.75)
    
    Returns:
        DataFrame with GWAS signal groups
    """
    # Check if we have any data to work with
    if gwas_coloc.empty and pair_coloc.empty:
        # Return empty DataFrame with expected columns
        return pd.DataFrame(columns=['signal_id', 'cluster_id', 'num_qtl_coloc', 'num_gwas_coloc', 'num_e_coloc', 'num_pc_coloc'])
    
    # Prepare GWAS co-localization data
    if not gwas_coloc.empty:
        gwas_coloc = _prepare_gwas_coloc_data(gwas_coloc)
        # Filter to significant co-localizations
        gwas_coloc_hits = gwas_coloc[gwas_coloc['PP.H4.abf'] > coloc_cutoff]
    else:
        gwas_coloc_hits = pd.DataFrame()
    
    # Filter pairwise co-localization data
    if not pair_coloc.empty:
        pair_coloc_hits = pair_coloc[pair_coloc['PP.H4.abf'] > coloc_cutoff].copy()
        # Create CS IDs for pairwise coloc data
        pair_coloc_hits['cs_id_1'] = pair_coloc_hits['qtl1_id'] + '_cs_' + pair_coloc_hits['idx1'].astype(int).astype(str)
        pair_coloc_hits['cs_id_2'] = pair_coloc_hits['qtl2_id'] + '_cs_' + pair_coloc_hits['idx2'].astype(int).astype(str)
    else:
        pair_coloc_hits = pd.DataFrame()

    # Create an undirected graph
    G = nx.Graph()
    
    # Add an edge for each gwas-qtl coloc
    if not gwas_coloc_hits.empty:
        for index, row in gwas_coloc_hits.iterrows():
            G.add_edge(row['gwas_cs_id'], row['qtl_cs_id'])

    # Add an edge for each qtl-qtl coloc
    if not pair_coloc_hits.empty:
        for index, row in pair_coloc_hits.iterrows():
            G.add_edge(row['cs_id_1'], row['cs_id_2'])

    # Get the connected components of the graph
    connected_components = list(nx.connected_components(G))
    
    # Generate underlying signal ids
    underlying_signals = ['-'.join(sorted(str(item) for item in component)) for component in connected_components]

    # Make a df
    underlying_signals = pd.DataFrame({'signal_id': underlying_signals})
    # Extract cluster_id from signal_id - extract from the first part before _cs_
    def extract_cluster_id(signal_id: str) -> Optional[str]:
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





def main() -> None:
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
            
            # Load pairwise co-localization data
            pair_coloc = _load_data_files(args.pair_coloc, logger, "pair_coloc")
            
            # Load SuSiE files
            try:
                pc_susie = pd.read_csv(args.pc_susie, sep='\t')
                if pc_susie.empty:
                    logger.warning(f'File is empty: {args.pc_susie}')
                    pc_susie = pd.DataFrame()
            except pd.errors.EmptyDataError:
                logger.warning(f'File is empty or has no valid data: {args.pc_susie}')
                pc_susie = pd.DataFrame()
                
            try:
                e_susie = pd.read_csv(args.e_susie, sep='\t')
                if e_susie.empty:
                    logger.warning(f'File is empty: {args.e_susie}')
                    e_susie = pd.DataFrame()
            except pd.errors.EmptyDataError:
                logger.warning(f'File is empty or has no valid data: {args.e_susie}')
                e_susie = pd.DataFrame()
            
            # Check if we have any data to work with
            if pair_coloc.empty and pc_susie.empty and e_susie.empty:
                logger.warning("No data found in any input files")
                # Create empty result with expected columns
                signal_groups = pd.DataFrame(columns=['signal_id', 'cluster_id', 'num_e_coloc', 'num_pc_coloc', 'lead_var_set'])
                logger.info(f'Saving empty results to {args.output}')
                signal_groups.to_csv(args.output, sep='\t', index=False)
                logger.info('QTL signal group analysis completed (no data found)')
            else:
                # Run QTL signal group analysis
                signal_groups = group_qtl_signals(pair_coloc, pc_susie, e_susie, args.coloc_cutoff)
                
                logger.info(f'Saving results to {args.output}')
                signal_groups.to_csv(args.output, sep='\t', index=False)
                
                logger.info('QTL signal group analysis completed successfully')
            
        elif args.mode == 'gwas':
            logger.info('Running GWAS signal groups analysis...')
            
            if not all([args.gwas_coloc, args.pair_coloc]):
                raise ValueError("GWAS mode requires --gwas-coloc and --pair-coloc arguments")
            
            # Load GWAS co-localization data
            gwas_coloc = _load_data_files(args.gwas_coloc, logger, "gwas_coloc")
            
            # Load pairwise co-localization data
            pair_coloc = _load_data_files(args.pair_coloc, logger, "pair_coloc")
            
            # Check if we have any data to work with
            if gwas_coloc.empty and pair_coloc.empty:
                logger.warning("No data found in either GWAS or pairwise co-localization files")
                # Create empty result with expected columns
                gwas_signal_groups = pd.DataFrame(columns=['signal_id', 'cluster_id', 'num_qtl_coloc', 'num_gwas_coloc', 'num_e_coloc', 'num_pc_coloc'])
                logger.info(f'Saving empty results to {args.output}')
                gwas_signal_groups.to_csv(args.output, sep='\t', index=False)
                logger.info('GWAS signal groups analysis completed (no data found)')
            else:
                # Run GWAS signal groups analysis
                gwas_signal_groups = group_gwas_signals(gwas_coloc, pair_coloc, args.coloc_cutoff)
                
                logger.info(f'Saving results to {args.output}')
                gwas_signal_groups.to_csv(args.output, sep='\t', index=False)
                
                logger.info('GWAS signal groups analysis completed successfully')
        
    except Exception as e:
        logger.error(f'Error in signal group analysis: {e}')
        sys.exit(1)


if __name__ == "__main__":
    main() 