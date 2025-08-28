#!/usr/bin/env python3
"""
Annotate Null Gene Clusters Script

This script generates null gene clusters of a specified size from the set of
expressed genes, excludes overlaps with provided real clusters,
adds the tissue identifier, annotates clusters using the same pipeline as
real clusters, and writes the annotated null clusters to output.
"""

import argparse
import logging
import os
import sys
from typing import Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

# Import shared functionality from annotate_clusters
try:
    from annotate_clusters import add_annotations
except ImportError:
    # Fallback when running directly
    sys.path.append(os.path.dirname(__file__))
    from annotate_clusters import add_annotations


def setup_logging(level: int = logging.INFO) -> None:
    """Configure root logger."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)],
    )


def load_expressed_genes(expression_path: str) -> pd.DataFrame:
    """Load expressed genes bed/tsv with required columns chr, start, end, gene_id."""
    genes = pd.read_table(expression_path).rename(columns={'#chr':'chr'})
    return genes[['chr', 'start', 'end', 'gene_id']]


def get_null_clusters(expressed_genes: pd.DataFrame, num_genes: int, cluster_df: Optional[pd.DataFrame]) -> pd.DataFrame:
    """
    Generate null clusters by sliding window across chromosomes and optionally
    exclude windows overlapping real clusters (two or more genes in common).
    """
    logger = logging.getLogger(__name__)
    logger.info(f'Generating null clusters of size {num_genes} genes')

    expressed_genes = expressed_genes.sort_values(['chr', 'start', 'end'])

    null_cluster_dfs = []
    for chr_id in range(1, 23):
        logger.debug(f'Processing chromosome {chr_id}')
        chr_subset = expressed_genes[expressed_genes['chr'] == f'chr{chr_id}']

        if len(chr_subset) < num_genes:
            logger.debug(f'chr{chr_id} has fewer than {num_genes} genes, skipping')
            continue

        cluster_ids = []
        for i in range(len(chr_subset) - num_genes + 1):
            cluster_transcripts = chr_subset['gene_id'].iloc[i:i+num_genes].astype(str).tolist()
            cluster_ids.append('_'.join(sorted(cluster_transcripts)))

        chr_null_df = pd.DataFrame({'cluster_id': cluster_ids, 'chr': f"chr{chr_id}"})

        if cluster_df is not None and 'cluster_id' in cluster_df.columns:
            try:
                real_transcripts = np.concatenate(cluster_df['cluster_id'].astype(str).str.split('_').values)
                chr_null_df['transcript_list'] = chr_null_df['cluster_id'].str.split('_')

                exploded = chr_null_df.explode('transcript_list')
                exploded['in_cluster'] = exploded['transcript_list'].isin(real_transcripts).astype(int)
                overlap_counts = exploded.groupby('cluster_id').agg({'in_cluster': 'sum'})
                excluded = overlap_counts[overlap_counts['in_cluster'] > 1].index
                chr_null_df = chr_null_df[~chr_null_df['cluster_id'].isin(excluded)]
                logger.debug(f'Excluded {len(excluded)} null clusters overlapping real clusters on chr{chr_id}')
            except Exception:
                pass

        null_cluster_dfs.append(chr_null_df[['cluster_id', 'chr']])

    null_df = pd.concat(null_cluster_dfs, ignore_index=True)
    null_df['num_genes'] = num_genes
    logger.info(f'Generated {len(null_df)} null clusters')
    return null_df[['cluster_id', 'chr', 'num_genes']]


def annotate_null_clusters_from_config(config: dict, tissue_id: str, num_genes: int) -> pd.DataFrame:
    """
    Construct paths directly from config using working_dir and specific path keys,
    generate and annotate null clusters, and return the result (no writing).
    """
    logger = logging.getLogger(__name__)
    wd = config['working_dir']

    clusters = f"{wd}/{config['clusters_dir']}/{tissue_id}.clusters.txt"
    expression = f"{wd}/{config['expression_dir']}/{tissue_id}.v8.normalized_expression.bed"
    covariates = f"{wd}/{config['covariates_dir']}/{tissue_id}.v8.covariates.txt"
    gencode = f"{wd}/{config['gencode_path']}"
    abc_path = f"{wd}/{config['abc_path']}"
    ctcf_match = f"{wd}/{config['ctcf_match_path']}"
    ctcf_dir = f"{wd}/{config['ctcf_dir']}"
    paralog = f"{wd}/{config['paralog_path']}"
    go = f"{wd}/{config['go_path']}"
    cross_map = f"{wd}/{config['cross_map_path']}"
    tad = f"{wd}/{config['tad_path']}"

    annotated = annotate_null_clusters(
        num_genes=num_genes,
        tissue_id=tissue_id,
        gencode_path=gencode,
        paralog_path=paralog,
        abc_path=abc_path,
        go_path=go,
        tad_path=tad,
        ctcf_match_path=ctcf_match,
        ctcf_dir=ctcf_dir,
        cross_map_path=cross_map,
        covariates_path=covariates,
        expression_path=expression,
        real_clusters_path=clusters,
    )

    return annotated


def annotate_null_clusters(
    num_genes: int,
    tissue_id: str,
    gencode_path: str,
    paralog_path: str,
    abc_path: str,
    go_path: str,
    tad_path: str,
    ctcf_match_path: str,
    ctcf_dir: str,
    cross_map_path: str,
    covariates_path: str,
    expression_path: str,
    real_clusters_path: str,
) -> pd.DataFrame:
    """Generate and annotate null clusters, excluding provided real clusters."""
    logger = logging.getLogger(__name__)
    expressed_genes = load_expressed_genes(expression_path)

    real_clusters_df = pd.read_table(real_clusters_path)
    null_clusters = get_null_clusters(expressed_genes, num_genes, real_clusters_df)
    null_clusters['tissue_id'] = tissue_id

    annotated_nulls = add_annotations(
        null_clusters,
        tissue_id,
        gencode_path,
        paralog_path,
        abc_path,
        go_path,
        tad_path,
        ctcf_match_path,
        ctcf_dir,
        cross_map_path,
        covariates_path,
        expression_path,
    )
    return annotated_nulls


def main() -> None:
    """Parse CLI args, generate and annotate null clusters, and write output."""
    setup_logging()
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description='Generate and annotate null gene clusters',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Keep ordering aligned with add_annotations, but replace clusters with cluster-size
    parser.add_argument('-t', '--tissue-id', required=True, help='GTEx tissue identifier')
    parser.add_argument('--cluster-size', type=int, required=True, help='Number of genes per null cluster')
    parser.add_argument('--gencode', required=True, help='GENCODE annotation file')
    parser.add_argument('--paralog', required=True, help='Paralog data file')
    parser.add_argument('--abc', required=True, help='ABC predictions file')
    parser.add_argument('--go', required=True, help='GO term file')
    parser.add_argument('--tad', required=True, help='TAD boundary file')
    parser.add_argument('--ctcf-match', required=True, help='CTCF-GTEx tissue matching file')
    parser.add_argument('--ctcf-dir', required=True, help='Directory containing CTCF binding site files')
    parser.add_argument('--cross-map', required=True, help='Cross-mappability file')
    parser.add_argument('-co', '--covariates', required=True, help='Covariates file')
    parser.add_argument('-e', '--expression', required=True, help='Normalized expression BED/TSV file')
    parser.add_argument('-c', '--clusters', required=True, help='Real cluster CSV/TSV to exclude overlaps')
    parser.add_argument('-o', '--output', required=True, help='Output file for annotated null clusters')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose logging')

    args = parser.parse_args()

    if args.verbose:
        setup_logging(logging.DEBUG)
        logger = logging.getLogger(__name__)

    try:
        # Validate inputs
        required_paths = {
            'gencode': args.gencode,
            'paralog': args.paralog,
            'abc': args.abc,
            'go': args.go,
            'tad': args.tad,
            'ctcf_match': args.ctcf_match,
            'cross_map': args.cross_map,
            'covariates': args.covariates,
            'expression': args.expression,
            'clusters': args.clusters,
        }
        missing = [k for k, p in required_paths.items() if not p or not os.path.exists(p)]
        if args.ctcf_dir and not os.path.isdir(args.ctcf_dir):
            missing.append('ctcf_dir')
        if missing:
            raise FileNotFoundError(f"Missing or non-existent required inputs: {', '.join(sorted(set(missing)))}")

        logger.info(f"Generating null clusters of size {args.cluster_size} for tissue {args.tissue_id}")

        annotated_nulls = annotate_null_clusters(
            num_genes=args.cluster_size,
            tissue_id=args.tissue_id,
            gencode_path=args.gencode,
            paralog_path=args.paralog,
            abc_path=args.abc,
            go_path=args.go,
            tad_path=args.tad,
            ctcf_match_path=args.ctcf_match,
            ctcf_dir=args.ctcf_dir,
            cross_map_path=args.cross_map,
            covariates_path=args.covariates,
            expression_path=args.expression,
            real_clusters_path=args.clusters,
        )

        logger.info(f"Writing annotated null clusters to {args.output}")
        annotated_nulls.to_csv(args.output, index=False, sep='\t')
        logger.info(f"Successfully wrote {len(annotated_nulls)} annotated null clusters to {args.output}")

    except Exception as e:
        logger.error(f"Failed to annotate null clusters: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()


