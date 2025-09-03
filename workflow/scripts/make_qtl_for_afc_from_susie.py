#!/usr/bin/env python3
"""Convert SuSiE credible set output to qtl_for_afc format.

Reads a SuSiE results table containing at least 'phenotype_id' and 'variant_id',
outputs a TSV with columns: pid, sid, sid_chr, sid_pos.

Usage:
  python make_qtl_for_afc_from_susie.py --susie <susie.txt> --output <out.tsv> [--min-pip 0.0]
"""

import argparse
import pandas as pd


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser()
	parser.add_argument("--susie", required=True, help="Path to SuSiE results (TSV)")
	parser.add_argument("--output", required=True, help="Path to output qtl_for_afc TSV")
	parser.add_argument("--qtl_type", required=True, help="Type of QTL to make (eqtl or pcqtl)")
	parser.add_argument("--min-pip", type=float, default=0.0, help="Minimum PIP to include (default 0.0)")

	return parser.parse_args()


def main() -> None:
	args = parse_args()

	df = pd.read_table(args.susie)
	if df.shape[0] == 0:
		out = pd.DataFrame(columns=["pid", "sid", "sid_chr", "sid_pos"])  # empty
	else:
		if "pip" in df.columns and args.min_pip > 0:
			df = df[df["pip"] >= args.min_pip]
		if args.qtl_type == "pcqtl":
			df['cluster_id'] = df['phenotype_id'].str.split("_pc").str[0]
		elif args.qtl_type == "eqtl":
			df['cluster_id'] = df['phenotype_id'].str.split("_e").str[0]
		df['pid'] = df['cluster_id'].str.split("_")
		df_explode = df.explode('pid')
		sub = df_explode[["pid", "variant_id"]].drop_duplicates().rename(columns={"variant_id": "sid"})
		parts = sub["sid"].astype(str).str.split("_", expand=True)
		if parts.shape[1] >= 2:
			sub["sid_chr"] = parts[0].str.replace("chr", "")
			sub["sid_pos"] = pd.to_numeric(parts[1], errors="coerce")
		out = sub

	out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
	main() 