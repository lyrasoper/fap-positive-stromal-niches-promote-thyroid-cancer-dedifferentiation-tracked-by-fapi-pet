#!/usr/bin/env python3
# Cancer Research submission - figure code release
# Builds: Figure 1H (preranked Hallmark GSEA).

import argparse
from pathlib import Path

import pandas as pd
import gseapy as gp


def main() -> None:
    parser = argparse.ArgumentParser(description="Run preranked Hallmark GSEA for Fig1H.")
    parser.add_argument("--rnk", required=True, help="CSV with gene and score columns.")
    parser.add_argument("--gmt", required=True, help="Hallmark GMT file.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rnk = pd.read_csv(args.rnk)
    if not {"gene", "score"}.issubset(rnk.columns):
      raise ValueError("Ranking file must contain 'gene' and 'score' columns.")
    rnk = rnk.dropna(subset=["gene", "score"]).drop_duplicates(subset=["gene"]).sort_values("score", ascending=False)

    pre_res = gp.prerank(
        rnk=rnk[["gene", "score"]],
        gene_sets=args.gmt,
        min_size=5,
        max_size=500,
        permutation_num=1000,
        seed=1,
        outdir=None,
        no_plot=True,
        threads=4,
    )

    res = pre_res.res2d.copy()
    res.to_csv(outdir / "Fig1H_gsea_full_results.csv", index=False)

    positive = res.copy()
    positive = positive[positive["NES"].notna() & positive["FDR q-val"].notna()]
    positive = positive[positive["NES"] > 0].sort_values(["FDR q-val", "NES"], ascending=[True, False])
    if positive.empty:
        positive = res[res["NES"].notna() & (res["NES"] > 0)].sort_values("NES", ascending=False)
    top = positive.head(10).copy()
    top["Pathway"] = (
        top["Term"]
        .str.replace("^HALLMARK_", "", regex=True)
        .str.replace("_", " ", regex=False)
        .str.title()
    )
    top = top.loc[:, ["Pathway", "NES", "FDR q-val", "NOM p-val"]]
    top = top.rename(columns={"FDR q-val": "FDR_qval", "NOM p-val": "nominal_pval"})
    top.to_csv(outdir / "Fig1H_top_positive_hallmark_pathways.csv", index=False)


if __name__ == "__main__":
    main()
