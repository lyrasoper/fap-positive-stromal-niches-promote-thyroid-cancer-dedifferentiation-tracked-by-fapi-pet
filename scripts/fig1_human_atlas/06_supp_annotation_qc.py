#!/usr/bin/env python3
# Cancer Research submission - figure code release
# Builds: Supplementary figures for Figure 1 (annotation/QC support).

"""Redesign Supplementary Fig. 1 as annotation/QC support for main Fig. 1.

The goal is to reduce overlap with the biological claims in main Fig. 1. This
preview emphasizes atlas integration, lineage annotation and CAF annotation
robustness, while keeping sample-level fibroblast/FAP+ CAF summaries secondary.
"""

from __future__ import annotations
import os

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


ROOT = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>"))
BASE_SCRIPT = ROOT / "scripts/fig1_human_atlas/05_supp_atlas_preview.py"
SRC = ROOT / "outputs/supplementary_figures/source_data"
META = ROOT / "outputs/fig1_scRNA/source_data/fig1_scrna_full_metadata.csv"
OUT = ROOT / "outputs/supplementary_figures/nature_preview_20260517"
OUT.mkdir(parents=True, exist_ok=True)

STEM = "Supplementary_Fig_1_annotation_QC_support_layout_preview"
PNG_OUT = OUT / f"{STEM}.png"
PDF_OUT = OUT / f"{STEM}.pdf"
SVG_OUT = OUT / f"{STEM}.svg"

spec = importlib.util.spec_from_file_location("suppfig1_base", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(base)


TISSUE_ORDER = base.TISSUE_ORDER
CELL_ORDER = base.CELL_ORDER
CAF_ORDER = base.CAF_ORDER
STUDY_ORDER = base.STUDY_ORDER
TISSUE_COLORS = base.TISSUE_COLORS
CELL_COLORS = base.CELL_COLORS
CAF_COLORS = base.CAF_COLORS
STUDY_COLORS = base.STUDY_COLORS


def configure_style() -> None:
    sns.set_theme(style="white", context="paper", font_scale=0.72, rc={"axes.linewidth": 0.45})
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 5.8,
            "axes.labelsize": 5.7,
            "axes.titlesize": 6.4,
            "xtick.labelsize": 5.0,
            "ytick.labelsize": 5.0,
            "legend.fontsize": 4.9,
            "legend.title_fontsize": 5.0,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def clean_axis(ax: plt.Axes, *, grid: bool = False) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", length=2.2, width=0.45, pad=1.3)
    if grid:
        ax.grid(axis="y", color="#DADADA", lw=0.35, alpha=0.75)
        ax.set_axisbelow(True)


def panel_label(ax: plt.Axes, label: str, *, x: float = -0.13, y: float = 1.13) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.0,
        fontweight="bold",
        color="black",
    )


def plot_sample_structure(ax: plt.Axes) -> None:
    meta = pd.read_csv(META, usecols=["Sample_id", "TissueType", "Data.ID", "Celltype_short"])
    sample_df = meta.drop_duplicates("Sample_id")
    counts = (
        sample_df.groupby(["Data.ID", "TissueType"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    pivot = (
        counts.pivot(index="Data.ID", columns="TissueType", values="n")
        .reindex(STUDY_ORDER)
        .reindex(columns=TISSUE_ORDER)
        .fillna(0)
    )
    x = np.arange(len(STUDY_ORDER))
    bottom = np.zeros(len(STUDY_ORDER))
    for tissue in TISSUE_ORDER:
        vals = pivot[tissue].to_numpy()
        ax.bar(x, vals, bottom=bottom, color=TISSUE_COLORS[tissue], width=0.66, edgecolor="white", lw=0.25)
        bottom += vals
    for xpos, total in zip(x, bottom):
        ax.text(xpos, total + 0.35, f"{int(total)}", ha="center", va="bottom", fontsize=4.8)
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("GSE", "GSE\n") for s in STUDY_ORDER])
    ax.set_ylabel("Samples")
    ax.set_title("Atlas sample structure", loc="left", fontweight="bold", pad=2)
    clean_axis(ax, grid=True)
    handles = [mpl.patches.Patch(facecolor=TISSUE_COLORS[t], label=t) for t in TISSUE_ORDER]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.22), ncol=3, frameon=False,
              handlelength=0.9, handletextpad=0.25, columnspacing=0.55)


def plot_cohort_mixing(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Supplementary_Fig1H_celltype_cohort_mixing.csv")
    df = df[df["Celltype_short"].isin(CELL_ORDER)].copy()
    pivot = (
        df.pivot(index="Celltype_short", columns="Data.ID", values="prop")
        .reindex(CELL_ORDER)
        .reindex(columns=STUDY_ORDER)
        .fillna(0)
    )
    x = np.arange(len(CELL_ORDER))
    bottom = np.zeros(len(CELL_ORDER))
    for study in STUDY_ORDER:
        vals = pivot[study].to_numpy()
        ax.bar(x, vals, bottom=bottom, color=STUDY_COLORS[study], width=0.72, edgecolor="white", lw=0.20)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(CELL_ORDER, rotation=25, ha="right")
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1.0, decimals=0))
    ax.set_ylabel("Cells")
    ax.set_title("Cohort mixing by lineage", loc="left", fontweight="bold", pad=2)
    clean_axis(ax, grid=True)
    handles = [mpl.patches.Patch(facecolor=STUDY_COLORS[s], label=s) for s in STUDY_ORDER]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2, frameon=False,
              handlelength=0.9, handletextpad=0.25, columnspacing=0.55)


def plot_fraction_context(ax: plt.Axes, csv_name: str, y_col: str, title: str, ylabel: str) -> None:
    df = pd.read_csv(SRC / csv_name)
    df = df[df["TissueType"].isin(TISSUE_ORDER)].copy()
    df["TissueType"] = pd.Categorical(df["TissueType"], categories=TISSUE_ORDER, ordered=True)
    sns.boxplot(
        data=df,
        x="TissueType",
        y=y_col,
        order=TISSUE_ORDER,
        palette=TISSUE_COLORS,
        ax=ax,
        width=0.52,
        fliersize=0,
        linewidth=0.5,
        saturation=0.90,
        boxprops={"alpha": 0.82},
        medianprops={"color": "black", "linewidth": 0.6},
        whiskerprops={"linewidth": 0.5},
        capprops={"linewidth": 0.5},
    )
    sns.stripplot(
        data=df,
        x="TissueType",
        y=y_col,
        order=TISSUE_ORDER,
        ax=ax,
        color="black",
        size=1.9,
        jitter=0.17,
        alpha=0.52,
        linewidth=0,
    )
    ax.set_title(title, loc="left", fontweight="bold", pad=2)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, min(1.0, df[y_col].max() * 1.16 + 0.03))
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1.0, decimals=0))
    clean_axis(ax, grid=True)


def plot_fap_lineage_dot(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Supplementary_Fig1I_FAP_celltype_expression.csv")
    df = df[df["Celltype_short"].isin(CELL_ORDER)].copy()
    summary = (
        df.groupby("Celltype_short", observed=True)
        .agg(
            pct_expr=("FAP", lambda s: float((s > 0).mean() * 100)),
            mean_expr=("FAP", "mean"),
            n=("FAP", "size"),
        )
        .reindex(CELL_ORDER)
        .reset_index()
    )
    x = np.arange(len(summary))
    ax.vlines(x, 0, summary["pct_expr"], color="#B8B8B8", lw=0.65, zorder=1)
    sizes = 28 + summary["mean_expr"].fillna(0).clip(0, 0.40) / 0.40 * 210
    ax.scatter(
        x,
        summary["pct_expr"],
        s=sizes,
        c=[CELL_COLORS[c] for c in summary["Celltype_short"]],
        edgecolor="#333333",
        linewidth=0.35,
        zorder=3,
    )
    for xpos, pct in zip(x, summary["pct_expr"]):
        ax.text(xpos, pct + 1.2, f"{pct:.1f}%", ha="center", va="bottom", fontsize=4.7)
    ax.set_xticks(x)
    ax.set_xticklabels(summary["Celltype_short"], rotation=25, ha="right")
    ax.set_ylim(0, max(12, summary["pct_expr"].max() * 1.25))
    ax.set_ylabel("FAP+ cells")
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=100, decimals=0))
    ax.set_title("FAP lineage restriction", loc="left", fontweight="bold", pad=2)
    clean_axis(ax, grid=True)
    handles = [
        plt.scatter([], [], s=s, color="#777777", edgecolor="#333333", linewidth=0.35)
        for s in [45, 130, 230]
    ]
    ax.legend(
        handles,
        ["low", "mid", "high"],
        title="Mean expr.",
        loc="upper right",
        frameon=False,
        borderpad=0.1,
        labelspacing=0.25,
        handletextpad=0.35,
    )


def plot_resolution_scan_compact(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Supplementary_Fig1J_caf_resolution_assignments.csv")
    res_cols = [c for c in df.columns if c.startswith("RNA_snn_res.")]
    long = df.melt(id_vars="cell", value_vars=res_cols, var_name="resolution", value_name="cluster")
    long["resolution"] = long["resolution"].str.replace("RNA_snn_res.", "", regex=False)
    counts = (
        long.groupby(["resolution", "cluster"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    counts["prop"] = counts["n"] / counts.groupby("resolution")["n"].transform("sum")
    resolutions = sorted(counts["resolution"].unique(), key=lambda x: float(x))
    clusters = sorted(counts["cluster"].unique(), key=lambda x: int(x) if str(x).isdigit() else str(x))
    pivot = (
        counts.pivot(index="resolution", columns="cluster", values="prop")
        .reindex(resolutions)
        .reindex(columns=clusters)
        .fillna(0)
    )
    colors = sns.color_palette("tab20", n_colors=len(clusters))
    x = np.arange(len(resolutions))
    bottom = np.zeros(len(resolutions))
    for color, cluster in zip(colors, clusters):
        vals = pivot[cluster].to_numpy()
        ax.bar(x, vals, bottom=bottom, color=color, edgecolor="white", linewidth=0.12, width=0.72)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(resolutions)
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1.0, decimals=0))
    ax.set_title("CAF reclustering resolution scan", loc="left", fontweight="bold", pad=2)
    ax.set_xlabel("Seurat resolution")
    ax.set_ylabel("CAF cells")
    clean_axis(ax, grid=True)


def plot_marker_effects(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Supplementary_Fig1L_top_caf_markers.csv")
    df = (
        df[df["cluster"].isin(CAF_ORDER)]
        .sort_values(["cluster", "avg_log2FC"], ascending=[True, False])
        .groupby("cluster", observed=True)
        .head(4)
        .copy()
    )
    rows = []
    y = 0
    for cluster in CAF_ORDER:
        sub = df[df["cluster"] == cluster].sort_values("avg_log2FC", ascending=True)
        for _, row in sub.iterrows():
            rows.append({**row.to_dict(), "y": y})
            y += 1
        y += 0.65
    plot_df = pd.DataFrame(rows)
    sizes = 25 + (plot_df["pct.1"].clip(0, 1) * 145)
    for cluster in CAF_ORDER:
        sub = plot_df[plot_df["cluster"] == cluster]
        ax.scatter(
            sub["avg_log2FC"],
            sub["y"],
            s=sizes.loc[sub.index],
            color=CAF_COLORS[cluster],
            edgecolor="#333333",
            linewidth=0.25,
            alpha=0.92,
            label=cluster,
        )
    ax.set_yticks(plot_df["y"])
    ax.set_yticklabels(plot_df["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("Average log2 fold change")
    ax.set_ylabel("")
    ax.set_title("Top CAF-subtype marker effect sizes", loc="left", fontweight="bold", pad=2)
    clean_axis(ax, grid=True)
    ax.grid(axis="x", color="#E6E6E6", lw=0.35, alpha=0.7)
    ax.set_xlim(0, max(10.2, plot_df["avg_log2FC"].max() + 0.4))

    for cluster in CAF_ORDER:
        sub = plot_df[plot_df["cluster"] == cluster]
        if sub.empty:
            continue
        ymid = (sub["y"].min() + sub["y"].max()) / 2
        ax.text(
            ax.get_xlim()[1] + 0.12,
            ymid,
            cluster,
            ha="left",
            va="center",
            fontsize=5.2,
            color=CAF_COLORS[cluster],
            fontweight="bold",
            clip_on=False,
        )
    handles = [
        plt.scatter([], [], s=s, color="#777777", edgecolor="#333333", linewidth=0.25)
        for s in [25 + 0.15 * 145, 25 + 0.35 * 145, 25 + 0.60 * 145]
    ]
    ax.legend(
        handles,
        ["15%", "35%", "60%"],
        title="pct.1",
        loc="lower right",
        frameon=False,
        borderpad=0.1,
        labelspacing=0.25,
        handletextpad=0.45,
    )


def build_figure() -> None:
    configure_style()

    fig = plt.figure(figsize=(12.0, 15.8), dpi=300)
    gs = GridSpec(
        nrows=6,
        ncols=6,
        figure=fig,
        height_ratios=[0.88, 1.05, 1.05, 0.92, 0.98, 1.35],
        hspace=0.92,
        wspace=0.62,
    )

    ax_a = fig.add_subplot(gs[0, 0:2])
    ax_b = fig.add_subplot(gs[0, 2:4])
    ax_c = fig.add_subplot(gs[0, 4:6])
    ax_d = fig.add_subplot(gs[1, 0:6])
    ax_e = fig.add_subplot(gs[2, 0:2])
    ax_f = fig.add_subplot(gs[2, 2:6])
    ax_g = fig.add_subplot(gs[3, 0:2])
    ax_h = fig.add_subplot(gs[3, 2:4])
    ax_i = fig.add_subplot(gs[3, 4:6])
    ax_j = fig.add_subplot(gs[4, 0:6])
    ax_k = fig.add_subplot(gs[5, 0:6])

    plot_sample_structure(ax_a)
    panel_label(ax_a, "a")

    base.plot_panel_f_qc(ax_b)
    ax_b.set_title("Sample-level QC profile", loc="left", fontweight="bold", pad=2)
    panel_label(ax_b, "b")

    plot_cohort_mixing(ax_c)
    panel_label(ax_c, "c")

    lineage_dot = pd.read_csv(SRC / "Supplementary_Fig1C_major_lineage_dotplot_data.csv")
    lineage_features = [
        "TG",
        "TPO",
        "PAX8",
        "FAP",
        "COL1A1",
        "PECAM1",
        "KDR",
        "LST1",
        "C1QC",
        "CD3D",
        "NKG7",
        "MS4A1",
        "CD79A",
    ]
    base.dotplot(
        ax_d,
        lineage_dot,
        "features.plot",
        "id",
        lineage_features,
        CELL_ORDER,
        "Canonical lineage markers",
        "RdBu_r",
        "Scaled exp.",
    )
    panel_label(ax_d, "d", x=-0.055, y=1.14)

    plot_fap_lineage_dot(ax_e)
    panel_label(ax_e, "e")

    caf_dot = pd.read_csv(SRC / "Supplementary_Fig1F_caf_dotplot_data.csv")
    caf_features = ["FAP", "CXCL12", "IL6", "COL1A1", "POSTN", "PLVAP", "RAMP2", "RGS5", "ACTA2", "SFRP5", "OGN"]
    base.dotplot(
        ax_f,
        caf_dot,
        "features.plot",
        "id",
        caf_features,
        CAF_ORDER,
        "CAF-subtype marker panel",
        "PuOr_r",
        None,
    )
    ax_f.set_yticklabels(["FAP+ infCAF", "ecmCAF", "EndMT", "myoCAF", "adiCAF"])
    panel_label(ax_f, "f")

    plot_fraction_context(
        ax_g,
        "Supplementary_Fig1D_sample_level_fibro_fraction.csv",
        "fraction",
        "Fibroblast fraction by sample",
        "Fraction",
    )
    panel_label(ax_g, "g")

    plot_fraction_context(
        ax_h,
        "Supplementary_Fig1E_sample_level_fap_infCAF_fraction.csv",
        "fraction",
        "FAP+ infCAF fraction by sample",
        "Fraction",
    )
    panel_label(ax_h, "h")

    plot_resolution_scan_compact(ax_i)
    panel_label(ax_i, "i")

    base.plot_panel_h(ax_j)
    ax_j.set_title("CAF subtype composition by sample", loc="left", fontweight="bold", pad=2)
    panel_label(ax_j, "j", x=-0.055, y=1.15)

    plot_marker_effects(ax_k)
    panel_label(ax_k, "k", x=-0.055, y=1.12)

    fig.subplots_adjust(top=0.975, bottom=0.055, left=0.075, right=0.90)
    fig.savefig(SVG_OUT, bbox_inches="tight")
    fig.savefig(PDF_OUT, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    build_figure()
    print(SVG_OUT)
    print(PDF_OUT)
    print(PNG_OUT)
