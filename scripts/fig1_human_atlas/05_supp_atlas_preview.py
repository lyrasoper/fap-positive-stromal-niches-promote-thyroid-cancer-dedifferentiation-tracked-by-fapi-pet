#!/usr/bin/env python3
# Cancer Research submission - figure code release
# Builds: Supplementary figures for Figure 1 (layout preview).

"""Build a Nature-style Python preview redesign for Supplementary Fig. 1.

This script intentionally writes to a preview folder and does not overwrite the
current submission package.
"""

from __future__ import annotations
import os

from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec


ROOT = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>"))
SRC = ROOT / "outputs/supplementary_figures/source_data"
META = ROOT / "outputs/fig1_scRNA/source_data/fig1_scrna_full_metadata.csv"
OUT = ROOT / "outputs/supplementary_figures/nature_preview_20260517"
OUT.mkdir(parents=True, exist_ok=True)

PNG_OUT = OUT / "Supplementary_Fig_1_harmonized_Nature_style_no_header_preview.png"
PDF_OUT = OUT / "Supplementary_Fig_1_harmonized_Nature_style_no_header_preview.pdf"


TISSUE_ORDER = ["NT", "PTC", "ATC"]
CELL_ORDER = ["Thyro", "Fibro", "Endo", "Myeloid", "T/NK", "B"]
CAF_ORDER = ["FAP+ infCAF", "ecmCAF", "EndMT CAF", "RGS15+ myoCAF", "adiCAF"]
STUDY_ORDER = ["GSE148673", "GSE184362", "GSE193581", "GSE210347"]

TISSUE_COLORS = {"NT": "#4C78A8", "PTC": "#59A14F", "ATC": "#E15759"}
CELL_COLORS = {
    "Thyro": "#4E79A7",
    "Fibro": "#E15759",
    "Endo": "#76B7B2",
    "Myeloid": "#B07AA1",
    "T/NK": "#F28E2B",
    "B": "#EDC948",
}
CAF_COLORS = {
    "FAP+ infCAF": "#D55E00",
    "ecmCAF": "#CC79A7",
    "EndMT CAF": "#009E73",
    "RGS15+ myoCAF": "#0072B2",
    "adiCAF": "#E69F00",
}
STUDY_COLORS = {
    "GSE148673": "#3B6EA8",
    "GSE184362": "#E17C05",
    "GSE193581": "#5E9F52",
    "GSE210347": "#B94E48",
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 5.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.45,
            "axes.labelsize": 5.9,
            "axes.titlesize": 6.4,
            "xtick.labelsize": 5.2,
            "ytick.labelsize": 5.2,
            "legend.fontsize": 5.0,
            "legend.title_fontsize": 5.2,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )
    sns.set_theme(style="white", rc={"axes.linewidth": 0.45})


def clean_axis(ax: plt.Axes, *, grid: bool = False) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", length=2.4, width=0.45, pad=1.5)
    if grid:
        ax.grid(axis="y", color="#D9D9D9", lw=0.35, alpha=0.7)
        ax.set_axisbelow(True)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.16,
        1.15,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.8,
        fontweight="bold",
    )


def format_p(p: float) -> str:
    if not np.isfinite(p):
        return "P=NA"
    if p < 1e-4:
        return f"P={p:.1e}"
    if p < 0.01:
        return f"P={p:.3f}"
    return f"P={p:.2f}"


def add_bracket(ax: plt.Axes, x1: int, x2: int, y: float, text: str) -> None:
    ylim = ax.get_ylim()
    h = (ylim[1] - ylim[0]) * 0.025
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=0.45)
    ax.text((x1 + x2) / 2, y + h * 1.25, text, ha="center", va="bottom", fontsize=5.2)


def dotplot(
    ax: plt.Axes,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    x_order: list[str],
    y_order: list[str],
    title: str,
    cmap: str,
    cbar_label: str | None = "Scaled exp.",
) -> None:
    data = df.copy()
    data[x_col] = pd.Categorical(data[x_col], categories=x_order, ordered=True)
    data[y_col] = pd.Categorical(data[y_col], categories=y_order, ordered=True)
    data = data.dropna(subset=[x_col, y_col]).sort_values([y_col, x_col])

    x_codes = data[x_col].cat.codes.to_numpy()
    y_codes = data[y_col].cat.codes.to_numpy()
    sizes = np.clip(data["pct.exp"].to_numpy(), 0, 100) * 0.95 + 5
    norm = Normalize(vmin=-2.0, vmax=2.0)
    sc = ax.scatter(
        x_codes,
        y_codes,
        s=sizes,
        c=data["avg.exp.scaled"],
        cmap=cmap,
        norm=norm,
        edgecolor="0.25",
        linewidth=0.18,
    )
    ax.set_xticks(range(len(x_order)))
    ax.set_xticklabels(x_order, rotation=45, ha="right")
    ax.set_yticks(range(len(y_order)))
    ax.set_yticklabels(y_order)
    ax.set_xlim(-0.55, len(x_order) - 0.45)
    ax.set_ylim(-0.55, len(y_order) - 0.45)
    ax.invert_yaxis()
    ax.set_title(title, loc="left", fontweight="bold", pad=2)
    ax.set_xlabel("")
    ax.set_ylabel("")
    clean_axis(ax)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=0)
    ax.grid(color="#ECECEC", lw=0.35)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.035, pad=0.015)
    cbar.outline.set_linewidth(0.35)
    cbar.ax.tick_params(length=1.8, width=0.35, labelsize=5.2)
    if cbar_label:
        cbar.set_label(cbar_label, fontsize=4.9, labelpad=1)


def plot_box_fraction(
    ax: plt.Axes,
    df: pd.DataFrame,
    y_col: str,
    title: str,
    ylabel: str,
    ylim_pad: float = 0.11,
) -> None:
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
        linewidth=0.55,
        saturation=0.95,
        boxprops={"alpha": 0.88},
        medianprops={"color": "black", "linewidth": 0.65},
        whiskerprops={"linewidth": 0.55},
        capprops={"linewidth": 0.55},
    )
    sns.stripplot(
        data=df,
        x="TissueType",
        y=y_col,
        order=TISSUE_ORDER,
        ax=ax,
        color="black",
        size=2.0,
        jitter=0.18,
        alpha=0.62,
        linewidth=0,
    )
    ax.set_title(title, loc="left", fontweight="bold", pad=2)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, min(1.0, df[y_col].max() * (1.0 + ylim_pad) + 0.04))
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1.0, decimals=0))
    clean_axis(ax, grid=True)

    ymax = ax.get_ylim()[1]
    comparisons = [("NT", "ATC"), ("PTC", "ATC")]
    levels = [ymax * 0.82, ymax * 0.91]
    for (left, right), y in zip(comparisons, levels):
        a = df.loc[df["TissueType"] == left, y_col].to_numpy()
        b = df.loc[df["TissueType"] == right, y_col].to_numpy()
        if len(a) > 0 and len(b) > 0:
            p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
            add_bracket(ax, TISSUE_ORDER.index(left), TISSUE_ORDER.index(right), y, format_p(p))


def plot_panel_a(ax: plt.Axes) -> None:
    meta = pd.read_csv(META, usecols=["TissueType", "Celltype_short"])
    meta = meta[meta["TissueType"].isin(TISSUE_ORDER) & meta["Celltype_short"].isin(CELL_ORDER)]
    comp = (
        meta.groupby(["TissueType", "Celltype_short"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    comp["prop"] = comp["n"] / comp.groupby("TissueType")["n"].transform("sum")
    pivot = comp.pivot(index="TissueType", columns="Celltype_short", values="prop")
    pivot = pivot.reindex(TISSUE_ORDER).reindex(columns=CELL_ORDER).fillna(0)

    bottom = np.zeros(len(TISSUE_ORDER))
    x = np.arange(len(TISSUE_ORDER))
    for cell in CELL_ORDER:
        vals = pivot[cell].to_numpy()
        ax.bar(x, vals, bottom=bottom, color=CELL_COLORS[cell], width=0.66, edgecolor="white", lw=0.25)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(TISSUE_ORDER)
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1.0, decimals=0))
    ax.set_title("Cell composition by disease state", loc="left", fontweight="bold", pad=2)
    ax.set_ylabel("Cells")
    clean_axis(ax, grid=True)
    handles = [mpl.patches.Patch(facecolor=CELL_COLORS[c], label=c) for c in CELL_ORDER]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.23),
        ncol=3,
        frameon=False,
        handlelength=0.9,
        handletextpad=0.3,
        columnspacing=0.6,
    )


def plot_panel_d(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Supplementary_Fig1I_FAP_celltype_expression.csv")
    df = df[df["Celltype_short"].isin(CELL_ORDER)].copy()
    pieces = []
    for _, group in df.groupby("Celltype_short", observed=True):
        pieces.append(group.sample(n=min(3500, len(group)), random_state=17))
    sampled = pd.concat(pieces, ignore_index=True)
    sampled["Celltype_short"] = pd.Categorical(sampled["Celltype_short"], categories=CELL_ORDER, ordered=True)
    sns.violinplot(
        data=sampled,
        x="Celltype_short",
        y="FAP",
        order=CELL_ORDER,
        palette=CELL_COLORS,
        inner=None,
        cut=0,
        linewidth=0.35,
        saturation=0.95,
        ax=ax,
    )
    sns.boxplot(
        data=sampled,
        x="Celltype_short",
        y="FAP",
        order=CELL_ORDER,
        width=0.16,
        color="white",
        fliersize=0,
        linewidth=0.45,
        medianprops={"color": "black", "linewidth": 0.55},
        ax=ax,
    )
    ax.set_title("FAP expression across lineages", loc="left", fontweight="bold", pad=2)
    ax.set_xlabel("")
    ax.set_ylabel("FAP")
    ax.set_ylim(0, np.nanpercentile(sampled["FAP"], 99.7) * 1.08)
    ax.tick_params(axis="x", rotation=30)
    clean_axis(ax, grid=True)


def plot_panel_e(ax: plt.Axes) -> None:
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
        ax.bar(x, vals, bottom=bottom, color=STUDY_COLORS[study], width=0.74, edgecolor="white", lw=0.25)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(CELL_ORDER, rotation=30, ha="right")
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1.0, decimals=0))
    ax.set_title("Cohort mixing across lineages", loc="left", fontweight="bold", pad=2)
    ax.set_ylabel("Cells")
    clean_axis(ax, grid=True)
    handles = [mpl.patches.Patch(facecolor=STUDY_COLORS[s], label=s) for s in STUDY_ORDER]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.26),
        ncol=2,
        frameon=False,
        handlelength=0.9,
        handletextpad=0.3,
        columnspacing=0.6,
    )


def plot_panel_f_qc(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Supplementary_Fig1G_sample_level_qc_summary.csv")
    df = df[df["TissueType"].isin(TISSUE_ORDER)].copy()
    metric_order = ["nFeature_RNA", "nCount_RNA (log10)", "percent.mt (%)"]
    metric_labels = ["Detected\ngenes", "UMI count\n(log10)", "Mitochondrial\nreads"]
    df["metric"] = pd.Categorical(df["metric"], categories=metric_order, ordered=True)
    df = df.dropna(subset=["metric"])
    df["z"] = df.groupby("metric", observed=True)["median"].transform(
        lambda s: (s - s.mean()) / (s.std(ddof=0) + 1e-9)
    )
    sns.boxplot(
        data=df,
        x="metric",
        y="z",
        hue="TissueType",
        order=metric_order,
        hue_order=TISSUE_ORDER,
        palette=TISSUE_COLORS,
        ax=ax,
        width=0.62,
        fliersize=0,
        linewidth=0.45,
        saturation=0.9,
        medianprops={"color": "black", "linewidth": 0.55},
    )
    sns.stripplot(
        data=df,
        x="metric",
        y="z",
        hue="TissueType",
        order=metric_order,
        hue_order=TISSUE_ORDER,
        dodge=True,
        ax=ax,
        color="black",
        size=1.6,
        alpha=0.48,
        linewidth=0,
        legend=False,
    )
    ax.axhline(0, color="#777777", lw=0.45, ls="--", zorder=0)
    ax.set_xticklabels(metric_labels)
    ax.set_title("Sample-level QC profile", loc="left", fontweight="bold", pad=2)
    ax.set_xlabel("")
    ax.set_ylabel("z-scored median")
    clean_axis(ax, grid=True)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[:3],
        labels[:3],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.26),
        ncol=3,
        frameon=False,
        handlelength=0.9,
        handletextpad=0.3,
        columnspacing=0.6,
    )


def plot_panel_i_resolution(ax: plt.Axes) -> None:
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
        ax.bar(x, vals, bottom=bottom, color=color, edgecolor="white", linewidth=0.15, width=0.72)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(resolutions)
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1.0, decimals=0))
    ax.set_title("CAF reclustering resolution scan", loc="left", fontweight="bold", pad=2)
    ax.set_xlabel("Seurat resolution")
    ax.set_ylabel("CAF cells")
    clean_axis(ax, grid=True)
    unique_counts = long.groupby("resolution", observed=True)["cluster"].nunique().reindex(resolutions)
    for xpos, n in enumerate(unique_counts):
        ax.text(xpos, 1.025, f"{int(n)}", ha="center", va="bottom", fontsize=5.0)
    ax.text(-0.45, 1.025, "clusters", ha="right", va="bottom", fontsize=5.0, color="#555555")


def plot_panel_h(ax: plt.Axes) -> None:
    df = pd.read_csv(SRC / "Supplementary_Fig1K_per_sample_caf_composition.csv")
    df = df[df["CAF_clusters"].isin(CAF_ORDER)].copy()
    fap_prop = (
        df[df["CAF_clusters"] == "FAP+ infCAF"]
        .set_index("Sample_id")["prop"]
        .to_dict()
    )
    sample_meta = df[["Sample_id", "TissueType"]].drop_duplicates()
    sample_meta["TissueType"] = pd.Categorical(sample_meta["TissueType"], categories=TISSUE_ORDER, ordered=True)
    sample_meta["fap_prop"] = sample_meta["Sample_id"].map(fap_prop).fillna(0)
    sample_meta = sample_meta.sort_values(["TissueType", "fap_prop", "Sample_id"])
    samples = sample_meta["Sample_id"].tolist()
    pivot = (
        df.pivot(index="Sample_id", columns="CAF_clusters", values="prop")
        .reindex(samples)
        .reindex(columns=CAF_ORDER)
        .fillna(0)
    )
    x = np.arange(len(samples))
    bottom = np.zeros(len(samples))
    for caf in CAF_ORDER:
        vals = pivot[caf].to_numpy()
        ax.bar(x, vals, bottom=bottom, color=CAF_COLORS[caf], width=0.88, edgecolor="white", lw=0.12)
        bottom += vals

    boundaries = []
    labels = []
    start = 0
    for tissue in TISSUE_ORDER:
        n = int((sample_meta["TissueType"] == tissue).sum())
        if n:
            labels.append((start + n / 2 - 0.5, tissue))
            start += n
            boundaries.append(start - 0.5)
    for b in boundaries[:-1]:
        ax.axvline(b, color="black", lw=0.45)
    for xpos, label in labels:
        ax.text(xpos, -0.13, label, ha="center", va="top", transform=ax.get_xaxis_transform(), fontsize=6.0)

    ax.set_xlim(-0.6, len(samples) - 0.4)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1.0, decimals=0))
    ax.set_title("Per-sample CAF subtype composition", loc="left", fontweight="bold", pad=2)
    ax.set_ylabel("CAF cells")
    clean_axis(ax, grid=True)
    handles = [mpl.patches.Patch(facecolor=CAF_COLORS[c], label=c) for c in CAF_ORDER]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.24),
        ncol=5,
        frameon=False,
        handlelength=0.9,
        handletextpad=0.3,
        columnspacing=0.8,
    )


def build_figure() -> None:
    configure_style()

    fig = plt.figure(figsize=(12.0, 16.2), dpi=300)
    gs = GridSpec(
        nrows=5,
        ncols=6,
        figure=fig,
        height_ratios=[1.02, 1.02, 1.08, 0.86, 1.02],
        hspace=0.98,
        wspace=0.58,
    )

    ax_a = fig.add_subplot(gs[0, 0:2])
    ax_b = fig.add_subplot(gs[0, 2:6])
    ax_c = fig.add_subplot(gs[1, 0:2])
    ax_d = fig.add_subplot(gs[1, 2:4])
    ax_e = fig.add_subplot(gs[1, 4:6])
    ax_f = fig.add_subplot(gs[2, 0:2])
    ax_g = fig.add_subplot(gs[2, 2:6])
    ax_h = fig.add_subplot(gs[3, 0:2])
    ax_i = fig.add_subplot(gs[3, 2:6])
    ax_j = fig.add_subplot(gs[4, 0:6])

    plot_panel_a(ax_a)
    panel_label(ax_a, "a")

    lineage_dot = pd.read_csv(SRC / "Supplementary_Fig1C_major_lineage_dotplot_data.csv")
    lineage_features = ["TG", "TPO", "PAX8", "FAP", "COL1A1", "PECAM1", "KDR", "LST1", "C1QC", "CD3D", "NKG7", "MS4A1", "CD79A"]
    dotplot(
        ax_b,
        lineage_dot,
        "features.plot",
        "id",
        lineage_features,
        CELL_ORDER,
        "Canonical lineage markers",
        "RdBu_r",
        "Scaled exp.",
    )
    panel_label(ax_b, "b")

    fibro = pd.read_csv(SRC / "Supplementary_Fig1D_sample_level_fibro_fraction.csv")
    plot_box_fraction(ax_c, fibro, "fraction", "Fibroblast enrichment", "Fraction")
    panel_label(ax_c, "c")

    plot_panel_d(ax_d)
    panel_label(ax_d, "d")

    plot_panel_e(ax_e)
    panel_label(ax_e, "e")

    plot_panel_f_qc(ax_f)
    panel_label(ax_f, "f")

    caf_dot = pd.read_csv(SRC / "Supplementary_Fig1F_caf_dotplot_data.csv")
    caf_features = ["FAP", "CXCL12", "IL6", "COL1A1", "POSTN", "PLVAP", "RAMP2", "RGS5", "ACTA2", "SFRP5", "OGN"]
    dotplot(
        ax_g,
        caf_dot,
        "features.plot",
        "id",
        caf_features,
        CAF_ORDER,
        "CAF subtype markers",
        "PuOr_r",
        None,
    )
    ax_g.set_yticklabels(["FAP+ infCAF", "ecmCAF", "EndMT", "myoCAF", "adiCAF"])
    panel_label(ax_g, "g")

    fap_inf = pd.read_csv(SRC / "Supplementary_Fig1E_sample_level_fap_infCAF_fraction.csv")
    plot_box_fraction(ax_h, fap_inf, "fraction", "FAP+ infCAF expansion", "Fraction")
    panel_label(ax_h, "h")

    plot_panel_i_resolution(ax_i)
    panel_label(ax_i, "i")

    plot_panel_h(ax_j)
    panel_label(ax_j, "j")

    fig.subplots_adjust(top=0.975, bottom=0.075, left=0.075, right=0.985)
    fig.savefig(PNG_OUT, dpi=300)
    fig.savefig(PDF_OUT)
    plt.close(fig)


if __name__ == "__main__":
    build_figure()
    print(PNG_OUT)
    print(PDF_OUT)
