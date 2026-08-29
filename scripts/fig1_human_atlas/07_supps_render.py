# Cancer Research submission - figure code release
# Builds: Supplementary-figure renderer for Figures 1-2 and the mouse
#          source panel of Supplementary Fig. S21. Internal "supp1-6" names
#          are this repository's own build order, not submission numbering.

"""Render Supps 1-6 in matplotlib (paper_A style), replacing the ggplot
versions produced by the corresponding R scripts.

All panel data is read from the pre-aggregated CSVs that the R scripts
exported into outputs/supplementary_figures/source_data/. Supp 6 (mouse
Seurat-only) requires an extra R bridge that exports the aggregates we
need; for this first pass we render what's available from CSV.

Usage:
    python3 scripts/fig1_human_atlas/07_supps_render.py              # all
    python3 scripts/fig1_human_atlas/07_supps_render.py supp3 supp5  # subset
"""
from __future__ import annotations
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).parent.parent / "_shared"))
from paper_A_style import (  # noqa: E402
    apply_pa_rc, panel_tag, style_axes, save, pa, pa_cmap_div, COL1, COL2,
)

ROOT = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>"))
SRC = ROOT / "outputs/supplementary_figures/source_data"
OUTDIR = ROOT / "outputs/supplementary_figures/paper_A_matplotlib"
OUTDIR.mkdir(parents=True, exist_ok=True)

apply_pa_rc()


# ── shared utilities ────────────────────────────────────────────────────────
def fmt_p(p: float) -> str:
    if p < 1e-100:
        return "P < 10$^{-100}$"
    if p < 0.001:
        exp = int(np.floor(np.log10(p)))
        return f"P < 10$^{{{exp+1}}}$"
    return f"P = {p:.3f}"


def fig_title(fig, n: int, text: str) -> None:
    fig.suptitle(f"Supplementary Fig. {n} | {text}",
                 fontsize=8, fontweight="bold",
                 y=0.995, x=0.04, ha="left")


def cbar_from(ax, mappable, label: str, *, fraction=0.04, pad=0.04,
              label_top=False):
    """Attach a slim colorbar. label_top=True places the label horizontally
    above the bar instead of as a rotated side label — use it when the panel
    sits to the left of another panel, so the rotated label can't collide
    with the neighbour's y-axis label.
    """
    cb = ax.figure.colorbar(mappable, ax=ax, fraction=fraction, pad=pad)
    cb.ax.tick_params(labelsize=5.2, width=0.4, length=2)
    if label_top:
        cb.ax.set_title(label, fontsize=5.8, pad=2)
    else:
        cb.set_label(label, fontsize=6)
    cb.outline.set_linewidth(0.4)
    return cb


def dotplot_size_legend(ax, sizes=(25, 50, 75, 100), max_pct=100,
                        max_size=55, base_size=4, loc=(1.18, 1.0)):
    """Anchored size legend for Seurat-style dot plots."""
    handles = [plt.scatter([], [], s=(p / max_pct) * max_size + base_size,
                           c="lightgrey", edgecolor="black", linewidth=0.2,
                           label=f"{p}%") for p in sizes]
    leg = ax.legend(handles=handles, title="% expressed",
                    loc="upper left", bbox_to_anchor=loc,
                    fontsize=4.6, title_fontsize=5,
                    handletextpad=0.3, borderpad=0.3,
                    labelspacing=0.6, frameon=False)
    ax.add_artist(leg)


def boxplot_by(ax, df, group_col: str, value_col: str, order, palette,
               *, jitter=True, rng=None):
    rng = rng or np.random.default_rng(42)
    data = [df.loc[df[group_col] == g, value_col].dropna().values for g in order]
    bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", linewidth=0.9),
                    whiskerprops=dict(color="black", linewidth=0.55),
                    capprops=dict(color="black", linewidth=0.55),
                    boxprops=dict(linewidth=0.55))
    for box, g in zip(bp["boxes"], order):
        box.set_facecolor(palette.get(g, pa["scatter_pt"]))
        box.set_alpha(0.85)
    if jitter:
        for i, g in enumerate(order, 1):
            v = data[i-1]
            x = i + (rng.random(len(v)) - 0.5) * 0.25
            ax.scatter(x, v, s=6, facecolor="black", alpha=0.5,
                       edgecolor="none")
    ax.set_xticks(range(1, len(order)+1))
    ax.set_xticklabels(order, rotation=18, ha="right")
    style_axes(ax)


def roc_panel(ax, df, model_col="model", fpr_col="fpr", tpr_col="tpr",
              palette_cycle=None):
    palette_cycle = palette_cycle or [
        pa["mode"]["PAX8-loss"], pa["mode"]["SERPINE1-gain"],
        pa["mode"]["Dual-axis"], pa["integrin"]["avb3"],
        pa["integrin"]["a5"], "#7CBF60",
    ]
    for i, (m, sub) in enumerate(df.groupby(model_col, sort=False)):
        sub = sub.sort_values(fpr_col)
        ax.plot(sub[fpr_col], sub[tpr_col],
                color=palette_cycle[i % len(palette_cycle)],
                linewidth=1.0, label=m)
    ax.plot([0, 1], [0, 1], ls="--", color="grey", linewidth=0.5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(loc="lower right", fontsize=5.2,
              handletextpad=0.3, borderpad=0.2, labelspacing=0.25)
    style_axes(ax)


# Palettes shared across Supps
TISSUE_PAL = pa["tissue_human"]
TISSUE_ORDER = ["NT", "PTC", "ATC"]


# ═══════════════════════════════════════════════════════════════════════════
# Supp 1 — Additional validation of the single-cell atlas + CAF annotation
# ═══════════════════════════════════════════════════════════════════════════
def supp1() -> None:
    """Six panels of EXTENDED validation that do NOT duplicate main Fig 1:
       a — Major lineage marker dot-plot (broad marker set; Fig 1 D only shows FAP/TDS)
       b — CAF subtype marker dot-plot (Fig 1 E shows UMAP; this shows marker specificity)
       c — Sample-level fibroblast fraction (per-sample variance; Fig 1 C only shows averages)
       d — Per-sample CAF composition with tissue stripe (Fig 1 G shows pooled averages)
       e — Cross-cohort cell-type mixing (QC, not in Fig 1)
       f — Top CAF subtype markers heatmap (not in Fig 1)
    """
    fig = plt.figure(figsize=(COL2, COL2 * 1.05))
    gs = gridspec.GridSpec(
        4, 6, figure=fig,
        height_ratios=[0.95, 0.95, 1.05, 1.00],
        hspace=1.05, wspace=1.30,
        left=0.07, right=0.96, top=0.955, bottom=0.12,
    )

    # ── (a) Major lineage marker dot-plot — row 1, full width ──
    ax_a = fig.add_subplot(gs[0, :])
    dp = pd.read_csv(SRC / "Supplementary_Fig1C_major_lineage_dotplot_data.csv")
    feats = list(pd.Index(dp["features.plot"]).unique())
    ids = list(pd.Index(dp["id"]).unique())
    Xc = {f: i for i, f in enumerate(feats)}
    Yc = {g: j for j, g in enumerate(ids)}
    sc = ax_a.scatter(
        dp["features.plot"].map(Xc), dp["id"].map(Yc),
        s=(dp["pct.exp"]/dp["pct.exp"].max())*65 + 4,
        c=dp["avg.exp.scaled"], cmap=pa_cmap_div(),
        vmin=-2, vmax=2, edgecolor="black", linewidth=0.2,
    )
    ax_a.set_xticks(range(len(feats)))
    ax_a.set_xticklabels(feats, rotation=45, ha="right", fontsize=6.0)
    ax_a.set_yticks(range(len(ids)))
    ax_a.set_yticklabels(ids, fontsize=6.2)
    ax_a.set_xlim(-0.5, len(feats) - 0.5)
    ax_a.set_ylim(-0.5, len(ids) - 0.5)
    cbar_from(ax_a, sc, "scaled avg expr")
    panel_tag(ax_a, "a", "Extended lineage marker validation across all cell types",
              dx=-12, dy=4)

    # ── (b) CAF subtype marker dot-plot — row 2, cols 0-3 ──
    ax_b = fig.add_subplot(gs[1, :4])
    caf = pd.read_csv(SRC / "Supplementary_Fig1F_caf_dotplot_data.csv")
    feats2 = list(pd.Index(caf["features.plot"]).unique())
    ids2 = list(pd.Index(caf["id"]).unique())
    Xc2 = {f: i for i, f in enumerate(feats2)}
    Yc2 = {g: j for j, g in enumerate(ids2)}
    sc2 = ax_b.scatter(
        caf["features.plot"].map(Xc2), caf["id"].map(Yc2),
        s=(caf["pct.exp"]/caf["pct.exp"].max())*65 + 4,
        c=caf["avg.exp.scaled"], cmap=pa_cmap_div(),
        vmin=-2, vmax=2, edgecolor="black", linewidth=0.2,
    )
    ax_b.set_xticks(range(len(feats2)))
    ax_b.set_xticklabels(feats2, rotation=45, ha="right", fontsize=6.0)
    ax_b.set_yticks(range(len(ids2)))
    ax_b.set_yticklabels(ids2, fontsize=6.2)
    ax_b.set_xlim(-0.5, len(feats2) - 0.5)
    ax_b.set_ylim(-0.5, len(ids2) - 0.5)
    cbar_from(ax_b, sc2, "scaled avg expr")
    panel_tag(ax_b, "b", "CAF subtype marker specificity", dx=-12, dy=4)

    # ── (c) Sample-level fibroblast fraction (boxplot, per-sample variance) ──
    ax_c = fig.add_subplot(gs[1, 4:])
    fib = pd.read_csv(SRC / "Supplementary_Fig1D_sample_level_fibro_fraction.csv")
    boxplot_by(ax_c, fib, "TissueType", "fraction", TISSUE_ORDER, TISSUE_PAL)
    ax_c.set_ylabel("Fibroblast fraction")
    panel_tag(ax_c, "c", "Per-sample fibroblast abundance", dx=-10, dy=4)

    # ── (d) Per-sample CAF composition — row 3, full width ──
    ax_d = fig.add_subplot(gs[2, :])
    comp = pd.read_csv(SRC / "Supplementary_Fig1K_per_sample_caf_composition.csv")
    sample_tissue = comp.drop_duplicates("Sample_id")[["Sample_id", "TissueType"]].copy()
    sample_tissue["TissueType"] = pd.Categorical(
        sample_tissue["TissueType"], categories=TISSUE_ORDER, ordered=True)
    sample_tissue = sample_tissue.sort_values(["TissueType", "Sample_id"]).reset_index(drop=True)
    samples = sample_tissue["Sample_id"].tolist()
    tissues = sample_tissue["TissueType"].astype(str).tolist()
    cafs = list(pd.Index(comp["CAF_clusters"]).unique())
    M2 = comp.pivot_table(index="Sample_id", columns="CAF_clusters",
                          values="prop", fill_value=0).reindex(samples)
    bottoms = np.zeros(len(samples))
    caf_pal = pa["caf"]
    caf_fallback = ["#E8A341", "#7E5CAB", "#16A085", "#3262A0", "#7CBF60", "#C03840"]
    for ci, c in enumerate(cafs):
        col = caf_pal.get(c, caf_fallback[ci % len(caf_fallback)])
        ax_d.bar(range(len(samples)), M2[c].values, bottom=bottoms,
                 color=col, edgecolor="white", linewidth=0.15, label=c)
        bottoms += M2[c].values
    ax_d.set_ylim(0, 1.04)
    ax_d.set_ylabel("CAF composition")
    stripe_height = 0.035
    for i, t in enumerate(tissues):
        ax_d.add_patch(Rectangle((i - 0.4, -stripe_height * 1.5),
                                 0.8, stripe_height,
                                 facecolor=TISSUE_PAL.get(t, "#999999"),
                                 edgecolor="none", clip_on=False, zorder=3))
    ax_d.set_xticks([])
    from collections import Counter
    tcount = Counter(tissues)
    cum = 0
    for t in TISSUE_ORDER:
        n_t = tcount.get(t, 0)
        if n_t == 0:
            continue
        mid = cum + n_t / 2 - 0.5
        ax_d.text(mid, -stripe_height * 4.5, f"{t} (n={n_t})",
                  ha="center", va="top", fontsize=6.5, fontweight="bold",
                  color=TISSUE_PAL.get(t, "#1A1A1A"))
        cum += n_t
    ax_d.set_xlim(-0.6, len(samples) - 0.4)
    ax_d.legend(fontsize=5.0, loc="center left", bbox_to_anchor=(1.005, 0.5),
                handletextpad=0.3, borderpad=0.2, labelspacing=0.25,
                title="CAF subtype", title_fontsize=5.4)
    style_axes(ax_d)
    panel_tag(ax_d, "d", "Per-sample CAF composition (inter-patient variability)",
              dx=-12, dy=4)

    # ── (e) Cross-cohort cell-type mixing — row 4, cols 0-2 ──
    ax_e = fig.add_subplot(gs[3, :3])
    mix = pd.read_csv(SRC / "Supplementary_Fig1H_celltype_cohort_mixing.csv")
    cells = list(pd.Index(mix["Celltype_short"]).unique())
    datasets = list(pd.Index(mix["Data.ID"]).unique())
    M = mix.pivot_table(index="Celltype_short", columns="Data.ID",
                        values="prop", fill_value=0).reindex(index=cells, columns=datasets)
    DATASET_PAL = ["#3262A0", "#7E5CAB", "#E8A341", "#1F6E3A",
                   "#C03840", "#5786C4", "#7CBF60", "#9E9E9E"]
    bottoms = np.zeros(len(cells))
    for j, d in enumerate(datasets):
        ax_e.bar(range(len(cells)), M[d].values, bottom=bottoms,
                 color=DATASET_PAL[j % len(DATASET_PAL)],
                 edgecolor="white", linewidth=0.3, label=d)
        bottoms += M[d].values
    ax_e.set_xticks(range(len(cells)))
    ax_e.set_xticklabels(cells, rotation=45, ha="right", fontsize=6.0)
    ax_e.set_ylim(0, 1)
    ax_e.set_ylabel("Proportion")
    # legend horizontal BELOW the panel — keeps it out of panel f's y-axis
    # labels, which sit on panel f's left edge in the same row.
    ax_e.legend(fontsize=4.6, loc="upper center", bbox_to_anchor=(0.5, -0.42),
                ncols=min(len(datasets), 4), handletextpad=0.3, borderpad=0.3,
                columnspacing=0.9, labelspacing=0.3,
                title="Dataset", title_fontsize=5.0, frameon=False)
    style_axes(ax_e)
    panel_tag(ax_e, "e", "Cross-cohort cell-type mixing (QC)", dx=-12, dy=4)

    # ── (f) Top CAF subtype markers heatmap — row 4, cols 3-5 ──
    # Transposed layout: ~20 marker genes on the x-axis (rotated 90°, each
    # gene gets a full column of room) and the 5 CAF subtypes on the y-axis.
    # A wide-and-short heatmap is legible where ~20 gene *rows* crammed into
    # a quarter-panel height collapse into an unreadable label blob. Genes
    # stay grouped by their cluster so the block-diagonal structure still
    # reads — it just runs left-to-right instead of top-to-bottom.
    ax_f = fig.add_subplot(gs[3, 3:])
    mk = pd.read_csv(SRC / "Supplementary_Fig1L_top_caf_markers.csv")
    mk_top = (mk.sort_values("avg_log2FC", ascending=False)
                .groupby("cluster").head(4)
                .reset_index(drop=True))
    cluster_order = list(mk_top["cluster"].drop_duplicates())
    mk_top["cluster"] = pd.Categorical(mk_top["cluster"], categories=cluster_order)
    mk_top = mk_top.sort_values(["cluster", "avg_log2FC"], ascending=[True, False])
    gene_order = mk_top["gene"].tolist()
    pv = mk_top.pivot_table(index="gene", columns="cluster",
                            values="avg_log2FC", fill_value=0, observed=False)
    # reindex to cluster-grouped gene order, then transpose → clusters × genes
    pv = pv.reindex(index=gene_order, columns=cluster_order).T
    vmax = float(np.abs(pv.values).max())
    im = ax_f.imshow(pv.values, cmap=pa_cmap_div(), vmin=-vmax, vmax=vmax,
                     aspect="auto", interpolation="nearest")
    ax_f.set_xticks(range(pv.shape[1]))
    ax_f.set_xticklabels(pv.columns, rotation=90, fontsize=5.4)
    ax_f.set_yticks(range(pv.shape[0]))
    ax_f.set_yticklabels(pv.index, fontsize=6.2)
    for ii in range(pv.shape[0]):
        for jj in range(pv.shape[1]):
            v = pv.values[ii, jj]
            if abs(v) >= 1.0:
                ax_f.text(jj, ii, f"{v:.1f}", ha="center", va="center",
                          fontsize=4.4,
                          color="white" if abs(v) > vmax * 0.6 else "#1A1A1A")
    cbar_from(ax_f, im, "log2FC")
    panel_tag(ax_f, "f", "Top CAF subtype-specific markers (top-4 by log2FC)",
              dx=-12, dy=4)

    fig_title(fig, 1, "Extended validation of the single-cell atlas and CAF annotation (orthogonal to Fig. 1)")
    save(fig, "Supplementary_Fig_1_Additional_validation_of_the_single_cell_atlas_and_CAF_annotation", OUTDIR)
    print("  ✓ Supp 1 rendered (6-panel, de-conflicted from Fig 1)")


# ═══════════════════════════════════════════════════════════════════════════
# Supp 2 — Extended TCGA-THCA validation
# ═══════════════════════════════════════════════════════════════════════════
def supp2() -> None:
    fig = plt.figure(figsize=(COL2, COL2 * 1.25))
    gs = gridspec.GridSpec(4, 6, figure=fig,
                           hspace=0.85, wspace=1.70,
                           left=0.06, right=0.985, top=0.955, bottom=0.06)

    # (a) bubble: gene × spearman ρ with FAP
    ax_a = fig.add_subplot(gs[0, :3])
    bub = pd.read_csv(SRC / "Supplementary_Fig2A_bubble_source_data.csv")
    bub = bub.sort_values("Correlation")
    cols_a = [pa["mode"]["PAX8-loss"] if d == "Negative" else pa["mode"]["SERPINE1-gain"] for d in bub["Direction"]]
    ax_a.scatter(bub["Correlation"], range(len(bub)),
                      s=-np.log10(bub["P_value"].clip(lower=1e-100))*2 + 8,
                      c=cols_a, edgecolor="black", linewidth=0.3, alpha=0.9)
    ax_a.axvline(0, color="grey", lw=0.4, ls="--")
    ax_a.set_yticks(range(len(bub))); ax_a.set_yticklabels(bub["Gene"], fontsize=5.4)
    ax_a.set_xlabel("Spearman ρ vs FAP")
    style_axes(ax_a)
    panel_tag(ax_a, "a", "TCGA-THCA: marker-FAP correlation", dx=-12, dy=4)

    # (b) ECM genes vs FAP bar
    ax_b = fig.add_subplot(gs[0, 3:])
    ecm = pd.read_csv(SRC / "Supplementary_Fig2E_ECMgene_FAP_spearman.csv")
    ecm = ecm.sort_values("rho")
    ax_b.barh(range(len(ecm)), ecm["rho"], color=pa["mode"]["SERPINE1-gain"],
              edgecolor="black", linewidth=0.3)
    ax_b.set_yticks(range(len(ecm))); ax_b.set_yticklabels(ecm["Gene"], fontsize=5.4)
    ax_b.set_xlabel("Spearman ρ vs FAP")
    style_axes(ax_b)
    panel_tag(ax_b, "b", "ECM gene–FAP correlations", dx=-12, dy=4)

    # (c) FAP correlation with stromal ssGSEA modules — bar (not a sparse 3×3)
    ax_c = fig.add_subplot(gs[1, :3])
    mod = pd.read_csv(SRC / "Supplementary_Fig2F_FAP_module_ssGSEA.csv")
    if "FAP" in mod.columns:
        cor = (mod.drop(columns=["SAMPLE_ID"]).corr(method="spearman")
                  .loc["FAP"].drop("FAP").sort_values())
        bar_cols = [pa["mode"]["PAX8-loss"] if v < 0 else pa["mode"]["SERPINE1-gain"]
                    for v in cor.values]
        ax_c.barh(range(len(cor)), cor.values, color=bar_cols,
                  edgecolor="black", linewidth=0.3, height=0.62)
        ax_c.axvline(0, color="grey", lw=0.4, ls="--")
        ax_c.set_yticks(range(len(cor)))
        ax_c.set_yticklabels(cor.index, fontsize=6.2)
        ax_c.set_xlabel("Spearman ρ with FAP")
        ax_c.set_xlim(min(cor.min() * 1.25, -0.1), max(cor.max() * 1.25, 0.1))
        for i, v in enumerate(cor.values):
            ax_c.text(v + (0.015 if v >= 0 else -0.015), i, f"{v:.2f}",
                      va="center", ha="left" if v >= 0 else "right",
                      fontsize=5.4)
    style_axes(ax_c)
    panel_tag(ax_c, "c", "FAP correlation with stromal modules", dx=-12, dy=4)

    # (d) OS KM-style: FAP-high vs FAP-low — Kaplan-Meier via step
    ax_d = fig.add_subplot(gs[1, 3:])
    osd = pd.read_csv(SRC / "Supplementary_Fig2H_FAP_OS_input.csv")
    osd = osd.dropna(subset=["OS_MONTHS", "status", "FAP_group"])
    for grp, sub in osd.groupby("FAP_group"):
        sub = sub.sort_values("OS_MONTHS")
        n = len(sub)
        at_risk = np.arange(n, 0, -1)
        surv = np.cumprod(np.where(sub["status"].values > 0, (at_risk - 1) / at_risk, 1.0))
        color = pa["mode"]["SERPINE1-gain"] if "high" in grp.lower() else pa["mode"]["PAX8-loss"]
        ax_d.step(sub["OS_MONTHS"].values, surv, where="post",
                  color=color, lw=1.0, label=f"{grp} (n={n})")
    ax_d.set_xlabel("Overall survival (months)")
    ax_d.set_ylabel("Survival probability")
    ax_d.legend(fontsize=5.2, loc="lower left",
                handletextpad=0.3, borderpad=0.2)
    style_axes(ax_d)
    panel_tag(ax_d, "d", "Overall survival by FAP tertile (TCGA-THCA)", dx=-12, dy=4)

    # (e) per-histology FAP-Spearman heatmap — transposed (genes on x, wide+short)
    ax_e = fig.add_subplot(gs[2, :3])
    ph = pd.read_csv(SRC / "Supplementary_Fig2J_perhistology_FAP_spearman.csv")
    ph_num = ph.set_index("Gene").select_dtypes(include=[np.number]).T
    im = ax_e.imshow(ph_num.values, cmap=pa_cmap_div(), vmin=-1, vmax=1,
                     aspect="auto")
    ax_e.set_xticks(range(ph_num.shape[1]))
    ax_e.set_xticklabels(ph_num.columns, rotation=90, fontsize=4.8)
    ax_e.set_yticks(range(ph_num.shape[0]))
    ax_e.set_yticklabels(ph_num.index, fontsize=6.2)
    cbar_from(ax_e, im, "Spearman ρ", fraction=0.032, label_top=True)
    panel_tag(ax_e, "e", "Per-histology Spearman vs FAP", dx=-12, dy=4)

    # (f) alt CAF signatures by histology
    ax_f = fig.add_subplot(gs[2, 3:])
    alt = pd.read_csv(SRC / "Supplementary_Fig2K_alt_CAF_signatures.csv")
    sigs = list(pd.Index(alt["Signature"]).unique())
    hist_order = sorted(alt["Histology"].dropna().unique())
    for i, sig in enumerate(sigs):
        sub = alt[alt["Signature"] == sig]
        means = sub.groupby("Histology")["Score"].mean().reindex(hist_order)
        sems = sub.groupby("Histology")["Score"].sem().reindex(hist_order)
        x = np.arange(len(hist_order)) + i*0.18 - (len(sigs)-1)*0.09
        ax_f.bar(x, means.values, width=0.16, yerr=sems.values,
                 capsize=1.2, edgecolor="black", linewidth=0.3,
                 error_kw=dict(linewidth=0.4),
                 label=sig)
    ax_f.set_xticks(range(len(hist_order)))
    ax_f.set_xticklabels(hist_order, rotation=20, ha="right", fontsize=5.6)
    ax_f.set_ylabel("ssGSEA score")
    ax_f.legend(fontsize=4.8, loc="center left", bbox_to_anchor=(1.01, 0.5),
                handletextpad=0.3, borderpad=0.3, labelspacing=0.3,
                title="Signature", title_fontsize=5.2)
    style_axes(ax_f)
    panel_tag(ax_f, "f", "Alternative CAF signatures by histology", dx=-12, dy=4)

    # (g) ROC: tall-cell vs classical (FAP/module)
    ax_g = fig.add_subplot(gs[3, :3])
    roc = pd.read_csv(SRC / "Supplementary_Fig2L_TallCell_vs_Classical_ROC.csv")
    roc_panel(ax_g, roc)
    panel_tag(ax_g, "g", "Tall cell vs classical: FAP/module ROC", dx=-12, dy=4)

    # (h) BRAF-like vs RAS-like comparison — grouped boxplot by class × readout
    ax_h = fig.add_subplot(gs[3, 3:])
    bf = pd.read_csv(SRC / "Supplementary_Fig2M_BRAF_RAS_class.csv")
    classes = sorted(bf["Class"].dropna().unique())
    readouts = sorted(bf["readout"].dropna().unique())
    readout_pal = [pa["mode"]["SERPINE1-gain"], pa["mode"]["PAX8-loss"],
                   pa["mode"]["Dual-axis"], pa["integrin"]["avb3"]]
    n_r = len(readouts)
    width = 0.8 / n_r
    for i, r in enumerate(readouts):
        positions = [j + (i - (n_r-1)/2) * width for j in range(len(classes))]
        data_r = [bf.loc[(bf["Class"] == c) & (bf["readout"] == r), "score"].dropna().values
                  for c in classes]
        bp = ax_h.boxplot(data_r, positions=positions, widths=width*0.85,
                          patch_artist=True, showfliers=False,
                          medianprops=dict(color="black", linewidth=0.7),
                          whiskerprops=dict(linewidth=0.45),
                          capprops=dict(linewidth=0.45),
                          boxprops=dict(linewidth=0.45))
        for box in bp["boxes"]:
            box.set_facecolor(readout_pal[i % len(readout_pal)])
            box.set_alpha(0.85)
        # one legend handle per readout
        ax_h.bar(0, 0, color=readout_pal[i % len(readout_pal)],
                 edgecolor="black", linewidth=0.3, label=r)
    ax_h.set_xticks(range(len(classes)))
    ax_h.set_xticklabels(classes, rotation=15, ha="right", fontsize=5.6)
    ax_h.set_ylabel("Score (z)")
    ax_h.legend(fontsize=4.8, loc="upper right",
                handletextpad=0.3, borderpad=0.2, labelspacing=0.2)
    style_axes(ax_h)
    panel_tag(ax_h, "h", "BRAF-like vs RAS-like (TCGA-THCA)", dx=-12, dy=4)

    fig_title(fig, 2, "Extended TCGA-THCA validation")
    save(fig, "Supplementary_Fig_2_Extended_TCGA_THCA_validation", OUTDIR)
    print("  ✓ Supp 2 rendered")


# ═══════════════════════════════════════════════════════════════════════════
# Supp 3 — Extended GEO meta-cohort validation
# ═══════════════════════════════════════════════════════════════════════════
def supp3() -> None:
    fig = plt.figure(figsize=(COL2, COL2 * 1.1))
    gs = gridspec.GridSpec(3, 6, figure=fig,
                           hspace=0.85, wspace=1.70,
                           left=0.06, right=0.985, top=0.955, bottom=0.06)

    # (a) GEO bar — diff genes vs FAP
    ax_a = fig.add_subplot(gs[0, :3])
    diff = pd.read_csv(SRC / "Supplementary_Fig3A_GEO_FAP_diffgenes.csv").sort_values("rho")
    cols = [pa["mode"]["PAX8-loss"] if r < 0 else pa["mode"]["SERPINE1-gain"] for r in diff["rho"]]
    ax_a.barh(range(len(diff)), diff["rho"], color=cols, edgecolor="black", linewidth=0.3)
    ax_a.set_yticks(range(len(diff))); ax_a.set_yticklabels(diff["Gene"], fontsize=5.6)
    ax_a.set_xlabel("Spearman ρ vs FAP")
    style_axes(ax_a)
    panel_tag(ax_a, "a", "GEO: differentiation genes vs FAP", dx=-12, dy=4)

    # (b) ECM genes — bar
    ax_b = fig.add_subplot(gs[0, 3:])
    ecm = pd.read_csv(SRC / "Supplementary_Fig3G_GEO_FAP_ECMgenes.csv").sort_values("rho")
    ax_b.barh(range(len(ecm)), ecm["rho"], color=pa["mode"]["SERPINE1-gain"],
              edgecolor="black", linewidth=0.3)
    ax_b.set_yticks(range(len(ecm))); ax_b.set_yticklabels(ecm["Gene"], fontsize=5.6)
    ax_b.set_xlabel("Spearman ρ vs FAP")
    style_axes(ax_b)
    panel_tag(ax_b, "b", "GEO: ECM genes vs FAP", dx=-12, dy=4)

    # (c) ssGSEA: TDS vs FAP_CAF program, colored by group
    ax_c = fig.add_subplot(gs[1, :3])
    ss = pd.read_csv(SRC / "Supplementary_Fig3CF_ssGSEA_scores.csv")
    groups = sorted(ss["Group"].dropna().unique())
    grp_pal = {g: TISSUE_PAL.get(g, pa["scatter_pt"]) for g in groups}
    # fallbacks
    if not all(g in TISSUE_PAL for g in groups):
        cycle = [pa["mode"]["PAX8-loss"], pa["mode"]["Dual-axis"],
                 pa["mode"]["SERPINE1-gain"], pa["integrin"]["avb3"]]
        grp_pal = {g: cycle[i % len(cycle)] for i, g in enumerate(groups)}
    for g in groups:
        sub = ss[ss["Group"] == g]
        ax_c.scatter(sub["FAP_CAF_program"], sub["Differentiation_TDS"],
                     s=10, facecolor=grp_pal[g], edgecolor="black",
                     linewidth=0.3, alpha=0.85, label=g)
    ax_c.set_xlabel("FAP-CAF program (ssGSEA)")
    ax_c.set_ylabel("Differentiation TDS")
    ax_c.legend(fontsize=5.2, loc="upper right",
                handletextpad=0.3, borderpad=0.2, labelspacing=0.25)
    style_axes(ax_c)
    panel_tag(ax_c, "c", "TDS vs FAP-CAF program (GEO meta-cohort)", dx=-12, dy=4)

    # (d) per-cohort FAP — strip by GEO + group
    ax_d = fig.add_subplot(gs[1, 3:])
    pcf = pd.read_csv(SRC / "Supplementary_Fig3D_per_cohort_FAP.csv")
    geos = sorted(pcf["GEO"].dropna().unique())
    for j, g in enumerate(geos):
        sub = pcf[pcf["GEO"] == g]
        for k, grp in enumerate(sorted(sub["Group"].dropna().unique())):
            v = sub.loc[sub["Group"] == grp, "FAP"].values
            x = np.full(len(v), j) + (k - 0.5) * 0.3
            color = grp_pal.get(grp, pa["scatter_pt"])
            ax_d.scatter(x, v, s=8, facecolor=color, edgecolor="black",
                         linewidth=0.2, alpha=0.75,
                         label=grp if j == 0 else None)
    ax_d.set_xticks(range(len(geos)))
    ax_d.set_xticklabels(geos, rotation=20, ha="right", fontsize=5.6)
    ax_d.set_ylabel("FAP (log2 expr)")
    ax_d.legend(fontsize=5, loc="upper right",
                handletextpad=0.3, borderpad=0.2, labelspacing=0.2)
    style_axes(ax_d)
    panel_tag(ax_d, "d", "Per-cohort FAP expression", dx=-12, dy=4)

    # (e) per-cohort FAP–gene Spearman heatmap — transposed (genes on x, wide+short)
    ax_e = fig.add_subplot(gs[2, :3])
    pc = pd.read_csv(SRC / "Supplementary_Fig3H_percohort_FAP_spearman.csv")
    pc_num = pc.set_index("Gene").select_dtypes(include=[np.number]).T
    im = ax_e.imshow(pc_num.values, cmap=pa_cmap_div(), vmin=-1, vmax=1, aspect="auto")
    ax_e.set_xticks(range(pc_num.shape[1]))
    ax_e.set_xticklabels(pc_num.columns, rotation=90, fontsize=4.6)
    ax_e.set_yticks(range(pc_num.shape[0]))
    ax_e.set_yticklabels(pc_num.index, fontsize=6.2)
    cbar_from(ax_e, im, "Spearman ρ", fraction=0.032, label_top=True)
    panel_tag(ax_e, "e", "Per-cohort gene–FAP correlations", dx=-12, dy=4)

    # (f) NT vs ATC ROC
    ax_f = fig.add_subplot(gs[2, 3:])
    roc = pd.read_csv(SRC / "Supplementary_Fig3J_NT_vs_ATC_ROC.csv")
    roc_panel(ax_f, roc)
    panel_tag(ax_f, "f", "NT vs ATC discrimination ROC", dx=-12, dy=4)

    fig_title(fig, 3, "Extended GEO meta-cohort validation")
    save(fig, "Supplementary_Fig_3_Extended_GEO_meta_cohort_validation", OUTDIR)
    print("  ✓ Supp 3 rendered")


# ═══════════════════════════════════════════════════════════════════════════
# Supp 4 — Spatial transcriptomics cohort overview + cluster architecture
# ═══════════════════════════════════════════════════════════════════════════
def supp4() -> None:
    """Five panels of spatial-transcriptomics QC + architecture, orthogonal
    to Fig 2. Old panel d (section FAP vs TDS scatter) dropped — it duplicated
    Fig 2 D (spot-level) and Supp 5 a (per-section) at coarser resolution.
       a — Spot UMAP, Seurat clusters (centroid-labelled, no legend)
       b — Spot UMAP, FAP expression
       c — Section-level cluster fraction
       d — FAP-positive spot fraction vs TDS by tertile
       e — Section-level dominant-class composition
    """
    import matplotlib.patheffects as pe
    fig = plt.figure(figsize=(COL2, COL2 * 1.00))
    gs = gridspec.GridSpec(3, 6, figure=fig,
                           height_ratios=[1.30, 0.85, 0.95],
                           hspace=1.00, wspace=1.20,
                           left=0.06, right=0.96, top=0.945, bottom=0.075)

    sp = pd.read_csv(SRC / "Supplementary_Fig4AB_umap_spot_metadata.csv")
    sp_sub = sp.sample(min(15000, len(sp)), random_state=42)

    # (a) UMAP clusters — centroid number labels instead of cramped legend
    ax_a = fig.add_subplot(gs[0, :3])
    clusters = sorted(sp_sub["seurat_clusters"].dropna().unique())
    cycle = plt.cm.tab20(np.linspace(0, 1, max(len(clusters), 20)))
    for i, c in enumerate(clusters):
        sub = sp_sub[sp_sub["seurat_clusters"] == c]
        ax_a.scatter(sub["UMAP_1"], sub["UMAP_2"], s=1.2, c=[cycle[i]],
                     edgecolor="none", alpha=0.65, rasterized=True)
    for c in clusters:
        sub = sp_sub[sp_sub["seurat_clusters"] == c]
        cx, cy = sub["UMAP_1"].median(), sub["UMAP_2"].median()
        ax_a.text(cx, cy, str(int(c)), fontsize=5.6, fontweight="bold",
                  ha="center", va="center", color="#1A1A1A",
                  path_effects=[pe.withStroke(linewidth=1.7, foreground="white")])
    ax_a.set_xlabel("UMAP 1"); ax_a.set_ylabel("UMAP 2")
    style_axes(ax_a)
    panel_tag(ax_a, "a",
              f"Spot UMAP — Seurat clusters ({len(clusters)} clusters)",
              dx=-12, dy=4)

    # (b) UMAP colored by FAP expression
    ax_b = fig.add_subplot(gs[0, 3:])
    sc = ax_b.scatter(sp_sub["UMAP_1"], sp_sub["UMAP_2"],
                      s=1.2, c=sp_sub["FAP"], cmap="rocket_r",
                      vmin=0, vmax=sp_sub["FAP"].quantile(0.99),
                      edgecolor="none", alpha=0.85, rasterized=True)
    ax_b.set_xlabel("UMAP 1"); ax_b.set_ylabel("UMAP 2")
    cbar_from(ax_b, sc, "FAP (log1p)")
    style_axes(ax_b)
    panel_tag(ax_b, "b", "Spot UMAP — FAP expression", dx=-12, dy=4)

    # (c) section-level cluster fraction — full width, legend outside
    ax_c = fig.add_subplot(gs[1, :])
    cf = pd.read_csv(SRC / "Supplementary_Fig4D_section_cluster_fraction.csv")
    samples = sorted(cf["sample_id"].unique())
    cluster_ids = sorted(cf["seurat_clusters"].unique())
    cluster_pal = plt.cm.tab20(np.linspace(0, 1, max(len(cluster_ids), 20)))
    pivot_c = cf.pivot_table(index="sample_id", columns="seurat_clusters",
                             values="cluster_fraction", fill_value=0).reindex(samples)
    bottoms = np.zeros(len(samples))
    for i, c in enumerate(cluster_ids):
        ax_c.bar(range(len(samples)), pivot_c[c].values, bottom=bottoms,
                 color=cluster_pal[i], edgecolor="white", linewidth=0.3,
                 label=f"c{int(c)}")
        bottoms += pivot_c[c].values
    ax_c.set_xticks(range(len(samples)))
    ax_c.set_xticklabels(samples, rotation=45, ha="right", fontsize=5.8)
    ax_c.set_ylim(0, 1)
    ax_c.set_ylabel("Cluster fraction")
    ax_c.legend(fontsize=4.8, loc="center left", bbox_to_anchor=(1.005, 0.5),
                ncols=2, handletextpad=0.3, borderpad=0.3, labelspacing=0.25,
                title="Cluster", title_fontsize=5.2)
    style_axes(ax_c)
    panel_tag(ax_c, "c", "Section-level cluster fraction", dx=-12, dy=4)

    # (d) FAP-positive spot fraction vs TDS by tertile
    ax_d = fig.add_subplot(gs[2, :3])
    tt = pd.read_csv(SRC / "Supplementary_Fig4H_per_section_FAP_TDS_tertile.csv")
    tert_pal = {"low": pa["mode"]["PAX8-loss"], "mid": pa["mode"]["Dual-axis"],
                "high": pa["mode"]["SERPINE1-gain"]}
    for tert, sub in tt.groupby("TDS_tertile"):
        col = tert_pal.get(str(tert).lower(), pa["scatter_pt"])
        ax_d.scatter(sub["fap_pos_frac"], sub["median_TDS"],
                     s=24, facecolor=col, edgecolor="black", linewidth=0.3,
                     label=f"TDS {tert}")
    ax_d.set_xlabel("FAP-positive spot fraction")
    ax_d.set_ylabel("Median TDS")
    ax_d.legend(fontsize=5.2, loc="upper right",
                handletextpad=0.3, borderpad=0.2, labelspacing=0.25)
    style_axes(ax_d)
    panel_tag(ax_d, "d", "FAP-positive fraction vs TDS by tertile", dx=-12, dy=4)

    # (e) section-level dominant-class composition
    ax_e = fig.add_subplot(gs[2, 3:])
    mc = pd.read_csv(SRC / "Supplementary_Fig4G_section_module_composition.csv")
    mc = mc.dropna(subset=["dominant_class"])
    cls = sorted(mc["dominant_class"].astype(str).unique())
    samples_e = sorted(mc["sample_id"].unique())
    pv = mc.pivot_table(index="sample_id", columns="dominant_class",
                        values="spot_fraction", fill_value=0).reindex(samples_e)
    cls_pal = {"identity": pa["mode"]["PAX8-loss"],
               "intermediate": pa["mode"]["Dual-axis"],
               "invasion": pa["mode"]["SERPINE1-gain"]}
    bottoms = np.zeros(len(samples_e))
    for c in cls:
        col = cls_pal.get(c, "#888888")
        ax_e.bar(range(len(samples_e)), pv[c].values, bottom=bottoms,
                 color=col, edgecolor="white", linewidth=0.3, label=c)
        bottoms += pv[c].values
    ax_e.set_xticks(range(len(samples_e)))
    ax_e.set_xticklabels(samples_e, rotation=45, ha="right", fontsize=5.8)
    ax_e.set_ylim(0, 1)
    ax_e.set_ylabel("Spot fraction")
    ax_e.legend(fontsize=5.2, loc="center left", bbox_to_anchor=(1.005, 0.5),
                handletextpad=0.3, borderpad=0.2, labelspacing=0.25)
    style_axes(ax_e)
    panel_tag(ax_e, "e", "Section-level dominant-class composition", dx=-12, dy=4)

    fig_title(fig, 4, "Spatial transcriptomics cohort overview and marker-guided spot-cluster architecture")
    save(fig, "Supplementary_Fig_4_Spatial_transcriptomics_cohort_overview_and_marker_guided_spot_cluster_architecture", OUTDIR)
    print("  ✓ Supp 4 rendered (5-panel, section FAP-vs-TDS dropped)")


# ═══════════════════════════════════════════════════════════════════════════
# Supp 5 — Section-level spatial support for Fig 2 (spatial transcriptomics)
# ═══════════════════════════════════════════════════════════════════════════
def supp5() -> None:
    """Section-level spatial support for Fig 2 (orthogonal to the main figure):
    mean FAP (log2 CP10K) across the Visium sections. The tissue-level
    experimental (wet-lab) panels are out of scope for this bioinformatics-only
    repository and have been removed.
    """
    fig, ax_a = plt.subplots(figsize=(COL1, COL1 * 0.72))

    # section-level mean FAP — horizontal bar
    ss = pd.read_csv(SRC / "Supplementary_Fig5C_section_level_spatial_summary.csv").sort_values("mean_FAP_logCP10K")
    ax_a.barh(range(len(ss)), ss["mean_FAP_logCP10K"], color=pa["scatter_pt"],
              edgecolor="black", linewidth=0.3)
    ax_a.set_yticks(range(len(ss))); ax_a.set_yticklabels(ss["slice"], fontsize=6.0)
    ax_a.set_xlabel(r"Mean FAP (log$_2$CP10K)")
    ax_a.set_ylabel("Spatial section")
    style_axes(ax_a)

    fig_title(fig, 5, "Section-level mean FAP across the spatial transcriptomics cohort (support for Fig. 2)")
    save(fig, "Supplementary_Fig_5_Section_level_mean_FAP_across_the_spatial_cohort_support_for_Fig2", OUTDIR)
    print("  ✓ Supp 5 rendered (1-panel, spatial only; experimental panels removed)")


# ═══════════════════════════════════════════════════════════════════════════
# Supp 6 — Mouse scRNA + bulk validation (Seurat-only; uses R bridge)
# ═══════════════════════════════════════════════════════════════════════════
def supp6() -> None:
    """Supp 6 — Source overview + extended transcriptomic validation for
    the time-resolved mouse model. Reads the per-panel CSVs that the R
    builder (scripts/fig5_mouse_braf/02_supp_expanded.R) already exports to
    `outputs/supplementary_figures/source_data/supp_fig6/`.
    """
    SUB = SRC / "supp_fig6"
    tp_order = ["1 month", "2 months", "4 months"]
    tp_pal = {"1 month": "#E78B8B", "2 months": "#86C06B", "4 months": "#6BAED6"}
    lineage_pal = {"Epithelial": pa["mode"]["PAX8-loss"],
                   "Fibroblast": pa["mode"]["SERPINE1-gain"],
                   "Myeloid":    "#D55E00",
                   "Endothelial": pa["integrin"]["avb3"],
                   "T_cell":     pa["mode"]["Dual-axis"]}

    fig = plt.figure(figsize=(COL2, COL2 * 1.30))
    gs = gridspec.GridSpec(4, 6, figure=fig,
                           hspace=0.95, wspace=1.10,
                           left=0.06, right=0.985, top=0.955, bottom=0.05)

    # (a) mouse scRNA timepoint overview (bar)
    ax_a = fig.add_subplot(gs[0, :2])
    ov = pd.read_csv(SUB / "supp_fig6_panelA_mouse_scRNA_overview.csv")
    ov["Time point"] = pd.Categorical(ov["Time point"], categories=tp_order, ordered=True)
    ov = ov.sort_values("Time point")
    xpos = np.arange(len(ov))
    ax_a.bar(xpos - 0.18, ov["All cells"], width=0.32,
             color=pa["scatter_pt"], edgecolor="black", linewidth=0.3, label="All cells")
    ax_a.bar(xpos + 0.18, ov["Fibroblast/CAF cells"], width=0.32,
             color=pa["mode"]["SERPINE1-gain"], edgecolor="black", linewidth=0.3,
             label="Fibroblast/CAF")
    ax_a.set_xticks(xpos); ax_a.set_xticklabels(ov["Time point"], rotation=15, ha="right", fontsize=5.6)
    ax_a.set_ylabel("n cells")
    ax_a.legend(fontsize=5.2, loc="upper left",
                handletextpad=0.3, borderpad=0.2)
    style_axes(ax_a)
    panel_tag(ax_a, "a", "Mouse scRNA timepoint overview", dx=-10, dy=4)

    # (b) bulk cohort overview (table)
    ax_b = fig.add_subplot(gs[0, 2:])
    ax_b.axis("off")
    bk = pd.read_csv(SUB / "supp_fig6_panelB_bulk_mTC_overview.csv")
    tbl = ax_b.table(cellText=bk.values.tolist(),
                     colLabels=list(bk.columns),
                     loc="center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(5.2); tbl.scale(1, 1.15)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.3); cell.set_edgecolor("#9A9A9A")
        if r == 0:
            cell.set_facecolor("#EAEAEA"); cell.set_text_props(weight="bold")
    panel_tag(ax_b, "b", "Bulk mouse-TC cohorts", dx=-6, dy=4)

    # (c) GSE30427 radar — gene × WT/FTC/ATC
    ax_c = fig.add_subplot(gs[1, :3], projection="polar")
    rd = pd.read_csv(SUB / "supp_fig6_panelD_gse30427_radar_scaled.csv")
    genes = rd["gene"].tolist()
    grp_pal = {"WT": pa["tissue_human"]["NT"], "FTC": pa["mode"]["Dual-axis"],
               "ATC": pa["mode"]["SERPINE1-gain"]}
    angles = np.linspace(0, 2*np.pi, len(genes), endpoint=False).tolist()
    angles_closed = angles + [angles[0]]
    for g in ["WT", "FTC", "ATC"]:
        vals = rd[g].tolist() + [rd[g].iloc[0]]
        ax_c.plot(angles_closed, vals, color=grp_pal[g], linewidth=1.0, label=g)
        ax_c.fill(angles_closed, vals, color=grp_pal[g], alpha=0.18)
    ax_c.set_xticks(angles); ax_c.set_xticklabels(genes, fontsize=4.8)
    ax_c.set_ylim(0, 105)
    ax_c.tick_params(axis="y", labelsize=5)
    ax_c.legend(fontsize=5.2, loc="lower left", bbox_to_anchor=(-0.18, -0.05),
                handletextpad=0.3, borderpad=0.2)
    ax_c.set_title("GSE30427 differentiation gene radar (scaled 0-100)",
                   loc="left", fontsize=7, pad=6)
    # panel_tag on polar — use figure coords
    fig.text(gs[1, :3].get_position(fig).x0 - 0.005,
             gs[1, :3].get_position(fig).y1 - 0.005,
             "c", fontsize=9, fontweight="bold")

    # (d) QC by timepoint × metric (boxplot-like points across timepoints)
    ax_d = fig.add_subplot(gs[1, 3:])
    qc = pd.read_csv(SUB / "supp_fig6_panelE_qc_median.csv")
    qc["Time"] = pd.Categorical(qc["Time"], categories=tp_order, ordered=True)
    metrics_q = sorted(qc["metric"].unique())
    n_q = len(metrics_q)
    width_q = 0.8 / max(n_q, 1)
    metric_pal = [pa["mode"]["PAX8-loss"], pa["mode"]["SERPINE1-gain"],
                  pa["mode"]["Dual-axis"], pa["integrin"]["avb3"]]
    for i, m in enumerate(metrics_q):
        sub = qc[qc["metric"] == m].sort_values("Time")
        x = np.arange(len(sub)) + (i - (n_q-1)/2) * width_q
        ax_d.bar(x, sub["median"], width=width_q*0.9,
                 color=metric_pal[i % len(metric_pal)],
                 edgecolor="black", linewidth=0.3, label=m)
    ax_d.set_xticks(range(len(tp_order))); ax_d.set_xticklabels(tp_order, rotation=15, ha="right", fontsize=5.6)
    ax_d.set_ylabel("Median QC value")
    ax_d.legend(fontsize=4.6, loc="upper right",
                handletextpad=0.3, borderpad=0.2, labelspacing=0.2)
    style_axes(ax_d)
    panel_tag(ax_d, "d", "Per-timepoint QC medians", dx=-10, dy=4)

    # (e) lineage marker dotplot
    ax_e = fig.add_subplot(gs[2, :3])
    lin = pd.read_csv(SUB / "supp_fig6_panelF_lineage_dotplot.csv")
    feats = list(pd.Index(lin["features.plot"]).unique())
    ids = list(pd.Index(lin["id"]).unique())
    Xc = {f: i for i, f in enumerate(feats)}; Yc = {g: j for j, g in enumerate(ids)}
    sc = ax_e.scatter(lin["features.plot"].map(Xc), lin["id"].map(Yc),
                      s=(lin["pct.exp"]/lin["pct.exp"].max())*55 + 4,
                      c=lin["avg.exp.scaled"], cmap=pa_cmap_div(),
                      vmin=-2, vmax=2, edgecolor="black", linewidth=0.2)
    ax_e.set_xticks(range(len(feats))); ax_e.set_xticklabels(feats, rotation=45, ha="right", fontsize=5.4)
    ax_e.set_yticks(range(len(ids))); ax_e.set_yticklabels(ids, fontsize=5.6)
    cbar_from(ax_e, sc, "scaled avg expr")
    panel_tag(ax_e, "e", "Mouse lineage marker dot-plot", dx=-12, dy=4)

    # (f) CAF subtype dot-plot
    ax_f = fig.add_subplot(gs[2, 3:])
    cdp = pd.read_csv(SUB / "supp_fig6_panelG_CAF_dotplot.csv")
    feats2 = list(pd.Index(cdp["features.plot"]).unique())
    ids2 = list(pd.Index(cdp["id"]).unique())
    Xc2 = {f: i for i, f in enumerate(feats2)}; Yc2 = {g: j for j, g in enumerate(ids2)}
    sc2 = ax_f.scatter(cdp["features.plot"].map(Xc2), cdp["id"].map(Yc2),
                       s=(cdp["pct.exp"]/cdp["pct.exp"].max())*55 + 4,
                       c=cdp["avg.exp.scaled"], cmap=pa_cmap_div(),
                       vmin=-2, vmax=2, edgecolor="black", linewidth=0.2)
    ax_f.set_xticks(range(len(feats2))); ax_f.set_xticklabels(feats2, rotation=45, ha="right", fontsize=5.4)
    ax_f.set_yticks(range(len(ids2))); ax_f.set_yticklabels(ids2, fontsize=5.6)
    cbar_from(ax_f, sc2, "scaled avg expr")
    panel_tag(ax_f, "f", "Mouse CAF subtype dot-plot", dx=-12, dy=4)

    # (g) cell composition by lineage × timepoint (stacked bar)
    ax_g = fig.add_subplot(gs[3, :3])
    cc = pd.read_csv(SUB / "supp_fig6_panelI_cell_composition.csv")
    cc["Time"] = pd.Categorical(cc["Time"], categories=tp_order, ordered=True)
    lineages = sorted(cc["Lineage"].unique())
    M = cc.pivot_table(index="Time", columns="Lineage",
                       values="fraction", fill_value=0)
    M = M.reindex(tp_order)
    bottoms = np.zeros(len(tp_order))
    for lin_ in lineages:
        col = lineage_pal.get(lin_, "#999999")
        ax_g.bar(range(len(tp_order)), M[lin_].values, bottom=bottoms,
                 color=col, edgecolor="white", linewidth=0.3, label=lin_)
        bottoms += M[lin_].values
    ax_g.set_xticks(range(len(tp_order))); ax_g.set_xticklabels(tp_order, rotation=15, ha="right", fontsize=5.6)
    ax_g.set_ylabel("Fraction")
    ax_g.set_ylim(0, 1)
    ax_g.legend(fontsize=5.2, loc="center left", bbox_to_anchor=(1.01, 0.5),
                handletextpad=0.3, borderpad=0.2, labelspacing=0.2)
    style_axes(ax_g)
    panel_tag(ax_g, "g", "Cell composition × timepoint", dx=-12, dy=4)

    # (h) bulk ssGSEA — GSE30427 module × group
    ax_h = fig.add_subplot(gs[3, 3:])
    ks = pd.read_csv(SUB / "supp_fig6_panelK_GSE30427_ssGSEA.csv")
    groups_k = sorted(ks["Group"].dropna().unique())
    grp_pal_k = {g: tp_pal.get(g, pa["scatter_pt"]) for g in groups_k}
    # fallback Hi-contrast palette
    if not any(g in tp_pal for g in groups_k):
        cycle = [pa["tissue_human"]["NT"], pa["tissue_human"]["PTC"],
                 pa["tissue_human"]["ATC"], pa["mode"]["Dual-axis"]]
        grp_pal_k = {g: cycle[i % len(cycle)] for i, g in enumerate(groups_k)}
    modules = sorted(ks["Module"].unique())
    n_m = len(modules); width_h = 0.8 / max(n_m, 1)
    for i, m in enumerate(modules):
        for j, g in enumerate(groups_k):
            vals = ks.loc[(ks["Module"] == m) & (ks["Group"] == g), "Score"].dropna().values
            if len(vals) == 0:
                continue
            x_center = j + (i - (n_m-1)/2) * width_h
            ax_h.bar(x_center, vals.mean(), width=width_h*0.85,
                     yerr=vals.std(ddof=1)/np.sqrt(len(vals)) if len(vals) > 1 else 0,
                     capsize=1.2,
                     color=grp_pal_k.get(g, pa["scatter_pt"]),
                     edgecolor="black", linewidth=0.3,
                     error_kw=dict(linewidth=0.4),
                     label=m if j == 0 else None)
    ax_h.set_xticks(range(len(groups_k))); ax_h.set_xticklabels(groups_k, rotation=15, ha="right", fontsize=5.6)
    ax_h.set_ylabel("ssGSEA score")
    ax_h.legend(fontsize=4.6, loc="upper right",
                handletextpad=0.3, borderpad=0.2, labelspacing=0.2,
                ncols=2)
    style_axes(ax_h)
    panel_tag(ax_h, "h", "GSE30427 module ssGSEA by group", dx=-12, dy=4)

    fig_title(fig, 6,
              "Source overview and extended transcriptomic validation for the time-resolved mouse model")
    save(fig, "Supplementary_Fig_6_Source_overview_and_extended_transcriptomic_validation_for_the_time_resolved_mouse_model",
         OUTDIR)
    print("  ✓ Supp 6 rendered (real data, 8 panels)")


# ═══════════════════════════════════════════════════════════════════════════
SUPP_FUNCS = {
    "supp1": supp1, "supp2": supp2, "supp3": supp3,
    "supp4": supp4, "supp5": supp5, "supp6": supp6,
}


def main() -> None:
    args = sys.argv[1:]
    targets = list(SUPP_FUNCS) if not args else args
    print(f"paper_A matplotlib render → {OUTDIR.relative_to(ROOT)}")
    for name in targets:
        fn = SUPP_FUNCS.get(name)
        if fn is None:
            print(f"  ✗ unknown: {name}")
            continue
        print(f"  • {name}")
        try:
            fn()
        except Exception as e:
            import traceback
            print(f"  ✗ {name} FAILED: {e}")
            traceback.print_exc()


if __name__ == "__main__":
    main()
