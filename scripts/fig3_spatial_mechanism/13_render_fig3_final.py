#!/usr/bin/env python3
# Cancer Research submission - figure code release
# Builds: Figure 3, final submitted version. Entry point for the figure: it
#          monkeypatches panel A to color couplings by receiver and to italicize
#          the non-integrin receivers (CD44, SDC1, DDR1), matching the legend,
#          then runs the full render and all beautify passes.
#          Pipeline: build_fig3_panelA_paired -> fig3_full_render ->
#          fig3_beautify -> this script.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

"""Rebuild the whole Figure 3 with panel A's couplings coloured by receiver.

Monkeypatches two things on the existing pipeline; no source file is edited.
  1. full.draw_panel_A  -> per-receptor edge and label colours
  2. split_ecm_receptors -> italicise non-integrin receptors instead of
     recolouring them (the recolour would overwrite the new encoding)
"""
import os
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.text import Text

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "_shared"))

import build_fig3_panelA_paired as A
import fig3_full_render as full
import fig3_beautify as B

# Okabe-Ito, safe under deuteranopia/protanopia
BY_RECEPTOR = {
    "ITGA2": "#D55E00", "ITGB1": "#0072B2", "ITGA3": "#009E73",
    "CD44":  "#CC79A7", "SDC1":  "#E69F00", "ITGA5": "#9A9A9A",
    "DDR1":  "#6E6E6E", "CXCR4": "#0B7A75", "TGFBR2": "#B8860B",
}
SENDER_INK = "#4A5568"
NONINTEGRIN = {"CD44", "SDC1", "DDR1"}


def _labels(ax, genes, groups, cmap, italic=()):
    y = np.arange(len(genes))[::-1]
    ax.set_xlim(0, 1); ax.set_ylim(-0.5, len(genes) - 0.5)
    prev, start = None, None
    for i, g in enumerate(groups):
        if g != prev:
            if prev is not None:
                ax.axhspan(y[start] + 0.5, y[i-1] - 0.5,
                           color=A.MAINT_BAND if prev == "maintenance" else A.REWIR_BAND, zorder=0)
            start, prev = i, g
    if prev is not None:
        ax.axhspan(y[start] + 0.5, y[len(genes)-1] - 0.5,
                   color=A.MAINT_BAND if prev == "maintenance" else A.REWIR_BAND, zorder=0)
    for yi, g in zip(y, genes):
        ax.text(0.5, yi, g, ha="center", va="center", fontsize=8.5, fontweight="bold",
                style="italic" if g in italic else "normal",
                color=cmap.get(g, SENDER_INK), zorder=4)
    ax.set_xticks([]); ax.set_yticks([]); full.hide_spines(ax)


def draw_panel_A(fig, gs):
    """Per-receptor coupling colours. Same data, same layout as the original."""
    sender_df = pd.read_csv(A.SENDER, sep="\t")
    recv_df   = pd.read_csv(A.RECV,   sep="\t")
    liana_df  = pd.read_csv(A.LIANA,  sep="\t")
    sg = A.SENDER_MAINT + A.SENDER_REWIR
    rg = A.RECV_MAINT   + A.RECV_REWIR
    le = A.parse_liana_edges(liana_df, sg, rg, top_n=12)
    lset = {(l, r) for (l, r, _, _) in le}
    cur = [(l, r, pw, 0.30) for (l, r, pw) in A.CURATED_EDGES
           if (l, r) not in lset and l in sg and r in rg]
    edges = ([(l, r, pw, s, "LIANA") for (l, r, pw, s) in le]
           + [(l, r, pw, s, "curated") for (l, r, pw, s) in cur])

    gA = gridspec.GridSpecFromSubplotSpec(
        1, 5, subplot_spec=gs, width_ratios=[1.1, 1.15, 2.35, 1.15, 1.1], wspace=0.05)

    # sender bars (mirrored)
    ax = fig.add_subplot(gA[0])
    s = sender_df.set_index("gene").reindex(sg).reset_index()
    y = np.arange(len(sg))[::-1]
    dv, pv, pc = s["delta_log_norm"].values, s["p_adj_BH"].values, s["pct_grp1"].values
    xl = max(abs(dv)) * 1.4
    for i, (d, p, c) in enumerate(zip(dv, pv, pc)):
        ax.barh(y[i], -d, color=A.UP_COLOR if d >= 0 else A.DOWN_COLOR,
                edgecolor="black", linewidth=0.35, height=0.55, zorder=2)
        st = A.sig_stars(p)
        if st:
            ax.text(-d - 0.08, y[i], st, fontsize=6.5, fontweight="bold",
                    va="center", ha="right", color="#333", zorder=4)
        ax.scatter(-xl*1.03, y[i], s=35 + 230*c, facecolor="#E2A736",
                   edgecolor="#333", linewidth=0.35, clip_on=False, zorder=3)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlim(-xl*1.15, 0.25); ax.set_ylim(-0.5, len(sg)-0.5)
    ax.set_yticks([]); ax.invert_xaxis(); ax.tick_params(axis="x", labelsize=6.2)
    full.hide_spines(ax, keep=("bottom",)); ax.spines["bottom"].set_linewidth(0.4)
    ax.set_xlabel("Δ log-norm (FAP-CAF − other)", fontsize=6.3)
    ax.set_title("Sender", fontsize=9.5, fontweight="bold", loc="right",
                 color="#B42318", pad=3)

    _labels(fig.add_subplot(gA[1]), sg,
            ["maintenance"]*len(A.SENDER_MAINT) + ["rewiring"]*len(A.SENDER_REWIR),
            {g: SENDER_INK for g in sg})

    axm = fig.add_subplot(gA[2])
    axm.set_xlim(0, 1); axm.set_ylim(-0.5, max(len(sg), len(rg)) - 0.5); axm.axis("off")
    ys = {g: (len(sg)-1)-i for i, g in enumerate(sg)}
    yr = {g: (len(rg)-1)-i for i, g in enumerate(rg)}
    for (lig, rec, pw, st, src) in sorted(edges, key=lambda e: 0 if e[4] == "curated" else 1):
        if lig not in ys or rec not in yr:
            continue
        A.draw_connection(axm, ys[lig], yr[rec], BY_RECEPTOR.get(rec, "#4A5568"), st,
                          x_left=0.01, x_right=0.99,
                          linestyle="solid" if src == "LIANA" else "dashed")
    axm.text(0.5, -1.15,
             "solid = LIANA top-12  ·  dashed = curated-only  ·  "
             "thickness: −log10(rank)  ·  colour: receiver",
             ha="center", va="top", fontsize=6.0, style="italic", color="#666",
             transform=axm.transData, clip_on=False)

    _labels(fig.add_subplot(gA[3]), rg,
            ["maintenance"]*len(A.RECV_MAINT) + ["rewiring"]*len(A.RECV_REWIR),
            BY_RECEPTOR, italic=NONINTEGRIN)

    ax = fig.add_subplot(gA[4])
    r = recv_df.set_index("gene").reindex(rg).reset_index()
    y = np.arange(len(rg))[::-1]
    dv, pv, pc = r["delta_log_norm"].values, r["p_adj_BH"].values, r["pct_grp1"].values
    xl = max(abs(dv)) * 1.4
    for i, (d, p, c) in enumerate(zip(dv, pv, pc)):
        ax.barh(y[i], d, color=A.UP_COLOR if d >= 0 else A.DOWN_COLOR,
                edgecolor="black", linewidth=0.35, height=0.55, zorder=2)
        st = A.sig_stars(p)
        if st:
            off = 0.05 if d >= 0 else -0.05
            ax.text(d + off, y[i], st, fontsize=6.5, fontweight="bold", va="center",
                    ha="left" if d >= 0 else "right", color="#333", zorder=4)
        ax.scatter(xl*1.03, y[i], s=35 + 230*c, facecolor="#E2A736",
                   edgecolor="#333", linewidth=0.35, clip_on=False, zorder=3)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlim(-0.25, xl*1.15); ax.set_ylim(-0.5, len(rg)-0.5)
    ax.set_yticks([]); ax.tick_params(axis="x", labelsize=6.2)
    full.hide_spines(ax, keep=("bottom",)); ax.spines["bottom"].set_linewidth(0.4)
    ax.set_xlabel("Δ log-norm (terminal − lineage)", fontsize=6.3)
    ax.set_title("Receiver", fontsize=9.5, fontweight="bold", loc="left",
                 color="#C0392B", pad=3)


def italicise_nonintegrin(fig):
    n = 0
    for t in fig.findobj(Text):
        if t.get_text().strip() in B.NONINTEGRIN_RECEPTORS:
            t.set_style("italic"); n += 1
    return n



# ---- panel H: draw the per-section spread instead of writing it out ----
def draw_panel_H_persection(fig, gs, full_module):
    """One-line annotation; the 8 per-section fits carry the range visually.

    The pooled coefficient is descriptive (spots are nested within sections), so it
    is reported in the legend rather than given equal typographic weight on the panel.
    """
    import pandas as pd
    from scipy import stats
    from matplotlib.colors import LinearSegmentedColormap
    req = ["fap_niche_proximity_index", "collagen_sender_program_activity",
           "itga2_receiver_program_activity", "dediff_shift_index"]
    d = pd.read_csv(full_module.AUG_CTX, sep="\t").dropna(subset=req).copy()
    d = d[d["epi_valid_primary"] == True].copy()

    def s01(v):
        a = np.asarray(v, float)
        return (a - np.nanmin(a)) / max(np.nanmax(a) - np.nanmin(a), 1e-9)

    d["idx3"] = (s01(d[req[0]]) + s01(d[req[1]]) + s01(d[req[2]])) / 3.0
    x, y = d["idx3"].to_numpy(float), d[req[3]].to_numpy(float)
    rho_pool, _ = stats.spearmanr(x, y)
    secs = [(s, dd) for s, dd in d.groupby("sample") if len(dd) >= 30]
    per = np.asarray([stats.spearmanr(dd["idx3"], dd[req[3]])[0] for _, dd in secs])

    ax = fig.add_subplot(gs)
    cmap = LinearSegmentedColormap.from_list(
        "density_neutral_exact", ["#F4F6FA", "#B0C4D9", "#4F6C8F", "#1A2A44"])
    hb = ax.hexbin(x, y, gridsize=42, cmap=cmap, bins="log", mincnt=1,
                   linewidths=0, rasterized=True)
    for _, dd in secs:                                   # per-section fits
        sx = dd["idx3"].to_numpy(float); sy = dd[req[3]].to_numpy(float)
        s_, i_ = np.polyfit(sx, sy, 1)
        xs = np.linspace(np.nanmin(sx), np.nanmax(sx), 50)
        ax.plot(xs, s_ * xs + i_, color="#555555", lw=0.6, alpha=0.75, zorder=3)
    sl, ic = np.polyfit(x, y, 1)
    xx = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    ax.plot(xx, sl * xx + ic, color="#1A1A1A", lw=1.6, zorder=4)
    ax.text(0.03, 0.97,
            f"per-section $\\rho$ = {np.median(per):+.3f} ($n$ = {len(per)})",
            transform=ax.transAxes, fontsize=7.0, va="top", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.26", fc="white", ec="#888", lw=0.4))
    ax.axhline(0, color="#888", lw=0.4)
    ax.set_xlabel("ECM-integrin support", fontsize=7.0)
    ax.set_ylabel("Dediff shift", fontsize=7.0)
    full_module.hide_spines(ax, keep=("left", "bottom"))
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.5)
    cb = plt.colorbar(hb, ax=ax, fraction=0.048, pad=0.02)
    cb.set_label("log10 spot count", fontsize=6.2)
    cb.ax.tick_params(labelsize=5.6)
    return len(d), float(np.median(per)), float(per.min()), float(per.max())



# ---- panel B: four independent measurements, real deltas, equal-weight integration ----
_B_SPEC = [("ITGA2", "ECM/integrin", "ITGA2"), ("ITGB1", "ECM/integrin", "ITGB1"),
           ("CXCL12/\nCXCR4", "CXCL12-CXCR4", "CXCR4")]
_B_COLS = ["cell_delta_vs_noncore", "spatial_delta_seed_vs_other",
           "cell_delta_terminal_vs_lineage", "spatial_delta_terminal_vs_lineage"]
_B_ROWS = ["CAF\n· cell", "CAF\n· spatial", "Epi\n· cell", "Epi\n· spatial"]


def draw_panel_B_support(fig, gs):
    """Four measurements in Δ log-norm, plus a physically separated integrated index.

    Rows 1–4 print the real mean difference; fill is rescaled within each row so the
    leading receptor is 1. The integrated row is a unitless 0–1 index (unweighted mean
    of the four row-rescaled values) and is set apart so its different unit is visible.
    """
    from matplotlib.colors import LinearSegmentedColormap
    r = pd.read_csv(full.FIBRO_P.parent / "ranked_candidate_interactions.tsv", sep="\t")
    M = np.array([[r[(r.pathway == p) & (r.receptor == rc)][c].mean() for c in _B_COLS]
                  for _, p, rc in _B_SPEC])
    nint = [len(r[(r.pathway == p) & (r.receptor == rc)]) for _, p, rc in _B_SPEC]
    norm = M / M.max(axis=0)
    integ = norm.mean(axis=1)
    cmap = LinearSegmentedColormap.from_list("pa_seq_B", ["#FDF3EF", "#E8B4AD", "#C0392B"])
    xt = [f"{d}\n({k})" for (d, _, _), k in zip(_B_SPEC, nint)]

    gB = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs,
                                          height_ratios=[4, 1], hspace=0.42)
    ax = fig.add_subplot(gB[0])
    ax.imshow(norm.T, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks([]); ax.set_yticks(range(4))
    ax.set_yticklabels(_B_ROWS, fontsize=6.0, linespacing=1.25)
    for i in range(4):
        for j in range(3):
            ax.text(j, i, f"{M[j, i]:.2f}", ha="center", va="center", fontsize=6.3,
                    fontweight="bold", color="white" if norm[j, i] > 0.62 else "#3A2018")
    for xx in np.arange(0.5, 2.5): ax.axvline(xx, color="white", lw=0.9)
    for yy in np.arange(0.5, 3.5): ax.axhline(yy, color="white", lw=0.9)
    ax.set_ylabel("Δ log-norm", fontsize=6.0, labelpad=2)
    for sp in ax.spines.values(): sp.set_linewidth(0.45)
    ax.tick_params(length=0, pad=2)

    ax2 = fig.add_subplot(gB[1])
    ax2.imshow(integ.reshape(1, -1), cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax2.set_xticks(range(3)); ax2.set_xticklabels(xt, fontsize=6.0)
    ax2.set_yticks([0]); ax2.set_yticklabels(["Integrated\n(0–1)"], fontsize=6.0, linespacing=1.25)
    for j in range(3):
        ax2.text(j, 0, f"{integ[j]:.2f}", ha="center", va="center", fontsize=6.3,
                 fontweight="bold", color="white" if integ[j] > 0.62 else "#3A2018")
    for xx in np.arange(0.5, 2.5): ax2.axvline(xx, color="white", lw=0.9)
    for sp in ax2.spines.values(): sp.set_linewidth(0.45)
    ax2.tick_params(length=0, pad=2)


full.draw_panel_A = draw_panel_A
B.full = full
B.split_ecm_receptors = italicise_nonintegrin
B.draw_panel_H_persection = draw_panel_H_persection
B.moderate.draw_panel_B_support = draw_panel_B_support

OUTDIR = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>")) / "outputs" / "fig3"
OUTDIR.mkdir(parents=True, exist_ok=True)
out = OUTDIR / "Figure3_paperA_recolored"
n_spots, rmed = B.build(out)
print(f"\nrebuilt -> {out}.pdf/.png   spots={n_spots:,}  per-section rho={rmed:+.3f}")
