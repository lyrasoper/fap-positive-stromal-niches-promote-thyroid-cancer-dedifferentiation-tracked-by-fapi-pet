# Cancer Research submission - figure code release
# Builds: Figure 3 layout variant. Imported by fig3_beautify.py.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

"""
Moderately simplified NC main-text Fig. 7.

Compared with the full atlas version, this keeps the A-H evidence chain but
compresses repeated or atlas-like content:

  A  ITGA2-centred stromal-epithelial mechanism model
  B  pathway-level support
  C  LIANA consensus support
  D  downstream program convergence
  E  representative spatial maps (5 samples x 5 evidence rows)
  F  distance-dependent epithelial-state redistribution
  G  per-sample paired state shift
  H  spot-level ECM-integrin/dedifferentiation closure

The full 8-sample x 7-row spatial atlas and robustness panels should remain in
the supplement. This script does not overwrite any existing final figure.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats

import fig3_full_render as full


OUT = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>")) / "outputs" / "fig7_redesign"
AUG_CTX = OUT / "augmented_context_v3.tsv"
REP_SAMPLES = ["P12", "P17", "P26", "P44", "P83"]


plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.size": 7.6,
    "axes.labelsize": 7.4,
    "axes.titlesize": 8.8,
    "xtick.labelsize": 6.7,
    "ytick.labelsize": 6.7,
    "legend.fontsize": 6.7,
    "axes.linewidth": 0.58,
    "xtick.major.width": 0.52,
    "ytick.major.width": 0.52,
    "xtick.major.size": 2.3,
    "ytick.major.size": 2.3,
    "lines.linewidth": 1.0,
})


CMAP_SEQ = LinearSegmentedColormap.from_list(   # refined muted warm sequential 2026-06-02 (editorial, less candy-red)
    "nc_seq", ["#FBF5F1", "#EFD3C4", "#DBA088", "#C06F58", "#9C3D2F", "#611812"]
)
CMAP_ECM = LinearSegmentedColormap.from_list(
    "ecm", ["#FAF3F0", "#F1C27D", "#C07A39", "#6C3A16"]
)
CMAP_SUPPORT = LinearSegmentedColormap.from_list(   # same refined muted warm ramp as CMAP_SEQ (consistency)
    "support_red", ["#FBF5F1", "#EFD3C4", "#DBA088", "#C06F58", "#9C3D2F", "#611812"]
)
CMAP_INTEG = LinearSegmentedColormap.from_list(
    "integrin", ["#F6EEF3", "#C99BC3", "#7E3783", "#31114A"]
)
CMAP_LINEAGE = LinearSegmentedColormap.from_list(
    "lineage_teal", ["#F3F7F6", "#9ED3C8", "#2C7F62", "#14433A"]
)
CMAP_TERMINAL = LinearSegmentedColormap.from_list(
    "terminal_red", ["#FCE9E5", "#F29A7F", "#C93B26", "#3A0C09"]
)


def scale01(x):
    x = np.asarray(x, dtype=float)
    lo, hi = np.nanmin(x), np.nanmax(x)
    return (x - lo) / max(hi - lo, 1e-9)


def add_panel_label(fig, x, y, letter):
    fig.text(x, y, letter, fontsize=12, fontweight="bold", ha="left", va="top")


def add_scalebar(ax, length=0.18):
    x0, y0 = 0.045, 0.055
    ax.plot([x0, x0 + length], [y0, y0], color="black", lw=1.25,
            transform=ax.transAxes, solid_capstyle="butt")


def draw_panel_A_schematic(fig, gs):
    ax = fig.add_subplot(gs)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    accent = "#C65A43"
    accent_light = "#F1C8B8"
    neutral = "#5A5A5A"
    pale = "#EEF1F4"
    thyroid = "#6B8FA5"

    xs = [0.085, 0.285, 0.500, 0.700, 0.895]
    y = 0.55

    for x0, x1 in zip(xs[:-1], xs[1:]):
        ax.add_patch(
            patches.FancyArrowPatch(
                (x0 + 0.080, y), (x1 - 0.090, y),
                arrowstyle="-|>", mutation_scale=9.5, lw=0.9,
                color=accent, shrinkA=0, shrinkB=0,
                transform=ax.transAxes,
            )
        )

    # 1. FAP+ CAF niche: activated spindle cell.
    fibro = patches.Polygon(
        [
            (xs[0] - 0.070, y),
            (xs[0] - 0.022, y + 0.060),
            (xs[0] + 0.064, y + 0.030),
            (xs[0] + 0.075, y),
            (xs[0] + 0.020, y - 0.055),
            (xs[0] - 0.055, y - 0.030),
        ],
        closed=True, facecolor=accent_light, edgecolor=accent, lw=1.0,
        transform=ax.transAxes,
    )
    ax.add_patch(fibro)
    ax.add_patch(
        patches.Ellipse((xs[0] + 0.005, y), 0.040, 0.022, angle=12,
                        facecolor="#F8E7DF", edgecolor=accent, lw=0.65,
                        transform=ax.transAxes)
    )

    # 2. ECM/collagen sender fibrils.
    for k, off in enumerate([-0.034, -0.012, 0.012, 0.034]):
        t = np.linspace(-0.060, 0.060, 80)
        yy = y + off + 0.008 * np.sin((t + 0.06) * 42 + k)
        ax.plot(xs[1] + t, yy, color=accent if k < 2 else "#AD7A4C",
                lw=1.15, alpha=0.95, transform=ax.transAxes)
    ax.add_patch(
        patches.Ellipse((xs[1], y), 0.150, 0.100, facecolor="#FFF8F4",
                        edgecolor="#D5B5A7", lw=0.55, alpha=0.65,
                        transform=ax.transAxes)
    )

    # 3. Tumour epithelial receiver with ITGA2-centred integrins.
    ax.add_patch(
        patches.Ellipse((xs[2], y), 0.145, 0.122, facecolor="#F7F7F7",
                        edgecolor=neutral, lw=0.95, transform=ax.transAxes)
    )
    ax.add_patch(
        patches.Ellipse((xs[2] + 0.008, y - 0.006), 0.060, 0.042,
                        facecolor="#ECECF2", edgecolor="#8A8A9B", lw=0.65,
                        transform=ax.transAxes)
    )
    for dx in [-0.040, 0.000, 0.040]:
        base = (xs[2] + dx, y + 0.066)
        ax.plot([base[0], base[0]], [base[1] - 0.020, base[1] + 0.014],
                color=accent, lw=1.0, transform=ax.transAxes)
        ax.plot([base[0], base[0] - 0.016], [base[1] + 0.014, base[1] + 0.034],
                color=accent, lw=1.0, transform=ax.transAxes)
        ax.plot([base[0], base[0] + 0.016], [base[1] + 0.014, base[1] + 0.034],
                color=accent, lw=1.0, transform=ax.transAxes)

    # 4. FAK-ERK/MAPK signalling.
    for i, (lab, xx) in enumerate([("FAK", xs[3] - 0.036), ("ERK\nMAPK", xs[3] + 0.042)]):
        ax.add_patch(
            patches.Circle((xx, y), 0.043, facecolor="#FFF5F1",
                           edgecolor=accent, lw=0.95, transform=ax.transAxes)
        )
        ax.text(xx, y, lab, ha="center", va="center", fontsize=6.1,
                color="#4A2520", fontweight="bold", transform=ax.transAxes)
    ax.add_patch(
        patches.FancyArrowPatch((xs[3] - 0.002, y), (xs[3] + 0.020, y),
                                arrowstyle="-|>", mutation_scale=7.2,
                                lw=0.75, color=accent, transform=ax.transAxes)
    )
    ax.text(xs[3] + 0.087, y + 0.030, "↑", ha="center", va="center",
            fontsize=9.0, color=accent, fontweight="bold", transform=ax.transAxes)

    # 5. Dedifferentiation output: thyroid identity fades.
    ax.add_patch(
        patches.Ellipse((xs[4] - 0.018, y - 0.002), 0.122, 0.102,
                        facecolor="#F7F8F9", edgecolor=neutral, lw=0.9,
                        transform=ax.transAxes)
    )
    ax.add_patch(
        patches.Ellipse((xs[4] + 0.048, y + 0.018), 0.056, 0.030,
                        facecolor="none", edgecolor=thyroid, lw=0.85,
                        alpha=0.40, transform=ax.transAxes)
    )
    ax.add_patch(
        patches.Ellipse((xs[4] + 0.048, y - 0.020), 0.056, 0.030,
                        facecolor="none", edgecolor=thyroid, lw=0.85,
                        alpha=0.28, transform=ax.transAxes)
    )
    for lab, yy in zip(["PAX8", "NIS", "TG", "TPO"], [0.585, 0.558, 0.531, 0.504]):
        ax.text(xs[4] - 0.018, yy, lab, ha="center", va="center",
                fontsize=5.7, color="#8B98A3", alpha=0.58,
                transform=ax.transAxes)
    ax.text(xs[4] + 0.070, y - 0.046, "↓", ha="center", va="center",
            fontsize=9.0, color=accent, fontweight="bold", transform=ax.transAxes)

    labels = [
        ("FAP⁺ CAF\nniche", xs[0]),
        ("ECM/collagen\nsender program", xs[1]),
        ("ITGA2-centred\nintegrin receiver", xs[2]),
        ("FAK-ERK/MAPK\nsignalling ↑", xs[3]),
        ("Thyroid-lineage loss\n/ dedifferentiation", xs[4]),
    ]
    for lab, xx in labels:
        ax.text(xx, 0.265, lab, ha="center", va="top",
                fontsize=6.8, color="#222", fontweight="bold",
                linespacing=1.15, transform=ax.transAxes)

    sublabels = [
        "",
        "COL1A1 · COL1A2 · COL3A1\nCOL6A3 · POSTN · THBS2",
        "ITGA2 · ITGB1",
        "",
        "PAX8 · NIS · TG · TPO ↓",
    ]
    for lab, xx in zip(sublabels, xs):
        if lab:
            ax.text(xx, 0.108, lab, ha="center", va="top",
                    fontsize=5.8, color="#555", linespacing=1.12,
                    transform=ax.transAxes)

    ax.text(0.5, 0.925,
            "ITGA2-centred ECM-integrin axis linking FAP-high stroma to thyroid cancer dedifferentiation",
            ha="center", va="center", fontsize=8.8, fontweight="bold",
            color="#222", transform=ax.transAxes)


def draw_panel_B_support(fig, gs):
    # Integrin-receptor-resolved sender->receiver support, fully data-driven from the
    # ranked candidate-interaction table (no hardcoded values). Columns retain the
    # ITGA2-centred ECM-integrin axis and CXCL12-CXCR4 as a lower-weight context branch.
    # The former HGF-MET / "DCN-MET" column was removed 2026-06-02: HGF was deleted
    # project-wide (2026-05-12) and no DCN-MET interaction exists in the source data.
    ranked = pd.read_csv(full.FIBRO_P.parent / "ranked_candidate_interactions.tsv", sep="\t")
    col_spec = [   # (display label, pathway, receptor)
        ("ITGA2",          "ECM/integrin", "ITGA2"),
        ("ITGB1",          "ECM/integrin", "ITGB1"),
        ("CXCL12/\nCXCR4", "CXCL12-CXCR4", "CXCR4"),
    ]
    pathways = [c[0] for c in col_spec]
    sender_cell, sender_spat, recv_cell, recv_spat = [], [], [], []
    for _disp, _pth, _rec in col_spec:
        sub = ranked[(ranked.pathway == _pth) & (ranked.receptor == _rec)]
        sender_cell.append(sub["cell_delta_vs_noncore"].mean())
        sender_spat.append(sub["spatial_delta_seed_vs_other"].mean())
        recv_cell.append(sub["cell_delta_terminal_vs_lineage"].mean())
        recv_spat.append(sub["spatial_delta_terminal_vs_lineage"].mean())
    sender_cell = np.asarray(sender_cell); sender_spat = np.asarray(sender_spat)
    recv_cell = np.asarray(recv_cell); recv_spat = np.asarray(recv_spat)

    def row_max(x):  # proportional within-row scaling; keeps real positives off zero
        x = np.asarray(x, dtype=float)
        m = np.nanmax(x)
        return x / m if np.isfinite(m) and m > 0 else x

    fibro_support = (row_max(sender_cell) + row_max(sender_spat)) / 2
    recv_support = (row_max(recv_cell) + row_max(recv_spat)) / 2
    state_support = row_max(recv_spat)   # real terminal-vs-lineage spatial alignment
    mat = np.vstack([fibro_support, recv_support, state_support])
    integrated = mat.mean(axis=0)
    mat4 = np.vstack([mat, integrated])

    ax = fig.add_subplot(gs)
    ax.imshow(mat4, cmap=CMAP_SEQ, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(pathways)))
    ax.set_xticklabels(pathways, fontsize=6.2)
    ax.set_yticks(range(4))
    ax.set_yticklabels(["Fibro\nsender", "Epi\nreceiver", "State\nshift", "Integrated"], fontsize=6.2)
    for i in range(mat4.shape[0]):
        for j in range(mat4.shape[1]):
            v = mat4[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=6.3, fontweight="bold",
                    color="white" if v > 0.55 else "#3A2018")
    # thin white separators → tiled, editorial heatmap look
    for xx in np.arange(0.5, mat4.shape[1] - 0.5):
        ax.axvline(xx, color="white", lw=0.9)
    for yy in np.arange(0.5, mat4.shape[0] - 0.5):
        ax.axhline(yy, color="white", lw=0.9)
    ax.axhline(2.5, color="white", lw=1.8)   # bolder rule before the Integrated row
    for spine in ax.spines.values():
        spine.set_linewidth(0.45)
    ax.set_title("", fontsize=8.8, fontweight="bold", pad=3)


def draw_panel_C_liana(fig, gs):
    li = pd.read_csv(full.LIANA_F, sep="\t")
    li = li[~li["receptor.complex"].astype(str).str.contains("ITGA5", regex=False)]
    li = li.sort_values("aggregate_rank").head(8).copy()
    li["pair"] = (
        li["ligand.complex"].astype(str)
        + "→"
        + li["receptor.complex"].astype(str).str.replace("_", "/", regex=False)
    )
    # colour dots by RECEPTOR (adds information + harmonized warm sweep; ITGA2 = unified red)
    RECEPTOR_COLORS = {       # warm = integrins (ECM-red theme); cool = non-integrin receptors. Okabe-Ito, colour-blind safe
        "ITGA2": "#C0392B",   # red    — integrin (matches figure-wide ECM red)
        "ITGB1": "#E8843A",   # orange — integrin
        "ITGA3": "#E69F00",   # amber  — integrin
        "CD44":  "#0072B2",   # blue   — non-integrin
        "SDC1":  "#009E73",   # green  — non-integrin
    }
    def _recep_key(rc):
        for t in str(rc).split("_"):
            if t in RECEPTOR_COLORS:
                return t
        return str(rc).split("_")[0]
    li["recep"] = li["receptor.complex"].map(_recep_key)
    li = li.iloc[::-1].reset_index(drop=True)

    method_cols = [
        ("natmi.rank", "NATMI"),
        ("connectome.rank", "Conn."),
        ("logfc.rank", "logFC"),
        ("sca.rank", "SCsR"),
        ("cellphonedb.rank", "CPDB"),
    ]
    sub = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=gs, width_ratios=[1.10, 1.85], wspace=0.035
    )
    ax_lab = fig.add_subplot(sub[0])
    ax = fig.add_subplot(sub[1], sharey=ax_lab)

    y = np.arange(len(li))

    ax_lab.set_xlim(0, 1)
    ax_lab.set_ylim(-0.6, len(li) - 0.4)
    ax_lab.axis("off")
    for yi, label in zip(y, li["pair"].tolist()):
        ax_lab.text(
            0.98, yi, label,
            ha="right", va="center", fontsize=5.8, color="#222",
            transform=ax_lab.transData,
        )

    for xi, (col, _) in enumerate(method_cols):
        vals = li[col].astype(float).values
        nlog = -np.log10(np.maximum(vals, 1))
        mn, mx = np.nanmin(nlog), np.nanmax(nlog)
        sizes = 18 + (nlog - mn) / max(mx - mn, 1e-9) * 145
        colors = li["recep"].map(lambda r: RECEPTOR_COLORS.get(r, "#9E9E9E"))
        for yi in y:
            ax.scatter(xi, yi, s=sizes[yi], color=colors.iloc[yi],
                       edgecolor="black", linewidth=0.28, alpha=0.92, zorder=3)
    ax.set_xticks(range(len(method_cols)))
    ax.set_xticklabels([x[1] for x in method_cols], rotation=32, ha="right", fontsize=5.9)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0, pad=0)
    ax.set_xlim(-0.5, len(method_cols) - 0.5)
    ax.set_ylim(-0.6, len(li) - 0.4)
    ax.grid(axis="both", color="#EDEDED", linewidth=0.35, zorder=1)
    full.hide_spines(ax, keep=("left", "bottom"))
    ax.set_title("", fontsize=8.8, fontweight="bold", pad=3)


def draw_panel_D_programs(fig, gs):
    # Data-driven from the Hallmark sender/receiver program summary (no hardcoded scores; 2026-06-02).
    # Replaces a stale hardcoded IL6/JAK = 1.26 (source table value is 0.86).
    hall = pd.read_csv(full.HALL_P, sep="\t").set_index("hallmark")
    sel = [   # curated "shared downstream" set: display label -> hallmark key
        ("EMT",       "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
        ("IFN-γ",     "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
        ("Inflam.",   "HALLMARK_INFLAMMATORY_RESPONSE"),
        ("IL6/JAK",   "HALLMARK_IL6_JAK_STAT3_SIGNALING"),
        ("TNF/NF-κB", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
    ]
    programs = sorted([(lab, float(hall.loc[key, "mean_signed_score"])) for lab, key in sel],
                      key=lambda t: t[1], reverse=True)   # rank by score, descending
    ax = fig.add_subplot(gs)
    names = [x[0] for x in programs]
    vals = [x[1] for x in programs]
    y = np.arange(len(programs))[::-1]
    vmax = max(vals)
    bar_colors = [CMAP_SEQ(0.45 + 0.5 * (v / vmax)) for v in vals]   # value-graded (magnitude-encoded), refined ramp
    ax.barh(y, vals, color=bar_colors, edgecolor="white", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=6.3)
    ax.set_xlabel("Mean signed score", fontsize=6.7)
    for yy, v in zip(y, vals):
        ax.text(v + 0.05, yy, f"{v:.2f}", va="center", ha="left", fontsize=5.8)
    full.hide_spines(ax, keep=("left", "bottom"))
    ax.set_title("", fontsize=8.8, fontweight="bold", pad=3)


def draw_panel_E_spatial(fig, gs):
    d = pd.read_csv(AUG_CTX, sep="\t")
    d = d[d["sample"].isin(REP_SAMPLES)].copy()
    d["sample"] = pd.Categorical(d["sample"], categories=REP_SAMPLES, ordered=True)
    d["integrin_receiver_program_activity"] = d["itga2_receiver_program_activity"]
    d["local_ecm_integrin_support_index"] = (
        scale01(d["fap_niche_proximity_index"].values)
        + scale01(d["collagen_sender_program_activity"].values)
        + scale01(d["itga2_receiver_program_activity"].values)
    ) / 3.0

    rows = [
        ("fap_niche_proximity_index", CMAP_SEQ, "FAP-high niche\nproximity", "Niche index", "all"),
        ("local_ecm_integrin_support_index", CMAP_SUPPORT, "ECM-integrin\nsupport index", "Support", "all"),
        ("integrin_receiver_program_activity", CMAP_INTEG, "ITGA2\nreceiver", "Program", "epi"),
        ("epi_deconv_load_lineage_preserved_epithelial", CMAP_LINEAGE, "Lineage-\npreserved", "State load", "epi"),
        ("epi_deconv_load_terminal_dedifferentiated", CMAP_TERMINAL, "Terminal\ndedifferentiated", "State load", "epi"),
    ]

    g = gridspec.GridSpecFromSubplotSpec(
        len(rows), 7, subplot_spec=gs,
        width_ratios=[0.50, 1, 1, 1, 1, 1, 0.22],
        wspace=0.045, hspace=0.10,
    )
    for ri, (col, cmap, label, cbar_label, gate) in enumerate(rows):
        ax_label = fig.add_subplot(g[ri, 0])
        ax_label.axis("off")
        ax_label.text(0.98, 0.5, label, ha="right", va="center",
                      fontsize=7.2, fontweight="bold", color="#222",
                      transform=ax_label.transAxes)
        cap = float(np.nanquantile(d[col].values, 0.99))
        sc = None
        for ci, sample in enumerate(REP_SAMPLES):
            ax = fig.add_subplot(g[ri, ci + 1])
            ds = d[d["sample"] == sample]
            if gate == "fibro":
                fg = ds[ds["candidate_fibro_spot"] == True]
            elif gate == "epi":
                fg = ds[ds["epi_valid_primary"] == True]
            else:
                fg = ds
            ax.scatter(ds["x_plot"], ds["y_plot"], s=2.4, c="#E7EBF0",
                       edgecolors="none", alpha=0.86, rasterized=True)
            sc = ax.scatter(fg["x_plot"], fg["y_plot"], s=2.7,
                            c=fg[col].values, cmap=cmap, vmin=0, vmax=cap,
                            edgecolors="none", alpha=0.97, rasterized=True)
            ax.set_xlim(-0.55, 0.55)
            ax.set_ylim(0.55, -0.55)
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if ri == 0:
                ax.set_title(sample, fontsize=7.5, fontweight="bold", pad=2)
            if ri == len(rows) - 1 and ci == 0:
                add_scalebar(ax)

        ax_cb_host = fig.add_subplot(g[ri, 6])
        ax_cb_host.axis("off")
        cax = inset_axes(ax_cb_host, width="45%", height="80%", loc="center left",
                         bbox_to_anchor=(0.02, 0, 1, 1),
                         bbox_transform=ax_cb_host.transAxes, borderpad=0)
        cb = plt.colorbar(sc, cax=cax)
        cb.set_label(cbar_label, fontsize=6.2)
        cb.ax.tick_params(labelsize=5.5, width=0.40, length=1.8)
        cb.outline.set_linewidth(0.50)


def draw_panel_F_distance(fig, gs):
    r = pd.read_csv(full.RING, sep="\t")
    r = r[
        (r.threshold_label == "top20")
        & (r.gate == "primary")
        & (r.cohort == "formal_noP1")
    ]
    states_keep = ["Lineage-preserved epithelial", "Terminal dedifferentiated"]
    r = r[r.state_display.isin(states_keep)].copy()
    ring_levels = ["0-1", "1-2", "2-3", "3-5", ">5"]
    r = r[r["ring_bin"].isin(ring_levels)].copy()
    r["ring_bin"] = pd.Categorical(r["ring_bin"], categories=ring_levels, ordered=True)
    r = r.sort_values(["state_display", "ring_bin"])

    ax = fig.add_subplot(gs)
    curves = {}
    label_map = {
        "Lineage-preserved epithelial": "Lineage-preserved",
        "Terminal dedifferentiated": "Terminal dediff.",
    }
    for state in states_keep:
        sub = r[r.state_display == state]
        x = np.arange(len(sub))
        ax.plot(
            x, sub.mean_load.values, "-o",
            color=full.STATE_COLORS[state], lw=2.25, markersize=4.8,
            markeredgecolor="white", markeredgewidth=0.5,
            label=label_map[state], zorder=3,
        )
        se = sub.se_load.values
        ax.fill_between(
            x, sub.mean_load.values - 1.96 * se,
            sub.mean_load.values + 1.96 * se,
            color=full.STATE_COLORS[state], alpha=0.22, linewidth=0,
            zorder=1,
        )
        curves[state] = sub.mean_load.values

    lin_y = curves.get("Lineage-preserved epithelial")
    ter_y = curves.get("Terminal dedifferentiated")
    if lin_y is not None and ter_y is not None and len(lin_y) == len(ter_y):
        idx_cross = int(np.argmin(np.abs(lin_y - ter_y)))
        y_cross = (lin_y[idx_cross] + ter_y[idx_cross]) / 2
        ax.annotate(
            "state\ncrossover",
            xy=(idx_cross, y_cross),
            xytext=(idx_cross + 0.55, y_cross + 0.045),
            fontsize=6.0, color="#222", ha="left",
            arrowprops=dict(arrowstyle="->", lw=0.55, color="#444",
                            connectionstyle="arc3,rad=-0.18"),
        )

    ax.set_xticks(range(len(ring_levels)))
    ax.set_xticklabels(ring_levels, fontsize=6.8)
    ax.set_xlabel("Ring distance from FAP-high fibro seed (spot units)", fontsize=7.2)
    ax.set_ylabel("Mean epithelial state load", fontsize=7.2)
    ax.legend(loc="upper left", frameon=False, fontsize=6.4, handlelength=1.7)
    full.hide_spines(ax, keep=("left", "bottom"))
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.set_title("", fontsize=8.8, fontweight="bold", pad=3)


def draw_panel_G_sample_shift(fig, gs):
    eff = pd.read_csv(full.SAMPLE_EFF, sep="\t")
    eff = eff[
        (eff.threshold_label == "top20")
        & (eff.gate == "primary")
        & (eff.cohort == "formal_noP1")
    ]
    states = ["Lineage-preserved epithelial", "Terminal dedifferentiated"]
    eff = eff[eff.state_display.isin(states)]

    ax = fig.add_subplot(gs)
    lin_vals, ter_vals = [], []
    for sample in full.SAMPLES:
        sub = eff[eff["sample"] == sample]
        lv = sub[sub.state_display == states[0]]["neighborhood_diff"].values
        tv = sub[sub.state_display == states[1]]["neighborhood_diff"].values
        if len(lv) == 0 or len(tv) == 0:
            continue
        lin_vals.append(float(lv[0]))
        ter_vals.append(float(tv[0]))
        ax.plot([0, 1], [lv[0], tv[0]], color="#8A8A8A", lw=0.65, alpha=0.60, zorder=1)
        ax.scatter(0, lv[0], s=42, color=full.STATE_COLORS[states[0]],
                   edgecolor="white", linewidth=0.6, zorder=2)
        ax.scatter(1, tv[0], s=42, color=full.STATE_COLORS[states[1]],
                   edgecolor="white", linewidth=0.6, zorder=2)

    if lin_vals and ter_vals:
        ax.plot([-0.18, 0.18], [np.mean(lin_vals), np.mean(lin_vals)],
                color=full.STATE_COLORS[states[0]], lw=2.6, solid_capstyle="round", zorder=3)
        ax.plot([0.82, 1.18], [np.mean(ter_vals), np.mean(ter_vals)],
                color=full.STATE_COLORS[states[1]], lw=2.6, solid_capstyle="round", zorder=3)
        _, p = stats.wilcoxon(lin_vals, ter_vals)
        p_fmt = f"{p:.1e}" if p < 1e-2 else f"{p:.2g}"
        ax.text(0.5, 0.96, f"P = {p_fmt}; n = {len(lin_vals)}",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=6.4, color="#222",
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#888", lw=0.35))
    ax.axhline(0, color="#888", lw=0.50, ls="--")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Lineage-\npreserved", "Terminal\ndediff."], fontsize=6.7)
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylabel("State-load Δ", fontsize=7.2, labelpad=1)
    ax.set_title("", fontsize=8.8, fontweight="bold", pad=3)
    full.hide_spines(ax, keep=("left", "bottom"))


def main():
    fig = plt.figure(figsize=(7.2, 10.9), dpi=300)
    fig.patch.set_facecolor("white")

    outer = gridspec.GridSpec(
        4, 1, figure=fig,
        height_ratios=[1.45, 1.42, 4.15, 1.55],
        hspace=0.27, left=0.060, right=0.985,
        top=0.985, bottom=0.050,
    )

    full.draw_panel_A(fig, outer[0])   # reverted 2026-06-05 per user: gene-row sender->receiver Panel A (schematic kept as draw_panel_A_schematic, unused)

    g_mid = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=outer[1],
        width_ratios=[1.18, 1.35, 0.92], wspace=0.45,
    )
    draw_panel_B_support(fig, g_mid[0])
    draw_panel_C_liana(fig, g_mid[1])
    draw_panel_D_programs(fig, g_mid[2])

    draw_panel_E_spatial(fig, outer[2])

    g_bottom = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=outer[3],
        width_ratios=[1.18, 0.88, 1.15], wspace=0.34,
    )
    draw_panel_F_distance(fig, g_bottom[0])
    draw_panel_G_sample_shift(fig, g_bottom[1])
    full.draw_panel_E2(fig, g_bottom[2])

    add_panel_label(fig, 0.014, 0.982, "A")
    add_panel_label(fig, 0.014, 0.815, "B")
    add_panel_label(fig, 0.365, 0.815, "C")
    add_panel_label(fig, 0.735, 0.815, "D")
    add_panel_label(fig, 0.014, 0.665, "E")
    add_panel_label(fig, 0.014, 0.218, "F")
    add_panel_label(fig, 0.365, 0.218, "G")
    add_panel_label(fig, 0.610, 0.218, "H")

    out_base = OUT / "Figure_7_moderate_NC_preview"
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(out_base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(out_base.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(out_base.with_suffix(".tif"), dpi=600, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    print(f"Saved: {out_base.with_suffix('.pdf')}")
    print(f"Saved: {out_base.with_suffix('.svg')}")
    print(f"Saved: {out_base.with_suffix('.png')}")
    print(f"Saved: {out_base.with_suffix('.tif')}")


if __name__ == "__main__":
    main()
