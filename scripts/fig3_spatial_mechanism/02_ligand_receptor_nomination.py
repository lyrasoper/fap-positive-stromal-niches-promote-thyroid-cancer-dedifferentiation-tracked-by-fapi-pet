# Cancer Research submission - figure code release
# Builds: Figure 3 (LIANA/NicheNet ligand-receptor nomination and spatial atlas).
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

"""
Fig. 7 Nature-style polish — Python-only preview/final render.

This script keeps the current scientific content and source tables intact, but
retunes the visual hierarchy for a cleaner high-impact journal page: tighter
typography, less saturated colors, more room for the spatial atlas, shorter
labels, and publication-friendly PNG/PDF exports.

Layout:
  Row A  (nomination heatmaps)           — fibro sender | epi receiver
  Row B  (maintenance vs rewiring)       — 4×4 support + top-5 LR + shared
  Row C  (spatial atlas, 4 rows × 8 s)   — niche | collagen | ITGA2 | shift
  Row D  (quantitative closure)          — distance | sample-level | per-spot
"""
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy import stats

# ---------- project paths ----------
ROOT = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>"))
FIG7_BASE = ROOT / "outputs" / "fig7_redesign"
# build_fig3_panelA_paired.py (Panel A helper) is a sibling module in this scripts dir
_SCRIPT_DIR = str(Path(__file__).resolve().parent)
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

# ---------- Nature-style global config ----------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 6.8,
    "axes.labelsize": 6.8,
    "axes.titlesize": 8.2,
    "xtick.labelsize": 6.0,
    "ytick.labelsize": 6.0,
    "legend.fontsize": 6.0,
    "axes.linewidth": 0.45,
    "xtick.major.width": 0.45,
    "ytick.major.width": 0.45,
    "xtick.major.size": 2.2,
    "ytick.major.size": 2.2,
    "lines.linewidth": 0.9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})

# palettes
CMAP_SEQ = LinearSegmentedColormap.from_list(
    "nature_seq", ["#FEF7F3", "#F2C8B8", "#C85B4A", "#7D1F1A"])
CMAP_DIV = LinearSegmentedColormap.from_list(
    "nature_div", ["#4F6C8F", "#A8B8CD", "#FAF7F5", "#E9B2A3", "#B9473C"])
# lineage = teal, terminal = red (consistent across Fig 7, 8, Supp)
STATE_COLORS = {
    "Follicular-lineage high":     "#2C7F5B",
    "Lineage-preserved epithelial": "#2C7F5B",
    "Plasticity-high transition":  "#E89B3C",
    "Terminal dedifferentiated":   "#C0392B",
}
BRANCH_COLORS = {"maintenance": "#B42318", "rewiring": "#E89B3C"}
# Sample colors (8 formal)
SAMPLES = ["P12", "P17", "P26", "P32", "P44", "P57", "P83", "P98"]

# ---------- data paths ----------
RESULTS = Path((os.environ.get("EXTERNAL_DATA", "<EXTERNAL_DATA>") + "/scRNA_atlas/results"))
OUT     = FIG7_BASE / "nature_preview_20260518"
OUT.mkdir(parents=True, exist_ok=True)
AUG_CTX = FIG7_BASE / "augmented_context_v3.tsv"
RING    = RESULTS / "20260410_spatial_state_load_wholeslide/state_load_ring_mean_summary.tsv"
SAMPLE_EFF = RESULTS / "20260410_spatial_state_load_wholeslide/sample_state_load_effects.tsv"
RANKED  = RESULTS / "20260410_fap_fibro_epi_mechanism_prioritization/ranked_candidate_interactions.tsv"
FIBRO_P = RESULTS / "20260410_fap_fibro_epi_mechanism_prioritization/fibro_pathway_program_summary.tsv"
RECV_P  = RESULTS / "20260410_fap_fibro_epi_mechanism_prioritization/epithelial_receiver_program_summary.tsv"
HALL_P  = RESULTS / "20260410_hallmark_sender_receiver_programs/hallmark_priority_summary.tsv"
LIANA_F = RESULTS / "20260421_fig7_liana_nichenet/tables/liana_focused_FAP_to_terminal_top50.tsv"
NN_F    = RESULTS / "20260421_fig7_liana_nichenet/tables/nichenet_ligand_activity.tsv"

ECM_LIGANDS = {"COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL6A1","COL6A2",
               "COL6A3","COL5A1","COL5A2","COL15A1","POSTN","FN1","THBS1",
               "THBS2","VCAN","LAMA4","LAMB1","LAMC1","TNC","NID1","HSPG2",
               "BGN","DCN","LUM","SPARC","MMP14"}
ECM_RECEPTOR_TOKENS = ("ITG", "CD44", "SDC", "DDR", "LRP")

def _pathway_family(lig, rec):
    lig, rec = str(lig).upper(), str(rec).upper()
    if lig in ECM_LIGANDS and any(rec.startswith(t) or t in rec
                                  for t in ECM_RECEPTOR_TOKENS):
        return "ECM/integrin"
    if lig == "DCN" and "MET" in rec:
        return "DCN-MET"
    if "HGF" in lig or "MET" in rec:
        return "HGF-MET"
    if "CXCL" in lig or "CXCR" in rec:
        return "CXCL12-CXCR4"
    if "TGFB" in lig or "TGFBR" in rec:
        return "TGF-beta"
    if "WNT" in lig or "FZD" in rec or "ROR" in rec:
        return "WNT"
    return "Other"

LIANA_PATHWAY_COLORS = {
    "ECM/integrin":  "#B42318",
    "DCN-MET":       "#D98519",
    "HGF-MET":       "#E89B3C",
    "CXCL12-CXCR4":  "#0B7A75",
    "TGF-beta":      "#6F4E9E",
    "WNT":           "#8A8F96",
    "Other":         "#4A5568",
}


# ---------- helpers ----------
def panel_label_fig(fig, x, y, letter, fontsize=12):
    fig.text(x, y, letter, fontsize=fontsize, fontweight='bold')


def hide_spines(ax, keep=()):
    for s in ("top", "right", "left", "bottom"):
        if s not in keep:
            ax.spines[s].set_visible(False)


def add_scalebar(ax, length_data=0.10, label="", loc=(0.04, 0.05), lw=1.8, color="black"):
    """Add scale bar in axis data coords (normalized coords since spots are scaled)."""
    x0 = loc[0]; y0 = loc[1]
    ax.plot([x0, x0 + length_data], [y0, y0], color=color, lw=lw,
            transform=ax.transAxes, solid_capstyle='butt')
    if label:
        ax.text(x0 + length_data/2, y0 + 0.02, label, fontsize=5.8,
                ha='center', va='bottom', transform=ax.transAxes)


def _read_tsv_guarded(p, **kwargs):
    """Read an upstream intermediate TSV, failing loudly if it is missing.

    These tables are produced by upstream pipeline steps (the
    20260410_*/20260421_* result directories under EXTERNAL_DATA, and the
    augmented_context_v3.tsv built under PROJECT_ROOT); see README.
    """
    p = Path(p)
    if not p.exists():
        raise FileNotFoundError(
            f"Missing {p} — produced by an upstream step; see README.")
    return pd.read_csv(p, **kwargs)


def sig_mark(q):
    if pd.isna(q): return ""
    if q < 0.001: return "***"
    if q < 0.01:  return "**"
    if q < 0.05:  return "*"
    return ""


# ============================================================
#  Panel A — Gene-level nomination (fibro sender + epi receiver)
# ============================================================
def draw_panel_A(fig, gs):
    """Paired sender → receiver L-R nomination panel.
    Uses per-gene Wilcoxon-tested delta values from
    results/20260421_fig7_liana_nichenet/tables/A_{sender,receiver}_*_deltas.tsv
    and the LIANA top-12 pairs for the central coupling curves.
    """
    from build_fig3_panelA_paired import (
        SENDER, RECV, LIANA,
        SENDER_MAINT, SENDER_REWIR, RECV_MAINT, RECV_REWIR,
        CURATED_EDGES, PATHWAY_COLORS, UP_COLOR, DOWN_COLOR,
        sig_stars,
        parse_liana_edges, draw_connection, draw_side_labels,
        sender_pathway, recv_pathway,
    )

    sender_df = _read_tsv_guarded(SENDER, sep="\t")
    recv_df   = _read_tsv_guarded(RECV,   sep="\t")
    liana_df  = _read_tsv_guarded(LIANA,  sep="\t")
    sender_genes = SENDER_MAINT + SENDER_REWIR
    recv_genes   = RECV_MAINT   + RECV_REWIR

    liana_edges = parse_liana_edges(liana_df, sender_genes, recv_genes, top_n=12)
    liana_set = {(l, r) for (l, r, _, _) in liana_edges}
    curated_extra = [(l, r, pw, 0.30) for (l, r, pw) in CURATED_EDGES
                     if (l, r) not in liana_set
                     and l in sender_genes and r in recv_genes]
    all_edges = ([(l, r, pw, s, "LIANA")   for (l, r, pw, s) in liana_edges]
               + [(l, r, pw, s, "curated") for (l, r, pw, s) in curated_extra])

    # 5-column inner layout
    gA = gridspec.GridSpecFromSubplotSpec(
        1, 5, subplot_spec=gs,
        width_ratios=[1.1, 1.15, 2.35, 1.15, 1.1], wspace=0.05)

    # ---- left sender Δ + %expr ----
    ax_sL = fig.add_subplot(gA[0])
    s_df = sender_df.set_index("gene").reindex(sender_genes).reset_index()
    y_s = np.arange(len(sender_genes))[::-1]
    d_s = s_df["delta_log_norm"].values
    p_s = s_df["p_adj_BH"].values
    pc_s = s_df["pct_grp1"].values
    xlim_s = max(abs(d_s)) * 1.4
    for i, (d, p, pc) in enumerate(zip(d_s, p_s, pc_s)):
        yi = y_s[i]
        col = UP_COLOR if d >= 0 else DOWN_COLOR
        ax_sL.barh(yi, -d, color=col, edgecolor="black",
                   linewidth=0.35, height=0.55, zorder=2)
        st = sig_stars(p)
        if st:
            ax_sL.text(-d - 0.08, yi, st, fontsize=5.8, fontweight="bold",
                       va='center', ha='right', color="#333", zorder=4)
    for i, pc in enumerate(pc_s):
        ax_sL.scatter(-xlim_s*1.03, y_s[i], s=28 + 190 * pc,
                      facecolor="#E2A736", edgecolor="#333",
                      linewidth=0.35, clip_on=False, zorder=3)
    ax_sL.axvline(0, color="black", linewidth=0.5)
    ax_sL.set_xlim(-xlim_s*1.15, 0.25)
    ax_sL.set_ylim(-0.5, len(sender_genes) - 0.5)
    ax_sL.set_yticks([])
    ax_sL.invert_xaxis()
    ax_sL.tick_params(axis='x', labelsize=5.7)
    hide_spines(ax_sL, keep=("bottom",))
    ax_sL.spines['bottom'].set_linewidth(0.4)
    ax_sL.set_xlabel("Δ log-norm (FAP-CAF − other)", fontsize=5.9)
    ax_sL.set_title("Sender", fontsize=8.4, fontweight='bold', loc='right',
                    color="#B42318", pad=3)

    # ---- sender gene label column ----
    ax_sLab = fig.add_subplot(gA[1])
    sender_groups = ["maintenance"]*len(SENDER_MAINT) + ["rewiring"]*len(SENDER_REWIR)
    draw_side_labels(ax_sLab, sender_genes, sender_groups,
                     {g: sender_pathway(g) for g in sender_genes})

    # ---- middle connection area ----
    ax_mid = fig.add_subplot(gA[2])
    ax_mid.set_xlim(0, 1); ax_mid.set_ylim(-0.5, max(len(sender_genes),
                                                      len(recv_genes)) - 0.5)
    ax_mid.axis('off')
    y_send_map = {g: (len(sender_genes) - 1) - i for i, g in enumerate(sender_genes)}
    y_recv_map = {g: (len(recv_genes)   - 1) - i for i, g in enumerate(recv_genes)}
    for (lig, rec, pw, strength, src) in sorted(all_edges,
                                                  key=lambda e: 0 if e[4]=="curated" else 1):
        if lig not in y_send_map or rec not in y_recv_map:
            continue
        c = PATHWAY_COLORS.get(pw, "#4A5568")
        ls = "solid" if src == "LIANA" else "dashed"
        draw_connection(ax_mid, y_send_map[lig], y_recv_map[rec], c, strength,
                        x_left=0.01, x_right=0.99, linestyle=ls)

    # ---- receiver gene label column ----
    ax_rLab = fig.add_subplot(gA[3])
    recv_groups = ["maintenance"]*len(RECV_MAINT) + ["rewiring"]*len(RECV_REWIR)
    draw_side_labels(ax_rLab, recv_genes, recv_groups,
                     {g: recv_pathway(g) for g in recv_genes})

    # ---- right receiver Δ + %expr ----
    ax_rR = fig.add_subplot(gA[4])
    r_df = recv_df.set_index("gene").reindex(recv_genes).reset_index()
    y_r = np.arange(len(recv_genes))[::-1]
    d_r = r_df["delta_log_norm"].values
    p_r = r_df["p_adj_BH"].values
    pc_r = r_df["pct_grp1"].values
    xlim_r = max(abs(d_r)) * 1.4
    for i, (d, p, pc) in enumerate(zip(d_r, p_r, pc_r)):
        yi = y_r[i]
        col = UP_COLOR if d >= 0 else DOWN_COLOR
        ax_rR.barh(yi, d, color=col, edgecolor="black",
                   linewidth=0.35, height=0.55, zorder=2)
        st = sig_stars(p)
        if st:
            tx = d + (0.05 if d >= 0 else -0.05)
            ha = "left" if d >= 0 else "right"
            ax_rR.text(tx, yi, st, fontsize=5.8, fontweight="bold",
                       va='center', ha=ha, color="#333", zorder=4)
    for i, pc in enumerate(pc_r):
        ax_rR.scatter(xlim_r*1.03, y_r[i], s=28 + 190 * pc,
                      facecolor="#E2A736", edgecolor="#333",
                      linewidth=0.35, clip_on=False, zorder=3)
    ax_rR.axvline(0, color="black", linewidth=0.5)
    ax_rR.set_xlim(-0.25, xlim_r*1.15)
    ax_rR.set_ylim(-0.5, len(recv_genes) - 0.5)
    ax_rR.set_yticks([])
    ax_rR.tick_params(axis='x', labelsize=5.7)
    hide_spines(ax_rR, keep=("bottom",))
    ax_rR.spines['bottom'].set_linewidth(0.4)
    ax_rR.set_xlabel("Δ log-norm (terminal − lineage)", fontsize=5.9)
    ax_rR.set_title("Receiver", fontsize=8.4, fontweight='bold', loc='left',
                    color="#C0392B", pad=3)

    # small legend for connection style
    ax_mid.text(0.5, -1.15,
                "solid = LIANA top-12   ·   dashed = curated-only   ·   "
                "line thickness ∝ −log10(aggregate rank)",
                ha='center', va='top', fontsize=5.5, style='italic',
                color='#666', transform=ax_mid.transData, clip_on=False)


# ============================================================
#  Panel B — Support heatmap + Top-5 LR + Shared downstream programs
# ============================================================
def draw_panel_B(fig, gs):
    fibro = _read_tsv_guarded(FIBRO_P, sep="\t")
    epi = _read_tsv_guarded(RECV_P, sep="\t")
    hall = pd.read_csv(HALL_P, sep="\t") if HALL_P.exists() else None

    # ----- support matrix -----
    # NOTE: Fourth pathway is renamed to "HGF-MET / DCN-MET" so the B1 label
    # reflects that LIANA rescues this axis via DCN→MET while canonical
    # HGF-MET remains weak. The underlying pathway key still maps to the
    # HGF-MET entry in fibro_pathway_program_summary.tsv and
    # epithelial_receiver_program_summary.tsv.
    pathways_data  = ["ECM/integrin", "CXCL12-CXCR4", "TGF-beta", "HGF-MET"]
    pathways       = ["ECM/integrin", "CXCL12-CXCR4", "TGF-beta", "HGF-MET /\nDCN-MET"]
    def scale01(x):
        x = np.asarray(x, dtype=float)
        lo, hi = np.nanmin(x), np.nanmax(x)
        if not np.isfinite(hi - lo) or hi - lo == 0:
            return np.full_like(x, 0.5)
        return (x - lo) / (hi - lo)

    fibro = fibro[fibro.pathway.isin(pathways_data)].set_index("pathway").reindex(pathways_data)
    epi   = epi[epi.pathway.isin(pathways_data)].set_index("pathway").reindex(pathways_data)

    sender_cell  = scale01(fibro["cell_delta"].values)
    sender_spat  = scale01(fibro["spatial_delta"].values)
    recv_cell    = scale01(epi["cell_delta"].values)
    recv_spat    = scale01(epi["spatial_delta"].values)
    fibro_support = (sender_cell + sender_spat) / 2
    recv_support  = (recv_cell + recv_spat) / 2
    integrated   = (fibro_support + recv_support) / 2

    # heatmap: 2 rows (sender / receiver); integrated → separate bar
    support_mat = pd.DataFrame({
        "Fibro sender support":         fibro_support,
        "Epithelial receiver support":  recv_support,
    }, index=pathways)
    integrated_series = pd.Series(integrated, index=pathways,
                                   name="Integrated")

    gB = gridspec.GridSpecFromSubplotSpec(
        1, 4, subplot_spec=gs,
        width_ratios=[1.05, 0.20, 1.46, 1.18], wspace=0.82)

    # B1 — support heatmap (2 rows, 4 pathways)
    # Pathway display labels: shorter, single-line where possible
    pathway_short = {"ECM/integrin":      "ECM /\nintegrin",
                     "CXCL12-CXCR4":      "CXCL12 /\nCXCR4",
                     "TGF-beta":          "TGF-β",
                     "HGF-MET /\nDCN-MET": "HGF-MET /\nDCN-MET"}
    pathway_xtick = [pathway_short.get(p, p) for p in pathways]
    # Y-axis labels: shorter, capitalised
    y_short = {"Fibro sender support":         "Fibro\nsender",
               "Epithelial receiver support":  "Epi\nreceiver"}
    y_labels = [y_short.get(c, c) for c in support_mat.columns]

    axB1 = fig.add_subplot(gB[0])
    im = axB1.imshow(support_mat.values.T, cmap=CMAP_SEQ, vmin=0, vmax=1,
                     aspect="auto")
    axB1.set_xticks(range(support_mat.shape[0]))
    axB1.set_xticklabels(pathway_xtick, fontsize=6.3, rotation=0, ha='center')
    axB1.set_yticks(range(support_mat.shape[1]))
    axB1.set_yticklabels(y_labels, fontsize=6.3)
    for i in range(support_mat.shape[1]):
        for j in range(support_mat.shape[0]):
            v = support_mat.values[j, i]
            axB1.text(j, i, f"{v:.2f}", ha='center', va='center',
                      fontsize=6.4, fontweight='bold',
                      color="white" if v > 0.62 else "#222")
    axB1.add_patch(Rectangle((-0.5, -0.5), 1, support_mat.shape[1],
                              fill=False, edgecolor="#B42318", lw=1.6))
    axB1.set_title("Mechanism support index",
                   fontsize=7.6, fontweight='bold', pad=3)
    cb = plt.colorbar(im, ax=axB1, fraction=0.05, pad=0.04)
    cb.ax.tick_params(labelsize=6.2)
    cb.set_label("Support (0-1)", fontsize=6.1)
    for s in axB1.spines.values(): s.set_linewidth(0.5)

    # B1-bar — integrated priority score as a compact vertical bar
    axB1bar = fig.add_subplot(gB[1])
    xs = np.arange(len(pathways))
    bar_colors = ["#B42318" if p == "ECM/integrin" else "#BCC3CC" for p in pathways]
    axB1bar.bar(xs, integrated_series.values, color=bar_colors,
                edgecolor='black', linewidth=0.4, width=0.72)
    for x, v in zip(xs, integrated_series.values):
        axB1bar.text(x, v + 0.025, f"{v:.2f}", ha='center', va='bottom',
                     fontsize=6.8, fontweight='bold', color='#222')
    axB1bar.set_xticks([])
    axB1bar.set_ylim(0, max(integrated_series.values) * 1.32)
    axB1bar.set_ylabel("Integrated\nsupport", fontsize=5.7)
    axB1bar.tick_params(axis='y', labelsize=5.6)
    axB1bar.set_title("Integrated", fontsize=6.8, fontweight='bold', pad=3)
    for s in ('top', 'right'): axB1bar.spines[s].set_visible(False)
    axB1bar.spines['left'].set_linewidth(0.5)
    axB1bar.spines['bottom'].set_linewidth(0.5)

    # B2 — LIANA consensus top-N L-R dotplot (real data)
    #   replaces the earlier custom "Integrated priority" 5-bar.
    #   Columns: 5 LIANA methods (dots · size = -log10 method rank) +
    #            aggregate rank bar (−log10) + NicheNet AUPR bar.
    axB2 = fig.add_subplot(gB[2])
    N_LIANA = 20   # expanded from 10 → 20 to surface COL1A1→ITGA5 (rank 17) and ANGPTL2→ITGA5_ITGB1, cross-figure link to Fig 8 panel g
    try:
        li = pd.read_csv(LIANA_F, sep="\t")
        nn = pd.read_csv(NN_F, sep="\t")
        nn_map = dict(zip(nn["test_ligand"], nn["aupr_corrected"]))
        li = li.sort_values("aggregate_rank").head(N_LIANA).copy()
        li["lr_pair"]  = li["ligand.complex"] + "→" + li["receptor.complex"]
        li["pathway"]  = [_pathway_family(l, r)
                          for l, r in zip(li["ligand.complex"],
                                          li["receptor.complex"])]
        li["nn_aupr"]  = li["ligand.complex"].map(nn_map)
        li = li.iloc[::-1].reset_index(drop=True)   # top rank at top of plot

        method_cols = [("natmi.rank",       "NATMI"),
                       ("connectome.rank",  "Conn."),
                       ("logfc.rank",       "logFC"),
                       ("sca.rank",         "SCsR"),
                       ("cellphonedb.rank", "CPDB")]
        # use gridspec inside axB2: turn axB2 into a container
        axB2.set_axis_off()
        subB = gridspec.GridSpecFromSubplotSpec(
            1, 3, subplot_spec=gB[2],
            width_ratios=[2.9, 0.85, 0.75], wspace=0.08)
        axM = fig.add_subplot(subB[0])
        axA = fig.add_subplot(subB[1])
        axN = fig.add_subplot(subB[2])

        ys = np.arange(len(li))
        for xi, (col, _) in enumerate(method_cols):
            vals = li[col].astype(float).values
            nlog = -np.log10(np.maximum(vals, 1))
            mn, mx = np.nanmin(nlog), np.nanmax(nlog)
            sizes = 28 + (nlog - mn) / max(mx - mn, 1e-9) * 200
            colors = li["pathway"].map(lambda p: LIANA_PATHWAY_COLORS.get(p, "#4A5568"))
            for yi in ys:
                axM.scatter(xi, yi, s=sizes[yi], color=colors.iloc[yi],
                            edgecolor="black", linewidth=0.35, alpha=0.92,
                            zorder=3)
        axM.set_xticks(range(len(method_cols)))
        axM.set_xticklabels([c[1] for c in method_cols],
                            fontsize=6.3, rotation=30, ha='right')
        axM.set_yticks(ys)
        # Bold + dark-red emphasis on ITGA5-containing pairs (cross-figure link to Fig 8 panel g)
        lr_labels = li["lr_pair"].tolist()
        axM.set_yticklabels(lr_labels, fontsize=5.9)
        for tick_label, lr in zip(axM.get_yticklabels(), lr_labels):
            if "ITGA5" in lr:
                tick_label.set_fontweight('bold')
                tick_label.set_color('#7A1A1A')   # darker than #B42318 for emphasis
        axM.set_xlim(-0.95, len(method_cols) - 0.5)   # left pad so x=0 bubbles clear the L-R labels
        axM.set_ylim(-0.6, len(li) - 0.4)
        axM.grid(axis='both', color='#eee', linewidth=0.4, zorder=1)
        hide_spines(axM, keep=("left", "bottom"))
        axM.spines['left'].set_linewidth(0.5)
        axM.spines['bottom'].set_linewidth(0.5)
        axM.set_title("LIANA L-R consensus ranking",
                      fontsize=7.6, fontweight='bold', pad=3, loc='left')

        # aggregate rank bar
        agg  = li["aggregate_rank"].astype(float).values
        nlog = -np.log10(np.maximum(agg, 1e-7))
        bar_colors = li["pathway"].map(lambda p: LIANA_PATHWAY_COLORS.get(p, "#4A5568"))
        axA.barh(ys, nlog, color=bar_colors, edgecolor='black',
                 linewidth=0.3, height=0.72)
        for yi, (v, r) in enumerate(zip(nlog, agg)):
            lbl = f"{r:.1e}" if r < 1e-3 else f"{r:.2g}"
            axA.text(v + max(nlog)*0.02, yi, lbl,
                     va='center', ha='left', fontsize=5.8, color='#222')
        axA.set_yticks([])
        axA.set_ylim(-0.6, len(li) - 0.4)
        axA.set_xlim(0, nlog.max() * 1.45)
        axA.set_xlabel("agg rank", fontsize=5.8)
        axA.tick_params(axis='x', labelsize=5.2)
        hide_spines(axA, keep=("left", "bottom"))
        axA.spines['left'].set_linewidth(0.5)
        axA.spines['bottom'].set_linewidth(0.5)
        axA.grid(axis='x', color='#eee', linewidth=0.4, zorder=0)

        # NicheNet AUPR
        nn_vals = li["nn_aupr"].astype(float).values
        nn_filled = np.where(np.isnan(nn_vals), 0, nn_vals)
        axN.barh(ys, nn_filled, color=bar_colors, edgecolor='black',
                 linewidth=0.3, height=0.72, alpha=0.9)
        for yi, v in enumerate(nn_vals):
            if not np.isnan(v):
                axN.text(v + 0.003, yi, f"{v:.2f}",
                         va='center', ha='left', fontsize=5.8, color='#222')
        axN.set_yticks([])
        axN.set_ylim(-0.6, len(li) - 0.4)
        xmax = max(0.001, np.nanmax(nn_filled))
        axN.set_xlim(0, xmax * 1.45)
        axN.set_xlabel("AUPR_c", fontsize=5.8)
        axN.tick_params(axis='x', labelsize=5.2)
        hide_spines(axN, keep=("left", "bottom"))
        axN.spines['left'].set_linewidth(0.5)
        axN.spines['bottom'].set_linewidth(0.5)
        axN.grid(axis='x', color='#eee', linewidth=0.4, zorder=0)
    except Exception as e:
        axB2.text(0.5, 0.5, f"LIANA data not available:\n{e}",
                  ha='center', va='center', fontsize=7, transform=axB2.transAxes)
        axB2.set_axis_off()

    # B3 — shared downstream programs
    # Curated top-5 tumor-relevant Hallmark programs by epithelium terminal-vs-
    # lineage FDR. Previous version was missing the top-3 FDR hit
    # "INTERFERON GAMMA RESPONSE" (FDR=1e-5, same tier as TNFA and INFLAM);
    # we now include it and keep IL6-JAK-STAT3 for the downstream-convergence
    # narrative because it's strongest in the spatial terminal-rich contrast.
    axB3 = fig.add_subplot(gB[3])
    programs_default = [
        ("EMT", 2.70),                                 # 5/5 contrast support
        ("Inflam. response",                  1.41),   # FDR=1e-5 ***
        ("IFN-γ response",                    2.36),   # FDR=1e-5
        ("TNFα/NF-κB",                        0.92),   # FDR=4e-6
        ("IL6/JAK-STAT3",                     1.26),   # spatial terminal-rich
    ]
    programs = programs_default
    if hall is not None and "hallmark" in hall.columns:
        try:
            # allow override from hall table if the curated set matches
            score_col = next((c for c in ["mean_signed_score", "signed_score",
                                          "priority_score", "score"]
                              if c in hall.columns), None)
            if score_col is not None:
                wanted = {p[0].upper() for p in programs_default}
                hh = hall.copy()
                hh["label_up"] = hh["hallmark"].str.replace("_", " ").str.replace(
                    "HALLMARK ", "", regex=False).str.strip().str.upper()
                sub = hh[hh["label_up"].isin(wanted)]
                if len(sub) == len(programs_default):
                    programs = list(zip(sub["label_up"].str.title(),
                                        sub[score_col].values))
        except Exception:
            programs = programs_default

    names = [p[0] for p in programs]
    vals  = [p[1] for p in programs]
    y = np.arange(len(programs))[::-1]
    bar_colors = ["#C94F3D" if v > 0 else "#4F6C8F" for v in vals]
    axB3.barh(y, vals, color=bar_colors, edgecolor="black", linewidth=0.4, alpha=0.92)
    axB3.set_yticks(y)
    axB3.set_yticklabels(names, fontsize=6.2)
    axB3.set_xlabel("Mean signed score", fontsize=6.4)
    axB3.axvline(0, color='#888', lw=0.4)
    for j, v in zip(y, vals):
        axB3.text(v + 0.05 if v >= 0 else v - 0.05, j, f"{v:.2f}",
                  va='center', ha='left' if v >= 0 else 'right', fontsize=5.8)
    hide_spines(axB3, keep=("left", "bottom"))
    axB3.spines['left'].set_linewidth(0.5); axB3.spines['bottom'].set_linewidth(0.5)
    axB3.set_title("Shared downstream programs",
                   fontsize=7.6, fontweight='bold', pad=3)


# ============================================================
#  Panel C — Spatial atlas (4 rows × 8 samples)
# ============================================================
def draw_panel_C(fig, gs):
    d = _read_tsv_guarded(AUG_CTX, sep="\t")
    d = d[d["sample"].isin(SAMPLES)].copy()

    # Compose the local ECM-integrin interaction support index
    # (extended 2026-05-26 to include ITGA5 alongside ITGA2 — ITGA5 is the
    # sole integrin in Fig 8 v13 ECM/focal-adhesion sentinel panel and the
    # strongest cross-cohort proteomic bridge among integrins)
    def _scale01(s):
        s = np.asarray(s, dtype=float)
        lo, hi = np.nanmin(s), np.nanmax(s)
        return (s - lo) / max(hi - lo, 1e-9)
    d["local_ecm_integrin_support_index"] = (
        _scale01(d["fap_niche_proximity_index"].values)
        + _scale01(d["collagen_sender_program_activity"].values)
        + _scale01(d["itga2_receiver_program_activity"].values)
        + _scale01(d["itga5_receiver_program_activity"].values)
    ) / 4.0

    # dedicated label column + 8 sample columns + reserved colorbar column
    # Tighter spacing + slimmer non-sample columns so individual spatial
    # maps are visibly larger.
    N_ROWS = 7   # added ITGA5 receiver row (was 6)
    nc_gridspec = gridspec.GridSpecFromSubplotSpec(
        N_ROWS, 11, subplot_spec=gs, wspace=0.028, hspace=0.085,
        width_ratios=[0.36] + [1]*8 + [0.025, 0.36])

    # 6 rows: niche + mechanism composite + full 4-state
    # State-specific palettes use one hue per state for instant
    # recognition across Figs 7/8 and Supp.
    CMAP_LIN   = LinearSegmentedColormap.from_list(
        "lin", ["#F3F7F6", "#9ED3C8", "#2C7F62", "#14433A"])    # teals
    CMAP_TERM  = LinearSegmentedColormap.from_list(
        "term", ["#FCE9E5", "#F29A7F", "#C93B26", "#3A0C09"])     # deep red (terminal endpoint)

    # 6-row layout: niche + two mechanism components + composite support +
    # two state-outcome rows (lineage DOWN, terminal UP).
    # Palette retuned (option C): red/tan/PURPLE/red/teal/DEEP-RED so Row C
    # shows a clear visual hierarchy (4 distinct hues + a 2-step red scale).
    CMAP_ECM    = LinearSegmentedColormap.from_list(
        "ecm", ["#FAF3F0", "#F1C27D", "#C07A39", "#6C3A16"])      # warm tan  (collagen)
    CMAP_INTEG  = LinearSegmentedColormap.from_list(
        "integ", ["#F6EEF3", "#C99BC3", "#7E3783", "#31114A"])    # purple    (ITGA2 receiver)
    CMAP_INTEG2 = LinearSegmentedColormap.from_list(
        "integ2", ["#FCEEE0", "#F2B36A", "#C9601F", "#5A2604"])   # deep orange (ITGA5 receiver — distinguishable from ITGA2 purple; signals cross-figure link to Fig 8)
    rows = [
        ("fap_niche_proximity_index",             CMAP_SEQ,
            "C1", "FAP-high niche\nproximity",
            "Niche index",       "seq", "all"),
        ("collagen_sender_program_activity",      CMAP_ECM,
            "C2", "Collagen sender\nprogram activity",
            "Program activity",  "seq", "fibro"),
        ("itga2_receiver_program_activity",       CMAP_INTEG,
            "C3", "ITGA2 receiver\nprogram activity",
            "Program activity",  "seq", "epi"),
        ("itga5_receiver_program_activity",       CMAP_INTEG2,
            "C4", "ITGA5 receiver\nprogram activity ★",  # ★ signals cross-figure link to Fig 8 panel g
            "Program activity",  "seq", "epi"),
        ("local_ecm_integrin_support_index",      CMAP_SEQ,
            "C5", "ECM–integrin\nsupport index",
            "Support (0-1)",     "seq", "all"),
        ("epi_deconv_load_lineage_preserved_epithelial",  CMAP_LIN,
            "C6", "Lineage-\npreserved",
            "State load",        "seq", "epi"),
        ("epi_deconv_load_terminal_dedifferentiated",     CMAP_TERM,
            "C7", "Terminal\ndedifferentiated",
            "State load",        "seq", "epi"),
    ]

    for ri, (col, cmap, tag, short_label, cbar_label, kind, gate) in enumerate(rows):
        # ---- dedicated label column (horizontal, right-aligned) ----
        ax_lab = fig.add_subplot(nc_gridspec[ri, 0])
        ax_lab.set_xticks([]); ax_lab.set_yticks([])
        for sp in ax_lab.spines.values(): sp.set_visible(False)
        ax_lab.set_xlim(0, 1); ax_lab.set_ylim(0, 1)
        ax_lab.text(0.98, 0.5, short_label,
                    transform=ax_lab.transAxes,
                    ha='right', va='center',
                    fontsize=6.9, fontweight='bold', color="#222")

        # compute robust caps
        if kind == "div":
            cap = float(np.nanquantile(np.abs(d[col].values), 0.99))
            vmin, vmax = -cap, cap
        else:
            cap = float(np.nanquantile(d[col].values, 0.99))
            vmin, vmax = 0.0, cap
        for ci, s in enumerate(SAMPLES):
            # sample axes are now in columns 1..8 (column 0 is the label)
            ax = fig.add_subplot(nc_gridspec[ri, ci + 1])
            ds = d[d["sample"] == s]
            # apply gate
            if gate == "fibro":
                ds_fg = ds[ds["candidate_fibro_spot"] == True]
            elif gate == "epi":
                ds_fg = ds[ds["epi_valid_primary"] == True]
            else:
                ds_fg = ds
            # background
            ax.scatter(ds["x_plot"], ds["y_plot"], s=1.28, c="#E8ECF1",
                       edgecolors='none', alpha=0.72)
            # foreground
            sc = ax.scatter(ds_fg["x_plot"], ds_fg["y_plot"], s=1.42,
                            c=ds_fg[col].values, cmap=cmap,
                            vmin=vmin, vmax=vmax, edgecolors='none', alpha=0.94)
            ax.set_xlim(-0.55, 0.55)
            ax.set_ylim(0.55, -0.55)  # invert Y as in R version
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if ri == 0:
                ax.set_title(s, fontsize=7.0, fontweight='bold', pad=1.6)
            if ri == N_ROWS - 1 and ci == 0:
                add_scalebar(ax, length_data=0.18, label="", lw=1.5)
        # row-level shared colorbar in the reserved rightmost column cell,
        # placed cleanly outside the sample grid. Column 9 is the visual gap.
        cax_host = fig.add_subplot(nc_gridspec[ri, 10])
        cax_host.set_xticks([]); cax_host.set_yticks([])
        for s in cax_host.spines.values(): s.set_visible(False)
        # place a thin colorbar at the left edge of the host cell so it reads
        # as a per-row annotation rather than a separate panel
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cax = inset_axes(cax_host, width="22%", height="80%",
                         loc='center left',
                         bbox_to_anchor=(0.0, 0, 1, 1),
                         bbox_transform=cax_host.transAxes,
                         borderpad=0)
        cb = plt.colorbar(sc, cax=cax)
        cb.set_label(cbar_label, fontsize=5.8)
        cb.ax.tick_params(labelsize=5.1)
        cb.outline.set_linewidth(0.4)


# ============================================================
#  Panel D — Distance-dependent redistribution  (2 curves + 95% CI)
# ============================================================
def draw_panel_D(fig, gs):
    r = _read_tsv_guarded(RING, sep="\t")
    # filter to top20 / primary / formal_noP1 / the two focus states
    r = r[(r.threshold_label == "top20") & (r.gate == "primary")
          & (r.cohort == "formal_noP1")]
    states_keep = ["Lineage-preserved epithelial", "Terminal dedifferentiated"]
    r = r[r.state_display.isin(states_keep)]
    # ring order (matches upstream categorical values)
    ring_levels = ["0-1", "1-2", "2-3", "3-5", ">5"]
    r = r[r["ring_bin"].isin(ring_levels)].copy()
    r["ring_bin"] = pd.Categorical(r["ring_bin"], categories=ring_levels, ordered=True)
    r = r.sort_values(["state_display", "ring_bin"])

    ax = fig.add_subplot(gs)
    # store curves to compute the crossover later
    curves = {}
    for state in states_keep:
        sub = r[r.state_display == state]
        x = np.arange(len(sub))
        ax.plot(x, sub.mean_load.values, '-o',
                color=STATE_COLORS[state], lw=2.0, markersize=4.4,
                markeredgecolor='white', markeredgewidth=0.5,
                label=state, zorder=3)
        se = sub.se_load.values
        ax.fill_between(x, sub.mean_load.values - 1.96*se,
                        sub.mean_load.values + 1.96*se,
                        color=STATE_COLORS[state], alpha=0.22, linewidth=0,
                        zorder=1)
        curves[state] = sub.mean_load.values

    # crossover annotation — find the ring bin closest to the two-curve intersection
    lin_y = curves.get("Lineage-preserved epithelial")
    ter_y = curves.get("Terminal dedifferentiated")
    if lin_y is not None and ter_y is not None and len(lin_y) == len(ter_y):
        diff = lin_y - ter_y
        idx_cross = int(np.argmin(np.abs(diff)))
        x_cross = idx_cross
        y_cross = (lin_y[idx_cross] + ter_y[idx_cross]) / 2
        ax.annotate("lineage ↔ terminal\ncrossover",
                    xy=(x_cross, y_cross),
                    xytext=(x_cross + 0.9, y_cross + 0.06),
                    fontsize=5.7, color="#222",
                    arrowprops=dict(arrowstyle="->", lw=0.55, color="#444",
                                    connectionstyle="arc3,rad=-0.2"))

    ax.set_xticks(range(len(ring_levels)))
    ax.set_xticklabels(ring_levels, fontsize=6.2)
    ax.set_xlabel("Ring distance from FAP-high fibro seed (spot units)",
                  fontsize=6.4)
    ax.set_ylabel("Mean epithelial state load", fontsize=6.4)
    ax.legend(loc='best', frameon=False, fontsize=5.8)
    hide_spines(ax, keep=("left", "bottom"))
    ax.spines['left'].set_linewidth(0.5); ax.spines['bottom'].set_linewidth(0.5)
    ax.set_title("Distance-dependent redistribution",
                 fontsize=8.0, fontweight='bold', pad=3)
    ax.axhline(0, color="#888", lw=0.4)


# ============================================================
#  Panel E1 — Sample-level paired dumbbell (lineage vs terminal)
# ============================================================
def draw_panel_E1(fig, gs):
    eff = _read_tsv_guarded(SAMPLE_EFF, sep="\t")
    eff = eff[(eff.threshold_label == "top20") & (eff.gate == "primary")
              & (eff.cohort == "formal_noP1")]
    eff = eff[eff.state_display.isin(
        ["Lineage-preserved epithelial", "Terminal dedifferentiated"])]

    ax = fig.add_subplot(gs)
    samples = sorted(eff["sample"].unique())
    lin_vals, ter_vals = [], []
    for s in samples:
        sub = eff[eff["sample"] == s]
        lv = sub[sub.state_display == "Lineage-preserved epithelial"]["neighborhood_diff"].values
        tv = sub[sub.state_display == "Terminal dedifferentiated"]["neighborhood_diff"].values
        if not len(lv) or not len(tv):
            continue
        lin_vals.append(lv[0]); ter_vals.append(tv[0])
        # connecting line
        ax.plot([0, 1], [lv[0], tv[0]],
                color="#888", lw=0.6, alpha=0.55, zorder=1)
        # points
        ax.scatter(0, lv[0], s=46,
                   color=STATE_COLORS["Lineage-preserved epithelial"],
                   edgecolor="white", linewidth=0.6, zorder=2)
        ax.scatter(1, tv[0], s=46,
                   color=STATE_COLORS["Terminal dedifferentiated"],
                   edgecolor="white", linewidth=0.6, zorder=2)
        # sample labels on OUTSIDE (left of lineage, right of terminal)
        ax.text(-0.07, lv[0], s, fontsize=6.3, va='center', ha='right',
                color='#444')
        ax.text(1.07, tv[0], s, fontsize=6.3, va='center', ha='left',
                color='#444')
    # mean lines for each state
    if lin_vals and ter_vals:
        lin_mean = float(np.mean(lin_vals))
        ter_mean = float(np.mean(ter_vals))
        ax.plot([-0.18, 0.18], [lin_mean, lin_mean],
                color=STATE_COLORS["Lineage-preserved epithelial"],
                lw=2.6, zorder=3, solid_capstyle='round')
        ax.plot([0.82, 1.18], [ter_mean, ter_mean],
                color=STATE_COLORS["Terminal dedifferentiated"],
                lw=2.6, zorder=3, solid_capstyle='round')
    ax.axhline(0, color='#888', lw=0.5, ls='--')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Lineage-\npreserved", "Terminal\ndedifferentiated"],
                       fontsize=6.1)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel("Per-sample Δ", fontsize=6.4)
    ax.set_title("Per-sample paired shift",
                 fontsize=8.0, fontweight='bold', pad=3)
    hide_spines(ax, keep=("left", "bottom"))
    ax.spines['left'].set_linewidth(0.5); ax.spines['bottom'].set_linewidth(0.5)
    # paired test annotation — aggressive simplification
    if len(lin_vals) > 2:
        try:
            _, p = stats.wilcoxon(lin_vals, ter_vals)
            p_fmt = f"{p:.1e}" if p < 1e-2 else f"{p:.2g}"
            ax.text(0.5, 0.97,
                    f"P = {p_fmt} ({len(lin_vals)}/{len(lin_vals)})",
                    transform=ax.transAxes, ha='center', va='top',
                    fontsize=5.8, color='#222',
                    bbox=dict(boxstyle="round,pad=0.26", fc="white",
                              ec="#888", lw=0.35))
        except Exception:
            pass


# ============================================================
#  Panel E2 — Per-spot scatter (hexbin)
# ============================================================
def draw_panel_E2(fig, gs):
    """Per-spot ECM-integrin interaction support × dedifferentiation shift.
    Uses the composite ECM-integrin support index (min-max mean of niche
    proximity + collagen sender + ITGA2 + ITGA5 receivers). Updated
    2026-05-26 to include ITGA5 alongside ITGA2 (cross-figure link to Fig 8
    v13 ECM/focal-adhesion sentinel panel)."""
    d = _read_tsv_guarded(AUG_CTX, sep="\t")
    # Compose the composite support index (same formula as Panel C)
    def _scale01(s):
        s = np.asarray(s, dtype=float)
        lo, hi = np.nanmin(s), np.nanmax(s)
        return (s - lo) / max(hi - lo, 1e-9)
    d = d.dropna(subset=["fap_niche_proximity_index",
                         "collagen_sender_program_activity",
                         "itga2_receiver_program_activity",
                         "itga5_receiver_program_activity",
                         "dediff_shift_index"])
    # Restrict to epi-valid spots (matches the F3 R-pipeline scope)
    if "epi_valid_primary" in d.columns:
        d = d[d["epi_valid_primary"] == True]
    d["local_ecm_integrin_support_index"] = (
        _scale01(d["fap_niche_proximity_index"].values)
        + _scale01(d["collagen_sender_program_activity"].values)
        + _scale01(d["itga2_receiver_program_activity"].values)
        + _scale01(d["itga5_receiver_program_activity"].values)
    ) / 4.0
    # Per-sample ρ (computed on the full epi-valid set BEFORE down-sampling)
    per_sample_rhos = []
    if "sample" in d.columns:
        for s, dd in d.groupby("sample"):
            if len(dd) >= 30:
                try:
                    r_s, _ = stats.spearmanr(
                        dd["local_ecm_integrin_support_index"].values,
                        dd["dediff_shift_index"].values,
                        nan_policy='omit')
                    if np.isfinite(r_s):
                        per_sample_rhos.append(r_s)
                except Exception:
                    pass

    x = d["local_ecm_integrin_support_index"].values
    y = d["dediff_shift_index"].values
    n_sample = 25000
    if len(x) > n_sample:
        ix = np.random.default_rng(42).choice(len(x), n_sample, replace=False)
        x = x[ix]; y = y[ix]
    rho, p = stats.spearmanr(x, y, nan_policy='omit')

    ax = fig.add_subplot(gs)
    # Neutral blue-grey density colormap so the hexbin doesn't visually
    # compete with Row C's red-terminal / C6 palette.
    density_cmap = LinearSegmentedColormap.from_list(
        "density_neutral",
        ["#F4F6FA", "#B0C4D9", "#4F6C8F", "#1A2A44"])
    hb = ax.hexbin(x, y, gridsize=42, cmap=density_cmap, bins='log',
                   mincnt=1, linewidths=0)
    m, b = np.polyfit(x, y, 1)
    xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    ax.plot(xs, m*xs + b, color="#1A1A1A", lw=1.6)

    # Aggressive simplification: title-only ρ, statistical detail lives in legend.
    stat_text = f"ρ = {rho:+.2f}"
    ax.text(0.03, 0.97, stat_text,
            transform=ax.transAxes, fontsize=6.7, va='top', fontweight='bold',
            bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="#888", lw=0.4))
    ax.axhline(0, color='#888', lw=0.4)
    ax.set_xlabel("ECM-integrin support", fontsize=6.4)
    ax.set_ylabel("Dediff shift", fontsize=6.4)
    hide_spines(ax, keep=("left", "bottom"))
    ax.spines['left'].set_linewidth(0.5); ax.spines['bottom'].set_linewidth(0.5)
    ax.set_title("Spot-level closure",
                 fontsize=8.0, fontweight='bold', pad=3)
    cb = plt.colorbar(hb, ax=ax, fraction=0.048, pad=0.02)
    cb.set_label("log10 spot count", fontsize=5.8)
    cb.ax.tick_params(labelsize=5.1)


# ============================================================
#  Compose
# ============================================================
def main():
    fig = plt.figure(figsize=(8.27, 9.2), dpi=170)   # 210mm × 234mm — legible compromise; production rescales
    fig.patch.set_facecolor("white")

    outer = gridspec.GridSpec(
        4, 1, figure=fig,
        height_ratios=[1.18, 1.65, 4.95, 1.30],
        hspace=0.45, left=0.057, right=0.972,
        top=0.986, bottom=0.043)

    draw_panel_A(fig, outer[0])
    draw_panel_B(fig, outer[1])
    draw_panel_C(fig, outer[2])

    gD = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=outer[3],
        width_ratios=[1.12, 0.90, 1.28], wspace=0.32)
    draw_panel_D(fig, gD[0])
    draw_panel_E1(fig, gD[1])
    draw_panel_E2(fig, gD[2])

    # Top-level NC-style panel letters A–H aligned to each visual unit.
    # Row B gridspec  (wr=[1.00, 0.28, 1.80, 0.85], wspace=0.42):
    #   B  → col 0+1 (support heatmap + Integrated bar)     x ≈ 0.020
    #   C  → col 2   (LIANA consensus)                       x ≈ 0.420
    #   D  → col 3   (shared downstream)                     x ≈ 0.800
    # Row D gridspec  (wr=[1.1, 1.0, 1.25], wspace=0.35):
    #   F  → col 0   (distance curves)                       x ≈ 0.020
    #   G  → col 1   (sample-level shift)                    x ≈ 0.330
    #   H  → col 2   (per-spot closure)                      x ≈ 0.635
    panel_label_fig(fig, 0.019, 0.978, "A", fontsize=11.2)
    panel_label_fig(fig, 0.019, 0.808, "B", fontsize=11.2)
    panel_label_fig(fig, 0.398, 0.808, "C", fontsize=11.2)
    panel_label_fig(fig, 0.765, 0.808, "D", fontsize=11.2)
    panel_label_fig(fig, 0.019, 0.600, "E", fontsize=11.2)
    panel_label_fig(fig, 0.019, 0.137, "F", fontsize=11.2)
    panel_label_fig(fig, 0.341, 0.137, "G", fontsize=11.2)
    panel_label_fig(fig, 0.617, 0.137, "H", fontsize=11.2)

    # NB: no suptitle on the figure — the figure caption is provided
    # separately in the manuscript legend (NC submission convention).

    out_png = OUT / "Figure_7_Nature_style_preview.png"
    out_pdf = OUT / "Figure_7_Nature_style_preview.pdf"
    out_svg = OUT / "Figure_7_Nature_style_preview.svg"
    fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
    fig.savefig(out_svg, bbox_inches='tight', facecolor='white')
    print("Saved:", out_png)
    print("Saved:", out_pdf)
    print("Saved:", out_svg)


if __name__ == "__main__":
    main()
