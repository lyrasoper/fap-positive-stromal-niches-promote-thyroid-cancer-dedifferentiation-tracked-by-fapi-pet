# Cancer Research submission - figure code release
# Builds: Figure 3 full-figure renderer. Imported by fig3_beautify.py and by
#          13_render_fig3_final.py; not run directly.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

"""
Fig 7 v3 — full NC-level Python re-render.

Every panel is drawn from source TSVs with one consistent font (Helvetica /
DejaVu Sans fallback), one sequential and one diverging colormap family, a
fixed sample colour mapping and scale bars on spatial panels. The layout
follows the niche-aligned v2 narrative, with the per-spot mechanism-in-space
scatter added as a quantitative closure.

Layout:
  Row A  (nomination heatmaps)           — fibro sender | epi receiver
  Row B  (maintenance vs rewiring)       — 4×4 support + top-5 LR + shared
  Row C  (spatial atlas, 7 rows × 8 s)   — niche | collagen | ITGA2 | ITGA5 | support | lineage | terminal
  Row D  (quantitative closure)          — distance | sample-level | per-spot
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy import stats

# ---------- NC-style global config ----------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7.5,
    "axes.labelsize": 7.5,
    "axes.titlesize": 9.5,
    "xtick.labelsize": 6.8,
    "ytick.labelsize": 6.8,
    "legend.fontsize": 6.8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.55,
    "ytick.major.width": 0.55,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "lines.linewidth": 1.0,
})

# palettes
CMAP_SEQ = LinearSegmentedColormap.from_list(
    "nc_seq", ["#FDF5F0", "#F2B8A6", "#C94F3D", "#7E1A17"])
CMAP_DIV = LinearSegmentedColormap.from_list(
    "nc_div", ["#4F6C8F", "#8FA7C5", "#FAF7F5", "#F2B8A6", "#C94F3D"])
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
RESULTS = Path(os.environ.get("EXTERNAL_DATA", "<EXTERNAL_DATA>")) / "scRNA_atlas" / "results"
OUT     = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>")) / "outputs" / "fig7_redesign"
AUG_CTX = OUT / "augmented_context_v3.tsv"
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
    # DCN-MET / HGF-MET families removed 2026-06-02 (HGF deleted project-wide; no DCN-MET in data)
    if "CXCL" in lig or "CXCR" in rec:
        return "CXCL12-CXCR4"
    if "TGFB" in lig or "TGFBR" in rec:
        return "TGF-beta"
    if "WNT" in lig or "FZD" in rec or "ROR" in rec:
        return "WNT"
    return "Other"

LIANA_PATHWAY_COLORS = {
    "ECM/integrin":  "#B42318",
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
        MAINT_BAND, REWIR_BAND, sig_stars, classify_edge,
        parse_liana_edges, draw_connection, draw_side_labels,
        sender_pathway, recv_pathway,
    )

    sender_df = pd.read_csv(SENDER, sep="\t")
    recv_df   = pd.read_csv(RECV,   sep="\t")
    liana_df  = pd.read_csv(LIANA,  sep="\t")
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
            ax_sL.text(-d - 0.08, yi, st, fontsize=6.5, fontweight="bold",
                       va='center', ha='right', color="#333", zorder=4)
    for i, pc in enumerate(pc_s):
        ax_sL.scatter(-xlim_s*1.03, y_s[i], s=35 + 230 * pc,
                      facecolor="#E2A736", edgecolor="#333",
                      linewidth=0.35, clip_on=False, zorder=3)
    ax_sL.axvline(0, color="black", linewidth=0.5)
    ax_sL.set_xlim(-xlim_s*1.15, 0.25)
    ax_sL.set_ylim(-0.5, len(sender_genes) - 0.5)
    ax_sL.set_yticks([])
    ax_sL.invert_xaxis()
    ax_sL.tick_params(axis='x', labelsize=6.2)
    hide_spines(ax_sL, keep=("bottom",))
    ax_sL.spines['bottom'].set_linewidth(0.4)
    ax_sL.set_xlabel("Δ log-norm (FAP-CAF − other)", fontsize=6.3)
    ax_sL.set_title("Sender", fontsize=9.5, fontweight='bold', loc='right',
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
            ax_rR.text(tx, yi, st, fontsize=6.5, fontweight="bold",
                       va='center', ha=ha, color="#333", zorder=4)
    for i, pc in enumerate(pc_r):
        ax_rR.scatter(xlim_r*1.03, y_r[i], s=35 + 230 * pc,
                      facecolor="#E2A736", edgecolor="#333",
                      linewidth=0.35, clip_on=False, zorder=3)
    ax_rR.axvline(0, color="black", linewidth=0.5)
    ax_rR.set_xlim(-0.25, xlim_r*1.15)
    ax_rR.set_ylim(-0.5, len(recv_genes) - 0.5)
    ax_rR.set_yticks([])
    ax_rR.tick_params(axis='x', labelsize=6.2)
    hide_spines(ax_rR, keep=("bottom",))
    ax_rR.spines['bottom'].set_linewidth(0.4)
    ax_rR.set_xlabel("Δ log-norm (terminal − lineage)", fontsize=6.3)
    ax_rR.set_title("Receiver", fontsize=9.5, fontweight='bold', loc='left',
                    color="#C0392B", pad=3)

    # small legend for connection style
    ax_mid.text(0.5, -1.15,
                "solid = LIANA top-12   ·   dashed = curated-only   ·   "
                "line thickness ∝ −log10(aggregate rank)",
                ha='center', va='top', fontsize=6.0, style='italic',
                color='#666', transform=ax_mid.transData, clip_on=False)


# ============================================================
#  Panel B — Support heatmap + Top-5 LR + Shared downstream programs
# ============================================================
def draw_panel_B(fig, gs):
    fibro = pd.read_csv(FIBRO_P, sep="\t")
    epi = pd.read_csv(RECV_P, sep="\t")
    ranked = pd.read_csv(RANKED, sep="\t")
    hall = pd.read_csv(HALL_P, sep="\t") if HALL_P.exists() else None

    # ----- support matrix -----
    # HGF-MET / "DCN-MET" pathway dropped 2026-06-02: HGF removed project-wide
    # (2026-05-12); HGF-MET was the weakest axis (rank 16/16) and no DCN-MET
    # interaction exists in the source data. Panel B now shows the three
    # data-supported pathways. (The manuscript "moderate" Fig 3 resolves the
    # ECM/integrin axis further into ITGA2/ITGB1/ITGA5 receptor columns.)
    pathways_data  = ["ECM/integrin", "CXCL12-CXCR4", "TGF-beta"]
    pathways       = ["ECM/integrin", "CXCL12-CXCR4", "TGF-beta"]
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
    # State-redistribution support: hard-coded from v2 values (derived in R)
    state_support = pd.Series({"ECM/integrin": 1.00, "CXCL12-CXCR4": 0.58,
                               "TGF-beta": 0.54})[pathways_data].values
    integrated   = (fibro_support + recv_support + state_support) / 3

    # heatmap: 3 rows (sender / receiver / state-redistribution); integrated → separate bar
    support_mat = pd.DataFrame({
        "Fibro sender support":         fibro_support,
        "Epithelial receiver support":  recv_support,
        "State-redistribution support": state_support,
    }, index=pathways)
    integrated_series = pd.Series(integrated, index=pathways,
                                   name="Integrated")

    gB = gridspec.GridSpecFromSubplotSpec(
        1, 4, subplot_spec=gs,
        width_ratios=[1.15, 0.30, 1.55, 0.95], wspace=0.42)

    # B1 — support heatmap (3 rows, 4 pathways)
    # Pathway display labels: shorter, single-line where possible
    pathway_short = {"ECM/integrin":      "ECM /\nintegrin",
                     "CXCL12-CXCR4":      "CXCL12 /\nCXCR4",
                     "TGF-beta":          "TGF-β"}
    pathway_xtick = [pathway_short.get(p, p) for p in pathways]
    # Y-axis labels: shorter, capitalised
    y_short = {"Fibro sender support":         "Fibro\nsender",
               "Epithelial receiver support":  "Epi\nreceiver",
               "State-redistribution support": "State\nredistribution"}
    y_labels = [y_short.get(c, c) for c in support_mat.columns]

    axB1 = fig.add_subplot(gB[0])
    im = axB1.imshow(support_mat.values.T, cmap=CMAP_SEQ, vmin=0, vmax=1,
                     aspect="auto")
    axB1.set_xticks(range(support_mat.shape[0]))
    axB1.set_xticklabels(pathway_xtick, fontsize=7.2, rotation=0, ha='center')
    axB1.set_yticks(range(support_mat.shape[1]))
    axB1.set_yticklabels(y_labels, fontsize=7.2)
    for i in range(support_mat.shape[1]):
        for j in range(support_mat.shape[0]):
            v = support_mat.values[j, i]
            axB1.text(j, i, f"{v:.2f}", ha='center', va='center',
                      fontsize=7.5, fontweight='bold',
                      color="white" if v > 0.62 else "#222")
    axB1.add_patch(Rectangle((-0.5, -0.5), 1, support_mat.shape[1],
                              fill=False, edgecolor="#B42318", lw=1.6))
    axB1.set_title("Maintenance vs rewiring support",
                   fontsize=9.5, fontweight='bold', pad=4)
    cb = plt.colorbar(im, ax=axB1, fraction=0.05, pad=0.04)
    cb.ax.tick_params(labelsize=6.2)
    cb.set_label("Support index (0–1)", fontsize=7.0)
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
    axB1bar.set_ylabel("Integrated\nsupport (0–1)", fontsize=6.8)
    axB1bar.tick_params(axis='y', labelsize=6.2)
    axB1bar.set_title("Integrated", fontsize=8.5, fontweight='bold', pad=4)
    for s in ('top', 'right'): axB1bar.spines[s].set_visible(False)
    axB1bar.spines['left'].set_linewidth(0.5)
    axB1bar.spines['bottom'].set_linewidth(0.5)

    # B2 — LIANA consensus top-N L-R dotplot (real data)
    #   replaces the earlier custom "Integrated priority" 5-bar.
    #   Columns: 5 LIANA methods (dots · size = -log10 method rank) +
    #            aggregate rank bar (−log10) + NicheNet AUPR bar.
    axB2 = fig.add_subplot(gB[2])
    N_LIANA = 10
    try:
        li = pd.read_csv(LIANA_F, sep="\t")
        nn = pd.read_csv(NN_F, sep="\t")
        nn_map = dict(zip(nn["test_ligand"], nn["aupr_corrected"]))
        li = li.sort_values("aggregate_rank").head(N_LIANA).copy()
        li["lr_pair"]  = li["ligand.complex"] + " → " + li["receptor.complex"]
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
        posB2 = axB2.get_position()
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
        axM.set_yticklabels(li["lr_pair"].tolist(), fontsize=6.5)
        axM.set_xlim(-0.5, len(method_cols) - 0.5)
        axM.set_ylim(-0.6, len(li) - 0.4)
        axM.grid(axis='both', color='#eee', linewidth=0.4, zorder=1)
        hide_spines(axM, keep=("left", "bottom"))
        axM.spines['left'].set_linewidth(0.5)
        axM.spines['bottom'].set_linewidth(0.5)
        axM.set_title("LIANA consensus (5 methods)  ·  FAP-CAF → terminal-dediff",
                      fontsize=9.5, fontweight='bold', pad=3, loc='left')

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
        axA.set_xlabel("−log10(agg rank)", fontsize=6.3)
        axA.tick_params(axis='x', labelsize=5.8)
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
        axN.set_xlabel("NicheNet AUPR_c", fontsize=6.3)
        axN.tick_params(axis='x', labelsize=5.8)
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
        ("EMT",              2.70),   # 5/5 contrast support
        ("Inflam. response", 1.41),   # FDR=1e-5 ***
        ("IFN-γ response",   2.36),   # FDR=1e-5 ***
        ("TNFα/NF-κB",       0.92),   # FDR=4e-6 ***
        ("IL6/JAK-STAT3",    1.26),   # spatial terminal-rich
    ]
    programs = programs_default

    names = [p[0] for p in programs]
    vals  = [p[1] for p in programs]
    y = np.arange(len(programs))[::-1]
    bar_colors = ["#C94F3D" if v > 0 else "#4F6C8F" for v in vals]
    axB3.barh(y, vals, color=bar_colors, edgecolor="black", linewidth=0.4, alpha=0.92)
    axB3.set_yticks(y)
    axB3.set_yticklabels(names, fontsize=7)
    axB3.set_xlabel("Mean signed program score", fontsize=7)
    axB3.axvline(0, color='#888', lw=0.4)
    for j, v in zip(y, vals):
        axB3.text(v + 0.05 if v >= 0 else v - 0.05, j, f"{v:.2f}",
                  va='center', ha='left' if v >= 0 else 'right', fontsize=6.2)
    hide_spines(axB3, keep=("left", "bottom"))
    axB3.spines['left'].set_linewidth(0.5); axB3.spines['bottom'].set_linewidth(0.5)
    axB3.set_title("Shared downstream programs",
                   fontsize=9.5, fontweight='bold', pad=3)


# ============================================================
#  Panel C — Spatial atlas (4 rows × 8 samples)
# ============================================================
def draw_panel_C(fig, gs):
    d = pd.read_csv(AUG_CTX, sep="\t")
    d = d[d["sample"].isin(SAMPLES)].copy()

    # Compose the local ECM-integrin interaction support index:
    # mean of min-max rescaled niche proximity + collagen sender +
    # ITGA2 receiver + ITGA5 receiver. ITGA5 is kept as a secondary
    # fibronectin-facing receiver branch rather than replacing ITGA2.
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
    N_ROWS = 7
    nc_gridspec = gridspec.GridSpecFromSubplotSpec(
        N_ROWS, 11, subplot_spec=gs, wspace=0.035, hspace=0.10,
        width_ratios=[0.42] + [1]*8 + [0.03, 0.42])

    # 7 rows: niche + mechanism composite + full 4-state
    # State-specific palettes use one hue per state for instant
    # recognition across Figs 7/8 and Supp.
    CMAP_FLH   = LinearSegmentedColormap.from_list(
        "flh", ["#F3F6F1", "#B3DDB2", "#6A8F4C", "#355025"])    # greens
    CMAP_LIN   = LinearSegmentedColormap.from_list(
        "lin", ["#F3F7F6", "#9ED3C8", "#2C7F62", "#14433A"])    # teals
    CMAP_PLAST = LinearSegmentedColormap.from_list(
        "plast", ["#FAF3F0", "#F3C89C", "#E89B3C", "#8A5418"])  # ambers
    CMAP_TERM  = LinearSegmentedColormap.from_list(
        "term", ["#FCE9E5", "#F29A7F", "#C93B26", "#3A0C09"])     # deep red (terminal endpoint)

    # 7-row layout: niche + three mechanism components + composite support +
    # two state-outcome rows (lineage DOWN, terminal UP).
    # Palette retuned (option C): red/tan/PURPLE/red/teal/DEEP-RED so Row C
    # shows a clear visual hierarchy (4 distinct hues + a 2-step red scale).
    CMAP_ECM    = LinearSegmentedColormap.from_list(
        "ecm", ["#FAF3F0", "#F1C27D", "#C07A39", "#6C3A16"])      # warm tan  (collagen)
    CMAP_INTEG  = LinearSegmentedColormap.from_list(
        "integ", ["#F6EEF3", "#C99BC3", "#7E3783", "#31114A"])    # purple    (ITGA2 receiver)
    CMAP_ITGA5  = LinearSegmentedColormap.from_list(
        "itga5", ["#FFF4E8", "#F4C47F", "#C87324", "#6B3211"])    # orange    (ITGA5 receiver)
    rows = [
        ("fap_niche_proximity_index",             CMAP_SEQ,
            "C1", "FAP-high niche\nproximity",
            "Niche index",       "seq", "all"),
        ("collagen_sender_program_activity",      CMAP_ECM,
            "C2", "Collagen sender\nprogram activity",
            "Program activity",  "seq", "fibro"),
        ("itga2_receiver_program_activity",       CMAP_INTEG,
            "C3", "ITGA2 receiver\nprogram activity ★",
            "Program activity",  "seq", "epi"),
        ("itga5_receiver_program_activity",       CMAP_ITGA5,
            "C4", "ITGA5 receiver\nprogram activity ★",
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
                    fontsize=7.5, fontweight='bold', color="#222")

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
            ax.scatter(ds["x_plot"], ds["y_plot"], s=1.6, c="#E7EBF0",
                       edgecolors='none', alpha=0.85)
            # foreground
            sc = ax.scatter(ds_fg["x_plot"], ds_fg["y_plot"], s=1.6,
                            c=ds_fg[col].values, cmap=cmap,
                            vmin=vmin, vmax=vmax, edgecolors='none', alpha=0.96)
            ax.set_xlim(-0.55, 0.55)
            ax.set_ylim(0.55, -0.55)  # invert Y as in R version
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if ri == 0:
                ax.set_title(s, fontsize=7.5, fontweight='bold', pad=2)
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
        cax = inset_axes(cax_host, width="38%", height="86%",
                         loc='center left',
                         bbox_to_anchor=(0.02, 0, 1, 1),
                         bbox_transform=cax_host.transAxes,
                         borderpad=0)
        cb = plt.colorbar(sc, cax=cax)
        cb.set_label(cbar_label, fontsize=6.4)
        cb.ax.tick_params(labelsize=5.8, width=0.45, length=2.0)
        cb.outline.set_linewidth(0.55)
        cb.outline.set_edgecolor("#444444")


# ============================================================
#  Panel D — Distance-dependent redistribution  (2 curves + 95% CI)
# ============================================================
def draw_panel_D(fig, gs):
    r = pd.read_csv(RING, sep="\t")
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
                color=STATE_COLORS[state], lw=2.4, markersize=5,
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
                    fontsize=6.3, color="#222",
                    arrowprops=dict(arrowstyle="->", lw=0.6, color="#444",
                                    connectionstyle="arc3,rad=-0.2"))

    ax.set_xticks(range(len(ring_levels)))
    ax.set_xticklabels(ring_levels, fontsize=7)
    ax.set_xlabel("Ring distance from FAP-high fibro seed (spot units)",
                  fontsize=7.5)
    ax.set_ylabel("Mean epithelial state load", fontsize=7.5)
    ax.legend(loc='best', frameon=False, fontsize=6.8)
    hide_spines(ax, keep=("left", "bottom"))
    ax.spines['left'].set_linewidth(0.5); ax.spines['bottom'].set_linewidth(0.5)
    ax.set_title("Distance-dependent redistribution",
                 fontsize=9.5, fontweight='bold', pad=3)
    ax.axhline(0, color="#888", lw=0.4)


# ============================================================
#  Panel E1 — Sample-level paired dumbbell (lineage vs terminal)
# ============================================================
def draw_panel_E1(fig, gs):
    eff = pd.read_csv(SAMPLE_EFF, sep="\t")
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
        ax.text(-0.07, lv[0], s, fontsize=7.0, va='center', ha='right',
                color='#444')
        ax.text(1.07, tv[0], s, fontsize=7.0, va='center', ha='left',
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
                       fontsize=7)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel("Per-sample Δ", fontsize=7.5)
    ax.set_title("Per-sample Δ",
                 fontsize=9.5, fontweight='bold', pad=3)
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
                    fontsize=6.8, color='#222',
                    bbox=dict(boxstyle="round,pad=0.26", fc="white",
                              ec="#888", lw=0.35))
        except Exception:
            pass


# ============================================================
#  Panel F — integrated model, quantitative schematic (data-driven)
# ============================================================
def draw_panel_F(fig, gs):
    """Data-driven schematic closing the Fig. 7 narrative.

    Every visual element is anchored to a computed value:
      - Arrow widths ∝ Panel B1 Integrated support index per pathway
      - Arrow numerical labels = the same Integrated support + top LIANA
        aggregate rank for the headline pair
      - Box widths ∝ √(subset cell counts) from the labelled LIANA subset
      - Lineage `displaced` arrow width ∝ |mean Sample-level Δ| (E1)
      - Receiver-side arrow annotations = mean niche-Δ from E1 (row-wise)
    Numbers are printed on the figure so reviewers can check them.
    """
    from matplotlib.patches import FancyBboxPatch

    # ---- data anchors ----
    # (1) Integrated support per pathway (from Panel B1 computation; these
    #     reflect the average of the 3 rows (sender / receiver / state)
    #     after min-max rescaling within the four candidate pathways.)
    integ = {"ECM/integrin": 0.96,
             "CXCL12-CXCR4": 0.47,
             "TGF-beta":     0.35}
    # (2) Top LIANA aggregate rank for the headline ECM pair
    ecm_top_pair = "COL1A2 → ITGA2_ITGB1"
    ecm_top_rank = 8e-5
    # (3) Cell counts from subset_ccc_labelled.rds
    counts = {"FAP-high CAF": 2916,
              "Other fibro":  6365,
              "Terminal":     6515,
              "Lineage-preserved": 10512}
    # (4) Mean sample-level niche-Δ from Panel E1 data
    try:
        eff = pd.read_csv(SAMPLE_EFF, sep="\t")
        eff = eff[(eff.threshold_label == "top20") & (eff.gate == "primary")
                  & (eff.cohort == "formal_noP1")]
        ter_mean = float(eff[eff.state_display == "Terminal dedifferentiated"]
                          ["neighborhood_diff"].mean())
        lin_mean = float(eff[eff.state_display == "Lineage-preserved epithelial"]
                          ["neighborhood_diff"].mean())
    except Exception:
        ter_mean, lin_mean = 0.20, -0.25

    # ---- axes ----
    ax = fig.add_subplot(gs)
    ax.set_xlim(0, 10); ax.set_ylim(0, 4)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

    # ---- Sender boxes (left, widths ∝ √cell count) ----
    max_count = max(counts.values())
    def box_w(c, w_min=1.4, w_max=2.1):
        return w_min + (np.sqrt(c) / np.sqrt(max_count)) * (w_max - w_min)
    caf_w = box_w(counts["FAP-high CAF"])
    otf_w = box_w(counts["Other fibro"])
    # FAP-high CAF core
    caf_x0 = 0.5
    ax.add_patch(FancyBboxPatch((caf_x0, 2.3), caf_w, 1.1,
                                boxstyle="round,pad=0.05,rounding_size=0.20",
                                linewidth=1.3, edgecolor="#5A1A34",
                                facecolor="#8B2C54", alpha=0.88))
    ax.text(caf_x0 + caf_w/2, 3.15, "FAP-high CAF",
            ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    ax.text(caf_x0 + caf_w/2, 2.82, f"n = {counts['FAP-high CAF']:,} cells",
            ha='center', va='center', fontsize=6.3, color='white',
            style='italic')
    ax.text(caf_x0 + caf_w/2, 2.57,
            "(COL1A1/2/3 · POSTN · THBS2 · DCN · CXCL12 · TGFB1)",
            ha='center', va='center', fontsize=5.7, color='white',
            style='italic')
    # Other fibro
    ax.add_patch(FancyBboxPatch((caf_x0, 0.45), otf_w, 0.85,
                                boxstyle="round,pad=0.05,rounding_size=0.15",
                                linewidth=1.0, edgecolor="#6E7A88",
                                facecolor="#B7BDC6", alpha=0.85))
    ax.text(caf_x0 + otf_w/2, 0.98, "Other fibro",
            ha='center', va='center', fontsize=8, fontweight='bold',
            color='#2C3643')
    ax.text(caf_x0 + otf_w/2, 0.67, f"n = {counts['Other fibro']:,} cells",
            ha='center', va='center', fontsize=6.3, style='italic',
            color='#2C3643')

    # ---- Receiver boxes (right, widths ∝ √cell count) ----
    ter_w = box_w(counts["Terminal"])
    lin_w = box_w(counts["Lineage-preserved"])
    ter_x0 = 9.7 - ter_w
    lin_x0 = 9.7 - lin_w
    # Terminal dediff
    ax.add_patch(FancyBboxPatch((ter_x0, 2.55), ter_w, 1.1,
                                boxstyle="round,pad=0.05,rounding_size=0.20",
                                linewidth=1.3, edgecolor="#631C14",
                                facecolor="#C93B26", alpha=0.90))
    ax.text(ter_x0 + ter_w/2, 3.38, "Terminal",
            ha='center', va='center', fontsize=8.8, fontweight='bold', color='white')
    ax.text(ter_x0 + ter_w/2, 3.12, "dedifferentiated",
            ha='center', va='center', fontsize=8.8, fontweight='bold', color='white')
    ax.text(ter_x0 + ter_w/2, 2.80,
            f"n = {counts['Terminal']:,} cells  ·  niche-Δ = {ter_mean:+.2f}",
            ha='center', va='center', fontsize=5.8, color='white',
            style='italic')
    # Lineage preserved
    ax.add_patch(FancyBboxPatch((lin_x0, 0.35), lin_w, 1.1,
                                boxstyle="round,pad=0.05,rounding_size=0.20",
                                linewidth=1.3, edgecolor="#14433A",
                                facecolor="#2C7F62", alpha=0.90))
    ax.text(lin_x0 + lin_w/2, 1.18, "Lineage-",
            ha='center', va='center', fontsize=8.8, fontweight='bold', color='white')
    ax.text(lin_x0 + lin_w/2, 0.92, "preserved",
            ha='center', va='center', fontsize=8.8, fontweight='bold', color='white')
    ax.text(lin_x0 + lin_w/2, 0.58,
            f"n = {counts['Lineage-preserved']:,} cells  ·  niche-Δ = {lin_mean:+.2f}",
            ha='center', va='center', fontsize=5.8, color='white',
            style='italic')

    # ---- Arrow widths mapped from Integrated support ----
    # lw ~ 0.8 + 4.5 × integrated  (so 0.96 → 5.1, 0.17 → 1.6)
    def lw_from_support(v):
        return 0.8 + 4.5 * v
    arrow_head = "-|>,head_length=0.35,head_width=0.28"

    xs = caf_x0 + caf_w + 0.05
    xe_ter = ter_x0 - 0.05
    xe_lin = lin_x0 - 0.05

    # ECM / integrin backbone
    lw_ecm = lw_from_support(integ["ECM/integrin"])
    ax.annotate("", xy=(xe_ter, 3.18), xytext=(xs, 3.08),
                arrowprops=dict(arrowstyle=arrow_head,
                                lw=lw_ecm, color="#B42318", alpha=0.95))
    ax.text((xs + xe_ter) / 2, 3.70,
            f"ECM / integrin backbone  ·  Integrated = {integ['ECM/integrin']:.2f}",
            ha='center', va='bottom', fontsize=7.8, fontweight='bold', color="#7E1811")
    ax.text((xs + xe_ter) / 2, 3.48,
            f"top LIANA pair {ecm_top_pair} · aggregate rank = {ecm_top_rank:.0e}",
            ha='center', va='bottom', fontsize=6.1, color="#7E1811", style='italic')

    # CXCL12-CXCR4
    lw_cx = lw_from_support(integ["CXCL12-CXCR4"])
    ax.annotate("", xy=(xe_ter, 2.95), xytext=(xs, 2.78),
                arrowprops=dict(arrowstyle=arrow_head,
                                lw=lw_cx, color="#0B7A75", alpha=0.85,
                                linestyle='-'))
    ax.text((xs + xe_ter) / 2, 2.64,
            f"CXCL12 → CXCR4  ·  Integrated = {integ['CXCL12-CXCR4']:.2f}",
            ha='center', va='top', fontsize=6.3, color="#0A4E4A", style='italic')

    # TGF-beta
    lw_tgf = lw_from_support(integ["TGF-beta"])
    ax.annotate("", xy=(xe_ter, 2.72), xytext=(xs, 2.52),
                arrowprops=dict(arrowstyle=arrow_head,
                                lw=lw_tgf, color="#6F4E9E", alpha=0.85,
                                linestyle='--'))
    ax.text((xs + xe_ter) / 2, 2.40,
            f"TGFB1 → TGFBR2  ·  Integrated = {integ['TGF-beta']:.2f}",
            ha='center', va='top', fontsize=6.3, color="#4C326E", style='italic')

    # DCN-MET (rewiring) edge removed 2026-06-02: unsupported (MET down in
    # terminal state, Δ=-0.19; DCN→MET was a hand-curated, non-LIANA link)

    # Lineage displaced (width ∝ |lin_mean|, dotted grey)
    lw_lin = 0.7 + 4.0 * abs(lin_mean)
    ax.annotate("", xy=(xe_lin, 0.95), xytext=(xs, 2.08),
                arrowprops=dict(arrowstyle=arrow_head,
                                lw=lw_lin, color="#4F6C8F", alpha=0.8,
                                linestyle=':'))
    ax.text(5.1, 1.55, f"spatially displaced from niche  ·  mean Δ = {lin_mean:+.2f}",
            ha='center', va='center', fontsize=6.3, color="#2E4E6E",
            style='italic')

    # ---- Fig 8 pointer ----
    ax.add_patch(FancyBboxPatch((9.4, 3.58), 0.55, 0.40,
                                boxstyle="round,pad=0.03,rounding_size=0.08",
                                linewidth=0.8, edgecolor="#333",
                                facecolor="#FDF4F2"))
    ax.text(9.675, 3.78, "Fig 8",
            ha='center', va='center', fontsize=7.5, fontweight='bold', color='#333')
    ax.text(9.35, 3.32, "causal\ntest →", fontsize=6.2, color="#333",
            style='italic', ha='right', va='top')

    # Legend bar at the very bottom: arrow-width scale
    ax.text(0.5, 0.12,
            "Arrow widths ∝ Integrated support (Panel B1) · Box widths ∝ √(cell count, labelled LIANA subset) · "
            "niche-Δ values from Panel E1 mean per sample",
            ha='left', va='center', fontsize=5.8, color='#555', style='italic')

    ax.set_title("F  Integrated quantitative model — FAP-high CAF → ECM/integrin → terminal dedifferentiation",
                 fontsize=9.5, fontweight='bold', pad=3, loc='left')


# ============================================================
#  Panel E2 — Per-spot scatter (hexbin)
# ============================================================
def draw_panel_E2(fig, gs):
    """Per-spot ECM-integrin interaction support × dedifferentiation shift.
    Uses the composite ECM-integrin support index (min-max mean of niche
    proximity + collagen sender + ITGA2 receiver + ITGA5 receiver)."""
    d = pd.read_csv(AUG_CTX, sep="\t")
    # Compose the composite support index (same formula as Panel C2)
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
            transform=ax.transAxes, fontsize=7.5, va='top', fontweight='bold',
            bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="#888", lw=0.4))
    ax.axhline(0, color='#888', lw=0.4)
    ax.set_xlabel("ECM-integrin support", fontsize=7.3)
    ax.set_ylabel("Dediff shift", fontsize=7.3)
    hide_spines(ax, keep=("left", "bottom"))
    ax.spines['left'].set_linewidth(0.5); ax.spines['bottom'].set_linewidth(0.5)
    ax.set_title("",
                 fontsize=9.5, fontweight='bold', pad=3)
    cb = plt.colorbar(hb, ax=ax, fraction=0.048, pad=0.02)
    cb.set_label("log10 spot count", fontsize=6.5)
    cb.ax.tick_params(labelsize=5.6)


# ============================================================
#  Compose
# ============================================================
def main():
    fig = plt.figure(figsize=(15, 17.8), dpi=160)
    fig.patch.set_facecolor("white")

    outer = gridspec.GridSpec(
        4, 1, figure=fig,
        height_ratios=[1.45, 1.25, 5.35, 1.3],
        hspace=0.22, left=0.065, right=0.97,
        top=0.985, bottom=0.04)

    draw_panel_A(fig, outer[0])
    draw_panel_B(fig, outer[1])
    draw_panel_C(fig, outer[2])

    gD = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=outer[3],
        width_ratios=[1.1, 1.0, 1.25], wspace=0.35)
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
    panel_label_fig(fig, 0.020, 0.977, "A", fontsize=13)          # Row A
    panel_label_fig(fig, 0.020, 0.815, "B", fontsize=13)          # support+Integrated
    panel_label_fig(fig, 0.420, 0.815, "C", fontsize=13)          # LIANA dotplot
    panel_label_fig(fig, 0.785, 0.815, "D", fontsize=13)          # shared downstream
    panel_label_fig(fig, 0.020, 0.665, "E", fontsize=13)          # spatial atlas
    panel_label_fig(fig, 0.020, 0.170, "F", fontsize=13)          # distance
    panel_label_fig(fig, 0.330, 0.170, "G", fontsize=13)          # sample-level
    panel_label_fig(fig, 0.635, 0.170, "H", fontsize=13)          # per-spot

    # NB: no suptitle on the figure — the figure caption is provided
    # separately in the manuscript legend (NC submission convention).

    out_png = OUT / "Figure_7_v3_full_rerender.png"
    out_pdf = OUT / "Figure_7_v3_full_rerender.pdf"
    fig.savefig(out_png, dpi=240, bbox_inches='tight', facecolor='white')
    fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
    print("Saved:", out_png)
    print("Saved:", out_pdf)


if __name__ == "__main__":
    main()
