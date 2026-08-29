# Cancer Research submission - figure code release
# Builds: Figure 3A (paired sender-receiver nomination). Helper module
#          imported by 02_ligand_receptor_nomination.py.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

"""
Fig 7 Panel A — paired sender ↔ receiver L-R nomination (heavy-duty NC version)

Left block  = FAP-high CAF sender ligands
Right block = terminal-dediff receiver receptors
Middle      = curved connection lines, one per LIANA top-ranked ligand-
              receptor pair. Line thickness ∝ -log10(LIANA aggregate_rank);
              line colour = pathway family (ECM-integrin red,
              CXCL12-CXCR4 teal, TGF-β violet). DCN-MET axis dropped
              2026-06-02 (MET down in terminal state; DCN recoloured to ECM).

Each gene row carries:
  · A horizontal Δ-bar (log-normalized expression difference vs the
    paired contrast, from `A_sender_ligand_deltas.tsv` or
    `A_receiver_receptor_deltas.tsv`)
  · Asterisks encoding Wilcoxon BH-adjusted p-value
  · A dot indicating % expressing in the on-side group

This is a standalone panel for preview and for later merging into
Fig 7 v3 Panel A via its draw_panel_A() replacement.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.path import Path as MPLPath
from matplotlib.patches import PathPatch

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size":   8,
    "axes.linewidth": 0.55,
})

# ---------- paths (env-var driven; see README) ----------
EXTERNAL_DATA = Path(os.environ.get("EXTERNAL_DATA", "<EXTERNAL_DATA>"))
PROJECT_ROOT  = Path(os.environ.get("PROJECT_ROOT",  "<PROJECT_ROOT>"))
SRC_TBL = EXTERNAL_DATA / "fig7_liana_nichenet" / "tables"
SENDER = SRC_TBL / "A_sender_ligand_deltas.tsv"
RECV   = SRC_TBL / "A_receiver_receptor_deltas.tsv"
LIANA  = SRC_TBL / "liana_focused_FAP_to_terminal_top50.tsv"
OUT    = PROJECT_ROOT / "outputs" / "fig7_redesign"

# ---------- palette ----------
PATHWAY_COLORS = {       # harmonized to paper_A `pa` palette 2026-06-02
    "ECM/integrin":  "#C0392B",   # pa canonical "invasion/high" red (unifies all ECM reds figure-wide)
    "CXCL12-CXCR4":  "#0B7A75",   # teal
    "TGF-beta":      "#6A3D9A",   # pa intermediate purple
    "Other":         "#9E9E9E",
}
UP_COLOR   = "#C0392B"    # delta > 0  (unified red)
DOWN_COLOR = "#2C6FBB"    # delta < 0  (pa identity blue)
MAINT_BAND = "#F5E6E3"
REWIR_BAND = "#F1ECF7"

# ---------- gene plan ----------
SENDER_MAINT = ["COL1A2", "COL1A1", "COL3A1", "COL6A1", "POSTN", "THBS2", "DCN"]
SENDER_REWIR = ["CXCL12", "TGFB1"]
RECV_MAINT   = ["ITGA2", "ITGB1", "ITGA5", "ITGA3", "CD44", "SDC1", "DDR1"]
RECV_REWIR   = ["CXCR4", "TGFBR2"]

# Curated connections (always drawn; LIANA-derived ones overlaid by rank)
# ITGA5 link added to surface the receptor that ranks #11 in protein-level
# Fig 8 v13 ECM/focal-adhesion sentinel panel and #426/9534 in cross-cohort
# proteomics bridge score (highest among integrins).
CURATED_EDGES = [
    ("CXCL12", "CXCR4",  "CXCL12-CXCR4"),
    ("TGFB1",  "TGFBR2", "TGF-beta"),
    ("POSTN",  "ITGB1",  "ECM/integrin"),
    ("POSTN",  "ITGA3",  "ECM/integrin"),
    ("THBS2",  "ITGA2",  "ECM/integrin"),
    ("COL6A1", "ITGA3",  "ECM/integrin"),
    ("COL1A1", "ITGA5",  "ECM/integrin"),    # LIANA rank 17; cross-figure link to Fig 8 panel g
]


def sig_stars(p):
    if pd.isna(p):     return ""
    if p < 1e-10:      return "***"
    if p < 1e-3:       return "**"
    if p < 0.05:       return "*"
    return ""


def classify_edge(lig, rec):
    if lig == "CXCL12" and "CXCR" in rec:           return "CXCL12-CXCR4"
    if lig.startswith("TGFB") and rec.startswith("TGFBR"): return "TGF-beta"
    return "ECM/integrin"


def parse_liana_edges(liana_df, sender_genes, receiver_genes, top_n=12):
    """Return a list of (ligand, receptor, pathway, strength_0_to_1)."""
    edges = []
    df = liana_df.sort_values("aggregate_rank").copy()
    # Expand receptor complexes like ITGA2_ITGB1 into individual receptors
    for _, r in df.head(40).iterrows():  # scan a bit deeper to catch complexes
        lig = str(r["ligand.complex"])
        rec_c = str(r["receptor.complex"])
        if lig not in sender_genes:
            continue
        parts = rec_c.split("_")
        for p in parts:
            if p in receiver_genes:
                edges.append((lig, p, classify_edge(lig, p),
                              r["aggregate_rank"]))
    # dedup (keep best rank)
    best = {}
    for lig, rec, pw, rank in edges:
        key = (lig, rec)
        if key not in best or rank < best[key][3]:
            best[key] = (lig, rec, pw, rank)
    edges = list(best.values())[:top_n]
    # convert rank → strength (for line width / alpha)
    if len(edges):
        nlogs = [-np.log10(max(e[3], 1e-7)) for e in edges]
        mn, mx = min(nlogs), max(nlogs)
        def norm(v): return (v - mn) / max(mx - mn, 1e-9)
        edges = [(l, r, pw, norm(-np.log10(max(rk, 1e-7))))
                 for (l, r, pw, rk) in edges]
    return edges


def draw_connection(ax, y_left, y_right, color, strength,
                    x_left=0.03, x_right=0.97, linestyle="solid"):
    """Bezier curve between (x_left, y_left) and (x_right, y_right).
    linestyle = 'solid' for LIANA-supported edges, 'dashed' for curated-only."""
    if linestyle == "solid":
        alpha = 0.45 + 0.50 * strength
        lw    = 0.8  + 2.4  * strength
    else:   # dashed (curated)
        alpha = 0.55
        lw    = 0.9
    x_ctrl1 = x_left + (x_right - x_left) * 0.35
    x_ctrl2 = x_left + (x_right - x_left) * 0.65
    verts = [(x_left,  y_left),
             (x_ctrl1, y_left),
             (x_ctrl2, y_right),
             (x_right, y_right)]
    codes = [MPLPath.MOVETO, MPLPath.CURVE4, MPLPath.CURVE4, MPLPath.CURVE4]
    path = MPLPath(verts, codes)
    patch = PathPatch(path, facecolor="none",
                      edgecolor=color, lw=lw, alpha=alpha, zorder=3,
                      linestyle=("--" if linestyle == "dashed" else "-"))
    ax.add_patch(patch)


def draw_sidebar(ax, df, genes, side="left", title=""):
    """Draw: Δ-bar (with sign colour), p-stars, %expr dot."""
    df = df.set_index("gene").reindex(genes).reset_index()
    y = np.arange(len(genes))[::-1]   # top gene at top
    delta = df["delta_log_norm"].values
    pct   = df["pct_grp1"].values
    padj  = df["p_adj_BH"].values

    xmin = min(0, delta.min() * 1.1)
    xmax = max(0, delta.max() * 1.15)

    for i, (d, p, pc, gene) in enumerate(zip(delta, padj, pct, df["gene"])):
        yi = y[i]
        col = UP_COLOR if d >= 0 else DOWN_COLOR
        if side == "left":
            ax.barh(yi, d, color=col, edgecolor="black",
                    linewidth=0.35, height=0.55, zorder=2)
        else:
            ax.barh(yi, d, color=col, edgecolor="black",
                    linewidth=0.35, height=0.55, zorder=2)
        # significance stars
        stars = sig_stars(p)
        if stars:
            tx = d + (0.05 if d >= 0 else -0.05)
            ha = "left" if d >= 0 else "right"
            ax.text(tx, yi, stars, fontsize=7, fontweight="bold",
                    va="center", ha=ha, color="#333", zorder=4)
        # %expr dot (size)
        dot_x = xmax * 1.05 if side == "left" else xmin * 1.05
        ax.scatter(dot_x, yi,
                   s=35 + 230 * pc, facecolor="#4A5568",
                   edgecolor="white", linewidth=0.3,
                   clip_on=False, zorder=3)

    ax.axvline(0, color="black", linewidth=0.5, zorder=1)
    ax.set_xlim(xmin, xmax * 1.20)
    ax.set_ylim(-0.5, len(genes) - 0.5)
    ax.set_yticks([])
    ax.tick_params(axis='x', labelsize=6.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines['left'].set_linewidth(0.4)
    ax.spines['bottom'].set_linewidth(0.4)
    ax.set_xlabel("Δ  log-norm  (" + title + ")", fontsize=6.8)


def draw_side_labels(ax, genes, groups, pathway_map):
    """Gene-name column; background band by maintenance/rewiring."""
    y = np.arange(len(genes))[::-1]
    ax.set_xlim(0, 1); ax.set_ylim(-0.5, len(genes) - 0.5)
    # background bands
    prev_group = None
    band_start = None
    for i, g in enumerate(groups):
        if g != prev_group:
            if prev_group is not None and band_start is not None:
                band_color = MAINT_BAND if prev_group == "maintenance" else REWIR_BAND
                ax.axhspan(y[band_start] + 0.5,
                           y[i-1] - 0.5,
                           xmin=0, xmax=1,
                           color=band_color, zorder=0)
            band_start = i
            prev_group = g
    if prev_group is not None:
        band_color = MAINT_BAND if prev_group == "maintenance" else REWIR_BAND
        ax.axhspan(y[band_start] + 0.5, y[len(genes)-1] - 0.5,
                   xmin=0, xmax=1, color=band_color, zorder=0)

    for yi, gene, group in zip(y, genes, groups):
        fam = pathway_map.get(gene, "Other")
        c   = PATHWAY_COLORS.get(fam, "#4A5568")
        ax.text(0.5, yi, gene, ha='center', va='center',
                fontsize=8.5, fontweight='bold', color=c, zorder=4)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)


def sender_pathway(gene):
    if gene in ["COL1A1","COL1A2","COL3A1","COL6A1","POSTN","THBS2","DCN"]:
        return "ECM/integrin"   # DCN recoloured to ECM (CAF-up ECM proteoglycan); DCN-MET axis dropped 2026-06-02
    if gene == "CXCL12":  return "CXCL12-CXCR4"
    if gene == "TGFB1":   return "TGF-beta"
    return "Other"

def recv_pathway(gene):
    if gene in ["ITGA2","ITGB1","ITGA5","ITGA3","CD44","SDC1","DDR1"]:
        return "ECM/integrin"
    if gene == "CXCR4":   return "CXCL12-CXCR4"   # MET row + DCN-MET axis removed 2026-06-02 (MET down in terminal state)
    if gene == "TGFBR2":  return "TGF-beta"
    return "Other"


def main():
    sender = pd.read_csv(SENDER, sep="\t")
    recv   = pd.read_csv(RECV,   sep="\t")
    liana  = pd.read_csv(LIANA,  sep="\t")

    # restrict to genes we plot
    sender_genes = SENDER_MAINT + SENDER_REWIR
    recv_genes   = RECV_MAINT   + RECV_REWIR
    sender = sender[sender["gene"].isin(sender_genes)]
    recv   = recv[recv["gene"].isin(recv_genes)]

    # LIANA edges (plus curated)
    liana_edges = parse_liana_edges(liana, sender_genes, recv_genes, top_n=12)
    liana_set = {(l, r) for (l, r, _, _) in liana_edges}
    # Add curated edges that are NOT in LIANA top — drawn dashed
    curated_extra = [(l, r, pw, 0.30) for (l, r, pw) in CURATED_EDGES
                     if (l, r) not in liana_set
                     and l in sender_genes and r in recv_genes]
    # Tag each edge with its source so we can set linestyle
    # edge tuples become (lig, rec, pathway, strength, source)
    all_edges = ([(l, r, pw, s, "LIANA")   for (l, r, pw, s) in liana_edges]
               + [(l, r, pw, s, "curated") for (l, r, pw, s) in curated_extra])

    # Layout
    fig = plt.figure(figsize=(13, 5.6), dpi=180)
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(1, 5, figure=fig,
                           width_ratios=[1.1, 1.15, 2.40, 1.15, 1.1],
                           wspace=0.05, left=0.04, right=0.99,
                           top=0.87, bottom=0.18)

    # Left sender Δ bar (flipped so bars go left)
    ax_sL = fig.add_subplot(gs[0])
    sender_df = sender.set_index("gene").reindex(sender_genes).reset_index()
    y = np.arange(len(sender_genes))[::-1]
    delta_s = sender_df["delta_log_norm"].values
    pct_s   = sender_df["pct_grp1"].values
    padj_s  = sender_df["p_adj_BH"].values
    # sender bar drawn to the LEFT (so axis is inverted)
    for i, (d, p, pc) in enumerate(zip(delta_s, padj_s, pct_s)):
        yi = y[i]
        col = UP_COLOR if d >= 0 else DOWN_COLOR
        ax_sL.barh(yi, -d, color=col, edgecolor="black",
                   linewidth=0.35, height=0.55, zorder=2)
        stars = sig_stars(p)
        if stars:
            ax_sL.text(-d - 0.08, yi, stars, fontsize=7, fontweight="bold",
                        va='center', ha='right', color="#333", zorder=4)
    # %expr dot (far left edge)
    xlim_s = max(abs(delta_s)) * 1.4
    for i, pc in enumerate(pct_s):
        ax_sL.scatter(-xlim_s*1.03, y[i], s=35 + 230 * pc,
                       facecolor="#E2A736", edgecolor="#333",
                       linewidth=0.35, clip_on=False, zorder=3)
    ax_sL.axvline(0, color="black", linewidth=0.5)
    ax_sL.set_xlim(-xlim_s*1.15, 0.25)
    ax_sL.set_ylim(-0.5, len(sender_genes) - 0.5)
    ax_sL.set_yticks([])
    ax_sL.invert_xaxis()
    ax_sL.tick_params(axis='x', labelsize=6.5)
    for s_ in ("top", "right"): ax_sL.spines[s_].set_visible(False)
    ax_sL.spines['left'].set_visible(False)
    ax_sL.spines['bottom'].set_linewidth(0.4)
    ax_sL.set_xlabel("Δ log-norm  (FAP-CAF − other fibro)", fontsize=6.8)
    ax_sL.set_title("Sender", fontsize=8.5, fontweight='bold', loc='right',
                    color="#B42318", pad=4)

    # Sender gene label column
    ax_sLab = fig.add_subplot(gs[1])
    sender_groups = ["maintenance"]*len(SENDER_MAINT) + ["rewiring"]*len(SENDER_REWIR)
    draw_side_labels(ax_sLab, sender_genes, sender_groups,
                     {g: sender_pathway(g) for g in sender_genes})
    # little side bracket labels
    ax_sLab.text(0.03, len(sender_genes) - 0.5 - len(SENDER_MAINT)/2 + 0.5,
                 "maintenance", rotation=90, fontsize=6.3,
                 ha='center', va='center', color="#9E3B2F", fontweight='bold')
    ax_sLab.text(0.03, len(SENDER_REWIR)/2 - 0.5,
                 "rewiring", rotation=90, fontsize=6.3,
                 ha='center', va='center', color="#5E447F", fontweight='bold')

    # Middle connection area
    ax_mid = fig.add_subplot(gs[2])
    ax_mid.set_xlim(0, 1); ax_mid.set_ylim(-0.5, max(len(sender_genes),
                                                       len(recv_genes)) - 0.5)
    ax_mid.axis('off')
    # map gene → y in its sidebar
    y_send = {g: (len(sender_genes) - 1) - i for i, g in enumerate(sender_genes)}
    y_recv = {g: (len(recv_genes)   - 1) - i for i, g in enumerate(recv_genes)}
    # draw curated first (dashed, underneath) then LIANA on top (solid)
    def _edge_sort_key(e):
        return 0 if e[4] == "curated" else 1
    for (lig, rec, pw, strength, src) in sorted(all_edges, key=_edge_sort_key):
        if lig not in y_send or rec not in y_recv:
            continue
        c = PATHWAY_COLORS.get(pw, "#4A5568")
        ls = "solid" if src == "LIANA" else "dashed"
        draw_connection(ax_mid, y_send[lig], y_recv[rec], c, strength,
                        x_left=0.01, x_right=0.99, linestyle=ls)

    # Receiver gene label column
    ax_rLab = fig.add_subplot(gs[3])
    recv_groups = ["maintenance"]*len(RECV_MAINT) + ["rewiring"]*len(RECV_REWIR)
    draw_side_labels(ax_rLab, recv_genes, recv_groups,
                     {g: recv_pathway(g) for g in recv_genes})
    ax_rLab.text(0.97, len(recv_genes) - 0.5 - len(RECV_MAINT)/2 + 0.5,
                 "maintenance", rotation=90, fontsize=6.3,
                 ha='center', va='center', color="#9E3B2F", fontweight='bold')
    ax_rLab.text(0.97, len(RECV_REWIR)/2 - 0.5,
                 "rewiring", rotation=90, fontsize=6.3,
                 ha='center', va='center', color="#5E447F", fontweight='bold')

    # Right receiver Δ bar
    ax_rR = fig.add_subplot(gs[4])
    recv_df = recv.set_index("gene").reindex(recv_genes).reset_index()
    yR = np.arange(len(recv_genes))[::-1]
    delta_r = recv_df["delta_log_norm"].values
    pct_r   = recv_df["pct_grp1"].values
    padj_r  = recv_df["p_adj_BH"].values
    xlim_r = max(abs(delta_r)) * 1.4
    for i, (d, p, pc) in enumerate(zip(delta_r, padj_r, pct_r)):
        yi = yR[i]
        col = UP_COLOR if d >= 0 else DOWN_COLOR
        ax_rR.barh(yi, d, color=col, edgecolor="black",
                    linewidth=0.35, height=0.55, zorder=2)
        stars = sig_stars(p)
        if stars:
            tx = d + (0.05 if d >= 0 else -0.05)
            ha = "left" if d >= 0 else "right"
            ax_rR.text(tx, yi, stars, fontsize=7, fontweight="bold",
                        va='center', ha=ha, color="#333", zorder=4)
    for i, pc in enumerate(pct_r):
        ax_rR.scatter(xlim_r*1.03, yR[i], s=35 + 230 * pc,
                       facecolor="#E2A736", edgecolor="#333",
                       linewidth=0.35, clip_on=False, zorder=3)
    ax_rR.axvline(0, color="black", linewidth=0.5)
    ax_rR.set_xlim(-0.25, xlim_r*1.15)
    ax_rR.set_ylim(-0.5, len(recv_genes) - 0.5)
    ax_rR.set_yticks([])
    ax_rR.tick_params(axis='x', labelsize=6.5)
    for s_ in ("top", "right"): ax_rR.spines[s_].set_visible(False)
    ax_rR.spines['right'].set_visible(False)
    ax_rR.spines['bottom'].set_linewidth(0.4)
    ax_rR.set_xlabel("Δ log-norm  (terminal − lineage)", fontsize=6.8)
    ax_rR.set_title("Receiver", fontsize=8.5, fontweight='bold', loc='left',
                    color="#C0392B", pad=4)

    # Title
    fig.suptitle("Figure 7A  ·  Paired sender → receiver nomination  ·  "
                 "FAP-high CAF ↔ terminal-dedifferentiated epithelium",
                 fontsize=10.5, fontweight='bold', y=0.96)

    # Legends (bottom)
    from matplotlib.lines import Line2D
    pw_legend = [Patch(facecolor=c, edgecolor='none', label=pw)
                 for pw, c in PATHWAY_COLORS.items() if pw != "Other"]
    sign_legend = [Patch(facecolor=UP_COLOR, edgecolor='black', lw=0.3, label='Δ > 0'),
                   Patch(facecolor=DOWN_COLOR, edgecolor='black', lw=0.3, label='Δ < 0')]
    sz_legend = [Line2D([0],[0], marker='o', linestyle='',
                        markerfacecolor="#E2A736", markeredgecolor="#333",
                        markersize=np.sqrt(35+230*pc),
                        label=f"{int(pc*100)}% expr")
                 for pc in (0.10, 0.50, 0.90)]
    fig.legend(handles=pw_legend + sign_legend,
               loc='lower center', bbox_to_anchor=(0.32, 0.01),
               ncol=len(pw_legend)+len(sign_legend),
               fontsize=6.5, frameon=False)
    fig.legend(handles=sz_legend,
               loc='lower center', bbox_to_anchor=(0.78, 0.01),
               ncol=3, fontsize=6.5, frameon=False,
               title="% expressing in on-side group",
               title_fontsize=6.5)
    # connection-strength legend
    fig.text(0.50, 0.09,
             "line thickness ∝ −log10(LIANA aggregate rank)",
             ha='center', fontsize=6.3, style='italic', color='#555')
    fig.text(0.50, 0.065,
             "dashed = curated pair not in LIANA top-12   |   "
             "significance: *** P_adj < 1e-10,  ** < 1e-3,  * < 0.05",
             ha='center', fontsize=6.0, color='#666')

    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / "Figure_7_A_paired.png"
    pdf = OUT / "Figure_7_A_paired.pdf"
    fig.savefig(png, dpi=220, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf, bbox_inches='tight', facecolor='white')
    print("Saved:", png)
    print("Saved:", pdf)

    # print the edges for reference
    print("\n=== Connection list ===")
    for lig, rec, pw, strength, src in all_edges:
        print(f"  {lig:8s} → {rec:8s}  [{pw:14s}]  strength={strength:.2f}  ({src})")


if __name__ == "__main__":
    main()
