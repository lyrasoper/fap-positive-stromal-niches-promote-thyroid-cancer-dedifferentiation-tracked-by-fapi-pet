"""paper_A style module (Python / matplotlib).

Drop-in style baseline distilled from
  <EXTERNAL_DATA>/xenium/sota_2026_pipeline/scripts/84_paper_A_figures.py
+ extras (panel_tag, style_axes, save).

Usage:
    from paper_A_style import (apply_pa_rc, panel_tag, style_axes, save,
                               pa, MM, COL1, COL2)
    apply_pa_rc()
    fig = plt.figure(figsize=(COL2, COL2 * 0.6))
    ...
    save(fig, "Supplementary_Fig_X_…", OUTDIR)
"""
from __future__ import annotations
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns  # noqa: F401  (registers rocket / mako colormaps)
from matplotlib.colors import LinearSegmentedColormap

# ── canvas constants ────────────────────────────────────────────────────────
MM = 1.0 / 25.4
COL1 = 88 * MM     # single-column (Nature/NC)
COL2 = 180 * MM    # double-column


# ── palette ────────────────────────────────────────────────────────────────
# Mirrors paper_A_style.R; keys aligned with reference 84_paper_A_figures.py
# where possible (MODE_COLOR / INTEG_COLOR / CLASS_COLOR).
pa = {
    # 4-tone mode axis (PAX8-loss → SERPINE1-gain)
    "mode": {
        "PAX8-loss":     "#2C6FBB",
        "Weak":          "#9E9E9E",
        "Dual-axis":     "#6A3D9A",
        "SERPINE1-gain": "#C0392B",
    },
    # gene-class gradients
    "identity":     ["#08306B", "#08519C", "#2171B5", "#4292C6", "#6BAED6"],
    "intermediate": ["#3F007D", "#6A51A3", "#9970AB"],
    "invasion":     ["#67000D", "#A50F15", "#CB181D", "#EF3B2C", "#FB6A4A"],
    # short 3-tone class palette used for trend lines / annotations
    "class": {
        "identity":     "#2C6FBB",
        "intermediate": "#6A3D9A",
        "invasion":     "#C0392B",
    },
    # differentiation tier (WD → EMT)
    "diff_state": {
        "WD":   "#1B7F4B",
        "PTC":  "#86C06B",
        "PD":   "#F0A30A",
        "Stem": "#6A3D9A",
        "EMT":  "#C0392B",
    },
    # CAF subtypes
    "caf": {
        "myCAF":      "#7B3FA0",
        "iCAF":       "#F0A30A",
        "vCAF":       "#1F5C8B",
        "pCAF":       "#2E8B43",
        "CAF_unspec": "#BDBDBD",
    },
    # Okabe-Ito (color-blind safe) — for integrin / contrast pairs
    "integrin": {"avb3": "#009E73", "a5": "#D55E00"},
    # human tissue cohort
    "tissue_human": {
        "NT":  "#6BAED6",
        "PTC": "#F0A30A",
        "ATC": "#C0392B",
    },
    # mouse tissue cohort
    "tissue_mouse": {
        "mNT":   "#6BAED6",
        "mPTC":  "#F0A30A",
        "mATC":  "#C0392B",
        "mmATC": "#67000D",
    },
    # FAP-CAF rebadge
    "fap_caf": {
        "FAP+ infCAF":   "#F0A30A",
        "ecmCAF":        "#6A3D9A",
        "EndMT CAF":     "#009E73",
        "RGS15+ myoCAF": "#2C6FBB",
        "adiCAF":        "#86C06B",
    },
    # scatter / trend / diverging anchors
    "scatter_pt":  "#2C6FBB",
    "trend_line":  "#C0392B",
    "div_low":     "#2C6FBB",
    "div_mid":     "#FFFFFF",
    "div_high":    "#C0392B",
    # neutrals
    "axis_gray":    "#1A1A1A",
    "grid_gray":    "#EAEAEA",
    "annot_gray":   "#8A8A8A",
    "panel_border": "#9A9A9A",
}


# ── matplotlib rcParams (the 80% of "looking publishable") ─────────────────
def apply_pa_rc() -> None:
    """Set matplotlib rcParams to paper_A defaults: 7-pt Arial/Helvetica,
    thin axes, vector PDF (fonttype 42), 400 dpi raster.
    """
    mpl.rcParams.update({
        "font.family":          "sans-serif",
        "font.sans-serif":      ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
        "font.size":            7,
        "axes.labelsize":       7,
        "axes.titlesize":       7,
        "xtick.labelsize":      6,
        "ytick.labelsize":      6,
        "legend.fontsize":      5.6,
        "axes.linewidth":       0.6,
        "lines.linewidth":      1.1,
        "xtick.major.width":    0.6,
        "ytick.major.width":    0.6,
        "xtick.major.size":     2.4,
        "ytick.major.size":     2.4,
        "pdf.fonttype":         42,   # Illustrator-editable text
        "ps.fonttype":          42,
        "savefig.dpi":          400,
        "figure.dpi":           200,
        "savefig.bbox":         "tight",
        "savefig.facecolor":    "white",
        "axes.edgecolor":       "#1A1A1A",
        "axes.spines.top":      False,
        "axes.spines.right":    False,
        "legend.frameon":       False,
    })


# ── helpers (panel_tag / style_axes / save) ────────────────────────────────
def panel_tag(ax, letter: str, subtitle: str = "",
              dx: float = -26, dy: float = 7) -> None:
    """Bold lowercase letter at top-left of panel + optional left-aligned subtitle."""
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction",
                xytext=(dx, dy), textcoords="offset points",
                fontsize=9, fontweight="bold", va="bottom", ha="left")
    if subtitle:
        ax.set_title(subtitle, loc="left", fontsize=7, pad=5)


def style_axes(ax) -> None:
    """Remove top + right spines (Tufte-style minimalism)."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save(fig, name: str, outdir: Path | str = ".") -> None:
    """Save once as PDF (vector) + PNG (preview)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{name}.pdf")
    fig.savefig(outdir / f"{name}.png")
    plt.close(fig)


# ── diverging colormap (RdBu-flavoured, anchored at div_low/mid/high) ──────
def pa_cmap_div(name: str = "pa_div") -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list(
        name, [pa["div_low"], pa["div_mid"], pa["div_high"]]
    )


pa_diverging = pa_cmap_div  # alias

__all__ = [
    "apply_pa_rc", "panel_tag", "style_axes", "save",
    "pa", "pa_cmap_div", "pa_diverging",
    "MM", "COL1", "COL2",
]
