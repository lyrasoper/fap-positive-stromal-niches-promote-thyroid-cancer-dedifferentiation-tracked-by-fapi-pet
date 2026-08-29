# Cancer Research submission - figure code release
# Builds: Figure 2 (spatial panels).

"""Rebuild the SPATIAL-TRANSCRIPTOMICS panels of Fig 2 in matplotlib paper_A
style. Outputs vector PDFs.

Spatial-only version for the public bioinformatics repository: the experimental
(wet-lab) quantification panels have been removed; only spatial-transcriptomics
panels are built here.

Panels built:
  Panel B — TDS16_z spatial maps on tissue (8-slice tile grid)
  Panel C — FAP_logCP10K spatial maps on tissue
  Panel D — TDS_z vs FAP_logCP10K across all spatial spots
  Panel E — dose-response per sample (FAP bin -> TDS bin median)

Run: python3 scripts/fig2_spatial/03_panels_paperA.py
"""
from __future__ import annotations
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent.parent / "_shared"))
from paper_A_style import (  # noqa: E402
    apply_pa_rc, style_axes, save, pa, MM,
)

ROOT = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>"))
SRC = ROOT / "outputs/figure2_rebuilt/source_data"
OUT = ROOT / "outputs/fig2_panels_paperA"
OUT.mkdir(parents=True, exist_ok=True)

apply_pa_rc()


# ── shared palette ─────────────────────────────────────────────────────────
GROUP_PAL = {
    "NT":       "#7BA6C9",   # cool blue
    "PTC":      "#E8A341",   # amber
    "DDTC":     "#7E5CAB",   # violet (intermediate / dual-axis)
    "PDTC/ATC": "#C03840",   # coral red
    "PDTC":     "#D04E50",   # if split later
    "ATC":      "#7D1B1F",   # if split later
}
GROUP_ORDER_4 = ["NT", "PTC", "DDTC", "PDTC/ATC"]
SLICE_PAL = ["#3262A0", "#7E5CAB", "#E8A341", "#1F6E3A",
             "#C03840", "#5786C4", "#7CBF60", "#9E9E9E"]


def fmt_p_italic(p: float) -> str:
    if p < 1e-100:
        return r"$P$ < 10$^{-100}$"
    if p < 0.001:
        exp = int(np.floor(np.log10(p)))
        return rf"$P$ < 10$^{{{exp+1}}}$"
    return rf"$P$ = {p:.3f}"


# ── Panel D — TDS_z vs FAP_logCP10K hexbin density + LOESS-like ────────────
def panel_D_all_spots() -> None:
    p = SRC / "Fig2_spatial_spot_level.csv"
    if not p.exists():
        raise FileNotFoundError(
            f"Missing {p} — produced by an upstream step; see README.")
    df = pd.read_csv(p)
    x = df["FAP_logCP10K"].values
    y = df["TDS16_z"].values
    rho, p = spearmanr(x, y)

    fig, ax = plt.subplots(figsize=(85 * MM, 60 * MM))
    # hexbin density (much better than 36k scatter)
    hb = ax.hexbin(x, y, gridsize=60, cmap="rocket_r",
                   mincnt=1, linewidths=0, rasterized=True)
    cb = fig.colorbar(hb, ax=ax, fraction=0.045, pad=0.03)
    cb.ax.tick_params(labelsize=5.6, width=0.4, length=2)
    cb.set_label("n spots", fontsize=6)
    cb.outline.set_linewidth(0.4)

    # smoothed running mean
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    win = 600
    yr = pd.Series(ys).rolling(win, center=True, min_periods=win // 3).mean().values
    ax.plot(xs, yr, color=pa["trend_line"], linewidth=1.4, alpha=0.95)

    ax.set_xlabel(r"FAP (log$_2$CP10K)")
    ax.set_ylabel(r"TDS16 ($z$-score)")
    ax.set_title(rf"All spots   $\rho$ = {rho:.3f},  {fmt_p_italic(p)}",
                 fontsize=7, loc="left", pad=4)
    style_axes(ax)
    save(fig, "Fig2_Panel_D_FAP_vs_TDS_allspots_paperA", OUT)
    print(f"  ✓ Panel D — TDS vs FAP all spots (n={len(df):,}, ρ={rho:.3f})")


# ── Panel E — Dose-response per sample (binned) ────────────────────────────
def panel_E_dose_response() -> None:
    p_bs = SRC / "Fig2_spatial_binned_by_sample.csv"
    p_ba = SRC / "Fig2_spatial_binned_all.csv"
    for p in (p_bs, p_ba):
        if not p.exists():
            raise FileNotFoundError(
                f"Missing {p} — produced by an upstream step; see README.")
    bs = pd.read_csv(p_bs)
    ba = pd.read_csv(p_ba)
    slices = sorted(bs["slice"].unique())

    fig, ax = plt.subplots(figsize=(85 * MM, 60 * MM))
    for i, s in enumerate(slices):
        sub = bs[bs["slice"] == s].sort_values("FAP_bin")
        ax.plot(sub["FAP_bin"].values, sub["y_med"].values,
                color=SLICE_PAL[i % len(SLICE_PAL)],
                linewidth=0.7, alpha=0.85, marker="o", markersize=2.2,
                label=s)
    # overall line (bold)
    ax.plot(ba["FAP_bin"].values, ba["y_med"].values,
            color="black", linewidth=1.8, marker="o", markersize=3.0,
            label="All samples", zorder=10)
    ax.axhline(0, color="grey", linestyle="--", linewidth=0.4, alpha=0.6)
    ax.set_xlabel("FAP decile bin (low → high)")
    ax.set_ylabel(r"TDS16 ($z$-score, bin median)")
    ax.set_xticks(range(1, 11))
    leg = ax.legend(fontsize=5.0, loc="center left",
                    bbox_to_anchor=(1.005, 0.5),
                    handletextpad=0.3, borderpad=0.3, labelspacing=0.3,
                    title="Sample", title_fontsize=5.4)
    leg.get_frame().set_linewidth(0)
    style_axes(ax)
    ax.set_title("Dose-response per sample (10 FAP-bin median)",
                 fontsize=7, loc="left", pad=4)
    save(fig, "Fig2_Panel_E_dose_response_per_sample_paperA", OUT)
    print(f"  ✓ Panel E — Dose-response per sample (n={len(slices)} slices)")


# NOTE: the experimental (wet-lab) panels were removed from this
# spatial-only version.


# ── Panel B / C — spatial heatmap on tissue (8-slice tile grid) ────────────
# Reads outputs/figure2_rebuilt/source_data/Fig2_spots_with_coords.csv
# produced by scripts/fig2_spatial/01_extract_spot_xy.R.

SLICE_ORDER = ["P12", "P17", "P26", "P32", "P44", "P57", "P83", "P98"]


def _tile_grid_panel(df, value_col: str, cmap: str | object,
                     cbar_label: str, out_stem: str,
                     vmin, vmax, zero_threshold: float = 0.0,
                     spot_size_signal: float = 3.5) -> None:
    """8-tile (4 cols × 2 rows) spatial scatter. Two-layer rendering:
        layer 1: all spots in light grey (tissue outline, even where gene
                 not expressed)
        layer 2: spots with value > zero_threshold, colored + larger + on top
    Shared colorbar at right.
    """
    slices = [s for s in SLICE_ORDER if s in df["slice"].unique()]
    ncols, nrows = 4, 2
    fig = plt.figure(figsize=(170 * MM, 90 * MM))
    gs = plt.matplotlib.gridspec.GridSpec(
        nrows, ncols + 1, figure=fig,
        width_ratios=[1, 1, 1, 1, 0.08],
        wspace=0.20, hspace=0.30,
        left=0.02, right=0.96, top=0.96, bottom=0.04,
    )
    last_sc = None
    for k, s in enumerate(slices):
        r, c = divmod(k, ncols)
        ax = fig.add_subplot(gs[r, c])
        sub = df[df["slice"] == s]
        # Layer 1: all spots, light grey (tissue outline)
        ax.scatter(sub["x"], sub["y"], c="#D8D8D8",
                   s=1.6, edgecolor="none", marker=".",
                   rasterized=True, zorder=1)
        # Layer 2: signal-bearing spots, colored, larger, on top
        signal = sub[sub[value_col] > zero_threshold].sort_values(value_col)
        if len(signal) > 0:
            sc = ax.scatter(signal["x"], signal["y"], c=signal[value_col],
                            cmap=cmap, vmin=vmin, vmax=vmax,
                            s=spot_size_signal, edgecolor="none", marker=".",
                            rasterized=True, zorder=2)
            last_sc = sc
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.text(0.02, 0.98, s, transform=ax.transAxes,
                fontsize=7, fontweight="bold", va="top", ha="left",
                color="#1A1A1A",
                bbox=dict(facecolor="white", edgecolor="none",
                          alpha=0.85, pad=1.0))
    # shared colorbar (use last signal mappable; if all empty, sample one)
    if last_sc is None:
        # fallback: dummy scatter to define cmap
        last_sc = plt.cm.ScalarMappable(cmap=cmap,
                                        norm=plt.Normalize(vmin, vmax))
    cax = fig.add_subplot(gs[:, -1])
    cb = fig.colorbar(last_sc, cax=cax)
    cb.ax.tick_params(labelsize=5.6, width=0.4, length=2)
    cb.set_label(cbar_label, fontsize=6)
    cb.outline.set_linewidth(0.4)
    save(fig, out_stem, OUT)


def panel_B_TDS_on_tissue() -> None:
    coord_csv = SRC / "Fig2_spots_with_coords.csv"
    if not coord_csv.exists():
        print(f"  ⚠ Panel B skipped — missing {coord_csv.name} "
              f"(run scripts/fig2_spatial/01_extract_spot_xy.R first)")
        return
    df = pd.read_csv(coord_csv)
    # TDS_z is dense — every spot has a meaningful value (not sparse).
    # Set threshold below the min so all spots are colored (no grey layer).
    _tile_grid_panel(
        df, "TDS16_z",
        cmap="RdBu_r", cbar_label=r"TDS16 ($z$-score)",
        out_stem="Fig2_Panel_B_TDS_on_tissue_paperA",
        vmin=-1.5, vmax=1.5,
        zero_threshold=-1e9,    # all spots
        spot_size_signal=3.0,
    )
    n_spots = len(df)
    n_slices = df["slice"].nunique()
    print(f"  ✓ Panel B — TDS16_z on tissue ({n_slices} slices, {n_spots:,} spots)")


def panel_C_FAP_on_tissue() -> None:
    coord_csv = SRC / "Fig2_spots_with_coords.csv"
    if not coord_csv.exists():
        print(f"  ⚠ Panel C skipped — missing {coord_csv.name}")
        return
    df = pd.read_csv(coord_csv)
    # FAP is sparse (62% zero) — two-layer: grey for zero, rocket for signal.
    vmax_c = float(np.percentile(df["FAP_logCP10K"], 99))
    _tile_grid_panel(
        df, "FAP_logCP10K",
        cmap="rocket_r", cbar_label=r"FAP (log$_2$CP10K)",
        out_stem="Fig2_Panel_C_FAP_on_tissue_paperA",
        vmin=0, vmax=vmax_c,
        zero_threshold=0.0,
        spot_size_signal=4.5,    # signal spots bigger to pop on sparse panels
    )
    n_signal = (df["FAP_logCP10K"] > 0).sum()
    print(f"  ✓ Panel C — FAP_logCP10K on tissue "
          f"({n_signal:,}/{len(df):,} FAP+ spots, vmax={vmax_c:.2f})")


def main() -> None:
    print(f"paper_A Fig 2 spatial panel rebuild → {OUT.relative_to(ROOT)}")
    panel_B_TDS_on_tissue()
    panel_C_FAP_on_tissue()
    panel_D_all_spots()
    panel_E_dose_response()
    print()
    print("Done. Drop into Illustrator:")
    for f in sorted(OUT.glob("*.pdf")):
        print(f"  {f.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
