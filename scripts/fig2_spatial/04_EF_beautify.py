# Cancer Research submission - figure code release
# Builds: Figure 2 (panel rebuild).

"""Fig 2 panels E & F (current lettering) — beautify v2 preview.

Improves on scripts/fig2_spatial/03_panels_paperA.py:
  Panel E "All spots"  (rebuild Panel D):
      - hexbin with LOG color norm  -> the 62%-FAP=0 stripe no longer washes
        out the cloud; real density structure becomes visible
      - trend line cased in white so it reads over both pale & dark hexbins
  Panel F "Dose-response" (rebuild Panel E):
      - per-sample lines recede (alpha/​thin); bold black "All samples" leads
      - cross-sample IQR ribbon (q25-q75 of the 8 sample medians per bin)

Renders into outputs/fig2_panels_paperA/_v2_preview/ (PDF + PNG).
Run: python3 scripts/fig2_spatial/04_EF_beautify.py
"""
from __future__ import annotations
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LogNorm
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent.parent / "_shared"))
from paper_A_style import apply_pa_rc, style_axes, save, MM  # noqa: E402

ROOT = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>"))
SRC = ROOT / "outputs/figure2_rebuilt/source_data"
OUT = ROOT / "outputs/fig2_panels_paperA/_v2_preview"
OUT.mkdir(parents=True, exist_ok=True)

apply_pa_rc()
SLICE_PAL = ["#3262A0", "#7E5CAB", "#E8A341", "#1F6E3A",
             "#C03840", "#5786C4", "#7CBF60", "#9E9E9E"]


def fmt_p_italic(p: float) -> str:
    if p < 1e-100:
        return r"$P$ < 10$^{-100}$"
    if p < 0.001:
        exp = int(np.floor(np.log10(p)))
        return rf"$P$ < 10$^{{{exp+1}}}$"
    return rf"$P$ = {p:.3f}"


def panel_E_allspots_logdensity() -> None:
    df = pd.read_csv(SRC / "Fig2_spatial_spot_level.csv")
    x = df["FAP_logCP10K"].values
    y = df["TDS16_z"].values
    rho, p = spearmanr(x, y)

    fig, ax = plt.subplots(figsize=(85 * MM, 60 * MM))
    hb = ax.hexbin(x, y, gridsize=58, cmap="rocket_r",
                   norm=LogNorm(vmin=1, vmax=None),
                   mincnt=1, linewidths=0, rasterized=True)
    cb = fig.colorbar(hb, ax=ax, fraction=0.045, pad=0.03)
    cb.ax.tick_params(labelsize=5.6, width=0.4, length=2)
    cb.set_label("spots per hex (log)", fontsize=6)
    cb.outline.set_linewidth(0.4)

    # smoothed running mean, cased in white so it reads on pale & dark hexes
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    win = 600
    yr = pd.Series(ys).rolling(win, center=True, min_periods=win // 3).mean().values
    ax.plot(xs, yr, color="#1A4E8A", linewidth=1.5, alpha=0.98, zorder=6,
            path_effects=[pe.withStroke(linewidth=2.8, foreground="white")])

    ax.set_xlabel(r"FAP (logCP10K)")
    ax.set_ylabel(r"TDS16 ($z$-score)")
    ax.set_title(rf"All spots   $\rho$ = {rho:.3f},  {fmt_p_italic(p)}"
                 rf"   ($n$ = {len(df):,})",
                 fontsize=7, loc="left", pad=4)
    style_axes(ax)
    save(fig, "Fig2_E_allspots_logdensity_v2", OUT)
    print(f"  ✓ E v2 — log-density hexbin (n={len(df):,}, ρ={rho:.3f})")


def panel_F_dose_response_v2() -> None:
    bs = pd.read_csv(SRC / "Fig2_spatial_binned_by_sample.csv")
    ba = pd.read_csv(SRC / "Fig2_spatial_binned_all.csv")
    slices = sorted(bs["slice"].unique())

    # cross-sample IQR per FAP_bin (spread of the 8 sample medians)
    band = (bs.groupby("FAP_bin")["y_med"]
              .agg(q25=lambda s: np.percentile(s, 25),
                   q75=lambda s: np.percentile(s, 75)).reset_index())

    fig, ax = plt.subplots(figsize=(85 * MM, 60 * MM))
    ax.axhline(0, color="grey", linestyle="--", linewidth=0.4, alpha=0.6)
    ax.fill_between(band["FAP_bin"], band["q25"], band["q75"],
                    color="#B9C4D0", alpha=0.45, linewidth=0, zorder=1,
                    label="Cross-sample IQR")
    for i, s in enumerate(slices):
        sub = bs[bs["slice"] == s].sort_values("FAP_bin")
        ax.plot(sub["FAP_bin"].values, sub["y_med"].values,
                color=SLICE_PAL[i % len(SLICE_PAL)],
                linewidth=0.55, alpha=0.5, marker="o", markersize=1.7,
                zorder=3, label=s)
    ax.plot(ba["FAP_bin"].values, ba["y_med"].values,
            color="black", linewidth=1.9, marker="o", markersize=3.0,
            zorder=10, label="All samples")
    ax.set_xlabel("FAP decile bin (low → high)")
    ax.set_ylabel(r"TDS16 ($z$-score, bin median)")
    ax.set_xticks(range(1, 11))
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles, labels, fontsize=5.0, loc="center left",
                    bbox_to_anchor=(1.005, 0.5), handletextpad=0.3,
                    borderpad=0.3, labelspacing=0.3,
                    title="Sample", title_fontsize=5.4)
    leg.get_frame().set_linewidth(0)
    ax.set_title("Dose–response across FAP deciles",
                 fontsize=7, loc="left", pad=4)
    style_axes(ax)
    save(fig, "Fig2_F_dose_response_v2", OUT)
    print(f"  ✓ F v2 — dose-response + IQR ribbon ({len(slices)} slices)")


if __name__ == "__main__":
    print(f"Fig 2 E/F beautify v2 -> {OUT.relative_to(ROOT)}")
    panel_E_allspots_logdensity()
    panel_F_dose_response_v2()
    for f in sorted(OUT.glob("*.png")):
        print("   ", f.relative_to(ROOT))
