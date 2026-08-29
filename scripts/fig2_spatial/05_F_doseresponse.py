# Cancer Research submission - figure code release
# Builds: Figure 2 (FAP-TDS dose-response panel).

"""Fig 2 Panel F — 'honest dose-response' v3 preview.

Problem: FAP is 62% zero, so plain deciles put 6 of 10 bins inside the
FAP=0 mass -> the left half of the curve is not a dose gradient and its
wiggle is an artifact of ordering zeros.

Fix: decile ONLY the FAP+ spots (real gradient), and show the FAP-
spots as a single reference baseline. Now all 10 bins span rising FAP.

Renders into outputs/fig2_panels_paperA/_v2_preview/.
Run: python3 scripts/fig2_spatial/05_F_doseresponse.py
"""
from __future__ import annotations
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent.parent / "_shared"))
from paper_A_style import apply_pa_rc, style_axes, save, MM  # noqa: E402

ROOT = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>"))
SRC = ROOT / "outputs/figure2_rebuilt/source_data"
OUT = ROOT / "outputs/fig2_panels_paperA/_v2_preview"
SLICE_PAL = ["#3262A0", "#7E5CAB", "#E8A341", "#1F6E3A",
             "#C03840", "#5786C4", "#7CBF60", "#9E9E9E"]
apply_pa_rc()


def decile(sub: pd.DataFrame, k: int = 10) -> pd.DataFrame:
    pos = sub[sub["FAP_logCP10K"] > 0].copy()
    if len(pos) < k:
        return pd.DataFrame(columns=["bin", "x_med", "y_med"])
    pos["bin"] = pd.qcut(pos["FAP_logCP10K"].rank(method="first"),
                         k, labels=range(1, k + 1)).astype(int)
    g = pos.groupby("bin").agg(x_med=("FAP_logCP10K", "median"),
                               y_med=("TDS16_z", "median")).reset_index()
    return g


def main() -> None:
    df = pd.read_csv(SRC / "Fig2_spatial_spot_level.csv")
    slices = sorted(df["slice"].unique())

    per = {s: decile(df[df["slice"] == s]) for s in slices}
    agg = decile(df)                                   # pooled FAP+ deciles
    base_all = df.loc[df["FAP_logCP10K"] == 0, "TDS16_z"].median()
    base_per = [df.loc[(df["slice"] == s) & (df["FAP_logCP10K"] == 0),
                       "TDS16_z"].median() for s in slices]

    # cross-sample IQR per decile index
    rows = [g.assign(slice=s) for s, g in per.items() if len(g)]
    allper = pd.concat(rows) if rows else pd.DataFrame(columns=["bin", "y_med"])
    band = (allper.groupby("bin")["y_med"]
            .agg(q25=lambda s: np.percentile(s, 25),
                 q75=lambda s: np.percentile(s, 75)).reset_index())

    fig, ax = plt.subplots(figsize=(88 * MM, 60 * MM))
    ax.axhline(0, color="grey", linestyle="--", linewidth=0.4, alpha=0.6)

    # FAP- baseline reference (the 62% zero mass) as a shaded band + line
    b_lo, b_hi = np.nanpercentile(base_per, 25), np.nanpercentile(base_per, 75)
    ax.axhspan(b_lo, b_hi, color="#9E9E9E", alpha=0.16, linewidth=0, zorder=0)
    ax.axhline(base_all, color="#6E6E6E", linewidth=0.9, linestyle=(0, (4, 2)),
               zorder=2)
    ax.text(10.3, base_all, "FAP$^-$ baseline\n(62% of spots)", fontsize=5.0,
            va="center", ha="left", color="#5A5A5A")

    if len(band):
        ax.fill_between(band["bin"], band["q25"], band["q75"],
                        color="#B9C4D0", alpha=0.45, linewidth=0, zorder=1,
                        label="Cross-sample IQR")
    for i, s in enumerate(slices):
        g = per[s]
        if len(g):
            ax.plot(g["bin"], g["y_med"], color=SLICE_PAL[i % len(SLICE_PAL)],
                    linewidth=0.55, alpha=0.5, marker="o", markersize=1.7,
                    zorder=3, label=s)
    if len(agg):
        ax.plot(agg["bin"], agg["y_med"], color="black", linewidth=1.9,
                marker="o", markersize=3.0, zorder=10, label="All FAP$^+$ spots")

    ax.set_xlabel(r"FAP$^+$ decile  (low $\rightarrow$ high FAP)")
    ax.set_ylabel(r"TDS16 ($z$-score, bin median)")
    ax.set_xticks(range(1, 11))
    ax.set_xlim(0.5, 10.5)
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles, labels, fontsize=5.0, loc="center left",
                    bbox_to_anchor=(1.16, 0.5), handletextpad=0.3,
                    borderpad=0.3, labelspacing=0.3,
                    title="Sample", title_fontsize=5.4)
    leg.get_frame().set_linewidth(0)
    ax.set_title("Dose-response within FAP$^+$ spots", fontsize=7,
                 loc="left", pad=4)
    style_axes(ax)
    save(fig, "Fig2_F_honest_doseresponse_v3", OUT)

    # quick monotonic-trend check on aggregate
    if len(agg):
        from scipy.stats import spearmanr
        r, p = spearmanr(agg["bin"], agg["y_med"])
        print(f"  ✓ F v3 — FAP+ deciles; aggregate Spearman ρ(bin,TDS)={r:.2f} P={p:.1e}")
        print(f"    FAP- baseline TDS median = {base_all:.2f}; "
              f"FAP+ decile1={agg['y_med'].iloc[0]:.2f} -> decile10={agg['y_med'].iloc[-1]:.2f}")


if __name__ == "__main__":
    main()
