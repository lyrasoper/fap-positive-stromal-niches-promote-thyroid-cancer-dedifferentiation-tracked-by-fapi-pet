# Cancer Research submission - figure code release
# Builds: Shared revision helpers. fig3_beautify.py loads this module for its
#          module-level matplotlib rcParams (spines, line widths, font
#          embedding), which affect the rendered Figure 3. Its Figure 1 and
#          Figure 2 helper functions are not called by the Figure 3 pipeline.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

"""Python-only revisions for manuscript Figures 1--3.

Outputs
-------
* Fig. 1: corrects panel H to the Hallmark pathways present in the deposited
  source-data workbook. This is deliberately labelled a *partial* revision:
  panels B--H of the legacy single-cell figure cannot be regenerated without
  a Python-readable, duplicate-free cell-level object.
* Fig. 2: regenerates panels G and H from the 96-core plotting table and the
  statistics extracted from the deposited source-data workbook, then embeds
  source-consistent vector panels into the current full figure. These outputs
  remain provisional until the underlying FAP and TPO ROI measurements have
  been audited.
* Fig. 3: regenerates the full moderate main-text figure from its TSV/CSV
  sources and reports the exact epithelial-valid spot count in panel H.

The selected and exclusive graphical backend is Python/matplotlib. No R
process or graphics device is invoked.
"""

from __future__ import annotations

import json
import sys
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from scipy import stats


ROOT = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>"))
OUT = ROOT / ".codex_review" / "python_figure_revisions" / "figs_1_3"
CURRENT = ROOT / ".codex_review" / "legend_figure_assets"
FIG1_INSPECT = ROOT / ".codex_review" / "fig1_inspect.ndjson"
FIG2_INSPECT = ROOT / ".codex_review" / "fig2_inspect.ndjson"
FIG2_POINTS = Path(os.environ.get("EXTERNAL_DATA", "<EXTERNAL_DATA>")) / "fig2" / "Fig2_96P_FAP_IHCV6.csv"
FIG3_CODE = Path(__file__).resolve().parent


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.55,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "legend.frameon": False,
    }
)


def _records(path: Path) -> list[dict]:
    records: list[dict] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith("{"):
            records.append(json.loads(line))
    return records


def _sheet_record(path: Path, sheet: str, address: str | None = None) -> dict:
    matches = [r for r in _records(path) if r.get("sheet") == sheet]
    if address is not None:
        matches = [r for r in matches if r.get("address") == address]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {sheet}!{address} record, found {len(matches)}")
    return matches[0]


def _current_image(index: int) -> Path:
    matches = sorted(CURRENT.glob(f"image{index}.*"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one current image{index}, found {matches}")
    return matches[0]


def _save_panel_bundle(fig: plt.Figure, stem: Path, dpi: int = 600) -> None:
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )


def _save_full_bundle(fig: plt.Figure, stem: Path, dpi: int = 300) -> None:
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=dpi, facecolor="white")
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=dpi,
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )


def load_fig1h() -> pd.DataFrame:
    rec = _sheet_record(FIG1_INSPECT, "Fig1H_top_hallmark", "A1:D11")
    table = rec["preview"]
    df = pd.DataFrame(table[1:], columns=table[0])
    for col in ["NES", "FDR_qval", "nominal_pval"]:
        df[col] = pd.to_numeric(df[col])
    return df


def draw_fig1h(ax: plt.Axes, data: pd.DataFrame, fontsize: float = 6.2) -> None:
    label_map = {
        "Tnfa Signaling Via Nfkb": "TNFα signaling via NF-κB",
        "Interferon Gamma Response": "Interferon-γ response",
        "Inflammatory Response": "Inflammatory response",
        "Allograft Rejection": "Allograft rejection",
        "Kras Signaling Up": "KRAS signaling up",
        "Interferon Alpha Response": "Interferon-α response",
        "Hypoxia": "Hypoxia",
        "Reactive Oxygen Species Pathway": "Reactive oxygen species",
        "Complement": "Complement",
        "Glycolysis": "Glycolysis",
    }
    y = np.arange(len(data))
    values = data["NES"].to_numpy(float)
    ax.barh(y, values, height=0.70, color="#ED7D19", edgecolor="white", linewidth=0.3)
    for yy, (term, val) in enumerate(zip(data["Pathway"], values)):
        ax.text(val + 0.045, yy, label_map.get(term, term), va="center", ha="left",
                fontsize=fontsize, color="black")
    ax.set_yticks([])
    ax.invert_yaxis()
    # Preserve enough right margin for the longest pathway label in the
    # original-size composite; 2.65 clipped the TNFα/NF-κB label.
    ax.set_xlim(0, 3.0)
    ax.set_xlabel("NES", fontsize=fontsize + 0.6, fontweight="bold", labelpad=2)
    ax.tick_params(axis="x", labelsize=fontsize, length=2.3, pad=1.5)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_linewidth(0.75)
    ax.spines["bottom"].set_linewidth(0.75)
    ax.set_title("GSEA | FAP+ infCAF vs other CAFs", fontsize=fontsize + 1.1,
                 fontweight="bold", loc="left", pad=3)


def build_figure1_partial() -> None:
    data = load_fig1h()

    # Standalone editable panel.
    fig, ax = plt.subplots(figsize=(4.25, 2.65), dpi=300)
    draw_fig1h(ax, data, fontsize=7.1)
    fig.text(0.005, 0.985, "H", ha="left", va="top", fontsize=10, fontweight="bold")
    fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.18)
    _save_panel_bundle(fig, OUT / "Figure_1H_corrected_source_supported", dpi=600)
    plt.close(fig)

    # Full-figure partial revision: preserve every other pixel and replace H only.
    with Image.open(_current_image(1)) as source:
        base = source.convert("RGB")
    width, height = base.size
    # The replacement rectangle starts to the right of the panel-G legend.
    x0, y0, x1, y1 = 790, 1210, 1785, 1885
    clean = base.copy()
    ImageDraw.Draw(clean).rectangle((x0, y0, x1, y1), fill="white")
    full = plt.figure(figsize=(width / 300, height / 300), dpi=300)
    bg = full.add_axes([0, 0, 1, 1])
    bg.imshow(clean)
    bg.axis("off")
    rect = [825 / width, 1 - 1845 / height, (1740 - 825) / width, (1845 - 1265) / height]
    axh = full.add_axes(rect)
    draw_fig1h(axh, data, fontsize=4.4)
    # Retain the original panel label H (outside the replaced rectangle).
    _save_full_bundle(full, OUT / "Figure_1_partial_revision_H_only", dpi=300)
    plt.close(full)


def load_fig2_stats() -> tuple[float, dict[tuple[str, str], str], float, float]:
    kw = _sheet_record(FIG2_INSPECT, "Fig2G_statistics", "A1:G2")["preview"]
    dunn = _sheet_record(FIG2_INSPECT, "Fig2G_statistics", "A4:G11")["preview"]
    corr = _sheet_record(FIG2_INSPECT, "Fig2H_statistics", "A1:E2")["preview"]
    kw_df = pd.DataFrame(kw[1:], columns=kw[0])
    dunn_df = pd.DataFrame(dunn[2:], columns=dunn[1])
    corr_df = pd.DataFrame(corr[1:], columns=corr[0])
    adjacent = {
        (row.group1, row.group2): row.signif for row in dunn_df.itertuples(index=False)
        if (row.group1, row.group2) in {
            ("NT", "PTC"), ("PTC", "DDTC"), ("DDTC", "PDTC/ATC")
        }
    }
    return (
        float(kw_df.loc[0, "p_value"]),
        adjacent,
        float(corr_df.loc[0, "rho"]),
        float(corr_df.loc[0, "p_value"]),
    )


GROUP_ORDER = ["NT", "PTC", "DDTC", "PDTC/ATC"]
GROUP_COLORS = ["#8FB1C9", "#E6AF55", "#DE7A31", "#C65A52"]


def _add_bracket(ax: plt.Axes, x1: float, x2: float, y: float, h: float, text: str,
                 fontsize: float) -> None:
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="#333", lw=0.55, clip_on=False)
    ax.text((x1 + x2) / 2, y + h + 0.006, text, ha="center", va="bottom",
            fontsize=fontsize, fontweight="bold")


def draw_fig2g(ax: plt.Axes, df: pd.DataFrame, kw_p: float,
               adjacent: dict[tuple[str, str], str], fontsize: float = 5.2) -> None:
    values = [df.loc[df["Group"] == g, "FAP/DAPI"].to_numpy(float) for g in GROUP_ORDER]
    box = ax.boxplot(
        values,
        positions=np.arange(4),
        widths=0.58,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#222", "linewidth": 0.8},
        whiskerprops={"color": "#555", "linewidth": 0.55},
        capprops={"color": "#555", "linewidth": 0.55},
        boxprops={"edgecolor": "#555", "linewidth": 0.55},
    )
    for patch, color in zip(box["boxes"], GROUP_COLORS):
        patch.set_facecolor(color)
        patch.set_alpha(0.95)

    rng = np.random.default_rng(20260806)
    for i, vals in enumerate(values):
        x = i + rng.uniform(-0.13, 0.13, size=len(vals))
        ax.scatter(x, vals, s=9, c="#666", alpha=0.72, edgecolors="none", zorder=3)

    _add_bracket(ax, 0, 1, 0.655, 0.010, adjacent[("NT", "PTC")], fontsize)
    _add_bracket(ax, 1, 2, 0.700, 0.010, adjacent[("PTC", "DDTC")], fontsize)
    _add_bracket(ax, 2, 3, 0.745, 0.010, adjacent[("DDTC", "PDTC/ATC")], fontsize)
    counts = [len(x) for x in values]
    ax.set_xticks(range(4))
    ax.set_xticklabels([f"{g}\nn={n}" for g, n in zip(GROUP_ORDER, counts)], fontsize=fontsize)
    ax.set_ylabel("FAP/DAPI area fraction", fontsize=fontsize + 0.2, labelpad=2)
    ax.set_ylim(-0.01, 0.795)
    ax.set_yticks(np.arange(0, 0.8, 0.1))
    ax.tick_params(axis="y", labelsize=fontsize, length=2.2, pad=1.5)
    ax.tick_params(axis="x", length=2.2, pad=2)
    ax.text(0.985, 0.035,
            r"Kruskal–Wallis $P$ = 6.14 × 10$^{-16}$"
            "\nAdjacent Dunn–Holm tests",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=fontsize - 0.2,
            color="#555")
    ax.text(-0.20, 1.03, "G", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=fontsize + 2.5, fontweight="bold")
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)


def draw_fig2h(ax: plt.Axes, df: pd.DataFrame, rho: float, p_value: float,
               fontsize: float = 5.2) -> None:
    x = df["FAP/DAPI"].to_numpy(float)
    y = df["TPO/DAPI"].to_numpy(float)
    ax.scatter(x, y, s=10, c="#606060", alpha=0.88, edgecolors="none")
    slope, intercept = np.polyfit(x, y, 1)
    xx = np.linspace(x.min(), x.max(), 100)
    ax.plot(xx, slope * xx + intercept, color="#222", lw=0.75, ls="--")
    ax.set_xlabel("FAP/DAPI area fraction", fontsize=fontsize + 0.2, labelpad=2)
    ax.set_ylabel("TPO/DAPI (mean intensity)", fontsize=fontsize + 0.2, labelpad=2)
    ax.tick_params(labelsize=fontsize, length=2.2, pad=1.5)
    ax.text(
        0.985,
        0.985,
        "Spearman ρ = −0.531\n" + r"$P$ = 2.69 × 10$^{-8}$ ($n$ = 96)",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=fontsize,
    )
    ax.text(-0.20, 1.03, "H", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=fontsize + 2.5, fontweight="bold")
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)


def build_figure2() -> None:
    df = pd.read_csv(FIG2_POINTS)
    kw_p, adjacent, rho, p_value = load_fig2_stats()
    calc_kw = stats.kruskal(*[df.loc[df.Group == g, "FAP/DAPI"] for g in GROUP_ORDER]).pvalue
    calc_rho, calc_p = stats.spearmanr(df["FAP/DAPI"], df["TPO/DAPI"])
    if not np.isclose(calc_kw, kw_p, rtol=1e-3):
        raise RuntimeError((calc_kw, kw_p))
    if not (np.isclose(calc_rho, rho, atol=5e-4) and np.isclose(calc_p, p_value, rtol=1e-10)):
        raise RuntimeError((calc_rho, calc_p, rho, p_value))

    # Standalone vector panels.
    figg, axg = plt.subplots(figsize=(3.10, 2.55), dpi=300)
    draw_fig2g(axg, df, kw_p, adjacent, fontsize=6.4)
    figg.subplots_adjust(left=0.19, right=0.99, top=0.88, bottom=0.20)
    _save_panel_bundle(
        figg,
        OUT / "Figure_2G_source_consistent_PROVISIONAL_ROI_AUDIT",
        dpi=600,
    )
    plt.close(figg)

    figh, axh = plt.subplots(figsize=(3.10, 2.30), dpi=300)
    draw_fig2h(axh, df, rho, p_value, fontsize=6.4)
    figh.subplots_adjust(left=0.20, right=0.98, top=0.88, bottom=0.20)
    _save_panel_bundle(
        figh,
        OUT / "Figure_2H_source_consistent_PROVISIONAL_ROI_AUDIT",
        dpi=600,
    )
    plt.close(figh)

    # Full raster+vector composite, preserving the current layout and all other panels.
    with Image.open(_current_image(2)) as source:
        base = source.convert("RGB")
    width, height = base.size
    clean = base.copy()
    ImageDraw.Draw(clean).rectangle((925, 1665, width, height), fill="white")
    fullfig = plt.figure(figsize=(width / 300, height / 300), dpi=300)
    bg = fullfig.add_axes([0, 0, 1, 1])
    bg.imshow(clean)
    bg.axis("off")

    # Axes coordinates are tied to original pixels so the overall figure size is unchanged.
    rect_g = [1030 / width, 1 - 2160 / height, (1595 - 1030) / width, (2160 - 1725) / height]
    rect_h = [1030 / width, 1 - 2570 / height, (1565 - 1030) / width, (2570 - 2240) / height]
    axg = fullfig.add_axes(rect_g)
    axh = fullfig.add_axes(rect_h)
    draw_fig2g(axg, df, kw_p, adjacent, fontsize=3.9)
    draw_fig2h(axh, df, rho, p_value, fontsize=3.9)
    _save_full_bundle(
        fullfig,
        OUT / "Figure_2_revised_GH_PROVISIONAL_ROI_AUDIT",
        dpi=300,
    )
    plt.close(fullfig)


def draw_figure3_panel_h_exact(fig: plt.Figure, gs, full_module) -> tuple[int, float, float]:
    d = pd.read_csv(full_module.AUG_CTX, sep="\t")
    required = [
        "fap_niche_proximity_index",
        "collagen_sender_program_activity",
        "itga2_receiver_program_activity",
        "itga5_receiver_program_activity",
        "dediff_shift_index",
    ]
    d = d.dropna(subset=required).copy()
    d = d[d["epi_valid_primary"] == True].copy()

    def scale01(values: pd.Series) -> np.ndarray:
        arr = values.to_numpy(float)
        return (arr - np.nanmin(arr)) / max(np.nanmax(arr) - np.nanmin(arr), 1e-9)

    d["local_ecm_integrin_support_index"] = (
        scale01(d[required[0]])
        + scale01(d[required[1]])
        + scale01(d[required[2]])
        + scale01(d[required[3]])
    ) / 4.0
    x = d["local_ecm_integrin_support_index"].to_numpy(float)
    y = d["dediff_shift_index"].to_numpy(float)
    rho, p_value = stats.spearmanr(x, y)

    ax = fig.add_subplot(gs)
    density_cmap = LinearSegmentedColormap.from_list(
        "density_neutral_exact",
        ["#F4F6FA", "#B0C4D9", "#4F6C8F", "#1A2A44"],
    )
    hb = ax.hexbin(x, y, gridsize=42, cmap=density_cmap, bins="log",
                   mincnt=1, linewidths=0, rasterized=True)
    slope, intercept = np.polyfit(x, y, 1)
    xx = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    ax.plot(xx, slope * xx + intercept, color="#1A1A1A", lw=1.6)
    ax.text(
        0.03,
        0.97,
        f"ρ = {rho:+.2f}\n$n$ = {len(d):,} spots",
        transform=ax.transAxes,
        fontsize=7.1,
        va="top",
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.26", fc="white", ec="#888", lw=0.4),
    )
    ax.axhline(0, color="#888", lw=0.4)
    ax.set_xlabel("ECM-integrin support", fontsize=7.3)
    ax.set_ylabel("Dediff shift", fontsize=7.3)
    full_module.hide_spines(ax, keep=("left", "bottom"))
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    cb = plt.colorbar(hb, ax=ax, fraction=0.048, pad=0.02)
    cb.set_label("log10 spot count", fontsize=6.5)
    cb.ax.tick_params(labelsize=5.6)
    return len(d), float(rho), float(p_value)


def build_figure3() -> tuple[int, float, float]:
    sys.path.insert(0, str(FIG3_CODE))
    import fig3_full_render as full  # type: ignore
    import fig3_moderate_variant as moderate  # type: ignore

    fig = plt.figure(figsize=(7.2, 10.9), dpi=300)
    fig.patch.set_facecolor("white")
    outer = gridspec.GridSpec(
        4,
        1,
        figure=fig,
        height_ratios=[1.45, 1.42, 4.15, 1.55],
        hspace=0.27,
        left=0.060,
        right=0.985,
        top=0.985,
        bottom=0.050,
    )
    full.draw_panel_A(fig, outer[0])
    g_mid = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=outer[1], width_ratios=[1.18, 1.35, 0.92], wspace=0.45
    )
    moderate.draw_panel_B_support(fig, g_mid[0])
    moderate.draw_panel_C_liana(fig, g_mid[1])
    moderate.draw_panel_D_programs(fig, g_mid[2])
    moderate.draw_panel_E_spatial(fig, outer[2])
    g_bottom = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=outer[3], width_ratios=[1.18, 0.88, 1.15], wspace=0.34
    )
    moderate.draw_panel_F_distance(fig, g_bottom[0])
    moderate.draw_panel_G_sample_shift(fig, g_bottom[1])
    n_spots, rho, p_value = draw_figure3_panel_h_exact(fig, g_bottom[2], full)

    moderate.add_panel_label(fig, 0.014, 0.982, "A")
    moderate.add_panel_label(fig, 0.014, 0.815, "B")
    moderate.add_panel_label(fig, 0.365, 0.815, "C")
    moderate.add_panel_label(fig, 0.735, 0.815, "D")
    moderate.add_panel_label(fig, 0.014, 0.665, "E")
    moderate.add_panel_label(fig, 0.014, 0.218, "F")
    moderate.add_panel_label(fig, 0.365, 0.218, "G")
    moderate.add_panel_label(fig, 0.610, 0.218, "H")
    _save_panel_bundle(
        fig,
        OUT / "Figure_3_revised_exact_spot_count_SOURCE_DATA_UPDATE_REQUIRED",
        dpi=600,
    )
    plt.close(fig)
    return n_spots, rho, p_value


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    build_figure1_partial()
    build_figure2()
    n_spots, rho, p_value = build_figure3()
    print(f"Outputs written to {OUT}")
    print(
        "Figure 3H local analysis: "
        f"n={n_spots:,}, rho={rho:+.12f}, scipy_p={p_value!r}; "
        "use the conservative manuscript bound P < 0.0001."
    )


if __name__ == "__main__":
    main()
