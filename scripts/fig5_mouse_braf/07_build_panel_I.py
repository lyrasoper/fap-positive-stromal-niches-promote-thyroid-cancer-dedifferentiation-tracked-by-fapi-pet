# Cancer Research submission - figure code release
# Builds: Figure 5I panel. Reads the tables written by 05_recompute_fig5I_ssgsea.R.

"""重绘 Figure 5 面板 I（2026-08-24）。

背景
----
原面板 Y 轴标注 "ssGSEA score"，但归档中唯一可复现的打分脚本
（111_bulk_mTC/scripts/plot_rank_signature_per_dataset.py）用的是秩百分位均值，
其数值与图面不符；EMT / TGFβ 两个模块的逐样本分数从未归档。
经比对，原面板确为 GSVA-ssGSEA：以 TDS16 重算的 Differentiation_TDS 与图面
逐点吻合（r = 1.00000，RMSE 0.005，误差量级等于像素读数偏差）；
FAP_CAF 用 14 基因 FAP/ECM 表也高度吻合（r = 0.992）；EMT / TGFβ 无任何
候选基因集能复现，故改用可引用的 MSigDB Hallmark 定义整体重算。

绘图风格严格沿用原面板实测参数（见下方常量），字号仍为 Fig 5 的 7/8/9 + 12 阶梯。
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

import os
IO = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>")) / "outputs" / "fig5I"
SCORES = IO / "fig5I_ssGSEA_scores_recomputed.csv"   # written by 05_recompute_fig5I_ssgsea.R
STATS = IO / "fig5I_pairwise_stats.csv"              # written by 05_recompute_fig5I_ssgsea.R
OUT_PDF = IO / "Fig5_panelI_rebuilt.pdf"
OUT_PNG = IO / "Fig5_panelI_rebuilt.png"

# ── 原面板实测几何（183 mm 成品 PDF 的页面坐标，pt）────────────────────────
PAGE_W = 518.7401733398438
REG_Y0, REG_Y1 = 486.0, 652.0          # 面板 I 所在的空白带之间
REG_H = REG_Y1 - REG_Y0

AX_X0 = [28.2, 153.2, 278.3, 403.3]     # 四个子图绘图区左边界
AX_W = 102.7
AX_TOP, AX_BOT = 517.8, 618.7           # 绘图区上下边界（页面 y，向下为正）
STRIP_TOP, STRIP_BOT = 503.3, 517.8     # facet 标题框

LETTER_XY = (10.0, 493.9)               # 面板字母左上角
TITLE_Y = 489.3                         # 面板标题顶端
YLAB_X = 6.06                           # Y 轴标题（旋转 90°）中心 x
# 3.8 会让字形框伸到页面左缘外（bbox x0 = −1.378 pt，导出 TIFF 最左有墨像素列 = 0，
# 裁切边零余量）。原面板的标题 bbox x0 = +0.88 pt、留白 2.28 pt，故右移 2.26 pt。
XTICK_Y = 620.6                         # x 刻度标签顶端

LW_AXIS = 1.053                         # 轴线 / 箱线 / 须
LW_MEDIAN = 2.107                       # 中位线
LW_STRIP = 2.107                        # facet 框
LW_DOT = 0.700                          # 散点描边
LW_BRACKET = 0.579                      # 显著性括号
TICK_LEN = 2.7
BOX_W = 11.0                            # 箱宽
DOT_D = 4.3                             # 散点直径

FS_TICK = 8.0
FS_STRIP = 8.0
FS_YLAB = 8.0
FS_TITLE = 9.0
FS_STAR = 9.0
FS_NS = 7.0
FS_LETTER = 12.0

GROUPS = ["mNT", "mPTC", "mATC", "remATC"]
COLORS = {"mNT": "#0072B5", "mPTC": "#E18727",
          "mATC": "#BC3C29", "remATC": "#20854E"}
MODULES = ["Differentiation_TDS", "EMT_program", "TGFB_program", "FAP_CAF_program"]

# 显著性括号：(组1, 组2, 层级)。层级 0 最低。
# 1、4 组的比较互不重叠，可同层，比原图省一层空间。
BRACKETS = [("mNT", "mPTC", 0), ("mATC", "remATC", 0),
            ("mPTC", "mATC", 1), ("mPTC", "remATC", 2)]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.bbox": "standard", "savefig.pad_inches": 0.0,
    "axes.linewidth": LW_AXIS,
})


def sym(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def page_to_fig(y: float) -> float:
    """页面 y（向下）→ figure 分数坐标（向上）。"""
    return (REG_H - (y - REG_Y0)) / REG_H


def group_x(k: int, i: int) -> float:
    """第 k 个子图、第 i 组（0 基）的页面 x。ggplot 离散轴默认扩展 0.6。"""
    return AX_X0[k] + AX_W * ((i + 1) - 0.4) / 4.2


def main() -> None:
    df = pd.read_csv(SCORES)
    df = df[df["gene_set"] == "primary"]
    st = pd.read_csv(STATS)
    st = st[st["gene_set"] == "primary"]

    fig = plt.figure(figsize=(PAGE_W / 72.0, REG_H / 72.0), dpi=600)

    rng = np.random.default_rng(20260824)          # 抖动可复现
    log: dict = {"modules": {}}

    for k, mod in enumerate(MODULES):
        d = df[df["module"] == mod]
        vals = {g: d.loc[d["group"] == g, "score"].to_numpy() for g in GROUPS}
        lo = min(v.min() for v in vals.values())
        hi = max(v.max() for v in vals.values())
        span = hi - lo
        gap = 0.070 * span
        lvl = [hi + gap * (1.0 + 1.10 * L) for L in range(3)]
        ylim = (lo - 0.09 * span, lvl[2] + 2.30 * gap)

        ax = fig.add_axes([AX_X0[k] / PAGE_W, page_to_fig(AX_BOT),
                           AX_W / PAGE_W, (AX_BOT - AX_TOP) / REG_H])
        ax.set_xlim(0.4, 4.6)
        ax.set_ylim(*ylim)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_linewidth(LW_AXIS)
        ax.tick_params(axis="y", length=TICK_LEN, width=LW_AXIS,
                       labelsize=FS_TICK, pad=2.0)
        ax.tick_params(axis="x", length=0, pad=1.5)
        ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))

        pt2data = (ylim[1] - ylim[0]) / (AX_BOT - AX_TOP)   # 1 pt = ? 数据单位
        pt2x = (4.6 - 0.4) / AX_W

        for i, g in enumerate(GROUPS):
            v = vals[g]
            x = (i + 1)
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            iqr = q3 - q1
            wlo = v[v >= q1 - 1.5 * iqr].min()
            whi = v[v <= q3 + 1.5 * iqr].max()
            hw = BOX_W * pt2x / 2.0
            ax.plot([x, x], [wlo, q1], color="k", lw=LW_AXIS, zorder=2,
                    solid_capstyle="butt")
            ax.plot([x, x], [q3, whi], color="k", lw=LW_AXIS, zorder=2,
                    solid_capstyle="butt")
            ax.add_patch(Rectangle((x - hw, q1), 2 * hw, q3 - q1,
                                   facecolor=COLORS[g], edgecolor="k",
                                   lw=LW_AXIS, zorder=3))
            ax.plot([x - hw, x + hw], [med, med], color="k", lw=LW_MEDIAN,
                    zorder=4, solid_capstyle="butt")
            jit = rng.uniform(-0.18, 0.18, size=v.size)
            ax.scatter(x + jit, v, s=(DOT_D) ** 2, facecolor=COLORS[g],
                       edgecolor="k", linewidths=LW_DOT, zorder=5, clip_on=False)

        # ── 显著性括号 ────────────────────────────────────────────────────
        bar = []
        for g1, g2, L in BRACKETS:
            row = st[(st["module"] == mod)
                     & (st["comparison"] == f"{g1} vs {g2}")].iloc[0]
            s = sym(row["p_holm"])
            x1, x2 = GROUPS.index(g1) + 1, GROUPS.index(g2) + 1
            y = lvl[L]
            tick = 2.1 * pt2data
            ax.plot([x1, x1, x2, x2], [y - tick, y, y, y - tick],
                    color="k", lw=LW_BRACKET, zorder=6,
                    solid_capstyle="butt", clip_on=False)
            ax.text((x1 + x2) / 2, y + 0.4 * pt2data, s, ha="center", va="bottom",
                    fontsize=FS_STAR if s != "ns" else FS_NS, color="k",
                    zorder=6, clip_on=False)
            bar.append({"comparison": f"{g1} vs {g2}", "p_raw": float(row["p_raw"]),
                        "p_holm": float(row["p_holm"]), "symbol": s})
        log["modules"][mod] = {
            "ylim": [float(ylim[0]), float(ylim[1])],
            "medians": {g: float(np.median(vals[g])) for g in GROUPS},
            "brackets": bar,
        }

        # x 刻度标签（45°，Arial Bold 8 pt）
        ax.set_xticks([])
        for i, g in enumerate(GROUPS):
            fig.text(group_x(k, i) / PAGE_W, page_to_fig(XTICK_Y), g,
                     ha="right", va="top", rotation=45, rotation_mode="anchor",
                     fontsize=FS_TICK, fontweight="bold", family="Arial")

        # facet 标题框
        fig.patches.append(Rectangle(
            (AX_X0[k] / PAGE_W, page_to_fig(STRIP_BOT)),
            AX_W / PAGE_W, (STRIP_BOT - STRIP_TOP) / REG_H,
            transform=fig.transFigure, facecolor="white", edgecolor="k",
            lw=LW_STRIP, zorder=10, clip_on=False))
        fig.text((AX_X0[k] + AX_W / 2) / PAGE_W,
                 page_to_fig((STRIP_TOP + STRIP_BOT) / 2 + 2.8), mod,
                 ha="center", va="baseline", fontsize=FS_STRIP,
                 fontweight="bold", family="Arial", zorder=11)

    # 面板字母、标题、Y 轴标题
    fig.text(LETTER_XY[0] / PAGE_W, page_to_fig(LETTER_XY[1]), "I",
             ha="left", va="top", fontsize=FS_LETTER, fontweight="bold",
             family="Arial")
    fig.text((AX_X0[0] + (AX_X0[3] + AX_W)) / 2 / PAGE_W, page_to_fig(TITLE_Y),
             "Module scores (ssGSEA, GSE118022)", ha="center", va="top",
             fontsize=FS_TITLE, fontweight="bold", family="Arial")
    fig.text(YLAB_X / PAGE_W, page_to_fig((AX_TOP + AX_BOT) / 2),
             "ssGSEA score", ha="center", va="center", rotation=90,
             fontsize=FS_YLAB, family="Arial")

    fig.savefig(OUT_PDF, format="pdf")
    fig.savefig(OUT_PNG, format="png", dpi=600)
    (IO / "panelI_build_log.json").write_text(
        json.dumps(log, ensure_ascii=False, indent=1), encoding="utf8")
    print(f"输出: {OUT_PDF}")
    print(f"      {OUT_PNG}")
    for m, v in log["modules"].items():
        marks = " ".join(f"{b['comparison']}={b['symbol']}" for b in v["brackets"])
        print(f"  {m:20s} ylim=[{v['ylim'][0]:+.3f},{v['ylim'][1]:+.3f}]  {marks}")


if __name__ == "__main__":
    main()
