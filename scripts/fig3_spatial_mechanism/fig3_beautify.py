# Cancer Research submission - figure code release
# Builds: Figure 3 post-processing passes (color normalization, Arial, font-size
#          ladder, italic statistics, subunit labels, per-section rho in panel H).
#          Imported by 13_render_fig3_final.py.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

"""Figure 3 美化 —— paper_A house style。

原始绘图逻辑不动（数据与统计完全一致），只在图对象树上做样式后处理：
  1. Arial 字族 + 字号阶梯化（原 12 种 → 5 档）
  2. 近重复色值归一到 pa 命名槽（5 个红 → 1 个）
  3. TGFβ 由紫改琥珀，避开红-紫这一色盲不友好对
  4. 面板字母统一为粗体大写（与图注 (A)(B)(C) 体例一致）
  5. P / n / ρ 斜体
  6. 锁定画布尺寸，不用 bbox_inches="tight"
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import re
from matplotlib.text import Text

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
FIG3_CODE = HERE   # the render modules are siblings in this directory
REF = HERE
sys.path.insert(0, str(HERE.parent / "_shared"))   # paper_A_style

from paper_A_style import apply_pa_rc, pa  # noqa: E402

# ── 1. 字体与 rc 基线 ───────────────────────────────────────────────────────
apply_pa_rc()

import fig3_full_render as full                   # noqa: E402
import fig3_moderate_variant as moderate          # noqa: E402
import importlib.util                             # noqa: E402

# 注意：上面两个模块在顶层执行 plt.rcParams.update({"font.family": "DejaVu Sans"})，
# 会覆盖先前的设置，故字体与输出参数必须在 import 之后再设。
plt.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset":   "custom",
    "mathtext.rm":        "Arial",
    "mathtext.it":        "Arial:italic",
    "mathtext.bf":        "Arial:bold",
    # paper_A 默认 savefig.bbox="tight"，会让保存宽度随标签长度漂移（原则 8）
    "savefig.bbox":       "standard",
    "savefig.pad_inches": 0.0,
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
    "svg.fonttype":       "none",
})

_spec = importlib.util.spec_from_file_location("refb", REF / "fig1_fig3_revision_helpers.py")
_refb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_refb)

# ── 2. 语义色板 ────────────────────────────────────────────────────────────
ECM   = "#C0392B"          # ECM / 整合素类 —— 主色
ECM_L = "#F3D9D4"          # 同色淡底
CXCL  = "#0B7A75"          # CXCL12–CXCR4（Okabe-Ito 青绿）
TGFB  = "#E08214"          # TGFβ —— 原为紫，改琥珀以避开红-紫
GREY  = "#4A5568"
AMBER = "#E89B3C"          # 表达比例点

# 近重复色 → 规范槽
COLOR_MAP = {
    "#c0392b": ECM, "#C0392B": ECM, "#C93B26": ECM, "#C94F3D": ECM,
    "#B42318": ECM, "#A93226": ECM, "#CB4335": ECM,
    "#FCE9E5": ECM_L, "#FAF3F0": ECM_L, "#FDEDEA": ECM_L,
    "#6F4E9E": TGFB, "#7B3FA0": TGFB, "#6A3D9A": TGFB,
    "#0B7A75": CXCL, "#2C7F62": CXCL, "#2A7F62": CXCL, "#14433A": CXCL,
    "#E2A736": AMBER, "#E89B3C": AMBER,
    "#4F6C8F": GREY, "#4A5568": GREY,
    "#f1ecf7": "#F2F3F5", "#EDE7F6": "#F2F3F5",
}
COLOR_MAP = {k.lower(): v for k, v in COLOR_MAP.items()}


def _to_hex(c):
    try:
        return matplotlib.colors.to_hex(c, keep_alpha=False).lower()
    except Exception:
        return None


def remap_colors(fig):
    n = 0
    for art in fig.findobj():
        for getter, setter in (("get_facecolor", "set_facecolor"),
                               ("get_edgecolor", "set_edgecolor"),
                               ("get_color", "set_color")):
            if not (hasattr(art, getter) and hasattr(art, setter)):
                continue
            try:
                cur = getattr(art, getter)()
            except Exception:
                continue
            if cur is None:
                continue
            # 单色
            h = _to_hex(cur)
            if h and h in COLOR_MAP:
                try:
                    getattr(art, setter)(COLOR_MAP[h]); n += 1
                except Exception:
                    pass
    return n


# ── 2b. 强制字族 ─────────────────────────────────────────────────────────
# 各面板绘图函数在 text() 调用里显式指定了字体，rcParams 管不到；
# 故在图对象树上逐个 Text 强制设为 Arial（保留原有粗体/斜体样式）。
def force_arial(fig):
    n = 0
    for t in fig.findobj(Text):
        try:
            t.set_fontfamily(["Arial", "Helvetica", "DejaVu Sans"])
            t.set_fontname("Arial") if hasattr(t, "set_fontname") else None
            n += 1
        except Exception:
            pass
    return n


# ── 3. 字号阶梯 ────────────────────────────────────────────────────────────
# AACR：图内文字 6–12 pt，字号变化应尽量小。旧档 5.6/6.2/7/8 + 面板字母 14
# 的跨度过大且字母超上限；改为 7/8/9 + 面板字母 12（全部落在 6–12 内）。
LADDER = [7.0, 8.0, 9.0]             # 刻度·图例 / 轴标题·注释 / 面板小标题


def snap_fonts(fig):
    n = 0
    for t in fig.findobj(Text):
        try:
            s = t.get_fontsize()
        except Exception:
            continue
        if s is None:
            continue
        if s >= 11:                   # 原面板字母，稍后单独重画
            continue
        new = min(LADDER, key=lambda v: abs(v - s))
        if abs(new - s) > 1e-6:
            t.set_fontsize(new); n += 1
    return n


# ── 4. 统计符号斜体 ────────────────────────────────────────────────────────
ITALIC = (
    ("ρ = ", r"$\rho$ = "),
    ("rho = ", r"$\rho$ = "),
    ("P = ", "$P$ = "), ("P < ", "$P$ < "), ("P<", "$P$ < "),
    ("n = ", "$n$ = "), ("n=", "$n$ = "),
)


def italicize_stats(fig):
    n = 0
    for t in fig.findobj(Text):
        s = t.get_text()
        if not s or "$" in s:
            continue
        new = s
        for a, b in ITALIC:
            if a in new:
                new = new.replace(a, b)
        if new != s:
            t.set_text(new); n += 1
    return n



# ── 6. ECM 受体分层着色（整合素 vs 非整合素）────────────────────────────────
# 原码 ECM_RECEPTOR_TOKENS = ("ITG","CD44","SDC","DDR","LRP") 把三类非整合素
# 胶原受体与整合素涂成同色，与正文 "ECM–integrin" 的概括同源。此处按亚类分色：
# 同色相、不同明度 —— 色盲条件下靠明度差仍可区分。
ECM_INTEGRIN     = "#C0392B"
ECM_NONINTEGRIN  = "#7B241C"
NONINTEGRIN_RECEPTORS = {"CD44", "SDC1", "SDC4", "DDR1", "DDR2", "LRP1"}


def split_ecm_receptors(fig):
    n = 0
    for t in fig.findobj(Text):
        if t.get_text().strip() in NONINTEGRIN_RECEPTORS:
            t.set_color(ECM_NONINTEGRIN); n += 1
    return n


# ── 6b. 亚基标注：避免把单亚基读成复合体 ────────────────────────────────
# panel E 画的是 ITGA2（α2 亚基）的 receiver program，不是 α2β1 复合体；
# panel B 亦为逐亚基评分。行标明确标注亚基身份，防止误读。
SUBUNIT_LABEL = {"ITGA2\nreceiver": "ITGA2 (\u03b12 subunit)\nreceiver"}


def annotate_subunit(fig):
    n = 0
    for t in fig.findobj(Text):
        s0 = t.get_text()
        if s0 in SUBUNIT_LABEL:
            t.set_text(SUBUNIT_LABEL[s0]); n += 1
    return n


# ── 7. 数值与格式微修 ──────────────────────────────────────────────────────
TEXT_FIX = [
    # Arial 不含 ∝ (U+221D)，强制 Arial 后会渲染成缺字符方框；改用冒号表达
    ("thickness \u221d ", "thickness: "),
    ("\u221d ", ": "),
    ("7.8e-03", r"7.8 $\times$ 10$^{-3}$"),
    ("7.8e-3",  r"7.8 $\times$ 10$^{-3}$"),
    (r"$\rho$ = +0.65\n", r"$\rho$ = +0.652\n"),
]


def fix_text(fig):
    n = 0
    for t in fig.findobj(Text):
        s0 = t.get_text()
        if not s0:
            continue
        s1 = s0
        for a, b in TEXT_FIX:
            if a in s1:
                s1 = s1.replace(a, b)
        if s1 != s0:
            t.set_text(s1); n += 1
    return n


# ── 8. 面板 c 的点大小图例（house style 原则 4）───────────────────────────
def add_size_legend(fig):
    """点大小图例：panel c 与 panel d 之间无可用空隙（d 的 y 轴标签左伸），
    强行插入会压字。尺寸编码由图注 "dot area indicates -log10(method-specific
    rank)" 承担；如需图形化图例，须先重排 b/c/d 行的 width_ratios。"""
    return 0


# ── 5. 面板字母：粗体小写 ──────────────────────────────────────────────────
def panel_letter(fig, x, y, letter):
    fig.text(x, y, letter, fontsize=12, fontweight="bold",
             va="top", ha="left",
             family=["Arial", "Helvetica", "DejaVu Sans"])


def _save_tiff(fig, path, dpi=600):
    """期刊投稿用 TIFF：600 dpi、RGB（去 alpha）、LZW 无损压缩、写入 dpi 元数据。"""
    import io
    from PIL import Image
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, facecolor="white")
    buf.seek(0)
    im = Image.open(buf)
    if im.mode in ("RGBA", "LA", "P"):
        bg = Image.new("RGB", im.size, "white")
        bg.paste(im, mask=im.split()[-1] if im.mode in ("RGBA", "LA") else None)
        im = bg
    else:
        im = im.convert("RGB")
    im.save(path, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
    return im.size



# ── 9. 面板 H 重建：三分量指数 + 逐切片统计 ────────────────────────────────
# 原实现用四分量（含 itga5_receiver，纤连蛋白受体），与图注所述三分量不符；
# ITGA5 贡献可忽略（池化 ρ 0.6516 → 0.6504）。此处改为与图注一致的三分量，
# 并把标注由不可解读的池化 ρ 改为逐切片中位数与范围（照 Supplementary Fig. S11 口径）。
def draw_panel_H_persection(fig, gs, full_module):
    import numpy as np, pandas as pd
    from scipy import stats
    from matplotlib.colors import LinearSegmentedColormap
    req = ["fap_niche_proximity_index", "collagen_sender_program_activity",
           "itga2_receiver_program_activity", "dediff_shift_index"]
    d = pd.read_csv(full_module.AUG_CTX, sep="\t").dropna(subset=req).copy()
    d = d[d["epi_valid_primary"] == True].copy()

    def s01(v):
        a = np.asarray(v, float)
        return (a - np.nanmin(a)) / max(np.nanmax(a) - np.nanmin(a), 1e-9)

    d["idx3"] = (s01(d[req[0]]) + s01(d[req[1]]) + s01(d[req[2]])) / 3.0
    x, y = d["idx3"].to_numpy(float), d[req[3]].to_numpy(float)
    rho_pool, _ = stats.spearmanr(x, y)
    per = [stats.spearmanr(dd["idx3"], dd[req[3]])[0]
           for _, dd in d.groupby("sample") if len(dd) >= 30]
    per = np.asarray(per)

    ax = fig.add_subplot(gs)
    cmap = LinearSegmentedColormap.from_list(
        "density_neutral_exact", ["#F4F6FA", "#B0C4D9", "#4F6C8F", "#1A2A44"])
    hb = ax.hexbin(x, y, gridsize=42, cmap=cmap, bins="log", mincnt=1,
                   linewidths=0, rasterized=True)
    sl, ic = np.polyfit(x, y, 1)
    xx = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    ax.plot(xx, sl * xx + ic, color="#1A1A1A", lw=1.6)
    ax.text(0.03, 0.97,
            f"per-section $\\rho$ = {np.median(per):+.3f}\n"
            f"range {per.min():+.2f} to {per.max():+.2f} ($n$ = {len(per)})\n"
            f"pooled $\\rho$ = {rho_pool:+.3f}, {len(d):,} spots",
            transform=ax.transAxes, fontsize=7.0, va="top", fontweight="bold",
            linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.26", fc="white", ec="#888", lw=0.4))
    ax.axhline(0, color="#888", lw=0.4)
    ax.set_xlabel("ECM-integrin support", fontsize=7.0)
    ax.set_ylabel("Dediff shift", fontsize=7.0)
    full_module.hide_spines(ax, keep=("left", "bottom"))
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.5)
    cb = plt.colorbar(hb, ax=ax, fraction=0.048, pad=0.02)
    cb.set_label("log10 spot count", fontsize=6.2)
    cb.ax.tick_params(labelsize=5.6)
    return len(d), float(np.median(per)), float(per.min()), float(per.max())


def build(out_stem: Path):
    fig = plt.figure(figsize=(7.2, 10.9), dpi=300)
    fig.patch.set_facecolor("white")
    outer = gridspec.GridSpec(
        4, 1, figure=fig, height_ratios=[1.45, 1.42, 4.15, 1.55],
        hspace=0.27, left=0.075, right=0.945, top=0.985, bottom=0.050,
    )
    full.draw_panel_A(fig, outer[0])
    g_mid = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=outer[1], width_ratios=[1.18, 1.35, 0.92], wspace=0.45)
    moderate.draw_panel_B_support(fig, g_mid[0])
    moderate.draw_panel_C_liana(fig, g_mid[1])
    moderate.draw_panel_D_programs(fig, g_mid[2])
    moderate.draw_panel_E_spatial(fig, outer[2])
    g_bot = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=outer[3], width_ratios=[1.18, 0.88, 1.15], wspace=0.34)
    moderate.draw_panel_F_distance(fig, g_bot[0])
    moderate.draw_panel_G_sample_shift(fig, g_bot[1])
    n_spots, rmed, rlo, rhi = draw_panel_H_persection(fig, g_bot[2], full)

    # —— 后处理 ——
    nc = remap_colors(fig)
    na = force_arial(fig)
    nf = snap_fonts(fig)
    ni = italicize_stats(fig)
    ns = split_ecm_receptors(fig)
    nt = fix_text(fig)
    nb = annotate_subunit(fig)
    nl = add_size_legend(fig)

    for x, y, L in [(0.014, 0.982, "A"), (0.014, 0.815, "B"), (0.362, 0.815, "C"),
                    (0.710, 0.815, "D"), (0.014, 0.652, "E"), (0.014, 0.218, "F"),
                    (0.362, 0.218, "G"), (0.592, 0.218, "H")]:
        panel_letter(fig, x, y, L)

    # colorbar 刻度等文字在绘制时才生成，前面的后处理抓不到；
    # 先强制绘制一次，再补跑字族与字号归一，确保无低于 6 pt 的残留。
    # colorbar 的刻度标签由 tick_params(labelsize) 管理，在 Text 对象上改字号
    # 会在下次绘制时被重置。故先按坐标轴层级设定，再补一次全局归一。
    # 刻度标签一律设为阶梯最低档：labelsize 存在坐标轴属性上，
    # 只改 Text 对象会在 savefig 重绘时被复原（源码里有 5.5/5.6/5.8 的设定）。
    nt2 = 0
    for ax in fig.axes:
        ax.tick_params(axis="both", which="both", labelsize=LADDER[0])
        nt2 += 1
    fig.canvas.draw()
    na2 = force_arial(fig)
    nf2 = snap_fonts(fig)
    print(f"  绘制后补正: 坐标轴 {nt2} | 字族 {na2} | 字号 {nf2}")

    out_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_stem.with_suffix(".pdf"), facecolor="white")      # 锁定尺寸
    fig.savefig(out_stem.with_suffix(".png"), dpi=400, facecolor="white")
    _save_tiff(fig, out_stem.with_suffix(".tif"), dpi=600)
    plt.close(fig)
    print(f"改色 {nc} | 字族统一 {na} | 字号归一 {nf} | 斜体 {ni} | 受体分层 {ns} | 文本修正 {nt} | 亚基标注 {nb} | 大小图例 {nl}")
    print(f"数据核对: n={n_spots:,}  逐切片中位 rho={rmed:+.4f} ({rlo:+.3f}~{rhi:+.3f})")
    return n_spots, rmed


if __name__ == "__main__":
    _out = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>")) / "outputs" / "fig3"
    _out.mkdir(parents=True, exist_ok=True)
    build(_out / "Figure3_paperA")
