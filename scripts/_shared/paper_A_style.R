# paper_A style module
# Distilled from <EXTERNAL_DATA>/xenium/sota_2026_pipeline/outputs/figures/paper_A/
# Goal: a disciplined ~10-colour palette + minimal sans-serif theme so that
# every Supp-Fig panel in the FAP project shares one visual language.
#
# Usage in an R script:
#   source("scripts/_shared/paper_A_style.R")
#   common_theme <- pa_theme()
#   ggplot(...) + scale_fill_manual(values = pa$mode)
#
# Designed to drop in next to the existing `common_theme`/`save_combo` blocks
# in scripts 66-71 without touching their data logic.

# ─── palette ─────────────────────────────────────────────────────────────────
pa <- list(
  # primary 4-tone "mode" axis  (paper_A Fig 3/4 modes)
  mode = c(
    "PAX8-loss"     = "#3262A0",  # cool blue   — identity loss / NT-like
    "Weak"          = "#8E8E8E",  # neutral gray — unchanged / control
    "Dual-axis"     = "#7E5CAB",  # violet      — intermediate / co-altered
    "SERPINE1-gain" = "#C03840"   # coral red   — invasion / EMT
  ),

  # identity-loss gene gradient (paper_A Fig 2a/2d blues)
  identity = c("#1F4E89", "#2A6BAA", "#3262A0", "#3A75B5", "#5786C4"),

  # intermediate (paper_A Fig 2b purples)
  intermediate = c("#4D2D8C", "#7B4DA8", "#9B7BC1"),

  # invasion / EMT (paper_A Fig 2c reds)
  invasion = c("#7D1B1F", "#A82832", "#C03840", "#D04E50", "#E26568"),

  # differentiation tiers (paper_A Fig 3c stacked-bar)
  diff_state = c(
    "WD"   = "#1F6E3A",  # deep green
    "PTC"  = "#7CBF60",  # light green
    "PD"   = "#E8A341",  # amber
    "Stem" = "#7E5CAB",  # purple
    "EMT"  = "#C03840"   # red
  ),

  # CAF subtypes (paper_A FigS1)
  caf = c(
    "myCAF"        = "#7E5CAB",
    "iCAF"         = "#E8A341",
    "vCAF"         = "#2C5E9D",
    "pCAF"         = "#1F6E3A",
    "CAF_unspec"   = "#CCCCCC"
  ),

  # integrin contrast (paper_A Fig 4a/4d)
  integrin = c("avb3" = "#16A085", "a5" = "#E67E22"),

  # disease-state axis for the FAP human cohort (NT → PTC → ATC).
  # Re-mapped to mirror paper_A's blue→amber→red progression so it carries
  # the same "low → high dedifferentiation" semantics as the mode axis.
  tissue_human = c(
    "NT"  = "#7BA6C9",  # calm cool blue
    "PTC" = "#E8A341",  # amber (mid-severity)
    "ATC" = "#C03840"   # coral red (severe)
  ),

  # mouse cohort tissue (mNT, mPTC, mATC, mmATC) — same progression + dark cap
  tissue_mouse = c(
    "mNT"   = "#7BA6C9",
    "mPTC"  = "#E8A341",
    "mATC"  = "#C03840",
    "mmATC" = "#7D1B1F"
  ),

  # FAP-CAF rebadge (existing FAP/CAF subtype names → paper_A palette slots)
  fap_caf = c(
    "FAP+ infCAF"   = "#E8A341",   # amber  (≈ iCAF)
    "ecmCAF"        = "#7E5CAB",   # violet (≈ myCAF)
    "EndMT CAF"     = "#16A085",   # teal
    "RGS15+ myoCAF" = "#3262A0",   # deep blue
    "adiCAF"        = "#7CBF60"    # light green
  ),

  # scatter / trend defaults
  scatter_pt   = "#5786C4",   # azure
  trend_line   = "#C03840",   # coral

  # diverging heatmap anchors (paper_A heatmaps, RdBu-flavoured)
  div_low      = "#2C5E9D",
  div_mid      = "#F5F5F5",
  div_high     = "#C03840",
  supp_ref     = "#7030A0",

  # auxiliary grays
  axis_gray    = "#1A1A1A",
  grid_gray    = "#EAEAEA",
  annot_gray   = "#8A8A8A",
  panel_border = "#9A9A9A"
)


# ─── theme ───────────────────────────────────────────────────────────────────
pa_theme <- function(base_size = 10, base_family = "Helvetica",
                     border = TRUE) {
  # paper_A panels: thin axes, light horizontal gridlines, no chartjunk.
  t <- ggplot2::theme_classic(base_size = base_size,
                              base_family = base_family) +
    ggplot2::theme(
      text          = ggplot2::element_text(face = "plain",
                                            color = pa$axis_gray),
      axis.text     = ggplot2::element_text(color = pa$axis_gray),
      axis.title    = ggplot2::element_text(face = "plain"),
      axis.line     = ggplot2::element_line(color = pa$axis_gray,
                                            linewidth = 0.6),
      axis.ticks    = ggplot2::element_line(color = pa$axis_gray,
                                            linewidth = 0.6),
      panel.grid.major.y = ggplot2::element_line(color = pa$grid_gray,
                                                 linewidth = 0.3),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.text   = ggplot2::element_text(size = base_size - 2),
      legend.title  = ggplot2::element_text(size = base_size - 1,
                                            face = "plain"),
      legend.key    = ggplot2::element_blank(),
      plot.title    = ggplot2::element_text(face = "bold",
                                            hjust = 0.5,
                                            size = base_size + 1),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "plain")
    )
  if (border) {
    t <- t + ggplot2::theme(
      panel.border = ggplot2::element_rect(color = pa$panel_border,
                                           fill = NA,
                                           linewidth = 0.45),
      axis.line = ggplot2::element_blank()
    )
  }
  t
}


# ─── ggplot scale shortcuts ─────────────────────────────────────────────────
pa_scale_fill_mode    <- function(...)
  ggplot2::scale_fill_manual(values = pa$mode, ...)
pa_scale_color_mode   <- function(...)
  ggplot2::scale_color_manual(values = pa$mode, ...)
pa_scale_fill_diff    <- function(...)
  ggplot2::scale_fill_manual(values = pa$diff_state, ...)
pa_scale_fill_caf     <- function(...)
  ggplot2::scale_fill_manual(values = pa$caf, ...)
pa_scale_fill_tissue  <- function(...)
  ggplot2::scale_fill_manual(values = pa$tissue_human, ...)
pa_scale_color_tissue <- function(...)
  ggplot2::scale_color_manual(values = pa$tissue_human, ...)
pa_scale_color_div    <- function(midpoint = 0, limits = c(-1, 1),
                                  name = "Spearman\nρ", ...) {
  ggplot2::scale_color_gradient2(low = pa$div_low, mid = pa$div_mid,
                                 high = pa$div_high, midpoint = midpoint,
                                 limits = limits, name = name, ...)
}
pa_scale_fill_div     <- function(midpoint = 0, limits = c(-1, 1),
                                  name = "Spearman\nρ", ...) {
  ggplot2::scale_fill_gradient2(low = pa$div_low, mid = pa$div_mid,
                                high = pa$div_high, midpoint = midpoint,
                                limits = limits, name = name, ...)
}

invisible(NULL)
