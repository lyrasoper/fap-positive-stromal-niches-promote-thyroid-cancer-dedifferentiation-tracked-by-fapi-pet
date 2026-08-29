## ============================================================================
## _paperA_supp_theme.R  —  shared house-style foundation for the Fig 4
## supplementary bioinformatics set (Supplementary Fig. 14-20).
##
## Purpose: ONE sourced foundation so every Supp 14-20 vector rebuild matches
## the main Fig 4 panels (build_fig4_bioinformatics_panels_paperA.R) EXACTLY,
## and so the audit's cross-figure standardizations live in one place:
##   - off-white #F5F5F5 diverging midpoint everywhere
##   - Arial(=Helvetica via quartz) 7pt body, 5.5pt floor
##   - panel tags 10pt bold black lowercase (theme plot.tag = base+3)
##   - one standard colorbar key dimension
##   - one canonical n=4 / no-BH-FDR caveat string
##   - vector PDF (quartz, UTF-8 safe) + 600dpi ragg PNG via save_vec()
##
## Tokens are copied verbatim from build_fig4_bioinformatics_panels_paperA.R
## (the canonical house reference) — DO NOT diverge them.
## Source this at the top of each rebuild:  source("scripts/_shared/_paperA_supp_theme.R")
## ============================================================================

## resolve the same package library the house scripts use
.libPaths(c(file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_omics/.Rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(tibble)
  library(patchwork); library(ggrepel); library(ggtext); library(scales); library(grid)
})

## ---- palette (paper_A, verbatim) ------------------------------------------
col_wt   <- "#3B6B9C"; col_fap  <- "#C9493A"; col_grey <- "#B7B7B7"
col_up   <- col_fap;   col_dn   <- col_wt
col_purple <- "#4B1E70"; col_gold <- "#E2B200"
col_mid  <- "#F5F5F5"                       # off-white diverging midpoint (NEVER pure white)
col_na   <- "#EEEEEE"
col_lineage <- "#C9493A"; col_stromal <- "#3B6B9C"
grp_levels <- c("WT host", "FAP-deficient host")
grp_cols <- setNames(c(col_wt, col_fap), grp_levels)
## cell-type / CAF-subtype palette (verbatim from the house heatmaps)
ct_cols <- c("myCAF"="#C9493A","iCAF"="#E08F3A","apCAF"="#5E4A7B",
             "Thyroid_tumor"="#3B6B9C","Immune"="#88B0A7","Endothelial"="#B7B7B7")
## driver-class palette for the human-cohort supps (16/19/20) — same blue->red spine
driver_cols <- c("BRAF V600E"="#C9493A","RAS-driver"="#3B6B9C","Other"="#B7B7B7",
                 "BRAF"="#C9493A","RAS"="#3B6B9C")
## RNA vs protein layer encoding — ONE convention across 14/17/18/20
layer_cols <- c("RNA"="#6E8CB3","Protein"="#264653")
layer_levels <- c("RNA","Protein")

## ---- theme (paper_A, verbatim) --------------------------------------------
base_font_size <- 7
axis_min_size  <- 5.5      # hard legibility floor; nothing prints smaller
key_w <- unit(2.5, "mm")   # standard colorbar key width  (audit: standardize)
key_h <- unit(5.0, "mm")   # standard colorbar key height (audit: standardize)

theme_nc <- function(base = base_font_size){
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(face="plain", size = base + 1, hjust = 0),
      plot.subtitle = element_text(size = base - 1.2, colour = "grey35"),
      plot.tag      = element_text(face = "bold", size = base + 3, hjust = 0),  # 10pt bold black
      axis.title    = element_text(size = base),
      axis.text     = element_text(size = base - 0.5, colour = "black"),
      axis.line     = element_line(linewidth = 0.3),
      axis.ticks    = element_line(linewidth = 0.3),
      strip.background = element_blank(),
      strip.text    = element_text(size = base, face = "plain"),
      legend.key.size = unit(3.2, "mm"),
      legend.text   = element_text(size = base - 0.5),
      legend.title  = element_text(size = base - 0.5, face = "plain"),
      legend.background = element_blank(), legend.box.background = element_blank(),
      panel.grid    = element_blank(),
      plot.margin   = margin(3, 4, 3, 4)
    )
}

## ---- helpers (paper_A, verbatim) ------------------------------------------
tag_group <- function(df) df %>% mutate(group = factor(ifelse(grepl("^WT", sample), "WT host", "FAP-deficient host"), levels = grp_levels))
p_label <- function(p){ if (is.na(p)) return("ns"); if (p < 0.001) sprintf("italic(P) == %.1e", p) else sprintf("italic(P) == %.3f", p) }
safe_z <- function(x){ sx <- sd(x, na.rm=TRUE); mx <- mean(x, na.rm=TRUE); if (!is.finite(sx) || sx==0) return(rep(0,length(x))); (x-mx)/sx }
to_mgi <- function(x){ x <- toupper(x); paste0(substr(x,1,1), tolower(substr(x,2,nchar(x)))) }

## ---- standardized diverging fill (off-white midpoint everywhere) ----------
scale_fill_div <- function(limits = c(-2.2, 2.2), name = "z", ...){
  scale_fill_gradient2(low = col_dn, mid = col_mid, high = col_up, midpoint = 0,
                       limits = limits, oob = scales::squish, na.value = col_na,
                       name = name, ...)
}
## standardized colorbar guide for ggplot heatmaps
guide_cbar <- function() guides(fill = guide_colourbar(barwidth = key_w, barheight = key_h,
                                                       frame.colour = NA, ticks.colour = "grey40"))

## ---- vector export: quartz PDF (UTF-8 safe) + 600dpi ragg PNG --------------
## w,h in mm. Produces <dir>/<name>.pdf (editable vector) and <dir>/<name>.png.
save_vec <- function(p, name, dir, w, h){
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  pdf_p <- file.path(dir, paste0(name, ".pdf"))
  png_p <- file.path(dir, paste0(name, ".png"))
  if (capabilities("cairo")) grDevices::cairo_pdf(pdf_p, width = w/25.4, height = h/25.4) else if (capabilities("aqua")) grDevices::quartz(file = pdf_p, type = "pdf", width = w/25.4, height = h/25.4) else grDevices::pdf(pdf_p, width = w/25.4, height = h/25.4)
  print(p); dev.off()
  ggsave(png_p, p, width = w, height = h, units = "mm", dpi = 600, device = ragg::agg_png, limitsize = FALSE)
  cat("  wrote", name, sprintf("(%gx%g mm)  -> pdf+png\n", w, h))
  invisible(list(pdf = pdf_p, png = png_p))
}

cat("[_paperA_supp_theme.R] loaded — Helvetica 7pt, off-white midpoint, save_vec() ready\n")
