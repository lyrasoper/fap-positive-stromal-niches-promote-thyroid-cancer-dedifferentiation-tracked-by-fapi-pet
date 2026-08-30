#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S20 - combined cross-cohort PROGENy + GSEA figure.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.
## REBUILD Supp 20 (was a PIL vstack of FigS17 a-f + GSEA panel g) as ONE vector patchwork.
## Sources both builders, regrabs panels, reassembles titleless at 180mm, exports quartz vector PDF.

OUTDIR <- file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/_fig4supp_audit_20260620/rebuilds")

## --- source the two panel builders (they build p_A..p_F, p_g in the global env;
##     their own exports are harmless side effects we ignore) ---
suppressWarnings(suppressMessages(
  source(file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "scripts/fig4_bulk_omics/build_figS17_progeny_cross_cohort.R"))
))
suppressWarnings(suppressMessages(
  source(file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "scripts/fig4_bulk_figures/build_supp20_gsea_panel_g_20260616.R"))
))

suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(ragg) })

## panel g needs its tag (the standalone script omits it)
p_g <- p_g + labs(tag = "g")

## a-f layout (verbatim from FigS17) + 2 full-width rows for g; NO plot_annotation title
## (supp title lives in the SI caption, per house convention).
design2 <- "
AAABBBCCC
AAABBBCCC
AAABBBCCC
DDDDEEEEE
DDDDEEEEE
DDDDEEEEE
FFFFFFFFF
FFFFFFFFF
GGGGGGGGG
GGGGGGGGG
"

combined <- p_A + p_B + p_C + p_D + p_E + p_F + p_g +
  plot_layout(design = design2) &
  theme(plot.tag.position = c(0.01, 1.02))

W_MM <- 180; H_MM <- 274   # ratio 0.657 = embedded Supp20 (5669x8637)
pdf_path <- file.path(OUTDIR, "Supp20_progeny_gsea_vector.pdf")
png_path <- file.path(OUTDIR, "Supp20_progeny_gsea_vector.png")
grDevices::quartz(file = pdf_path, type = "pdf", width = W_MM/25.4, height = H_MM/25.4, family = "Helvetica")
print(combined); dev.off()
ggsave(png_path, combined, width = W_MM, height = H_MM, units = "mm", dpi = 600, device = ragg::agg_png)
cat("Supp20 combined vector written:\n  ", pdf_path, "\n  ", png_path, "\n")
