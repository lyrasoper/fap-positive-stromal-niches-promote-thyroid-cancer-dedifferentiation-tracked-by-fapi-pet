#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S15 (MOFA2 factors and ligand activities for
#          host Fap deficiency).
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

## REBUILD Supp15_mofa (src build_fig8_supp_mofa_nichenet.R) beautify 2026-06-20: quartz vector PDF + 180mm (was 220x240). Data identical. -> rebuilds/
## Fig. S8 | Unbiased multi-omics integration and stroma->tumor
## causal ligand prioritization.
##
## Combines MOFA2 Factor1 (cross-modality) with NicheNet stroma-lost
## and stroma-released ligand panels to support Fig. 8.

.libPaths(c(file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_omics/.Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
})

proj     <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_omics")
bulk_out <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
mofa_dir <- file.path(bulk_out, "mofa2")
nn_dir   <- file.path(bulk_out, "nichenet")
out_dir  <- file.path(bulk_out, "fig8_supp_panel")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- palette -----------------------------------------------------------
col_wt  <- "#3B6B9C"; col_fap <- "#C9493A"
col_up  <- col_fap;   col_dn  <- col_wt
col_grey <- "#B7B7B7"; col_purple <- "#4B1E70"
grp_levels <- c("WT host", "FAP-deficient host")
grp_cols <- setNames(c(col_wt, col_fap), grp_levels)

base_font_size <- 7
theme_nc <- function(base = base_font_size){
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(size = base + 1, hjust = 0),
      plot.subtitle = element_text(size = base - 0.5, colour = "grey30"),
      plot.tag      = element_text(face = "bold", size = base + 3, hjust = 0),
      axis.title    = element_text(size = base),
      axis.text     = element_text(size = base - 0.5, colour = "black"),
      axis.line     = element_line(linewidth = 0.6),
      axis.ticks    = element_line(linewidth = 0.6),
      strip.background = element_blank(),
      strip.text    = element_text(size = base),
      legend.key.size = unit(3.2, "mm"),
      legend.text   = element_text(size = base - 0.5),
      legend.title  = element_text(size = base - 0.5),
      panel.grid    = element_blank(),
      plot.margin   = margin(3, 4, 3, 4)
    )
}

tag_group <- function(df){
  df %>% mutate(group = factor(ifelse(grepl("^WT", sample), "WT host", "FAP-deficient host"),
                               levels = grp_levels))
}

p_label <- function(p){
  if (is.na(p)) return("ns")
  if (p < 0.001) sprintf("P = %.1e", p) else sprintf("P = %.3f", p)
}

## ---- load --------------------------------------------------------------
## NOTE: the MOFA2/NicheNet/dediff/BRS TSVs below are produced by upstream
## bulk-omics steps (MOFA2 fit, NicheNet ligand-activity, dediff-axis and
## TCGA-BRS scoring); see README. Guard each read.
read_tsv_guarded <- function(p, ...){
  if (!file.exists(p))
    stop(sprintf("Missing %s — produced by an upstream step; see README.", p))
  read_tsv(p, ...)
}

var_expl   <- read_tsv_guarded(file.path(mofa_dir, "variance_explained_per_factor.tsv"),
                       show_col_types = FALSE)
factor_df  <- read_tsv_guarded(file.path(mofa_dir, "factor_scores_per_sample.tsv"),
                       show_col_types = FALSE)
factor_stats <- read_tsv_guarded(file.path(mofa_dir, "factor_group_stats.tsv"),
                         show_col_types = FALSE)
factor_cor <- read_tsv_guarded(file.path(mofa_dir, "factor_phenotype_correlations.tsv"),
                       show_col_types = FALSE)

dediff_scores <- read_tsv_guarded(file.path(bulk_out, "dediff_axes/integrated/dediff_axes_scores_per_sample.tsv"),
                          show_col_types = FALSE)
brs_scores    <- read_tsv_guarded(file.path(bulk_out, "tcga_brs_official/bulk_rna/bulk_official_tcga_brs_scores.tsv"),
                          show_col_types = FALSE)

nn_down <- read_tsv_guarded(file.path(nn_dir, "primary_stromalost_ligand_ranking.tsv"),
                    show_col_types = FALSE)
nn_all <- read_tsv_guarded(file.path(nn_dir, "ligand_activities_all.tsv"),
                   show_col_types = FALSE)
nn_up <- nn_all %>%
  filter(analysis == "stroma-up ligands -> released tumor program") %>%
  arrange(desc(aupr_corrected))

## ============================================================
## Panel a — MOFA2 variance explained heatmap (compact)
## ============================================================
var_expl <- var_expl %>%
  mutate(factor = factor(factor, levels = paste0("Factor", 1:10)),
         view   = factor(view, levels = c("RNA", "Proteomics")))

p_a <- ggplot(var_expl, aes(view, factor, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.1f%%", pct),
                colour = pct > 18),
            size = (base_font_size - 0.8)/.pt, show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = col_purple,
                      name = "Var.\nexpl.",
                      labels = function(x) sprintf("%.0f%%", x)) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       tag = "a") +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key.height = unit(5, "mm"),
        legend.key.width  = unit(2.6, "mm"))

## ============================================================
## Panel b — Factor1 score WT vs FAP-KO boxplot (+ strip)
## ============================================================
fac1 <- factor_df %>% filter(factor == "Factor1") %>%
  mutate(group = factor(group, levels = grp_levels))
f1_stats <- factor_stats %>% filter(factor == "Factor1")
f1_lbl   <- sprintf("Delta = %+.2f\nWilcox %s",
                    f1_stats$delta_fap_minus_wt, p_label(f1_stats$wilcox_p))

p_b <- ggplot(fac1, aes(group, value, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 1.0,
              shape = 21, stroke = 0.25, colour = "black") +
  annotate("text", x = 1.5, y = max(fac1$value) + 0.3, label = f1_lbl,
           size = (base_font_size - 1.2)/.pt, lineheight = 1.05) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_x_discrete(labels = c("WT host" = "WT", "FAP-deficient host" = "FAP-KO")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "MOFA2 Factor1 score",
       tag = "b") +
  theme_nc()

## ============================================================
## Panel c — MOFA2 Factor1 vs TDS scatter (the KEY biology-link plot)
## ============================================================
tds_per_sample <- dediff_scores %>%
  filter(signature == "TDS / lineage identity", dataset == "RNA") %>%
  transmute(sample, tds_z = score_z)

fac1_tds <- fac1 %>% left_join(tds_per_sample, by = "sample") %>%
  filter(!is.na(tds_z))

r_spearman <- suppressWarnings(cor(fac1_tds$value, fac1_tds$tds_z, method = "spearman"))
r_pearson  <- suppressWarnings(cor(fac1_tds$value, fac1_tds$tds_z, method = "pearson"))
pct_p      <- suppressWarnings(cor.test(fac1_tds$value, fac1_tds$tds_z,
                                        method = "spearman", exact = FALSE)$p.value)

p_c <- ggplot(fac1_tds, aes(value, tds_z)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey45", fill = "grey80", linewidth = 0.4, alpha = 0.25) +
  geom_point(aes(fill = group), shape = 21, size = 2.4,
             stroke = 0.3, colour = "black") +
  geom_text_repel(aes(label = sub("FAP_BPC", "KO", sub("WT_BPC", "WT", sample)), colour = group),
                  size = (base_font_size - 1.3)/.pt,
                  segment.size = 0.2, min.segment.length = 0,
                  box.padding = 0.4, point.padding = 0.3,
                  max.overlaps = Inf, seed = 15, show.legend = FALSE) +
  annotate("text",
           x = max(fac1_tds$value), y = min(fac1_tds$tds_z),
           hjust = 1, vjust = 0,
           label = sprintf("Spearman %.2f (P = %.3f)\nPearson %.2f",
                           r_spearman, pct_p, r_pearson),
           size = (base_font_size - 1.5)/.pt, lineheight = 1.05,
           colour = "grey20") +
  scale_fill_manual(values = grp_cols, name = NULL) +
  scale_colour_manual(values = grp_cols, guide = "none") +
  labs(x = "MOFA2 Factor1 score (unbiased)",
       y = "16-gene TDS lineage z-score",
       tag = "c") +
  theme_nc() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"))

## ============================================================
## Panel d — Factor1 correlation ladder with all dediff axes + BRS
##
## NOTE: per the 2026-05-12 HGF removal pass, the "HGF-MET" entry is
## relabelled "MET (auxiliary)" on display — the receptor MET is kept as
## an auxiliary axis, but the HGF ligand axis is no longer claimed.
## Data row in factor_cor still carries the literal "HGF-MET" key so we
## remap it via a display-label vector below.
## ============================================================
focus_phen <- c("TDS / lineage identity", "BRAF-RAS score",
                "ERK score", "HGF-MET", "YAP/TAZ",
                "ECM-integrin-FAK", "Stemness / plasticity",
                "Notch signaling", "Hedgehog signaling",
                "TCGA-BRS (raw)")

display_lbl <- setNames(focus_phen, focus_phen)
display_lbl["HGF-MET"] <- "MET (auxiliary)"

f1_cor <- factor_cor %>% filter(factor == "Factor1", phen %in% focus_phen) %>%
  mutate(phen_disp = display_lbl[phen],
         phen_lbl = factor(phen_disp,
                            levels = rev(unname(display_lbl))),
         sig_tag = case_when(p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ ""))

p_d <- ggplot(f1_cor, aes(rho, phen_lbl, fill = rho)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_col(colour = "black", linewidth = 0.25) +
  geom_text(aes(label = sig_tag,
                x = rho + sign(rho) * 0.05),
            size = (base_font_size + 0.5)/.pt) +
  scale_fill_gradient2(low = col_dn, mid = "#F5F5F5", high = col_up,
                       midpoint = 0, limits = c(-1, 1),
                       oob = squish, guide = "none") +
  scale_x_continuous(limits = c(-1.15, 1.15),
                     breaks = seq(-1, 1, 0.5),
                     expand = expansion(mult = c(0.02, 0.08))) +
  ## significance-glyph KEY kept on-plot as a caption (not a descriptive subtitle)
  annotate("text", x = -1.15, y = 0.7, hjust = 0, vjust = 1,
           label = "* P < 0.05\n. P < 0.1",
           size = (base_font_size - 1.0)/.pt, lineheight = 1.05,
           colour = "grey20") +
  labs(x = "Spearman rho (Factor1 vs phenotype)",
       y = NULL,
       tag = "d") +
  theme_nc()

## ============================================================
## Panel e — NicheNet top 20 stroma-lost ligands (down in FAP-KO)
## ============================================================
top_dn <- nn_down %>% slice_head(n = 20) %>%
  mutate(test_ligand = factor(test_ligand, levels = rev(test_ligand)),
         prot_avail = ifelse(!is.na(prot_log2fc), "yes", "no"))

p_e <- ggplot(top_dn, aes(aupr_corrected, test_ligand)) +
  geom_col(aes(fill = rna_log2fc), colour = "black", linewidth = 0.2) +
  geom_point(data = subset(top_dn, prot_avail == "yes"),
             aes(x = aupr_corrected + max(aupr_corrected) * 0.04),
             shape = 8, size = 1.3, colour = "grey20") +
  scale_fill_gradient2(low = col_dn, mid = "#F5F5F5", high = col_up,
                       midpoint = 0, limits = c(-4, 4), oob = squish,
                       name = "RNA\nlog2FC") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  ## glyph KEY for the star marker kept on-plot
  annotate("text",
           x = max(top_dn$aupr_corrected), y = 1, hjust = 1, vjust = 0,
           label = "* = protein also measured",
           size = (base_font_size - 1.0)/.pt, colour = "grey20") +
  labs(x = "NicheNet ligand activity (AUPR corrected)",
       y = NULL,
       tag = "e") +
  theme_nc() +
  theme(axis.text.y = element_text(face = "italic",
                                   size = base_font_size - 0.5),
        legend.key.height = unit(6, "mm"),
        legend.key.width  = unit(2.4, "mm"))

## ============================================================
## Panel f — NicheNet top released ligands (UP in FAP-KO -> redifferentiation)
## ============================================================
top_up <- nn_up %>% slice_head(n = 20) %>%
  mutate(test_ligand = factor(test_ligand, levels = rev(test_ligand)))

# mark Wnt family for visual emphasis
top_up <- top_up %>%
  mutate(family = case_when(
    grepl("^Wnt", test_ligand) ~ "Wnt family",
    grepl("^Fgf", test_ligand) ~ "Fgf family",
    grepl("^Cxcl", test_ligand) ~ "Cxcl family",
    grepl("^Sema", test_ligand) ~ "Semaphorin",
    TRUE ~ "other"))

fam_cols <- c(`Wnt family` = "#C9493A", `Fgf family` = "#E2B200",
              `Cxcl family` = "#4B1E70", `Semaphorin` = "#165A42",
              other = "grey70")

p_f <- ggplot(top_up, aes(aupr_corrected, test_ligand, fill = family)) +
  geom_col(colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = fam_cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  labs(x = "NicheNet ligand activity (AUPR corrected)",
       y = NULL,
       tag = "f") +
  theme_nc() +
  theme(axis.text.y = element_text(face = "italic",
                                   size = base_font_size - 0.5),
        legend.position = "top", legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"))

## ============================================================
## assemble 6-panel Fig S8
## ============================================================
design <- "
AABBCCCC
AABBCCCC
DDDDCCCC
DDDDCCCC
EEEEFFFF
EEEEFFFF
EEEEFFFF
"

fig_supp <- p_a + p_b + p_c + p_d + p_e + p_f +
  plot_layout(design = design) &
  theme(plot.tag.position = c(0.01, 1.02))

## STRICT-NATURE: strip every panel title/subtitle globally (safety net over per-panel edits)
fig_supp <- fig_supp & labs(title = NULL, subtitle = NULL)

pdf_path <- file.path(out_dir, "FigS8_MOFA2_NicheNet_combined.pdf")
png_path <- file.path(out_dir, "FigS8_MOFA2_NicheNet_combined.png")

## --- vector export at TRUE SI print size (180mm); aspect preserved from 220x240
if (capabilities("cairo")) grDevices::cairo_pdf(pdf_path, width = 180/25.4, height = 196/25.4) else if (capabilities("aqua")) grDevices::quartz(file = pdf_path, type = "pdf", width = 180/25.4, height = 196/25.4) else grDevices::pdf(pdf_path, width = 180/25.4, height = 196/25.4)
print(fig_supp); dev.off()
ggsave(png_path, fig_supp, width = 180, height = 196, units = "mm", dpi = 600, device = ragg::agg_png)

cat("Fig S8 written:\n  ", pdf_path, "\n  ", png_path, "\n")
