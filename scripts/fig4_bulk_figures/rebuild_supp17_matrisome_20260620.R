#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S17 - matrisome/NNLS composition figure.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.
## REBUILD Supp17 matrisome (FigS10/fig_s10 from build_figS10_S11_supp.R) beautify 2026-06-20: quartz vector + 180mm. FigS11 export left intact (=Supp18 a-g, handled separately). -> rebuilds/
## Fig. S10 & Fig. S11 — Supplementary figures for WT vs FAP-deficient host
## xenograft stromal analysis.
##
## Fig S10 : Standardized stromal/CAF/Hallmark signature reference framework
##           (NABA matrisome + Elyada/Dominguez/Kieffer/Bartoschek CAF + Hallmark)
##           Rationale: the peer-reviewed backbone that motivated the B1
##           scaffold-gene analysis, kept for method transparency.
##
## Fig S11 : Full mechanism-refinement details — B1 (scaffold genes
##           across 3 signatures) + B2 (NNLS cell-fraction deconvolution
##           and ULM activity) combined into one supplementary page.


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales); library(grid)
})

proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
out_root <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
fig_dir <- file.path(out_root, "fig_S10_S11_supp")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- palette / theme ---------------------------------------------------
col_wt <- "#3B6B9C"; col_fap <- "#C9493A"
col_up <- col_fap;   col_dn <- col_wt
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
      axis.line     = element_line(linewidth = 0.3),
      axis.ticks    = element_line(linewidth = 0.3),
      strip.background = element_blank(),
      strip.text    = element_text(size = base),
      legend.key.size = unit(3.2, "mm"),
      legend.text   = element_text(size = base - 0.5),
      legend.title  = element_text(size = base - 0.5),
      panel.grid    = element_blank(),
      plot.margin   = margin(3, 4, 3, 4)
    )
}
p_label <- function(p){
  if (is.na(p)) return("ns")
  if (p < 0.001) sprintf("P = %.1e", p) else sprintf("P = %.3f", p)
}

## =======================================================================
## Load pre-computed tables
## =======================================================================
sig_stats   <- read_tsv(file.path(out_root, "naba_caf_signatures/signature_stats_by_group.tsv"),
                        show_col_types = FALSE)
sig_scores  <- read_tsv(file.path(out_root, "naba_caf_signatures/signature_scores_per_sample.tsv"),
                        show_col_types = FALSE)
core_tbl    <- read_tsv(file.path(out_root, "ecm_mycaf_core_genes/core_gene_stats.tsv"),
                        show_col_types = FALSE)
overlap_tbl <- read_tsv(file.path(out_root, "ecm_mycaf_core_genes/signature_overlap_table.tsv"),
                        show_col_types = FALSE)
nnls_df     <- read_tsv(file.path(out_root, "caf_deconvolution/nnls_cell_fractions.tsv"),
                        show_col_types = FALSE)
ulm_df      <- read_tsv(file.path(out_root, "caf_deconvolution/ulm_activity_scores.tsv"),
                        show_col_types = FALSE)
decon_stats <- read_tsv(file.path(out_root, "caf_deconvolution/deconvolution_group_stats.tsv"),
                        show_col_types = FALSE)

## =======================================================================
## FIG S10 — Standardized signature reference framework
## =======================================================================
heatmap_delta <- function(df, title_str, tag_str, subtitle_str = NULL){
  if (is.null(subtitle_str))
    subtitle_str <- "* P<0.05, . P<0.1 (Wilcoxon)"
  dset <- df %>%
    mutate(sig_tag = case_when(wilcox_p < 0.05 ~ "*",
                               wilcox_p < 0.1  ~ ".",
                               TRUE ~ ""),
           dataset = factor(dataset, levels = c("RNA", "Proteomics")))
  ord <- df %>% filter(dataset == "RNA") %>%
    arrange(delta_fap_minus_wt) %>% pull(signature)
  dset <- dset %>% mutate(signature = factor(signature, levels = rev(ord)),
                          label_disp = sprintf("%s%.2f",
                                                sig_tag, delta_fap_minus_wt))
  max_abs <- max(abs(dset$delta_fap_minus_wt), na.rm = TRUE)
  ggplot(dset, aes(dataset, signature, fill = delta_fap_minus_wt)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = label_disp,
                  colour = abs(delta_fap_minus_wt) > max_abs * 0.6),
              size = (base_font_size - 1.5)/.pt, show.legend = FALSE) +
    scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                         midpoint = 0, limits = c(-max_abs, max_abs),
                         oob = squish, name = "Delta\nFAP-KO\nminus WT") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, title = title_str, tag = tag_str,
         subtitle = subtitle_str) +
    theme_nc() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_text(size = base_font_size - 0.5),
          legend.key.height = unit(8, "mm"),
          legend.key.width  = unit(2.6, "mm"))
}

s10_a <- heatmap_delta(
  sig_stats %>% filter(collection == "NABA matrisome"),
  "NABA matrisome categories", "a",
  "Naba et al. 2012, MSigDB C2:CGP")

s10_b <- heatmap_delta(
  sig_stats %>% filter(collection == "CAF subtype (published)"),
  "CAF subtype signatures (published)", "b",
  "Elyada 2019 | Dominguez 2020 | Kieffer 2020 | Bartoschek 2018")

s10_c <- heatmap_delta(
  sig_stats %>% filter(collection == "MSigDB Hallmark"),
  "MSigDB Hallmark (fibrosis / EMT / inflammation context)", "c",
  "Liberzon et al. 2015 Cell Syst")

## Spotlight — top 4 signatures per-sample score
spotlight_sigs <- sig_stats %>%
  filter(!is.na(wilcox_p)) %>%
  group_by(collection) %>% arrange(wilcox_p) %>%
  slice_head(n = 2) %>% pull(signature) %>% unique() %>% head(6)
spot_df <- sig_scores %>%
  filter(signature %in% spotlight_sigs) %>%
  mutate(signature = factor(signature, levels = spotlight_sigs),
         dataset = factor(dataset, levels = c("RNA", "Proteomics")),
         group = factor(group, levels = grp_levels))

p_lab_d <- sig_stats %>% filter(signature %in% spotlight_sigs) %>%
  mutate(signature = factor(signature, levels = spotlight_sigs),
         dataset = factor(dataset, levels = c("RNA", "Proteomics")),
         label = vapply(wilcox_p, p_label, character(1)))

s10_d <- ggplot(spot_df, aes(group, score, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.8,
              shape = 21, stroke = 0.2, colour = "black") +
  geom_text(data = p_lab_d, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.3,
            size = (base_font_size - 1.5)/.pt) +
  facet_grid(signature ~ dataset, scales = "free_y", switch = "y") +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_x_discrete(labels = c("WT host" = "WT",
                              "FAP-deficient host" = "KO")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.25))) +
  labs(x = NULL, y = "ULM activity score",
       title = "Top standardized signatures per sample",
       subtitle = "Most-significant two signatures per collection",
       tag = "d") +
  theme_nc() +
  theme(strip.text.y.left = element_text(angle = 0,
                                         size = base_font_size - 1.5,
                                         face = "plain"),
        strip.placement = "outside",
        panel.spacing.y = unit(1.5, "mm"))

s10_design <- "
AABB
AABB
AABB
CCDD
CCDD
CCDD
"
fig_s10 <- s10_a + s10_b + s10_c + s10_d +
  plot_layout(design = s10_design) +
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    theme = theme(
      plot.title = element_text(face = "bold", size = base_font_size + 1),
      plot.subtitle = element_text(size = base_font_size - 0.5,
                                    colour = "grey30"))) &
  theme(plot.tag.position = c(0.01, 1.02))

## --- Supp17 vector export (FigS10=fig_s10) 180mm; was 220x240
grDevices::quartz(file = file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/_fig4supp_audit_20260620/rebuilds/Supp17_matrisome_vector.pdf"), type = "pdf", width = 180/25.4, height = 196/25.4, family = "Helvetica")
print(fig_s10); dev.off()
ggsave(file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/_fig4supp_audit_20260620/rebuilds/Supp17_matrisome_vector.png"), fig_s10, width = 180, height = 196, units = "mm", dpi = 600, device = ragg::agg_png)

cat("Fig S10 written.\n")

## =======================================================================
## FIG S11 — Full mechanism-refinement details (B1 + B2)
## =======================================================================

## ---- Panel a: 3-signature membership -----------------------------------
panel_a_df <- overlap_tbl %>%
  mutate(gene = factor(gene, levels = rev(overlap_tbl$gene)),
         group = case_when(
           n_sig == 3 ~ "3-way (scaffold core)",
           n_sig == 2 ~ "2-way (extended core)",
           TRUE ~ "unique")) %>%
  pivot_longer(c(Elyada, Dominguez, Kieffer),
               names_to = "signature", values_to = "present") %>%
  mutate(signature = factor(signature,
                            levels = c("Elyada", "Dominguez", "Kieffer")))

s11_a <- ggplot(panel_a_df, aes(signature, gene, fill = present)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = subset(panel_a_df, present),
             aes(colour = group), size = 1.1, shape = 16) +
  scale_fill_manual(values = c(`TRUE` = "#F0E4D7", `FALSE` = "white"),
                    guide = "none") +
  scale_colour_manual(values = c("3-way (scaffold core)" = col_fap,
                                  "2-way (extended core)" = col_purple,
                                  "unique" = col_grey),
                      name = NULL) +
  scale_x_discrete(position = "top",
                   labels = c("Elyada\n2019","Dominguez\n2020","Kieffer\n2020")) +
  labs(x = NULL, y = NULL,
       title = "Signature membership",
       subtitle = sprintf("%d / %d genes in 3 / 2+ of 3 signatures",
                          sum(overlap_tbl$n_sig == 3),
                          sum(overlap_tbl$n_sig >= 2)),
       tag = "a") +
  theme_nc() +
  theme(axis.text.y = element_text(face = "italic",
                                   size = base_font_size - 1),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(2.8, "mm"))

## ---- Panel b: scaffold gene forest (12 genes x 2 modalities) ------------
core12 <- core_tbl %>% filter(n_sig >= 2) %>%
  arrange(rna_log2fc) %>%
  mutate(gene_mgi = factor(gene_mgi, levels = gene_mgi))

b_df <- core12 %>%
  select(gene_mgi, n_sig,
         RNA = rna_log2fc, Protein = prot_log2fc) %>%
  pivot_longer(c(RNA, Protein), names_to = "modality",
               values_to = "log2fc") %>%
  mutate(modality = factor(modality, levels = c("RNA","Protein")))

tag_df_b <- bind_rows(
  core12 %>% select(gene_mgi, rna_q) %>% mutate(
    modality = "RNA",
    tag = case_when(rna_q < 0.05 ~ "*", rna_q < 0.1 ~ ".", TRUE ~ "")),
  core12 %>% select(gene_mgi, prot_p) %>% mutate(
    modality = "Protein",
    tag = case_when(prot_p < 0.05 ~ "*", prot_p < 0.1 ~ ".", TRUE ~ ""))
) %>%
  left_join(b_df %>% select(gene_mgi, modality, log2fc),
            by = c("gene_mgi","modality")) %>%
  mutate(modality = factor(modality, levels = c("RNA","Protein")))

s11_b <- ggplot(b_df, aes(log2fc, gene_mgi, fill = log2fc)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
  geom_col(colour = "black", linewidth = 0.2) +
  geom_text(data = tag_df_b,
            aes(x = log2fc + sign(log2fc) * 0.15,
                y = gene_mgi, label = tag),
            inherit.aes = FALSE,
            size = (base_font_size + 0.5)/.pt) +
  facet_wrap(~ modality, nrow = 1) +
  scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                       midpoint = 0, limits = c(-2.2, 2.2), oob = squish,
                       guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.15))) +
  labs(x = "log2 FC (FAP-KO vs WT)", y = NULL,
       title = "Scaffold genes (2+ signatures)",
       subtitle = "* FDR<0.05 or P<0.05; . <0.1",
       tag = "b") +
  theme_nc() +
  theme(axis.text.y = element_text(face = "italic",
                                   size = base_font_size - 0.5),
        panel.spacing.x = unit(2, "mm"))

## ---- Panel c: RNA vs Protein scatter with sign-concordance -------------
sc_df <- core_tbl %>% filter(n_sig >= 2,
                              !is.na(rna_log2fc), !is.na(prot_log2fc)) %>%
  mutate(class = case_when(n_sig == 3 ~ "3-way core",
                           n_sig == 2 ~ "2-way extended"))
bn <- binom.test(sum(sc_df$sign_concordant), nrow(sc_df),
                 p = 0.5, alternative = "greater")
rho_sc <- suppressWarnings(cor(sc_df$rna_log2fc, sc_df$prot_log2fc,
                                method = "spearman"))
max_abs_sc <- max(abs(c(sc_df$rna_log2fc, sc_df$prot_log2fc)),
                  na.rm = TRUE) * 1.05

s11_c <- ggplot(sc_df, aes(rna_log2fc, prot_log2fc)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              linewidth = 0.25, colour = "grey70") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey30", fill = "grey85", linewidth = 0.3) +
  geom_point(aes(fill = class), shape = 21, size = 2.2,
             stroke = 0.25, colour = "black") +
  geom_text_repel(aes(label = gene_mgi, colour = class),
                  size = (base_font_size - 1.5)/.pt,
                  min.segment.length = 0, segment.size = 0.2,
                  max.overlaps = Inf, box.padding = 0.3,
                  fontface = "italic", show.legend = FALSE) +
  scale_fill_manual(values = c("3-way core" = col_fap,
                               "2-way extended" = col_purple),
                    name = NULL) +
  scale_colour_manual(values = c("3-way core" = col_fap,
                                  "2-way extended" = col_purple)) +
  coord_cartesian(xlim = c(-max_abs_sc, max_abs_sc),
                  ylim = c(-max_abs_sc, max_abs_sc)) +
  labs(x = "RNA log2FC (FAP-KO vs WT)",
       y = "Protein log2FC (FAP-KO vs WT)",
       title = "Cross-modality agreement",
       subtitle = sprintf("Spearman %.2f | %d / %d sign-concordant (binomial P = %.2g)",
                          rho_sc, sum(sc_df$sign_concordant),
                          nrow(sc_df), bn$p.value),
       tag = "c") +
  theme_nc() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"))

## ---- Panel d: NNLS stacked bar per sample ------------------------------
samp_order <- c("WT_BPC1","WT_BPC2","WT_BPC3","WT_BPC4",
                "FAP_BPC1","FAP_BPC2","FAP_BPC3","FAP_BPC4")
ct_cols <- c("myCAF"        = "#C9493A",
             "iCAF"         = "#E08F3A",
             "apCAF"        = "#5E4A7B",
             "Thyroid_tumor"= "#3B6B9C",
             "Immune"       = "#88B0A7",
             "Endothelial"  = "#B7B7B7")
ct_levels <- c("myCAF","iCAF","apCAF","Thyroid_tumor","Immune","Endothelial")

s11_d <- ggplot(nnls_df %>%
                  mutate(sample = factor(sample, levels = samp_order),
                         celltype = factor(celltype, levels = ct_levels)),
                aes(sample, fraction, fill = celltype)) +
  geom_col(colour = "white", linewidth = 0.2, width = 0.9) +
  scale_fill_manual(values = ct_cols, name = NULL) +
  scale_y_continuous(expand = c(0, 0), labels = label_percent()) +
  labs(x = NULL, y = "Estimated cell fraction",
       title = "NNLS deconvolution",
       subtitle = "Marker-based, reference-free",
       tag = "d") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1,
                                   colour = ifelse(grepl("^WT", samp_order),
                                                   col_wt, col_fap)),
        legend.position = "right",
        legend.key.size = unit(2.4, "mm"))

## ---- Panel e: NNLS fraction box by cell type ---------------------------
box_df <- nnls_df %>%
  mutate(celltype = factor(celltype, levels = ct_levels),
         group = factor(group, levels = grp_levels))
box_stat <- decon_stats %>% filter(method == "NNLS fraction") %>%
  mutate(celltype = factor(celltype, levels = ct_levels),
         label = vapply(wilcox_p, p_label, character(1)))

s11_e <- ggplot(box_df, aes(group, fraction, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.7,
              shape = 21, stroke = 0.15, colour = "black") +
  geom_text(data = box_stat, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.3,
            size = (base_font_size - 1.8)/.pt) +
  facet_wrap(~ celltype, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_x_discrete(labels = c("WT host" = "WT",
                              "FAP-deficient host" = "KO")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)),
                     labels = label_percent(accuracy = 1)) +
  labs(x = NULL, y = "Fraction",
       title = "NNLS per-cell-type fraction by group",
       tag = "e") +
  theme_nc() +
  theme(strip.text = element_text(size = base_font_size - 1),
        panel.spacing.x = unit(1.5, "mm"))

## ---- Panel f: Cross-method consistency scatter -------------------------
cross_df <- decon_stats %>%
  select(celltype, delta_fap_minus_wt, wilcox_p, method) %>%
  pivot_wider(names_from = method,
              values_from = c(delta_fap_minus_wt, wilcox_p),
              names_sep = "_") %>%
  rename(delta_NNLS = `delta_fap_minus_wt_NNLS fraction`,
         delta_ULM  = `delta_fap_minus_wt_ULM activity`,
         p_NNLS = `wilcox_p_NNLS fraction`,
         p_ULM  = `wilcox_p_ULM activity`) %>%
  mutate(sig_any = (!is.na(p_NNLS) & p_NNLS < 0.1) |
                    (!is.na(p_ULM)  & p_ULM  < 0.1))

s11_f <- ggplot(cross_df, aes(delta_NNLS, delta_ULM)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_point(aes(fill = celltype, size = sig_any),
             shape = 21, stroke = 0.3, colour = "black") +
  geom_text_repel(aes(label = celltype, colour = celltype),
                  size = (base_font_size - 1)/.pt,
                  fontface = "bold",
                  min.segment.length = 0, segment.size = 0.25,
                  box.padding = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = ct_cols, guide = "none") +
  scale_colour_manual(values = ct_cols, guide = "none") +
  scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2.5),
                    labels = c(`TRUE` = "P<0.1 (any method)",
                               `FALSE` = "ns"),
                    name = NULL) +
  labs(x = "Delta (FAP-KO - WT) - NNLS fraction",
       y = "Delta (FAP-KO - WT) - ULM activity (z)",
       title = "Cross-method consistency",
       subtitle = "NNLS fraction vs orthogonal ULM activity",
       tag = "f") +
  theme_nc() +
  theme(legend.position = "top")

## ---- Panel g: ULM activity per cell type -------------------------------
ulm_stat <- decon_stats %>% filter(method == "ULM activity") %>%
  mutate(celltype = factor(celltype, levels = ct_levels),
         label = vapply(wilcox_p, p_label, character(1)))

s11_g <- ggplot(ulm_df %>%
                  mutate(celltype = factor(celltype, levels = ct_levels),
                         group = factor(group, levels = grp_levels)),
                aes(group, score, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.7,
              shape = 21, stroke = 0.15, colour = "black") +
  geom_text(data = ulm_stat, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.3,
            size = (base_font_size - 1.8)/.pt) +
  facet_wrap(~ celltype, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_x_discrete(labels = c("WT host" = "WT",
                              "FAP-deficient host" = "KO")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "ULM activity (z)",
       title = "ULM activity (orthogonal estimator)",
       tag = "g") +
  theme_nc() +
  theme(strip.text = element_text(size = base_font_size - 1),
        panel.spacing.x = unit(1.5, "mm"))

## ---- S11 assemble ------------------------------------------------------
s11_design <- "
AAABBBCCCC
AAABBBCCCC
AAABBBCCCC
DDDDDDFFFF
DDDDDDFFFF
EEEEEEEEEE
EEEEEEEEEE
GGGGGGGGGG
GGGGGGGGGG
"

fig_s11 <- s11_a + s11_b + s11_c + s11_d + s11_e + s11_f + s11_g +
  plot_layout(design = s11_design) +
  plot_annotation(
    title = "Fig. S11. Full mechanism-refinement evidence (scaffold gene + CAF deconvolution).",
    subtitle = "B1 scaffold-gene analysis (a-c) + B2 cell-type deconvolution (d-g). Supports Fig. 8 Panel e.",
    theme = theme(
      plot.title = element_text(face = "bold", size = base_font_size + 1),
      plot.subtitle = element_text(size = base_font_size - 0.5,
                                    colour = "grey30"))) &
  theme(plot.tag.position = c(0.01, 1.02))

ggsave(file.path(fig_dir, "FigS11_mechanism_refinement_full.pdf"),
       fig_s11, width = 230, height = 310, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(file.path(fig_dir, "FigS11_mechanism_refinement_full.png"),
       fig_s11, width = 230, height = 310, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("Fig S11 written.\n")

cat("\nAll supplementary figures in:\n  ", fig_dir, "\n")
