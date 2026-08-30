#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S19 - thyroid-lineage promoter methylation in TCGA-THCA.
## Fig. S14 — Thyroid lineage-TF promoter methylation in TCGA-THCA.
##
## Goal: provide a methylation-level explanation for why
##   (a) terminal-state signature in Fig S12c did not fully reverse
##       in FAP-KO mouse (epigenetic lock hypothesis), and
##   (b) Factor1 signature anti-correlates with TDS in TCGA
##       (BRAF-like tumours are methylated / lineage-silenced).
##
## Method:
##   * TCGA-THCA HumanMethylation450K beta matrix from GDC/Xena
##   * Extract probes mapped to PAX8, NKX2-1, FOXE1, TG, TPO,
##     SLC5A5, TSHR, DIO1, DIO2, DUOX1, DUOX2 (thyroid lineage TFs
##     + iodide handling + thyroid-hormone synthesis)
##   * Per-sample per-gene mean promoter beta
##   * Integrate with Factor1 score + BRAF/RAS status + TDS from
##     prior TCGA analysis
##
## Panels:
##   a  Per-gene promoter methylation: driver class (BRAF/RAS/other)
##   b  PAX8 promoter methylation vs PAX8 mRNA scatter (TCGA)
##   c  PAX8 methylation vs Factor1 signature scatter
##   d  Multi-gene methylation heatmap across 3 driver groups
##   e  Methylation vs TDS scatter (inverse relationship)
##   f  Methylation x dedifferentiation ladder (TCGA + Landa)


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
})

proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
out_root <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
fig_dir <- file.path(out_root, "fig_S14_pax8_methylation")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

meth_dir <- file.path(proj, "reference/tcga_thca_methylation")

col_wt   <- "#3B6B9C"; col_fap  <- "#C9493A"; col_grey <- "#B7B7B7"
col_up   <- col_fap;   col_dn   <- col_wt
col_purple <- "#4B1E70"
col_gold <- "#E2B200"

base_font_size <- 7
theme_nc <- function(base = base_font_size){
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(size = base + 1, hjust = 0, face = "plain"),
      plot.subtitle = element_text(size = base - 0.5, colour = "grey30"),
      plot.tag      = element_text(face = "bold", size = base + 3, hjust = 0),
      axis.title    = element_text(size = base),
      axis.text     = element_text(size = base - 0.5, colour = "black"),
      axis.line     = element_line(linewidth = 0.6),
      axis.ticks    = element_line(linewidth = 0.6),
      strip.background = element_blank(),
      strip.text    = element_text(size = base, face = "plain"),
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

## ============================================================
## 1. Load beta matrix (thyroid-lineage probes only, ~305 probes)
## ============================================================
beta_mat_raw <- read_tsv(file.path(meth_dir, "thca_methylation_thyroid_probes.tsv"),
                         show_col_types = FALSE)
cat(sprintf("Beta matrix: %d probes x %d samples\n",
            nrow(beta_mat_raw), ncol(beta_mat_raw) - 1))

## probe annotation
probe_map <- read_tsv(file.path(meth_dir, "thyroid_lineage_probes.tsv"),
                      col_names = c("probe_id","gene","chrom","chromStart","chromEnd"),
                      show_col_types = FALSE)

## ============================================================
## 2. Build per-sample per-gene mean beta
## ============================================================
# shape: rows = probes, cols = samples (TCGA short ID format)
beta_long <- beta_mat_raw %>%
  rename(probe_id = 1) %>%
  pivot_longer(-probe_id, names_to = "sample_full", values_to = "beta") %>%
  mutate(sample = substr(sample_full, 1, 15)) %>%  # truncate to -01A etc
  left_join(probe_map %>% dplyr::select(probe_id, gene), by = "probe_id") %>%
  filter(!is.na(gene), !is.na(beta))

gene_beta <- beta_long %>%
  group_by(sample, gene) %>%
  summarise(mean_beta = mean(beta, na.rm = TRUE),
            n_probes = n(), .groups = "drop")

gene_beta_wide <- gene_beta %>%
  dplyr::select(-n_probes) %>%
  pivot_wider(names_from = gene, values_from = mean_beta,
              names_prefix = "meth_")

cat(sprintf("Per-sample methylation table: %d samples, %d genes\n",
            nrow(gene_beta_wide), sum(grepl("^meth_", colnames(gene_beta_wide)))))

## ============================================================
## 3. Integrate with Factor1 + TCGA annotations
## ============================================================
tcga_master <- read_tsv(file.path(out_root,
                                   "factor1_tcga_projection/factor1_tcga_merged_annot.tsv"),
                        show_col_types = FALSE)

# keep only primary tumour (01) — already in tcga_master
# Filter methylation to primary tumour (sample code 01)
primary_samples <- gene_beta_wide %>%
  mutate(sample_type = substr(sample, 14, 15)) %>%
  filter(sample_type == "01")

merged <- tcga_master %>%
  inner_join(primary_samples %>% dplyr::select(-sample_type),
             by = "sample")
cat(sprintf("Samples with methylation + Factor1 + annotation: %d\n",
            nrow(merged)))

# Add PAX8 mRNA (human gene symbol) from raw TCGA read - we need to re-read
# Instead get it from the factor1_tcga_scores.tsv? No — we need raw PAX8
# expression. Let's quickly rebuild log2-CPM just for PAX8 from the rds.
cat("Loading TCGA-THCA counts for PAX8 mRNA ...\n")

suppressPackageStartupMessages({
  library(edgeR); library(org.Hs.eg.db); library(AnnotationDbi)
})
tcga_counts <- readRDS(file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_public/TCGA-THCA.rds"))
ens_ids <- sub("\\..*$", "", rownames(tcga_counts))
sym_map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
                                                  keys = unique(ens_ids),
                                                  columns = "SYMBOL",
                                                  keytype = "ENSEMBL")) %>%
  as_tibble() %>% distinct(ENSEMBL, .keep_all = TRUE)
rowsym <- sym_map$SYMBOL[match(ens_ids, sym_map$ENSEMBL)]
keep <- !is.na(rowsym)
tcga_counts <- tcga_counts[keep, ]; rowsym <- rowsym[keep]
dup_syms <- rowsym[duplicated(rowsym)] %>% unique()
if (length(dup_syms) > 0){
  var_r <- apply(tcga_counts, 1, var)
  rank_within <- ave(-var_r, rowsym, FUN = rank)
  tcga_counts <- tcga_counts[rank_within == 1, , drop = FALSE]
  rowsym <- rowsym[rank_within == 1]
}
rownames(tcga_counts) <- rowsym
tcga_counts <- tcga_counts[, substr(colnames(tcga_counts), 14, 15) == "01"]
colnames(tcga_counts) <- substr(colnames(tcga_counts), 1, 15)
dge <- edgeR::DGEList(counts = tcga_counts)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
tcga_log <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

thyroid_genes <- c("PAX8","NKX2-1","FOXE1","TG","TPO","SLC5A5","TSHR",
                   "DIO1","DIO2","DUOX1","DUOX2")
thyroid_in <- intersect(thyroid_genes, rownames(tcga_log))
rna_df <- as.data.frame(t(tcga_log[thyroid_in, , drop = FALSE])) %>%
  tibble::rownames_to_column("sample")
colnames(rna_df)[-1] <- paste0("rna_", colnames(rna_df)[-1])

merged <- merged %>% inner_join(rna_df, by = "sample")

write_tsv(merged, file.path(fig_dir, "pax8_methylation_merged.tsv"))

## ============================================================
## 4. Per-gene methylation by driver class (stats)
## ============================================================
driver_levels <- c("BRAF V600E", "RAS driver", "other/none")
merged <- merged %>%
  mutate(brs_group = factor(brs_group, levels = driver_levels))

driver_stats <- list()
genes_to_test <- c("PAX8","NKX2-1","FOXE1","TG","TPO","SLC5A5","TSHR",
                   "DIO1","DIO2","DUOX1","DUOX2")
for (g in genes_to_test){
  col_name <- paste0("meth_", g)
  if (!col_name %in% colnames(merged)) next
  df <- merged %>% dplyr::select(brs_group, beta = all_of(col_name)) %>%
    filter(!is.na(beta))
  braf <- df %>% filter(brs_group == "BRAF V600E") %>% pull(beta)
  ras  <- df %>% filter(brs_group == "RAS driver") %>% pull(beta)
  if (length(braf) >= 3 && length(ras) >= 3){
    wp <- suppressWarnings(wilcox.test(braf, ras)$p.value)
  } else wp <- NA
  driver_stats[[g]] <- tibble(
    gene = g,
    n_BRAF = length(braf),
    n_RAS  = length(ras),
    mean_BRAF = mean(braf, na.rm = TRUE),
    mean_RAS  = mean(ras,  na.rm = TRUE),
    delta_BRAF_minus_RAS = mean_BRAF - mean_RAS,
    wilcox_p = wp)
}
driver_stats_tbl <- bind_rows(driver_stats)
write_tsv(driver_stats_tbl, file.path(fig_dir, "promoter_methylation_by_driver_stats.tsv"))

## ============================================================
## Panels
## ============================================================
## Panel a — promoter methylation by driver for all 11 genes
box_long <- merged %>%
  dplyr::select(sample, brs_group, starts_with("meth_")) %>%
  pivot_longer(-c(sample, brs_group),
               names_to = "gene", values_to = "beta") %>%
  mutate(gene = sub("^meth_", "", gene),
         gene = factor(gene, levels = genes_to_test))

stat_a <- driver_stats_tbl %>%
  transmute(gene = factor(gene, levels = genes_to_test),
            lbl = vapply(wilcox_p, p_label, character(1)))

p_A <- ggplot(box_long %>% filter(!is.na(beta)),
              aes(brs_group, beta, fill = brs_group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.3,
              shape = 21, stroke = 0.1, colour = "black", alpha = 0.6) +
  geom_text(data = stat_a, aes(x = 1.5, y = Inf, label = lbl),
            inherit.aes = FALSE, vjust = 1.3,
            size = (base_font_size - 1.8)/.pt) +
  facet_wrap(~ gene, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c("BRAF V600E" = col_wt,
                               "RAS driver" = col_fap,
                               "other/none" = col_grey),
                    guide = "none") +
  scale_x_discrete(labels = c("BRAF V600E" = "BRAF",
                              "RAS driver" = "RAS",
                              "other/none" = "other")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = NULL, y = "Mean promoter methylation (beta)",
       title = "Thyroid lineage TF promoter methylation by driver class",
       subtitle = "Wilcoxon BRAF vs RAS per gene",
       tag = "a") +
  theme_nc() +
  theme(strip.text = element_text(face = "italic",
                                   size = base_font_size - 0.5))

## Panel b — PAX8 methylation vs PAX8 mRNA scatter
if ("meth_PAX8" %in% colnames(merged) && "rna_PAX8" %in% colnames(merged)){
  scatter_pax8 <- merged %>%
    dplyr::select(sample, brs_group, meth_PAX8, rna_PAX8) %>%
    filter(!is.na(meth_PAX8), !is.na(rna_PAX8))
  cor_test_pax8 <- cor.test(scatter_pax8$meth_PAX8, scatter_pax8$rna_PAX8,
                             method = "spearman", exact = FALSE)
  p_B <- ggplot(scatter_pax8, aes(meth_PAX8, rna_PAX8)) +
    geom_point(aes(fill = brs_group), shape = 21,
               size = 1.0, stroke = 0.2, colour = "black", alpha = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "grey25", fill = "grey85", linewidth = 0.4) +
    scale_fill_manual(values = c("BRAF V600E" = col_wt,
                                 "RAS driver" = col_fap,
                                 "other/none" = col_grey),
                      name = NULL) +
    labs(x = "PAX8 promoter methylation (beta)",
         y = "PAX8 mRNA (log2 CPM)",
         title = "PAX8 methylation vs PAX8 mRNA",
         subtitle = sprintf("Spearman rho = %.2f, %s",
                            cor_test_pax8$estimate,
                            p_label(cor_test_pax8$p.value))) +
    theme_nc() +
    theme(plot.title = element_text(size = base_font_size, face = "bold",
                                     hjust = 0),
          plot.subtitle = element_text(size = base_font_size - 1.5),
          legend.position = "top", legend.direction = "horizontal",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(1, "mm"))
} else {
  p_B <- ggplot() + theme_void() +
    annotate("text", x = 1, y = 1, label = "PAX8 methylation missing") +
    labs(tag = "b")
}

## Panel c — PAX8 methylation vs Factor1 signature
if ("meth_PAX8" %in% colnames(merged)){
  scatter_pf <- merged %>%
    dplyr::select(sample, brs_group, meth_PAX8, score_weighted) %>%
    filter(!is.na(meth_PAX8), !is.na(score_weighted))
  cor_test_pf <- cor.test(scatter_pf$meth_PAX8, scatter_pf$score_weighted,
                           method = "spearman", exact = FALSE)
  p_C <- ggplot(scatter_pf, aes(meth_PAX8, score_weighted)) +
    geom_point(aes(fill = brs_group), shape = 21,
               size = 1.0, stroke = 0.2, colour = "black", alpha = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "grey25", fill = "grey85", linewidth = 0.4) +
    scale_fill_manual(values = c("BRAF V600E" = col_wt,
                                 "RAS driver" = col_fap,
                                 "other/none" = col_grey),
                      guide = "none") +
    labs(x = "PAX8 promoter methylation (beta)",
         y = "Factor1 signature score",
         title = "PAX8 methylation tracks Factor1 axis",
         subtitle = sprintf("Spearman rho = %.2f, %s",
                            cor_test_pf$estimate,
                            p_label(cor_test_pf$p.value))) +
    theme_nc() +
    theme(plot.title = element_text(size = base_font_size, face = "bold",
                                     hjust = 0),
          plot.subtitle = element_text(size = base_font_size - 1.5))
} else {
  p_C <- ggplot() + theme_void() +
    annotate("text", x = 1, y = 1, label = "PAX8 methylation missing") +
    labs(tag = "c")
}

## Panel d — Multi-gene methylation heatmap by driver (group-mean z per gene)
hm_df <- merged %>%
  dplyr::select(brs_group, starts_with("meth_")) %>%
  pivot_longer(-brs_group, names_to = "gene", values_to = "beta") %>%
  mutate(gene = sub("^meth_", "", gene)) %>%
  filter(!is.na(beta)) %>%
  group_by(gene, brs_group) %>%
  summarise(mean_beta = mean(beta), .groups = "drop") %>%
  group_by(gene) %>%
  mutate(z = scale(mean_beta)[, 1]) %>% ungroup() %>%
  mutate(gene = factor(gene, levels = genes_to_test),
         brs_group = factor(brs_group,
                            levels = c("BRAF V600E", "other/none", "RAS driver")))

p_D <- ggplot(hm_df, aes(brs_group, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", mean_beta),
                colour = abs(z) > 0.5),
            size = (base_font_size - 1.8)/.pt, show.legend = FALSE) +
  scale_fill_gradient2(low = col_dn, mid = "#F5F5F5", high = col_up,
                       midpoint = 0, oob = squish,
                       name = "z (group mean)") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  labs(x = NULL, y = NULL,
       title = "Thyroid-lineage methylation by driver (group mean)",
       subtitle = "higher beta = more methylated = silenced",
       tag = "d") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key.height = unit(6, "mm"),
        legend.key.width  = unit(2.4, "mm"))

## Panel e — PAX8 methylation vs TDS
if ("meth_PAX8" %in% colnames(merged)){
  scatter_tds <- merged %>%
    dplyr::select(sample, brs_group, meth_PAX8, tds) %>%
    filter(!is.na(meth_PAX8), !is.na(tds))
  cor_test_tds <- cor.test(scatter_tds$meth_PAX8, scatter_tds$tds,
                            method = "spearman", exact = FALSE)
  p_E <- ggplot(scatter_tds, aes(meth_PAX8, tds)) +
    geom_point(aes(fill = brs_group), shape = 21,
               size = 1.0, stroke = 0.2, colour = "black", alpha = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "grey25", fill = "grey85", linewidth = 0.4) +
    scale_fill_manual(values = c("BRAF V600E" = col_wt,
                                 "RAS driver" = col_fap,
                                 "other/none" = col_grey),
                      guide = "none") +
    labs(x = "PAX8 promoter methylation (beta)",
         y = "TDS (16-gene thyroid differentiation)",
         title = "Methylation anti-correlates with TDS",
         subtitle = sprintf("Spearman rho = %.2f, %s",
                            cor_test_tds$estimate,
                            p_label(cor_test_tds$p.value))) +
    theme_nc() +
    theme(plot.title = element_text(size = base_font_size, face = "bold",
                                     hjust = 0),
          plot.subtitle = element_text(size = base_font_size - 1.5))
} else {
  p_E <- ggplot() + theme_void() + labs(tag = "e")
}

## Panel f — mean 11-gene promoter methylation composite score by driver
merged <- merged %>%
  mutate(mean_lineage_meth = rowMeans(
    dplyr::select(., starts_with("meth_")), na.rm = TRUE))

comp_df <- merged %>%
  dplyr::select(brs_group, mean_lineage_meth) %>%
  filter(!is.na(mean_lineage_meth))
comp_braf <- comp_df %>% filter(brs_group == "BRAF V600E") %>% pull(mean_lineage_meth)
comp_ras  <- comp_df %>% filter(brs_group == "RAS driver") %>% pull(mean_lineage_meth)
wp_comp <- suppressWarnings(wilcox.test(comp_braf, comp_ras)$p.value)

p_F <- ggplot(comp_df, aes(brs_group, mean_lineage_meth, fill = brs_group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.35,
              shape = 21, stroke = 0.12, colour = "black", alpha = 0.75) +
  annotate("text", x = 1.5, y = max(comp_df$mean_lineage_meth) * 1.02,
           label = sprintf("BRAF vs RAS %s", p_label(wp_comp)),
           size = (base_font_size - 1.2)/.pt) +
  scale_fill_manual(values = c("BRAF V600E" = col_wt,
                               "RAS driver" = col_fap,
                               "other/none" = col_grey),
                    guide = "none") +
  scale_x_discrete(labels = c("BRAF V600E" = "BRAF",
                              "RAS driver" = "RAS",
                              "other/none" = "other")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "Mean promoter beta\n(11 thyroid-lineage genes)",
       title = "Composite lineage-silencing score",
       subtitle = sprintf("BRAF-like tumours show elevated lineage-TF methylation (n=%d/%d/%d)",
                          length(comp_braf), length(comp_ras),
                          sum(comp_df$brs_group == "other/none"))) +
  theme_nc() +
  theme(plot.title = element_text(size = base_font_size, face = "bold",
                                   hjust = 0),
        plot.subtitle = element_text(size = base_font_size - 1.5))

## ============================================================
## assemble Fig S14
## ============================================================
design <- "
AAAAAA
AAAAAA
AAAAAA
BBBCCC
BBBCCC
DDEEFF
DDEEFF
"

## In-figure plot title + subtitle stripped per scripts-122/123 convention
## (Supplementary figure title lives in the SI docx caption, not on the figure).
## The old title "Fig. S14. ..." also referenced stale pre-2026-05-12 numbering.
fig_s14 <- p_A + p_B + p_C + p_D + p_E + p_F +
  plot_layout(design = design) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(fig_dir, "FigS14_pax8_methylation.pdf")
png_path <- file.path(fig_dir, "FigS14_pax8_methylation.png")

ggsave(pdf_path, fig_s14, width = 240, height = 310, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(png_path, fig_s14, width = 240, height = 310, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("\nFig S14 written:\n  ", pdf_path, "\n  ", png_path, "\n")

cat("\n==== Per-gene driver stats ====\n")
print(driver_stats_tbl %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))

cat(sprintf("\nComposite 11-gene methylation: BRAF mean = %.3f, RAS mean = %.3f, P = %.2e\n",
            mean(comp_braf, na.rm = TRUE),
            mean(comp_ras, na.rm = TRUE),
            wp_comp))
