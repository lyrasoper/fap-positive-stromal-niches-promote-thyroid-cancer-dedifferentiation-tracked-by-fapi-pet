#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S17 - stromal ECM marker overview heatmap.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

script_path <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
}

outdir <- file.path(root, "outputs", "wt_vs_fapko_bulk_omics_20260407", "integrated")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

marker_sets <- list(
  "FAP-rich stroma" = c("FAP", "PDGFRA", "PDGFRB", "THY1", "PDPN"),
  "collagen-rich stroma" = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1"),
  "collagen/fibrotic ECM program" = c("COL1A1", "COL1A2", "FN1", "POSTN", "SPARC", "LOX"),
  "CAF/ECM fibrosis program" = c("FAP", "ACTA2", "TAGLN", "MYLK", "CTHRC1", "TGFBR1"),
  "fibroproliferative niche" = c("FAP", "POSTN", "TGFB1", "TGFBR1", "PDGFRB", "ITGA11"),
  "stromal matrix output" = c("DCN", "LUM", "BGN", "VCAN", "FN1", "THBS2"),
  "ECM remodeling axis" = c("MMP2", "MMP14", "TIMP2", "LOX", "LOXL2", "SERPINE1")
)

facet_levels <- names(marker_sets)

marker_map <- bind_rows(lapply(seq_along(marker_sets), function(i) {
  category <- facet_levels[i]
  genes <- unname(marker_sets[[i]])
  tibble(
    category = category,
    gene = genes,
    gene_order = seq_along(genes)
  )
}))

read_rna_diff <- function(path) {
  read.delim(path, check.names = FALSE) |>
    transmute(
      gene = toupper(gene_name),
      dataset = "RNA",
      log2_fc = as.numeric(`log2(fc)`),
      p_value = as.numeric(pval),
      q_value = as.numeric(qval)
    ) |>
    distinct(gene, .keep_all = TRUE)
}

read_protein_diff <- function(path) {
  read.delim(path, check.names = FALSE) |>
    transmute(
      gene = toupper(gene_name),
      dataset = "Proteomics",
      log2_fc = log2(as.numeric(FC)),
      p_value = as.numeric(p_value),
      q_value = as.numeric(q_value)
    ) |>
    distinct(gene, .keep_all = TRUE)
}

diff_df <- bind_rows(
  read_rna_diff(file.path(root, "001_mRNA_Summary", "3.DiffExpression", "FAP_BPCVSWT_BPC", "FAP_BPCVSWT_BPC_Gene_differential_expression.txt")),
  read_protein_diff(file.path(root, "002_DIA_Summary", "04.DiffExp", "COND1", "FAP_BPCVSWT_BPC", "FAP_BPCVSWT_BPC_diff_annotation.txt"))
)

plot_df <- marker_map |>
  tidyr::crossing(dataset = c("RNA", "Proteomics")) |>
  left_join(diff_df, by = c("gene", "dataset")) |>
  mutate(
    category = factor(category, levels = facet_levels),
    dataset = factor(dataset, levels = c("RNA", "Proteomics")),
    gene_facet = paste(category, gene, sep = "___"),
    sig_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.05 ~ "**",
      p_value < 0.10 ~ "*",
      TRUE ~ ""
    )
  ) |>
  mutate(
    gene_facet = factor(
      gene_facet,
      levels = rev(unique(paste(marker_map$category, marker_map$gene, sep = "___")))
    )
  ) |>
  ungroup()

wide_df <- plot_df |>
  transmute(category, gene = sub("^.*___", "", as.character(gene_facet)), dataset, log2_fc, p_value, q_value) |>
  pivot_wider(
    names_from = dataset,
    values_from = c(log2_fc, p_value, q_value),
    names_vary = "slowest"
  ) |>
  arrange(category, desc(gene))

write.table(
  wide_df,
  file.path(outdir, "19_stromal_ecm_marker_overview.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  marker_map,
  file.path(outdir, "19_stromal_ecm_marker_categories.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p <- ggplot(plot_df, aes(x = dataset, y = gene_facet, fill = log2_fc)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sig_label), size = 3.4) +
  scale_fill_gradient2(
    low = "#2563EB",
    mid = "white",
    high = "#DC2626",
    midpoint = 0,
    na.value = "grey88",
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    title = "Expanded stromal / ECM marker overview",
    subtitle = "* p < 0.10, ** p < 0.05; grey = not detected in that omics layer",
    x = NULL,
    y = NULL,
    fill = "log2FC\n(FAP/WT)"
  ) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey40"),
    strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 0),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9),
    panel.spacing.y = unit(0.3, "lines")
  )

p_compact <- ggplot(plot_df, aes(x = dataset, y = gene_facet, fill = log2_fc)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = sig_label), size = 2.9) +
  scale_fill_gradient2(
    low = "#5E81AC",
    mid = "#FBFBFA",
    high = "#D08770",
    midpoint = 0,
    na.value = "grey92",
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  labs(
    x = NULL,
    y = NULL,
    fill = "log2FC"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey96", colour = "grey55", linewidth = 0.5),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 8.1),
    axis.text.x = element_text(face = "bold", size = 8.5),
    axis.text.y = element_text(face = "bold", size = 7.8),
    legend.title = element_text(size = 8.2, face = "bold"),
    legend.text = element_text(size = 7.6),
    legend.key.height = unit(0.42, "in"),
    legend.key.width = unit(0.16, "in"),
    panel.spacing.y = unit(0.16, "lines"),
    plot.margin = margin(3, 4, 3, 4)
  )

ggsave(
  file.path(outdir, "19_stromal_ecm_marker_overview_heatmap.png"),
  p,
  width = 6.7,
  height = 12.5,
  dpi = 300
)
ggsave(
  file.path(outdir, "19_stromal_ecm_marker_overview_heatmap.pdf"),
  p,
  width = 6.7,
  height = 12.5
)

ggsave(
  file.path(outdir, "19_stromal_ecm_marker_overview_heatmap_compact.png"),
  p_compact,
  width = 5.9,
  height = 10.8,
  dpi = 300
)
ggsave(
  file.path(outdir, "19_stromal_ecm_marker_overview_heatmap_compact.pdf"),
  p_compact,
  width = 5.9,
  height = 10.8
)

message("Expanded stromal / ECM marker overview saved to: ", outdir)
