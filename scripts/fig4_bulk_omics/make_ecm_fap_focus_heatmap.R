#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Figure 4E - ECM/FAP focus heatmap.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
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

focus_markers <- c(
  "FAP",
  "ACTA2",
  "COL1A1",
  "COL1A2",
  "FN1",
  "POSTN",
  "SPARC",
  "TGFB1",
  "TGFBR1"
)

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

plot_df <- bind_rows(
  read_rna_diff(file.path(root, "001_mRNA_Summary", "3.DiffExpression", "FAP_BPCVSWT_BPC", "FAP_BPCVSWT_BPC_Gene_differential_expression.txt")),
  read_protein_diff(file.path(root, "002_DIA_Summary", "04.DiffExp", "COND1", "FAP_BPCVSWT_BPC", "FAP_BPCVSWT_BPC_diff_annotation.txt"))
) |>
  filter(gene %in% focus_markers) |>
  mutate(
    gene = factor(gene, levels = rev(focus_markers)),
    dataset = factor(dataset, levels = c("RNA", "Proteomics")),
    sig_label = case_when(
      q_value < 0.05 ~ "**",
      q_value < 0.25 ~ "*",
      TRUE ~ ""
    )
  )

write.table(
  plot_df,
  file.path(outdir, "18_ecm_fap_focused_marker_changes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p <- ggplot(plot_df, aes(x = dataset, y = gene, fill = log2_fc)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sig_label), size = 4.1) +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#DC2626",
    midpoint = 0,
    na.value = "grey88"
  ) +
  labs(
    title = "Focused ECM / collagen / FAP marker changes",
    subtitle = "* q < 0.25, ** q < 0.05",
    x = NULL,
    y = NULL,
    fill = "log2FC\n(FAP/WT)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

ggsave(file.path(outdir, "18_ecm_fap_focused_marker_heatmap.png"), p, width = 4.8, height = 5.8, dpi = 300)
ggsave(file.path(outdir, "18_ecm_fap_focused_marker_heatmap.pdf"), p, width = 4.8, height = 5.8)

message("Focused ECM / collagen / FAP heatmap saved to: ", outdir)
