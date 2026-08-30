#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Figure 4D (MAPK/ERK output row) - PROGENy pathway footprints.

script_path <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
}

project_lib <- file.path(root, ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggrepel)
  library(ggplot2)
  library(pheatmap)
  library(progeny)
  library(RColorBrewer)
  library(tidyr)
})

theme_set(
  theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      axis.text = element_text(colour = "black")
    )
)

outdir <- file.path(root, "outputs", "wt_vs_fapko_bulk_omics_20260407", "progeny")
rna_dir <- file.path(outdir, "rna")
proteomics_dir <- file.path(outdir, "proteomics")
integrated_dir <- file.path(outdir, "integrated")
dir.create(rna_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proteomics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(integrated_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("WT host", "FAP-deficient host")
group_colors <- c("WT host" = "#4C78A8", "FAP-deficient host" = "#E45756")
pathway_levels <- c(
  "Androgen", "EGFR", "Estrogen", "Hypoxia", "JAK-STAT", "MAPK", "NFkB",
  "p53", "PI3K", "TGFb", "TNFa", "Trail", "VEGF", "WNT"
)

label_group <- function(x) {
  clean_x <- sub("^FPKM\\.", "", x)
  ifelse(grepl("^FAP_BPC", clean_x), "FAP-deficient host", "WT host")
}

clean_sample_name <- function(x) {
  sub("^FPKM\\.", "", x)
}

collapse_by_gene <- function(df, gene_col, sample_cols) {
  df |>
    group_by(.data[[gene_col]]) |>
    summarise(across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

read_sample_info <- function(path) {
  info <- read.delim(path, check.names = FALSE, comment.char = "")
  names(info)[1] <- "RawID"
  info
}

read_rna_expression <- function(path) {
  df <- read.delim(path, check.names = FALSE)
  sample_cols <- grep("^FPKM\\.", names(df), value = TRUE)
  df2 <- df |>
    select(gene_name, all_of(sample_cols)) |>
    mutate(across(all_of(sample_cols), as.numeric)) |>
    filter(!is.na(gene_name), gene_name != "")
  collapse_by_gene(df2, "gene_name", sample_cols)
}

read_protein_expression <- function(path, sample_info) {
  df <- read.delim(path, check.names = FALSE)
  sample_cols <- grep("^0208JLYAstral_DIA_zyl_P4_", names(df), value = TRUE)
  raw_ids <- sub("^0208JLYAstral_DIA_zyl_P4_", "", sample_cols)
  mapped_names <- sample_info$SampleID[match(raw_ids, sample_info$RawID)]
  rename_map <- setNames(sample_cols, mapped_names)
  df2 <- df |>
    mutate(gene_name = sub(";.*$", "", Genes)) |>
    select(gene_name, all_of(sample_cols)) |>
    rename(!!!rename_map) |>
    mutate(across(-gene_name, as.numeric)) |>
    filter(!is.na(gene_name), gene_name != "")
  collapse_by_gene(df2, "gene_name", mapped_names)
}

prepare_matrix <- function(expr_df, gene_col, sample_cols) {
  mat <- as.matrix(expr_df[, sample_cols, drop = FALSE])
  rownames(mat) <- expr_df[[gene_col]]
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  mat <- log2(mat + 1)
  # PROGENy cannot score pathways when input genes still contain missing values.
  mat[rowSums(!is.finite(mat)) == 0, , drop = FALSE]
}

safe_scale_cols <- function(mat) {
  out <- scale(mat)
  out[!is.finite(out)] <- 0
  out
}

write_table_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

p_to_symbol <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  }
  if (p_value < 1e-4) {
    return("****")
  }
  if (p_value < 1e-3) {
    return("***")
  }
  if (p_value < 1e-2) {
    return("**")
  }
  if (p_value < 0.05) {
    return("*")
  }
  "ns"
}

plot_boxplots <- function(score_df, stats_df, dataset_name, path_png, path_pdf) {
  ann <- stats_df |>
    mutate(
      sig_label = vapply(p_value, p_to_symbol, character(1))
    ) |>
    select(pathway, sig_label)

  bounds <- score_df |>
    group_by(pathway) |>
    summarise(
      y_min = min(score_z, na.rm = TRUE),
      y_max = max(score_z, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      span = pmax(y_max - y_min, 0.2),
      bracket_y = y_max + span * 0.12,
      bracket_top = bracket_y + span * 0.06,
      label_y = bracket_top + span * 0.05
    )

  ann <- left_join(ann, bounds, by = "pathway")

  plot_colors <- c("WT host" = "#D47A6A", "FAP-deficient host" = "#46B58A")

  p <- ggplot(score_df, aes(x = group, y = score_z)) +
    geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.56, alpha = 0.92, linewidth = 0.7, colour = "#2C2C2C") +
    geom_point(
      aes(fill = group),
      shape = 21,
      size = 2.1,
      stroke = 0.35,
      colour = "#2C2C2C",
      position = position_jitter(width = 0.10, height = 0)
    ) +
    geom_segment(
      data = ann,
      inherit.aes = FALSE,
      aes(x = 1, xend = 1, y = bracket_y, yend = bracket_top),
      linewidth = 0.6,
      colour = "#2C2C2C"
    ) +
    geom_segment(
      data = ann,
      inherit.aes = FALSE,
      aes(x = 1, xend = 2, y = bracket_top, yend = bracket_top),
      linewidth = 0.6,
      colour = "#2C2C2C"
    ) +
    geom_segment(
      data = ann,
      inherit.aes = FALSE,
      aes(x = 2, xend = 2, y = bracket_top, yend = bracket_y),
      linewidth = 0.6,
      colour = "#2C2C2C"
    ) +
    geom_text(
      data = ann,
      inherit.aes = FALSE,
      aes(x = 1.5, y = label_y, label = sig_label),
      size = 3.1,
      fontface = "bold"
    ) +
    facet_wrap(~ pathway, scales = "free_y", ncol = 7) +
    scale_fill_manual(values = plot_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.24))) +
    labs(
      title = paste0(dataset_name, " PROGENy"),
      x = NULL,
      y = "PROGENy pathway score z"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1, face = "bold", size = 8.8),
      axis.text.y = element_text(size = 8.2),
      axis.title.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0, size = 12),
      strip.text = element_text(face = "bold", size = 8.7, margin = margin(4, 5, 4, 5)),
      strip.background = element_rect(fill = "white", colour = "#2C2C2C", linewidth = 0.8),
      panel.spacing = unit(0.65, "lines"),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.6, colour = "#2C2C2C")
    )

  ggsave(path_png, p, width = 20, height = 7.5, dpi = 300)
  ggsave(path_pdf, p, width = 20, height = 7.5)
}

plot_heatmap_dual <- function(mat, annotation_col, title, path_png, path_pdf, width = 7.8, height = 5.8, cluster_rows = FALSE, cluster_cols = FALSE) {
  pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(101)
  ann_colors <- list(Group = group_colors)

  png(path_png, width = width, height = height, units = "in", res = 300)
  pheatmap(
    mat,
    color = pal,
    breaks = seq(-2.5, 2.5, length.out = 102),
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    main = title,
    fontsize_row = 10,
    fontsize_col = 9,
    border_color = "white"
  )
  dev.off()

  pdf(path_pdf, width = width, height = height)
  pheatmap(
    mat,
    color = pal,
    breaks = seq(-2.5, 2.5, length.out = 102),
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    main = title,
    fontsize_row = 10,
    fontsize_col = 9,
    border_color = "white"
  )
  dev.off()
}

plot_delta_concordance <- function(delta_df, path_png, path_pdf) {
  rho <- if (stats::sd(delta_df$rna_delta_z, na.rm = TRUE) == 0 || stats::sd(delta_df$proteomics_delta_z, na.rm = TRUE) == 0) {
    NA_real_
  } else {
    cor(delta_df$rna_delta_z, delta_df$proteomics_delta_z, method = "spearman")
  }

  p <- ggplot(delta_df, aes(rna_delta_z, proteomics_delta_z, label = pathway)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(size = 3.2, colour = "#7C3AED", alpha = 0.85) +
    ggrepel::geom_text_repel(size = 3.4, show.legend = FALSE) +
    annotate(
      "text",
      x = min(delta_df$rna_delta_z, na.rm = TRUE),
      y = max(delta_df$proteomics_delta_z, na.rm = TRUE),
      label = paste0("Spearman rho = ", ifelse(is.na(rho), "NA", sprintf("%.2f", rho))),
      hjust = 0,
      vjust = 1,
      size = 3.5
    ) +
    labs(
      title = "PROGENy delta concordance across RNA and proteomics",
      subtitle = "Delta = mean z-score in FAP-deficient host minus WT host",
      x = "RNA delta z",
      y = "Proteomics delta z"
    )

  ggsave(path_png, p, width = 6.8, height = 5.8, dpi = 300)
  ggsave(path_pdf, p, width = 6.8, height = 5.8)
}

format_top_pathways <- function(stats_df, direction = c("up", "down"), n = 5) {
  direction <- match.arg(direction)
  ordered <- if (direction == "up") {
    stats_df |>
      arrange(desc(delta_fap_minus_wt_z), q_value)
  } else {
    stats_df |>
      arrange(delta_fap_minus_wt_z, q_value)
  }

  paste(
    ordered |>
      slice_head(n = n) |>
      transmute(label = paste0(pathway, "(delta=", sprintf("%.2f", delta_fap_minus_wt_z), ", q=", formatC(q_value, format = "e", digits = 1), ")")) |>
      pull(label),
    collapse = "; "
  )
}

run_progeny_modality <- function(expr_df, gene_col, sample_cols, dataset_name, output_dir) {
  expr_mat <- prepare_matrix(expr_df, gene_col, sample_cols)
  raw_scores <- as.data.frame(progeny::progeny(expr_mat, scale = TRUE, organism = "Mouse", top = 100, perm = 1, z_scores = FALSE))
  raw_scores <- raw_scores[, pathway_levels, drop = FALSE]
  z_scores <- as.data.frame(safe_scale_cols(as.matrix(raw_scores)))
  z_scores <- z_scores[, pathway_levels, drop = FALSE]

  sample_meta <- tibble(
    sample_raw = rownames(raw_scores),
    sample = clean_sample_name(sample_raw),
    group = factor(label_group(sample_raw), levels = group_levels)
  ) |>
    arrange(group, sample)

  raw_scores <- raw_scores[sample_meta$sample_raw, , drop = FALSE]
  z_scores <- z_scores[sample_meta$sample_raw, , drop = FALSE]

  raw_long <- raw_scores |>
    tibble::rownames_to_column("sample_raw") |>
    mutate(sample = clean_sample_name(sample_raw), group = factor(label_group(sample_raw), levels = group_levels)) |>
    pivot_longer(all_of(pathway_levels), names_to = "pathway", values_to = "score_raw")

  z_long <- z_scores |>
    tibble::rownames_to_column("sample_raw") |>
    mutate(sample = clean_sample_name(sample_raw), group = factor(label_group(sample_raw), levels = group_levels)) |>
    pivot_longer(all_of(pathway_levels), names_to = "pathway", values_to = "score_z")

  score_df <- left_join(raw_long, z_long, by = c("sample_raw", "sample", "group", "pathway")) |>
    mutate(dataset = dataset_name)

  stats_df <- score_df |>
    group_by(pathway) |>
    summarise(
      wt_mean_raw = mean(score_raw[group == "WT host"], na.rm = TRUE),
      fap_mean_raw = mean(score_raw[group == "FAP-deficient host"], na.rm = TRUE),
      delta_fap_minus_wt_raw = fap_mean_raw - wt_mean_raw,
      wt_mean_z = mean(score_z[group == "WT host"], na.rm = TRUE),
      fap_mean_z = mean(score_z[group == "FAP-deficient host"], na.rm = TRUE),
      delta_fap_minus_wt_z = fap_mean_z - wt_mean_z,
      p_value = tryCatch(wilcox.test(score_raw ~ group, exact = FALSE)$p.value, error = function(e) NA_real_),
      .groups = "drop"
    ) |>
    mutate(
      q_value = p.adjust(p_value, method = "BH"),
      dataset = dataset_name
    )

  model_mat <- as.matrix(progeny::getModel(organism = "Mouse", top = 100))
  coverage_df <- tibble(
    pathway = colnames(model_mat),
    genes_in_model = colSums(model_mat != 0),
    genes_present = colSums(model_mat[rownames(model_mat) %in% rownames(expr_mat), , drop = FALSE] != 0),
    coverage_pct = 100 * genes_present / genes_in_model,
    dataset = dataset_name
  ) |>
    arrange(match(pathway, pathway_levels))

  sample_ann <- data.frame(Group = sample_meta$group)
  rownames(sample_ann) <- sample_meta$sample

  sample_heatmap_mat <- t(as.matrix(z_scores))
  colnames(sample_heatmap_mat) <- sample_meta$sample
  rownames(sample_heatmap_mat) <- pathway_levels

  group_mean_mat <- score_df |>
    group_by(pathway, group) |>
    summarise(mean_z = mean(score_z, na.rm = TRUE), .groups = "drop") |>
    mutate(group = factor(group, levels = group_levels)) |>
    arrange(match(pathway, pathway_levels), group) |>
    select(pathway, group, mean_z) |>
    pivot_wider(names_from = group, values_from = mean_z) |>
    as.data.frame()
  rownames(group_mean_mat) <- group_mean_mat$pathway
  group_mean_mat <- as.matrix(group_mean_mat[, group_levels, drop = FALSE])
  colnames(group_mean_mat) <- c("WT host", "FAP-deficient host")
  rownames(group_mean_mat) <- pathway_levels

  group_ann <- data.frame(Group = factor(group_levels, levels = group_levels))
  rownames(group_ann) <- colnames(group_mean_mat)

  plot_boxplots(
    score_df,
    stats_df,
    dataset_name,
    file.path(output_dir, paste0("01_", tolower(dataset_name), "_progeny_boxplots.png")),
    file.path(output_dir, paste0("01_", tolower(dataset_name), "_progeny_boxplots.pdf"))
  )
  plot_heatmap_dual(
    sample_heatmap_mat,
    sample_ann,
    paste0(dataset_name, ": sample-level PROGENy heatmap"),
    file.path(output_dir, paste0("02_", tolower(dataset_name), "_progeny_heatmap_samples.png")),
    file.path(output_dir, paste0("02_", tolower(dataset_name), "_progeny_heatmap_samples.pdf"))
  )
  plot_heatmap_dual(
    group_mean_mat,
    group_ann,
    paste0(dataset_name, ": group mean PROGENy heatmap"),
    file.path(output_dir, paste0("03_", tolower(dataset_name), "_progeny_heatmap_group_means.png")),
    file.path(output_dir, paste0("03_", tolower(dataset_name), "_progeny_heatmap_group_means.pdf")),
    width = 5.4,
    height = 5.8
  )

  write_table_tsv(score_df, file.path(output_dir, paste0(tolower(dataset_name), "_progeny_scores_per_sample.tsv")))
  write_table_tsv(stats_df, file.path(output_dir, paste0(tolower(dataset_name), "_progeny_pathway_stats.tsv")))
  write_table_tsv(coverage_df, file.path(output_dir, paste0(tolower(dataset_name), "_progeny_model_coverage.tsv")))

  list(
    score_df = score_df,
    stats_df = stats_df,
    coverage_df = coverage_df,
    group_mean_mat = group_mean_mat
  )
}

sample_info <- read_sample_info(file.path(root, "002_DIA_Summary", "Sample_Info.txt"))
rna_expr <- read_rna_expression(file.path(root, "001_mRNA_Summary", "1.GeneExpression", "1_genes_fpkm_expression.txt"))
rna_samples <- grep("^FPKM\\.", names(rna_expr), value = TRUE)
prot_expr <- read_protein_expression(file.path(root, "002_DIA_Summary", "01.RawData", "MatchResult", "protein.groups.intensity.txt"), sample_info)
prot_samples <- setdiff(names(prot_expr), "gene_name")

rna_res <- run_progeny_modality(rna_expr, "gene_name", rna_samples, "RNA", rna_dir)
prot_res <- run_progeny_modality(prot_expr, "gene_name", prot_samples, "Proteomics", proteomics_dir)

integrated_group_mean_df <- bind_rows(
  as.data.frame(rna_res$group_mean_mat) |>
    tibble::rownames_to_column("pathway") |>
    mutate(dataset = "RNA"),
  as.data.frame(prot_res$group_mean_mat) |>
    tibble::rownames_to_column("pathway") |>
    mutate(dataset = "Proteomics")
) |>
  pivot_longer(all_of(group_levels), names_to = "group", values_to = "mean_z") |>
  mutate(column_id = paste(dataset, group, sep = " | ")) |>
  select(pathway, column_id, mean_z) |>
  pivot_wider(names_from = column_id, values_from = mean_z) |>
  as.data.frame()

rownames(integrated_group_mean_df) <- integrated_group_mean_df$pathway
integrated_group_mean_mat <- as.matrix(integrated_group_mean_df[, c(
  "RNA | WT host",
  "RNA | FAP-deficient host",
  "Proteomics | WT host",
  "Proteomics | FAP-deficient host"
), drop = FALSE])
rownames(integrated_group_mean_mat) <- pathway_levels

integrated_ann <- data.frame(
  Group = factor(c("WT host", "FAP-deficient host", "WT host", "FAP-deficient host"), levels = group_levels)
)
rownames(integrated_ann) <- colnames(integrated_group_mean_mat)

plot_heatmap_dual(
  integrated_group_mean_mat,
  integrated_ann,
  "Integrated PROGENy group mean comparison",
  file.path(integrated_dir, "04_progeny_group_mean_comparison.png"),
  file.path(integrated_dir, "04_progeny_group_mean_comparison.pdf")
)

plot_heatmap_dual(
  integrated_group_mean_mat,
  integrated_ann,
  "Integrated PROGENy group mean comparison (row-clustered)",
  file.path(integrated_dir, "06_progeny_group_mean_comparison_row_clustered.png"),
  file.path(integrated_dir, "06_progeny_group_mean_comparison_row_clustered.pdf"),
  cluster_rows = TRUE,
  cluster_cols = FALSE
)

delta_df <- inner_join(
  rna_res$stats_df |>
    select(pathway, rna_delta_z = delta_fap_minus_wt_z, rna_q = q_value),
  prot_res$stats_df |>
    select(pathway, proteomics_delta_z = delta_fap_minus_wt_z, proteomics_q = q_value),
  by = "pathway"
) |>
  arrange(match(pathway, pathway_levels))

plot_delta_concordance(
  delta_df,
  file.path(integrated_dir, "05_progeny_delta_concordance.png"),
  file.path(integrated_dir, "05_progeny_delta_concordance.pdf")
)

write_table_tsv(
  as.data.frame(integrated_group_mean_mat) |>
    tibble::rownames_to_column("pathway"),
  file.path(integrated_dir, "progeny_group_mean_comparison.tsv")
)
write_table_tsv(delta_df, file.path(integrated_dir, "progeny_delta_concordance.tsv"))

summary_lines <- c(
  "WT host vs FAP-deficient host xenograft PROGENy summary",
  paste("Output directory:", outdir),
  paste("RNA outputs:", rna_dir),
  paste("Proteomics outputs:", proteomics_dir),
  paste("Integrated outputs:", integrated_dir),
  "",
  "Method:",
  "- Official PROGENy with Mouse model, scale=TRUE, top=100, z_scores=FALSE.",
  "- Input matrices were log2-transformed as log2(expression + 1) before PROGENy scoring.",
  "- Boxplots and heatmaps use modality-internal z-scores of the resulting pathway scores.",
  "",
  paste("Top RNA pathways higher in FAP-deficient host:", format_top_pathways(rna_res$stats_df, "up")),
  paste("Top RNA pathways higher in WT host:", format_top_pathways(rna_res$stats_df, "down")),
  paste("Top proteomics pathways higher in FAP-deficient host:", format_top_pathways(prot_res$stats_df, "up")),
  paste("Top proteomics pathways higher in WT host:", format_top_pathways(prot_res$stats_df, "down"))
)

writeLines(summary_lines, con = file.path(outdir, "summary_notes.txt"))

message("PROGENy analysis complete. Outputs written to: ", outdir)
