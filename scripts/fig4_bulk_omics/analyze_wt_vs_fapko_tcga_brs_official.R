#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Figure 4D (BRAF-RAS score row) - official TCGA-THCA BRS panel scoring.

script_path <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
})

theme_set(
  theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      axis.text = element_text(colour = "black")
    )
)

tcga_dir <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_public", "tcga_thca")
tcga_clin_path <- file.path(tcga_dir, "data_clinil_with_TCGA_2025.4.16.csv")
tcga_expr_path <- file.path(tcga_dir, "TCGA_mRNA482.csv")
panel_path <- file.path(root, "reference", "tcga_brs_71_gene_panel.tsv")

outdir <- file.path(root, "outputs", "wt_vs_fapko_bulk_omics_20260407", "tcga_brs_official")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
reference_dir <- file.path(outdir, "reference")
bulk_dir <- file.path(outdir, "bulk_rna")
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bulk_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("WT host", "FAP-deficient host")
tcga_plot_levels <- c("Braf-like", "Ras-like", "Unassigned")

clean_sample_name <- function(x) {
  sub("^FPKM\\.", "", x)
}

label_bulk_group <- function(x) {
  ifelse(grepl("^FPKM\\.FAP_BPC", x), "FAP-deficient host", "WT host")
}

read_tcga_panel <- function(path) {
  read.delim(path, check.names = FALSE) |>
    pull(symbol) |>
    unique() |>
    toupper()
}

read_tcga_expression <- function(path) {
  df <- read.csv(path, check.names = FALSE)
  sample_cols <- setdiff(names(df), c("", "Unnamed: 0", "Hugo_Symbol", "Entrez_Gene_Id"))
  df2 <- df |>
    select(Hugo_Symbol, all_of(sample_cols)) |>
    filter(!is.na(Hugo_Symbol), Hugo_Symbol != "") |>
    mutate(Hugo_Symbol = toupper(Hugo_Symbol))
  dup_idx <- duplicated(df2$Hugo_Symbol)
  df2 <- df2[!dup_idx, , drop = FALSE]
  mat <- as.matrix(df2[, sample_cols, drop = FALSE])
  rownames(mat) <- df2$Hugo_Symbol
  storage.mode(mat) <- "numeric"
  mat
}

read_bulk_expression <- function(path) {
  df <- read.delim(path, check.names = FALSE)
  sample_cols <- grep("^FPKM\\.", names(df), value = TRUE)
  df2 <- df |>
    select(gene_name, all_of(sample_cols)) |>
    filter(!is.na(gene_name), gene_name != "") |>
    mutate(
      gene_name = toupper(gene_name),
      across(all_of(sample_cols), as.numeric)
    )
  summarised <- df2 |>
    group_by(gene_name) |>
    summarise(across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  mat <- as.matrix(summarised[, sample_cols, drop = FALSE])
  rownames(mat) <- summarised$gene_name
  storage.mode(mat) <- "numeric"
  mat
}

score_against_centroids <- function(expr_mat, mean_vec, sd_vec, b_centroid, r_centroid, ref_range) {
  common_genes <- intersect(names(mean_vec), rownames(expr_mat))
  if (length(common_genes) < 10) {
    stop("Too few genes available for BRS scoring: ", length(common_genes))
  }

  log_mat <- log2(expr_mat[common_genes, , drop = FALSE] + 1)
  z_mat <- sweep(sweep(log_mat, 1, mean_vec[common_genes], "-"), 1, sd_vec[common_genes], "/")
  raw_score <- apply(
    z_mat,
    2,
    function(x) {
      mean((x - b_centroid[common_genes])^2) - mean((x - r_centroid[common_genes])^2)
    }
  )
  braf_distance <- apply(
    z_mat,
    2,
    function(x) sqrt(mean((x - b_centroid[common_genes])^2))
  )
  ras_distance <- apply(
    z_mat,
    2,
    function(x) sqrt(mean((x - r_centroid[common_genes])^2))
  )
  scaled_score <- (raw_score - ref_range[1]) / diff(ref_range) * 2 - 1

  list(
    genes_used = common_genes,
    z_mat = z_mat,
    raw_score = raw_score,
    scaled_score = scaled_score,
    braf_distance = braf_distance,
    ras_distance = ras_distance,
    braf_proximity = -braf_distance,
    ras_proximity = -ras_distance
  )
}

plot_reference_boxplot <- function(df, path_png, path_pdf) {
  plot_df <- df |>
    mutate(
      brs_label = ifelse(is.na(brs_label), "Unassigned", brs_label),
      brs_label = factor(brs_label, levels = tcga_plot_levels)
    )

  p <- ggplot(plot_df, aes(brs_label, brs_scaled, fill = brs_label)) +
    geom_boxplot(outlier.shape = NA, width = 0.58, alpha = 0.92, linewidth = 0.7, colour = "#2C2C2C") +
    geom_point(
      shape = 21,
      size = 1.8,
      stroke = 0.3,
      colour = "#2C2C2C",
      position = position_jitter(width = 0.12, height = 0)
    ) +
    geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.6, linetype = "dashed") +
    scale_fill_manual(values = c("Braf-like" = "#D47A6A", "Ras-like" = "#4C78A8", "Unassigned" = "#B8B8B8")) +
    labs(
      title = "TCGA reference BRS validation",
      subtitle = "Official 71-gene TCGA-BRS centroid framework; negative = BRAF-like, positive = RAS-like",
      x = NULL,
      y = "Scaled BRS"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(face = "bold", size = 9.2),
      axis.text.y = element_text(size = 8.4),
      axis.title.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.1)
    )

  ggsave(path_png, p, width = 6.6, height = 5.2, dpi = 320)
  ggsave(path_pdf, p, width = 6.6, height = 5.2)
}

plot_bulk_boxplot <- function(df, stats_df, path_png, path_pdf) {
  y_min <- min(df$brs_raw, na.rm = TRUE)
  y_max <- max(df$brs_raw, na.rm = TRUE)
  span <- max(y_max - y_min, 0.2)
  bracket_y <- y_max + span * 0.12
  bracket_top <- bracket_y + span * 0.08
  label_y <- bracket_top + span * 0.05

  sig_label <- if (is.na(stats_df$wilcox_p)) {
    "NA"
  } else if (stats_df$wilcox_p < 0.01) {
    "**"
  } else if (stats_df$wilcox_p < 0.05) {
    "*"
  } else if (stats_df$wilcox_p < 0.10) {
    "."
  } else {
    "ns"
  }

  p <- ggplot(df, aes(group, brs_raw, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.56, alpha = 0.92, linewidth = 0.7, colour = "#2C2C2C") +
    geom_point(
      shape = 21,
      size = 2.1,
      stroke = 0.35,
      colour = "#2C2C2C",
      position = position_jitter(width = 0.10, height = 0)
    ) +
    annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_top, linewidth = 0.6) +
    annotate("segment", x = 1, xend = 2, y = bracket_top, yend = bracket_top, linewidth = 0.6) +
    annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_top, linewidth = 0.6) +
    annotate("text", x = 1.5, y = label_y, label = sig_label, size = 3.2, fontface = "bold") +
    scale_fill_manual(values = c("WT host" = "#D47A6A", "FAP-deficient host" = "#46B58A")) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.24))) +
    labs(
      title = "WT host vs FAP-deficient host official TCGA-BRS",
      subtitle = "RNA-only projection onto TCGA THCA BRAF and RAS centroids; lower = more BRAF-like, higher = more RAS-like",
      x = NULL,
      y = "Raw centroid-distance BRS"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1, face = "bold", size = 8.8),
      axis.text.y = element_text(size = 8.2),
      axis.title.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.0)
    )

  ggsave(path_png, p, width = 5.8, height = 5.3, dpi = 320)
  ggsave(path_pdf, p, width = 5.8, height = 5.3)
}

compute_component_stats <- function(df) {
  tibble(
    metric = c("BRAF-centroid proximity", "RAS-centroid proximity"),
    wt_mean = c(
      mean(df$braf_proximity[df$group == "WT host"]),
      mean(df$ras_proximity[df$group == "WT host"])
    ),
    fap_mean = c(
      mean(df$braf_proximity[df$group == "FAP-deficient host"]),
      mean(df$ras_proximity[df$group == "FAP-deficient host"])
    ),
    delta_fap_minus_wt = fap_mean - wt_mean,
    wilcox_p = c(
      wilcox.test(braf_proximity ~ group, data = df, exact = FALSE)$p.value,
      wilcox.test(ras_proximity ~ group, data = df, exact = FALSE)$p.value
    ),
    welch_p = c(
      t.test(braf_proximity ~ group, data = df)$p.value,
      t.test(ras_proximity ~ group, data = df)$p.value
    )
  )
}

plot_component_boxplots <- function(df, stats_df, title_text, subtitle_text, path_png, path_pdf) {
  long_df <- df |>
    select(sample, group, braf_proximity, ras_proximity) |>
    pivot_longer(
      cols = c(braf_proximity, ras_proximity),
      names_to = "metric",
      values_to = "score"
    ) |>
    mutate(
      metric = recode(
        metric,
        braf_proximity = "BRAF-centroid proximity",
        ras_proximity = "RAS-centroid proximity"
      ),
      metric = factor(metric, levels = c("BRAF-centroid proximity", "RAS-centroid proximity"))
    )

  ann <- stats_df |>
    mutate(
      sig_label = case_when(
        is.na(wilcox_p) ~ "NA",
        wilcox_p < 0.01 ~ "**",
        wilcox_p < 0.05 ~ "*",
        wilcox_p < 0.10 ~ ".",
        TRUE ~ "ns"
      )
    )

  bounds <- long_df |>
    group_by(metric) |>
    summarise(
      y_min = min(score, na.rm = TRUE),
      y_max = max(score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      span = pmax(y_max - y_min, 0.2),
      bracket_y = y_max + span * 0.12,
      bracket_top = bracket_y + span * 0.08,
      label_y = bracket_top + span * 0.05
    )

  ann <- left_join(ann, bounds, by = "metric")

  p <- ggplot(long_df, aes(group, score, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.56, alpha = 0.92, linewidth = 0.7, colour = "#2C2C2C") +
    geom_point(
      shape = 21,
      size = 2.0,
      stroke = 0.35,
      colour = "#2C2C2C",
      position = position_jitter(width = 0.10, height = 0)
    ) +
    geom_segment(
      data = ann,
      inherit.aes = FALSE,
      aes(x = 1, xend = 1, y = bracket_y, yend = bracket_top),
      linewidth = 0.6
    ) +
    geom_segment(
      data = ann,
      inherit.aes = FALSE,
      aes(x = 1, xend = 2, y = bracket_top, yend = bracket_top),
      linewidth = 0.6
    ) +
    geom_segment(
      data = ann,
      inherit.aes = FALSE,
      aes(x = 2, xend = 2, y = bracket_y, yend = bracket_top),
      linewidth = 0.6
    ) +
    geom_text(
      data = ann,
      inherit.aes = FALSE,
      aes(x = 1.5, y = label_y, label = sig_label),
      size = 3.1,
      fontface = "bold"
    ) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("WT host" = "#D47A6A", "FAP-deficient host" = "#46B58A")) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.24))) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = NULL,
      y = "Centroid proximity\n(higher = closer)"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1, face = "bold", size = 8.8),
      axis.text.y = element_text(size = 8.2),
      axis.title.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.0),
      strip.text = element_text(face = "bold", size = 8.8),
      strip.background = element_rect(fill = "white", colour = "#2C2C2C", linewidth = 0.8)
    )

  ggsave(path_png, p, width = 8.8, height = 5.4, dpi = 320)
  ggsave(path_pdf, p, width = 8.8, height = 5.4)
}

tcga_panel <- read_tcga_panel(panel_path)
tcga_clin <- read.csv(tcga_clin_path, check.names = FALSE)
tcga_mat <- read_tcga_expression(tcga_expr_path)

common_tcga_samples <- intersect(colnames(tcga_mat), tcga_clin$SAMPLE_ID)
tcga_mat <- tcga_mat[, common_tcga_samples, drop = FALSE]
tcga_clin <- tcga_clin[match(common_tcga_samples, tcga_clin$SAMPLE_ID), , drop = FALSE]

train_hist <- c("Classical", "Follicular", "Tall Cell")
is_braf_train <- grepl("BRAF:p.V600E", tcga_clin$MUT_DRIVER_PROTEIN_CHANGE, fixed = TRUE) &
  tcga_clin$HISTOLOGICAL_TYPE %in% train_hist
is_ras_train <- grepl("HRAS:|KRAS:|NRAS:", tcga_clin$MUT_DRIVER_PROTEIN_CHANGE) &
  tcga_clin$HISTOLOGICAL_TYPE %in% train_hist

tcga_present <- intersect(tcga_panel, rownames(tcga_mat))
tcga_ref_expr <- log2(tcga_mat[tcga_present, is_braf_train | is_ras_train, drop = FALSE] + 1)
tcga_ref_mean <- rowMeans(tcga_ref_expr, na.rm = TRUE)
tcga_ref_sd <- apply(tcga_ref_expr, 1, sd, na.rm = TRUE)
tcga_ref_sd[!is.finite(tcga_ref_sd) | tcga_ref_sd == 0] <- 1

tcga_score_res <- score_against_centroids(
  expr_mat = tcga_mat,
  mean_vec = tcga_ref_mean,
  sd_vec = tcga_ref_sd,
  b_centroid = apply(
    sweep(sweep(log2(tcga_mat[tcga_present, is_braf_train, drop = FALSE] + 1), 1, tcga_ref_mean, "-"), 1, tcga_ref_sd, "/"),
    1,
    median,
    na.rm = TRUE
  ),
  r_centroid = apply(
    sweep(sweep(log2(tcga_mat[tcga_present, is_ras_train, drop = FALSE] + 1), 1, tcga_ref_mean, "-"), 1, tcga_ref_sd, "/"),
    1,
    median,
    na.rm = TRUE
  ),
  ref_range = c(0, 1)
)

# Recalculate with the true reference range now that raw scores are available.
tcga_b_centroid <- apply(tcga_score_res$z_mat[, is_braf_train, drop = FALSE], 1, median, na.rm = TRUE)
tcga_r_centroid <- apply(tcga_score_res$z_mat[, is_ras_train, drop = FALSE], 1, median, na.rm = TRUE)
tcga_raw <- apply(
  tcga_score_res$z_mat,
  2,
  function(x) mean((x - tcga_b_centroid)^2) - mean((x - tcga_r_centroid)^2)
)
tcga_ref_range <- range(tcga_raw[is_braf_train | is_ras_train], na.rm = TRUE)
tcga_scaled <- (tcga_raw - tcga_ref_range[1]) / diff(tcga_ref_range) * 2 - 1
tcga_braf_distance <- apply(
  tcga_score_res$z_mat,
  2,
  function(x) sqrt(mean((x - tcga_b_centroid)^2))
)
tcga_ras_distance <- apply(
  tcga_score_res$z_mat,
  2,
  function(x) sqrt(mean((x - tcga_r_centroid)^2))
)

tcga_reference_scores <- tibble(
  SAMPLE_ID = common_tcga_samples,
  brs_label = ifelse(is.na(tcga_clin$BRAFV600E_RAS), "Unassigned", tcga_clin$BRAFV600E_RAS),
  histological_type = tcga_clin$HISTOLOGICAL_TYPE,
  mut_driver = tcga_clin$MUT_DRIVER_PROTEIN_CHANGE,
  train_braf = is_braf_train,
  train_ras = is_ras_train,
  brs_raw = tcga_raw,
  brs_scaled = tcga_scaled,
  braf_distance = tcga_braf_distance,
  ras_distance = tcga_ras_distance,
  braf_proximity = -tcga_braf_distance,
  ras_proximity = -tcga_ras_distance,
  erk_score = tcga_clin$ERK_SCORE,
  tds = tcga_clin$DIFFERENTIATION_SCORE
)

reference_gene_status <- tibble(
  published_symbol = tcga_panel,
  present_in_tcga = tcga_panel %in% rownames(tcga_mat),
  present_in_bulk = NA
)

bulk_mat <- read_bulk_expression(file.path(root, "001_mRNA_Summary", "1.GeneExpression", "1_genes_fpkm_expression.txt"))
bulk_score_res <- score_against_centroids(
  expr_mat = bulk_mat,
  mean_vec = tcga_ref_mean,
  sd_vec = tcga_ref_sd,
  b_centroid = tcga_b_centroid,
  r_centroid = tcga_r_centroid,
  ref_range = tcga_ref_range
)

reference_gene_status$present_in_bulk <- reference_gene_status$published_symbol %in% rownames(bulk_mat)

bulk_scores <- tibble(
  sample_raw = names(bulk_score_res$raw_score),
  sample = clean_sample_name(names(bulk_score_res$raw_score)),
  group = factor(label_bulk_group(names(bulk_score_res$raw_score)), levels = group_levels),
  brs_raw = bulk_score_res$raw_score,
  brs_scaled_reference = bulk_score_res$scaled_score,
  braf_distance = bulk_score_res$braf_distance,
  ras_distance = bulk_score_res$ras_distance,
  braf_proximity = bulk_score_res$braf_proximity,
  ras_proximity = bulk_score_res$ras_proximity
)

bulk_stats <- tibble(
  genes_requested = length(tcga_panel),
  genes_present_in_reference = length(tcga_present),
  genes_present_in_bulk = length(bulk_score_res$genes_used),
  wt_mean_raw = mean(bulk_scores$brs_raw[bulk_scores$group == "WT host"]),
  fap_mean_raw = mean(bulk_scores$brs_raw[bulk_scores$group == "FAP-deficient host"]),
  delta_fap_minus_wt_raw = fap_mean_raw - wt_mean_raw,
  wt_mean_scaled_reference = mean(bulk_scores$brs_scaled_reference[bulk_scores$group == "WT host"]),
  fap_mean_scaled_reference = mean(bulk_scores$brs_scaled_reference[bulk_scores$group == "FAP-deficient host"]),
  delta_fap_minus_wt_scaled_reference = fap_mean_scaled_reference - wt_mean_scaled_reference,
  wilcox_p = wilcox.test(brs_raw ~ group, data = bulk_scores, exact = FALSE)$p.value,
  welch_p = t.test(brs_raw ~ group, data = bulk_scores)$p.value
)
bulk_component_stats <- compute_component_stats(bulk_scores)

plot_reference_boxplot(
  tcga_reference_scores,
  file.path(reference_dir, "01_tcga_reference_brs_boxplot.png"),
  file.path(reference_dir, "01_tcga_reference_brs_boxplot.pdf")
)
plot_bulk_boxplot(
  bulk_scores,
  bulk_stats,
  file.path(bulk_dir, "02_bulk_official_tcga_brs_boxplot.png"),
  file.path(bulk_dir, "02_bulk_official_tcga_brs_boxplot.pdf")
)
plot_component_boxplots(
  bulk_scores,
  bulk_component_stats,
  "WT host vs FAP-deficient host official centroid components",
  "Derived from the same official TCGA 71-gene centroid framework as BRS",
  file.path(bulk_dir, "03_bulk_official_braf_ras_component_boxplots.png"),
  file.path(bulk_dir, "03_bulk_official_braf_ras_component_boxplots.pdf")
)

write.table(reference_gene_status, file.path(reference_dir, "tcga_brs_gene_panel_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tcga_reference_scores, file.path(reference_dir, "tcga_reference_brs_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bulk_scores, file.path(bulk_dir, "bulk_official_tcga_brs_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bulk_stats, file.path(bulk_dir, "bulk_official_tcga_brs_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bulk_component_stats, file.path(bulk_dir, "bulk_official_braf_ras_component_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_lines <- c(
  "Official TCGA-BRS projection summary",
  paste("Output directory:", outdir),
  "",
  "Reference framework:",
  "- Panel: published 71-gene TCGA BRAF-RAS score panel.",
  "- Training labels: TCGA THCA samples with exact BRAF:p.V600E versus HRAS/KRAS/NRAS driver calls.",
  "- Histology restriction for centroid training: Classical, Follicular, and Tall Cell thyroid carcinoma.",
  "- Transformation: log2(expression + 1), then gene-wise z-scoring using the TCGA reference cohort.",
  "- Centroids: gene-wise medians of the BRAF V600E and RAS-mutant reference groups.",
  "- Raw BRS = mean squared distance to BRAF centroid minus mean squared distance to RAS centroid.",
  "- Negative raw BRS indicates BRAF-like state; positive raw BRS indicates RAS-like state.",
  "- Scaled BRS was linearly rescaled using the TCGA reference cohort range; external samples may fall outside [-1, 1].",
  "",
  paste0("Reference genes evaluable in TCGA: ", length(tcga_present), "/", length(tcga_panel)),
  paste0("Reference genes evaluable in bulk RNA: ", length(bulk_score_res$genes_used), "/", length(tcga_panel)),
  paste0("Bulk WT mean raw BRS: ", sprintf("%.3f", bulk_stats$wt_mean_raw)),
  paste0("Bulk FAP-deficient mean raw BRS: ", sprintf("%.3f", bulk_stats$fap_mean_raw)),
  paste0("Bulk delta (FAP - WT) raw BRS: ", sprintf("%.3f", bulk_stats$delta_fap_minus_wt_raw)),
  paste0("Bulk Wilcoxon p: ", formatC(bulk_stats$wilcox_p, format = "e", digits = 3)),
  paste0("Bulk Welch p: ", formatC(bulk_stats$welch_p, format = "e", digits = 3)),
  paste0("Bulk BRAF-centroid proximity delta (FAP - WT): ", sprintf("%.3f", bulk_component_stats$delta_fap_minus_wt[bulk_component_stats$metric == "BRAF-centroid proximity"])),
  paste0("Bulk RAS-centroid proximity delta (FAP - WT): ", sprintf("%.3f", bulk_component_stats$delta_fap_minus_wt[bulk_component_stats$metric == "RAS-centroid proximity"])),
  "",
  "Interpretation note:",
  "- This is a direct RNA-only projection onto TCGA thyroid reference centroids.",
  "- BRAF-like and RAS-like components are reported here as centroid proximities derived from the same official framework, not as separate published TCGA standalone scores.",
  "- Because the xenograft bulk RNA data are cross-study and mixed-compartment, relative group differences are more interpretable than absolute score magnitude."
)

writeLines(summary_lines, file.path(outdir, "summary_notes.txt"))

message("Official TCGA-BRS analysis completed. Results saved to: ", outdir)
