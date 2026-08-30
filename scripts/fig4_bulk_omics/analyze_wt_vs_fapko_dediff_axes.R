#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Figure 4D (lineage/dediff rows) - dedifferentiation axis scoring.

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
})

theme_set(
  theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      axis.text = element_text(colour = "black")
    )
)

outdir <- file.path(root, "outputs", "wt_vs_fapko_bulk_omics_20260407", "dediff_axes")
rna_dir <- file.path(outdir, "rna")
proteomics_dir <- file.path(outdir, "proteomics")
integrated_dir <- file.path(outdir, "integrated")
dir.create(rna_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proteomics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(integrated_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("WT host", "FAP-deficient host")

signature_defs <- list(
  "TDS / lineage identity" = c("DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1", "PAX8", "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR"),
  "ERK score" = c("DUSP4", "DUSP5", "DUSP6", "ETV1", "ETV4", "ETV5", "SPRY2", "SPRY4", "FOS", "FOSL1", "JUNB", "EGR1"),
  "Stemness / plasticity" = c("SOX2", "KLF4", "MYC", "NANOG", "POU5F1", "LIN28A", "PROM1", "ALDH1A1", "KIT", "ITGA6", "BMI1", "CD44"),
  "YAP/TAZ" = c("YAP1", "WWTR1", "CTGF", "CYR61", "ANKRD1", "AMOTL2", "AXL", "FGF1", "AJUBA", "BIRC5", "TAGLN", "SERPINE1"),
  "ECM-integrin-FAK" = c("FN1", "COL1A1", "COL1A2", "POSTN", "SPARC", "ITGA5", "ITGB1", "ITGAV", "PXN", "PTK2", "SRC", "TLN1", "VCL", "ZYX"),
  "HGF-MET" = c("HGF", "MET", "GAB1", "GRB2", "SHC1", "PTPN11", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "PIK3CA", "PIK3CB", "PIK3R1", "AKT1", "AKT2", "STAT3", "DUSP6", "ETV4", "ETV5"),
  "BRAF-like program" = c("DUSP4", "DUSP5", "DUSP6", "ETV1", "ETV4", "ETV5", "SPRY2", "SPRY4", "FOSL1", "CCND1", "MYC", "ERBB3"),
  "RAS-like lineage program" = c("FOXE1", "PAX8", "NKX2-1", "TG", "TPO", "TSHR", "SLC5A5", "DIO1", "DIO2", "GLIS3", "THRB"),
  "Notch signaling" = c("HES1", "HEY1", "HEY2", "NRARP", "DTX1", "NOTCH1", "NOTCH3", "JAG1"),
  "Hedgehog signaling" = c("GLI1", "GLI2", "PTCH1", "PTCH2", "HHIP", "SMO", "SUFU")
)

plot_signature_order <- c(
  "TDS / lineage identity",
  "Stemness / plasticity",
  "BRAF-RAS score",
  "ERK score",
  "YAP/TAZ",
  "ECM-integrin-FAK",
  "HGF-MET",
  "Notch signaling",
  "Hedgehog signaling"
)

signature_categories <- c(
  "TDS / lineage identity" = "Lineage / dediff",
  "Stemness / plasticity" = "Lineage / dediff",
  "BRAF-RAS score" = "Lineage / dediff",
  "ERK score" = "Signaling output",
  "YAP/TAZ" = "Mechanotransduction / niche",
  "ECM-integrin-FAK" = "Mechanotransduction / niche",
  "HGF-MET" = "RTK / stromal relay",
  "Notch signaling" = "Developmental signaling",
  "Hedgehog signaling" = "Developmental signaling"
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
    filter(!is.na(gene_name), gene_name != "") |>
    mutate(gene_name = toupper(gene_name))
  collapse_by_gene(df2, "gene_name", sample_cols)
}

read_protein_expression <- function(path, sample_info) {
  df <- read.delim(path, check.names = FALSE)
  sample_cols <- grep("^0208JLYAstral_DIA_zyl_P4_", names(df), value = TRUE)
  raw_ids <- sub("^0208JLYAstral_DIA_zyl_P4_", "", sample_cols)
  mapped_names <- sample_info$SampleID[match(raw_ids, sample_info$RawID)]
  rename_map <- setNames(sample_cols, mapped_names)
  df2 <- df |>
    mutate(gene_name = toupper(sub(";.*$", "", Genes))) |>
    select(gene_name, all_of(sample_cols)) |>
    rename(!!!rename_map) |>
    mutate(across(-gene_name, as.numeric)) |>
    filter(!is.na(gene_name), gene_name != "")
  collapse_by_gene(df2, "gene_name", mapped_names)
}

prepare_matrix <- function(expr_df, gene_col, sample_cols) {
  mat <- as.matrix(expr_df[, sample_cols, drop = FALSE])
  rownames(mat) <- toupper(expr_df[[gene_col]])
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  storage.mode(mat) <- "numeric"
  log2(mat + 1)
}

safe_scale_rows <- function(mat) {
  out <- t(scale(t(mat)))
  out[!is.finite(out)] <- 0
  out
}

compute_signature_matrix <- function(expr_mat, signatures, min_genes = 3) {
  z_mat <- safe_scale_rows(expr_mat)
  score_rows <- list()
  coverage_rows <- list()

  for (sig_name in names(signatures)) {
    genes <- unique(toupper(signatures[[sig_name]]))
    present <- intersect(genes, rownames(z_mat))
    missing <- setdiff(genes, rownames(z_mat))
    coverage_rows[[sig_name]] <- tibble(
      signature = sig_name,
      genes_requested = length(genes),
      genes_present = length(present),
      genes_missing = length(missing),
      coverage_pct = 100 * length(present) / length(genes),
      present_genes = paste(present, collapse = ","),
      missing_genes = paste(missing, collapse = ",")
    )
    if (length(present) < min_genes) {
      score_rows[[sig_name]] <- rep(NA_real_, ncol(z_mat))
    } else {
      score_rows[[sig_name]] <- colMeans(z_mat[present, , drop = FALSE], na.rm = TRUE)
    }
  }

  score_mat <- do.call(rbind, score_rows)
  rownames(score_mat) <- names(signatures)
  colnames(score_mat) <- colnames(z_mat)

  list(score_mat = score_mat, coverage_df = bind_rows(coverage_rows))
}

add_braf_ras_score <- function(score_mat, coverage_df, min_component_genes = 3) {
  braf_ok <- coverage_df$genes_present[match("BRAF-like program", coverage_df$signature)] >= min_component_genes
  ras_ok <- coverage_df$genes_present[match("RAS-like lineage program", coverage_df$signature)] >= min_component_genes
  if (isTRUE(braf_ok) && isTRUE(ras_ok)) {
    composite <- score_mat["BRAF-like program", ] - score_mat["RAS-like lineage program", ]
  } else {
    composite <- rep(NA_real_, ncol(score_mat))
  }
  score_mat <- rbind(score_mat, `BRAF-RAS score` = composite)

  coverage_df <- bind_rows(
    coverage_df,
    tibble(
      signature = "BRAF-RAS score",
      genes_requested = NA_integer_,
      genes_present = NA_integer_,
      genes_missing = NA_integer_,
      coverage_pct = NA_real_,
      present_genes = paste0(
        "BRAF-like present=",
        coverage_df$genes_present[match("BRAF-like program", coverage_df$signature)],
        "; RAS-like present=",
        coverage_df$genes_present[match("RAS-like lineage program", coverage_df$signature)]
      ),
      missing_genes = "Composite score = BRAF-like minus RAS-like"
    )
  )

  list(score_mat = score_mat, coverage_df = coverage_df)
}

make_long_scores <- function(score_mat, dataset_name, coverage_df) {
  z_scores <- t(scale(t(score_mat)))
  z_scores[!is.finite(z_scores)] <- 0

  raw_long <- as.data.frame(score_mat) |>
    tibble::rownames_to_column("signature") |>
    pivot_longer(-signature, names_to = "sample_raw", values_to = "score_raw")

  z_long <- as.data.frame(z_scores) |>
    tibble::rownames_to_column("signature") |>
    pivot_longer(-signature, names_to = "sample_raw", values_to = "score_z")

  left_join(raw_long, z_long, by = c("signature", "sample_raw")) |>
    mutate(
      sample = clean_sample_name(sample_raw),
      group = factor(label_group(sample_raw), levels = group_levels),
      dataset = dataset_name
    ) |>
    left_join(coverage_df |> select(signature, genes_requested, genes_present, coverage_pct, present_genes), by = "signature")
}

compute_stats <- function(score_df) {
  score_df |>
    group_by(signature) |>
    summarise(
      wt_mean_raw = mean(score_raw[group == "WT host"], na.rm = TRUE),
      fap_mean_raw = mean(score_raw[group == "FAP-deficient host"], na.rm = TRUE),
      delta_fap_minus_wt_raw = fap_mean_raw - wt_mean_raw,
      wt_mean_z = mean(score_z[group == "WT host"], na.rm = TRUE),
      fap_mean_z = mean(score_z[group == "FAP-deficient host"], na.rm = TRUE),
      delta_fap_minus_wt_z = fap_mean_z - wt_mean_z,
      p_value = tryCatch(wilcox.test(score_raw ~ group, exact = FALSE)$p.value, error = function(e) NA_real_),
      genes_requested = first(genes_requested),
      genes_present = first(genes_present),
      coverage_pct = first(coverage_pct),
      present_genes = first(present_genes),
      .groups = "drop"
    ) |>
    mutate(
      across(c(wt_mean_raw, fap_mean_raw, delta_fap_minus_wt_raw, wt_mean_z, fap_mean_z, delta_fap_minus_wt_z), ~ ifelse(is.nan(.x), NA_real_, .x)),
      q_value = p.adjust(p_value, method = "BH")
    )
}

plot_boxplots <- function(score_df, stats_df, dataset_name, path_png, path_pdf) {
  available <- score_df |>
    group_by(signature) |>
    summarise(any_value = any(is.finite(score_raw)), .groups = "drop") |>
    filter(any_value) |>
    pull(signature)

  plot_order <- intersect(plot_signature_order, available)
  if (length(plot_order) == 0) {
    return(invisible(NULL))
  }

  ann <- stats_df |>
    filter(signature %in% plot_order) |>
    mutate(
      sig_label = case_when(
        is.na(p_value) ~ "NA",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        p_value < 0.10 ~ ".",
        TRUE ~ "ns"
      )
    ) |>
    select(signature, sig_label)

  bounds <- score_df |>
    filter(signature %in% plot_order) |>
    group_by(signature) |>
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

  ann <- left_join(ann, bounds, by = "signature")

  plot_df <- score_df |>
    filter(signature %in% plot_order) |>
    mutate(signature = factor(signature, levels = plot_order))

  ann <- ann |>
    mutate(signature = factor(signature, levels = plot_order))

  plot_colors <- c("WT host" = "#D47A6A", "FAP-deficient host" = "#46B58A")

  p <- ggplot(plot_df, aes(x = group, y = score_z)) +
    geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.56, alpha = 0.92, linewidth = 0.7, colour = "#2C2C2C") +
    geom_point(
      aes(fill = group),
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
      size = 3.0,
      fontface = "bold"
    ) +
    facet_wrap(~ signature, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = plot_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.24))) +
    labs(
      title = paste0(dataset_name, " dedifferentiation-related axes"),
      subtitle = "Boxplots show modality-internal z-scores; . p < 0.10, * p < 0.05, ** p < 0.01",
      x = NULL,
      y = "Axis score z"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 36, hjust = 1, vjust = 1, face = "bold", size = 8.8),
      axis.text.y = element_text(size = 8.2),
      axis.title.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0, size = 12),
      plot.subtitle = element_text(size = 9.2),
      strip.text = element_text(face = "bold", size = 8.4, margin = margin(4, 5, 4, 5)),
      strip.background = element_rect(fill = "white", colour = "#2C2C2C", linewidth = 0.8),
      panel.spacing = unit(0.65, "lines"),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.6, colour = "#2C2C2C")
    )

  ggsave(path_png, p, width = 16.0, height = 7.6, dpi = 300)
  ggsave(path_pdf, p, width = 16.0, height = 7.6)
}

plot_overview_heatmap <- function(stats_df, path_png, path_pdf) {
  plot_df <- stats_df |>
    mutate(
      dataset = factor(dataset, levels = c("RNA", "Proteomics")),
      signature = factor(signature, levels = rev(plot_signature_order)),
      category = factor(signature_categories[as.character(signature)],
                        levels = rev(c("Lineage / dediff", "Signaling output", "Mechanotransduction / niche", "RTK / stromal relay", "Developmental signaling"))),
      sig_label = case_when(
        is.na(p_value) ~ "",
        p_value < 0.05 ~ "**",
        p_value < 0.10 ~ "*",
        TRUE ~ ""
      )
    )

  p <- ggplot(plot_df, aes(x = dataset, y = signature, fill = delta_fap_minus_wt_z)) +
    geom_tile(colour = "white", linewidth = 0.7) +
    geom_text(aes(label = sig_label), size = 3.2) +
    scale_fill_gradient2(
      low = "#5E81AC",
      mid = "#FBFBFA",
      high = "#D08770",
      midpoint = 0,
      na.value = "grey90",
      limits = c(-2.2, 2.2),
      oob = scales::squish
    ) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(
      title = "Dedifferentiation-related axis overview",
      subtitle = "Red = higher in FAP-deficient host; grey = insufficient coverage for scoring",
      x = NULL,
      y = NULL,
      fill = "delta z\n(FAP-WT)"
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey96", colour = "grey55", linewidth = 0.5),
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 7.4),
      axis.text.x = element_text(face = "bold", size = 8.8),
      axis.text.y = element_text(face = "bold", size = 8.2),
      legend.title = element_text(size = 8.4, face = "bold"),
      legend.text = element_text(size = 7.8),
      legend.key.height = unit(0.42, "in"),
      legend.key.width = unit(0.16, "in"),
      panel.spacing.y = unit(0.16, "lines")
    )

  ggsave(path_png, p, width = 5.8, height = 8.4, dpi = 300)
  ggsave(path_pdf, p, width = 5.8, height = 8.4)
}

plot_group_mean_heatmap <- function(score_df, path_png, path_pdf) {
  plot_df <- score_df |>
    group_by(dataset, signature, group) |>
    summarise(mean_z = mean(score_z, na.rm = TRUE), .groups = "drop") |>
    mutate(mean_z = ifelse(is.nan(mean_z), NA_real_, mean_z)) |>
    mutate(
      column_id = paste(dataset, group, sep = " | "),
      signature = factor(signature, levels = rev(plot_signature_order)),
      category = factor(signature_categories[as.character(signature)],
                        levels = rev(c("Lineage / dediff", "Signaling output", "Mechanotransduction / niche", "RTK / stromal relay", "Developmental signaling"))),
      column_id = factor(column_id, levels = c(
        "RNA | WT host",
        "RNA | FAP-deficient host",
        "Proteomics | WT host",
        "Proteomics | FAP-deficient host"
      ))
    )

  p <- ggplot(plot_df, aes(x = column_id, y = signature, fill = mean_z)) +
    geom_tile(colour = "white", linewidth = 0.7) +
    scale_fill_gradient2(
      low = "#2563EB",
      mid = "white",
      high = "#DC2626",
      midpoint = 0,
      na.value = "grey90",
      limits = c(-2.4, 2.4),
      oob = scales::squish
    ) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(
      title = "Dedifferentiation-related group mean comparison",
      x = NULL,
      y = NULL,
      fill = "mean z"
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey96", colour = "grey55", linewidth = 0.5),
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 7.4),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 8.3),
      axis.text.y = element_text(face = "bold", size = 8.2),
      legend.title = element_text(size = 8.4, face = "bold"),
      legend.text = element_text(size = 7.8),
      legend.key.height = unit(0.42, "in"),
      legend.key.width = unit(0.16, "in"),
      panel.spacing.y = unit(0.16, "lines")
    )

  ggsave(path_png, p, width = 6.2, height = 8.6, dpi = 300)
  ggsave(path_pdf, p, width = 6.2, height = 8.6)
}

plot_delta_concordance <- function(stats_df, path_png, path_pdf) {
  delta_df <- inner_join(
    stats_df |>
      filter(dataset == "RNA") |>
      select(signature, rna_delta_z = delta_fap_minus_wt_z, rna_q = q_value),
    stats_df |>
      filter(dataset == "Proteomics") |>
      select(signature, proteomics_delta_z = delta_fap_minus_wt_z, proteomics_q = q_value),
    by = "signature"
  ) |>
    filter(is.finite(rna_delta_z), is.finite(proteomics_delta_z))

  if (nrow(delta_df) < 2) {
    return(invisible(NULL))
  }

  rho <- cor(delta_df$rna_delta_z, delta_df$proteomics_delta_z, method = "spearman")
  label_df <- delta_df |>
    arrange(desc(abs(rna_delta_z) + abs(proteomics_delta_z))) |>
    slice_head(n = 6)

  p <- ggplot(delta_df, aes(rna_delta_z, proteomics_delta_z, label = signature)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(size = 3.0, colour = "#7C3AED", alpha = 0.85) +
    geom_text(data = label_df, size = 2.8, check_overlap = TRUE, nudge_y = 0.03) +
    annotate(
      "text",
      x = min(delta_df$rna_delta_z, na.rm = TRUE),
      y = max(delta_df$proteomics_delta_z, na.rm = TRUE),
      label = paste0("Spearman rho = ", sprintf("%.2f", rho)),
      hjust = 0,
      vjust = 1,
      size = 3.4
    ) +
    labs(
      title = "Dedifferentiation-axis concordance across RNA and proteomics",
      subtitle = "Only axes with evaluable scores in both modalities are shown",
      x = "RNA delta z",
      y = "Proteomics delta z"
    )

  ggsave(path_png, p, width = 6.8, height = 5.6, dpi = 300)
  ggsave(path_pdf, p, width = 6.8, height = 5.6)
}

format_top_axes <- function(stats_df, direction = c("up", "down"), n = 5) {
  direction <- match.arg(direction)
  ordered <- if (direction == "up") {
    stats_df |>
      arrange(desc(delta_fap_minus_wt_raw), q_value)
  } else {
    stats_df |>
      arrange(delta_fap_minus_wt_raw, q_value)
  }

  paste(
    ordered |>
      slice_head(n = n) |>
      transmute(label = paste0(signature, "(delta=", sprintf("%.2f", delta_fap_minus_wt_raw), ", q=", formatC(q_value, format = "e", digits = 1), ")")) |>
      pull(label),
    collapse = "; "
  )
}

sample_info <- read_sample_info(file.path(root, "002_DIA_Summary", "Sample_Info.txt"))
rna_expr <- read_rna_expression(file.path(root, "001_mRNA_Summary", "1.GeneExpression", "1_genes_fpkm_expression.txt"))
rna_samples <- grep("^FPKM\\.", names(rna_expr), value = TRUE)
prot_expr <- read_protein_expression(file.path(root, "002_DIA_Summary", "01.RawData", "MatchResult", "protein.groups.intensity.txt"), sample_info)
prot_samples <- setdiff(names(prot_expr), "gene_name")

rna_mat <- prepare_matrix(rna_expr, "gene_name", rna_samples)
prot_mat <- prepare_matrix(prot_expr, "gene_name", prot_samples)

rna_res <- compute_signature_matrix(rna_mat, signature_defs, min_genes = 3)
rna_res <- add_braf_ras_score(rna_res$score_mat, rna_res$coverage_df, min_component_genes = 3)
prot_res <- compute_signature_matrix(prot_mat, signature_defs, min_genes = 3)
prot_res <- add_braf_ras_score(prot_res$score_mat, prot_res$coverage_df, min_component_genes = 3)

rna_scores <- make_long_scores(rna_res$score_mat, "RNA", rna_res$coverage_df) |>
  filter(signature %in% plot_signature_order) |>
  mutate(signature = factor(signature, levels = plot_signature_order))

prot_scores <- make_long_scores(prot_res$score_mat, "Proteomics", prot_res$coverage_df) |>
  filter(signature %in% plot_signature_order) |>
  mutate(signature = factor(signature, levels = plot_signature_order))

all_scores <- bind_rows(rna_scores, prot_scores) |>
  mutate(
    signature = factor(as.character(signature), levels = plot_signature_order),
    dataset = factor(dataset, levels = c("RNA", "Proteomics"))
  )

rna_stats <- compute_stats(rna_scores) |>
  mutate(dataset = "RNA")
prot_stats <- compute_stats(prot_scores) |>
  mutate(dataset = "Proteomics")
all_stats <- bind_rows(rna_stats, prot_stats) |>
  mutate(dataset = factor(dataset, levels = c("RNA", "Proteomics")))

write.table(rna_scores, file.path(rna_dir, "rna_dediff_axes_scores_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(prot_scores, file.path(proteomics_dir, "proteomics_dediff_axes_scores_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rna_res$coverage_df, file.path(rna_dir, "rna_dediff_axes_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(prot_res$coverage_df, file.path(proteomics_dir, "proteomics_dediff_axes_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rna_stats, file.path(rna_dir, "rna_dediff_axes_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(prot_stats, file.path(proteomics_dir, "proteomics_dediff_axes_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_scores, file.path(integrated_dir, "dediff_axes_scores_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_stats, file.path(integrated_dir, "dediff_axes_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_boxplots(
  rna_scores,
  rna_stats,
  "RNA",
  file.path(rna_dir, "01_rna_dediff_axes_boxplots.png"),
  file.path(rna_dir, "01_rna_dediff_axes_boxplots.pdf")
)
plot_boxplots(
  prot_scores,
  prot_stats,
  "Proteomics",
  file.path(proteomics_dir, "01_proteomics_dediff_axes_boxplots.png"),
  file.path(proteomics_dir, "01_proteomics_dediff_axes_boxplots.pdf")
)
plot_overview_heatmap(
  all_stats,
  file.path(integrated_dir, "02_dediff_axes_overview_heatmap.png"),
  file.path(integrated_dir, "02_dediff_axes_overview_heatmap.pdf")
)
plot_group_mean_heatmap(
  all_scores,
  file.path(integrated_dir, "03_dediff_axes_group_mean_heatmap.png"),
  file.path(integrated_dir, "03_dediff_axes_group_mean_heatmap.pdf")
)
plot_delta_concordance(
  all_stats,
  file.path(integrated_dir, "04_dediff_axes_delta_concordance.png"),
  file.path(integrated_dir, "04_dediff_axes_delta_concordance.pdf")
)

summary_lines <- c(
  "WT host vs FAP-deficient host xenograft dedifferentiation-axis summary",
  paste("Output directory:", outdir),
  "",
  "Method:",
  "- Scores were computed from log2(expression + 1) matrices.",
  "- Each gene was row-z-scored within modality across samples.",
  "- Signature score = mean row-z of genes in the signature.",
  "- Default minimum coverage for scoring = 3 genes per signature.",
  "- BRAF-RAS score = BRAF-like score minus RAS-like lineage score; reported only when both components met coverage threshold.",
  "",
  paste("Top RNA axes higher in FAP-deficient host:", format_top_axes(rna_stats, "up")),
  paste("Top RNA axes higher in WT host:", format_top_axes(rna_stats, "down")),
  paste("Top proteomics axes higher in FAP-deficient host:", format_top_axes(prot_stats, "up")),
  paste("Top proteomics axes higher in WT host:", format_top_axes(prot_stats, "down")),
  "",
  "Coverage notes:",
  paste0("- RNA TDS / lineage identity coverage: ", rna_res$coverage_df$genes_present[match("TDS / lineage identity", rna_res$coverage_df$signature)], "/", rna_res$coverage_df$genes_requested[match("TDS / lineage identity", rna_res$coverage_df$signature)]),
  paste0("- Proteomics TDS / lineage identity coverage: ", prot_res$coverage_df$genes_present[match("TDS / lineage identity", prot_res$coverage_df$signature)], "/", prot_res$coverage_df$genes_requested[match("TDS / lineage identity", prot_res$coverage_df$signature)]),
  paste0("- Proteomics BRAF-like coverage: ", prot_res$coverage_df$genes_present[match("BRAF-like program", prot_res$coverage_df$signature)], "/", prot_res$coverage_df$genes_requested[match("BRAF-like program", prot_res$coverage_df$signature)]),
  paste0("- Proteomics RAS-like coverage: ", prot_res$coverage_df$genes_present[match("RAS-like lineage program", prot_res$coverage_df$signature)], "/", prot_res$coverage_df$genes_requested[match("RAS-like lineage program", prot_res$coverage_df$signature)]),
  paste0("- Proteomics Notch coverage: ", prot_res$coverage_df$genes_present[match("Notch signaling", prot_res$coverage_df$signature)], "/", prot_res$coverage_df$genes_requested[match("Notch signaling", prot_res$coverage_df$signature)]),
  paste0("- Proteomics Hedgehog coverage: ", prot_res$coverage_df$genes_present[match("Hedgehog signaling", prot_res$coverage_df$signature)], "/", prot_res$coverage_df$genes_requested[match("Hedgehog signaling", prot_res$coverage_df$signature)])
)

writeLines(summary_lines, file.path(outdir, "summary_notes.txt"))

message("Dedifferentiation-axis analysis completed. Results saved to: ", outdir)
