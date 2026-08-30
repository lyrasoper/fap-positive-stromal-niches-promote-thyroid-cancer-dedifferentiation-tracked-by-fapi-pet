#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Figure 4D and Supplementary Fig. S18 - Hallmark module scoring.

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
  library(ggplot2)
  library(pheatmap)
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

outdir <- file.path(root, "outputs", "wt_vs_fapko_bulk_omics_20260407", "hallmark")
rna_dir <- file.path(outdir, "rna")
proteomics_dir <- file.path(outdir, "proteomics")
integrated_dir <- file.path(outdir, "integrated")
dir.create(rna_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proteomics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(integrated_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("WT host", "FAP-deficient host")
group_colors <- c("WT host" = "#4C78A8", "FAP-deficient host" = "#E45756")

hallmark_gmt <- file.path(root, "reference", "mh.all.v2026.1.Mm.symbols.gmt")

focused_hallmarks <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_HYPOXIA",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_APOPTOSIS"
)

hallmark_theme_blocks <- list(
  `Stroma / EMT` = c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_ANGIOGENESIS",
    "HALLMARK_HYPOXIA",
    "HALLMARK_COAGULATION",
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_APICAL_SURFACE"
  ),
  `Inflammation / Immune` = c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_ALLOGRAFT_REJECTION"
  ),
  `Growth Signaling` = c(
    "HALLMARK_KRAS_SIGNALING_UP",
    "HALLMARK_KRAS_SIGNALING_DN",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    "HALLMARK_NOTCH_SIGNALING",
    "HALLMARK_HEDGEHOG_SIGNALING",
    "HALLMARK_ANDROGEN_RESPONSE",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE"
  ),
  `Proliferation / Stress` = c(
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MITOTIC_SPINDLE",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_UV_RESPONSE_UP",
    "HALLMARK_UV_RESPONSE_DN",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
  ),
  Metabolism = c(
    "HALLMARK_GLYCOLYSIS",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_FATTY_ACID_METABOLISM",
    "HALLMARK_ADIPOGENESIS",
    "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
    "HALLMARK_BILE_ACID_METABOLISM",
    "HALLMARK_PEROXISOME",
    "HALLMARK_HEME_METABOLISM",
    "HALLMARK_XENOBIOTIC_METABOLISM",
    "HALLMARK_PROTEIN_SECRETION"
  ),
  `Tissue Identity` = c(
    "HALLMARK_PANCREAS_BETA_CELLS",
    "HALLMARK_MYOGENESIS",
    "HALLMARK_SPERMATOGENESIS"
  )
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
  rownames(mat) <- expr_df[[gene_col]]
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  storage.mode(mat) <- "numeric"
  mat
}

parse_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parsed <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    list(name = parts[1], genes = unique(toupper(parts[-c(1, 2)])))
  })
  names(parsed) <- vapply(parsed, `[[`, character(1), "name")
  lapply(parsed, `[[`, "genes")
}

pretty_hallmark_name <- function(x) {
  pretty_map <- c(
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
    HALLMARK_TGF_BETA_SIGNALING = "TGF-beta",
    HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFa/NFkB",
    HALLMARK_IL6_JAK_STAT3_SIGNALING = "IL6-JAK-STAT3",
    HALLMARK_IL2_STAT5_SIGNALING = "IL2-STAT5",
    HALLMARK_PI3K_AKT_MTOR_SIGNALING = "PI3K-AKT-mTOR",
    HALLMARK_MTORC1_SIGNALING = "mTORC1",
    HALLMARK_WNT_BETA_CATENIN_SIGNALING = "WNT/beta-catenin",
    HALLMARK_KRAS_SIGNALING_UP = "KRAS up",
    HALLMARK_KRAS_SIGNALING_DN = "KRAS down",
    HALLMARK_E2F_TARGETS = "E2F targets",
    HALLMARK_G2M_CHECKPOINT = "G2M checkpoint",
    HALLMARK_P53_PATHWAY = "p53 pathway",
    HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = "ROS pathway",
    HALLMARK_INTERFERON_ALPHA_RESPONSE = "IFNa response",
    HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFNg response",
    HALLMARK_ESTROGEN_RESPONSE_EARLY = "Estrogen early",
    HALLMARK_ESTROGEN_RESPONSE_LATE = "Estrogen late",
    HALLMARK_PROTEIN_SECRETION = "Protein secretion"
  )
  mapped <- pretty_map[x]
  fallback <- gsub("_", " ", sub("^HALLMARK_", "", x))
  ifelse(is.na(mapped), fallback, mapped)
}

safe_scale_rows <- function(mat) {
  out <- t(scale(t(mat)))
  out[!is.finite(out)] <- 0
  out
}

rank_percentile_matrix <- function(mat) {
  rank_mat <- apply(mat, 2, function(x) {
    n_non_na <- sum(!is.na(x))
    if (n_non_na == 0) {
      return(rep(NA_real_, length(x)))
    }
    rank(x, ties.method = "average", na.last = "keep") / n_non_na
  })
  if (is.null(dim(rank_mat))) {
    rank_mat <- matrix(rank_mat, ncol = 1)
  }
  rownames(rank_mat) <- rownames(mat)
  colnames(rank_mat) <- colnames(mat)
  rank_mat
}

compute_signature_scores <- function(mat, signatures, min_genes = 4) {
  rank_mat <- rank_percentile_matrix(mat)
  coverage_rows <- vector("list", length(signatures))
  score_rows <- vector("list", length(signatures))
  names(score_rows) <- names(signatures)

  for (i in seq_along(signatures)) {
    sig_name <- names(signatures)[i]
    genes <- unique(signatures[[i]])
    present <- intersect(genes, rownames(mat))
    missing <- setdiff(genes, rownames(mat))
    sub_ranks <- rank_mat[present, , drop = FALSE]
    non_missing_counts <- colSums(!is.na(sub_ranks))
    score <- if (length(present) == 0) {
      rep(NA_real_, ncol(mat))
    } else {
      colMeans(sub_ranks, na.rm = TRUE) * 2 - 1
    }
    score[non_missing_counts < min_genes] <- NA_real_
    score_rows[[i]] <- score
    coverage_rows[[i]] <- data.frame(
      signature = sig_name,
      genes_requested = length(genes),
      genes_present = length(present),
      genes_missing = length(missing),
      coverage_pct = 100 * length(present) / length(genes),
      present_genes = paste(present, collapse = ","),
      missing_genes = paste(missing, collapse = ","),
      stringsAsFactors = FALSE
    )
  }

  score_mat <- do.call(rbind, score_rows)
  rownames(score_mat) <- names(signatures)
  colnames(score_mat) <- colnames(mat)
  list(
    score_mat = score_mat,
    coverage_df = bind_rows(coverage_rows)
  )
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

order_hallmarks <- function(sig_names) {
  block_order <- unlist(hallmark_theme_blocks, use.names = FALSE)
  c(intersect(block_order, sig_names), setdiff(sort(sig_names), block_order))
}

plot_boxplots <- function(score_df, stats_df, dataset_name, path_png, path_pdf) {
  ann <- stats_df |>
    filter(signature %in% focused_hallmarks) |>
    mutate(
      sig_label = vapply(p_value, p_to_symbol, character(1)),
      signature_label = pretty_hallmark_name(signature)
    ) |>
    select(signature, signature_label, sig_label)

  bounds <- score_df |>
    filter(signature %in% focused_hallmarks) |>
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
  plot_colors <- c("WT host" = "#D47A6A", "FAP-deficient host" = "#46B58A")

  plot_df <- score_df |>
    filter(signature %in% focused_hallmarks) |>
    mutate(
      signature = factor(signature, levels = focused_hallmarks),
      signature_label = pretty_hallmark_name(signature)
    )

  ann <- ann |>
    mutate(signature = factor(signature, levels = focused_hallmarks))

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
    facet_wrap(~ signature_label, scales = "free_y", ncol = 7) +
    scale_fill_manual(values = plot_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.24))) +
    labs(
      title = paste0(dataset_name, " Hallmark scores"),
      subtitle = "Focused Hallmarks, rank-based activity scores shown as modality-internal z-scores",
      x = NULL,
      y = "Hallmark score z"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1, face = "bold", size = 8.8),
      axis.text.y = element_text(size = 8.2),
      axis.title.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0, size = 12),
      plot.subtitle = element_text(size = 9.5),
      strip.text = element_text(face = "bold", size = 8.7, margin = margin(4, 5, 4, 5)),
      strip.background = element_rect(fill = "white", colour = "#2C2C2C", linewidth = 0.8),
      panel.spacing = unit(0.65, "lines"),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.6, colour = "#2C2C2C")
    )

  ggsave(path_png, p, width = 20.0, height = 7.6, dpi = 300)
  ggsave(path_pdf, p, width = 20.0, height = 7.6)
}

plot_heatmap_dual <- function(mat, annotation_col, title, path_png, path_pdf, width = 8.5, height = 10.0, cluster_rows = FALSE, cluster_cols = FALSE) {
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
    fontsize_row = 8,
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
    fontsize_row = 8,
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

  label_df <- bind_rows(
    delta_df |> arrange(desc(abs(rna_delta_z))) |> slice_head(n = 5),
    delta_df |> arrange(desc(abs(proteomics_delta_z))) |> slice_head(n = 5)
  ) |>
    distinct(signature, .keep_all = TRUE)

  p <- ggplot(delta_df, aes(rna_delta_z, proteomics_delta_z)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(size = 2.8, colour = "#7C3AED", alpha = 0.85) +
    geom_text(
      data = label_df,
      aes(label = pretty_hallmark_name(signature)),
      size = 2.8,
      check_overlap = TRUE,
      nudge_y = 0.03
    ) +
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
      title = "Hallmark delta concordance across RNA and proteomics",
      subtitle = "Delta = mean z-score in FAP-deficient host minus WT host",
      x = "RNA delta z",
      y = "Proteomics delta z"
    )

  ggsave(path_png, p, width = 6.8, height = 5.8, dpi = 300)
  ggsave(path_pdf, p, width = 6.8, height = 5.8)
}

format_top_signatures <- function(stats_df, direction = c("up", "down"), n = 6) {
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
      transmute(label = paste0(pretty_hallmark_name(signature), "(delta=", sprintf("%.2f", delta_fap_minus_wt_raw), ", q=", formatC(q_value, format = "e", digits = 1), ")")) |>
      pull(label),
    collapse = "; "
  )
}

run_hallmark_modality <- function(expr_df, gene_col, sample_cols, dataset_name, output_dir, signatures, signature_order) {
  expr_mat <- prepare_matrix(expr_df, gene_col, sample_cols)
  hallmark_res <- compute_signature_scores(expr_mat, signatures, min_genes = 4)
  raw_scores <- hallmark_res$score_mat[signature_order, , drop = FALSE]
  z_scores <- safe_scale_rows(raw_scores)

  sample_meta <- tibble(
    sample_raw = colnames(raw_scores),
    sample = clean_sample_name(colnames(raw_scores)),
    group = factor(label_group(colnames(raw_scores)), levels = group_levels)
  ) |>
    arrange(group, sample)

  raw_scores <- raw_scores[, sample_meta$sample_raw, drop = FALSE]
  z_scores <- z_scores[, sample_meta$sample_raw, drop = FALSE]

  raw_long <- as.data.frame(raw_scores) |>
    tibble::rownames_to_column("signature") |>
    pivot_longer(-signature, names_to = "sample_raw", values_to = "score_raw") |>
    mutate(
      sample = clean_sample_name(sample_raw),
      group = factor(label_group(sample_raw), levels = group_levels),
      dataset = dataset_name
    )

  z_long <- as.data.frame(z_scores) |>
    tibble::rownames_to_column("signature") |>
    pivot_longer(-signature, names_to = "sample_raw", values_to = "score_z")

  score_df <- left_join(raw_long, z_long, by = c("signature", "sample_raw")) |>
    mutate(signature_label = pretty_hallmark_name(signature))

  stats_df <- score_df |>
    group_by(signature) |>
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
      dataset = dataset_name,
      signature_label = pretty_hallmark_name(signature)
    ) |>
    arrange(match(signature, signature_order))

  coverage_df <- hallmark_res$coverage_df |>
    mutate(
      dataset = dataset_name,
      signature_label = pretty_hallmark_name(signature)
    ) |>
    arrange(match(signature, signature_order))

  sample_ann <- data.frame(Group = sample_meta$group)
  rownames(sample_ann) <- sample_meta$sample

  sample_heatmap_mat <- z_scores
  colnames(sample_heatmap_mat) <- sample_meta$sample
  rownames(sample_heatmap_mat) <- pretty_hallmark_name(rownames(sample_heatmap_mat))

  group_mean_mat <- score_df |>
    group_by(signature, group) |>
    summarise(mean_z = mean(score_z, na.rm = TRUE), .groups = "drop") |>
    mutate(group = factor(group, levels = group_levels)) |>
    arrange(match(signature, signature_order), group) |>
    select(signature, group, mean_z) |>
    pivot_wider(names_from = group, values_from = mean_z) |>
    as.data.frame()
  rownames(group_mean_mat) <- pretty_hallmark_name(group_mean_mat$signature)
  group_mean_mat <- as.matrix(group_mean_mat[, group_levels, drop = FALSE])

  group_ann <- data.frame(Group = factor(group_levels, levels = group_levels))
  rownames(group_ann) <- colnames(group_mean_mat)

  plot_boxplots(
    score_df,
    stats_df,
    dataset_name,
    file.path(output_dir, paste0("01_", tolower(dataset_name), "_hallmark_boxplots_focused.png")),
    file.path(output_dir, paste0("01_", tolower(dataset_name), "_hallmark_boxplots_focused.pdf"))
  )
  plot_heatmap_dual(
    sample_heatmap_mat,
    sample_ann,
    paste0(dataset_name, ": sample-level Hallmark heatmap"),
    file.path(output_dir, paste0("02_", tolower(dataset_name), "_hallmark_heatmap_samples.png")),
    file.path(output_dir, paste0("02_", tolower(dataset_name), "_hallmark_heatmap_samples.pdf"))
  )
  plot_heatmap_dual(
    group_mean_mat,
    group_ann,
    paste0(dataset_name, ": group mean Hallmark heatmap"),
    file.path(output_dir, paste0("03_", tolower(dataset_name), "_hallmark_heatmap_group_means.png")),
    file.path(output_dir, paste0("03_", tolower(dataset_name), "_hallmark_heatmap_group_means.pdf")),
    width = 5.6,
    height = 10.0
  )

  write_table_tsv(score_df, file.path(output_dir, paste0(tolower(dataset_name), "_hallmark_scores_per_sample.tsv")))
  write_table_tsv(stats_df, file.path(output_dir, paste0(tolower(dataset_name), "_hallmark_pathway_stats.tsv")))
  write_table_tsv(coverage_df, file.path(output_dir, paste0(tolower(dataset_name), "_hallmark_coverage.tsv")))

  list(
    score_df = score_df,
    stats_df = stats_df,
    coverage_df = coverage_df,
    group_mean_mat = group_mean_mat
  )
}

hallmark_signatures <- parse_gmt(hallmark_gmt)
signature_order <- order_hallmarks(names(hallmark_signatures))

sample_info <- read_sample_info(file.path(root, "002_DIA_Summary", "Sample_Info.txt"))
rna_expr <- read_rna_expression(file.path(root, "001_mRNA_Summary", "1.GeneExpression", "1_genes_fpkm_expression.txt"))
rna_samples <- grep("^FPKM\\.", names(rna_expr), value = TRUE)
prot_expr <- read_protein_expression(file.path(root, "002_DIA_Summary", "01.RawData", "MatchResult", "protein.groups.intensity.txt"), sample_info)
prot_samples <- setdiff(names(prot_expr), "gene_name")

rna_res <- run_hallmark_modality(rna_expr, "gene_name", rna_samples, "RNA", rna_dir, hallmark_signatures, signature_order)
prot_res <- run_hallmark_modality(prot_expr, "gene_name", prot_samples, "Proteomics", proteomics_dir, hallmark_signatures, signature_order)

integrated_group_mean_df <- bind_rows(
  as.data.frame(rna_res$group_mean_mat) |>
    tibble::rownames_to_column("signature_label") |>
    mutate(dataset = "RNA"),
  as.data.frame(prot_res$group_mean_mat) |>
    tibble::rownames_to_column("signature_label") |>
    mutate(dataset = "Proteomics")
) |>
  pivot_longer(all_of(group_levels), names_to = "group", values_to = "mean_z") |>
  mutate(column_id = paste(dataset, group, sep = " | ")) |>
  select(signature_label, column_id, mean_z) |>
  pivot_wider(names_from = column_id, values_from = mean_z) |>
  as.data.frame()

rownames(integrated_group_mean_df) <- integrated_group_mean_df$signature_label
integrated_group_mean_mat <- as.matrix(integrated_group_mean_df[, c(
  "RNA | WT host",
  "RNA | FAP-deficient host",
  "Proteomics | WT host",
  "Proteomics | FAP-deficient host"
), drop = FALSE])
integrated_group_mean_mat <- integrated_group_mean_mat[pretty_hallmark_name(signature_order), , drop = FALSE]

integrated_ann <- data.frame(
  Group = factor(c("WT host", "FAP-deficient host", "WT host", "FAP-deficient host"), levels = group_levels)
)
rownames(integrated_ann) <- colnames(integrated_group_mean_mat)

plot_heatmap_dual(
  integrated_group_mean_mat,
  integrated_ann,
  "Integrated Hallmark group mean comparison",
  file.path(integrated_dir, "04_hallmark_group_mean_comparison.png"),
  file.path(integrated_dir, "04_hallmark_group_mean_comparison.pdf")
)

plot_heatmap_dual(
  integrated_group_mean_mat,
  integrated_ann,
  "Integrated Hallmark group mean comparison (row-clustered)",
  file.path(integrated_dir, "06_hallmark_group_mean_comparison_row_clustered.png"),
  file.path(integrated_dir, "06_hallmark_group_mean_comparison_row_clustered.pdf"),
  cluster_rows = TRUE,
  cluster_cols = FALSE
)

delta_df <- inner_join(
  rna_res$stats_df |>
    select(signature, signature_label, rna_delta_z = delta_fap_minus_wt_z, rna_q = q_value),
  prot_res$stats_df |>
    select(signature, proteomics_delta_z = delta_fap_minus_wt_z, proteomics_q = q_value),
  by = "signature"
) |>
  arrange(match(signature, signature_order))

plot_delta_concordance(
  delta_df,
  file.path(integrated_dir, "05_hallmark_delta_concordance.png"),
  file.path(integrated_dir, "05_hallmark_delta_concordance.pdf")
)

write_table_tsv(
  as.data.frame(integrated_group_mean_mat) |>
    tibble::rownames_to_column("signature_label"),
  file.path(integrated_dir, "hallmark_group_mean_comparison.tsv")
)
write_table_tsv(delta_df, file.path(integrated_dir, "hallmark_delta_concordance.tsv"))

summary_lines <- c(
  "WT host vs FAP-deficient host xenograft Hallmark summary",
  paste("Output directory:", outdir),
  paste("RNA outputs:", rna_dir),
  paste("Proteomics outputs:", proteomics_dir),
  paste("Integrated outputs:", integrated_dir),
  "",
  "Method:",
  "- Official Mouse Hallmark GMT (mh.all.v2026.1.Mm.symbols.gmt).",
  "- Rank-based gene set score computed within each sample.",
  "- For each sample, genes were percentile-ranked within modality; score = mean percentile rank of hallmark genes scaled to [-1, 1].",
  "- Hallmark boxplots and heatmaps use modality-internal z-scores of the resulting rank-based scores.",
  "- Global hallmark coverage is reported separately for RNA and proteomics.",
  "",
  paste("Top RNA Hallmarks higher in FAP-deficient host:", format_top_signatures(rna_res$stats_df, "up")),
  paste("Top RNA Hallmarks higher in WT host:", format_top_signatures(rna_res$stats_df, "down")),
  paste("Top proteomics Hallmarks higher in FAP-deficient host:", format_top_signatures(prot_res$stats_df, "up")),
  paste("Top proteomics Hallmarks higher in WT host:", format_top_signatures(prot_res$stats_df, "down"))
)

writeLines(summary_lines, file.path(outdir, "summary_notes.txt"))

message("Hallmark analysis completed. Results saved to: ", outdir)
