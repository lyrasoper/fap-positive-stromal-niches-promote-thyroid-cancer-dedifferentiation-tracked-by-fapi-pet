#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Figure 4C-E upstream - differential RNA/DIA-protein testing, Fap-deficient vs wild-type host.

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggrepel)
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

script_path <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
}

outdir <- file.path(root, "outputs", "wt_vs_fapko_bulk_omics_20260407")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rna_dir <- file.path(outdir, "rna")
proteomics_dir <- file.path(outdir, "proteomics")
integrated_dir <- file.path(outdir, "integrated")
dir.create(rna_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proteomics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(integrated_dir, recursive = TRUE, showWarnings = FALSE)

label_group <- function(x) {
  clean_x <- sub("^FPKM\\.", "", x)
  ifelse(grepl("^FAP_BPC", clean_x), "FAP-deficient host", "WT host")
}

clean_pathway_name <- function(x) {
  x |>
    gsub("^r mmu [0-9]+ ", "", x = _) |>
    gsub("^[0-9]+ ", "", x = _) |>
    gsub("_", " ", x = _) |>
    gsub("\\s+", " ", x = _) |>
    trimws()
}

safe_scale_rows <- function(mat) {
  out <- t(scale(t(mat)))
  out[!is.finite(out)] <- 0
  out
}

prepare_matrix <- function(expr_df, gene_col, sample_cols, uppercase = FALSE) {
  mat <- as.matrix(expr_df[, sample_cols, drop = FALSE])
  rownames(mat) <- if (uppercase) toupper(expr_df[[gene_col]]) else expr_df[[gene_col]]
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  mat <- log2(mat + 1)
  mat
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
    mutate(gene_name = toupper(sub(";.*$", "", Genes))) |>
    select(gene_name, all_of(sample_cols)) |>
    rename(!!!rename_map) |>
    mutate(across(-gene_name, as.numeric)) |>
    filter(!is.na(gene_name), gene_name != "")
  collapse_by_gene(df2, "gene_name", mapped_names)
}

score_signatures <- function(expr_df, gene_col, sample_cols, signatures, dataset_name, uppercase = FALSE) {
  gene_names <- expr_df[[gene_col]]
  if (uppercase) {
    gene_names <- toupper(gene_names)
  }
  expr_mat <- as.matrix(expr_df[, sample_cols, drop = FALSE])
  rownames(expr_mat) <- gene_names
  expr_mat <- log2(expr_mat + 1)
  expr_mat <- safe_scale_rows(expr_mat)

  scores <- lapply(names(signatures), function(sig_name) {
    genes <- signatures[[sig_name]]
    if (uppercase) {
      genes <- toupper(genes)
    }
    present <- intersect(genes, rownames(expr_mat))
    if (length(present) == 0) {
      return(NULL)
    }
    tibble(
      dataset = dataset_name,
      signature = sig_name,
      sample = sample_cols,
      group = label_group(sample_cols),
      score = colMeans(expr_mat[present, , drop = FALSE], na.rm = TRUE),
      genes_present = length(present),
      genes_used = paste(present, collapse = ", ")
    )
  })

  bind_rows(scores)
}

read_rna_diff <- function(path) {
  read.delim(path, check.names = FALSE) |>
    transmute(
      gene = gene_name,
      dataset = "RNA",
      log2_fc = as.numeric(`log2(fc)`),
      p_value = as.numeric(pval),
      q_value = as.numeric(qval),
      regulation = regulation,
      significant = significant
    ) |>
    arrange(q_value, p_value, desc(abs(log2_fc))) |>
    distinct(gene, .keep_all = TRUE)
}

read_protein_diff <- function(path) {
  read.delim(path, check.names = FALSE) |>
    transmute(
      gene = toupper(gene_name),
      dataset = "Proteomics",
      log2_fc = log2(as.numeric(FC)),
      p_value = as.numeric(p_value),
      q_value = as.numeric(q_value),
      regulation = regulation,
      significant = Significant
    ) |>
    arrange(q_value, p_value, desc(abs(log2_fc))) |>
    distinct(gene, .keep_all = TRUE)
}

read_gsea <- function(path, dataset, database) {
  read.delim(path, check.names = FALSE) |>
    transmute(
      dataset = dataset,
      database = database,
      pathway = clean_pathway_name(NAME),
      NES = as.numeric(NES),
      p_value = as.numeric(NOM.pval),
      fdr = as.numeric(FDR.qval),
      direction = ifelse(NES > 0, "Higher in FAP-deficient host", "Higher in WT host")
    )
}

plot_signature_scores <- function(df, path_png, path_pdf) {
  pvals <- df |>
    group_by(dataset, signature) |>
    summarise(
      p_value = tryCatch(wilcox.test(score ~ group)$p.value, error = function(e) NA_real_),
      genes_present = first(genes_present),
      .groups = "drop"
    ) |>
    mutate(label = paste0("n=", genes_present, "\nP=", formatC(p_value, format = "e", digits = 1)))

  ymax <- df |>
    group_by(dataset, signature) |>
    summarise(y = max(score, na.rm = TRUE) + 0.25, .groups = "drop")

  ann <- left_join(pvals, ymax, by = c("dataset", "signature"))

  p <- ggplot(df, aes(x = group, y = score, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.15) +
    geom_point(position = position_jitter(width = 0.08), size = 2.4) +
    geom_text(
      data = ann,
      aes(x = 1.5, y = y, label = label),
      inherit.aes = FALSE,
      size = 3
    ) +
    facet_grid(dataset ~ signature, scales = "free_y") +
    scale_colour_manual(values = c("WT host" = "#4C78A8", "FAP-deficient host" = "#E45756")) +
    labs(
      title = "Curated pathway-signature shifts in WT vs FAP-deficient host xenografts",
      x = NULL,
      y = "Mean row-z score"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 25, hjust = 1)
    )

  ggsave(path_png, p, width = 13, height = 6.5, dpi = 300)
  ggsave(path_pdf, p, width = 13, height = 6.5)
}

plot_marker_heatmap <- function(df, path_png, path_pdf) {
  marker_order <- c(
    "PAX8", "NKX2-1", "TPO", "TG",
    "FAP", "ACTA2", "COL1A1", "COL1A2", "POSTN", "FN1", "SPARC", "TGFB1", "TGFBR1",
    "VIM", "MET", "HGF", "OSMR",
    "SPP1", "GPNMB", "C1QA", "C1QB", "CD274", "CXCL9", "CXCL10"
  )

  plot_df <- df |>
    mutate(
      marker = factor(gene, levels = rev(marker_order)),
      sig_label = case_when(
        q_value < 0.05 ~ "**",
        q_value < 0.25 ~ "*",
        TRUE ~ ""
      )
    ) |>
    filter(!is.na(marker))

  p <- ggplot(plot_df, aes(x = dataset, y = marker, fill = log2_fc)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sig_label), size = 4) +
    scale_fill_gradient2(low = "#3B82F6", mid = "white", high = "#DC2626", midpoint = 0) +
    labs(
      title = "Selected marker changes across RNA and proteomics",
      subtitle = "* q < 0.25, ** q < 0.05",
      x = NULL,
      y = NULL,
      fill = "log2FC\n(FAP/WT)"
    )

  ggsave(path_png, p, width = 5.2, height = 7.5, dpi = 300)
  ggsave(path_pdf, p, width = 5.2, height = 7.5)
}

plot_gsea_top <- function(df, dataset_name, path_png, path_pdf) {
  selected <- df |>
    filter(dataset == dataset_name) |>
    group_by(direction, pathway) |>
    slice_min(order_by = fdr, n = 1, with_ties = FALSE) |>
    ungroup() |>
    group_by(direction) |>
    arrange(fdr, desc(abs(NES))) |>
    slice_head(n = 10) |>
    ungroup() |>
    mutate(pathway_label = paste0(pathway, " [", database, "]")) |>
    group_by(direction) |>
    mutate(pathway_label = fct_reorder(pathway_label, NES)) |>
    ungroup()

  p <- ggplot(selected, aes(x = NES, y = pathway_label, size = -log10(pmax(fdr, 1e-6)), colour = database)) +
    geom_point(alpha = 0.9) +
    facet_wrap(~ direction, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = c("KEGG" = "#4C78A8", "Reactome" = "#72B7B2")) +
    labs(
      title = paste0(dataset_name, ": top ranked pathway changes"),
      subtitle = "Positive NES means higher in FAP-deficient host; negative NES means higher in WT host",
      x = "NES",
      y = NULL,
      size = "-log10(FDR)"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold")
    )

  ggsave(path_png, p, width = 10.5, height = 8.5, dpi = 300)
  ggsave(path_pdf, p, width = 10.5, height = 8.5)
}

plot_pca <- function(expr_df, gene_col, sample_cols, dataset_name, path_png, path_pdf, uppercase = FALSE) {
  mat <- prepare_matrix(expr_df, gene_col, sample_cols, uppercase = uppercase)
  vars <- apply(mat, 1, var, na.rm = TRUE)
  vars <- vars[is.finite(vars) & vars > 0]
  keep_n <- min(2000, length(vars))
  if (keep_n < 2) {
    return(invisible(NULL))
  }
  keep <- names(sort(vars, decreasing = TRUE))[seq_len(keep_n)]
  mat2 <- mat[keep, , drop = FALSE]
  mat2 <- mat2[rowSums(!is.finite(mat2)) == 0, , drop = FALSE]
  if (nrow(mat2) < 2) {
    return(invisible(NULL))
  }
  pca <- prcomp(t(mat2), center = TRUE, scale. = TRUE)
  pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pca_df <- tibble(
    sample = rownames(pca$x),
    group = label_group(rownames(pca$x)),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )

  p <- ggplot(pca_df, aes(PC1, PC2, colour = group, label = sample)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3.2, show.legend = FALSE) +
    scale_colour_manual(values = c("WT host" = "#4C78A8", "FAP-deficient host" = "#E45756")) +
    labs(
      title = paste0(dataset_name, ": PCA of sample-level expression"),
      x = paste0("PC1 (", sprintf("%.1f", pct[1]), "%)"),
      y = paste0("PC2 (", sprintf("%.1f", pct[2]), "%)"),
      colour = NULL
    ) +
    theme(legend.position = "top")

  ggsave(path_png, p, width = 6.5, height = 5.2, dpi = 300)
  ggsave(path_pdf, p, width = 6.5, height = 5.2)
}

plot_sample_corr <- function(expr_df, gene_col, sample_cols, dataset_name, path_png, path_pdf, uppercase = FALSE) {
  mat <- prepare_matrix(expr_df, gene_col, sample_cols, uppercase = uppercase)
  corr <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
  ann <- data.frame(Group = label_group(colnames(corr)))
  rownames(ann) <- colnames(corr)
  pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

  png(path_png, width = 6.5, height = 5.8, units = "in", res = 300)
  pheatmap(
    corr,
    color = pal,
    breaks = seq(-1, 1, length.out = 101),
    display_numbers = TRUE,
    number_format = "%.2f",
    annotation_col = ann,
    annotation_row = ann,
    main = paste0(dataset_name, ": sample correlation")
  )
  dev.off()

  pdf(path_pdf, width = 6.5, height = 5.8)
  pheatmap(
    corr,
    color = pal,
    breaks = seq(-1, 1, length.out = 101),
    display_numbers = TRUE,
    number_format = "%.2f",
    annotation_col = ann,
    annotation_row = ann,
    main = paste0(dataset_name, ": sample correlation")
  )
  dev.off()
}

plot_volcano <- function(diff_df, dataset_name, path_png, path_pdf, focus_genes = character()) {
  df <- diff_df |>
    filter(dataset == dataset_name) |>
    mutate(
      neg_log10_p = -log10(pmax(p_value, 1e-300)),
      sig_class = case_when(
        significant == "yes" & log2_fc > 0 ~ "Up in FAP-deficient host",
        significant == "yes" & log2_fc < 0 ~ "Up in WT host",
        TRUE ~ "NS"
      )
    )

  top_genes <- df |>
    filter(significant == "yes") |>
    arrange(p_value, desc(abs(log2_fc))) |>
    slice_head(n = 12) |>
    pull(gene)

  label_df <- df |>
    filter(toupper(gene) %in% unique(c(toupper(focus_genes), toupper(top_genes))))

  p <- ggplot(df, aes(log2_fc, neg_log10_p, colour = sig_class)) +
    geom_point(alpha = 0.8, size = 1.6) +
    geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60") +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 3,
      max.overlaps = 50,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = c(
        "Up in WT host" = "#3B82F6",
        "NS" = "grey75",
        "Up in FAP-deficient host" = "#DC2626"
      )
    ) +
    labs(
      title = paste0(dataset_name, ": differential expression volcano"),
      x = "log2FC (FAP-deficient / WT)",
      y = "-log10(P value)",
      colour = NULL
    ) +
    theme(legend.position = "top")

  ggsave(path_png, p, width = 7.5, height = 5.6, dpi = 300)
  ggsave(path_pdf, p, width = 7.5, height = 5.6)
}

plot_top_heatmap <- function(expr_df, gene_col, sample_cols, diff_df, dataset_name, path_png, path_pdf, uppercase = FALSE) {
  mat <- prepare_matrix(expr_df, gene_col, sample_cols, uppercase = uppercase)
  diff_sub <- diff_df |>
    filter(dataset == dataset_name, significant == "yes") |>
    arrange(p_value, desc(abs(log2_fc))) |>
    slice_head(n = 40)
  genes <- intersect(diff_sub$gene, rownames(mat))
  if (length(genes) < 10) {
    vars <- apply(mat, 1, var, na.rm = TRUE)
    genes <- names(sort(vars, decreasing = TRUE))[seq_len(min(40, length(vars)))]
  }
  plot_mat <- safe_scale_rows(mat[genes, , drop = FALSE])
  ann <- data.frame(Group = label_group(colnames(plot_mat)))
  rownames(ann) <- colnames(plot_mat)
  pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

  png(path_png, width = 7.2, height = 9, units = "in", res = 300)
  pheatmap(
    plot_mat,
    color = pal,
    annotation_col = ann,
    show_rownames = TRUE,
    fontsize_row = 7,
    main = paste0(dataset_name, ": top differential features")
  )
  dev.off()

  pdf(path_pdf, width = 7.2, height = 9)
  pheatmap(
    plot_mat,
    color = pal,
    annotation_col = ann,
    show_rownames = TRUE,
    fontsize_row = 7,
    main = paste0(dataset_name, ": top differential features")
  )
  dev.off()
}

plot_marker_boxplots <- function(expr_df, gene_col, sample_cols, markers, dataset_name, path_png, path_pdf, uppercase = FALSE) {
  mat <- prepare_matrix(expr_df, gene_col, sample_cols, uppercase = uppercase)
  row_map <- setNames(rownames(mat), toupper(rownames(mat)))
  present <- unname(row_map[intersect(toupper(markers), names(row_map))])
  df <- as.data.frame(mat[present, , drop = FALSE]) |>
    tibble::rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "sample", values_to = "log_expr") |>
    mutate(group = label_group(sample)) |>
    filter(is.finite(log_expr))

  p <- ggplot(df, aes(group, log_expr, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.2) +
    geom_point(position = position_jitter(width = 0.08), size = 2.1) +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = c("WT host" = "#4C78A8", "FAP-deficient host" = "#E45756")) +
    labs(
      title = paste0(dataset_name, ": selected marker expression"),
      x = NULL,
      y = "log2(expression + 1)"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 25, hjust = 1),
      strip.text = element_text(face = "bold", size = 9)
    )

  ggsave(path_png, p, width = 11, height = 8.5, dpi = 300)
  ggsave(path_pdf, p, width = 11, height = 8.5)
}

plot_fc_concordance <- function(rna_diff, prot_diff, path_png, path_pdf, label_genes = character()) {
  df <- inner_join(
    rna_diff |> mutate(gene = toupper(gene)) |> select(gene, rna_log2_fc = log2_fc, rna_p = p_value),
    prot_diff |> select(gene, prot_log2_fc = log2_fc, prot_p = p_value),
    by = "gene"
  ) |>
    mutate(
      highlight = gene %in% toupper(label_genes),
      corr_label = paste0("Spearman rho = ", sprintf("%.2f", cor(rna_log2_fc, prot_log2_fc, method = "spearman")))
    )

  label_df <- df |>
    filter(highlight) |>
    distinct(gene, .keep_all = TRUE)

  p <- ggplot(df, aes(rna_log2_fc, prot_log2_fc)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(alpha = 0.35, size = 1.4, colour = "#4B5563") +
    geom_point(data = label_df, colour = "#D97706", size = 2.2) +
    ggrepel::geom_text_repel(data = label_df, aes(label = gene), size = 3.2) +
    annotate("text", x = min(df$rna_log2_fc, na.rm = TRUE), y = max(df$prot_log2_fc, na.rm = TRUE),
             label = df$corr_label[1], hjust = 0, vjust = 1, size = 3.5) +
    labs(
      title = "RNA-proteomics fold-change concordance",
      x = "RNA log2FC (FAP-deficient / WT)",
      y = "Proteomics log2FC (FAP-deficient / WT)"
    )

  ggsave(path_png, p, width = 6.8, height = 5.8, dpi = 300)
  ggsave(path_pdf, p, width = 6.8, height = 5.8)
}

classify_tumor_axis <- function(pathway) {
  p <- tolower(pathway)
  case_when(
    grepl("cell cycle|dna replication|homologous recombination|fanconi|nucleotide excision repair|dna repair", p) ~ "Cell cycle / DNA repair",
    grepl("ecm|extracellular matrix|focal adhesion|integrin", p) ~ "ECM / adhesion",
    grepl("tgf beta", p) ~ "TGF-beta",
    grepl("wnt", p) ~ "WNT",
    grepl("notch", p) ~ "NOTCH",
    grepl("mapk|pi3k|akt|receptor tyrosine kinase|fgfr", p) ~ "RTK-MAPK / PI3K-AKT",
    grepl("vegf|vascular permeability|angiogenesis", p) ~ "Angiogenesis / VEGF",
    grepl("hif|hypoxia|glycolysis|pyruvate|central carbon metabolism|mtor", p) ~ "Hypoxia / metabolism",
    grepl("apoptosis|tp53|p53", p) ~ "Apoptosis / p53",
    TRUE ~ NA_character_
  )
}

plot_tumor_gsea_focus <- function(df, path_png, path_pdf) {
  plot_df <- df |>
    mutate(axis = classify_tumor_axis(pathway)) |>
    filter(!is.na(axis)) |>
    group_by(dataset, axis, direction) |>
    slice_min(order_by = fdr, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      axis = factor(axis, levels = c(
        "Cell cycle / DNA repair", "ECM / adhesion", "TGF-beta", "WNT", "NOTCH",
        "RTK-MAPK / PI3K-AKT", "Angiogenesis / VEGF", "Hypoxia / metabolism", "Apoptosis / p53"
      )),
      dataset = factor(dataset, levels = c("RNA", "Proteomics"))
    )

  p <- ggplot(plot_df, aes(dataset, axis, fill = NES, size = -log10(pmax(fdr, 1e-6)))) +
    geom_point(shape = 21, colour = "black", alpha = 0.9) +
    geom_text(aes(label = database), size = 2.8, vjust = 0.5) +
    scale_fill_gradient2(low = "#3B82F6", mid = "white", high = "#DC2626", midpoint = 0) +
    labs(
      title = "Tumor-focused pathway summary from GSEA",
      subtitle = "Each bubble shows the most significant pathway within a tumor-relevant axis; positive NES means higher in FAP-deficient host",
      x = NULL,
      y = NULL,
      fill = "NES",
      size = "-log10(FDR)"
    )

  ggsave(path_png, p, width = 8.4, height = 6.8, dpi = 300)
  ggsave(path_pdf, p, width = 8.4, height = 6.8)
}

write_summary <- function(path, summary_lines) {
  writeLines(summary_lines, con = path)
}

signatures <- list(
  "Thyroid differentiation" = c("Pax8", "Nkx2-1", "Tg", "Tpo", "Slc5a5", "Foxe1"),
  "CAF / ECM fibrosis" = c("Fap", "Col1a1", "Col1a2", "Col3a1", "Fn1", "Postn", "Sparc", "Acta2", "Tgfb1", "Tgfbr1"),
  "EMT / invasion" = c("Vim", "Cdh2", "Snai1", "Snai2", "Zeb1", "Zeb2", "Twist1", "Cd44", "Itga5", "Mmp2", "Mmp9"),
  "RTK / MAPK-PI3K core" = c("Hgf", "Met", "Gab1", "Grb2", "Pik3r1", "Pik3ca", "Pik3cb", "Akt1", "Akt2", "Ptpn11", "Stat3"),
  "Myeloid / complement" = c("C1qa", "C1qb", "C1qc", "Spp1", "Gpnmb", "Lyz2", "Tyrobp", "Fcgr3", "C3ar1", "Cxcl9", "Cxcl10", "Osmr")
)

tumor_signatures <- list(
  "Thyroid differentiation" = c("Pax8", "Nkx2-1", "Tg", "Tpo", "Slc5a5", "Foxe1"),
  "Cell cycle / proliferation" = c("Mki67", "Ccna2", "Ccnb1", "Cdk1", "Top2a", "Ube2c", "Pcna", "Birc5", "Mcm2", "Mcm5", "Rad51", "Brca1"),
  "EMT / invasion" = c("Vim", "Cdh2", "Snai1", "Snai2", "Zeb1", "Zeb2", "Twist1", "Itga5", "Fn1", "Postn", "Mmp2", "Mmp9"),
  "ECM / TGF-beta" = c("Fap", "Col1a1", "Col1a2", "Fn1", "Postn", "Sparc", "Acta2", "Tgfb1", "Tgfbr1", "Serpine1"),
  "RTK / MAPK" = c("Egfr", "Erbb2", "Kras", "Raf1", "Met", "Gab1", "Grb2", "Map2k1", "Map2k2", "Mapk1", "Mapk3", "Dusp6"),
  "PI3K / AKT / mTOR" = c("Pik3ca", "Pik3cb", "Pik3r1", "Akt1", "Akt2", "Mtor", "Rptor", "Rps6kb1", "Eif4ebp1", "Ptpn11"),
  "Hypoxia / glycolysis" = c("Hif1a", "Slc2a1", "Hk2", "Pfkp", "Aldoa", "Ldha", "Pdk1", "Pgk1", "Eno1"),
  "Angiogenesis / VEGF" = c("Vegfa", "Kdr", "Flt1", "Pecam1", "Cdh5", "Eng", "Angpt2"),
  "Apoptosis / p53" = c("Trp53", "Cdkn1a", "Bax", "Bak1", "Bcl2", "Bcl2l1", "Pmaip1", "Casp3", "Casp8")
)

sample_info <- read_sample_info(file.path(root, "002_DIA_Summary", "Sample_Info.txt"))

rna_expr <- read_rna_expression(file.path(root, "001_mRNA_Summary", "1.GeneExpression", "1_genes_fpkm_expression.txt"))
rna_samples <- grep("^FPKM\\.", names(rna_expr), value = TRUE)

prot_expr <- read_protein_expression(file.path(root, "002_DIA_Summary", "01.RawData", "MatchResult", "protein.groups.intensity.txt"), sample_info)
prot_samples <- setdiff(names(prot_expr), "gene_name")

rna_scores <- score_signatures(rna_expr, "gene_name", rna_samples, signatures, "RNA", uppercase = FALSE) |>
  mutate(sample = sub("^FPKM\\.", "", sample))
prot_scores <- score_signatures(prot_expr, "gene_name", prot_samples, signatures, "Proteomics", uppercase = TRUE)

rna_tumor_scores <- score_signatures(rna_expr, "gene_name", rna_samples, tumor_signatures, "RNA", uppercase = FALSE) |>
  mutate(sample = sub("^FPKM\\.", "", sample))
prot_tumor_scores <- score_signatures(prot_expr, "gene_name", prot_samples, tumor_signatures, "Proteomics", uppercase = TRUE)

signature_scores <- bind_rows(rna_scores, prot_scores) |>
  mutate(
    signature = factor(signature, levels = names(signatures)),
    dataset = factor(dataset, levels = c("RNA", "Proteomics")),
    group = factor(group, levels = c("WT host", "FAP-deficient host"))
  )

tumor_signature_scores <- bind_rows(rna_tumor_scores, prot_tumor_scores) |>
  mutate(
    signature = factor(signature, levels = names(tumor_signatures)),
    dataset = factor(dataset, levels = c("RNA", "Proteomics")),
    group = factor(group, levels = c("WT host", "FAP-deficient host"))
  )

rna_diff <- read_rna_diff(file.path(root, "001_mRNA_Summary", "3.DiffExpression", "FAP_BPCVSWT_BPC", "FAP_BPCVSWT_BPC_Gene_differential_expression.txt"))
prot_diff <- read_protein_diff(file.path(root, "002_DIA_Summary", "04.DiffExp", "COND1", "FAP_BPCVSWT_BPC", "FAP_BPCVSWT_BPC_diff_annotation.txt"))

selected_markers <- c(
  "PAX8", "NKX2-1", "TPO", "TG", "FAP", "ACTA2", "COL1A1", "COL1A2", "POSTN", "FN1", "SPARC",
  "TGFB1", "TGFBR1", "VIM", "MET", "HGF", "OSMR", "SPP1", "GPNMB", "C1QA", "C1QB", "CD274", "CXCL9", "CXCL10"
)

marker_panel_boxplot <- c("FAP", "PAX8", "TG", "TPO", "VIM", "MET", "COL1A1", "POSTN", "SPP1", "GPNMB", "C1QA", "OSMR")

marker_changes <- bind_rows(
  rna_diff |> mutate(gene = toupper(gene)),
  prot_diff
) |>
  filter(gene %in% selected_markers) |>
  distinct(gene, dataset, .keep_all = TRUE)

all_gsea <- bind_rows(
  read_gsea(file.path(root, "001_mRNA_Summary", "4.Gsea", "FAP_BPCVSWT_BPC", "KEGG", "FAP_BPCVSWT_BPC.Gsea.enrichment.KEGG.txt"), "RNA", "KEGG"),
  read_gsea(file.path(root, "001_mRNA_Summary", "4.Gsea", "FAP_BPCVSWT_BPC", "Reactome", "FAP_BPCVSWT_BPC.Gsea.enrichment.Reactome.txt"), "RNA", "Reactome"),
  read_gsea(file.path(root, "002_DIA_Summary", "06.GSEA", "COND1", "FAP_BPCVSWT_BPC", "kegg", "FAP_BPCVSWT_BPC.Gsea.enrichment.KEGG.txt"), "Proteomics", "KEGG"),
  read_gsea(file.path(root, "002_DIA_Summary", "06.GSEA", "COND1", "FAP_BPCVSWT_BPC", "reactome", "FAP_BPCVSWT_BPC.Gsea.enrichment.Reactome.txt"), "Proteomics", "Reactome")
)

plot_pca(
  rna_expr,
  "gene_name",
  rna_samples,
  "RNA",
  file.path(rna_dir, "01_rna_pca.png"),
  file.path(rna_dir, "01_rna_pca.pdf"),
  uppercase = FALSE
)

plot_sample_corr(
  rna_expr,
  "gene_name",
  rna_samples,
  "RNA",
  file.path(rna_dir, "02_rna_sample_correlation.png"),
  file.path(rna_dir, "02_rna_sample_correlation.pdf"),
  uppercase = FALSE
)

plot_volcano(
  rna_diff,
  "RNA",
  file.path(rna_dir, "03_rna_volcano.png"),
  file.path(rna_dir, "03_rna_volcano.pdf"),
  focus_genes = selected_markers
)

plot_top_heatmap(
  rna_expr,
  "gene_name",
  rna_samples,
  rna_diff,
  "RNA",
  file.path(rna_dir, "04_rna_top_diff_heatmap.png"),
  file.path(rna_dir, "04_rna_top_diff_heatmap.pdf"),
  uppercase = FALSE
)

plot_marker_boxplots(
  rna_expr,
  "gene_name",
  rna_samples,
  marker_panel_boxplot,
  "RNA",
  file.path(rna_dir, "05_rna_marker_boxplots.png"),
  file.path(rna_dir, "05_rna_marker_boxplots.pdf"),
  uppercase = FALSE
)

plot_signature_scores(
  signature_scores,
  file.path(integrated_dir, "06_curated_signature_scores.png"),
  file.path(integrated_dir, "06_curated_signature_scores.pdf")
)

plot_marker_heatmap(
  marker_changes,
  file.path(integrated_dir, "07_selected_marker_heatmap.png"),
  file.path(integrated_dir, "07_selected_marker_heatmap.pdf")
)

plot_gsea_top(
  all_gsea,
  "RNA",
  file.path(rna_dir, "08_rna_top_gsea_pathways.png"),
  file.path(rna_dir, "08_rna_top_gsea_pathways.pdf")
)

plot_pca(
  prot_expr,
  "gene_name",
  prot_samples,
  "Proteomics",
  file.path(proteomics_dir, "09_proteomics_pca.png"),
  file.path(proteomics_dir, "09_proteomics_pca.pdf"),
  uppercase = TRUE
)

plot_sample_corr(
  prot_expr,
  "gene_name",
  prot_samples,
  "Proteomics",
  file.path(proteomics_dir, "10_proteomics_sample_correlation.png"),
  file.path(proteomics_dir, "10_proteomics_sample_correlation.pdf"),
  uppercase = TRUE
)

plot_volcano(
  prot_diff,
  "Proteomics",
  file.path(proteomics_dir, "11_proteomics_volcano.png"),
  file.path(proteomics_dir, "11_proteomics_volcano.pdf"),
  focus_genes = selected_markers
)

plot_top_heatmap(
  prot_expr,
  "gene_name",
  prot_samples,
  prot_diff,
  "Proteomics",
  file.path(proteomics_dir, "12_proteomics_top_diff_heatmap.png"),
  file.path(proteomics_dir, "12_proteomics_top_diff_heatmap.pdf"),
  uppercase = TRUE
)

plot_marker_boxplots(
  prot_expr,
  "gene_name",
  prot_samples,
  marker_panel_boxplot,
  "Proteomics",
  file.path(proteomics_dir, "13_proteomics_marker_boxplots.png"),
  file.path(proteomics_dir, "13_proteomics_marker_boxplots.pdf"),
  uppercase = TRUE
)

plot_gsea_top(
  all_gsea,
  "Proteomics",
  file.path(proteomics_dir, "14_proteomics_top_gsea_pathways.png"),
  file.path(proteomics_dir, "14_proteomics_top_gsea_pathways.pdf")
)

plot_fc_concordance(
  rna_diff,
  prot_diff,
  file.path(integrated_dir, "15_rna_proteomics_fc_concordance.png"),
  file.path(integrated_dir, "15_rna_proteomics_fc_concordance.pdf"),
  label_genes = selected_markers
)

plot_signature_scores(
  tumor_signature_scores,
  file.path(integrated_dir, "16_tumor_signature_scores.png"),
  file.path(integrated_dir, "16_tumor_signature_scores.pdf")
)

plot_tumor_gsea_focus(
  all_gsea,
  file.path(integrated_dir, "17_tumor_gsea_focus.png"),
  file.path(integrated_dir, "17_tumor_gsea_focus.pdf")
)

signature_summary <- signature_scores |>
  group_by(dataset, signature, group) |>
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    genes_present = first(genes_present),
    .groups = "drop"
  ) |>
  mutate(group_key = ifelse(group == "FAP-deficient host", "fap_host", "wt_host")) |>
  select(dataset, signature, group_key, mean_score, sd_score, genes_present) |>
  pivot_wider(names_from = group_key, values_from = c(mean_score, sd_score), names_sep = "_") |>
  mutate(delta_fap_minus_wt = mean_score_fap_host - mean_score_wt_host)

top_rna_pathways <- all_gsea |>
  filter(dataset == "RNA") |>
  arrange(fdr, desc(abs(NES))) |>
  slice_head(n = 15)

top_prot_pathways <- all_gsea |>
  filter(dataset == "Proteomics") |>
  arrange(fdr, desc(abs(NES))) |>
  slice_head(n = 15)

write.table(
  filter(signature_scores, dataset == "RNA"),
  file.path(rna_dir, "rna_signature_scores_per_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  filter(tumor_signature_scores, dataset == "RNA"),
  file.path(rna_dir, "rna_tumor_signature_scores_per_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(rna_diff, file.path(rna_dir, "rna_differential_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filter(all_gsea, dataset == "RNA"), file.path(rna_dir, "rna_all_gsea_pathways.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top_rna_pathways, file.path(rna_dir, "top_rna_pathways.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(
  filter(signature_scores, dataset == "Proteomics"),
  file.path(proteomics_dir, "proteomics_signature_scores_per_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  filter(tumor_signature_scores, dataset == "Proteomics"),
  file.path(proteomics_dir, "proteomics_tumor_signature_scores_per_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(prot_diff, file.path(proteomics_dir, "proteomics_differential_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filter(all_gsea, dataset == "Proteomics"), file.path(proteomics_dir, "proteomics_all_gsea_pathways.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top_prot_pathways, file.path(proteomics_dir, "top_proteomics_pathways.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(signature_scores, file.path(integrated_dir, "signature_scores_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tumor_signature_scores, file.path(integrated_dir, "tumor_signature_scores_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(signature_summary, file.path(integrated_dir, "signature_scores_group_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(marker_changes, file.path(integrated_dir, "selected_marker_changes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_gsea, file.path(integrated_dir, "all_gsea_pathways.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

rna_stats <- rna_diff |>
  filter(significant == "yes") |>
  count(regulation) |>
  arrange(desc(regulation))

prot_stats <- prot_diff |>
  filter(significant == "yes") |>
  count(regulation) |>
  arrange(desc(regulation))

summary_lines <- c(
  "WT host vs FAP-deficient host xenograft bulk-omics summary",
  paste("Output directory:", outdir),
  paste("RNA outputs:", rna_dir),
  paste("Proteomics outputs:", proteomics_dir),
  paste("Integrated outputs:", integrated_dir),
  "",
  "Comparison direction:",
  "- Positive log2FC / positive NES: higher in FAP-deficient host xenografts.",
  "- Negative log2FC / negative NES: higher in WT host xenografts.",
  "",
  "RNA differential summary (vendor significant == yes):",
  paste(capture.output(print(rna_stats)), collapse = "\n"),
  "",
  "Proteomics differential summary (vendor significant == yes):",
  paste(capture.output(print(prot_stats)), collapse = "\n"),
  "",
  "Notes:",
  "- Fap is significantly lower in RNA from FAP-deficient host xenografts.",
  "- Pax8 is increased in RNA and shows a strong positive trend in proteomics.",
  "- Bulk GSEA and curated marker trends are not perfectly concordant for every stromal pathway; inspect both the pathway plots and the marker heatmap together.",
  "- Proteomics pathway FDR values are generally weaker than RNA, so interpret protein-level pathway shifts as supportive rather than definitive."
)

write_summary(file.path(outdir, "summary_notes.txt"), summary_lines)

message("Analysis complete. Outputs written to: ", outdir)
