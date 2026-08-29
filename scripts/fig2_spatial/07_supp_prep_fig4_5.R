#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary figures for Figure 2. The "fig4_5" in the filename
#          is this repository's internal supplement numbering.


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(ggrepel)
})

project_dir <- Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>")
spatial_rds <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "spatial_visium/Spatial_integrated9P.rds")
supp_fig_dir <- file.path(project_dir, "outputs", "supplementary_figures")
supp_source_dir <- file.path(supp_fig_dir, "source_data")
fig2_source_dir <- file.path(project_dir, "outputs", "figure2_rebuilt", "source_data")

dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_source_dir, recursive = TRUE, showWarnings = FALSE)

save_combo <- function(plot, stem, width, height) {
  ggsave(
    filename = file.path(supp_fig_dir, paste0(stem, ".pdf")),
    plot = plot,
    device = "pdf",
    width = width,
    height = height,
    units = "in"
  )
  ggsave(
    filename = file.path(supp_fig_dir, paste0(stem, ".png")),
    plot = plot,
    device = "png",
    width = width,
    height = height,
    units = "in",
    dpi = 600
  )
}

theme_supp <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(color = "black"),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

safe_get_counts <- function(obj, assay = "SCT", layer = "data") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

sample_palette <- c(
  "P1" = "#8C8C8C",
  "P12" = "#F28E2B",
  "P17" = "#E377C2",
  "P26" = "#59A14F",
  "P32" = "#00B894",
  "P44" = "#4E79A7",
  "P57" = "#5DA5DA",
  "P83" = "#8E63CE",
  "P98" = "#FF66B3"
)

broad_class_levels <- c(
  "Thyroid epithelial",
  "ECM/FAP stroma",
  "Myofibro/pericyte",
  "Myeloid",
  "Lymphoid",
  "Endothelial",
  "Mixed/low-signal"
)

main_fig2_sections <- c("P32", "P12", "P98", "P44", "P26", "P83", "P57", "P17")

message("Loading spatial object: ", spatial_rds)
if (!file.exists(spatial_rds)) {
  stop(sprintf("Missing %s — integrated spatial Visium object produced by an upstream step; see README.", spatial_rds))
}
obj <- readRDS(spatial_rds)

suffix <- sub(".*_(\\d+)$", "\\1", colnames(obj))
idx_to_sample <- setNames(names(obj@images), seq_along(names(obj@images)))
sample_id <- unname(idx_to_sample[suffix])
obj$sample_id <- sample_id

md <- obj@meta.data %>%
  mutate(
    spot = rownames(.),
    sample_id = factor(sample_id, levels = names(obj@images)),
    main_fig2_display = ifelse(sample_id %in% main_fig2_sections, "Yes", "No")
  )

emb <- Embeddings(obj, "umap")
umap_df <- data.frame(
  spot = rownames(emb),
  UMAP_1 = emb[, 1],
  UMAP_2 = emb[, 2],
  sample_id = md$sample_id,
  seurat_clusters = md$seurat_clusters,
  TDS_score = md$TDS_score,
  FAP = md$FAP,
  nCount_Spatial = md$nCount_Spatial,
  nFeature_Spatial = md$nFeature_Spatial,
  stringsAsFactors = FALSE
)

section_metadata <- tibble(
  sample_id = factor(names(obj@images), levels = names(obj@images)),
  spot_count_total = as.integer(table(md$sample_id)[names(obj@images)]),
  included_in_main_Fig2_display = ifelse(names(obj@images) %in% main_fig2_sections, "Yes", "No"),
  note = ifelse(
    names(obj@images) == "P1",
    "Retained in the integrated 9-section object but excluded from main Fig. 2 plotting.",
    "Displayed in main Fig. 2."
  )
)
write.csv(
  section_metadata,
  file.path(supp_source_dir, "Supplementary_Fig4_section_metadata.csv"),
  row.names = FALSE,
  quote = FALSE
)

section_summary <- md %>%
  group_by(sample_id) %>%
  summarise(
    spot_count_total = n(),
    median_nCount_Spatial = median(nCount_Spatial, na.rm = TRUE),
    median_nFeature_Spatial = median(nFeature_Spatial, na.rm = TRUE),
    median_TDS_score = median(TDS_score, na.rm = TRUE),
    mean_FAP = mean(FAP, na.rm = TRUE),
    FAP_positive_fraction = mean(FAP > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    included_in_main_Fig2_display = ifelse(as.character(sample_id) %in% main_fig2_sections, "Yes", "No")
  )

sample_order_by_tds <- section_summary %>%
  arrange(desc(median_TDS_score)) %>%
  pull(sample_id) %>%
  as.character()

mat <- safe_get_counts(obj, assay = "SCT", layer = "data")
cluster_factor <- obj$seurat_clusters

marker_sets <- list(
  "Thyroid epithelial" = c("EPCAM", "KRT19", "TG", "PAX8"),
  "ECM/FAP stroma" = c("FAP", "COL1A1", "COL1A2", "FN1", "POSTN", "CXCL12"),
  "Myofibro/pericyte" = c("ACTA2", "TAGLN", "RGS5", "PDGFRB"),
  "Myeloid" = c("PTPRC", "C1QC", "LST1"),
  "Lymphoid" = c("PTPRC", "CD3D", "NKG7", "MS4A1"),
  "Endothelial" = c("PECAM1", "KDR")
)
marker_sets <- lapply(marker_sets, intersect, rownames(mat))

cluster_levels <- levels(cluster_factor)
avg_expr <- sapply(cluster_levels, function(cl) Matrix::rowMeans(mat[, cluster_factor == cl, drop = FALSE]))
rownames(avg_expr) <- rownames(mat)
avg_expr_z <- t(scale(t(as.matrix(avg_expr))))
avg_expr_z[is.na(avg_expr_z)] <- 0
module_scores <- sapply(marker_sets, function(gs) colMeans(avg_expr_z[gs, , drop = FALSE]))

max_score <- apply(module_scores, 1, max)
second_score <- apply(module_scores, 1, function(x) sort(x, decreasing = TRUE)[2])
dominant_class <- apply(module_scores, 1, function(x) names(which.max(x))[1])
dominant_class[max_score < 0.4 | (max_score - second_score) < 0.12] <- "Mixed/low-signal"

cluster_module_scores <- data.frame(
  seurat_clusters = rownames(module_scores),
  module_scores,
  dominant_class = dominant_class,
  max_module_score = max_score,
  score_margin = max_score - second_score,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  mutate(
    dominant_class = factor(dominant_class, levels = broad_class_levels),
    cluster_num = as.numeric(as.character(seurat_clusters)),
    cluster_label = paste0("C", seurat_clusters)
  ) %>%
  arrange(dominant_class, cluster_num)

write.csv(
  cluster_module_scores,
  file.path(supp_source_dir, "Supplementary_Fig4_cluster_module_scores.csv"),
  row.names = FALSE,
  quote = FALSE
)

cluster_order <- cluster_module_scores$seurat_clusters

section_cluster_fraction <- md %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  count(sample_id, seurat_clusters, name = "spot_n") %>%
  group_by(sample_id) %>%
  mutate(cluster_fraction = spot_n / sum(spot_n)) %>%
  ungroup() %>%
  left_join(
    cluster_module_scores %>% select(seurat_clusters, dominant_class, cluster_label),
    by = "seurat_clusters"
  ) %>%
  mutate(
    sample_id = factor(sample_id, levels = sample_order_by_tds),
    seurat_clusters = factor(seurat_clusters, levels = cluster_order),
    dominant_class = factor(dominant_class, levels = broad_class_levels)
  ) %>%
  arrange(dominant_class, seurat_clusters, sample_id)

write.csv(
  section_cluster_fraction,
  file.path(supp_source_dir, "Supplementary_Fig4D_section_cluster_fraction.csv"),
  row.names = FALSE,
  quote = FALSE
)

section_metric_long <- section_summary %>%
  select(
    sample_id,
    spot_count_total,
    median_nCount_Spatial,
    median_nFeature_Spatial,
    median_TDS_score,
    mean_FAP,
    FAP_positive_fraction
  ) %>%
  mutate(sample_id = factor(sample_id, levels = sample_order_by_tds)) %>%
  pivot_longer(-sample_id, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  mutate(z_score = as.numeric(scale(value))) %>%
  ungroup()

metric_labels <- c(
  "spot_count_total" = "Spot count",
  "median_nCount_Spatial" = "Median UMI",
  "median_nFeature_Spatial" = "Median genes",
  "median_TDS_score" = "Median TDS score",
  "mean_FAP" = "Mean FAP",
  "FAP_positive_fraction" = "FAP-positive fraction"
)

write.csv(
  section_metric_long,
  file.path(supp_source_dir, "Supplementary_Fig4E_section_summary_metrics.csv"),
  row.names = FALSE,
  quote = FALSE
)

dot_features <- list(
  "Thyroid epithelial" = c("EPCAM", "KRT19", "TG", "PAX8"),
  "ECM/FAP stroma" = c("FAP", "COL1A1", "POSTN", "CXCL12"),
  "Myofibro/pericyte" = c("ACTA2", "TAGLN", "RGS5", "PDGFRB"),
  "Immune and endothelial" = c("PTPRC", "CD3D", "C1QC", "PECAM1", "KDR")
)

obj_dot <- obj
obj_dot@meta.data$FAP_meta <- obj_dot@meta.data$FAP
obj_dot@meta.data$FAP <- NULL
DefaultAssay(obj_dot) <- "SCT"

dot_data <- DotPlot(obj_dot, features = dot_features, group.by = "seurat_clusters")$data %>%
  mutate(
    features.plot = factor(features.plot, levels = rev(unique(features.plot))),
    id = factor(id, levels = cluster_order)
  )

write.csv(
  dot_data,
  file.path(supp_source_dir, "Supplementary_Fig4C_marker_dotplot_source.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  umap_df,
  file.path(supp_source_dir, "Supplementary_Fig4AB_umap_spot_metadata.csv"),
  row.names = FALSE,
  quote = FALSE
)

## Supplementary Fig. 4

p4_a <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = sample_id)) +
  geom_point(size = 0.12, alpha = 0.9) +
  scale_color_manual(values = sample_palette, drop = FALSE) +
  coord_equal() +
  labs(title = "Section identity", color = "Section") +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 11),
    legend.position = "right",
    legend.key.height = unit(0.35, "cm")
  )

p4_b <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = factor(seurat_clusters))) +
  geom_point(size = 0.12, alpha = 0.9) +
  scale_color_manual(values = scales::hue_pal()(length(unique(umap_df$seurat_clusters)))) +
  coord_equal() +
  labs(title = "Unsupervised spot clusters", color = "Cluster") +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 11),
    legend.position = "none"
  )

p4_c <- ggplot(dot_data, aes(id, features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0.3, 5.3), name = "% expressed") +
  scale_color_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, name = "Scaled\naverage") +
  labs(title = "Canonical marker landscape across spot clusters", x = "Spot cluster", y = NULL) +
  theme_supp(9) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

p4_d <- ggplot(section_cluster_fraction, aes(sample_id, cluster_label, fill = cluster_fraction)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2171B5", name = "Spot\nfraction") +
  labs(title = "Section-wise cluster composition", x = "Section (ordered by median TDS)", y = NULL) +
  theme_supp(9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7)
  )

p4_e <- ggplot(section_metric_long, aes(sample_id, factor(metric, levels = names(metric_labels)), fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#E15759", midpoint = 0, name = "Z-score") +
  scale_y_discrete(labels = metric_labels) +
  labs(title = "Section-level summary metrics", x = "Section (ordered by median TDS)", y = NULL) +
  theme_supp(9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  )

supp_fig4 <- (p4_a | p4_b) / p4_c / (p4_d | p4_e) +
  plot_annotation(
    title = "Supplementary Fig. 4 | Spatial transcriptomics cohort overview and marker-guided spot-cluster architecture",
    subtitle = "The integrated object contains 9 sections; P1 was retained in the object but excluded from the main Fig. 2 display workflow."
  )

save_combo(
  supp_fig4,
  "Supplementary_Fig_4_Spatial_transcriptomics_cohort_overview_and_marker_guided_spot_cluster_architecture",
  width = 12,
  height = 16
)

## Supplementary Fig. 5

spot_level_path <- file.path(fig2_source_dir, "Fig2_spatial_spot_level.csv")
if (!file.exists(spot_level_path)) {
  stop(sprintf("Missing %s — produced by an upstream Fig. 2 build step; see README.", spot_level_path))
}
spot_level <- read.csv(spot_level_path, stringsAsFactors = FALSE)
spot_level$slice <- factor(spot_level$slice, levels = main_fig2_sections)

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 6) {
    return(c(rho = NA_real_, p = NA_real_, n = sum(keep)))
  }
  ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(keep))
}

section_cor <- bind_rows(lapply(split(spot_level, spot_level$slice), function(x) {
  res <- safe_spearman(x$FAP_logCP10K, x$TDS16_z)
  data.frame(
    slice = as.character(x$slice[[1]]),
    rho = res[["rho"]],
    p = res[["p"]],
    n = res[["n"]],
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(
    slice = factor(slice, levels = main_fig2_sections),
    ci_low = pmax(rho - 1.96 / sqrt(pmax(n - 3, 1)), -1),
    ci_high = pmin(rho + 1.96 / sqrt(pmax(n - 3, 1)), 1)
  )

overall_cor <- safe_spearman(spot_level$FAP_logCP10K, spot_level$TDS16_z)

write.csv(
  section_cor,
  file.path(supp_source_dir, "Supplementary_Fig5A_per_section_spatial_correlations.csv"),
  row.names = FALSE,
  quote = FALSE
)

binned_by_sample_path <- file.path(fig2_source_dir, "Fig2_spatial_binned_by_sample.csv")
if (!file.exists(binned_by_sample_path)) {
  stop(sprintf("Missing %s — produced by an upstream Fig. 2 build step; see README.", binned_by_sample_path))
}
binned_by_sample <- read.csv(binned_by_sample_path, stringsAsFactors = FALSE) %>%
  mutate(slice = factor(slice, levels = main_fig2_sections))

write.csv(
  binned_by_sample,
  file.path(supp_source_dir, "Supplementary_Fig5B_spatial_binned_by_sample.csv"),
  row.names = FALSE,
  quote = FALSE
)

section_summary_fig2 <- spot_level %>%
  group_by(slice) %>%
  summarise(
    mean_FAP_logCP10K = mean(FAP_logCP10K, na.rm = TRUE),
    median_TDS16_z = median(TDS16_z, na.rm = TRUE),
    FAP_positive_fraction = mean(FAP_logCP10K > 0, na.rm = TRUE),
    n_spots = n(),
    .groups = "drop"
  ) %>%
  mutate(slice = factor(slice, levels = main_fig2_sections))

write.csv(
  section_summary_fig2,
  file.path(supp_source_dir, "Supplementary_Fig5C_section_level_spatial_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

p5_a <- ggplot(section_cor, aes(x = rho, y = slice, color = slice)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.15, linewidth = 0.55) +
  geom_point(size = 2.4) +
  geom_vline(xintercept = overall_cor[["rho"]], color = "black", linetype = 2, linewidth = 0.6) +
  scale_color_manual(values = sample_palette[names(sample_palette) %in% main_fig2_sections], drop = FALSE) +
  labs(
    title = "Per-section spot-level FAP-TDS associations",
    subtitle = paste0("Dashed line: overall rho = ", sprintf("%.3f", overall_cor[["rho"]])),
    x = "Spearman rho (FAP vs TDS)",
    y = NULL
  ) +
  theme_supp(9) +
  theme(legend.position = "none")

p5_b <- ggplot(binned_by_sample, aes(x = x_med, y = y_med, color = slice, group = slice)) +
  geom_path(linewidth = 0.6) +
  geom_point(size = 0.9) +
  facet_wrap(~slice, ncol = 4, scales = "free_x") +
  scale_color_manual(values = sample_palette[names(sample_palette) %in% main_fig2_sections], drop = FALSE) +
  labs(
    title = "Dose-response by section",
    subtitle = "Within-section medians across FAP bins (nbin = 10)",
    x = "Median FAP (logCP10K)",
    y = "Median TDS (z-score)"
  ) +
  theme_supp(8) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 7)
  )

section_level_cor <- safe_spearman(section_summary_fig2$mean_FAP_logCP10K, section_summary_fig2$median_TDS16_z)
p5_c <- ggplot(section_summary_fig2, aes(mean_FAP_logCP10K, median_TDS16_z, color = slice)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = 0.7) +
  ggrepel::geom_text_repel(aes(label = slice), size = 3.1, show.legend = FALSE, min.segment.length = 0) +
  scale_color_manual(values = sample_palette[names(sample_palette) %in% main_fig2_sections], drop = FALSE) +
  labs(
    title = "Section-level FAP-TDS summary",
    subtitle = paste0("rho = ", sprintf("%.2f", section_level_cor[["rho"]]), ", P = ", formatC(section_level_cor[["p"]], format = "e", digits = 2)),
    x = "Section mean FAP (logCP10K)",
    y = "Section median TDS (z-score)"
  ) +
  theme_supp(9) +
  theme(legend.position = "none")

supp_fig5 <- (p5_a | p5_b) / p5_c +
  plot_annotation(
    title = "Supplementary Fig. 5 | Section-level spatial support for the spatial analyses in Fig. 2",
    subtitle = "This figure extends the main text with per-section spatial FAP-TDS association summaries."
  )

save_combo(
  supp_fig5,
  "Supplementary_Fig_5_Section_level_spatial_support_for_the_spatial_analyses_in_Fig2",
  width = 14.5,
  height = 9.0
)

message("Saved Supplementary Fig. 4 and Supplementary Fig. 5 for Fig. 2 support.")
