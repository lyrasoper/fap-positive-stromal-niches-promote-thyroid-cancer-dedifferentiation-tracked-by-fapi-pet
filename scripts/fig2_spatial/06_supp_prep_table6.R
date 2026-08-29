# Cancer Research submission - figure code release
# Builds: Supplementary figures/tables for Figure 2.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(Matrix)
})

project_dir <- Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>")
spatial_rds <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "spatial_visium/Spatial_integrated9P.rds")
supp_fig_dir <- file.path(project_dir, "outputs/supplementary_figures")
supp_source_dir <- file.path(supp_fig_dir, "source_data")
supp_table_dir <- file.path(project_dir, "outputs/supplementary_tables")

dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_table_dir, recursive = TRUE, showWarnings = FALSE)

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

add_sheet <- function(wb, sheet_name, df) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df)
  hdr_style <- createStyle(textDecoration = "bold", halign = "center", wrapText = TRUE)
  addStyle(
    wb,
    sheet = sheet_name,
    style = hdr_style,
    rows = 1,
    cols = seq_len(ncol(df)),
    gridExpand = TRUE
  )
  freezePane(wb, sheet = sheet_name, firstActiveRow = 2)
  setColWidths(wb, sheet = sheet_name, cols = seq_len(ncol(df)), widths = "auto")
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

message("Loading spatial object: ", spatial_rds)
if (!file.exists(spatial_rds)) {
  stop(sprintf("Missing %s — integrated spatial Visium object produced by an upstream step; see README.", spatial_rds))
}
obj <- readRDS(spatial_rds)

suffix <- sub(".*_(\\d+)$", "\\1", colnames(obj))
idx_to_sample <- setNames(names(obj@images), seq_along(names(obj@images)))
sample_id <- unname(idx_to_sample[suffix])
obj$sample_id <- sample_id

main_fig2_sections <- c("P32", "P12", "P98", "P44", "P26", "P83", "P57", "P17")

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
  Senescence1 = md$Senescence1,
  HGF_MET1 = md$HGF_MET1,
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
    "Retained in the integrated 9-section object but excluded from main Fig. 2 plotting in the recovered workflow.",
    "Displayed in main Fig. 2."
  )
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
    median_Senescence1 = median(Senescence1, na.rm = TRUE),
    median_HGF_MET1 = median(HGF_MET1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    included_in_main_Fig2_display = ifelse(as.character(sample_id) %in% main_fig2_sections, "Yes", "No")
  )

sample_order_by_tds <- section_summary %>%
  arrange(desc(median_TDS_score)) %>%
  pull(sample_id) %>%
  as.character()

mat <- GetAssayData(obj, assay = "SCT", layer = "data")
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

cluster_order <- cluster_module_scores$seurat_clusters

section_cluster_fraction <- md %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  count(sample_id, seurat_clusters, name = "spot_n") %>%
  group_by(sample_id) %>%
  mutate(cluster_fraction = spot_n / sum(spot_n)) %>%
  ungroup() %>%
  left_join(
    cluster_module_scores %>%
      select(seurat_clusters, dominant_class, cluster_label),
    by = "seurat_clusters"
  ) %>%
  mutate(
    sample_id = factor(sample_id, levels = sample_order_by_tds),
    seurat_clusters = factor(seurat_clusters, levels = cluster_order),
    dominant_class = factor(dominant_class, levels = broad_class_levels)
  ) %>%
  arrange(dominant_class, seurat_clusters, sample_id)

section_metric_long <- section_summary %>%
  select(
    sample_id,
    spot_count_total,
    median_nCount_Spatial,
    median_nFeature_Spatial,
    median_TDS_score,
    mean_FAP,
    FAP_positive_fraction,
    median_Senescence1,
    median_HGF_MET1
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
  "FAP_positive_fraction" = "FAP-positive fraction",
  "median_Senescence1" = "Median senescence score",
  "median_HGF_MET1" = "Median HGF-MET score"
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
    id = factor(id, levels = cluster_order),
    features.plot = factor(features.plot, levels = rev(unique(features.plot)))
  )

for (old in Sys.glob(file.path(supp_source_dir, "Supplementary_Fig4*.csv"))) {
  file.remove(old)
}

write.csv(
  umap_df,
  file.path(supp_source_dir, "Supplementary_Fig4AB_umap_spot_metadata.csv"),
  row.names = FALSE
)
write.csv(
  dot_data,
  file.path(supp_source_dir, "Supplementary_Fig4C_marker_dotplot_source.csv"),
  row.names = FALSE
)
write.csv(
  section_cluster_fraction,
  file.path(supp_source_dir, "Supplementary_Fig4D_section_cluster_fraction.csv"),
  row.names = FALSE
)
write.csv(
  section_metric_long,
  file.path(supp_source_dir, "Supplementary_Fig4E_section_summary_metrics.csv"),
  row.names = FALSE
)
write.csv(
  cluster_module_scores,
  file.path(supp_source_dir, "Supplementary_Fig4_cluster_module_scores.csv"),
  row.names = FALSE
)

p_umap_sample <- ggplot(
  umap_df,
  aes(x = UMAP_1, y = UMAP_2, color = sample_id)
) +
  geom_point(size = 0.18, alpha = 0.7) +
  scale_color_manual(values = sample_palette[names(obj@images)], drop = FALSE) +
  coord_equal() +
  labs(title = "Section identity", color = "Section") +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(face = "bold", color = "black"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size = 12.5, face = "bold")
  )

p_umap_cluster <- ggplot(
  umap_df,
  aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)
) +
  geom_point(size = 0.18, alpha = 0.75) +
  coord_equal() +
  labs(title = "Unsupervised spot clusters", color = "Cluster") +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(face = "bold", color = "black"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 12.5, face = "bold")
  )

p_dot <- ggplot(dot_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
  scale_size(range = c(0.2, 5.8)) +
  labs(
    title = "Canonical marker landscape across spot clusters",
    x = "Spot cluster",
    y = NULL,
    size = "% expressed",
    color = "Scaled\naverage"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    plot.title = element_text(size = 12.5, face = "bold"),
    strip.background = element_blank(),
    strip.text.y = element_text(face = "bold")
  ) +
  facet_grid(features.plot ~ ., scales = "free_y", space = "free_y")

p_comp <- ggplot(
  section_cluster_fraction,
  aes(x = sample_id, y = cluster_label, fill = cluster_fraction)
) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradientn(
    colors = c("white", "#A5D8FF", "#4D96FF", "#124076"),
    limits = c(0, max(section_cluster_fraction$cluster_fraction))
  ) +
  facet_grid(dominant_class ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Section-wise cluster composition",
    x = "Section (ordered by median TDS)",
    y = NULL,
    fill = "Spot\nfraction"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(size = 12.5, face = "bold"),
    strip.background = element_blank(),
    strip.text.y = element_text(face = "bold"),
    panel.spacing.y = unit(0.15, "lines")
  )

p_metric <- ggplot(
  section_metric_long,
  aes(x = sample_id, y = factor(metric, levels = rev(names(metric_labels))), fill = z_score)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_y_discrete(labels = rev(metric_labels[names(metric_labels)])) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(
    title = "Section-level summary metrics",
    x = "Section (ordered by median TDS)",
    y = NULL,
    fill = "Z-score"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 12.5, face = "bold")
  )

supp4 <- (p_umap_sample | p_umap_cluster) / p_dot / (p_comp | p_metric) +
  plot_layout(heights = c(0.95, 1.35, 1.15), widths = c(1.55, 1)) +
  plot_annotation(
    title = "Supplementary Fig. 4 | Spatial transcriptomics cohort overview and marker-guided spot-cluster architecture",
    tag_levels = "A"
  ) &
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    plot.tag = element_text(size = 18, face = "bold")
  )

save_combo(
  supp4,
  "Supplementary_Fig_4_Spatial_transcriptomics_cohort_overview_and_marker_guided_spot_cluster_architecture",
  14,
  18
)

table6_path <- file.path(
  supp_table_dir,
  "Supplementary_Table_6_Spatial_transcriptomics_cohort_metadata_and_spot_cluster_summary.xlsx"
)

selection_note <- data.frame(
  Item = c(
    "Integrated object",
    "Main Fig. 2 display",
    "Reason P1 is absent from main Fig. 2"
  ),
  Detail = c(
    "The recovered Seurat object contains 9 sections (P1, P32, P12, P98, P44, P26, P83, P57, P17) and 40,544 spots.",
    "The recovered main-Fig. 2 plotting workflow uses 8 sections: P32, P12, P98, P44, P26, P83, P57, and P17.",
    "Recovered server-side plotting code explicitly removed P1 before generating the final display panels."
  ),
  stringsAsFactors = FALSE
)

manifest <- data.frame(
  File = c(
    "Supplementary_Fig_4_Spatial_transcriptomics_cohort_overview_and_marker_guided_spot_cluster_architecture.png",
    "Supplementary_Fig_4_Spatial_transcriptomics_cohort_overview_and_marker_guided_spot_cluster_architecture.pdf",
    "Supplementary_Fig4AB_umap_spot_metadata.csv",
    "Supplementary_Fig4C_marker_dotplot_source.csv",
    "Supplementary_Fig4D_section_cluster_fraction.csv",
    "Supplementary_Fig4E_section_summary_metrics.csv",
    "Supplementary_Fig4_cluster_module_scores.csv"
  ),
  Description = c(
    "Supplementary Fig. 4 raster figure.",
    "Supplementary Fig. 4 vector figure.",
    "Spot-level UMAP coordinates and per-spot metadata used for panels A-B.",
    "Dot-plot source data used for panel C.",
    "Per-section cluster fractions used for panel D.",
    "Section-level summary metrics used for panel E.",
    "Marker-guided module scores and dominant broad class per spot cluster."
  ),
  stringsAsFactors = FALSE
)

wb6 <- createWorkbook()
add_sheet(wb6, "section_metadata", section_metadata)
add_sheet(wb6, "section_summary_metrics", section_summary)
add_sheet(wb6, "cluster_module_scores", cluster_module_scores)
add_sheet(wb6, "section_cluster_fraction", section_cluster_fraction)
add_sheet(wb6, "selection_note", selection_note)
add_sheet(wb6, "source_file_manifest", manifest)
saveWorkbook(wb6, table6_path, overwrite = TRUE)

message("Supplementary Fig. 4 written to: ", supp_fig_dir)
message("Supplementary Table 6 written to: ", table6_path)
