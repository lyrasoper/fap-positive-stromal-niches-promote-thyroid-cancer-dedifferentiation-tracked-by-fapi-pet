# Cancer Research submission - figure code release
# Builds: Supplementary figures for Figure 1 (atlas QC and annotation).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(png)
  library(grid)
})

project_dir <- Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>")
input_rds <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/new_seurat_final_2026-02-05.rds")
caf_rds <- file.path(project_dir, "data/fig1_scRNA/caf_seurat_reclustered_res0.1_to_0.6.rds")

main_fig1_dir <- file.path(project_dir, "outputs/fig1")
main_scrna_dir <- file.path(project_dir, "outputs/fig1_scRNA")
supp_dir <- file.path(project_dir, "outputs/supplementary_figures")
supp_source_dir <- file.path(supp_dir, "source_data")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_source_dir, recursive = TRUE, showWarnings = FALSE)

celltype_cols <- c(
  "Fibro" = "#B0302A",
  "Thyro" = "#1F5AA6",
  "T/NK" = "#D76414",
  "Myeloid" = "#4B1F7A",
  "Endo" = "#008C9E",
  "B" = "#B88900"
)

tissue_cols <- c(
  "NT" = "#E78B8B",
  "PTC" = "#2FBF71",
  "ATC" = "#7BAAF7"
)

caf_cols <- c(
  "FAP+ infCAF" = "#D76414",
  "ecmCAF" = "#B0302A",
  "EndMT CAF" = "#008C9E",
  "RGS15+ myoCAF" = "#4B1F7A",
  "adiCAF" = "#B88900"
)

tissue_levels <- c("NT", "PTC", "ATC")

save_combo <- function(plot, stem, width, height) {
  ggsave(
    filename = file.path(supp_dir, paste0(stem, ".pdf")),
    plot = plot,
    device = "pdf",
    width = width,
    height = height,
    units = "in"
  )
  ggsave(
    filename = file.path(supp_dir, paste0(stem, ".png")),
    plot = plot,
    device = "png",
    width = width,
    height = height,
    units = "in",
    dpi = 600
  )
}

image_panel <- function(path, left = 0, right = 1, top = 0, bottom = 1) {
  # Panel PNGs are produced by upstream main-figure scripts (see README).
  if (!file.exists(path)) stop(sprintf("Missing %s — produced by an upstream main-figure step; see README.", path))
  img <- png::readPNG(path)
  h <- dim(img)[1]
  w <- dim(img)[2]
  row_idx <- seq.int(max(1, floor(top * h) + 1), min(h, ceiling(bottom * h)))
  col_idx <- seq.int(max(1, floor(left * w) + 1), min(w, ceiling(right * w)))
  cropped <- img[row_idx, col_idx, , drop = FALSE]
  patchwork::wrap_elements(full = grid::rasterGrob(cropped, interpolate = TRUE))
}

assign_common_meta <- function(obj) {
  obj$Celltype <- as.character(obj$Celltype)
  obj$Celltype[obj$Celltype == "Meyloid_cells"] <- "Myeloid_cells"
  obj$Celltype_short <- dplyr::recode(
    obj$Celltype,
    "Thyrocytes" = "Thyro",
    "Endo_cells" = "Endo",
    "Fibro_cells" = "Fibro",
    "Myeloid_cells" = "Myeloid",
    "T_NK" = "T/NK",
    "B_cells" = "B"
  )
  obj$Celltype_short <- factor(obj$Celltype_short, levels = c("Thyro", "T/NK", "Myeloid", "Fibro", "Endo", "B"))
  obj$TissueType <- dplyr::recode(obj$Tissue_label_detail8, "Normal" = "NT", "PTC" = "PTC", "ATC" = "ATC")
  obj$TissueType <- factor(obj$TissueType, levels = tissue_levels)
  obj
}

assign_caf_clusters <- function(caf) {
  cluster_to_caf <- c(
    "0" = "RGS15+ myoCAF",
    "1" = "FAP+ infCAF",
    "2" = "ecmCAF",
    "3" = "FAP+ infCAF",
    "4" = "adiCAF",
    "5" = "FAP+ infCAF",
    "6" = "RGS15+ myoCAF",
    "7" = "FAP+ infCAF",
    "8" = "EndMT CAF",
    "9" = "RGS15+ myoCAF",
    "10" = "EndMT CAF",
    "11" = "RGS15+ myoCAF",
    "12" = "RGS15+ myoCAF"
  )
  caf@meta.data$CAF_clusters <- factor(
    cluster_to_caf[as.character(caf$RNA_snn_res.0.1)],
    levels = c("FAP+ infCAF", "ecmCAF", "EndMT CAF", "RGS15+ myoCAF", "adiCAF")
  )
  caf@meta.data$TissueType <- factor(dplyr::recode(caf$Tissue_label_detail8, "Normal" = "NT", "PTC" = "PTC", "ATC" = "ATC"), levels = tissue_levels)
  caf
}

# caf_rds is produced by an upstream CAF reclustering step (see README).
if (!file.exists(input_rds)) stop(sprintf("Missing %s — provided as external single-cell atlas data; see README.", input_rds))
if (!file.exists(caf_rds)) stop(sprintf("Missing %s — produced by an upstream step; see README.", caf_rds))

obj <- assign_common_meta(readRDS(input_rds))
caf <- assign_caf_clusters(readRDS(caf_rds))

# Supplementary Fig. 1
supp1_a <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "TissueType",
  cols = tissue_cols,
  pt.size = 0.22,
  raster = TRUE
) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(face = "bold", color = "black"),
    legend.title = element_blank()
  )

study_cols <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")
supp1_b <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "Data.ID",
  cols = study_cols,
  pt.size = 0.22,
  raster = TRUE
) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(face = "bold", color = "black"),
    legend.title = element_blank()
  )

major_features <- c("TG", "TPO", "PAX8", "FAP", "COL1A1", "PECAM1", "KDR", "LST1", "C1QC", "CD3D", "NKG7", "MS4A1", "CD79A")
supp1_c <- DotPlot(
  obj,
  features = major_features,
  group.by = "Celltype_short",
  dot.scale = 6
) +
  RotatedAxis() +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )
write.csv(supp1_c$data, file.path(supp_source_dir, "Supplementary_Fig1C_major_lineage_dotplot_data.csv"), row.names = FALSE)

sample_fibro <- obj@meta.data %>%
  dplyr::count(Sample_id, TissueType, Celltype_short, name = "n") %>%
  group_by(Sample_id, TissueType) %>%
  mutate(total_cells = sum(n), fraction = n / total_cells) %>%
  ungroup() %>%
  filter(Celltype_short == "Fibro")
write.csv(sample_fibro, file.path(supp_source_dir, "Supplementary_Fig1D_sample_level_fibro_fraction.csv"), row.names = FALSE)

supp1_d <- ggplot(sample_fibro, aes(x = TissueType, y = fraction, fill = TissueType)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.16, size = 2, shape = 21, stroke = 0.7) +
  scale_fill_manual(values = tissue_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Fibroblast fraction") +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(face = "bold", color = "black"),
    legend.position = "none"
  )

sample_fap <- caf@meta.data %>%
  dplyr::count(Sample_id, TissueType, CAF_clusters, name = "n") %>%
  group_by(Sample_id, TissueType) %>%
  mutate(total_caf = sum(n), fraction = n / total_caf) %>%
  ungroup() %>%
  filter(CAF_clusters == "FAP+ infCAF")
write.csv(sample_fap, file.path(supp_source_dir, "Supplementary_Fig1E_sample_level_fap_infCAF_fraction.csv"), row.names = FALSE)

supp1_e <- ggplot(sample_fap, aes(x = TissueType, y = fraction, fill = TissueType)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.16, size = 2, shape = 21, stroke = 0.7) +
  scale_fill_manual(values = tissue_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "FAP+ infCAF fraction among CAFs") +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(face = "bold", color = "black"),
    legend.position = "none"
  )

caf_features <- c("FAP", "CXCL12", "IL6", "COL1A1", "POSTN", "PLVAP", "RAMP2", "RGS5", "ACTA2", "SFRP5", "OGN")
supp1_f <- DotPlot(
  caf,
  features = caf_features,
  group.by = "CAF_clusters",
  dot.scale = 6
) +
  RotatedAxis() +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold")
  )
write.csv(supp1_f$data, file.path(supp_source_dir, "Supplementary_Fig1F_caf_dotplot_data.csv"), row.names = FALSE)

design1 <- "
AB
CC
DE
FF
"
supp1 <- supp1_a + supp1_b + supp1_c + supp1_d + supp1_e + supp1_f +
  plot_layout(design = design1) +
  plot_annotation(
    title = "Supplementary Fig. 1 | Additional validation of the single-cell atlas and CAF annotation",
    tag_levels = "A"
  ) &
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    plot.tag = element_text(size = 18, face = "bold")
  )
save_combo(
  supp1,
  "Supplementary_Fig_1_Additional_validation_of_the_single_cell_atlas_and_CAF_annotation",
  12,
  16
)

# Supplementary Fig. 2
s2_a <- image_panel(file.path(main_fig1_dir, "FIG1I_FAP_TDSgenes_Spearman_bubble_recovered.png"), left = 0.03, right = 0.98, top = 0.02, bottom = 0.98)
s2_b <- image_panel(file.path(main_fig1_dir, "Fig1k_FAP_vs_TG_recovered.png"), left = 0.02, right = 0.98, top = 0.02, bottom = 0.98)
s2_c <- image_panel(file.path(main_fig1_dir, "Fig1L_FAP_vs_TPO_recovered.png"), left = 0.02, right = 0.98, top = 0.02, bottom = 0.98)
s2_d <- image_panel(file.path(main_fig1_dir, "Fig1M_FAP_vs_PAX8_recovered.png"), left = 0.02, right = 0.98, top = 0.02, bottom = 0.98)

supp2 <- (s2_a | s2_b) / (s2_c | s2_d) +
  plot_annotation(
    title = "Supplementary Fig. 2 | Extended TCGA-THCA validation",
    tag_levels = "A"
  ) &
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    plot.tag = element_text(size = 18, face = "bold")
  )
save_combo(
  supp2,
  "Supplementary_Fig_2_Extended_TCGA_THCA_validation",
  12,
  10
)

bubble_src <- file.path(main_fig1_dir, "FIG1I_FAP_TDSgenes_Spearman_results.csv")
if (!file.exists(bubble_src)) stop(sprintf("Missing %s — produced by an upstream main-figure step; see README.", bubble_src))
file.copy(
  bubble_src,
  file.path(supp_source_dir, "Supplementary_Fig2A_bubble_source_data.csv"),
  overwrite = TRUE
)

# Supplementary Fig. 3
# Split the recovered two-panel ssGSEA plot into clean left/right panels so the
# original main-figure letters (O/P) do not carry over into the supplementary layout.
s3_a <- image_panel(
  file.path(main_fig1_dir, "FIG1o_ssGSEA_TDS_vs_FAPCAF_boxplot_recovered_gseapy.png"),
  left = 0.08, right = 0.49, top = 0.06, bottom = 0.98
)
s3_b <- image_panel(
  file.path(main_fig1_dir, "FIG1o_ssGSEA_TDS_vs_FAPCAF_boxplot_recovered_gseapy.png"),
  left = 0.56, right = 0.98, top = 0.06, bottom = 0.98
)
s3_c <- image_panel(
  file.path(main_fig1_dir, "FIG1p_heatmap_TDS_vs_FAPECM_recovered.png"),
  left = 0.02, right = 0.98, top = 0.02, bottom = 0.98
)

supp3 <- ((s3_a | s3_b) / s3_c) +
  plot_layout(heights = c(1, 1.55)) +
  plot_annotation(
    title = "Supplementary Fig. 3 | Extended GEO meta-cohort validation",
    tag_levels = "A"
  ) &
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    plot.tag = element_text(size = 18, face = "bold")
  )
save_combo(
  supp3,
  "Supplementary_Fig_3_Extended_GEO_meta_cohort_validation",
  13,
  13
)

message("Supplementary Figs. 1-3 written to: ", supp_dir)
