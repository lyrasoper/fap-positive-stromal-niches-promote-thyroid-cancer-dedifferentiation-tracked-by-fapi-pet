# Cancer Research submission - figure code release
# Builds: Supplementary figures for Figure 1 (expanded atlas panels).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  library(patchwork)
  library(clustree)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "scripts/_shared/paper_A_style.R"))

project_dir <- Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>")
input_rds  <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/new_seurat_final_2026-02-05.rds")
caf_rds    <- file.path(project_dir, "data/fig1_scRNA/caf_seurat_reclustered_res0.1_to_0.6.rds")

supp_dir   <- file.path(project_dir, "outputs/supplementary_figures")
src_dir    <- file.path(supp_dir, "source_data")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir,  recursive = TRUE, showWarnings = FALSE)

tissue_levels <- c("NT", "PTC", "ATC")
tissue_cols   <- c("NT" = "#E78B8B", "PTC" = "#2FBF71", "ATC" = "#7BAAF7")
celltype_levels <- c("Thyro", "T/NK", "Myeloid", "Fibro", "Endo", "B")
study_cols    <- c("GSE148673" = "#1F77B4", "GSE184362" = "#FF7F0E",
                   "GSE193581" = "#2CA02C", "GSE210347" = "#D62728")
caf_cols <- c("FAP+ infCAF" = "#D76414", "ecmCAF" = "#B0302A",
              "EndMT CAF" = "#008C9E", "RGS15+ myoCAF" = "#4B1F7A",
              "adiCAF" = "#B88900")

common_theme <- theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(text = element_text(face = "plain", color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "plain"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "plain"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

save_combo <- function(plot, stem, width, height) {
  # Unicode-safe output: ragg PNG + quartz PDF. The base R pdf() device with the
  # built-in Helvetica AFM is Latin-1-only and mangles Greek (rho/sigma); cairo
  # is unavailable here (no XQuartz), so use ragg + macOS-native quartz instead.
  png_path <- file.path(supp_dir, paste0(stem, ".png"))
  ragg::agg_png(png_path, width = width, height = height, units = "in", res = 500)
  print(plot); dev.off()

  pdf_path <- file.path(supp_dir, paste0(stem, ".pdf"))
  if (capabilities("cairo")) grDevices::cairo_pdf(pdf_path, width = width, height = height) else if (capabilities("aqua")) grDevices::quartz(file = pdf_path, type = "pdf", width = width, height = height) else grDevices::pdf(pdf_path, width = width, height = height)
  print(plot); dev.off()
}

message("[1/5] Loading full scRNA object ...")
obj <- readRDS(input_rds)

obj$Celltype <- as.character(obj$Celltype)
obj$Celltype[obj$Celltype == "Meyloid_cells"] <- "Myeloid_cells"
obj$Celltype_short <- dplyr::recode(obj$Celltype,
  "Thyrocytes" = "Thyro", "Endo_cells" = "Endo", "Fibro_cells" = "Fibro",
  "Myeloid_cells" = "Myeloid", "T_NK" = "T/NK", "B_cells" = "B")
obj$Celltype_short <- factor(obj$Celltype_short, levels = celltype_levels)
obj$TissueType <- factor(dplyr::recode(obj$Tissue_label_detail8,
  "Normal" = "NT", "PTC" = "PTC", "ATC" = "ATC"), levels = tissue_levels)

if (!"percent.mt" %in% colnames(obj@meta.data)) {
  mt_pat <- "^MT-"; if (sum(grepl(mt_pat, rownames(obj))) < 3) mt_pat <- "^mt-"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pat)
}

message("[2/5] Panels A–C, I (full atlas) ...")
# A: per-tissue UMAPs (NT | PTC | ATC) — shows cluster-level changes across disease states
supp1_a <- DimPlot(obj, reduction = "umap", group.by = "TissueType",
                   split.by = "TissueType",
                   cols = tissue_cols, pt.size = 1.0, raster = TRUE,
                   raster.dpi = c(1200, 1200), ncol = 3) +
  common_theme +
  labs(title = "Per-tissue UMAPs (disease-state split)") +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_blank(),
        panel.spacing.x = unit(0.8, "lines"))

# B: per-cohort UMAPs — integration/batch check across the 4 public scRNA studies
supp1_b <- DimPlot(obj, reduction = "umap", group.by = "Data.ID",
                   split.by = "Data.ID",
                   cols = study_cols, pt.size = 1.0, raster = TRUE,
                   raster.dpi = c(1200, 1200), ncol = 4) +
  common_theme +
  labs(title = "Per-cohort UMAPs (source-study split, integration QC)") +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_blank(),
        panel.spacing.x = unit(0.8, "lines"))

major_features <- c("TG", "TPO", "PAX8", "FAP", "COL1A1", "PECAM1", "KDR",
                    "LST1", "C1QC", "CD3D", "NKG7", "MS4A1", "CD79A")
supp1_c <- DotPlot(obj, features = major_features, group.by = "Celltype_short",
                   dot.scale = 6,
                   cols = c("grey92", "#B0302A")) +
  RotatedAxis() + common_theme +
  labs(title = "Lineage marker expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())
write.csv(supp1_c$data, file.path(src_dir, "Supplementary_Fig1C_major_lineage_dotplot_data.csv"), row.names = FALSE)

# Panel I: FAP violin across 6 lineages
fap_vec <- FetchData(obj, vars = c("FAP", "Celltype_short", "TissueType"))
colnames(fap_vec)[1] <- "FAP"
write.csv(fap_vec, file.path(src_dir, "Supplementary_Fig1I_FAP_celltype_expression.csv"), row.names = FALSE)
supp1_i <- ggplot(fap_vec, aes(x = Celltype_short, y = FAP, fill = Celltype_short)) +
  geom_violin(scale = "width", trim = TRUE, width = 0.90, linewidth = 0.35, alpha = 0.88) +
  geom_boxplot(width = 0.10, outlier.shape = NA, fill = "white",
               linewidth = 0.35, colour = "grey20") +
  scale_fill_manual(values = c(Fibro = "#B0302A", Thyro = "#1F5AA6", `T/NK` = "#D76414",
                               Myeloid = "#4B1F7A", Endo = "#008C9E", B = "#B88900")) +
  labs(x = NULL, y = "FAP expression (log-norm)",
       title = "FAP across major lineages") +
  common_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1))

message("[3/5] Panels D, E, G (QC) + H (integration mixing) ...")
meta <- obj@meta.data %>% as_tibble(rownames = "cell")
sample_fibro <- meta %>%
  count(Sample_id, TissueType, Celltype_short, name = "n") %>%
  group_by(Sample_id, TissueType) %>% mutate(total = sum(n), fraction = n / total) %>%
  ungroup() %>% filter(Celltype_short == "Fibro")
write.csv(sample_fibro, file.path(src_dir, "Supplementary_Fig1D_sample_level_fibro_fraction.csv"), row.names = FALSE)
tissue_cmp <- list(c("NT","PTC"), c("PTC","ATC"), c("NT","ATC"))
supp1_d <- ggplot(sample_fibro, aes(TissueType, fraction, fill = TissueType)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.90, linewidth = 0.55) +
  geom_jitter(width = 0.16, size = 1.6, shape = 21, stroke = 0.35,
              colour = "black", alpha = 0.85) +
  stat_compare_means(comparisons = tissue_cmp, method = "wilcox.test",
                     label = "p.signif", size = 3.5, fontface = "bold",
                     tip.length = 0.01, bracket.size = 0.4) +
  scale_fill_manual(values = tissue_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Fibroblast fraction", title = "Per-sample fibroblast fraction") +
  common_theme + theme(legend.position = "none")

# Panel G: QC (nFeature, nCount, percent.mt) by TissueType
qc_df <- meta %>%
  select(Sample_id, TissueType, nFeature_RNA, nCount_RNA, percent.mt) %>%
  pivot_longer(cols = c(nFeature_RNA, nCount_RNA, percent.mt),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = dplyr::recode(metric,
          "nFeature_RNA" = "nFeature_RNA",
          "nCount_RNA"   = "nCount_RNA (log10)",
          "percent.mt"   = "percent.mt (%)"),
         value  = ifelse(metric == "nCount_RNA (log10)", log10(pmax(value, 1)), value))
write.csv(qc_df %>% group_by(Sample_id, TissueType, metric) %>%
            summarise(median = median(value), .groups = "drop"),
          file.path(src_dir, "Supplementary_Fig1G_sample_level_qc_summary.csv"), row.names = FALSE)
supp1_g <- ggplot(qc_df, aes(TissueType, value, fill = TissueType)) +
  geom_violin(scale = "width", trim = TRUE, width = 0.90, linewidth = 0.30, alpha = 0.88) +
  geom_boxplot(width = 0.10, outlier.shape = NA, fill = "white",
               linewidth = 0.30, colour = "grey20") +
  scale_fill_manual(values = tissue_cols) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL, title = "Per-cell QC by disease state") +
  common_theme +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_blank(),
        panel.spacing = unit(0.8, "lines"))

# Panel H: Data.ID mixing per celltype (batch-effect QC)
mix_df <- meta %>% count(Celltype_short, Data.ID) %>%
  group_by(Celltype_short) %>% mutate(prop = n / sum(n)) %>% ungroup()
write.csv(mix_df, file.path(src_dir, "Supplementary_Fig1H_celltype_cohort_mixing.csv"), row.names = FALSE)
supp1_h <- ggplot(mix_df, aes(Celltype_short, prop, fill = Data.ID)) +
  geom_col(width = 0.78, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = study_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Cohort composition",
       title = "Cohort mixing across lineages (integration QC)") +
  common_theme +
  theme(legend.title = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.35, "cm"),
        axis.text.x = element_text(angle = 25, hjust = 1))

rm(obj); gc(verbose = FALSE)

message("[4/5] CAF panels E, F, J, K, L ...")
caf <- readRDS(caf_rds)
cluster_to_caf <- c("0" = "RGS15+ myoCAF", "1" = "FAP+ infCAF", "2" = "ecmCAF",
                    "3" = "FAP+ infCAF", "4" = "adiCAF", "5" = "FAP+ infCAF",
                    "6" = "RGS15+ myoCAF", "7" = "FAP+ infCAF", "8" = "EndMT CAF",
                    "9" = "RGS15+ myoCAF", "10" = "EndMT CAF",
                    "11" = "RGS15+ myoCAF", "12" = "RGS15+ myoCAF")
caf@meta.data$CAF_clusters <- factor(cluster_to_caf[as.character(caf$RNA_snn_res.0.1)],
                           levels = c("FAP+ infCAF", "ecmCAF", "EndMT CAF",
                                      "RGS15+ myoCAF", "adiCAF"))
caf@meta.data$TissueType <- factor(dplyr::recode(caf$Tissue_label_detail8,
  "Normal" = "NT", "PTC" = "PTC", "ATC" = "ATC"), levels = tissue_levels)

sample_fap <- caf@meta.data %>%
  count(Sample_id, TissueType, CAF_clusters, name = "n") %>%
  group_by(Sample_id, TissueType) %>% mutate(total = sum(n), fraction = n / total) %>%
  ungroup() %>% filter(CAF_clusters == "FAP+ infCAF")
write.csv(sample_fap, file.path(src_dir, "Supplementary_Fig1E_sample_level_fap_infCAF_fraction.csv"), row.names = FALSE)
supp1_e <- ggplot(sample_fap, aes(TissueType, fraction, fill = TissueType)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.90, linewidth = 0.55) +
  geom_jitter(width = 0.16, size = 1.6, shape = 21, stroke = 0.35,
              colour = "black", alpha = 0.85) +
  stat_compare_means(comparisons = tissue_cmp, method = "wilcox.test",
                     label = "p.signif", size = 3.5, fontface = "bold",
                     tip.length = 0.01, bracket.size = 0.4) +
  scale_fill_manual(values = tissue_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "FAP+ infCAF fraction", title = "Per-sample FAP+ infCAF fraction") +
  common_theme + theme(legend.position = "none")

caf_features <- c("FAP", "CXCL12", "IL6", "COL1A1", "POSTN", "PLVAP",
                  "RAMP2", "RGS5", "ACTA2", "SFRP5", "OGN")
supp1_f <- DotPlot(caf, features = caf_features, group.by = "CAF_clusters",
                   dot.scale = 6, cols = c("grey92", "#4B1F7A")) +
  RotatedAxis() + common_theme +
  labs(title = "CAF subtype marker expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())
write.csv(supp1_f$data, file.path(src_dir, "Supplementary_Fig1F_caf_dotplot_data.csv"), row.names = FALSE)

# Panel J: clustree resolution stability
res_cols <- c("RNA_snn_res.0.1", "RNA_snn_res.0.2", "RNA_snn_res.0.3",
              "RNA_snn_res.0.4", "RNA_snn_res.0.5", "RNA_snn_res.0.6")
res_cols <- intersect(res_cols, colnames(caf@meta.data))
clust_meta <- caf@meta.data[, res_cols, drop = FALSE]
write.csv(clust_meta %>% tibble::rownames_to_column("cell"),
          file.path(src_dir, "Supplementary_Fig1J_caf_resolution_assignments.csv"), row.names = FALSE)
supp1_j <- clustree(clust_meta, prefix = "RNA_snn_res.",
                    node_colour = "sc3_stability", node_size_range = c(4, 10)) +
  labs(title = "CAF clustering resolution stability (res 0.1 \u2192 0.6)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

# Panel K: Per-sample CAF subtype composition
caf_sample <- caf@meta.data %>%
  count(Sample_id, TissueType, CAF_clusters) %>%
  group_by(Sample_id, TissueType) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  complete(tidyr::nesting(Sample_id, TissueType), CAF_clusters,
           fill = list(n = 0, prop = 0))
sample_order <- caf_sample %>% filter(CAF_clusters == "FAP+ infCAF") %>%
  arrange(TissueType, prop) %>% pull(Sample_id)
caf_sample$Sample_id <- factor(caf_sample$Sample_id, levels = sample_order)
write.csv(caf_sample, file.path(src_dir, "Supplementary_Fig1K_per_sample_caf_composition.csv"), row.names = FALSE)
supp1_k <- ggplot(caf_sample, aes(Sample_id, prop, fill = CAF_clusters)) +
  geom_col(width = 1, colour = "white", linewidth = 0.15) +
  facet_grid(~ TissueType, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = caf_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Samples (ordered by FAP+ infCAF proportion)",
       y = "CAF composition", fill = NULL,
       title = "Per-sample CAF subtype composition") +
  common_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        panel.spacing.x = unit(0.6, "lines"))

# Panel L: Top markers per CAF subtype (heatmap)
message("  finding CAF markers (this can take ~1-2 min) ...")
Idents(caf) <- caf$CAF_clusters
caf_markers <- FindAllMarkers(caf, only.pos = TRUE, min.pct = 0.15,
                              logfc.threshold = 0.25, verbose = FALSE)
top_markers <- caf_markers %>%
  group_by(cluster) %>% arrange(desc(avg_log2FC)) %>%
  slice_head(n = 8) %>% ungroup()
top_genes <- unique(top_markers$gene)
write.csv(top_markers, file.path(src_dir, "Supplementary_Fig1L_top_caf_markers.csv"), row.names = FALSE)

# averaged scaled expression per CAF subtype
avg_mat <- AverageExpression(caf, features = top_genes, group.by = "CAF_clusters",
                             assays = "RNA", slot = "data")[["RNA"]]
avg_mat_s <- t(scale(t(as.matrix(avg_mat))))
# order columns by CAF factor level
ord_cols <- c("FAP+ infCAF", "ecmCAF", "EndMT CAF", "RGS15+ myoCAF", "adiCAF")
avg_mat_s <- avg_mat_s[, intersect(ord_cols, colnames(avg_mat_s))]
# order rows by which cluster they mark (from top_markers)
row_order <- top_markers %>% filter(gene %in% rownames(avg_mat_s)) %>%
  mutate(cluster = factor(cluster, levels = ord_cols)) %>%
  arrange(cluster, desc(avg_log2FC)) %>% pull(gene) %>% unique()
avg_mat_s <- avg_mat_s[row_order, , drop = FALSE]

col_fun <- colorRamp2(c(-2, 0, 2), c("#2C5E9D", "white", "#C03840"))
row_anno <- rowAnnotation(
  Cluster = top_markers %>% filter(gene %in% row_order) %>%
    distinct(gene, .keep_all = TRUE) %>%
    arrange(match(gene, row_order)) %>% pull(cluster) %>% as.character(),
  col = list(Cluster = caf_cols),
  show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"))
ht_l <- Heatmap(avg_mat_s, name = "z-score",
                col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
                row_names_gp = gpar(fontface = "bold", fontsize = 7),
                column_names_gp = gpar(fontface = "bold", fontsize = 9),
                column_names_rot = 30, left_annotation = row_anno,
                heatmap_legend_param = list(direction = "vertical",
                                            legend_height = unit(3, "cm"),
                                            title_position = "topleft"),
                column_title = "Top markers across CAF subtypes",
                column_title_gp = gpar(fontface = "bold", fontsize = 11))
# render to grob so patchwork can embed
supp1_l <- patchwork::wrap_elements(full = grid.grabExpr(draw(ht_l,
            heatmap_legend_side = "right", annotation_legend_side = "right")))

message("[5/5] Assembling Supp Fig 1 ...")
# Six-row layout:
#   Row 1: A (3 per-tissue UMAPs) | B (4 per-cohort UMAPs) on same row
#   Row 2: I (FAP violin) | D | E (per-sample fractions)
#   Rows 3-6: C/G, F/H, K, J/L
design1 <- "
AAAAAAABBBBBBBB
IIIIIDDDDDEEEEE
CCCCCCGGGGGGGGG
FFFFFFHHHHHHHHH
KKKKKKKKKKKKKKK
JJJJJJLLLLLLLLL
"
supp1 <- supp1_a + supp1_b + supp1_c + supp1_d + supp1_e + supp1_f +
         supp1_g + supp1_h + supp1_i + supp1_j + supp1_k + supp1_l +
  plot_layout(design = design1,
              heights = c(1.0, 1.0, 1.0, 1.0, 1.1, 1.8)) +
  plot_annotation(
    title = "Supplementary Fig. 1 | Additional validation of the single-cell atlas and CAF annotation",
    tag_levels = "a"
  ) &
  theme(plot.title = element_text(size = 14, face = "bold", family = "Helvetica"),
        plot.tag   = element_text(size = 14, face = "bold", family = "Helvetica"))

save_combo(supp1, "Supplementary_Fig_1_Additional_validation_of_the_single_cell_atlas_and_CAF_annotation",
           width = 20, height = 26)

message("DONE Supp Fig 1")
