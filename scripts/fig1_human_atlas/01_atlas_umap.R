# Cancer Research submission - figure code release
# Builds: Figure 1 (human single-cell atlas).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(MASS)
})

project_dir <- Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>")
input_rds <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/new_seurat_final_2026-02-05.rds")
gmt_file <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "spatial_visium/reference/h.all.v2026.1.Hs.symbols.gmt")
python_bin <- Sys.getenv("PYTHON_BIN", unset = file.path(project_dir, ".conda/envs/fig1o/bin/python"))
gsea_script <- file.path(project_dir, "scripts/fig1_human_atlas/02_gsea_fig1h_prerank.py")

outdir <- file.path(project_dir, "outputs/fig1_scRNA")
panel_dir <- file.path(outdir, "panels")
source_dir <- file.path(outdir, "source_data")
data_dir <- file.path(project_dir, "data/fig1_scRNA")
recovered_dir <- file.path(data_dir, "recovered_code_reference")
caf_rds <- file.path(data_dir, "caf_seurat_reclustered_res0.1_to_0.6.rds")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(recovered_dir, recursive = TRUE, showWarnings = FALSE)

celltype_cols <- c(
  "Fibro" = "#B0302A",
  "Thyro" = "#1F5AA6",
  "T/NK" = "#D76414",
  "Myeloid" = "#4B1F7A",
  "Endo" = "#008C9E",
  "B" = "#B88900"
)

caf_cols <- c(
  "FAP+ infCAF" = "#D76414",
  "ecmCAF" = "#B0302A",
  "EndMT CAF" = "#008C9E",
  "RGS15+ myoCAF" = "#4B1F7A",
  "adiCAF" = "#B88900"
)

tissue_levels <- c("NT", "PTC", "ATC")

save_plot <- function(plot, stem, width, height) {
  ggsave(
    filename = file.path(panel_dir, paste0(stem, ".pdf")),
    plot = plot,
    device = "pdf",
    width = width,
    height = height,
    units = "in"
  )
  ggsave(
    filename = file.path(panel_dir, paste0(stem, ".png")),
    plot = plot,
    device = "png",
    width = width,
    height = height,
    units = "in",
    dpi = 600
  )
}

plot_density_like_reference <- function(seu, feature, reduction = "umap", grid_n = 200,
                                        q_cap = 0.99, dens_cap = 0.995, sample_n = 50000,
                                        pt.size = 0.32, trans = "sqrt") {
  um <- as.data.frame(Embeddings(seu, reduction)[, 1:2, drop = FALSE])
  colnames(um) <- c("UMAP_1", "UMAP_2")
  expr <- FetchData(seu, vars = feature)[[1]]
  df <- cbind(um, value = expr)

  w <- pmin(df$value, quantile(df$value, q_cap, na.rm = TRUE))
  w <- pmax(w, 0)
  w[is.na(w)] <- 0
  prob <- w
  if (sum(prob) == 0) prob <- rep(1, length(prob))
  prob <- prob / sum(prob)

  set.seed(1)
  idx <- sample(seq_len(nrow(df)), size = min(sample_n, nrow(df)), replace = TRUE, prob = prob)
  kd <- MASS::kde2d(
    x = df$UMAP_1[idx],
    y = df$UMAP_2[idx],
    n = grid_n,
    lims = c(range(df$UMAP_1), range(df$UMAP_2))
  )

  ix <- findInterval(df$UMAP_1, kd$x, all.inside = TRUE)
  iy <- findInterval(df$UMAP_2, kd$y, all.inside = TRUE)
  df$density <- kd$z[cbind(ix, iy)]
  cap <- quantile(df$density, dens_cap, na.rm = TRUE)
  df$density2 <- pmin(df$density, cap)

  ggplot(df, aes(UMAP_1, UMAP_2, color = density2)) +
    geom_point(size = pt.size) +
    scale_color_gradientn(
      colors = c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"),
      name = "Density",
      trans = trans,
      labels = scales::number_format(accuracy = 0.01)
    ) +
    labs(title = feature, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 10) +
    theme(
      text = element_text(face = "bold", color = "black"),
      plot.title = element_text(hjust = 0)
    )
}

cluster_signature_scores <- function(avg_expr) {
  marker_sets <- list(
    "RGS15+ myoCAF" = c("RGS5", "ACTA2", "TAGLN", "MYLK", "CNN1", "DES"),
    "adiCAF" = c("SFRP5", "ANGPTL7", "CLDN1", "TENM2", "ALDH1A2", "IGF1", "OGN", "C7", "CFD"),
    "FAP+ infCAF" = c("FAP", "TGFB1", "SPP1", "HGF", "TNC", "CXCL12", "CXCL14", "IL6", "CCL2", "MMP2", "MMP9", "DPP4"),
    "ecmCAF" = c("COL1A1", "COL3A1", "FN1", "POSTN", "COL5A1", "SPARC"),
    "EndMT CAF" = c("PLVAP", "RAMP2", "CLDN5", "FLT1", "RGCC", "RBP7", "PECAM1", "EMCN")
  )

  present <- intersect(unique(unlist(marker_sets)), rownames(avg_expr))
  z_mat <- t(scale(t(as.matrix(avg_expr[present, , drop = FALSE]))))
  z_mat[is.na(z_mat)] <- 0

  score_mat <- sapply(names(marker_sets), function(label) {
    genes <- intersect(marker_sets[[label]], rownames(z_mat))
    colMeans(z_mat[genes, , drop = FALSE])
  })

  score_df <- as.data.frame(score_mat)
  score_df$cluster <- sub("^g", "", rownames(score_df))
  score_df$best_label <- apply(score_mat, 1, function(x) names(which.max(x))[1])
  score_df
}

recluster_caf <- function(obj) {
  caf <- subset(obj, subset = Celltype == "Fibro_cells")
  caf[["percent.mt"]] <- PercentageFeatureSet(caf, pattern = "^MT-")
  set.seed(1)
  caf <- ScaleData(caf, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
  caf <- FindVariableFeatures(caf, nfeatures = 4000, verbose = FALSE)
  caf <- RunPCA(caf, npcs = 50, verbose = FALSE)
  caf <- FindNeighbors(caf, reduction = "pca", dims = 1:50, verbose = FALSE)
  caf <- FindClusters(caf, resolution = seq(0.1, 0.6, 0.1), verbose = FALSE)
  caf <- RunUMAP(caf, reduction = "pca", dims = 1:50, seed.use = 1, verbose = FALSE)
  caf
}

plot_panel_a <- function(meta) {
  counts <- meta %>%
    distinct(Sample_id, TissueType) %>%
    count(TissueType, name = "n") %>%
    complete(TissueType = tissue_levels, fill = list(n = 0))

  labels <- c(
    sprintf("NT (n = %s)", counts$n[counts$TissueType == "NT"]),
    sprintf("PTC (n = %s)", counts$n[counts$TissueType == "PTC"]),
    sprintf("ATC (n = %s)", counts$n[counts$TissueType == "ATC"])
  )

  ggplot() +
    annotate("text", x = 0, y = 3.35, label = "Public TC\nsingle-cell cohort", hjust = 0, size = 5.2, fontface = "bold") +
    annotate("text", x = 0.15, y = c(2.3, 1.55, 0.8), label = labels, hjust = 0, size = 6.3, fontface = "bold") +
    coord_cartesian(xlim = c(0, 3.3), ylim = c(0.4, 3.7), clip = "off") +
    theme_void()
}

obj <- readRDS(input_rds)
obj$Celltype <- as.character(obj$Celltype)
obj$Celltype[obj$Celltype == "Meyloid_cells"] <- "Myeloid_cells"
obj$Celltype_short <- recode(
  obj$Celltype,
  "Thyrocytes" = "Thyro",
  "Endo_cells" = "Endo",
  "Fibro_cells" = "Fibro",
  "Myeloid_cells" = "Myeloid",
  "T_NK" = "T/NK",
  "B_cells" = "B"
)
obj$Celltype_short <- factor(obj$Celltype_short, levels = c("Thyro", "T/NK", "Myeloid", "Fibro", "Endo", "B"))
obj$TissueType <- recode(obj$Tissue_label_detail8, "Normal" = "NT", "PTC" = "PTC", "ATC" = "ATC")
obj$TissueType <- factor(obj$TissueType, levels = tissue_levels)

meta_full <- obj@meta.data %>%
  tibble::rownames_to_column("cell_barcode")
write.csv(meta_full, file.path(source_dir, "fig1_scrna_full_metadata.csv"), row.names = FALSE)

panel_a <- plot_panel_a(obj@meta.data)
save_plot(panel_a, "Fig1A_cohort_overview", 3.2, 2.4)

panel_b <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "Celltype_short",
  cols = celltype_cols[c("Thyro", "T/NK", "Myeloid", "Fibro", "Endo", "B")],
  pt.size = 0.22,
  alpha = 1,
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  raster = TRUE
) +
  labs(title = NULL) +
  theme_classic(base_size = 12) +
  theme(text = element_text(face = "bold", color = "black"))
save_plot(panel_b, "Fig1B_major_lineage_umap", 4.6, 3.7)

df_prop <- obj@meta.data %>%
  count(TissueType, Celltype_short, name = "n") %>%
  complete(TissueType = tissue_levels, Celltype_short = names(celltype_cols), fill = list(n = 0)) %>%
  mutate(
    TissueType = factor(TissueType, levels = tissue_levels),
    Celltype_short = factor(Celltype_short, levels = c("Fibro", "Thyro", "T/NK", "Myeloid", "Endo", "B"))
  ) %>%
  group_by(TissueType) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
write.csv(df_prop, file.path(source_dir, "Fig1C_celltype_proportions.csv"), row.names = FALSE)

fibro_df <- df_prop %>%
  filter(Celltype_short == "Fibro") %>%
  mutate(label = percent(prop, accuracy = 0.1), y_text = pmin(prop + 0.08, 0.98))

panel_c <- ggplot(df_prop, aes(x = TissueType, y = prop, fill = Celltype_short)) +
  geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
  geom_errorbar(
    data = fibro_df,
    aes(x = TissueType, ymin = prop, ymax = prop),
    inherit.aes = FALSE,
    width = 0.82,
    linewidth = 1.4,
    color = "white"
  ) +
  geom_errorbar(
    data = fibro_df,
    aes(x = TissueType, ymin = prop, ymax = prop),
    inherit.aes = FALSE,
    width = 0.82,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_line(
    data = fibro_df,
    aes(x = TissueType, group = 1, y = prop),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_point(
    data = fibro_df,
    aes(x = TissueType, y = prop),
    inherit.aes = FALSE,
    shape = 21,
    size = 3,
    stroke = 0.8,
    fill = "black",
    color = "black"
  ) +
  geom_text(
    data = fibro_df,
    aes(x = TissueType, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 3.3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = celltype_cols[c("Fibro", "Thyro", "T/NK", "Myeloid", "Endo", "B")], breaks = c("Fibro", "Thyro", "T/NK", "Myeloid", "Endo", "B"), drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0.1)) +
  labs(title = "Cell Type Proportion", x = NULL, y = "Proportion") +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    plot.margin = margin(5.5, 16, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")
save_plot(panel_c, "Fig1C_celltype_proportion", 3.6, 3.2)

panel_d <- DotPlot(
  obj,
  features = c("FAP", "TDS_Score1"),
  group.by = "TissueType",
  dot.scale = 10,
  scale.by = "size",
  col.min = -1,
  col.max = 0.5
) +
  scale_color_gradientn(
    colors = c("#330066", "#336699", "#66CC66", "#FFCC33"),
    limits = c(-1, 0.5),
    oob = scales::squish
  ) +
  scale_size(range = c(2.5, 10)) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  )
panel_d$data$features.plot <- factor(
  recode(panel_d$data$features.plot, "FAP" = "FAP", "TDS_Score1" = "TDS Score"),
  levels = c("FAP", "TDS Score")
)
save_plot(panel_d, "Fig1D_fap_tds_dotplot", 3.3, 3.2)

caf_seurat <- if (file.exists(caf_rds)) {
  readRDS(caf_rds)
} else {
  tmp <- recluster_caf(obj)
  saveRDS(tmp, caf_rds)
  tmp
}

avg_expr <- AggregateExpression(caf_seurat, group.by = "RNA_snn_res.0.1", assays = "RNA", slot = "data", verbose = FALSE)$RNA
score_df <- cluster_signature_scores(avg_expr)
write.csv(score_df, file.path(source_dir, "Fig1_caf_cluster_signature_scores.csv"), row.names = FALSE)

# Manual subtype mapping derived from recovered code and current marker signatures.
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

caf_seurat@meta.data$CAF_clusters <- factor(
  cluster_to_caf[as.character(caf_seurat$RNA_snn_res.0.1)],
  levels = c("ecmCAF", "FAP+ infCAF", "EndMT CAF", "RGS15+ myoCAF", "adiCAF")
)
caf_seurat@meta.data$CAF_clusters_prop <- factor(
  cluster_to_caf[as.character(caf_seurat$RNA_snn_res.0.1)],
  levels = c("FAP+ infCAF", "ecmCAF", "EndMT CAF", "RGS15+ myoCAF", "adiCAF")
)
write.csv(
  caf_seurat@meta.data %>% tibble::rownames_to_column("cell_barcode"),
  file.path(source_dir, "Fig1_caf_reclustered_metadata.csv"),
  row.names = FALSE
)

panel_e <- DimPlot(
  caf_seurat,
  reduction = "umap",
  group.by = "CAF_clusters",
  cols = caf_cols[c("ecmCAF", "FAP+ infCAF", "EndMT CAF", "RGS15+ myoCAF", "adiCAF")],
  pt.size = 0.55,
  alpha = 1,
  label = TRUE,
  repel = TRUE,
  label.size = 3.8,
  raster = TRUE
) +
  labs(title = "CAF Subtypes") +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
save_plot(panel_e, "Fig1E_caf_umap", 4.5, 3.2)

panel_f <- plot_density_like_reference(caf_seurat, "FAP") +
  theme(
    text = element_text(face = "bold", color = "black"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8)
  )
save_plot(panel_f, "Fig1F_fap_density", 3.6, 3.2)

df_caf_prop <- caf_seurat@meta.data %>%
  transmute(TissueType = recode(Tissue_label_detail8, "Normal" = "NT", "PTC" = "PTC", "ATC" = "ATC"),
            CAF_clusters = CAF_clusters_prop) %>%
  count(TissueType, CAF_clusters, name = "n") %>%
  complete(TissueType = tissue_levels, CAF_clusters = levels(caf_seurat$CAF_clusters_prop), fill = list(n = 0)) %>%
  mutate(TissueType = factor(TissueType, levels = tissue_levels)) %>%
  group_by(TissueType) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
write.csv(df_caf_prop, file.path(source_dir, "Fig1G_caf_proportions.csv"), row.names = FALSE)

fap_df <- df_caf_prop %>%
  filter(CAF_clusters == "FAP+ infCAF") %>%
  mutate(label = percent(prop, accuracy = 0.1))

panel_g <- ggplot(df_caf_prop, aes(x = TissueType, y = prop, fill = CAF_clusters)) +
  geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
  geom_line(
    data = fap_df,
    aes(x = TissueType, group = 1, y = prop),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_point(
    data = fap_df,
    aes(x = TissueType, y = prop),
    inherit.aes = FALSE,
    size = 2.4,
    color = "black"
  ) +
  geom_text(
    data = fap_df,
    aes(x = TissueType, y = prop, label = label),
    inherit.aes = FALSE,
    nudge_y = -0.07,
    size = 3.3,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = caf_cols[c("FAP+ infCAF", "ecmCAF", "EndMT CAF", "RGS15+ myoCAF", "adiCAF")], drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0.08)) +
  labs(title = NULL, x = NULL, y = "Proportion") +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )
save_plot(panel_g, "Fig1G_caf_subtype_proportion", 3.6, 3.2)

Idents(caf_seurat) <- "CAF_clusters"
deg_file <- file.path(source_dir, "Fig1H_deg_fap_infCAF_vs_other.csv")
rank_file <- file.path(source_dir, "Fig1H_prerank_gene_list.csv")
if (file.exists(deg_file) && file.exists(rank_file)) {
  deg_res <- read.csv(deg_file, check.names = FALSE)
  logfc_col <- if ("avg_log2FC" %in% colnames(deg_res)) "avg_log2FC" else "avg_logFC"
} else {
  deg_res <- FindMarkers(
    caf_seurat,
    ident.1 = "FAP+ infCAF",
    ident.2 = setdiff(levels(caf_seurat$CAF_clusters), "FAP+ infCAF"),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.1,
    return.thresh = 1
  )
  logfc_col <- if ("avg_log2FC" %in% colnames(deg_res)) "avg_log2FC" else "avg_logFC"
  deg_res$gene <- rownames(deg_res)
  deg_res <- deg_res %>% arrange(desc(.data[[logfc_col]]))
  write.csv(deg_res, deg_file, row.names = FALSE)
  write.csv(
    deg_res %>% dplyr::select(gene, score = all_of(logfc_col)),
    rank_file,
    row.names = FALSE
  )
}

fig1h_top_file <- file.path(source_dir, "Fig1H_top_positive_hallmark_pathways.csv")
if (!file.exists(fig1h_top_file) && file.exists(python_bin) && file.exists(gsea_script) && file.exists(gmt_file)) {
  system2(
    python_bin,
    args = c(
      shQuote(gsea_script),
      "--rnk", shQuote(file.path(source_dir, "Fig1H_prerank_gene_list.csv")),
      "--gmt", shQuote(gmt_file),
      "--outdir", shQuote(source_dir)
    )
  )
}

if (!file.exists(fig1h_top_file)) {
  warning("Fig1H GSEA results not found; run fig1_human_atlas/02_gsea_fig1h_prerank.py with PYTHON_BIN set. Using an empty placeholder for panel H.")
  panel_h <- ggplot() + theme_void() +
    annotate("text", x = 0, y = 0, label = "Fig 1H (GSEA)\nnot generated", size = 3)
} else {
plot_df_h <- read.csv(fig1h_top_file, check.names = FALSE)
plot_df_h$Pathway <- factor(plot_df_h$Pathway, levels = rev(plot_df_h$Pathway))
panel_h <- ggplot(plot_df_h, aes(x = NES, y = Pathway)) +
  geom_vline(xintercept = 0, linewidth = 0.8, color = "black") +
  geom_col(fill = "#E67E22", width = 0.75) +
  scale_y_discrete(position = "right") +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(face = "bold", color = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y.right = element_text(size = 9, face = "bold"),
    axis.text.y.left = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(5.5, 2, 5.5, 5.5)
  ) +
  labs(title = "GSEA | FAP+ Fibroblast", x = "NES", y = NULL) +
  coord_cartesian(clip = "off")
save_plot(panel_h, "Fig1H_gsea_fap_infCAF", 5.2, 3.0)
}

panel_b_combined <- panel_b
panel_c_combined <- panel_c + theme(legend.position = "none")
panel_e_combined <- panel_e + theme(legend.position = "none")

design <- "
AABBC
DDEEF
GGHHH
"
fig_ah <- panel_a + panel_b_combined + panel_c_combined + panel_d + panel_e_combined + panel_f + panel_g + panel_h +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(
  filename = file.path(outdir, "Figure1_AH_scRNA_combined.pdf"),
  plot = fig_ah,
  device = "pdf",
  width = 13,
  height = 10.5,
  units = "in"
)
ggsave(
  filename = file.path(outdir, "Figure1_AH_scRNA_combined.png"),
  plot = fig_ah,
  device = "png",
  width = 13,
  height = 10.5,
  units = "in",
  dpi = 600
)

message("Figure 1 A-H outputs written to: ", outdir)
