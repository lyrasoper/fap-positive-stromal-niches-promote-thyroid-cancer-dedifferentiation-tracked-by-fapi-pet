# Cancer Research submission - figure code release
# Builds: Figure 2 (spatial-transcriptomics panels only; the wet-lab
#          quantification panels of Figure 2 are not built here).

options(stringsAsFactors = FALSE)

required_pkgs <- c(
  "Seurat", "Matrix", "dplyr", "ggplot2", "patchwork",
  "readxl", "ggpubr", "rstatix"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    ". Please install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readxl)
  library(ggpubr)
  library(rstatix)
})

get_project_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    script_path <- gsub("~\\+~", " ", script_path)
    return(dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

save_plot <- function(plot_obj, stem, out_dir, width, height) {
  ggsave(
    filename = file.path(out_dir, paste0(stem, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
  ggsave(
    filename = file.path(out_dir, paste0(stem, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    bg = "white"
  )
}

theme_fig2 <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black", face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(color = "black"),
      plot.margin = margin(6, 6, 6, 6)
    )
}

theme_spatial_clean <- function() {
  theme_void(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 9, color = "black"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-300) return("<1e-300")
  if (p < 0.001) return(format(p, scientific = TRUE, digits = 2))
  sprintf("%.3f", p)
}

calc_logcp10k_gene <- function(obj, gene, assay = "Spatial") {
  count_layers <- grep("^counts", Layers(obj[[assay]]), value = TRUE)
  if (length(count_layers) == 0) stop("No counts.* layers found in assay: ", assay)

  out <- rep(NA_real_, ncol(obj))
  names(out) <- colnames(obj)

  for (lyr in count_layers) {
    m <- LayerData(obj[[assay]], layer = lyr)
    if (!gene %in% rownames(m)) next
    lib <- Matrix::colSums(m)
    lib[lib == 0] <- NA_real_
    v <- as.numeric(m[gene, ])
    out[colnames(m)] <- log1p((v / lib) * 1e4)
  }

  out
}

calc_tds16 <- function(obj, genes, assay = "Spatial") {
  count_layers <- grep("^counts", Layers(obj[[assay]]), value = TRUE)
  if (length(count_layers) == 0) stop("No counts.* layers found in assay: ", assay)

  tds <- rep(NA_real_, ncol(obj))
  names(tds) <- colnames(obj)
  matched_union <- character(0)

  for (lyr in count_layers) {
    m <- LayerData(obj[[assay]], layer = lyr)
    genes_use <- intersect(genes, rownames(m))
    matched_union <- union(matched_union, genes_use)
    if (length(genes_use) < 3) next

    lib <- Matrix::colSums(m)
    lib[lib == 0] <- NA_real_
    norm <- t(t(m) / lib) * 1e4
    norm <- log1p(norm)

    tds[colnames(m)] <- Matrix::colMeans(norm[genes_use, , drop = FALSE])
  }

  list(tds = tds, matched_union = matched_union)
}

project_dir <- get_project_dir()
data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "outputs", "figure2_rebuilt")
panel_dir <- file.path(output_dir, "panels")
source_dir <- file.path(output_dir, "source_data")
reference_dir <- file.path(output_dir, "reference")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)

spatial_rds_path <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "spatial_visium/Spatial_integrated9P.rds")

if (!file.exists(spatial_rds_path)) {
  stop(sprintf("Missing %s — integrated Visium object produced by an upstream step; see README.", spatial_rds_path))
}

message("Loading spatial object from: ", spatial_rds_path)
obj <- readRDS(spatial_rds_path)
DefaultAssay(obj) <- "Spatial"

sample_order <- c("P32", "P12", "P98", "P44", "P26", "P83", "P57", "P17")
available_images <- intersect(sample_order, names(obj@images))

spot2img <- rep(NA_character_, ncol(obj))
names(spot2img) <- colnames(obj)
for (im in available_images) {
  spot2img[Cells(obj[[im]])] <- im
}
obj$image_id <- spot2img

tds16_genes <- c(
  "TG", "TPO", "PAX8", "DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1",
  "GLIS3", "NKX2-1", "SLC26A4", "SLC5A5", "SLC5A8", "THRA", "THRB", "TSHR"
)

obj$FAP_logCP10K <- calc_logcp10k_gene(obj, "FAP", assay = "Spatial")
tds_res <- calc_tds16(obj, tds16_genes, assay = "Spatial")
obj$TDS16 <- tds_res$tds
obj$TDS16_z <- as.numeric(scale(obj$TDS16))

spatial_df <- obj@meta.data %>%
  mutate(
    spot = rownames(.),
    slice = factor(image_id, levels = sample_order)
  ) %>%
  filter(!is.na(slice), is.finite(FAP_logCP10K), is.finite(TDS16_z))

write.csv(
  spatial_df %>% select(spot, slice, FAP_logCP10K, TDS16, TDS16_z),
  file.path(source_dir, "Fig2_spatial_spot_level.csv"),
  row.names = FALSE
)

sample_palette <- c(
  "P12" = "#F28E2B",
  "P17" = "#E377C2",
  "P26" = "#59A14F",
  "P32" = "#00B894",
  "P44" = "#4E79A7",
  "P57" = "#5DA5DA",
  "P83" = "#8E63CE",
  "P98" = "#FF66B3"
)
sample_palette <- sample_palette[intersect(names(sample_palette), as.character(sample_order))]

add_scale_fap <- function(p) {
  p +
    scale_fill_gradientn(
      colors = c("#2C7BB6", "#FFFFBF", "#D7191C"),
      na.value = "grey90"
    ) +
    labs(fill = "FAP\n(logCP10K)") +
    theme_spatial_clean()
}

add_scale_tds <- function(p) {
  p +
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      na.value = "grey90",
      limits = c(-2, 2)
    ) +
    labs(fill = "TDS\n(z-score)") +
    theme_spatial_clean()
}

plist_fap <- SpatialFeaturePlot(
  obj,
  features = "FAP_logCP10K",
  images = available_images,
  alpha = c(0.15, 1),
  min.cutoff = "q2",
  max.cutoff = "q80",
  pt.size.factor = 6,
  combine = FALSE
)
plist_fap <- Map(function(p, nm) add_scale_fap(p + ggtitle(nm)), plist_fap, available_images)
fig2_b <- wrap_plots(plist_fap, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
save_plot(fig2_b, "FIG2C_FAP_logCP10K_rebuilt", panel_dir, width = 9, height = 13)

plist_tds <- SpatialFeaturePlot(
  obj,
  features = "TDS16_z",
  images = available_images,
  alpha = c(0.15, 1),
  min.cutoff = -2,
  max.cutoff = 2,
  pt.size.factor = 6,
  combine = FALSE
)
plist_tds <- Map(function(p, nm) add_scale_tds(p + ggtitle(nm)), plist_tds, available_images)
fig2_c <- wrap_plots(plist_tds, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
save_plot(fig2_c, "FIG2B_TDS16_z_rebuilt", panel_dir, width = 9, height = 13)

cor_res <- suppressWarnings(cor.test(
  spatial_df$FAP_logCP10K,
  spatial_df$TDS16_z,
  method = "spearman"
))

fig2_e_anno <- sprintf(
  "Spearman rho = %.3f\np = %s",
  unname(cor_res$estimate),
  format_p(cor_res$p.value)
)

fig2_e <- ggplot(spatial_df, aes(x = FAP_logCP10K, y = TDS16_z)) +
  geom_point(size = 0.25, alpha = 0.14, color = "grey55") +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.3, color = "black", span = 0.9) +
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = fig2_e_anno,
    hjust = 1.05,
    vjust = 1.08,
    size = 3.4,
    fill = "white",
    color = "black"
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "All spots",
    x = "FAP (logCP10K)",
    y = "TDS (z-score)"
  ) +
  theme_fig2(12)
save_plot(fig2_e, "FIG2E_FAP_vs_TDS_all_spots_rebuilt", panel_dir, width = 5.8, height = 5.2)

nbin <- 10
df_bin_slice <- spatial_df %>%
  group_by(slice) %>%
  mutate(FAP_bin = dplyr::ntile(FAP_logCP10K, nbin)) %>%
  group_by(slice, FAP_bin) %>%
  summarise(
    n = n(),
    x_med = median(FAP_logCP10K),
    y_med = median(TDS16_z),
    .groups = "drop"
  )
df_bin_all <- spatial_df %>%
  mutate(FAP_bin = dplyr::ntile(FAP_logCP10K, nbin)) %>%
  group_by(FAP_bin) %>%
  summarise(
    x_med = median(FAP_logCP10K),
    y_med = median(TDS16_z),
    .groups = "drop"
  )

write.csv(df_bin_slice, file.path(source_dir, "Fig2_spatial_binned_by_sample.csv"), row.names = FALSE)
write.csv(df_bin_all, file.path(source_dir, "Fig2_spatial_binned_all.csv"), row.names = FALSE)

fig2_f <- ggplot(
  df_bin_slice,
  aes(x = x_med, y = y_med, group = slice, color = slice)
) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 1.5, alpha = 0.95) +
  geom_line(
    data = df_bin_all,
    aes(x = x_med, y = y_med),
    inherit.aes = FALSE,
    linewidth = 1.8,
    color = "black"
  ) +
  geom_point(
    data = df_bin_all,
    aes(x = x_med, y = y_med),
    inherit.aes = FALSE,
    size = 2.0,
    color = "black"
  ) +
  scale_color_manual(values = sample_palette, drop = FALSE) +
  labs(
    title = sprintf("Dose-response per sample (nbin=%d)", nbin),
    x = "FAP (logCP10K) - bin median",
    y = "TDS (z-score) - bin median",
    color = "Sample"
  ) +
  theme_fig2(12) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.1, alpha = 1)))
save_plot(fig2_f, "FIG2F_dose_response_per_sample_rebuilt", panel_dir, width = 7.2, height = 5.2)

# NOTE: the experimental (wet-lab) panels were removed from this
# spatial-only version for the public bioinformatics repository.

fig2_partial <- (
  wrap_elements(full = fig2_c) | wrap_elements(full = fig2_b)
) / (
  fig2_e | fig2_f
) +
  plot_annotation(tag_levels = list(c("B", "C", "E", "F"))) &
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.01, 0.99)
  )
save_plot(fig2_partial, "FIGURE2_partial_spatial_validation_rebuilt", output_dir, width = 17, height = 22)

reference_candidates <- c(
  file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "AAAAA/Afigure/FIG2/FIG2 A-E.tif"),
  file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "A/FIG2/Fig2/FIG2B_FAP_logCP10K_300dpi.png"),
  file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "A/FIG2/Fig2/FIG2C_TDS16_z_300dpi.png"),
  file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "A/FIG2/Fig2/FIG2D_FAP_vs_TDS_300dpi.png")
)
for (ref in reference_candidates[file.exists(reference_candidates)]) {
  file.copy(ref, file.path(reference_dir, basename(ref)), overwrite = TRUE)
}

writeLines(
  c(
    "Figure 2 rebuild summary (spatial transcriptomics only)",
    paste0("Spatial object: ", spatial_rds_path),
    paste0("Spatial spots used: ", nrow(spatial_df)),
    paste0("Matched TDS16 genes: ", paste(tds_res$matched_union, collapse = ", ")),
    paste0("Fig2E spatial FAP vs TDS rho = ", sprintf("%.4f", unname(cor_res$estimate))),
    "Panels regenerated here: B, C, E, F (spatial transcriptomics only).",
    "Experimental (wet-lab) and image panels are not part of this spatial-only version."
  ),
  con = file.path(output_dir, "FIGURE2_rebuild_notes.txt")
)

message("Figure 2 rebuilt outputs written to: ", output_dir)
