# Cancer Research submission - figure code release
# Builds: Figure 3 (single-cell reference panels).
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

# Single-cell panels for the Fig 3 spatial-mechanism figure (human scRNA:
# epithelial / fibroblast / immune / endothelial objects + fibroblast subset).
# Output files retain legacy "Fig3B-H" panel labels from an earlier manuscript
# version. Reads <EXTERNAL_DATA>/... Seurat objects; see ../../README.md.

options(stringsAsFactors = FALSE)

required_pkgs <- c("Seurat", "dplyr", "ggplot2", "patchwork", "ggsci", "tibble", "scales", "Nebulosa")
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
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggsci)
  library(tibble)
  library(scales)
  library(Nebulosa)
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

theme_afigure <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )
}

save_plot <- function(plot_obj, stem, out_dir, width = 6, height = 5) {
  ggsave(
    filename = file.path(out_dir, paste0(stem, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 320,
    bg = "white",
    device = grDevices::png
  )
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  tryCatch(
    {
      grDevices::pdf(pdf_path, width = width, height = height, bg = "white", useDingbats = FALSE)
      print(plot_obj)
      grDevices::dev.off()
    },
    error = function(e) {
      if (names(grDevices::dev.cur()) != "null device") {
        grDevices::dev.off()
      }
      message("PDF export skipped for ", stem, ": ", conditionMessage(e))
    }
  )
}

project_dir <- Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>")
external_data <- Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>")
data_dir <- file.path(external_data, "fig3")
out_dir <- file.path(project_dir, "outputs", "fig3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

full_obj_path <- file.path(data_dir, "epi_fib_imm_endo_seurat_object.rds")
fibro_obj_path <- file.path(data_dir, "fibroblast_seurat_object.rds")
gsea_top10_path <- file.path(data_dir, "Fig3H_Top10_posNES_table.csv")

required_files <- c(full_obj_path, fibro_obj_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

scRNA_ha1 <- readRDS(full_obj_path)
fibroblast <- readRDS(fibro_obj_path)

if (!"celltype" %in% colnames(scRNA_ha1[[]])) {
  stop("The full Seurat object does not contain a 'celltype' column.")
}
if (!"orig.ident" %in% colnames(scRNA_ha1[[]])) {
  stop("The full Seurat object does not contain an 'orig.ident' column.")
}
if (!"RNA_snn_res.2" %in% colnames(fibroblast[[]])) {
  stop("The fibroblast Seurat object does not contain 'RNA_snn_res.2'.")
}

if (!"pca" %in% names(scRNA_ha1@reductions)) {
  set.seed(1)
  scRNA_ha1 <- RunPCA(scRNA_ha1, verbose = FALSE)
}
if (!"umap" %in% names(fibroblast@reductions)) {
  set.seed(1)
  fibroblast <- RunUMAP(fibroblast, reduction = "pca", dims = 1:30)
}

celltype_order <- c("Fibroblast", "Epithelial", "Myeloid", "Endothelial", "T_cell", "NK_cells", "B_cell", "Neutrophils", "other")
celltype_palette <- c(
  "Fibroblast" = "#bf362c",
  "Epithelial" = "#2e64ae",
  "Myeloid" = "#5a2d91",
  "Endothelial" = "#1698a6",
  "T_cell" = "#df6f10",
  "NK_cells" = "#ea7d11",
  "B_cell" = "#c79400",
  "Neutrophils" = "#36a85a",
  "other" = "#8d8d8d"
)

scRNA_ha1$celltype <- factor(scRNA_ha1$celltype, levels = celltype_order)
sample_order <- c("mPTC_1month", "mPTC_2month", "mPTC_4month")
scRNA_ha1$orig.ident <- factor(scRNA_ha1$orig.ident, levels = sample_order)

p_fig3b <- DimPlot(
  scRNA_ha1,
  reduction = "pca",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  scale_colour_manual(values = celltype_palette, drop = FALSE) +
  labs(title = "celltype", x = "PC_1", y = "PC_2") +
  theme_afigure(base_size = 18)
save_plot(p_fig3b, "Fig3B_Celltype_PCA_recovered", out_dir, width = 8.6, height = 6.8)

cell_prop_df <- as.data.frame(prop.table(table(scRNA_ha1$celltype, scRNA_ha1$orig.ident), margin = 2))
colnames(cell_prop_df) <- c("celltype", "Sample", "Proportion")
cell_prop_df$celltype <- factor(cell_prop_df$celltype, levels = celltype_order)
cell_prop_df$Sample <- factor(cell_prop_df$Sample, levels = sample_order)

fibro_line_df <- cell_prop_df %>% filter(celltype == "Fibroblast")

p_fig3c <- ggplot(cell_prop_df, aes(x = Sample, y = Proportion, fill = celltype)) +
  geom_col(width = 0.82, colour = "white") +
  geom_line(
    data = fibro_line_df,
    aes(x = Sample, y = Proportion, group = 1),
    inherit.aes = FALSE,
    linewidth = 1.2,
    colour = "black"
  ) +
  geom_point(
    data = fibro_line_df,
    aes(x = Sample, y = Proportion),
    inherit.aes = FALSE,
    size = 5.5,
    colour = "black"
  ) +
  geom_text(
    data = fibro_line_df,
    aes(x = Sample, y = Proportion, label = percent(Proportion, accuracy = 0.1)),
    inherit.aes = FALSE,
    vjust = -1.2,
    fontface = "bold",
    size = 6
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = celltype_palette, drop = FALSE) +
  labs(
    title = "Cell Type Proportion in Each Sample",
    x = "Sample",
    y = "Proportion"
  ) +
  theme_afigure(base_size = 17)
save_plot(p_fig3c, "Fig3C_Celltype_Proportion_FibroHighlighted_recovered", out_dir, width = 8.2, height = 6.8)

tds_genes_mouse <- c("Tg", "Tpo", "Pax8", "Dio1", "Dio2", "Duox1", "Duox2", "Foxe1", "Glis3", "Nkx2-1", "Slc26a4", "Slc5a5", "Slc5a8", "Thra", "Thrb", "Tshr")
tds_genes_mouse <- intersect(tds_genes_mouse, rownames(scRNA_ha1))

if (!"TDS_Score1" %in% colnames(scRNA_ha1[[]])) {
  set.seed(1)
  scRNA_ha1 <- AddModuleScore(scRNA_ha1, features = list(tds_genes_mouse), name = "TDS_Score")
}

p_fig3d <- DotPlot(
  scRNA_ha1,
  features = c("Fap", "TDS_Score1"),
  group.by = "orig.ident"
) +
  scale_colour_gradientn(colours = c("#431878", "#336699", "#66CC66", "#FFD43B")) +
  labs(x = "Features", y = "Identity") +
  theme_afigure(base_size = 16)
save_plot(p_fig3d, "Fig3D_DotPlot_FAP_TDS_Score_recovered", out_dir, width = 5.8, height = 5.8)

fibroblast$Fib_clu <- NA_character_
fibroblast$Fib_clu[fibroblast$RNA_snn_res.2 %in% c(11, 13)] <- "Ramp2+ endmtCAF"
fibroblast$Fib_clu[fibroblast$RNA_snn_res.2 %in% c(6, 17, 19, 22)] <- "Rgs5+ myoCAF"
fibroblast$Fib_clu[fibroblast$RNA_snn_res.2 == 20] <- "Sox10+ pnCAF"
fibroblast$Fib_clu[fibroblast$RNA_snn_res.2 %in% c(3, 5, 7, 8, 10, 14, 16)] <- "FAP+ inflaCAF"
fibroblast$Fib_clu[is.na(fibroblast$Fib_clu)] <- "FAP+ inflaCAF"

fibro_order <- c("FAP+ inflaCAF", "Ramp2+ endmtCAF", "Rgs5+ myoCAF", "Sox10+ pnCAF")
fibro_palette <- c(
  "FAP+ inflaCAF" = "#e16c0d",
  "Ramp2+ endmtCAF" = "#1a99aa",
  "Rgs5+ myoCAF" = "#5c2d91",
  "Sox10+ pnCAF" = "#bf9000"
)

fibroblast$Fib_clu <- factor(fibroblast$Fib_clu, levels = fibro_order)
fibroblast$orig.ident <- factor(fibroblast$orig.ident, levels = sample_order)

p_fig3e <- DimPlot(
  fibroblast,
  reduction = "umap",
  group.by = "Fib_clu",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  scale_colour_manual(values = fibro_palette, drop = FALSE) +
  labs(x = "umap_1", y = "umap_2") +
  theme_afigure(base_size = 18)
save_plot(p_fig3e, "Fig3E_mouseCAF_Subtypes_UMAP_recovered", out_dir, width = 8.0, height = 6.2)

p_fig3f <- FeaturePlot(
  fibroblast,
  features = "Fap",
  reduction = "umap",
  raster = FALSE
) +
  scale_colour_gradientn(colours = c("#f7fbff", "#6baed6", "#08306b")) +
  labs(x = "umap_1", y = "umap_2") +
  theme_afigure(base_size = 17)
save_plot(p_fig3f, "Fig3F_Fap_Density_mouseCAF_recovered", out_dir, width = 7.2, height = 6.0)

fibro_prop_df <- as.data.frame(prop.table(table(fibroblast$Fib_clu, fibroblast$orig.ident), margin = 2))
colnames(fibro_prop_df) <- c("Fib_clu", "Sample", "Proportion")
fibro_prop_df$Fib_clu <- factor(fibro_prop_df$Fib_clu, levels = fibro_order)
fibro_prop_df$Sample <- factor(fibro_prop_df$Sample, levels = sample_order)

fap_line_df <- fibro_prop_df %>% filter(Fib_clu == "FAP+ inflaCAF")

p_fig3g <- ggplot(fibro_prop_df, aes(x = Sample, y = Proportion, fill = Fib_clu)) +
  geom_col(width = 0.82, colour = "white") +
  geom_line(
    data = fap_line_df,
    aes(x = Sample, y = Proportion, group = 1),
    inherit.aes = FALSE,
    linewidth = 1.2,
    colour = "black"
  ) +
  geom_point(
    data = fap_line_df,
    aes(x = Sample, y = Proportion),
    inherit.aes = FALSE,
    size = 5.5,
    colour = "black"
  ) +
  geom_text(
    data = fap_line_df,
    aes(x = Sample, y = Proportion, label = percent(Proportion, accuracy = 0.1)),
    inherit.aes = FALSE,
    vjust = -1.2,
    fontface = "bold",
    size = 6
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = fibro_palette, drop = FALSE) +
  labs(x = "Sample", y = "Proportion") +
  theme_afigure(base_size = 17)
save_plot(p_fig3g, "Fig3G_mouseCAF_Subtype_Proportion_FAPlabel_recovered", out_dir, width = 8.0, height = 6.2)

if (file.exists(gsea_top10_path)) {
  gsea_top10 <- read.csv(gsea_top10_path, check.names = FALSE) %>%
    mutate(Pathway = factor(Pathway, levels = rev(Pathway)))

  p_fig3h <- ggplot(gsea_top10, aes(x = Pathway, y = NES)) +
    geom_col(fill = "#f1871a") +
    coord_flip() +
    labs(title = "GSEA | FAP+ CAF (Mouse)", x = NULL, y = "NES") +
    theme_afigure(base_size = 17)
  save_plot(p_fig3h, "Fig3H_GSEA_FAPposCAF_recovered", out_dir, width = 10.5, height = 5.8)
}

message("Recovered outputs written to: ", out_dir)
