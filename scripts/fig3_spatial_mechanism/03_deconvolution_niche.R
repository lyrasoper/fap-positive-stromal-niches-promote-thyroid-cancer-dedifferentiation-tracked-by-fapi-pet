#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Figure 3 / Supplementary Fig. S7 (RCTD deconvolution and niche states).
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.


suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(png)
  library(grid)
})

project_root <- Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>")
outdir <- file.path(project_root, "outputs", "fig7_draft")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  scrna = file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/new_seurat_final_2026-02-05.rds"),
  fib_ann = file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/results/candidate_state_labels_refined/fibro_candidate_annotations.csv"),
  thy_ann = file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/results/candidate_state_labels_refined/thyrocyte_candidate_annotations.csv"),
  my_ann = file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/results/targeted_reclustering_markers/Meyloid_cells/Meyloid_cells_refined_cell_annotations.csv"),
  B = file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "spatial_visium/outputs/rctd_formal_all9_fast_collection/figures/rctd_core_states_noP1_refined.png"),
  C = file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "spatial_visium/outputs/fap_caf_tumor_hallmarks_noP1/selected_tumor_hallmarks_barplot.png"),
  D = file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "spatial_visium/outputs/focused_tam_transfer_all9/figures/fap_tam_context_with_c1qc_noP1.png"),
  E = file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "spatial_visium/results/12_tri_niche_fap_spp1_dediff/tri_niche_summary_figure.png"),
  F_data = file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "spatial_visium/results/12_tri_niche_fap_spp1_dediff/malignant_tri_high_models.csv")
)

missing_files <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_files)) {
  stop("Missing input files: ", paste(missing_files, collapse = ", "))
}

theme_fig <- function(base_size = 10.5) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 0.5, hjust = 0),
      plot.subtitle = element_text(size = base_size - 1.5, hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

image_panel <- function(path, title = NULL, subtitle = NULL) {
  img <- png::readPNG(path)
  g <- rasterGrob(img, interpolate = TRUE)
  ggdraw() +
    draw_grob(g, 0, 0, 1, 1) +
    {
      if (!is.null(title)) {
        draw_label(title, x = 0.02, y = 0.985, hjust = 0, vjust = 1, fontface = "bold", size = 10.5)
      }
    } +
    {
      if (!is.null(subtitle)) {
        draw_label(subtitle, x = 0.02, y = 0.935, hjust = 0, vjust = 1, size = 8.3)
      }
    }
}

add_panel_label <- function(plot_obj, label) {
  ggdraw(plot_obj) +
    draw_label(label, x = 0.0, y = 1.02, hjust = 0, vjust = 1, fontface = "bold", size = 18)
}

build_panel_a <- function() {
  obj <- readRDS(paths$scrna)
  emb <- as.data.frame(Embeddings(obj, "umap"))
  emb$cell_id <- rownames(emb)
  colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")

  fib <- read.csv(paths$fib_ann, stringsAsFactors = FALSE)
  thy <- read.csv(paths$thy_ann, stringsAsFactors = FALSE)
  my <- read.csv(paths$my_ann, stringsAsFactors = FALSE)

  emb$ref_state <- "Other cells"
  emb$ref_state[match(fib$cell_id[fib$fibro_state == "FAP_CAF_core"], emb$cell_id, nomatch = 0)] <- "FAP_CAF_core"
  emb$ref_state[match(thy$cell_id[thy$epithelial_state == "dediff_candidate"], emb$cell_id, nomatch = 0)] <- "dediff_candidate"
  emb$ref_state[match(my$cell_id[my$myeloid_state == "SPP1_GPNMB_scar_like_TAM"], emb$cell_id, nomatch = 0)] <- "SPP1_GPNMB_scar_like_TAM"
  emb$ref_state[match(my$cell_id[my$myeloid_state == "C1QC_C3_CX3CR1_resident_macrophage"], emb$cell_id, nomatch = 0)] <- "C1QC_C3_CX3CR1_resident_macrophage"

  set.seed(1)
  other_ids <- which(emb$ref_state == "Other cells")
  keep_other <- if (length(other_ids) > 60000) sample(other_ids, 60000) else other_ids
  keep_idx <- c(keep_other, which(emb$ref_state != "Other cells"))
  plot_df <- emb[keep_idx, , drop = FALSE]

  plot_df$ref_state <- factor(
    plot_df$ref_state,
    levels = c(
      "Other cells",
      "FAP_CAF_core",
      "dediff_candidate",
      "SPP1_GPNMB_scar_like_TAM",
      "C1QC_C3_CX3CR1_resident_macrophage"
    )
  )

  pal <- c(
    "Other cells" = "#d1d5db",
    "FAP_CAF_core" = "#e07a5f",
    "dediff_candidate" = "#3d405b",
    "SPP1_GPNMB_scar_like_TAM" = "#81b29a",
    "C1QC_C3_CX3CR1_resident_macrophage" = "#6d597a"
  )

  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = ref_state)) +
    geom_point(
      data = plot_df[plot_df$ref_state == "Other cells", , drop = FALSE],
      size = 0.18,
      alpha = 0.10,
      stroke = 0
    ) +
    geom_point(
      data = plot_df[plot_df$ref_state != "Other cells", , drop = FALSE],
      size = 0.42,
      alpha = 0.88,
      stroke = 0
    ) +
    scale_color_manual(
      values = pal,
      labels = c(
        "Other cells",
        "FAP+ CAF",
        "Dediff tumor",
        "SPP1/GPNMB TAM",
        "C1QC/CX3CR1 macrophage"
      )
    ) +
    labs(
      title = "scRNA reference states used for spatial projection",
      subtitle = "Key sender, receiver, and myeloid niche states defined in the integrated single-cell atlas",
      x = "UMAP 1",
      y = "UMAP 2",
      color = NULL
    ) +
    theme_fig(10.2) +
    theme(
      legend.position = c(0.995, 0.50),
      legend.justification = c(1, 0.5),
      legend.background = element_rect(fill = scales::alpha("white", 0.75), color = NA),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      axis.title = element_text(size = 9),
      legend.text = element_text(size = 8.5),
      legend.key.height = unit(0.38, "cm"),
      legend.key.width = unit(0.28, "cm"),
      plot.title = element_text(size = 10.4, face = "bold"),
      plot.subtitle = element_text(size = 8.1)
    )

  ggsave(file.path(outdir, "Fig7A_scrna_reference_states_v7.png"), p, width = 7.1, height = 5.4, dpi = 320, bg = "white")
  ggsave(file.path(outdir, "Fig7A_scrna_reference_states_v7.pdf"), p, width = 7.1, height = 5.4, bg = "white")
  p
}

build_panel_f <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  keep <- c("TDS", "EMT score", "IL6/JAK/STAT3 score", "CD274", "Checkpoint score")
  df <- df[df$outcome_label %in% keep, , drop = FALSE]
  order_vec <- c("TDS", "EMT score", "IL6/JAK/STAT3 score", "CD274", "Checkpoint score")
  df$outcome_label <- factor(df$outcome_label, levels = rev(order_vec))
  df$direction <- ifelse(df$beta > 0, "Higher in tri-high", "Lower in tri-high")

  p <- ggplot(df, aes(x = beta, y = outcome_label)) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.5) +
    geom_errorbar(
      aes(xmin = ci_low, xmax = ci_high),
      orientation = "y",
      width = 0.18,
      linewidth = 0.65,
      color = "grey25"
    ) +
    geom_point(
      aes(fill = direction),
      shape = 21,
      size = 3.6,
      stroke = 0.45,
      color = "black"
    ) +
    scale_fill_manual(values = c("Higher in tri-high" = "#d76a6a", "Lower in tri-high" = "#4c78a8")) +
    labs(
      title = "Tri-high niches mark a more dedifferentiated malignant phenotype",
      subtitle = "Within malignant spots: TDS decreases while EMT, IL6/JAK/STAT3, CD274, and checkpoint score increase",
      x = "Beta in tri-high malignant spots",
      y = NULL
    ) +
    theme_fig(10.4) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(face = "bold"),
      axis.title.x = element_text(size = 9.0),
      plot.title = element_text(size = 10.2, face = "bold"),
      plot.subtitle = element_text(size = 8.0),
      axis.text.x = element_text(size = 8.3)
    )

  ggsave(file.path(outdir, "Fig7F_tri_high_malignant_phenotype_v7.png"), p, width = 6.2, height = 4.1, dpi = 320, bg = "white")
  ggsave(file.path(outdir, "Fig7F_tri_high_malignant_phenotype_v7.pdf"), p, width = 6.2, height = 4.1, bg = "white")
  p
}

panel_a <- build_panel_a()

panel_b <- image_panel(
  paths$B,
  title = "Formal RCTD spatial state maps",
  subtitle = "FAP+ CAF, dedifferentiated tumor, differentiated tumor, and SPP1-high myeloid states across design-matched sections"
)

panel_c <- image_panel(
  paths$C,
  title = "Tumor-related hallmarks near FAP+ CAF",
  subtitle = "Nearby tumor-proxy spots show enrichment of EMT, IL6/JAK/STAT3, inflammatory, hypoxia, and cell-cycle programs"
)

panel_d <- image_panel(
  paths$D,
  title = "FAP-associated TAM spatial states",
  subtitle = "Focused transfer highlights SPP1/GPNMB TAM and C1QC/CX3CR1 macrophage context around FAP-rich dedifferentiated regions"
)

panel_e <- image_panel(
  paths$E,
  title = "Tri-niche hotspot architecture",
  subtitle = "Projected fibroblast, myeloid, and dedifferentiated tumor states form a spatially coupled niche"
)

panel_f <- build_panel_f(paths$F_data)

row1 <- plot_grid(
  add_panel_label(panel_a, "A"),
  add_panel_label(panel_b, "B"),
  ncol = 2,
  rel_widths = c(0.86, 1.14)
)

row2 <- plot_grid(
  add_panel_label(panel_c, "C"),
  add_panel_label(panel_d, "D"),
  ncol = 2,
  rel_widths = c(0.95, 1.05)
)

row3 <- plot_grid(
  add_panel_label(panel_e, "E"),
  add_panel_label(panel_f, "F"),
  ncol = 2,
  rel_widths = c(1.08, 0.92)
)

title_block <- ggdraw() +
  draw_label(
    "Figure 7 draft v8 | Single-cell state projection reveals a fibro-myeloid niche linked to thyroid cancer dedifferentiation",
    x = 0, y = 0.84, hjust = 0, vjust = 1, fontface = "bold", size = 14.8
  ) +
  draw_label(
    "scRNA-defined states were projected onto spatial transcriptomics, revealing FAP+ CAF-rich niches associated with myeloid context, hallmark pathway activation, and tri-high malignant phenotypes.",
    x = 0, y = 0.34, hjust = 0, vjust = 1, size = 8.9
  )

fig <- plot_grid(
  title_block,
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(0.09, 0.31, 0.30, 0.30)
)

outfile_png <- file.path(outdir, "Figure7_deconvolution_mechanism_draft_v8.png")
outfile_pdf <- file.path(outdir, "Figure7_deconvolution_mechanism_draft_v8.pdf")

ggsave(outfile_png, fig, width = 17, height = 19, dpi = 320, bg = "white", limitsize = FALSE)
ggsave(outfile_pdf, fig, width = 17, height = 19, bg = "white", limitsize = FALSE)

note_lines <- c(
  "Figure 7 draft v8 adds a cleaner hallmark pathway panel to the main figure.",
  "Panel A: enlarged scRNA reference-state dots with legend moved into the plotting area.",
  "Panel B: formal RCTD spatial state maps (noP1 refined).",
  "Panel C: selected tumor-related hallmarks enriched in tumor-proxy spots near FAP+ CAF.",
  "Panel D: focused Seurat transfer TAM context with SPP1 and C1QC macrophage states.",
  "Panel E: tri-niche summary built from projected spatial states.",
  "Panel F: tri-high malignant phenotype built from malignant_tri_high_models.csv.",
  "Recommended supplementary materials for this version: near-vs-far hallmark heatmap, spatial LR network, scRNA pathway ranking, receiver-state mechanism heatmap, and full correlation tables."
)
writeLines(note_lines, file.path(outdir, "Figure7_deconvolution_mechanism_draft_v8_notes.txt"))

message("Saved: ", outfile_png)
message("Saved: ", outfile_pdf)
