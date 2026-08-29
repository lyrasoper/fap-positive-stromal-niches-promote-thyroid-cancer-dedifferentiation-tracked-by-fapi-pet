#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S9 (niche hallmark panels). The "supp11_15"
#          directory names below are earlier numbering.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.


suppressPackageStartupMessages({
  library(ggplot2)
})

base_dir <- file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "adobe_iiiu")
out_dir <- file.path(base_dir, "supp11_15_remade_vector_panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source_pkg <- file.path(base_dir, "supp11_15_source_based_for_ai")

save_pdf_png <- function(plot_obj, stem, width, height, dpi = 300) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  ggsave(pdf_path, plot = plot_obj, width = width, height = height, device = "pdf", bg = "white")
  ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = dpi, bg = "white")
  invisible(c(pdf = pdf_path, png = png_path))
}

# -------------------------
# Supp Fig. 11A remake
# -------------------------
supp11_src <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "spatial_visium/outputs/fap_caf_tumor_hallmarks_per_sample/per_sample_hallmark_stats_formal.csv")
supp11_df <- read.csv(supp11_src, stringsAsFactors = FALSE, check.names = FALSE)

sample_order <- c("P32", "P12", "P98", "P44", "P26", "P83", "P57", "P17")
hallmark_order <- c(
  "APICAL SURFACE",
  "TNFA SIGNALING VIA NFKB",
  "IL6 JAK STAT3 SIGNALING",
  "E2F TARGETS",
  "G2M CHECKPOINT",
  "INFLAMMATORY RESPONSE",
  "APICAL JUNCTION",
  "HYPOXIA"
)

supp11_sub <- supp11_df[supp11_df$sample %in% sample_order & supp11_df$label %in% hallmark_order, ]
supp11_long <- rbind(
  data.frame(
    sample = supp11_sub$sample,
    hallmark = supp11_sub$label,
    group = "Far",
    score = supp11_sub$mean_score_far,
    stringsAsFactors = FALSE
  ),
  data.frame(
    sample = supp11_sub$sample,
    hallmark = supp11_sub$label,
    group = "Near",
    score = supp11_sub$mean_score_near,
    stringsAsFactors = FALSE
  )
)
supp11_long$sample <- factor(supp11_long$sample, levels = sample_order)
supp11_long$group <- factor(supp11_long$group, levels = c("Far", "Near"))
supp11_long$hallmark <- factor(supp11_long$hallmark, levels = rev(hallmark_order))

p11 <- ggplot(supp11_long, aes(x = group, y = hallmark, fill = score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~ sample, nrow = 2) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Mean z-score") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

save_pdf_png(p11, "Supp11A_remade_vector", width = 9.4, height = 5.2)

# -------------------------
# Supp Fig. 13A remake
# -------------------------
supp13a_src <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/topic_service_reassessment_20260327/pathway_ranking_fap_dediff_v1/tables/pathway_ranking_integrated_scores.csv")
supp13a_df <- read.csv(supp13a_src, stringsAsFactors = FALSE, check.names = FALSE)
supp13a_df <- supp13a_df[order(supp13a_df$integrated_score, decreasing = TRUE), ]
supp13a_df$mechanism_label <- factor(supp13a_df$mechanism_label, levels = rev(supp13a_df$mechanism_label))

p13a <- ggplot(supp13a_df, aes(x = integrated_score, y = mechanism_label, fill = integrated_score)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.3f", integrated_score)), hjust = -0.12, size = 3.5) +
  scale_fill_gradient(low = "#dbe7f0", high = "#0b4f6c") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.2))) +
  labs(x = "Integrated pathway support score", y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

save_pdf_png(p13a, "Supp13A_remade_vector", width = 8.5, height = 4.8)

# -------------------------
# Supp Fig. 15A remake
# -------------------------
supp15a_src <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/topic_service_reassessment_20260327/stat3_pathway_activation_v1/tables/stat3_module_subtype_ranking_filtered_n100.csv")
supp15a_df <- read.csv(supp15a_src, stringsAsFactors = FALSE, check.names = FALSE)
supp15a_df <- head(supp15a_df[order(supp15a_df$mean_stat3_module_score, decreasing = TRUE), ], 15)
supp15a_df$subtype_full <- factor(supp15a_df$subtype_full, levels = rev(supp15a_df$subtype_full))

lineage_cols <- c(
  "Fibro" = "#f8766d",
  "Myeloid" = "#7CAE00",
  "T_NK" = "#00BFC4",
  "Thyrocyte" = "#C77CFF"
)

p15a <- ggplot(supp15a_df, aes(x = mean_stat3_module_score, y = subtype_full, fill = lineage)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", mean_stat3_module_score)), hjust = -0.12, size = 3.4) +
  scale_fill_manual(values = lineage_cols, drop = FALSE) +
  xlim(min(0, min(supp15a_df$mean_stat3_module_score) * 1.05), max(supp15a_df$mean_stat3_module_score) * 1.25) +
  labs(x = "Mean STAT3 module score", y = NULL, fill = "lineage") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

save_pdf_png(p15a, "Supp15A_remade_vector", width = 9, height = 6)

# -------------------------
# Supp Fig. 15B remake
# -------------------------
supp15b_src <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/topic_service_reassessment_20260327/bgn_tlr4_stat3_cd274_focus/tables/sample_level_targeted_axis_metrics.csv")
supp15b_df <- read.csv(supp15b_src, stringsAsFactors = FALSE, check.names = FALSE)

ligand_cols <- c(
  "fibro__FAP_CAF_core__BGN_mean",
  "fibro__FAP_CAF_core__VCAN_mean",
  "fibro__FAP_CAF_core__POSTN_mean",
  "fibro__FAP_CAF_core__CXCL12_mean"
)
ycol <- "thyro__dediff_candidate__CD274_mean"

ligand_labels <- c(
  "fibro__FAP_CAF_core__BGN_mean" = "FAP+ CAF BGN",
  "fibro__FAP_CAF_core__VCAN_mean" = "FAP+ CAF VCAN",
  "fibro__FAP_CAF_core__POSTN_mean" = "FAP+ CAF POSTN",
  "fibro__FAP_CAF_core__CXCL12_mean" = "FAP+ CAF CXCL12"
)

long_parts <- lapply(ligand_cols, function(col_name) {
  data.frame(
    ligand_metric = ligand_labels[[col_name]],
    x = supp15b_df[[col_name]],
    y = supp15b_df[[ycol]],
    Histology = supp15b_df$Histology,
    stringsAsFactors = FALSE
  )
})
supp15b_long <- do.call(rbind, long_parts)
supp15b_long <- supp15b_long[!is.na(supp15b_long$x) & !is.na(supp15b_long$y), ]
supp15b_long$ligand_metric <- factor(supp15b_long$ligand_metric, levels = unname(ligand_labels))
supp15b_long$Histology <- factor(supp15b_long$Histology, levels = c("ATC", "Normal", "PTC"))

hist_cols <- c("ATC" = "#F8766D", "Normal" = "#00BA38", "PTC" = "#619CFF")

p15b <- ggplot(supp15b_long, aes(x = x, y = y, color = Histology)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "#444444") +
  facet_wrap(~ ligand_metric, scales = "free_x") +
  scale_color_manual(values = hist_cols, drop = FALSE) +
  labs(
    x = "Sample-level mean expression in FAP+ CAF",
    y = "Sample-level mean CD274 in dediff epithelial"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "#ebebeb"),
    panel.grid.minor = element_blank()
  )

save_pdf_png(p15b, "Supp15B_remade_vector", width = 10.5, height = 4.5)

# -------------------------
# Copy remade outputs into source-based package
# -------------------------
copy_if_exists <- function(from_stem, subdir, panel_name) {
  for (ext in c("pdf", "png")) {
    from <- file.path(out_dir, paste0(from_stem, ".", ext))
    to <- file.path(source_pkg, subdir, paste0(panel_name, "_remade_vector.", ext))
    if (file.exists(from)) {
      file.copy(from, to, overwrite = TRUE)
    }
  }
}

copy_if_exists("Supp11A_remade_vector", "SuppFig11", "Panel_A")
copy_if_exists("Supp13A_remade_vector", "SuppFig13", "Panel_A")
copy_if_exists("Supp15A_remade_vector", "SuppFig15", "Panel_A")
copy_if_exists("Supp15B_remade_vector", "SuppFig15", "Panel_B")

manifest_lines <- c(
  "Remade vector-style panels",
  "Supp11A_remade_vector.pdf/png",
  "Supp13A_remade_vector.pdf/png",
  "Supp15A_remade_vector.pdf/png",
  "Supp15B_remade_vector.pdf/png",
  "",
  "These were regenerated from source CSV tables to provide cleaner Illustrator-editable assets.",
  "Supp11A uses per-sample hallmark stats.",
  "Supp13A uses integrated pathway ranking scores.",
  "Supp15A uses STAT3 subtype ranking scores.",
  "Supp15B uses sample-level targeted axis metrics."
)
writeLines(manifest_lines, file.path(out_dir, "README.txt"))
