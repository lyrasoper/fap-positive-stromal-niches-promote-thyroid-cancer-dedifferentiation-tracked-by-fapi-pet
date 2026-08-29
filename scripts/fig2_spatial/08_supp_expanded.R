# Cancer Research submission - figure code release
# Builds: Supplementary figures for Figure 2 (expanded panels).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "scripts/_shared/paper_A_style.R"))

project_dir <- Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>")
src_dir     <- file.path(project_dir, "outputs/supplementary_figures/source_data")
supp_dir    <- file.path(project_dir, "outputs/supplementary_figures")
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

common_theme <- theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(text = element_text(face = "plain", color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "plain"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "plain"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        panel.border = element_rect(colour = "grey30", fill = NA,
                                    linewidth = 0.45),
        axis.line = element_blank())

save_combo <- function(plot, stem, width, height) {
  # Unicode-safe output: ragg PNG + macOS quartz PDF. The base R pdf() device
  # with the built-in Helvetica AFM is Latin-1-only and mangles Greek (rho/
  # sigma); cairo is unavailable here (no XQuartz).
  png_path <- file.path(supp_dir, paste0(stem, ".png"))
  ragg::agg_png(png_path, width = width, height = height, units = "in", res = 500)
  print(plot); dev.off()

  pdf_path <- file.path(supp_dir, paste0(stem, ".pdf"))
  if (capabilities("cairo")) grDevices::cairo_pdf(pdf_path, width = width, height = height) else if (capabilities("aqua")) grDevices::quartz(file = pdf_path, type = "pdf", width = width, height = height) else grDevices::pdf(pdf_path, width = width, height = height)
  print(plot); dev.off()
}

# ---------------- Load source data ----------------
# Source CSVs are produced by an upstream Supp Fig 4 source-data export step; see README.
message("[1/4] Loading Supp Fig 4 source CSVs ...")
read_csv_guarded <- function(p) {
  if (!file.exists(p))
    stop(sprintf("Missing %s — produced by an upstream step; see README.", p))
  read_csv(p, show_col_types = FALSE)
}
umap_df     <- read_csv_guarded(file.path(src_dir, "Supplementary_Fig4AB_umap_spot_metadata.csv"))
marker_df   <- read_csv_guarded(file.path(src_dir, "Supplementary_Fig4C_marker_dotplot_source.csv"))
cluster_frac<- read_csv_guarded(file.path(src_dir, "Supplementary_Fig4D_section_cluster_fraction.csv"))
section_sum <- read_csv_guarded(file.path(src_dir, "Supplementary_Fig4E_section_summary_metrics.csv"))
module_sc   <- read_csv_guarded(file.path(src_dir, "Supplementary_Fig4_cluster_module_scores.csv"))

sample_levels <- c("P1","P32","P12","P98","P44","P26","P83","P57","P17")
main_fig2_sections <- c("P12","P17","P26","P32","P44","P57","P83","P98")
umap_df$sample_id <- factor(umap_df$sample_id, levels = sample_levels)
sample_cols <- setNames(
  c("grey60","#6699CC","#E69F00","#009E73","#CC79A7","#56B4E9","#C03840","#F0E442","#C03840"),
  sample_levels)

# ---------------- Panel A: UMAP by section ----------------
message("[2/4] Rebuilding A-E ...")
p4_a <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = sample_id)) +
  geom_point(size = 0.35, alpha = 0.75) +
  scale_color_manual(values = sample_cols, drop = FALSE, name = "Section") +
  labs(title = "Section identity") +
  common_theme +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.key.size = unit(0.35, "cm"))

# ---------------- Panel B: UMAP by spot cluster ----------------
n_clust <- length(unique(umap_df$seurat_clusters))
p4_b <- ggplot(umap_df, aes(UMAP_1, UMAP_2,
                            color = factor(seurat_clusters))) +
  geom_point(size = 0.35, alpha = 0.75) +
  scale_color_manual(values = colorRampPalette(
    c("#D62728","#E69F00","#009E73","#1F77B4","#6A4C93","#C03840",
      "#56B4E9","#C03840","#8C564B","#2CA02C","#AA4488","#FFC857"))(n_clust),
    name = "Cluster") +
  labs(title = "Unsupervised spot clusters") +
  common_theme +
  guides(color = "none")

# ---------------- Panel C: canonical marker DotPlot ----------------
# Columns in Fig4C CSV: avg.exp, pct.exp, features.plot, id, avg.exp.scaled, feature.groups
marker_df$id <- factor(marker_df$id,
                       levels = as.character(sort(unique(
                         suppressWarnings(as.integer(marker_df$id))))))
p4_c <- ggplot(marker_df, aes(id, features.plot,
                              size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient2(low = "#2C5E9D", mid = "white", high = "#C03840",
                        midpoint = 0, name = "Scaled\naverage") +
  scale_size(range = c(0.5, 6), name = "% expressed") +
  labs(title = "Canonical marker landscape across spot clusters",
       x = "Spot cluster", y = NULL) +
  common_theme +
  theme(axis.text.x = element_text(angle = 0),
        axis.text.y = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"))

# ---------------- Panel D: section-wise cluster composition ----------------
cluster_frac$sample_id <- factor(cluster_frac$sample_id, levels = sample_levels)
p4_d <- ggplot(cluster_frac,
               aes(factor(seurat_clusters), sample_id,
                   fill = cluster_fraction)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#1F5AA6", name = "Spot\nfraction") +
  labs(title = "Section-wise cluster composition",
       x = "Spot cluster", y = NULL) +
  common_theme +
  theme(axis.text.y = element_text(size = 8),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.45))

# ---------------- Panel E: section-level summary metrics heatmap ----------------
section_sum$sample_id <- factor(section_sum$sample_id, levels = sample_levels)
metric_levels <- c("FAP_positive_fraction","mean_FAP","median_TDS_score",
                   "median_nFeature_Spatial","median_nCount_Spatial","spot_count_total")
section_sum$metric <- factor(section_sum$metric, levels = metric_levels)
p4_e <- ggplot(section_sum, aes(sample_id, metric, fill = z_score)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2C5E9D", mid = "white", high = "#C03840",
                       midpoint = 0, name = "Z-score") +
  labs(title = "Section-level summary metrics",
       x = "Section (ordered by median TDS)", y = NULL) +
  common_theme +
  theme(axis.text.y = element_text(size = 8))

# ---------------- Panel F (NEW): section-level FAP+ fraction vs median TDS scatter ----------------
message("[3/4] Building F, G, H, I ...")
sec_wide <- section_sum %>%
  select(sample_id, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value)
write.csv(sec_wide, file.path(src_dir, "Supplementary_Fig4F_section_FAP_vs_TDS.csv"),
          row.names = FALSE)

f_cor <- with(sec_wide, suppressWarnings(cor.test(FAP_positive_fraction,
                                                   median_TDS_score,
                                                   method = "spearman",
                                                   exact = FALSE)))
p4_f <- ggplot(sec_wide, aes(FAP_positive_fraction, median_TDS_score, color = sample_id)) +
  geom_point(size = 3.8, show.legend = FALSE) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey30",
              linewidth = 0.7) +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 3,
                           show.legend = FALSE, min.segment.length = 0,
                           fontface = "bold") +
  scale_color_manual(values = sample_cols, drop = FALSE) +
  labs(x = "Section FAP-positive spot fraction",
       y = "Section median TDS (z-score)",
       title = sprintf(
         "Section-level FAP+ fraction vs TDS (Spearman \u03c1 = %.2f, P = %s)",
         f_cor$estimate,
         ifelse(f_cor$p.value < 0.001, "< 0.001",
                formatC(f_cor$p.value, format = "e", digits = 2)))) +
  common_theme

# ---------------- Panel G (NEW): per-section dominant module composition ----------------
# Fig4D already has dominant_class per cluster.  Aggregate per sample.
cluster_frac_mod <- cluster_frac %>%
  group_by(sample_id, dominant_class) %>%
  summarise(spot_fraction = sum(cluster_fraction, na.rm = TRUE), .groups = "drop")
cluster_frac_mod$dominant_class <- factor(
  cluster_frac_mod$dominant_class,
  levels = c("Thyroid epithelial","ECM/FAP stroma","Myofibro/pericyte",
             "Myeloid","Lymphoid","Endothelial","Unassigned"))
mod_pal <- c(
  "Thyroid epithelial" = "#1F5AA6",
  "ECM/FAP stroma"     = "#B0302A",
  "Myofibro/pericyte"  = "#8B5FBF",
  "Myeloid"            = "#D76414",
  "Lymphoid"           = "#B88900",
  "Endothelial"        = "#008C9E",
  "Unassigned"         = "grey70")
write.csv(cluster_frac_mod,
          file.path(src_dir, "Supplementary_Fig4G_section_module_composition.csv"),
          row.names = FALSE)
p4_g <- ggplot(cluster_frac_mod,
               aes(sample_id, spot_fraction, fill = dominant_class)) +
  geom_col(width = 0.78, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = mod_pal, name = NULL,
                    na.translate = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = "Section (ordered by median TDS)",
       y = "Spot composition",
       title = "Per-section dominant module composition") +
  common_theme +
  theme(legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        axis.text.x = element_text(size = 9))

# ---------------- Panel H (NEW): FAP+ spot fraction distribution by section grade ----------------
# use spot-level UMAP_df to compute per-section FAP positivity
spot_fap <- umap_df %>%
  mutate(FAP_positive = FAP > 0) %>%
  group_by(sample_id) %>%
  summarise(fap_pos_frac = mean(FAP_positive, na.rm = TRUE),
            n_spots = n(),
            median_TDS = median(TDS_score, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(TDS_tertile = cut(median_TDS,
                           breaks = quantile(median_TDS, c(0, 1/3, 2/3, 1),
                                             na.rm = TRUE),
                           labels = c("Low TDS", "Mid TDS", "High TDS"),
                           include.lowest = TRUE))
write.csv(spot_fap, file.path(src_dir, "Supplementary_Fig4H_per_section_FAP_TDS_tertile.csv"),
          row.names = FALSE)
p4_h <- ggplot(spot_fap, aes(TDS_tertile, fap_pos_frac, fill = TDS_tertile)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.55, alpha = 0.9) +
  geom_jitter(aes(color = sample_id), width = 0.15, size = 2.5,
              show.legend = FALSE) +
  scale_fill_manual(values = c("Low TDS" = "#7BAAF7",
                               "Mid TDS" = "#2FBF71",
                               "High TDS" = "#E78B8B")) +
  scale_color_manual(values = sample_cols, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Section FAP+ spot fraction",
       title = "FAP positivity is enriched in low-TDS sections") +
  common_theme +
  theme(legend.position = "none")

# ---------------- Panel I (NEW): spot-level FAP vs TDS density ----------------
p4_i <- ggplot(umap_df, aes(FAP, TDS_score)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon",
                  contour_var = "ndensity") +
  geom_smooth(method = "lm", se = FALSE, colour = "#B0302A", linewidth = 0.7) +
  scale_fill_gradient(low = "white", high = "#1F5AA6", guide = "none") +
  labs(x = "FAP (log-norm)", y = "TDS score (z)",
       title = "Spot-level joint density (all 9 sections pooled)") +
  common_theme

# ---------------- Assemble ----------------
message("[4/4] Assemble Supp Fig 4 (9 panels) ...")
# 7-row layout:
#   Row 1: A | B (2 UMAPs)
#   Row 2-3: C (marker dotplot, tall for readability)
#   Row 4: D (section composition tile) | E (metrics heatmap)
#   Row 5: F | G (section scatter | module stack)
#   Row 6: H | I (TDS tertile | joint density)
design4 <- "
AAAABBBB
CCCCCCCC
CCCCCCCC
DDDDEEEE
FFFFGGGG
HHHHIIII
"
supp4 <- p4_a + p4_b + p4_c + p4_d + p4_e + p4_f + p4_g + p4_h + p4_i +
  plot_layout(design = design4,
              heights = c(1.0, 1.0, 1.0, 1.0, 1.1, 1.1)) +
  plot_annotation(
    title = "Supplementary Fig. 4 | Spatial transcriptomics cohort overview and marker-guided spot-cluster architecture",
    tag_levels = "a"
  ) &
  theme(plot.title = element_text(size = 13.5, face = "bold", family = "Helvetica"),
        plot.tag   = element_text(size = 14, face = "bold", family = "Helvetica"))

save_combo(supp4,
  "Supplementary_Fig_4_Spatial_transcriptomics_cohort_overview_and_marker_guided_spot_cluster_architecture",
  width = 14, height = 21)
message("DONE Supp Fig 4")
