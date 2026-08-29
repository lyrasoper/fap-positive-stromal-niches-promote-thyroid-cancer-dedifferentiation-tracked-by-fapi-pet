# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S21 (expanded mouse panels).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  library(patchwork)
  library(clustree)
  library(GSVA)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(gridExtra)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "scripts/_shared/paper_A_style.R"))

project_dir   <- Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>")
fig3_data_dir <- file.path(project_dir, "data/fig3")
full_obj_path <- file.path(fig3_data_dir, "epi_fib_imm_endo_seurat_object.rds")
fib_obj_path  <- file.path(fig3_data_dir, "fibroblast_seurat_object.rds")

bulk_dir      <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "mouse_bulk/mTC")

supp_dir      <- file.path(project_dir, "outputs/supplementary_figures")
src_dir       <- file.path(supp_dir, "source_data")
sub_src       <- file.path(src_dir, "supp_fig6")
dir.create(sub_src, recursive = TRUE, showWarnings = FALSE)

common_theme <- theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(text = element_text(face = "plain", color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "plain"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "plain"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.45),
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

tp_levels   <- c("1 month", "2 months", "4 months")
# paper_A 3-tone progression: 1->2->4 month mirrors the low->high
# dedifferentiation axis (cool blue -> amber -> coral red).
tp_palette  <- c("1 month"  = "#7BA6C9",
                 "2 months" = "#E8A341",
                 "4 months" = "#C03840")

lineage_palette <- c(Epithelial = "#1F5AA6", Fibroblast = "#B0302A",
                     Myeloid    = "#D76414", Endothelial = "#008C9E",
                     T_cell     = "#B88900")
caf_palette <- c("FAP+ inflaCAF"   = "#D76414",
                 "Ramp2+ endmtCAF" = "#008C9E",
                 "Rgs5+ myoCAF"    = "#4B1F7A",
                 "Sox10+ pnCAF"    = "#B88900")

# =============================================================
# Load mouse scRNA objects
# =============================================================
message("[1/6] Loading mouse scRNA objects ...")
obj <- readRDS(full_obj_path)
fib <- readRDS(fib_obj_path)

# derive time point from orig.ident (expected like mPTC_1month / 2month / 4month)
oi <- unique(as.character(obj$orig.ident))
message("  orig.ident values: ", paste(oi, collapse = ", "))
obj$Time <- dplyr::case_when(
  grepl("1[._ ]?month|_1m", obj$orig.ident, ignore.case = TRUE) ~ "1 month",
  grepl("2[._ ]?month|_2m", obj$orig.ident, ignore.case = TRUE) ~ "2 months",
  grepl("4[._ ]?month|_4m", obj$orig.ident, ignore.case = TRUE) ~ "4 months",
  TRUE ~ NA_character_)
if (all(is.na(obj$Time))) {
  # fallback: treat the 3 ordered levels in orig.ident as 1/2/4 months
  ord <- sort(unique(as.character(obj$orig.ident)))
  if (length(ord) == 3) {
    obj$Time <- dplyr::case_when(
      obj$orig.ident == ord[1] ~ "1 month",
      obj$orig.ident == ord[2] ~ "2 months",
      obj$orig.ident == ord[3] ~ "4 months")
  }
}
obj$Time <- factor(obj$Time, levels = tp_levels)
fib$Time <- factor(dplyr::case_when(
  grepl("1[._ ]?month|_1m", fib$orig.ident, ignore.case = TRUE) ~ "1 month",
  grepl("2[._ ]?month|_2m", fib$orig.ident, ignore.case = TRUE) ~ "2 months",
  grepl("4[._ ]?month|_4m", fib$orig.ident, ignore.case = TRUE) ~ "4 months"),
  levels = tp_levels)

# unify lineage column
obj$Lineage <- factor(obj$celltype,
                      levels = c("Epithelial","Fibroblast","Myeloid",
                                 "Endothelial","T_cell"))

# =============================================================
# Panel A: Mouse scRNA source overview table
# =============================================================
message("[2/6] Panels A, B, E, F, I ...")
scrna_overview <- as_tibble(obj@meta.data) %>%
  count(Time) %>%
  rename(`Time point` = Time, `All cells` = n) %>%
  left_join(
    as_tibble(fib@meta.data) %>% count(Time) %>%
      rename(`Time point` = Time, `Fibroblast/CAF cells` = n),
    by = "Time point")
write.csv(scrna_overview,
          file.path(sub_src, "supp_fig6_panelA_mouse_scRNA_overview.csv"),
          row.names = FALSE)

p6_a <- ggtexttable(scrna_overview, rows = NULL,
                    theme = ttheme("classic", base_size = 10,
                                   padding = unit(c(4, 3), "mm"))) %>%
  tab_add_title(text = "Mouse scRNA-seq source used for Fig. 3",
                face = "bold", size = 11, padding = unit(1, "lines"))
p6_a <- wrap_elements(full = p6_a)

# =============================================================
# Panel B: Independent mouse bulk transcriptomic cohorts
# =============================================================
bulk_overview <- tibble(
  cohort            = c("GSE55933", "GSE118022", "GSE30427"),
  n                 = c(10, 24, 28),
  `Displayed groups` = c("mPTC \u00d7 5, mATC \u00d7 5",
                        "mNT \u00d7 6, mPTC \u00d7 4, mATC \u00d7 5, remATC \u00d7 9",
                        "WT \u00d7 18, FTC \u00d7 5, ATC \u00d7 5"),
  `Primary use`     = c("Reference cohort",
                        "Fig. 3I / Supp. Fig. 6C",
                        "Supp. Fig. 6D"))
write.csv(bulk_overview,
          file.path(sub_src, "supp_fig6_panelB_bulk_mTC_overview.csv"),
          row.names = FALSE)

p6_b <- ggtexttable(bulk_overview, rows = NULL,
                    theme = ttheme("classic", base_size = 10,
                                   padding = unit(c(4, 3), "mm"))) %>%
  tab_add_title(text = "Independent mouse bulk transcriptomic cohorts",
                face = "bold", size = 11, padding = unit(1, "lines"))
p6_b <- wrap_elements(full = p6_b)

# =============================================================
# Panel E: Mouse scRNA QC by time point
# =============================================================
qc_df <- as_tibble(obj@meta.data) %>%
  select(orig.ident, Time, nFeature_RNA, nCount_RNA, percent.mt) %>%
  pivot_longer(cols = c(nFeature_RNA, nCount_RNA, percent.mt),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = dplyr::recode(metric,
          nFeature_RNA = "nFeature_RNA",
          nCount_RNA   = "nCount_RNA (log10)",
          percent.mt   = "percent.mt (%)"),
         value = ifelse(metric == "nCount_RNA (log10)",
                        log10(pmax(value, 1)), value))
write.csv(qc_df %>% group_by(Time, metric) %>%
            summarise(median = median(value), .groups = "drop"),
          file.path(sub_src, "supp_fig6_panelE_qc_median.csv"),
          row.names = FALSE)

p6_e <- ggplot(qc_df, aes(Time, value, fill = Time)) +
  geom_violin(scale = "width", trim = TRUE, width = 0.9,
              linewidth = 0.3, alpha = 0.88) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white",
               linewidth = 0.3, colour = "grey20") +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = tp_palette) +
  labs(x = NULL, y = NULL,
       title = "Per-cell QC metrics across time points") +
  common_theme +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10))

# =============================================================
# Panel F: Mouse major-lineage marker DotPlot
# =============================================================
m_lineage_features <- c("Epcam","Krt18","Krt19","Tg","Pax8","Fap","Col1a1",
                         "Pdgfra","Acta2","Ptprc","C1qa","Cd3d","Pecam1",
                         "Kdr","Vwf")
p6_f <- DotPlot(obj, features = m_lineage_features, group.by = "Lineage",
                dot.scale = 6, cols = c("grey92", "#C03840")) +
  RotatedAxis() + common_theme +
  labs(title = "Mouse lineage marker expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())
write.csv(p6_f$data,
          file.path(sub_src, "supp_fig6_panelF_lineage_dotplot.csv"),
          row.names = FALSE)

# =============================================================
# Panel G: Mouse CAF subtype marker DotPlot
# =============================================================
caf_subtypes <- c("FAP+ inflaCAF","Ramp2+ endmtCAF","Rgs5+ myoCAF","Sox10+ pnCAF")
# CAF identity comes from celltype_detail1 in the full object
caf_mask <- obj$celltype_detail1 %in% caf_subtypes
caf_obj  <- subset(obj, cells = colnames(obj)[caf_mask])
Idents(caf_obj) <- factor(caf_obj$celltype_detail1, levels = caf_subtypes)
m_caf_features <- c("Fap","Cxcl12","Il6","Postn","Col1a1","Ramp2","Plvap",
                    "Rgs5","Acta2","Tagln","Sox10","Mpz")
p6_g <- DotPlot(caf_obj, features = m_caf_features, dot.scale = 6,
                cols = c("grey92", "#7E5CAB")) +
  RotatedAxis() + common_theme +
  labs(title = "Mouse CAF subtype marker expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())
write.csv(p6_g$data,
          file.path(sub_src, "supp_fig6_panelG_CAF_dotplot.csv"),
          row.names = FALSE)

# =============================================================
# Panel H: CAF clustering resolution stability (clustree)
# =============================================================
res_cols_fib <- grep("^RNA_snn_res\\.", colnames(fib@meta.data), value = TRUE)
# pick a reasonable ladder
chosen_res <- intersect(c("RNA_snn_res.0.1","RNA_snn_res.0.2","RNA_snn_res.0.3",
                          "RNA_snn_res.0.4","RNA_snn_res.0.5","RNA_snn_res.0.6"),
                        res_cols_fib)
if (length(chosen_res) < 2 && length(res_cols_fib) >= 2) {
  chosen_res <- res_cols_fib[1:min(6, length(res_cols_fib))]
}
clust_meta <- fib@meta.data[, chosen_res, drop = FALSE]
write.csv(tibble::rownames_to_column(clust_meta, "cell"),
          file.path(sub_src, "supp_fig6_panelH_CAF_resolution_assignments.csv"),
          row.names = FALSE)
p6_h <- clustree(clust_meta, prefix = "RNA_snn_res.",
                 node_colour = "sc3_stability",
                 node_size_range = c(4, 10)) +
  labs(title = "CAF clustering resolution stability (mouse fibroblasts)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

# =============================================================
# Panel I: Per-sample cell composition across time points
# =============================================================
cell_comp <- as_tibble(obj@meta.data) %>%
  count(Time, Lineage) %>%
  group_by(Time) %>%
  mutate(fraction = n / sum(n)) %>% ungroup()
write.csv(cell_comp,
          file.path(sub_src, "supp_fig6_panelI_cell_composition.csv"),
          row.names = FALSE)
p6_i <- ggplot(cell_comp, aes(Time, fraction, fill = Lineage)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = lineage_palette, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Cell composition",
       title = "Cell composition by time point (mouse scRNA)") +
  common_theme +
  theme(legend.position = "right",
        legend.key.size = unit(0.4, "cm"))

# =============================================================
# Panel C: GSE118022 marker heatmap + Panel J: ssGSEA boxplots
# =============================================================
message("[3/6] GSE118022 bulk analysis ...")
read_bulk_matrix <- function(csv_path) {
  mat <- read.csv(csv_path, check.names = FALSE, row.names = 1)
  as.matrix(mat)
}
gse118 <- read_bulk_matrix(file.path(bulk_dir, "GSE118022_expr.csv"))
# Wait — GSE118022_expr.csv is the METADATA. The matrix is GSE118022_groupinf.csv
# (which is actually the expression matrix, confusingly named).
gse118_mat  <- read_bulk_matrix(file.path(bulk_dir, "GSE118022_groupinf.csv"))
gse118_meta <- read.csv(file.path(bulk_dir, "GSE118022_expr.csv"),
                        check.names = FALSE) %>%
  rename(SAMPLE_ID = sample, Group = group)

common_118 <- intersect(gse118_meta$SAMPLE_ID, colnames(gse118_mat))
gse118_meta <- gse118_meta %>% filter(SAMPLE_ID %in% common_118)
gse118_mat  <- gse118_mat[, gse118_meta$SAMPLE_ID, drop = FALSE]
message(sprintf("  GSE118022 aligned samples: %d", nrow(gse118_meta)))

m_diff_genes <- c("Tg","Tpo","Slc5a5","Nkx2-1","Pax8","Foxe1","Dio1","Dio2",
                  "Duox2","Iyd","Tshr","Slc26a4")
m_ecm_genes  <- c("Fap","Col1a1","Col1a2","Col3a1","Postn","Fn1","Sparc",
                  "Acta2","Lox","Dcn","Lum","Tagln","Ctgf","Serpine1")
m_diff_genes <- intersect(m_diff_genes, rownames(gse118_mat))
m_ecm_genes  <- intersect(m_ecm_genes,  rownames(gse118_mat))

# ---- Panel C: heatmap ----
group_order_118 <- c("mNT","mPTC","mATC","remATC")
gse118_meta$Group <- factor(gse118_meta$Group, levels = group_order_118)
ord_meta <- gse118_meta %>% arrange(Group, SAMPLE_ID)
heat_genes <- c(m_diff_genes, m_ecm_genes)
heat_mat   <- gse118_mat[heat_genes, ord_meta$SAMPLE_ID, drop = FALSE]
heat_mat_z <- t(scale(t(heat_mat)))

col_anno_c <- HeatmapAnnotation(
  Group = ord_meta$Group,
  col = list(Group = c("mNT" = "#7BA6C9", "mPTC" = "#E8A341",
                       "mATC" = "#C03840", "remATC" = "#7D1B1F")),
  show_annotation_name = FALSE,
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(Group = list(labels_gp = gpar(fontface = "bold"))))
row_anno_c <- rowAnnotation(
  Module = c(rep("Differentiation", length(m_diff_genes)),
             rep("FAP/ECM",         length(m_ecm_genes))),
  col = list(Module = c("Differentiation" = "#2C5E9D", "FAP/ECM" = "#C03840")),
  annotation_legend_param = list(Module = list(labels_gp = gpar(fontface = "bold"))),
  show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"))
col_fun_c <- colorRamp2(c(-2, 0, 2), c("#2C5E9D", "white", "#C03840"))

ht_c <- Heatmap(heat_mat_z, name = "z-score", col = col_fun_c,
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_column_names = FALSE,
                top_annotation = col_anno_c, left_annotation = row_anno_c,
                row_names_gp = gpar(fontface = "bold", fontsize = 8),
                row_split = factor(
                  c(rep("Differentiation", length(m_diff_genes)),
                    rep("FAP/ECM",         length(m_ecm_genes))),
                  levels = c("Differentiation","FAP/ECM")),
                row_title_gp = gpar(fontface = "bold", fontsize = 10),
                row_gap = unit(2, "mm"),
                column_title = "GSE118022: differentiation vs FAP/ECM genes",
                column_title_gp = gpar(fontface = "bold", fontsize = 11),
                column_split = ord_meta$Group,
                heatmap_legend_param = list(direction = "vertical",
                                            title_position = "topleft",
                                            legend_height = unit(3, "cm")))
p6_c <- wrap_elements(full = grid.grabExpr(draw(ht_c,
            heatmap_legend_side = "right", annotation_legend_side = "right",
            merge_legend = TRUE)))

# ---- Panel J: GSE118022 ssGSEA boxplots ----
message("  panel J: GSE118022 ssGSEA ...")
ssgsea_input_118 <- as.matrix(gse118_mat[
  unique(c(m_diff_genes, m_ecm_genes)), gse118_meta$SAMPLE_ID, drop = FALSE])
storage.mode(ssgsea_input_118) <- "numeric"
sspar_118 <- ssgseaParam(exprData = ssgsea_input_118,
                         geneSets = list(Differentiation_TDS = m_diff_genes,
                                         FAP_CAF_program     = m_ecm_genes),
                         normalize = FALSE, verbose = FALSE)
ssmat_118 <- gsva(sspar_118, verbose = FALSE)
ssdf_118  <- as.data.frame(t(ssmat_118)) %>%
  tibble::rownames_to_column("SAMPLE_ID") %>%
  left_join(gse118_meta, by = "SAMPLE_ID") %>%
  pivot_longer(cols = c("Differentiation_TDS","FAP_CAF_program"),
               names_to = "Module", values_to = "Score")
write.csv(ssdf_118,
          file.path(sub_src, "supp_fig6_panelJ_GSE118022_ssGSEA.csv"),
          row.names = FALSE)

cmp_118 <- list(c("mNT","mPTC"), c("mPTC","mATC"),
                c("mATC","remATC"), c("mNT","mATC"))
p6_j <- ggplot(ssdf_118, aes(Group, Score, fill = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.55, alpha = 0.9) +
  geom_jitter(width = 0.15, size = 1.1, shape = 21, stroke = 0.3, alpha = 0.75) +
  stat_compare_means(comparisons = cmp_118, method = "wilcox.test",
                     label = "p.signif", size = 2.8, fontface = "bold",
                     tip.length = 0.01, bracket.size = 0.4) +
  facet_wrap(~Module, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("mNT" = "#74C69D", "mPTC" = "#F4A261",
                               "mATC" = "#D62828", "remATC" = "#E76F51")) +
  labs(x = NULL, y = "ssGSEA score",
       title = "GSE118022 ssGSEA scores across Braf-driven mouse disease states") +
  common_theme +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 20, hjust = 1),
        panel.spacing = unit(0.8, "lines"))

# =============================================================
# Panel K: GSE30427 ssGSEA + Panel D: radar (reuse)
# =============================================================
message("[4/6] GSE30427 bulk analysis ...")
gse30 <- read.csv(file.path(bulk_dir, "GSE30427_expr.csv"),
                  check.names = FALSE, row.names = 1)
gse30_mat <- as.matrix(gse30)
# read groupinf via python-compatible skip (tolerate jagged rows)
raw <- readLines(file.path(bulk_dir, "GSE30427_groupinf.csv"))
hdr <- strsplit(raw[1], ",")[[1]]
gse30_meta_raw <- do.call(rbind, lapply(raw[-1], function(r) {
  f <- strsplit(r, ",", fixed = TRUE)[[1]]
  if (length(f) < 2) return(NULL)
  data.frame(title         = f[1],
             geo_accession = f[2], stringsAsFactors = FALSE)
})) %>% as_tibble()
gse30_meta <- gse30_meta_raw %>%
  mutate(Group = dplyr::case_when(
    grepl("wild type", title, ignore.case = TRUE)  ~ "WT",
    grepl("follicular", title, ignore.case = TRUE) ~ "FTC",
    grepl("anaplastic", title, ignore.case = TRUE) ~ "ATC",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(Group)) %>%
  rename(SAMPLE_ID = geo_accession) %>%
  mutate(Group = factor(Group, levels = c("WT","FTC","ATC")))
common_30 <- intersect(gse30_meta$SAMPLE_ID, colnames(gse30_mat))
gse30_meta <- gse30_meta %>% filter(SAMPLE_ID %in% common_30)
gse30_mat  <- gse30_mat[, gse30_meta$SAMPLE_ID, drop = FALSE]
message(sprintf("  GSE30427 aligned samples: %d", nrow(gse30_meta)))

m_diff_30 <- intersect(c("Tg","Tpo","Slc5a5","Nkx2-1","Pax8","Foxe1","Dio1",
                         "Dio2","Duox2","Iyd","Tshr","Slc26a4"), rownames(gse30_mat))
m_ecm_30  <- intersect(c("Fap","Col1a1","Col1a2","Col3a1","Postn","Fn1","Sparc",
                         "Acta2","Lox","Dcn","Lum","Tagln","Ctgf","Serpine1"),
                       rownames(gse30_mat))

ssgsea_input_30 <- as.matrix(gse30_mat[
  unique(c(m_diff_30, m_ecm_30)), gse30_meta$SAMPLE_ID, drop = FALSE])
storage.mode(ssgsea_input_30) <- "numeric"
sspar_30 <- ssgseaParam(exprData = ssgsea_input_30,
                        geneSets = list(Differentiation_TDS = m_diff_30,
                                        FAP_CAF_program     = m_ecm_30),
                        normalize = FALSE, verbose = FALSE)
ssmat_30 <- gsva(sspar_30, verbose = FALSE)
ssdf_30  <- as.data.frame(t(ssmat_30)) %>%
  tibble::rownames_to_column("SAMPLE_ID") %>%
  left_join(gse30_meta, by = "SAMPLE_ID") %>%
  pivot_longer(cols = c("Differentiation_TDS","FAP_CAF_program"),
               names_to = "Module", values_to = "Score")
write.csv(ssdf_30,
          file.path(sub_src, "supp_fig6_panelK_GSE30427_ssGSEA.csv"),
          row.names = FALSE)

cmp_30 <- list(c("WT","FTC"), c("FTC","ATC"), c("WT","ATC"))
p6_k <- ggplot(ssdf_30, aes(Group, Score, fill = Group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.55, alpha = 0.9) +
  geom_jitter(width = 0.15, size = 1.1, shape = 21, stroke = 0.3, alpha = 0.75) +
  stat_compare_means(comparisons = cmp_30, method = "wilcox.test",
                     label = "p.signif", size = 2.8, fontface = "bold",
                     tip.length = 0.01, bracket.size = 0.4) +
  facet_wrap(~Module, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("WT" = "#7BA6C9", "FTC" = "#E8A341",
                               "ATC" = "#C03840")) +
  labs(x = NULL, y = "ssGSEA score",
       title = "GSE30427 ssGSEA scores (WT vs FTC vs ATC)") +
  common_theme +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        panel.spacing = unit(0.8, "lines"))

# ---- Panel D: radar (reuse scaled csv) ----
message("[5/6] Panel D radar ...")
radar_path <- file.path(sub_src, "supp_fig6_panelD_gse30427_radar_scaled.csv")
if (file.exists(radar_path)) {
  radar_df <- read_csv(radar_path, show_col_types = FALSE)
  radar_long <- radar_df %>% pivot_longer(-gene, names_to = "Group", values_to = "Score") %>%
    mutate(Group = factor(Group, levels = c("WT","FTC","ATC")))
  radar_long$gene <- factor(radar_long$gene, levels = unique(radar_long$gene))
  p6_d <- ggplot(radar_long, aes(x = gene, y = Score, color = Group, group = Group)) +
    geom_polygon(aes(fill = Group), alpha = 0.18, linewidth = 0.6) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.6) +
    coord_polar() +
    scale_color_manual(values = c("WT" = "#7BA6C9", "FTC" = "#E8A341", "ATC" = "#C03840")) +
    scale_fill_manual(values = c("WT" = "#7BA6C9", "FTC" = "#E8A341", "ATC" = "#C03840")) +
    labs(title = "GSE30427 differentiation loss vs FAP/ECM activation",
         x = NULL, y = NULL, color = NULL, fill = NULL) +
    common_theme +
    theme(axis.text.x = element_text(size = 7, face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          legend.position = "top")
} else {
  stop(sprintf("Missing %s — produced by an upstream GSE30427 prep step; see README.", radar_path))
}

# =============================================================
# Assemble
# =============================================================
message("[6/6] Assembling Supp Fig 6 (7 panels, 4 rows) ...")
# Panels A, B (tables) removed — cohort counts now in Supplementary Table 7.
# Two more dropped 2026-05-20 to de-conflict with main Fig 3:
#   - cell composition by time point (p6_i) — duplicated Fig 3 C
#   - GSE118022 ssGSEA boxplots      (p6_j) — duplicated Fig 3 I
# Remaining 7 figure panels re-lettered a-g:
#   A = QC (was E)          B = lineage DotPlot (was F)
#   C = CAF DotPlot (was G)  D = clustree (was H)
#   E = GSE118022 heatmap (was C)
#   F = GSE30427 ssGSEA (was K)  G = GSE30427 radar (was D)
design6 <- "
AAAAAABBBBBB
CCCCCCDDDDDD
EEEEEEEEEEEE
FFFFFFFGGGGG
"
# Order matches alphabetical areas in design6:
#   A B C D E F G
supp6 <- p6_e + p6_f + p6_g + p6_h + p6_c + p6_k + p6_d +
  plot_layout(design = design6,
              heights = c(1.0, 1.0, 1.2, 1.15)) +
  plot_annotation(
    title = "Supplementary Fig. 6 | Source overview and extended transcriptomic validation for the time-resolved mouse model",
    caption = "Mouse scRNA-seq and bulk mTC cohort tables are provided in Supplementary Table 7.",
    tag_levels = "a"
  ) &
  theme(plot.title = element_text(size = 13.5, face = "bold", family = "Helvetica"),
        plot.tag   = element_text(size = 14, face = "bold", family = "Helvetica"),
        plot.caption = element_text(size = 9, face = "italic",
                                    hjust = 0, colour = "grey30"))

save_combo(supp6,
  "Supplementary_Fig_6_Source_overview_and_extended_transcriptomic_validation_for_the_time_resolved_mouse_model",
  width = 16, height = 20)

message("DONE Supp Fig 6")
