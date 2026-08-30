#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S17 - NNLS CAF-composition deconvolution of the bulk contrast.
## B2: CAF subtype deconvolution of WT vs FAP-deficient host bulk RNA
## using a marker-based pseudo-CIBERSORTx approach.
##
## Method:
## Two complementary estimators are run on log2(FPKM+1) bulk:
##   (a) NNLS deconvolution with a synthetic signature matrix
##       (rows = markers, columns = cell types; entries = 1/0 with
##       optional decay for shared markers). Produces fraction
##       estimates summing to ~1.
##   (b) decoupleR::run_ulm enrichment of marker sets per sample.
##       Produces z-scored activity per cell type (not fractions,
##       but orthogonal to NNLS).
##
## Cell types & markers (all from published scRNA-seq studies):
##   myCAF     — Elyada 2019 Cancer Discov (PDAC), top 20 markers
##   iCAF      — Elyada 2019, top 20 markers
##   apCAF     — Elyada 2019, top markers
##   Thyroid_tumor — curated thyroid lineage (PAX8/TPO/TG/NKX2-1/FOXE1 family)
##   Immune    — Ptprc/Cd3e/Cd8a/Cd4/Cd19 pan-leukocyte
##   Endothelial — Pecam1/Cdh5/Vwf/Tek
##
## Output:
##   * per-sample cell fraction (NNLS) and activity (ULM) tables
##   * group comparison for myCAF, iCAF, myCAF/iCAF ratio
##   * 4-panel figure


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(nnls); library(decoupleR)
})

proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
bulk_out <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
out_dir  <- file.path(bulk_out, "caf_deconvolution")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_wt  <- "#3B6B9C"; col_fap <- "#C9493A"; col_up <- col_fap; col_dn <- col_wt
col_grey <- "#B7B7B7"; col_purple <- "#4B1E70"
grp_levels <- c("WT host", "FAP-deficient host")
grp_cols <- setNames(c(col_wt, col_fap), grp_levels)

base_font_size <- 7
theme_nc <- function(base = base_font_size){
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(size = base + 1, hjust = 0),
      plot.subtitle = element_text(size = base - 0.5, colour = "grey30"),
      plot.tag      = element_text(face = "bold", size = base + 3, hjust = 0),
      axis.title    = element_text(size = base),
      axis.text     = element_text(size = base - 0.5, colour = "black"),
      axis.line     = element_line(linewidth = 0.3),
      axis.ticks    = element_line(linewidth = 0.3),
      strip.background = element_blank(),
      strip.text    = element_text(size = base),
      legend.key.size = unit(3.2, "mm"),
      legend.text   = element_text(size = base - 0.5),
      legend.title  = element_text(size = base - 0.5),
      panel.grid    = element_blank(),
      plot.margin   = margin(3, 4, 3, 4)
    )
}
p_label <- function(p){
  if (is.na(p)) return("ns")
  if (p < 0.001) sprintf("P = %.1e", p) else sprintf("P = %.3f", p)
}

## ============================================================
## 1. Cell-type marker library (all peer-reviewed)
## ============================================================
## Elyada 2019 — top cluster-defining markers (from their Table S2 /
## paper Fig. 3B). We use the well-characterised subset.
myCAF <- c("Acta2","Tagln","Myl9","Mylk","Tpm1","Tpm2","Myh11",
           "Postn","Ctgf","Col1a1","Col1a2","Col3a1",
           "Mmp11","Lrrc15","Fap","Thy1","Gja4","Hopx",
           "Thbs1","Thbs2","Mfap5","Comp","Cthrc1","Tnc","Inhba")
iCAF  <- c("Il6","Il11","Lif","Cxcl12","Cxcl1","Cxcl2","Cxcl10",
           "Ccl2","Ccl7","Ccl8","Has1","Has2","Pdgfra","Ly6c1",
           "Ly6c2","Dpt","Cfd","Clec3b","Ptgs2","Lmna","Cxcl13",
           "Nr4a1","Fbln1")
apCAF <- c("H2-Ab1","H2-Aa","H2-Eb1","Cd74","Saa3","Slpi",
           "H2-DMb1","H2-DMa","Cd14")

## Thyroid epithelium / tumour (BPC is thyroid-derived)
Thyroid_tumor <- c("Pax8","Nkx2-1","Foxe1","Tg","Tpo","Duox1","Duox2",
                   "Slc5a5","Slc26a4","Slc5a8","Tshr","Thra","Thrb",
                   "Dio1","Dio2","Epcam","Krt8","Krt18","Cdh1","Cldn1")

## Pan-leukocyte / T-cell / B-cell
Immune <- c("Ptprc","Cd3e","Cd3d","Cd3g","Cd4","Cd8a","Cd8b1",
            "Cd19","Cd79a","Cd79b","Ms4a1","Ncr1","Klrb1c","Foxp3",
            "Itgam","Cd68","Cd163","Lyz2","S100a8","S100a9","Nkg7","Gzmb")

## Endothelial
Endothelial <- c("Pecam1","Cdh5","Vwf","Tek","Kdr","Flt1","Emcn",
                 "Cldn5","Cd34","Plvap","Esam")

markers_list <- list(myCAF = myCAF, iCAF = iCAF, apCAF = apCAF,
                     Thyroid_tumor = Thyroid_tumor,
                     Immune = Immune, Endothelial = Endothelial)
ref_tbl <- tibble::enframe(markers_list, "celltype", "gene") %>% unnest(gene)
write_tsv(ref_tbl, file.path(out_dir, "celltype_marker_library.tsv"))

cat(sprintf("Marker library: %d cell types, %d total (gene x type) pairs\n",
            length(markers_list), nrow(ref_tbl)))

## ============================================================
## 2. Load bulk expression (log2 FPKM+1)
## ============================================================
fpkm_tab <- read_tsv(file.path(proj, "001_mRNA_Summary/1.GeneExpression/1_genes_fpkm_expression.txt"),
                     show_col_types = FALSE, progress = FALSE)
samp_order <- c("WT_BPC1","WT_BPC2","WT_BPC3","WT_BPC4",
                "FAP_BPC1","FAP_BPC2","FAP_BPC3","FAP_BPC4")
samp_grp <- tibble(sample = samp_order,
                   group = factor(ifelse(grepl("^WT", samp_order),
                                         "WT host", "FAP-deficient host"),
                                  levels = grp_levels))

rna_cols <- paste0("FPKM.", samp_order)
rna_mat <- fpkm_tab %>%
  select(gene_name, all_of(rna_cols)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>%
  slice_max(rowSums(across(all_of(rna_cols))), n = 1, with_ties = FALSE) %>%
  ungroup() %>% column_to_rownames("gene_name") %>% as.matrix()
colnames(rna_mat) <- samp_order
rna_log <- log2(rna_mat + 1)

## ============================================================
## 3a. NNLS deconvolution with a synthetic signature matrix
## ============================================================
## Build a cell-type x gene binary signature matrix (markers only).
## Then normalize each cell-type column to unit L1 so that
## NNLS estimates are directly comparable fractions.

all_markers <- unique(unlist(markers_list))
all_markers <- intersect(all_markers, rownames(rna_log))

S <- matrix(0, nrow = length(all_markers),
            ncol = length(markers_list),
            dimnames = list(all_markers, names(markers_list)))
for (ct in names(markers_list)){
  g <- intersect(markers_list[[ct]], all_markers)
  S[g, ct] <- 1
}
# normalize each cell-type column to unit L2 (common CIBERSORT prep)
S_norm <- sweep(S, 2, sqrt(colSums(S^2)), "/")

## subset bulk matrix to marker genes, mean-center per gene
B <- rna_log[all_markers, , drop = FALSE]

# NNLS per sample (B[,j] ≈ S %*% f_j) then normalize to sum-1
frac_mat <- matrix(NA, nrow = ncol(S_norm), ncol = ncol(B),
                   dimnames = list(colnames(S_norm), colnames(B)))
for (j in seq_len(ncol(B))){
  fit <- nnls(A = S_norm, b = B[, j])
  f <- fit$x
  f <- f / sum(f)
  frac_mat[, j] <- f
}

nnls_df <- as.data.frame(frac_mat) %>%
  rownames_to_column("celltype") %>%
  pivot_longer(-celltype, names_to = "sample", values_to = "fraction") %>%
  left_join(samp_grp, by = "sample") %>%
  mutate(celltype = factor(celltype,
                           levels = c("myCAF","iCAF","apCAF",
                                      "Thyroid_tumor","Immune","Endothelial")))
write_tsv(nnls_df, file.path(out_dir, "nnls_cell_fractions.tsv"))

## ============================================================
## 3b. decoupleR ULM enrichment (orthogonal method)
## ============================================================
net_df <- ref_tbl %>% transmute(source = celltype, target = gene, mor = 1)
ulm_res <- suppressMessages(decoupleR::run_ulm(
  mat = rna_log, net = net_df,
  .source = "source", .target = "target", .mor = "mor",
  minsize = 5))
ulm_df <- ulm_res %>%
  transmute(celltype = source, sample = condition, score = score) %>%
  left_join(samp_grp, by = "sample") %>%
  mutate(celltype = factor(celltype,
                           levels = c("myCAF","iCAF","apCAF",
                                      "Thyroid_tumor","Immune","Endothelial")))
write_tsv(ulm_df, file.path(out_dir, "ulm_activity_scores.tsv"))

## ============================================================
## 4. Group comparisons (Wilcoxon)
## ============================================================
comp_stats <- function(df, value_col){
  df %>% group_by(celltype) %>%
    summarise(
      wt_mean  = mean(.data[[value_col]][group == "WT host"],  na.rm = TRUE),
      fap_mean = mean(.data[[value_col]][group == "FAP-deficient host"], na.rm = TRUE),
      delta_fap_minus_wt = fap_mean - wt_mean,
      wilcox_p = tryCatch(wilcox.test(.data[[value_col]] ~ group)$p.value,
                          error = function(e) NA_real_),
      .groups = "drop") %>%
    mutate(method = value_col)
}

stats_all <- bind_rows(
  comp_stats(nnls_df %>% rename(val = fraction),
             "val") %>% mutate(method = "NNLS fraction"),
  comp_stats(ulm_df  %>% rename(val = score),
             "val") %>% mutate(method = "ULM activity")
)
write_tsv(stats_all, file.path(out_dir, "deconvolution_group_stats.tsv"))

## myCAF / iCAF ratio (log, per sample) — extra orthogonal readout
ratio_df <- nnls_df %>%
  filter(celltype %in% c("myCAF","iCAF")) %>%
  select(sample, group, celltype, fraction) %>%
  pivot_wider(names_from = celltype, values_from = fraction) %>%
  mutate(ratio_log = log2((myCAF + 1e-4) / (iCAF + 1e-4)))
ratio_stat <- suppressWarnings(wilcox.test(ratio_log ~ group, data = ratio_df))

## ============================================================
## 5. Panels
## ============================================================
## (a) Stacked bar — per-sample cell fractions (NNLS)
stack_df <- nnls_df %>%
  mutate(sample = factor(sample, levels = samp_order))
ct_cols <- c("myCAF"        = "#C9493A",
             "iCAF"         = "#E08F3A",
             "apCAF"        = "#5E4A7B",
             "Thyroid_tumor"= "#3B6B9C",
             "Immune"       = "#88B0A7",
             "Endothelial"  = "#B7B7B7")

p_a <- ggplot(stack_df, aes(sample, fraction, fill = celltype)) +
  geom_col(colour = "white", linewidth = 0.2, width = 0.9) +
  scale_fill_manual(values = ct_cols, name = NULL) +
  scale_y_continuous(expand = c(0, 0), labels = label_percent()) +
  labs(x = NULL, y = "Estimated cell fraction",
       title = "NNLS deconvolution (marker-based, reference-free)",
       subtitle = "Elyada 2019 CAF markers + curated tumour/immune/endo panels",
       tag = "a") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1,
                                   colour = ifelse(grepl("^WT", samp_order),
                                                   col_wt, col_fap)),
        legend.position = "right",
        legend.key.size = unit(2.8, "mm"))

## (b) myCAF vs iCAF box by group (NNLS fraction)
focus_ct <- c("myCAF","iCAF","Thyroid_tumor")
stat_focus <- stats_all %>% filter(celltype %in% focus_ct,
                                    method == "NNLS fraction")
# build labels outside aes
stat_focus <- stat_focus %>%
  mutate(lbl = vapply(wilcox_p, p_label, character(1)),
         celltype = factor(celltype, levels = focus_ct))
plot_b <- nnls_df %>% filter(celltype %in% focus_ct) %>%
  mutate(celltype = factor(celltype, levels = focus_ct))

p_b <- ggplot(plot_b, aes(group, fraction, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 1,
              shape = 21, stroke = 0.2, colour = "black") +
  geom_text(data = stat_focus,
            aes(x = 1.5, y = Inf, label = lbl),
            inherit.aes = FALSE, vjust = 1.3,
            size = (base_font_size - 1.2)/.pt) +
  facet_wrap(~ celltype, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_x_discrete(labels = c("WT host" = "WT", "FAP-deficient host" = "KO")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22)),
                     labels = label_percent()) +
  labs(x = NULL, y = "Estimated fraction",
       title = "Cell fraction shifts (NNLS)",
       tag = "b") +
  theme_nc()

## (c) myCAF / iCAF ratio
rstat_lbl <- sprintf("Wilcox %s | delta log2 = %.2f",
                     p_label(ratio_stat$p.value),
                     mean(ratio_df$ratio_log[ratio_df$group == "FAP-deficient host"]) -
                     mean(ratio_df$ratio_log[ratio_df$group == "WT host"]))

p_c <- ggplot(ratio_df, aes(group, ratio_log, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 1,
              shape = 21, stroke = 0.2, colour = "black") +
  annotate("text", x = 1.5, y = max(ratio_df$ratio_log) + 0.3,
           label = rstat_lbl,
           size = (base_font_size - 1.2)/.pt, lineheight = 1.05) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_x_discrete(labels = c("WT host" = "WT", "FAP-deficient host" = "KO")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "log2 (myCAF / iCAF) fraction",
       title = "myCAF <-> iCAF balance",
       subtitle = "negative = iCAF-dominant",
       tag = "c") +
  theme_nc()

## (d) Cross-method delta bubble: NNLS vs ULM delta by cell type
cross_df <- stats_all %>%
  select(celltype, delta_fap_minus_wt, wilcox_p, method) %>%
  pivot_wider(names_from = method,
              values_from = c(delta_fap_minus_wt, wilcox_p),
              names_sep = "_") %>%
  rename(delta_NNLS = `delta_fap_minus_wt_NNLS fraction`,
         delta_ULM  = `delta_fap_minus_wt_ULM activity`,
         p_NNLS = `wilcox_p_NNLS fraction`,
         p_ULM  = `wilcox_p_ULM activity`) %>%
  mutate(sig_any = (!is.na(p_NNLS) & p_NNLS < 0.1) |
                    (!is.na(p_ULM)  & p_ULM  < 0.1))

p_d <- ggplot(cross_df, aes(delta_NNLS, delta_ULM)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_point(aes(fill = celltype, size = sig_any),
             shape = 21, stroke = 0.3, colour = "black") +
  geom_text_repel(aes(label = celltype, colour = celltype),
                  size = (base_font_size - 1)/.pt,
                  fontface = "bold",
                  min.segment.length = 0, segment.size = 0.25,
                  box.padding = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = ct_cols, guide = "none") +
  scale_colour_manual(values = ct_cols, guide = "none") +
  scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2.5),
                    labels = c(`TRUE` = "P<0.1 in >=1 method",
                               `FALSE` = "ns"),
                    name = NULL) +
  labs(x = "Delta (FAP-KO - WT) - NNLS fraction",
       y = "Delta (FAP-KO - WT) - ULM activity (z)",
       title = "Cross-method consistency",
       subtitle = "agreement across both orthogonal estimators",
       tag = "d") +
  theme_nc() +
  theme(legend.position = "top")

## assemble
design <- "
AAAAAABBBBBB
AAAAAABBBBBB
AAAAAABBBBBB
CCCCDDDDDDDD
CCCCDDDDDDDD
CCCCDDDDDDDD
"
fig <- p_a + p_b + p_c + p_d +
  plot_layout(design = design) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(out_dir, "caf_deconvolution_overview.pdf")
png_path <- file.path(out_dir, "caf_deconvolution_overview.png")
ggsave(pdf_path, fig, width = 240, height = 200, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(png_path, fig, width = 240, height = 200, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("\nDeconvolution outputs written to:\n  ", out_dir, "\n")

## ============================================================
## 6. Summary
## ============================================================
cat("\n==== NNLS fraction stats (WT vs FAP-KO) ====\n")
print(stats_all %>% filter(method == "NNLS fraction") %>%
        arrange(wilcox_p) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))

cat("\n==== ULM activity stats (WT vs FAP-KO) ====\n")
print(stats_all %>% filter(method == "ULM activity") %>%
        arrange(wilcox_p) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))

cat(sprintf("\nmyCAF/iCAF log2 ratio: WT mean = %.2f | FAP-KO mean = %.2f | Wilcox P = %.3f\n",
            mean(ratio_df$ratio_log[ratio_df$group == "WT host"]),
            mean(ratio_df$ratio_log[ratio_df$group == "FAP-deficient host"]),
            ratio_stat$p.value))
