#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S16 - external-cohort validation panels.
## Fig. S9 | External validation of the MOFA2 Factor1 signature in
## TCGA-THCA (BRAF/RAS discrimination, survival) and the Landa GSE76039
## PDTC/ATC cohort (dedifferentiation-grade ladder).


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(survival); library(survminer)
})

proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
bulk_out <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
out_dir  <- file.path(bulk_out, "fig9_external_validation_panel")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- palette ------------------------------------------------------------
col_wt  <- "#3B6B9C"; col_fap <- "#C9493A"; col_up <- col_fap; col_dn <- col_wt
col_grey <- "#B7B7B7"; col_purple <- "#4B1E70"

lad_cols <- c(
  "TCGA RAS-like"  = "#3B6B9C",
  "TCGA other"     = "#B7B7B7",
  "TCGA BRAF-like" = "#5E4A7B",
  "Landa PDTC"     = "#E08F3A",
  "Landa ATC"      = "#C9493A"
)
driver_cols <- c("BRAF V600E" = "#5E4A7B",   # match ladder BRAF-like (purple); consistent driver encoding across panels
                 "RAS driver" = "#3B6B9C",   # match ladder RAS-like (blue)
                 "other/none" = col_grey)

base_font_size <- 7
theme_nc <- function(base = base_font_size){
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(size = base + 1, hjust = 0),
      plot.subtitle = element_text(size = base - 0.5, colour = "grey30"),
      plot.tag      = element_text(face = "bold", size = base + 3, hjust = 0),
      axis.title    = element_text(size = base),
      axis.text     = element_text(size = base - 0.5, colour = "black"),
      axis.line     = element_line(linewidth = 0.6),
      axis.ticks    = element_line(linewidth = 0.6),
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
## Load pre-computed data
## ============================================================
tcga <- read_tsv(file.path(bulk_out,
                           "factor1_tcga_projection/factor1_tcga_merged_annot.tsv"),
                 show_col_types = FALSE)
combined <- read_tsv(file.path(bulk_out,
                               "factor1_gse76039_projection/factor1_combined_scores.tsv"),
                     show_col_types = FALSE)
surv_master <- read_tsv(file.path(bulk_out,
                                  "factor1_tcga_survival/factor1_tcga_survival_master.tsv"),
                        show_col_types = FALSE)
surv_stats <- read_tsv(file.path(bulk_out,
                                 "factor1_tcga_survival/factor1_tcga_survival_stats.tsv"),
                       show_col_types = FALSE)

## lineage gene expression — load from the saved GSE77 run by replaying
## only the lineage part (re-run minimal processing here)
lineage_tbl_path <- file.path(bulk_out, "fig9_external_validation_panel/lineage_z_cache.tsv")
if (!file.exists(lineage_tbl_path)){
  suppressPackageStartupMessages({
    library(edgeR); library(org.Hs.eg.db); library(AnnotationDbi)
  })
  tcga_counts <- readRDS(file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_public/TCGA-THCA.rds"))
  ens_ids <- sub("\\..*$", "", rownames(tcga_counts))
  sym_map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
                                                    keys = unique(ens_ids),
                                                    columns = "SYMBOL",
                                                    keytype = "ENSEMBL")) %>%
    as_tibble() %>% distinct(ENSEMBL, .keep_all = TRUE)
  rowsym <- sym_map$SYMBOL[match(ens_ids, sym_map$ENSEMBL)]
  keep <- !is.na(rowsym)
  tcga_counts <- tcga_counts[keep, ]; rowsym <- rowsym[keep]
  dup_syms <- rowsym[duplicated(rowsym)] %>% unique()
  if (length(dup_syms) > 0){
    var_r <- apply(tcga_counts, 1, var)
    rank_within <- ave(-var_r, rowsym, FUN = rank)
    tcga_counts <- tcga_counts[rank_within == 1, , drop = FALSE]
    rowsym <- rowsym[rank_within == 1]
  }
  rownames(tcga_counts) <- rowsym
  sample_type <- substr(colnames(tcga_counts), 14, 15)
  tcga_counts <- tcga_counts[, sample_type == "01"]
  colnames(tcga_counts) <- substr(colnames(tcga_counts), 1, 15)
  dge <- edgeR::DGEList(counts = tcga_counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  tcga_logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  gse_exp <- read.csv(file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_public/GSE76039_EXP.csv"),
                      row.names = 1, check.names = FALSE)
  gse_mat <- as.matrix(gse_exp)

  lineage_genes <- c("PAX8", "NKX2-1", "TG", "TPO", "SLC5A5", "FOXE1", "DUOX2")
  lineage_all <- intersect(intersect(lineage_genes, rownames(tcga_logcpm)),
                           rownames(gse_mat))
  tcga_ann <- tcga %>%
    transmute(sample, ladder = case_when(
      brs_group == "BRAF V600E" ~ "TCGA BRAF-like",
      brs_group == "RAS driver" ~ "TCGA RAS-like",
      TRUE ~ "TCGA other"))
  gse_score <- read_tsv(file.path(bulk_out,
                                  "factor1_gse76039_projection/factor1_gse76039_scores.tsv"),
                        show_col_types = FALSE)
  tcga_z <- t(scale(t(tcga_logcpm[lineage_all, intersect(colnames(tcga_logcpm), tcga_ann$sample), drop = FALSE])))
  gse_z  <- t(scale(t(gse_mat[lineage_all, , drop = FALSE])))
  lineage_df <- bind_rows(
    as.data.frame(tcga_z) %>% rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "sample", values_to = "z") %>%
      left_join(tcga_ann, by = "sample"),
    as.data.frame(gse_z) %>% rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "sample", values_to = "z") %>%
      left_join(gse_score %>% transmute(sample, ladder = paste0("Landa ", as.character(class))),
                by = "sample")
  ) %>% filter(!is.na(ladder))
  write_tsv(lineage_df, lineage_tbl_path)
} else {
  lineage_df <- read_tsv(lineage_tbl_path, show_col_types = FALSE)
}

lineage_df <- lineage_df %>%
  mutate(ladder = factor(ladder,
                         levels = c("TCGA RAS-like", "TCGA other",
                                    "TCGA BRAF-like", "Landa PDTC",
                                    "Landa ATC")),
         gene = factor(gene, levels = c("PAX8","NKX2-1","TG","TPO",
                                        "SLC5A5","FOXE1","DUOX2")))

## ============================================================
## Panel a — TCGA Factor1 score by driver class
## ============================================================
tcga <- tcga %>%
  mutate(brs_group = factor(brs_group,
                            levels = c("BRAF V600E","RAS driver","other/none")))
braf_ras <- tcga %>% filter(brs_group %in% c("BRAF V600E","RAS driver"))
wt_br <- suppressWarnings(wilcox.test(score_weighted ~ brs_group,
                                      data = braf_ras))

p_a <- ggplot(tcga, aes(brs_group, score_weighted, fill = brs_group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.4,
              shape = 21, stroke = 0.15, colour = "black", alpha = 0.7) +
  annotate("text", x = 1.5, y = max(tcga$score_weighted) + 0.03,
           label = sprintf("BRAF vs RAS %s", p_label(wt_br$p.value)),
           size = (base_font_size - 1.2)/.pt) +
  scale_fill_manual(values = driver_cols, guide = "none") +
  labs(x = NULL, y = "Factor1 signature score",
       title = "TCGA-THCA by driver class",
       subtitle = sprintf("n = %d primary tumors", nrow(tcga)),
       tag = "a") +
  theme_nc()

## ============================================================
## Panel b — Factor1 vs BRS scatter
## ============================================================
cor_brs <- suppressWarnings(cor.test(tcga$score_weighted, tcga$brs_raw,
                                     method = "spearman", exact = FALSE))
lab_b <- sprintf("Spearman %.2f | %s",
                 cor_brs$estimate, p_label(cor_brs$p.value))
p_b <- ggplot(tcga, aes(score_weighted, brs_raw)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_point(aes(fill = brs_group), shape = 21, size = 1.1,
             stroke = 0.18, colour = "black", alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey25", fill = "grey85", linewidth = 0.35) +
  scale_fill_manual(values = driver_cols, name = NULL) +
  labs(x = "Factor1 signature score",
       y = "Official 71-gene BRS (raw)",
       title = "Factor1 inversely tracks BRS",
       subtitle = lab_b,
       tag = "b") +
  theme_nc() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"))

## ============================================================
## Panel c — 5-tier dedifferentiation ladder
## ============================================================
ladder <- combined %>%
  mutate(ladder = case_when(
    cohort == "TCGA-THCA (DTC)" & subgroup == "TCGA RAS"   ~ "TCGA RAS-like",
    cohort == "TCGA-THCA (DTC)" & subgroup == "TCGA other" ~ "TCGA other",
    cohort == "TCGA-THCA (DTC)" & subgroup == "TCGA BRAF"  ~ "TCGA BRAF-like",
    cohort == "GSE76039 (PDTC)" ~ "Landa PDTC",
    cohort == "GSE76039 (ATC)"  ~ "Landa ATC"
  ),
  ladder = factor(ladder,
                  levels = c("TCGA RAS-like", "TCGA other",
                             "TCGA BRAF-like", "Landa PDTC", "Landa ATC")))

aov_ld <- aov(score_weighted ~ ladder, data = ladder)
aov_p_ld <- summary(aov_ld)[[1]]$`Pr(>F)`[1]
kt <- suppressWarnings(cor.test(as.numeric(ladder$ladder),
                                ladder$score_weighted,
                                method = "kendall"))

p_c <- ggplot(ladder, aes(ladder, score_weighted, fill = ladder)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.5,
              shape = 21, stroke = 0.15, colour = "black", alpha = 0.7) +
  annotate("text", x = 3, y = Inf,
           label = sprintf("ANOVA %s | Kendall tau %.2f (%s)",
                           p_label(aov_p_ld),
                           kt$estimate, p_label(kt$p.value)),
           vjust = 1.4, size = (base_font_size - 1.2)/.pt) +
  scale_fill_manual(values = lad_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.14))) +
  labs(x = NULL, y = "Factor1 signature score",
       title = "Dedifferentiation-grade ladder",
       subtitle = "TCGA-THCA split by driver + Landa PDTC/ATC",
       tag = "c") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1,
                                   size = base_font_size - 0.5))

## ============================================================
## Panel d — Thyroid lineage genes across ladder
## ============================================================
p_d <- ggplot(lineage_df, aes(ladder, z, fill = ladder)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.25, alpha = 0.85) +
  facet_wrap(~ gene, nrow = 1) +
  scale_fill_manual(values = lad_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(x = NULL, y = "z (per cohort)",
       title = "Thyroid lineage genes collapse along the ladder",
       tag = "d") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1,
                                   size = base_font_size - 1.5),
        panel.spacing.x = unit(2, "mm"))

## ============================================================
## Panel e — TCGA PFI KM by Factor1 tertile
## ============================================================
master <- surv_master %>%
  filter(!is.na(PFI), !is.na(PFI.time)) %>%
  mutate(score_tertile = factor(score_tertile, levels = c("low","mid","high")))

sfit_pfi <- survfit(Surv(PFI.time, PFI) ~ score_tertile, data = master)
diff_pfi <- survdiff(Surv(PFI.time, PFI) ~ score_tertile, data = master)
lr_pfi <- 1 - pchisq(diff_pfi$chisq, df = length(diff_pfi$n) - 1)

km_pfi <- ggsurvplot(
  sfit_pfi, data = master,
  risk.table = TRUE, conf.int = FALSE,
  palette = c(col_wt, "grey60", col_fap),
  legend.labs = c("low","mid","high"),
  legend.title = "Factor1 tertile",
  pval = sprintf("Log-rank %s", p_label(lr_pfi)),
  pval.size = 2.4,
  risk.table.height = 0.28,
  ggtheme = theme_nc() +
    theme(legend.position = "top", legend.direction = "horizontal"),
  font.main = 7, font.x = 6, font.y = 6,
  font.tickslab = 6, font.legend = 6,
  tables.theme = theme_nc() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6))
)
km_pfi$plot <- km_pfi$plot +
  labs(x = "PFI time (days)",
       y = "Progression-free probability",
       title = "TCGA-THCA PFI by Factor1 tertile",
       tag = "e")
km_grid <- cowplot::plot_grid(km_pfi$plot, km_pfi$table, ncol = 1,
                              rel_heights = c(3, 1.2))

## ============================================================
## Panel f — Cox HR forest (OS / DSS / PFI; univariate + adj)
## ============================================================
forest_tbl <- surv_stats %>%
  transmute(outcome, n, events,
            HR_uni = cox_HR, lo_uni = cox_HR_lo, hi_uni = cox_HR_hi,
            p_uni  = cox_p,
            HR_mv  = cox_mv_HR_score, p_mv = cox_mv_p_score) %>%
  pivot_longer(cols = -c(outcome, n, events),
               names_to = c(".value", "model"),
               names_pattern = "([a-zA-Z_]+)_(uni|mv)") %>%
  mutate(model = recode(model, uni = "univariate",
                        mv = "age + stage adj"),
         model = factor(model, levels = c("univariate", "age + stage adj")),
         row_lbl = paste0(outcome, " | ", model),
         row_lbl = factor(row_lbl, levels = rev(unique(row_lbl))),
         p_str = vapply(p, p_label, character(1)),
         panel_lbl = sprintf("HR %.2f (%.2f-%.2f), %s",
                             HR, lo, hi, p_str))

p_f <- ggplot(forest_tbl, aes(HR, row_lbl)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             linewidth = 0.3, colour = "grey45") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.4,
                colour = "grey25", orientation = "y") +
  geom_point(aes(fill = log2(HR)), shape = 21, size = 2.2,
             stroke = 0.3, colour = "black") +
  geom_text(aes(x = hi, label = panel_lbl),
            hjust = -0.08, size = (base_font_size - 1.5)/.pt) +
  scale_fill_gradient2(low = col_wt, mid = "#F5F5F5", high = col_fap,
                       midpoint = 0, guide = "none") +
  scale_x_log10(labels = function(x) format(x, digits = 2),
                expand = expansion(mult = c(0.05, 0.45))) +
  labs(x = "Cox HR per 1 unit Factor1 score (log scale)",
       y = NULL,
       title = "Factor1 - survival in TCGA-THCA",
       subtitle = "PFI most informative; OS/DSS events few",
       tag = "f") +
  theme_nc()

## ============================================================
## assemble
## ============================================================
design <- "
AABBCC
AABBCC
DDDDDD
DDDDDD
EEFFFF
EEFFFF
EEFFFF
"

## Panel e is km_grid (cowplot) — need to wrap.
## We'll use wrap_plots with patchwork-compatible grob via cowplot::plot_grid().
km_wrap <- patchwork::wrap_elements(full = km_grid)

fig <- p_a + p_b + p_c + p_d + km_wrap + p_f +
  plot_layout(design = design) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(out_dir, "FigS9_factor1_external_validation.pdf")
png_path <- file.path(out_dir, "FigS9_factor1_external_validation.png")

ggsave(pdf_path, fig, width = 220, height = 280, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(png_path, fig, width = 220, height = 280, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("Fig S9 written:\n  ", pdf_path, "\n  ", png_path, "\n")
