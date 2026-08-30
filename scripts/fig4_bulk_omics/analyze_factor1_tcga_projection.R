#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S16 - projection of the mouse-derived Factor 1 onto TCGA-THCA.
## Step 1: Project the MOFA2 Factor1 "FAP-KO unbiased signature"
## onto TCGA-THCA bulk RNA, and test association with BRAF vs RAS
## driver status, histology, TDS and official 71-gene BRS.
##
## Data:
##   * MOFA2 Factor1 RNA weights (mouse MGI symbols)
##   * nichenetr mouse->human ortholog converter
##   * TCGA-THCA raw counts (Ensembl IDs, 60660 x 571)
##   * Pre-computed TCGA BRS reference (BRAF/RAS/hist/TDS/BRS)
##
## Method:
##   * Normalize TCGA counts to log2(CPM+1), mean-center features
##   * Build Factor1 signature = top 100 + / top 100 - loading genes
##   * Per-sample score = weighted sum (mean over top genes, z by weight sign)
##   * Complementary ssGSEA-style score using decoupleR::run_ulm
##   * Stratify by BRAF-like / RAS-like class; Wilcoxon test
##   * Histology x score forest (Classical PTC / FV-PTC / Tall cell)
##   * Correlate with TDS, BRS, ERK score


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(org.Hs.eg.db); library(AnnotationDbi)
  library(nichenetr); library(edgeR); library(decoupleR)
})

proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
bulk_out <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
mofa_dir <- file.path(bulk_out, "mofa2")
out_dir  <- file.path(bulk_out, "factor1_tcga_projection")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_wt  <- "#3B6B9C"; col_fap <- "#C9493A"; col_up <- col_fap; col_dn <- col_wt
col_purple <- "#4B1E70"; col_grey <- "#B7B7B7"

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
## 1. Build Factor1 signature (mouse -> human)
## ============================================================
weights <- read_tsv(file.path(mofa_dir, "factor_weights_all.tsv"),
                    show_col_types = FALSE) %>%
  filter(factor == "Factor1", view == "RNA") %>%
  mutate(gene_mouse = sub("_(RNA|Proteomics)$", "", feature))

## Convert mouse -> human symbols
hs_map <- suppressWarnings(convert_mouse_to_human_symbols(weights$gene_mouse))
weights$gene_human <- unname(hs_map)

weights <- weights %>% filter(!is.na(gene_human), gene_human != "") %>%
  distinct(gene_human, .keep_all = TRUE)
cat(sprintf("Factor1 RNA weights translatable to human: %d / %d\n",
            nrow(weights), length(hs_map)))

## Top-K signature (symmetric)
top_k <- 100
sig_up   <- weights %>% arrange(desc(value)) %>% slice_head(n = top_k)
sig_down <- weights %>% arrange(value)       %>% slice_head(n = top_k)
signature_df <- bind_rows(
  sig_up   %>% transmute(gene_human, weight = value, direction = "up_in_FAPKO"),
  sig_down %>% transmute(gene_human, weight = value, direction = "down_in_FAPKO")
)
write_tsv(signature_df, file.path(out_dir, "factor1_signature_top100_human.tsv"))

## ============================================================
## 2. Load TCGA-THCA counts, annotate symbols, normalize
## ============================================================
cat("Loading TCGA-THCA count matrix ...\n")
tcga_counts <- readRDS(file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_public/TCGA-THCA.rds"))
cat(sprintf("Counts: %d x %d\n", nrow(tcga_counts), ncol(tcga_counts)))

## Ensembl version stripping
ens_ids <- sub("\\..*$", "", rownames(tcga_counts))
sym_map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
                                                  keys = unique(ens_ids),
                                                  columns = c("SYMBOL"),
                                                  keytype = "ENSEMBL")) %>%
  as_tibble() %>%
  distinct(ENSEMBL, .keep_all = TRUE)
rowsym <- sym_map$SYMBOL[match(ens_ids, sym_map$ENSEMBL)]
stopifnot(length(rowsym) == nrow(tcga_counts))

## Collapse duplicate symbols by max variance
keep <- !is.na(rowsym)
tcga_counts <- tcga_counts[keep, ]
rowsym <- rowsym[keep]

dup_syms <- rowsym[duplicated(rowsym)] %>% unique()
if (length(dup_syms) > 0){
  var_r <- apply(tcga_counts, 1, var)
  rank_within <- ave(-var_r, rowsym, FUN = rank)
  keep_r <- rank_within == 1
  tcga_counts <- tcga_counts[keep_r, , drop = FALSE]
  rowsym <- rowsym[keep_r]
}
rownames(tcga_counts) <- rowsym

## Filter: keep primary tumors (TCGA sample type 01) only
barcode <- colnames(tcga_counts)
sample_type <- substr(barcode, 14, 15)
tumor_keep <- sample_type %in% c("01")  # primary tumor
tcga_counts <- tcga_counts[, tumor_keep]
cat(sprintf("Tumor samples retained: %d\n", ncol(tcga_counts)))

## Build short sample id (e.g. TCGA-BJ-A0YZ-01) to match BRS ref
short_id <- substr(colnames(tcga_counts), 1, 15)
colnames(tcga_counts) <- short_id

## Normalization: log2(CPM+1) then per-gene z-score
dge <- edgeR::DGEList(counts = tcga_counts)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

## ============================================================
## 3. Score Factor1 signature on TCGA
## ============================================================
## (a) weighted mean using Factor1 weights
sig_expr <- log_cpm[rownames(log_cpm) %in% signature_df$gene_human, , drop = FALSE]
z_mat <- t(scale(t(sig_expr)))                   # gene-wise z
w_vec <- signature_df$weight[match(rownames(z_mat), signature_df$gene_human)]

score_weighted <- colSums(w_vec * z_mat, na.rm = TRUE) / sum(abs(w_vec))

## (b) decoupleR ULM score (using top-100+ vs top-100-)
net_df <- signature_df %>%
  transmute(source = "Factor1_FAPKO",
            target = gene_human,
            mor    = sign(weight))
ulm_res <- decoupleR::run_ulm(mat = log_cpm,
                              net = net_df,
                              .source = "source",
                              .target = "target",
                              .mor    = "mor",
                              minsize = 20)
score_ulm <- ulm_res %>% filter(source == "Factor1_FAPKO") %>%
  transmute(sample = condition, score_ulm = score)

## merge
score_df <- tibble(sample = colnames(log_cpm), score_weighted = score_weighted) %>%
  left_join(score_ulm, by = "sample")

write_tsv(score_df, file.path(out_dir, "factor1_tcga_scores.tsv"))

## ============================================================
## 4. Merge with TCGA BRS reference & driver / histology
## ============================================================
brs_ref <- read_tsv(file.path(bulk_out, "tcga_brs_official/reference/tcga_reference_brs_scores.tsv"),
                    show_col_types = FALSE)
cat("BRS ref samples:", nrow(brs_ref), "\n")

annot <- brs_ref %>%
  transmute(sample = SAMPLE_ID, brs_label, histological_type, mut_driver,
            train_braf, train_ras, brs_raw, tds, erk_score)

merged <- score_df %>% inner_join(annot, by = "sample")
cat(sprintf("Scored + annotated TCGA samples: %d\n", nrow(merged)))

## set analysis groups
merged <- merged %>%
  mutate(brs_group = case_when(
           train_braf == TRUE ~ "BRAF V600E",
           train_ras  == TRUE ~ "RAS driver",
           TRUE ~ "other/none"),
         brs_group = factor(brs_group,
                            levels = c("BRAF V600E", "RAS driver", "other/none")),
         brs_label_f = factor(brs_label, levels = c("Braf-like", "Ras-like")),
         histology = factor(histological_type,
                            levels = c("Classical", "Follicular", "Tall Cell")))

write_tsv(merged, file.path(out_dir, "factor1_tcga_merged_annot.tsv"))

## ============================================================
## 5. Stats
## ============================================================
## BRAF V600E vs RAS driver
braf_ras <- merged %>% filter(brs_group %in% c("BRAF V600E", "RAS driver"))
wt_brafras <- suppressWarnings(wilcox.test(score_weighted ~ brs_group,
                                           data = braf_ras))
## Braf-like vs Ras-like (centroid-assigned)
brs_lab <- merged %>% filter(!is.na(brs_label_f))
wt_bl <- suppressWarnings(wilcox.test(score_weighted ~ brs_label_f,
                                      data = brs_lab))
## correlations
cor_tds  <- suppressWarnings(cor.test(merged$score_weighted, merged$tds,
                                      method = "spearman", exact = FALSE))
cor_brs  <- suppressWarnings(cor.test(merged$score_weighted, merged$brs_raw,
                                      method = "spearman", exact = FALSE))
cor_erk  <- suppressWarnings(cor.test(merged$score_weighted, merged$erk_score,
                                      method = "spearman", exact = FALSE))

cat("\n==== summary ====\n")
cat(sprintf("N (scored & annotated): %d\n", nrow(merged)))
cat(sprintf("BRAF V600E vs RAS driver Wilcoxon P = %.2e\n", wt_brafras$p.value))
cat(sprintf("Braf-like vs Ras-like centroid    P = %.2e\n", wt_bl$p.value))
cat(sprintf("Spearman (score vs TDS):  rho = %.3f, P = %.1e\n",
            cor_tds$estimate, cor_tds$p.value))
cat(sprintf("Spearman (score vs BRS):  rho = %.3f, P = %.1e\n",
            cor_brs$estimate, cor_brs$p.value))
cat(sprintf("Spearman (score vs ERK):  rho = %.3f, P = %.1e\n",
            cor_erk$estimate, cor_erk$p.value))

## ============================================================
## 6. Visualize (4-panel)
## ============================================================
## Panel a — Factor1 score by driver (BRAF/RAS/other)
p_a <- ggplot(merged, aes(brs_group, score_weighted, fill = brs_group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.4,
              shape = 21, stroke = 0.15, colour = "black", alpha = 0.75) +
  annotate("text", x = 1.5, y = max(merged$score_weighted, na.rm = TRUE) + 0.04,
           label = sprintf("BRAF vs RAS %s", p_label(wt_brafras$p.value)),
           size = (base_font_size - 1.5)/.pt) +
  scale_fill_manual(values = c("BRAF V600E" = col_wt,
                               "RAS driver" = col_fap,
                               "other/none" = col_grey),
                    guide = "none") +
  labs(x = NULL, y = "Factor1 weighted score",
       title = "Factor1 score by driver class",
       subtitle = sprintf("TCGA-THCA primary tumors, n = %d", nrow(merged)),
       tag = "a") +
  theme_nc()

## Panel b — score vs TDS scatter (color by driver)
lab_b <- sprintf("Spearman %.2f | %s",
                 cor_tds$estimate, p_label(cor_tds$p.value))
p_b <- ggplot(merged, aes(score_weighted, tds)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_point(aes(fill = brs_group), shape = 21, size = 1.1,
             stroke = 0.2, colour = "black", alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey25", fill = "grey85", linewidth = 0.35) +
  scale_fill_manual(values = c("BRAF V600E" = col_wt,
                               "RAS driver" = col_fap,
                               "other/none" = col_grey),
                    name = NULL) +
  labs(x = "Factor1 weighted score",
       y = "TDS (16-gene thyroid differentiation)",
       title = "Factor1 tracks TDS in TCGA-THCA",
       subtitle = lab_b,
       tag = "b") +
  theme_nc() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"))

## Panel c — score vs BRS scatter
lab_c <- sprintf("Spearman %.2f | %s",
                 cor_brs$estimate, p_label(cor_brs$p.value))
p_c <- ggplot(merged, aes(score_weighted, brs_raw)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_point(aes(fill = brs_group), shape = 21, size = 1.1,
             stroke = 0.2, colour = "black", alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey25", fill = "grey85", linewidth = 0.35) +
  scale_fill_manual(values = c("BRAF V600E" = col_wt,
                               "RAS driver" = col_fap,
                               "other/none" = col_grey),
                    guide = "none") +
  labs(x = "Factor1 weighted score",
       y = "Official 71-gene BRS (raw)",
       title = "Factor1 inversely tracks BRS",
       subtitle = lab_c,
       tag = "c") +
  theme_nc()

## Panel d — histology facet of Factor1 score
hist_stat <- merged %>% filter(!is.na(histology)) %>%
  group_by(histology) %>%
  summarise(n = n(),
            mean_score = mean(score_weighted),
            .groups = "drop") %>%
  mutate(label = sprintf("n=%d", n))

aov_res <- aov(score_weighted ~ histology, data = merged %>% filter(!is.na(histology)))
aov_p <- summary(aov_res)[[1]]$`Pr(>F)`[1]

p_d <- ggplot(merged %>% filter(!is.na(histology)),
              aes(histology, score_weighted, fill = histology)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.4, shape = 21,
              stroke = 0.15, colour = "black", alpha = 0.7) +
  geom_text(data = hist_stat,
            aes(x = histology, y = -Inf, label = label),
            inherit.aes = FALSE, vjust = -0.6,
            size = (base_font_size - 1.5)/.pt) +
  annotate("text", x = 2, y = Inf,
           label = sprintf("ANOVA %s", p_label(aov_p)),
           vjust = 1.5, size = (base_font_size - 1.5)/.pt) +
  scale_fill_manual(values = c("Classical" = "#88B0A7",
                               "Follicular" = "#D9A66A",
                               "Tall Cell" = "#8B5E83"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.20))) +
  labs(x = NULL, y = "Factor1 weighted score",
       title = "Factor1 across TCGA-THCA histology",
       tag = "d") +
  theme_nc()

## Panel e — top signature genes histology split (heatmap top 20)
top_sig <- signature_df %>% arrange(desc(abs(weight))) %>% slice_head(n = 24) %>%
  pull(gene_human) %>% intersect(rownames(log_cpm))

hm_tbl <- log_cpm[top_sig, , drop = FALSE] %>%
  as.data.frame() %>% tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  left_join(merged %>% dplyr::select(sample, brs_group), by = "sample") %>%
  group_by(gene, brs_group) %>%
  summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
  group_by(gene) %>% mutate(z = scale(mean_expr)[, 1]) %>% ungroup() %>%
  left_join(signature_df %>% dplyr::select(gene = gene_human, direction, weight),
            by = "gene") %>%
  mutate(gene = factor(gene,
                       levels = signature_df %>% filter(gene_human %in% top_sig) %>%
                         arrange(weight) %>% pull(gene_human)),
         brs_group = factor(brs_group,
                            levels = c("BRAF V600E", "other/none", "RAS driver")))

p_e <- ggplot(hm_tbl, aes(brs_group, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                       midpoint = 0, oob = squish,
                       limits = c(-1.2, 1.2), name = "z\n(mean)") +
  labs(x = NULL, y = NULL,
       title = "Top Factor1 genes by driver (group z of mean log-CPM)",
       tag = "e") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size = base_font_size - 1.5,
                                   face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key.height = unit(7, "mm"),
        legend.key.width  = unit(2.6, "mm"))

## Panel f — signature weight top +/- bar chart
sig_top20 <- signature_df %>% arrange(desc(abs(weight))) %>%
  slice_head(n = 24) %>%
  mutate(gene_human = factor(gene_human,
                             levels = unique(gene_human[order(weight)])),
         dir = ifelse(weight > 0, "up in FAP-KO", "down in FAP-KO"))

p_f <- ggplot(sig_top20, aes(weight, gene_human, fill = weight)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_col(colour = "black", linewidth = 0.2) +
  scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                       midpoint = 0, guide = "none") +
  labs(x = "Factor1 weight (mouse -> human)",
       y = NULL,
       title = "Signature backbone (top 24 |w|)",
       tag = "f") +
  theme_nc() +
  theme(axis.text.y = element_text(face = "italic",
                                   size = base_font_size - 1.5))

## assemble
design <- "
AAABBB
AAABBB
CCCDDD
CCCDDD
EEEFFF
EEEFFF
"
fig <- p_a + p_b + p_c + p_d + p_e + p_f +
  plot_layout(design = design) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(out_dir, "factor1_tcga_projection_overview.pdf")
png_path <- file.path(out_dir, "factor1_tcga_projection_overview.png")
ggsave(pdf_path, fig, width = 220, height = 260, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(png_path, fig, width = 220, height = 260, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("\nTCGA projection outputs written to:\n  ", out_dir, "\n")

## summary stats table
summary_tbl <- tibble(
  comparison = c("BRAF V600E vs RAS driver (Wilcox)",
                 "Braf-like vs Ras-like centroid (Wilcox)",
                 "Factor1 vs TDS (Spearman)",
                 "Factor1 vs BRS (Spearman)",
                 "Factor1 vs ERK (Spearman)",
                 "Factor1 vs histology (ANOVA)"),
  statistic = c(wt_brafras$statistic, wt_bl$statistic,
                cor_tds$estimate, cor_brs$estimate, cor_erk$estimate,
                aov_p),
  p_value = c(wt_brafras$p.value, wt_bl$p.value,
              cor_tds$p.value, cor_brs$p.value, cor_erk$p.value, aov_p)
)
write_tsv(summary_tbl, file.path(out_dir, "factor1_tcga_summary_stats.tsv"))
