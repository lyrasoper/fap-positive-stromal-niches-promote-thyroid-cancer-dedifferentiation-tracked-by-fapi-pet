#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S16 - projection of Factor 1 onto GSE76039.
## Step 2: Project the MOFA2 Factor1 signature onto the Landa et al.
## GSE76039 cohort of poorly-differentiated (PDTC) and anaplastic (ATC)
## thyroid carcinomas, and combine with TCGA-THCA (DTC) for a
## dedifferentiation-grade ladder.
##
## Data:
##   * GSE76039_EXP.csv  (HG-U133 Plus 2, gcRMA-normalized, gene symbol)
##   * GSE76039_cli.csv  (20 ATC + 17 PDTC)
##   * MOFA2 Factor1 signature (top 100+ / top 100- human symbols)
##   * TCGA-THCA scored output from the previous step (for ladder)


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(decoupleR)
})

proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
bulk_out <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
sig_path <- file.path(bulk_out, "factor1_tcga_projection/factor1_signature_top100_human.tsv")
tcga_scored_path <- file.path(bulk_out, "factor1_tcga_projection/factor1_tcga_merged_annot.tsv")
out_dir <- file.path(bulk_out, "factor1_gse76039_projection")
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
## 1. Load signature & GSE76039
## ============================================================
signature_df <- read_tsv(sig_path, show_col_types = FALSE)

gse_exp <- read.csv(file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_public/GSE76039_EXP.csv"),
                    row.names = 1, check.names = FALSE)
gse_mat <- as.matrix(gse_exp)
cat(sprintf("GSE76039 expression: %d genes x %d samples\n",
            nrow(gse_mat), ncol(gse_mat)))

gse_cli <- read.csv(file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_public/GSE76039_cli.csv"),
                    stringsAsFactors = FALSE)
gse_cli <- gse_cli %>%
  transmute(sample = geo_accession,
            source = source_name_ch1,
            tumor_type = .data[["tumor.type.ch1"]],
            gender = .data[["gender.ch1"]],
            tissue = .data[["tissue.ch1"]]) %>%
  mutate(class = case_when(
           grepl("^Anaplastic", source) ~ "ATC",
           grepl("^Poorly",     source) ~ "PDTC",
           TRUE ~ "other"),
         class = factor(class, levels = c("PDTC", "ATC")))
stopifnot(all(gse_cli$sample %in% colnames(gse_mat)))
gse_mat <- gse_mat[, gse_cli$sample, drop = FALSE]

## ============================================================
## 2. Factor1 score on GSE76039
## ============================================================
common_sig <- intersect(signature_df$gene_human, rownames(gse_mat))
cat(sprintf("Signature genes present in GSE76039: %d / %d\n",
            length(common_sig), nrow(signature_df)))

sig_expr <- gse_mat[common_sig, , drop = FALSE]
z_mat <- t(scale(t(sig_expr)))
w_vec <- signature_df$weight[match(rownames(z_mat), signature_df$gene_human)]
score_weighted <- colSums(w_vec * z_mat, na.rm = TRUE) / sum(abs(w_vec))

net_df <- signature_df %>%
  transmute(source = "Factor1_FAPKO",
            target = gene_human,
            mor    = sign(weight))
ulm_res <- decoupleR::run_ulm(mat = gse_mat, net = net_df,
                              .source = "source", .target = "target",
                              .mor = "mor", minsize = 20)
score_ulm <- ulm_res %>% filter(source == "Factor1_FAPKO") %>%
  transmute(sample = condition, score_ulm = score)

gse_score <- tibble(sample = colnames(gse_mat),
                    score_weighted = score_weighted) %>%
  left_join(score_ulm, by = "sample") %>%
  left_join(gse_cli, by = "sample")

write_tsv(gse_score, file.path(out_dir, "factor1_gse76039_scores.tsv"))

## ============================================================
## 3. Stats within GSE76039 and combined with TCGA-THCA
## ============================================================
wt_pa <- suppressWarnings(wilcox.test(score_weighted ~ class,
                                      data = gse_score))
cat(sprintf("\nPDTC vs ATC Wilcoxon P = %.3e\n", wt_pa$p.value))

tcga_scored <- read_tsv(tcga_scored_path, show_col_types = FALSE) %>%
  transmute(sample, score_weighted, cohort = "TCGA-THCA (DTC)",
            subgroup = case_when(
              brs_group == "BRAF V600E" ~ "TCGA BRAF",
              brs_group == "RAS driver" ~ "TCGA RAS",
              TRUE ~ "TCGA other"
            ))
gse_long <- gse_score %>%
  transmute(sample, score_weighted,
            cohort = paste0("GSE76039 (", class, ")"),
            subgroup = as.character(class))
combined <- bind_rows(tcga_scored, gse_long) %>%
  mutate(cohort = factor(cohort,
                         levels = c("TCGA-THCA (DTC)",
                                    "GSE76039 (PDTC)",
                                    "GSE76039 (ATC)")),
         subgroup = factor(subgroup,
                           levels = c("TCGA RAS","TCGA other","TCGA BRAF",
                                      "PDTC","ATC")))

## TCGA (split by RAS vs BRAF) + PDTC + ATC as 5-level ladder
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

## ANOVA across ladder
aov_ld <- aov(score_weighted ~ ladder, data = ladder)
aov_p_ld <- summary(aov_ld)[[1]]$`Pr(>F)`[1]

## Pairwise Wilcoxon with BH
pw <- pairwise.wilcox.test(ladder$score_weighted, ladder$ladder,
                           p.adjust.method = "BH")

## trend test (Jonckheere-Terpstra): use Kendall tau across ordered ladder
kt <- suppressWarnings(cor.test(as.numeric(ladder$ladder),
                                ladder$score_weighted,
                                method = "kendall"))

write_tsv(combined, file.path(out_dir, "factor1_combined_scores.tsv"))
write_tsv(as_tibble(pw$p.value, rownames = "grp1") %>%
            pivot_longer(-grp1, names_to = "grp2", values_to = "padj") %>%
            filter(!is.na(padj)),
          file.path(out_dir, "factor1_ladder_pairwise_padj.tsv"))

## ============================================================
## 4. Visualization
## ============================================================
lad_cols <- c(
  "TCGA RAS-like"  = "#3B6B9C",
  "TCGA other"     = "#B7B7B7",
  "TCGA BRAF-like" = "#5E4A7B",
  "Landa PDTC"     = "#E08F3A",
  "Landa ATC"      = "#C9493A"
)

## Panel a — ladder boxplot (DTC → PDTC → ATC)
p_a <- ggplot(ladder, aes(ladder, score_weighted, fill = ladder)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.55,
              shape = 21, stroke = 0.15, colour = "black", alpha = 0.7) +
  annotate("text", x = 3, y = Inf,
           label = sprintf("ANOVA %s | Kendall tau %.2f (%s)",
                           p_label(aov_p_ld),
                           kt$estimate, p_label(kt$p.value)),
           vjust = 1.5, size = (base_font_size - 1)/.pt) +
  scale_fill_manual(values = lad_cols, guide = "none") +
  labs(x = NULL, y = "Factor1 signature score",
       title = "Factor1 score across dedifferentiation grade",
       subtitle = "TCGA-THCA (DTC) split by driver, Landa PDTC/ATC",
       tag = "a") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1,
                                   size = base_font_size - 0.5))

## Panel b — PDTC vs ATC detail
p_b <- ggplot(gse_score, aes(class, score_weighted, fill = class)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.8,
              shape = 21, stroke = 0.2, colour = "black") +
  annotate("text", x = 1.5, y = max(gse_score$score_weighted) + 0.02,
           label = sprintf("PDTC vs ATC %s", p_label(wt_pa$p.value)),
           size = (base_font_size - 1)/.pt) +
  scale_fill_manual(values = c("PDTC" = "#E08F3A", "ATC" = "#C9493A"),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL, y = "Factor1 signature score",
       title = "GSE76039 Landa cohort",
       subtitle = sprintf("n = %d PDTC, n = %d ATC",
                          sum(gse_score$class == "PDTC"),
                          sum(gse_score$class == "ATC")),
       tag = "b") +
  theme_nc()

## Panel c — Factor1 signature gene heatmap across ladder (top 25 genes,
##            rows scaled across all samples; columns ordered by ladder)
top_sig <- signature_df %>% arrange(desc(abs(weight))) %>%
  slice_head(n = 25) %>% pull(gene_human)

## Combined expression matrix: standardize per gene and per cohort,
## then align on shared genes.
# Build a long df of (sample, gene, z_expr) per cohort, then concat
# TCGA expression matrix reload
tcga_counts <- readRDS(file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_public/TCGA-THCA.rds"))
library(edgeR); library(org.Hs.eg.db); library(AnnotationDbi)
ens_ids <- sub("\\..*$", "", rownames(tcga_counts))
sym_map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
                                                  keys = unique(ens_ids),
                                                  columns = c("SYMBOL"),
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
barcode <- colnames(tcga_counts); short_id <- substr(barcode, 1, 15)
sample_type <- substr(barcode, 14, 15)
tcga_counts <- tcga_counts[, sample_type == "01"]
colnames(tcga_counts) <- substr(colnames(tcga_counts), 1, 15)
dge <- edgeR::DGEList(counts = tcga_counts)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
tcga_logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

# Common genes
tcga_ann <- read_tsv(tcga_scored_path, show_col_types = FALSE) %>%
  transmute(sample, ladder = case_when(
    brs_group == "BRAF V600E" ~ "TCGA BRAF-like",
    brs_group == "RAS driver" ~ "TCGA RAS-like",
    TRUE ~ "TCGA other"
  ))
tcga_samples <- intersect(colnames(tcga_logcpm), tcga_ann$sample)

# Select matched top_sig genes
g_tcga <- intersect(top_sig, rownames(tcga_logcpm))
g_gse  <- intersect(top_sig, rownames(gse_mat))
g_all  <- intersect(g_tcga, g_gse)

# Build per-sample z_expr (standardize per gene WITHIN each cohort)
tcga_sub <- tcga_logcpm[g_all, tcga_samples, drop = FALSE]
gse_sub  <- gse_mat[g_all, , drop = FALSE]
tcga_z <- t(scale(t(tcga_sub)))
gse_z  <- t(scale(t(gse_sub)))
z_all  <- cbind(tcga_z, gse_z)

long_z <- as.data.frame(z_all) %>% rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "z") %>%
  left_join(bind_rows(tcga_ann,
                      gse_score %>%
                        transmute(sample, ladder = paste0("Landa ", as.character(class)))),
            by = "sample") %>%
  filter(!is.na(ladder)) %>%
  mutate(ladder = factor(ladder,
                         levels = c("TCGA RAS-like","TCGA other",
                                    "TCGA BRAF-like","Landa PDTC","Landa ATC"))) %>%
  group_by(gene, ladder) %>% summarise(z_mean = mean(z, na.rm = TRUE), .groups = "drop") %>%
  left_join(signature_df %>% dplyr::select(gene = gene_human, weight),
            by = "gene") %>%
  mutate(gene = factor(gene,
                       levels = signature_df %>% filter(gene_human %in% g_all) %>%
                         arrange(weight) %>% pull(gene_human)))

p_c <- ggplot(long_z, aes(ladder, gene, fill = z_mean)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                       midpoint = 0, oob = squish,
                       limits = c(-1.1, 1.1), name = "z\n(cohort)") +
  labs(x = NULL, y = NULL,
       title = "Top Factor1 genes across ladder",
       tag = "c") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1,
                                   size = base_font_size - 0.5),
        axis.text.y = element_text(size = base_font_size - 1.5,
                                   face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key.height = unit(6, "mm"),
        legend.key.width  = unit(2.6, "mm"))

## Panel d — PAX8 / TG / TPO / NIS mean across ladder (redifferentiation readout)
lineage_genes <- c("PAX8", "NKX2-1", "TG", "TPO", "SLC5A5", "FOXE1", "DUOX2")
lineage_tcga <- intersect(lineage_genes, rownames(tcga_logcpm))
lineage_gse  <- intersect(lineage_genes, rownames(gse_mat))
lineage_all  <- intersect(lineage_tcga, lineage_gse)

lineage_z <- cbind(
  t(scale(t(tcga_logcpm[lineage_all, tcga_samples, drop = FALSE]))),
  t(scale(t(gse_mat[lineage_all, , drop = FALSE])))
)
lin_long <- as.data.frame(lineage_z) %>% rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "z") %>%
  left_join(bind_rows(tcga_ann,
                      gse_score %>% transmute(sample,
                                              ladder = paste0("Landa ", as.character(class)))),
            by = "sample") %>%
  filter(!is.na(ladder)) %>%
  mutate(ladder = factor(ladder,
                         levels = c("TCGA RAS-like","TCGA other",
                                    "TCGA BRAF-like","Landa PDTC","Landa ATC")),
         gene = factor(gene, levels = lineage_all))

p_d <- ggplot(lin_long, aes(ladder, z, fill = ladder)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.25, alpha = 0.85) +
  facet_wrap(~ gene, nrow = 1) +
  scale_fill_manual(values = lad_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(x = NULL, y = "z (per cohort)",
       title = "Thyroid lineage genes across ladder",
       tag = "d") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1,
                                   size = base_font_size - 1.2),
        panel.spacing.x = unit(2, "mm"))

## assemble
design <- "
AABB
AABB
CCCC
CCCC
DDDD
"
fig <- p_a + p_b + p_c + p_d +
  plot_layout(design = design) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(out_dir, "factor1_gse76039_ladder_overview.pdf")
png_path <- file.path(out_dir, "factor1_gse76039_ladder_overview.png")
ggsave(pdf_path, fig, width = 220, height = 240, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(png_path, fig, width = 220, height = 240, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("\nGSE76039 + ladder outputs written to:\n  ", out_dir, "\n")

## summary block
cat("\n==== summary ====\n")
cat(sprintf("GSE76039: %d genes in signature overlap\n", length(common_sig)))
cat(sprintf("PDTC vs ATC Wilcoxon P = %.3e\n", wt_pa$p.value))
cat(sprintf("Ladder ANOVA P = %.3e\n", aov_p_ld))
cat(sprintf("Ladder Kendall tau = %.2f (P = %.3e)\n",
            kt$estimate, kt$p.value))
cat("Pairwise adj p-values (BH):\n")
print(pw$p.value)
