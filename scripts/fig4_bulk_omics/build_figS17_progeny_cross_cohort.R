#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S20 - cross-cohort PROGENy comparison.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.
## Fig. S17 — PROGENy 14-pathway cross-cohort activity.
##
## Purpose: provide an orthogonal, peer-reviewed pathway-level
## overview across
##   (i)   Fig 8 mouse WT vs FAP-KO
##   (ii)  TCGA-THCA n=481 (BRAF / RAS / other)
##   (iii) Landa GSE76039 PDTC / ATC n=37
##
## Method: decoupleR::run_mlm with PROGENy top-500 weighted targets,
## run separately per organism (human / mouse) and per dataset.
##
## Panels:
##   a  Fig 8 mouse — 14-pathway activity delta (FAP-KO - WT)
##   b  TCGA-THCA — pathway activity BRAF vs RAS (delta + P)
##   c  Landa GSE76039 — PDTC vs ATC delta
##   d  Cross-cohort consistency heatmap (3 datasets x 14 pathways)
##   e  Factor1 x pathway activity correlation (TCGA)
##   f  5-tier dedifferentiation ladder for MAPK/WNT/TGFb/JAK-STAT


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(decoupleR)
})

proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
out_root <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
fig_dir <- file.path(out_root, "fig_S17_progeny_cross_cohort")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

col_wt   <- "#3B6B9C"; col_fap  <- "#C9493A"; col_grey <- "#B7B7B7"
col_up   <- col_fap;   col_dn   <- col_wt
col_purple <- "#4B1E70"
grp_levels <- c("WT host", "FAP-deficient host")
grp_cols <- setNames(c(col_wt, col_fap), grp_levels)

base_font_size <- 7
theme_nc <- function(base = base_font_size){
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(size = base + 1, hjust = 0, face = "plain"),
      plot.subtitle = element_text(size = base - 0.5, colour = "grey30"),
      plot.tag      = element_text(face = "bold", size = base + 3, hjust = 0),
      axis.title    = element_text(size = base),
      axis.text     = element_text(size = base - 0.5, colour = "black"),
      axis.line     = element_line(linewidth = 0.3),
      axis.ticks    = element_line(linewidth = 0.3),
      strip.background = element_blank(),
      strip.text    = element_text(size = base, face = "plain"),
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
## 1. Load PROGENy networks (human + mouse)
## ============================================================
cat("Loading PROGENy networks ...\n")
net_h <- suppressMessages(decoupleR::get_progeny(organism = "human", top = 500))
net_m <- suppressMessages(decoupleR::get_progeny(organism = "mouse", top = 500))

## pathway order
pw_levels <- c("MAPK","PI3K","EGFR","VEGF","JAK-STAT","NFkB","TNFa",
               "TGFb","WNT","Hypoxia","p53","Trail","Estrogen","Androgen")
## display relabel only (data keys/factor levels stay as-is)
pw_lab <- function(x){ x <- gsub("NFkB","NF-κB",x); x <- gsub("TNFa","TNF-α",x); x <- gsub("TGFb","TGF-β",x); x }

## ============================================================
## 2. Score Fig 8 mouse bulk RNA (WT vs FAP-KO)
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
  dplyr::select(gene_name, all_of(rna_cols)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>%
  slice_max(rowSums(across(all_of(rna_cols))), n = 1, with_ties = FALSE) %>%
  ungroup() %>% column_to_rownames("gene_name") %>% as.matrix()
colnames(rna_mat) <- samp_order
rna_log <- log2(rna_mat + 1)

cat("Running PROGENy ULM on mouse ...\n")
mouse_ulm <- suppressMessages(decoupleR::run_mlm(
  mat = rna_log, net = net_m,
  .source = "source", .target = "target", .mor = "weight",
  minsize = 10))

mouse_scores <- mouse_ulm %>%
  transmute(sample = condition, pathway = source, score = score) %>%
  left_join(samp_grp, by = "sample")

mouse_stats <- mouse_scores %>% group_by(pathway) %>%
  summarise(wt_mean = mean(score[group == "WT host"]),
            fap_mean = mean(score[group == "FAP-deficient host"]),
            delta_fap_minus_wt = fap_mean - wt_mean,
            wilcox_p = tryCatch(wilcox.test(score ~ group)$p.value,
                                error = function(e) NA_real_),
            .groups = "drop") %>%
  mutate(pathway = factor(pathway, levels = pw_levels))

write_tsv(mouse_stats, file.path(fig_dir, "mouse_progeny_stats.tsv"))

## ============================================================
## 3. Score TCGA-THCA
## ============================================================
suppressPackageStartupMessages({
  library(edgeR); library(org.Hs.eg.db); library(AnnotationDbi)
})
cat("Loading TCGA-THCA ...\n")
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
tcga_counts <- tcga_counts[, substr(colnames(tcga_counts), 14, 15) == "01"]
colnames(tcga_counts) <- substr(colnames(tcga_counts), 1, 15)
dge <- edgeR::DGEList(counts = tcga_counts)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
tcga_log <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

cat("Running PROGENy mlm on TCGA ...\n")
tcga_ulm <- suppressMessages(decoupleR::run_mlm(
  mat = tcga_log, net = net_h,
  .source = "source", .target = "target", .mor = "weight",
  minsize = 10))

tcga_scores <- tcga_ulm %>%
  transmute(sample = condition, pathway = source, score = score)

## merge with TCGA annotation + Factor1
tcga_master <- read_tsv(file.path(out_root,
                                  "factor1_tcga_projection/factor1_tcga_merged_annot.tsv"),
                        show_col_types = FALSE) %>%
  dplyr::select(sample, score_weighted, brs_group, tds, brs_raw)
tcga_wide <- tcga_scores %>%
  pivot_wider(names_from = pathway, values_from = score) %>%
  inner_join(tcga_master, by = "sample")
cat(sprintf("TCGA merged: %d samples\n", nrow(tcga_wide)))

## BRAF vs RAS stats
tcga_bra <- tcga_wide %>% filter(brs_group %in% c("BRAF V600E","RAS driver"))
tcga_stats <- tibble(pathway = pw_levels) %>%
  rowwise() %>%
  mutate(wt_like_mean_BRAF = mean(tcga_bra[[pathway]][tcga_bra$brs_group == "BRAF V600E"]),
         wt_like_mean_RAS  = mean(tcga_bra[[pathway]][tcga_bra$brs_group == "RAS driver"]),
         delta_BRAF_minus_RAS = wt_like_mean_BRAF - wt_like_mean_RAS,
         wilcox_p = tryCatch(
           wilcox.test(tcga_bra[[pathway]] ~ tcga_bra$brs_group)$p.value,
           error = function(e) NA_real_)) %>%
  ungroup() %>%
  mutate(pathway = factor(pathway, levels = pw_levels))

write_tsv(tcga_stats, file.path(fig_dir, "tcga_progeny_stats.tsv"))

## Factor1 correlations per pathway
factor1_cor <- tibble(pathway = pw_levels) %>%
  rowwise() %>%
  mutate(rho = cor(tcga_wide[[pathway]], tcga_wide$score_weighted,
                    method = "spearman", use = "pairwise.complete.obs"),
         p = cor.test(tcga_wide[[pathway]], tcga_wide$score_weighted,
                     method = "spearman", exact = FALSE,
                     use = "pairwise.complete.obs")$p.value) %>%
  ungroup() %>%
  mutate(pathway = factor(pathway, levels = pw_levels))
write_tsv(factor1_cor, file.path(fig_dir, "tcga_factor1_pathway_correlations.tsv"))

## ============================================================
## 4. Score Landa GSE76039
## ============================================================
gse_exp <- read.csv(file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_public/GSE76039_EXP.csv"),
                    row.names = 1, check.names = FALSE)
gse_mat <- as.matrix(gse_exp)
cat("Running PROGENy mlm on Landa ...\n")
gse_ulm <- suppressMessages(decoupleR::run_mlm(
  mat = gse_mat, net = net_h,
  .source = "source", .target = "target", .mor = "weight",
  minsize = 10))
gse_scores <- gse_ulm %>%
  transmute(sample = condition, pathway = source, score = score)

landa_master <- read_tsv(file.path(out_root,
                                   "factor1_gse76039_projection/factor1_gse76039_scores.tsv"),
                         show_col_types = FALSE) %>%
  dplyr::select(sample, class)
gse_wide <- gse_scores %>%
  pivot_wider(names_from = pathway, values_from = score) %>%
  inner_join(landa_master, by = "sample") %>%
  mutate(class = factor(class, levels = c("PDTC","ATC")))

landa_stats <- tibble(pathway = pw_levels) %>%
  rowwise() %>%
  mutate(pdtc_mean = mean(gse_wide[[pathway]][gse_wide$class == "PDTC"]),
         atc_mean  = mean(gse_wide[[pathway]][gse_wide$class == "ATC"]),
         delta_ATC_minus_PDTC = atc_mean - pdtc_mean,
         wilcox_p = tryCatch(
           wilcox.test(gse_wide[[pathway]] ~ gse_wide$class)$p.value,
           error = function(e) NA_real_)) %>%
  ungroup() %>%
  mutate(pathway = factor(pathway, levels = pw_levels))
write_tsv(landa_stats, file.path(fig_dir, "landa_progeny_stats.tsv"))

## ============================================================
## 5. Cross-cohort consistency table
##    We compare direction of pathway changes:
##    Mouse: FAP-KO minus WT
##    TCGA:  BRAF minus RAS
##    Landa: ATC minus PDTC
## ============================================================
cross_tbl <- tibble(pathway = pw_levels) %>%
  left_join(mouse_stats %>% transmute(pathway, mouse_delta = delta_fap_minus_wt,
                                       mouse_p = wilcox_p),
            by = "pathway") %>%
  left_join(tcga_stats %>% transmute(pathway, tcga_delta = delta_BRAF_minus_RAS,
                                      tcga_p = wilcox_p),
            by = "pathway") %>%
  left_join(landa_stats %>% transmute(pathway, landa_delta = delta_ATC_minus_PDTC,
                                       landa_p = wilcox_p),
            by = "pathway") %>%
  mutate(pathway = factor(pathway, levels = pw_levels))
write_tsv(cross_tbl, file.path(fig_dir, "cross_cohort_progeny_table.tsv"))

cat("\n==== Mouse FAP-KO vs WT (top 5 significant) ====\n")
print(mouse_stats %>% arrange(wilcox_p) %>% head(5) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))
cat("\n==== TCGA BRAF vs RAS (top 5 significant) ====\n")
print(tcga_stats %>% arrange(wilcox_p) %>% head(5) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))
cat("\n==== Landa ATC vs PDTC (top 5 significant) ====\n")
print(landa_stats %>% arrange(wilcox_p) %>% head(5) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))
cat("\n==== Factor1 pathway correlations (TCGA, top 5) ====\n")
print(factor1_cor %>% arrange(p) %>% head(5) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))

## ============================================================
## Panels
## ============================================================
## Panel a — Mouse Fig 8 delta (FAP-KO - WT) bar
a_df <- mouse_stats %>%
  mutate(pathway = factor(pathway, levels = rev(pw_levels)),
         tag = case_when(wilcox_p < 0.05 ~ "*",
                         wilcox_p < 0.1  ~ ".",
                         TRUE ~ ""))
max_a <- max(abs(a_df$delta_fap_minus_wt)) * 1.1
p_A <- ggplot(a_df, aes(delta_fap_minus_wt, pathway, fill = delta_fap_minus_wt)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_col(colour = "black", linewidth = 0.2) +
  geom_text(aes(x = delta_fap_minus_wt + sign(delta_fap_minus_wt) * max_a * 0.05,
                label = tag),
            size = (base_font_size + 2.5)/.pt, fontface = "bold") +
  scale_fill_gradient2(low = col_dn, mid = "#F5F5F5", high = col_up,
                       midpoint = 0, limits = c(-max_a, max_a),
                       oob = squish, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_discrete(labels = pw_lab) +
  labs(x = "Delta activity (FAP-KO − WT)", y = NULL,
       title = "a  Mouse: PROGENy deltas (Fig 8)",
       subtitle = "* P<0.05, . <0.1 (Wilcoxon)",
       tag = "a") +
  theme_nc() +
  theme(axis.text.y = element_text(size = base_font_size - 0.5))

## Panel b — TCGA BRAF-RAS delta bar
b_df <- tcga_stats %>%
  mutate(pathway = factor(pathway, levels = rev(pw_levels)),
         tag = case_when(wilcox_p < 0.05 ~ "*",
                         wilcox_p < 0.1  ~ ".",
                         TRUE ~ ""))
max_b <- max(abs(b_df$delta_BRAF_minus_RAS)) * 1.1
p_B <- ggplot(b_df, aes(delta_BRAF_minus_RAS, pathway, fill = delta_BRAF_minus_RAS)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_col(colour = "black", linewidth = 0.2) +
  geom_text(aes(x = delta_BRAF_minus_RAS + sign(delta_BRAF_minus_RAS) * max_b * 0.05,
                label = tag),
            size = (base_font_size + 2.5)/.pt, fontface = "bold") +
  scale_fill_gradient2(low = col_dn, mid = "#F5F5F5", high = col_up,
                       midpoint = 0, limits = c(-max_b, max_b),
                       oob = squish, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = "Delta activity (BRAF − RAS)", y = NULL,
       title = "b  TCGA-THCA: PROGENy deltas",
       subtitle = "n=213 BRAF vs 48 RAS",
       tag = "b") +
  theme_nc() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

## Panel c — Landa ATC-PDTC delta bar
c_df <- landa_stats %>%
  mutate(pathway = factor(pathway, levels = rev(pw_levels)),
         tag = case_when(wilcox_p < 0.05 ~ "*",
                         wilcox_p < 0.1  ~ ".",
                         TRUE ~ ""))
max_c <- max(abs(c_df$delta_ATC_minus_PDTC)) * 1.1
p_C <- ggplot(c_df, aes(delta_ATC_minus_PDTC, pathway, fill = delta_ATC_minus_PDTC)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_col(colour = "black", linewidth = 0.2) +
  geom_text(aes(x = delta_ATC_minus_PDTC + sign(delta_ATC_minus_PDTC) * max_c * 0.05,
                label = tag),
            size = (base_font_size + 2.5)/.pt, fontface = "bold") +
  scale_fill_gradient2(low = col_dn, mid = "#F5F5F5", high = col_up,
                       midpoint = 0, limits = c(-max_c, max_c),
                       oob = squish, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = "Delta activity (ATC − PDTC)", y = NULL,
       title = "c  Landa GSE76039: PROGENy deltas",
       subtitle = "n=20 ATC vs 17 PDTC",
       tag = "c") +
  theme_nc() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

## Panel d — Cross-cohort consistency heatmap
hm_df <- cross_tbl %>%
  pivot_longer(c(mouse_delta, tcga_delta, landa_delta),
               names_to = "cohort", values_to = "delta") %>%
  pivot_longer(c(mouse_p, tcga_p, landa_p),
               names_to = "cohort_p", values_to = "p") %>%
  filter(sub("_delta","", cohort) == sub("_p","", cohort_p)) %>%
  mutate(cohort = recode(cohort,
                          mouse_delta = "Mouse\n(FAP-KO - WT)",
                          tcga_delta  = "TCGA\n(BRAF - RAS)",
                          landa_delta = "Landa\n(ATC - PDTC)"),
         cohort = factor(cohort,
                         levels = c("Mouse\n(FAP-KO - WT)",
                                    "TCGA\n(BRAF - RAS)",
                                    "Landa\n(ATC - PDTC)")),
         pathway = factor(pathway, levels = rev(pw_levels)),
         tag = case_when(is.na(p) ~ "",
                         p < 0.05 ~ "*",
                         p < 0.1  ~ ".",
                         TRUE ~ ""),
         # normalize each cohort's delta range to [-1,1] for comparable heat
         max_c = ave(abs(delta), cohort, FUN = function(x) max(x, na.rm = TRUE)),
         delta_norm = delta / max_c)

p_D <- ggplot(hm_df, aes(cohort, pathway, fill = delta_norm)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%s%.2f", tag, delta),
                colour = abs(delta_norm) > 0.55),
            size = (base_font_size - 2)/.pt, show.legend = FALSE) +
  scale_fill_gradient2(low = col_dn, mid = "#F5F5F5", high = col_up,
                       midpoint = 0, limits = c(-1, 1), oob = squish,
                       name = "Normalized\ndelta") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(labels = pw_lab) +
  labs(x = NULL, y = NULL,
       title = "d  Cross-cohort PROGENy consistency",
       subtitle = "raw delta annotated; within-cohort [-1,1] normalized for color",
       tag = "d") +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = base_font_size - 0.5),
        axis.text.x.top = element_text(face = "bold",
                                        size = base_font_size - 0.5,
                                        lineheight = 0.9),
        legend.key.height = unit(6, "mm"),
        legend.key.width  = unit(2.6, "mm"))

## Panel e — Factor1 × pathway correlation in TCGA
e_df <- factor1_cor %>%
  mutate(pathway = factor(pathway, levels = rev(pw_levels)),
         tag = case_when(p < 0.001 ~ "***",
                         p < 0.01  ~ "**",
                         p < 0.05  ~ "*",
                         TRUE ~ ""))
p_E <- ggplot(e_df, aes(rho, pathway, fill = rho)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_col(colour = "black", linewidth = 0.2) +
  geom_text(aes(x = rho + sign(rho) * 0.03, label = tag),
            size = (base_font_size + 2.5)/.pt, fontface = "bold") +
  scale_fill_gradient2(low = col_dn, mid = "#F5F5F5", high = col_up,
                       midpoint = 0, limits = c(-1, 1), oob = squish,
                       guide = "none") +
  scale_x_continuous(limits = c(-0.1, 0.82), breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
  scale_y_discrete(labels = pw_lab) +
  labs(x = "Spearman rho (TCGA Factor1 × pathway)", y = NULL,
       title = "e  TCGA Factor1 x PROGENy pathway activity",
       subtitle = "* P<0.05, ** <0.01, *** <0.001",
       tag = "e") +
  theme_nc() +
  theme(axis.text.y = element_text(size = base_font_size - 0.5))

## Panel f — 4 focus pathways on 5-tier ladder (TCGA driver + Landa)
focus_pw <- c("MAPK","WNT","TGFb","JAK-STAT")

ladder_tcga <- tcga_wide %>%
  mutate(ladder = case_when(
    brs_group == "RAS driver" ~ "TCGA RAS",
    brs_group == "BRAF V600E" ~ "TCGA BRAF",
    TRUE ~ "TCGA other")) %>%
  dplyr::select(sample, ladder, all_of(focus_pw))
ladder_landa <- gse_wide %>%
  mutate(ladder = paste0("Landa ", as.character(class))) %>%
  dplyr::select(sample, ladder, all_of(focus_pw))

## within-cohort z-score normalization before combining (batch fix)
znorm <- function(df, cols){
  for (c in cols) if (c %in% colnames(df)) df[[c]] <- as.numeric(scale(df[[c]]))
  df
}
ladder_tcga  <- znorm(ladder_tcga,  focus_pw)
ladder_landa <- znorm(ladder_landa, focus_pw)

ladder_all <- bind_rows(ladder_tcga, ladder_landa) %>%
  mutate(ladder = factor(ladder,
                         levels = c("TCGA RAS", "TCGA other",
                                    "TCGA BRAF", "Landa PDTC", "Landa ATC")))

lad_cols <- c(
  "TCGA RAS"   = "#3B6B9C",
  "TCGA other" = "#B7B7B7",
  "TCGA BRAF"  = "#5E4A7B",
  "Landa PDTC" = "#E08F3A",
  "Landa ATC"  = "#C9493A"
)

ladder_long <- ladder_all %>%
  pivot_longer(-c(sample, ladder), names_to = "pathway", values_to = "z") %>%
  mutate(pathway = factor(pathway, levels = focus_pw))

kt_stats <- ladder_long %>%
  group_by(pathway) %>%
  summarise(tau = suppressWarnings(cor.test(as.numeric(ladder), z,
                                             method = "kendall")$estimate),
            p = suppressWarnings(cor.test(as.numeric(ladder), z,
                                           method = "kendall")$p.value),
            .groups = "drop") %>%
  mutate(pathway = factor(pathway, levels = focus_pw),
         label = sprintf("tau = %.2f\n%s", tau, vapply(p, p_label, character(1))))

p_F <- ggplot(ladder_long, aes(ladder, z, fill = ladder)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.10, height = 0, size = 0.22,
              shape = 16, colour = "grey30", alpha = 0.3) +
  geom_text(data = kt_stats, aes(x = 3, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.4,
            size = (base_font_size - 1.5)/.pt, lineheight = 1.0) +
  facet_wrap(~ pathway, nrow = 1, scales = "free_y", labeller = as_labeller(pw_lab)) +
  scale_fill_manual(values = lad_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "PROGENy activity (within-cohort z)",
       title = "f  Dedifferentiation-grade ladder - focus pathways",
       subtitle = "MAPK / WNT / TGFb / JAK-STAT",
       tag = "f") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1,
                                   size = base_font_size - 1),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        panel.spacing.x = unit(2, "mm"))

## ============================================================
## assemble
## ============================================================
design <- "
AAABBBCCC
AAABBBCCC
AAABBBCCC
DDDDEEEEE
DDDDEEEEE
DDDDEEEEE
FFFFFFFFF
FFFFFFFFF
"

fig_s17 <- p_A + p_B + p_C + p_D + p_E + p_F +
  plot_layout(design = design) +
  plot_annotation(
    title = "Fig. S17. PROGENy 14-pathway activity across mouse FAP-KO, TCGA-THCA, and Landa PDTC/ATC.",
    subtitle = "Orthogonal peer-reviewed pathway framework; confirms MAPK / JAK-STAT / TNFa directionality with Wnt/TGFb asymmetry.",
    theme = theme(plot.title = element_text(face = "bold",
                                             size = base_font_size + 1),
                  plot.subtitle = element_text(size = base_font_size - 0.5,
                                                colour = "grey30"))) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(fig_dir, "FigS17_progeny_cross_cohort.pdf")
png_path <- file.path(fig_dir, "FigS17_progeny_cross_cohort.png")

ggsave(pdf_path, fig_s17, width = 240, height = 290, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(png_path, fig_s17, width = 240, height = 290, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("\nFig S17 written:\n  ", pdf_path, "\n  ", png_path, "\n")
