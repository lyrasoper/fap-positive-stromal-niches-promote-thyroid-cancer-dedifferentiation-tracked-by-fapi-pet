#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Figure 4E and Supplementary Fig. S17 - ECM/myCAF core-gene selection.
## B1: Core-gene analysis of the ecm-myCAF / LRRC15+ myCAF subtype.
##
## Intersection of three independent published signatures:
##   Elyada 2019 myCAF (PDAC scRNA-seq)        — 16 genes
##   Dominguez 2020 LRRC15+ myCAF (pan-cancer)  — 18 genes
##   Kieffer 2020 ecm-myCAF (breast)           — 13 genes
##
## Goal: show that the 4-12 shared scaffold genes (FAP, LRRC15, POSTN,
## MMP11, COL1A1/2, COL3A1, MFAP5, COMP, TNC, THY1, ACTA2, TAGLN)
## concordantly collapse in FAP-KO at both RNA and protein levels,
## converting the bulk-signature trend into single-gene evidence.


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
})

proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
bulk_out <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
out_dir  <- file.path(bulk_out, "ecm_mycaf_core_genes")
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
to_mgi <- function(x){
  x <- toupper(x); paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
}

## ---- define the 3 myCAF signatures -------------------------------------
elyada <- c("Acta2","Tagln","Myl9","Tpm1","Tpm2","Mmp11","Postn","Hopx",
            "Thy1","Gja4","Ctgf","Thbs1","Col1a1","Col3a1","Lrrc15","Fap")
dominguez <- c("Lrrc15","Col10a1","Col11a1","Cthrc1","Inhba","Mmp11",
               "Mfap5","Postn","Sulf1","Comp","Fap","Ltbp2","Tagln",
               "Acta2","Tnc","Thbs2","Cilp","Fbn1")
kieffer <- c("Fap","Lrrc15","Col1a1","Col1a2","Col3a1","Mfap5","Postn",
             "Mmp11","Comp","Tnc","Pdpn","Thy1","S100a4")

sig_list <- list(Elyada_myCAF = elyada,
                 Dominguez_LRRC15_myCAF = dominguez,
                 Kieffer_ecm_myCAF = kieffer)

all_genes <- unique(unlist(sig_list))
membership <- sapply(sig_list, function(s) all_genes %in% s)
rownames(membership) <- all_genes
overlap_n <- rowSums(membership)

core3 <- all_genes[overlap_n == 3]       # in all 3
core2 <- all_genes[overlap_n >= 2]       # in >= 2
cat(sprintf("3-way intersection: %d genes: %s\n",
            length(core3), paste(core3, collapse = ", ")))
cat(sprintf("2+ signatures:     %d genes: %s\n",
            length(core2), paste(core2, collapse = ", ")))

overlap_tbl <- tibble(gene = all_genes,
                      n_sig = overlap_n,
                      Elyada    = membership[, "Elyada_myCAF"],
                      Dominguez = membership[, "Dominguez_LRRC15_myCAF"],
                      Kieffer   = membership[, "Kieffer_ecm_myCAF"]) %>%
  arrange(desc(n_sig), gene)
write_tsv(overlap_tbl, file.path(out_dir, "signature_overlap_table.tsv"))

## ---- load DE and expression --------------------------------------------
deg_rna  <- read_tsv(file.path(bulk_out, "rna/rna_differential_expression.tsv"),
                     show_col_types = FALSE) %>%
  mutate(gene_mgi = to_mgi(gene))
deg_prot <- read_tsv(file.path(bulk_out, "proteomics/proteomics_differential_expression.tsv"),
                     show_col_types = FALSE) %>%
  mutate(gene_mgi = to_mgi(gene))

fisher_combine <- function(p1, p2){
  p1 <- pmax(p1, 1e-300); p2 <- pmax(p2, 1e-300)
  pchisq(-2 * (log(p1) + log(p2)), df = 4, lower.tail = FALSE)
}

core_tbl <- tibble(gene_mgi = to_mgi(all_genes)) %>%
  left_join(deg_rna %>% transmute(gene_mgi,
                                   rna_log2fc  = log2_fc,
                                   rna_p       = p_value,
                                   rna_q       = q_value,
                                   rna_sig     = significant == "yes"),
            by = "gene_mgi") %>%
  left_join(deg_prot %>% transmute(gene_mgi,
                                    prot_log2fc = log2_fc,
                                    prot_p      = p_value,
                                    prot_q      = q_value,
                                    prot_sig    = significant == "yes"),
            by = "gene_mgi") %>%
  left_join(overlap_tbl %>% transmute(gene_mgi = to_mgi(gene),
                                       n_sig, Elyada, Dominguez, Kieffer),
            by = "gene_mgi") %>%
  mutate(sign_concordant = sign(rna_log2fc) == sign(prot_log2fc),
         fisher_p = fisher_combine(rna_p, prot_p))
core_tbl <- core_tbl %>% mutate(fisher_q = p.adjust(fisher_p, method = "BH"))

write_tsv(core_tbl, file.path(out_dir, "core_gene_stats.tsv"))

cat("\nPer-gene table (2+ signatures):\n")
print(core_tbl %>% filter(n_sig >= 2) %>%
        dplyr::select(gene_mgi, n_sig, rna_log2fc, rna_q, prot_log2fc, prot_p,
                      fisher_p, fisher_q, sign_concordant) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))

## ---- sample-level FPKM and protein loading ------------------------------
fpkm_tab <- read_tsv(file.path(proj, "001_mRNA_Summary/1.GeneExpression/1_genes_fpkm_expression.txt"),
                     show_col_types = FALSE, progress = FALSE)
prot_tab <- read_tsv(file.path(proj, "002_DIA_Summary/02.ProteinExp/protein_annotation_profile.txt"),
                     show_col_types = FALSE, progress = FALSE)

samp_order <- c("WT_BPC1","WT_BPC2","WT_BPC3","WT_BPC4",
                "FAP_BPC1","FAP_BPC2","FAP_BPC3","FAP_BPC4")

get_rna_mat <- function(genes){
  rows <- fpkm_tab %>% filter(gene_name %in% genes) %>%
    distinct(gene_name, .keep_all = TRUE)
  cols <- paste0("FPKM.", samp_order)
  mat <- as.matrix(rows[, cols]); rownames(mat) <- rows$gene_name
  colnames(mat) <- samp_order
  log2(mat + 1)
}
get_prot_mat <- function(genes){
  rows <- prot_tab %>% mutate(g = to_mgi(gene_name)) %>%
    filter(g %in% to_mgi(genes)) %>% distinct(g, .keep_all = TRUE)
  cols <- samp_order[samp_order %in% colnames(rows)]
  mat <- suppressWarnings(as.matrix(rows[, cols, drop = FALSE]))
  rownames(mat) <- rows$g
  mode(mat) <- "numeric"
  # impute NA with column min/2
  for (j in seq_len(ncol(mat))){
    mn <- suppressWarnings(min(mat[, j], na.rm = TRUE))
    if (is.infinite(mn)) mn <- 0
    mat[is.na(mat[, j]), j] <- mn / 2
  }
  log2(mat + 1)
}

rna_core_mat  <- get_rna_mat(core2)
prot_core_mat <- get_prot_mat(core2)

## ---- Panel A — Venn-style overlap bar ----------------------------------
panel_a_df <- overlap_tbl %>%
  mutate(gene = factor(gene, levels = rev(overlap_tbl$gene)),
         in_3 = n_sig == 3,
         group = case_when(
           n_sig == 3 ~ "3-way (scaffold core)",
           n_sig == 2 ~ "2-way (extended core)",
           TRUE ~ "unique")) %>%
  pivot_longer(c(Elyada, Dominguez, Kieffer),
               names_to = "signature", values_to = "present") %>%
  mutate(signature = factor(signature,
                            levels = c("Elyada", "Dominguez", "Kieffer")))

p_a <- ggplot(panel_a_df, aes(signature, gene, fill = present)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = subset(panel_a_df, present),
             aes(colour = group), size = 1.1, shape = 16) +
  scale_fill_manual(values = c(`TRUE` = "#F0E4D7", `FALSE` = "white"),
                    guide = "none") +
  scale_colour_manual(values = c("3-way (scaffold core)" = col_fap,
                                  "2-way (extended core)" = col_purple,
                                  "unique" = col_grey),
                     name = NULL) +
  scale_x_discrete(position = "top",
                   labels = c("Elyada\n2019","Dominguez\n2020","Kieffer\n2020")) +
  labs(x = NULL, y = NULL,
       title = "Signature membership",
       subtitle = sprintf("%d / %d genes shared by 3 / 2 of 3 signatures",
                          sum(overlap_tbl$n_sig == 3),
                          sum(overlap_tbl$n_sig >= 2)),
       tag = "a") +
  theme_nc() +
  theme(axis.text.y = element_text(face = "italic",
                                   size = base_font_size - 1),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(2.8, "mm"))

## ---- Panel B — Dual-modality log2FC forest for 2+ signatures genes ----
forest_df <- core_tbl %>% filter(n_sig >= 2) %>%
  arrange(n_sig, rna_log2fc) %>%
  mutate(gene_mgi = factor(gene_mgi, levels = rev(gene_mgi))) %>%
  dplyr::select(gene_mgi, n_sig,
                RNA = rna_log2fc, Proteomics = prot_log2fc,
                rna_q, prot_p) %>%
  pivot_longer(c(RNA, Proteomics), names_to = "modality", values_to = "log2fc") %>%
  mutate(modality = factor(modality, levels = c("RNA","Proteomics")))

forest_sig <- core_tbl %>% filter(n_sig >= 2) %>%
  mutate(gene_mgi = factor(gene_mgi,
                           levels = rev((core_tbl %>%
                                          filter(n_sig >= 2) %>%
                                          arrange(n_sig, rna_log2fc))$gene_mgi)),
         tag_rna = case_when(rna_q < 0.05 ~ "*", rna_q < 0.1 ~ ".", TRUE ~ ""),
         tag_prot = case_when(prot_p < 0.05 ~ "*", prot_p < 0.1 ~ ".", TRUE ~ ""))

max_abs <- max(abs(forest_df$log2fc), na.rm = TRUE)

p_b <- ggplot(forest_df, aes(log2fc, gene_mgi, fill = log2fc)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
  geom_col(colour = "black", linewidth = 0.2) +
  geom_text(data = forest_sig,
            aes(x = rna_log2fc + sign(rna_log2fc) * 0.15,
                y = gene_mgi, label = tag_rna),
            inherit.aes = FALSE,
            size = (base_font_size + 0.5)/.pt) +
  facet_wrap(~ modality, nrow = 1) +
  scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                       midpoint = 0, guide = "none",
                       limits = c(-max_abs, max_abs), oob = squish) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(x = "log2 FC (FAP-KO vs WT host)",
       y = NULL,
       title = "Scaffold genes shared by 2+ signatures",
       subtitle = "* FDR<0.05 (RNA) or P<0.05 (Prot), . <0.1",
       tag = "b") +
  theme_nc() +
  theme(axis.text.y = element_text(face = "italic",
                                   size = base_font_size - 0.5),
        panel.spacing.x = unit(2, "mm"))

## ---- Panel C — RNA vs Protein log2FC scatter (sign concordance) -------
sc_df <- core_tbl %>% filter(n_sig >= 2,
                               !is.na(rna_log2fc), !is.na(prot_log2fc)) %>%
  mutate(class = case_when(n_sig == 3 ~ "3-way core",
                           n_sig == 2 ~ "2-way extended"),
         label = gene_mgi)

# simple binomial test for concordance (one-sided, H0: p=0.5)
bn <- binom.test(sum(sc_df$sign_concordant), nrow(sc_df),
                 p = 0.5, alternative = "greater")
rho_sc <- suppressWarnings(cor(sc_df$rna_log2fc, sc_df$prot_log2fc,
                               method = "spearman"))

max_abs_sc <- max(abs(c(sc_df$rna_log2fc, sc_df$prot_log2fc)), na.rm = TRUE) * 1.05

p_c <- ggplot(sc_df, aes(rna_log2fc, prot_log2fc)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              linewidth = 0.25, colour = "grey70") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey30", fill = "grey85", linewidth = 0.3) +
  geom_point(aes(fill = class), shape = 21, size = 2.2,
             stroke = 0.25, colour = "black") +
  geom_text_repel(aes(label = label, colour = class),
                  size = (base_font_size - 1.5)/.pt,
                  min.segment.length = 0, segment.size = 0.2,
                  max.overlaps = Inf, box.padding = 0.3,
                  fontface = "italic", show.legend = FALSE) +
  scale_fill_manual(values = c("3-way core" = col_fap,
                               "2-way extended" = col_purple),
                    name = NULL) +
  scale_colour_manual(values = c("3-way core" = col_fap,
                                  "2-way extended" = col_purple)) +
  coord_cartesian(xlim = c(-max_abs_sc, max_abs_sc),
                  ylim = c(-max_abs_sc, max_abs_sc)) +
  labs(x = "RNA log2FC (FAP-KO vs WT)",
       y = "Protein log2FC (FAP-KO vs WT)",
       title = "Scaffold genes collapse on both modalities",
       subtitle = sprintf("Spearman %.2f | %d / %d sign-concordant (binomial one-sided P = %.2g)",
                          rho_sc, sum(sc_df$sign_concordant), nrow(sc_df),
                          bn$p.value),
       tag = "c") +
  theme_nc() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"))

## ---- Panel D — Per-sample heatmap for 3-way core genes ----------------
hm_genes <- core2
rna_z <- t(scale(t(rna_core_mat[intersect(rownames(rna_core_mat), hm_genes), ])))
prot_z <- t(scale(t(prot_core_mat[intersect(rownames(prot_core_mat), to_mgi(hm_genes)), ])))

row_order <- core_tbl %>% filter(gene_mgi %in% hm_genes) %>%
  arrange(n_sig, rna_log2fc) %>% pull(gene_mgi)

z_df <- bind_rows(
  as.data.frame(rna_z) %>% rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "z") %>%
    mutate(modality = "RNA"),
  as.data.frame(prot_z) %>% rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "z") %>%
    mutate(modality = "Protein")
) %>%
  mutate(gene = factor(gene, levels = rev(row_order)),
         modality = factor(modality, levels = c("RNA","Protein")),
         sample = factor(sample, levels = samp_order),
         group = factor(ifelse(grepl("^WT", sample), "WT host",
                               "FAP-deficient host"),
                        levels = grp_levels))

p_d <- ggplot(z_df, aes(sample, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  facet_grid(. ~ modality, scales = "free_x") +
  scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                       midpoint = 0, limits = c(-2.2, 2.2), oob = squish,
                       name = "z-score") +
  labs(x = NULL, y = NULL,
       title = "Per-sample z-score of scaffold genes",
       subtitle = "rows sorted by membership + RNA log2FC",
       tag = "d") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1,
                                   size = base_font_size - 1.5,
                                   colour = ifelse(samp_order %in%
                                                     samp_order[grepl("^WT", samp_order)],
                                                   col_wt, col_fap)),
        axis.text.y = element_text(face = "italic",
                                   size = base_font_size - 0.5),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key.height = unit(6, "mm"),
        legend.key.width  = unit(2.4, "mm"),
        panel.spacing.x = unit(2, "mm"),
        strip.text = element_text(face = "bold"))

## ---- Panel E — Fisher combined P ranking bar --------------------------
rank_df <- core_tbl %>% filter(n_sig >= 2, !is.na(fisher_p)) %>%
  arrange(fisher_p) %>%
  mutate(rank_order = row_number(),
         neg_log10p = -log10(pmax(fisher_p, 1e-300)),
         gene_mgi = factor(gene_mgi, levels = rev(gene_mgi)),
         sig_tag = case_when(fisher_q < 0.05 ~ "*",
                             fisher_q < 0.25 ~ ".",
                             TRUE ~ ""))

p_e <- ggplot(rank_df, aes(neg_log10p, gene_mgi, fill = rna_log2fc)) +
  geom_col(colour = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%sq = %.2g",
                                 sig_tag, fisher_q)),
            hjust = -0.1, size = (base_font_size - 1.5)/.pt) +
  scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                       midpoint = 0, oob = squish,
                       limits = c(-4, 4),
                       name = "RNA\nlog2FC") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.35))) +
  labs(x = expression(-log[10]~"Fisher combined P"),
       y = NULL,
       title = "Per-gene cross-modality evidence",
       subtitle = "BH-adj q shown; * q<0.05",
       tag = "e") +
  theme_nc() +
  theme(axis.text.y = element_text(face = "italic",
                                   size = base_font_size - 1),
        legend.key.height = unit(6, "mm"),
        legend.key.width  = unit(2.4, "mm"))

## ---- assemble ---------------------------------------------------------
design <- "
AAABBBCCCC
AAABBBCCCC
AAABBBCCCC
DDDDDDEEEE
DDDDDDEEEE
DDDDDDEEEE
"
fig <- p_a + p_b + p_c + p_d + p_e +
  plot_layout(design = design) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(out_dir, "ecm_mycaf_core_gene_overview.pdf")
png_path <- file.path(out_dir, "ecm_mycaf_core_gene_overview.png")
ggsave(pdf_path, fig, width = 240, height = 220, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(png_path, fig, width = 240, height = 220, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("\nCore-gene outputs written to:\n  ", out_dir, "\n")

## ---- short text summary -----------------------------------------------
cat("\n==== TAKEAWAY ====\n")
cat(sprintf("Total unique genes across 3 signatures: %d\n", length(all_genes)))
cat(sprintf("3-way intersection (scaffold core): %d genes: %s\n",
            length(core3), paste(core3, collapse = ", ")))
cat(sprintf("2+ signatures: %d genes\n", length(core2)))
cat(sprintf("Cross-modality sign-concordant (among 2+ core): %d / %d  (binomial one-sided P = %.2g)\n",
            sum(sc_df$sign_concordant), nrow(sc_df), bn$p.value))
cat(sprintf("Spearman (RNA vs Prot log2FC): %.2f\n", rho_sc))
cat("\nTop-ranked cross-modality scaffold genes:\n")
print(rank_df %>% slice_head(n = 8) %>%
        dplyr::select(gene_mgi, n_sig, rna_log2fc, prot_log2fc,
                      rna_q, prot_p, fisher_p, fisher_q) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))
