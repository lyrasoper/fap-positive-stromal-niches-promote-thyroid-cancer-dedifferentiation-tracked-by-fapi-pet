#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S14 (human spatial candidates vs tumors from
#          Fap-deficient hosts; decoupleR ULM scoring).
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

## REBUILD Supp14_crossecho (src build_fig7_fig8_cross_echo_v2_itga5_highlight.R) beautify 2026-06-20: quartz vector PDF + 180mm (was 240x280). Data identical. -> rebuilds/
## Fig. S12 — Fig 7 <-> Fig 8 cross-echo.
##
## Goal: explicitly bridge the HUMAN spatial transcriptomics (Fig 7)
## and the MOUSE multi-omics xenograft (Fig 8) evidence by:
##
##   A. Projecting Fig 7's top ECM-integrin sender-receiver L-R pairs
##      (COL1A1/2/3, POSTN, THBS2, VCAN, COL6A3 <-> ITGA2/ITGB1/ITGA5)
##      onto Fig 8 bulk RNA + DIA proteomics FAP-KO vs WT log2FC, and
##   B. Scoring Fig 7's 4 epithelial state signatures
##      (follicular_lineage_high, lineage_preserved_epithelial,
##       plasticity_high_transition, terminal_dedifferentiated)
##      via decoupleR::run_ulm on Fig 8 samples, showing directionality
##      WT vs FAP-KO.
##
## Source files (Fig 7):
##   ranked_candidate_interactions.tsv          — L-R pairs
##   epithelial_state_score_gene_sets.tsv       — state signatures

.libPaths(c(file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_omics/.Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(decoupleR)
})

proj <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "bulk_omics")
out_root <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
fig_dir <- file.path(out_root, "fig7_fig8_cross_echo")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

fig7_root <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/results")

col_wt   <- "#3B6B9C"; col_fap  <- "#C9493A"; col_grey <- "#B7B7B7"
col_up   <- col_fap;   col_dn   <- col_wt
col_purple <- "#4B1E70"
grp_levels <- c("WT host", "FAP-deficient host")
grp_cols <- setNames(c(col_wt, col_fap), grp_levels)

base_font_size <- 7
theme_nc <- function(base = base_font_size){
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(size = 7.5, hjust = 0, face = "bold",
                                   margin = margin(b = 1.5)),
      plot.subtitle = element_text(size = 6.5, colour = "grey30",
                                   margin = margin(b = 2)),
      plot.tag      = element_text(face = "bold", size = base + 3, hjust = 0),
      axis.title    = element_text(size = base),
      axis.text     = element_text(size = base - 0.5, colour = "black"),
      axis.line     = element_line(linewidth = 0.6),
      axis.ticks    = element_line(linewidth = 0.6),
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
to_mgi <- function(x){
  x <- toupper(x); paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
}

## Guard for upstream intermediate files (produced by earlier Fig 7/Fig 8 steps).
read_tsv_guard <- function(p, ...){
  if (!file.exists(p))
    stop(sprintf("Missing %s — produced by an upstream step; see README.", p))
  read_tsv(p, ...)
}

## ============================================================
## Load Fig 7 assets
## ============================================================
lr_all <- read_tsv_guard(
  file.path(fig7_root, "20260410_fap_fibro_epi_mechanism_prioritization",
            "ranked_candidate_interactions.tsv"),
  show_col_types = FALSE)

state_sets <- read_tsv_guard(
  file.path(fig7_root, "20260409_fap_epi_state_association",
            "epithelial_state_score_gene_sets.tsv"),
  show_col_types = FALSE)

## ============================================================
## Load Fig 8 data
## ============================================================
deg_rna  <- read_tsv_guard(file.path(out_root, "rna/rna_differential_expression.tsv"), show_col_types = FALSE)
deg_prot <- read_tsv_guard(file.path(out_root, "proteomics/proteomics_differential_expression.tsv"), show_col_types = FALSE)

fpkm_tab <- read_tsv_guard(file.path(proj, "001_mRNA_Summary/1.GeneExpression/1_genes_fpkm_expression.txt"),
                     show_col_types = FALSE, progress = FALSE)
prot_tab <- read_tsv_guard(file.path(proj, "002_DIA_Summary/02.ProteinExp/protein_annotation_profile.txt"),
                     show_col_types = FALSE, progress = FALSE)

samp_order <- c("WT_BPC1","WT_BPC2","WT_BPC3","WT_BPC4",
                "FAP_BPC1","FAP_BPC2","FAP_BPC3","FAP_BPC4")
samp_grp <- tibble(sample = samp_order,
                   group = factor(ifelse(grepl("^WT", samp_order),
                                         "WT host", "FAP-deficient host"),
                                  levels = grp_levels))

## =======================================================================
## Enhancement A — Fig 7 top L-R pair log2FC in Fig 8
## =======================================================================
## Fig 7 top ECM-integrin + CXCL12 pairs
lr_top <- lr_all %>%
  arrange(priority_rank) %>%
  filter(priority_rank <= 12) %>%
  transmute(priority_rank, ligand, receptor, pathway,
            pair_label = paste0(ligand, " -> ", receptor))

genes_to_test <- unique(c(lr_top$ligand, lr_top$receptor))

## RNA log2FC (mouse MGI symbol)
deg_rna_mgi <- deg_rna %>% mutate(gene_mgi = to_mgi(gene))
deg_prot_mgi <- deg_prot %>% mutate(gene_mgi = to_mgi(gene))

lookup_fc <- function(human_gene, dataset){
  mgi <- to_mgi(human_gene)
  if (dataset == "RNA"){
    hit <- deg_rna_mgi %>% filter(gene_mgi == mgi) %>% slice_head(n = 1)
  } else {
    hit <- deg_prot_mgi %>% filter(gene_mgi == mgi) %>% slice_head(n = 1)
  }
  if (nrow(hit) == 0) return(c(NA, NA))
  c(hit$log2_fc[1], hit$p_value[1])
}

lr_detail <- lr_top %>%
  rowwise() %>%
  mutate(
    ligand_rna_fc   = lookup_fc(ligand,   "RNA")[1],
    ligand_rna_p    = lookup_fc(ligand,   "RNA")[2],
    ligand_prot_fc  = lookup_fc(ligand,   "Protein")[1],
    ligand_prot_p   = lookup_fc(ligand,   "Protein")[2],
    receptor_rna_fc = lookup_fc(receptor, "RNA")[1],
    receptor_rna_p  = lookup_fc(receptor, "RNA")[2],
    receptor_prot_fc= lookup_fc(receptor, "Protein")[1],
    receptor_prot_p = lookup_fc(receptor, "Protein")[2]
  ) %>% ungroup()

write_tsv(lr_detail, file.path(fig_dir, "fig7_lr_pairs_in_fig8_stats.tsv"))

## Panel A1: ligand side (sender)
a_df <- lr_detail %>%
  transmute(pair_label, ligand, receptor, pathway,
            RNA = ligand_rna_fc, Protein = ligand_prot_fc,
            rna_p = ligand_rna_p, prot_p = ligand_prot_p) %>%
  distinct(ligand, .keep_all = TRUE) %>%
  arrange(RNA) %>%
  mutate(gene_ordered = factor(ligand, levels = ligand))

a_long <- a_df %>%
  pivot_longer(c(RNA, Protein), names_to = "modality", values_to = "log2fc") %>%
  mutate(modality = factor(modality, levels = c("RNA", "Protein")))

a_tag <- bind_rows(
  a_df %>% transmute(gene_ordered, log2fc = RNA, modality = "RNA",
                     p = rna_p),
  a_df %>% transmute(gene_ordered, log2fc = Protein, modality = "Protein",
                     p = prot_p)
) %>%
  mutate(modality = factor(modality, levels = c("RNA","Protein")),
         tag = case_when(is.na(p) ~ "NA",
                         p < 0.05 ~ "*",
                         p < 0.1  ~ ".",
                         TRUE ~ ""))

max_abs_a <- max(abs(a_long$log2fc), na.rm = TRUE) * 1.05

p_A1 <- ggplot(a_long, aes(log2fc, gene_ordered, fill = log2fc)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
  geom_col(colour = "black", linewidth = 0.2) +
  geom_text(data = a_tag %>% filter(!is.na(log2fc)),
            aes(x = log2fc + sign(log2fc) * max_abs_a * 0.05,
                y = gene_ordered, label = tag),
            inherit.aes = FALSE,
            size = (base_font_size + 0.5)/.pt) +
  facet_wrap(~ modality, nrow = 1) +
  scale_fill_gradient2(low = col_dn, mid = "#F5F5F5", high = col_up,
                       midpoint = 0, limits = c(-max_abs_a, max_abs_a),
                       oob = squish, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(x = "log2 FC (FAP-KO vs WT)", y = NULL,
       title = "Sender ligands reproduce in mouse FAP-KO",
       subtitle = "Fig 3 ECM/chemokine ligands · mouse RNA | protein",
       tag = "a") +
  theme_nc() +
  theme(axis.text.y = element_text(face = "italic",
                                   size = base_font_size - 0.3),
        panel.spacing.x = unit(2, "mm"))

## Panel A2: receptor side
b_df <- lr_detail %>%
  transmute(pair_label, ligand, receptor, pathway,
            RNA = receptor_rna_fc, Protein = receptor_prot_fc,
            rna_p = receptor_rna_p, prot_p = receptor_prot_p) %>%
  distinct(receptor, .keep_all = TRUE) %>%
  arrange(RNA) %>%
  mutate(gene_ordered = factor(receptor, levels = receptor))

## v2 — markdown-styled receptor labels: ★ ITGA2 (lead α2β1 receiver) in bold deep red
##      (mirrors Fig 4 g; ITGA5 demoted to secondary α5β1 per the α2β1-lead narrative)
recv_levels <- as.character(b_df$gene_ordered)
recv_md_lookup <- setNames(
  ifelse(recv_levels == "ITGA2",
         "<span style='color:#C92E2E;'>**★ ITGA2**</span>",
         paste0("*", recv_levels, "*")),
  recv_levels
)
b_df <- b_df %>%
  mutate(gene_md = factor(recv_md_lookup[as.character(gene_ordered)],
                          levels = unname(recv_md_lookup)))

b_long <- b_df %>%
  pivot_longer(c(RNA, Protein), names_to = "modality", values_to = "log2fc") %>%
  mutate(modality = factor(modality, levels = c("RNA", "Protein")))

b_tag <- bind_rows(
  b_df %>% transmute(gene_md, log2fc = RNA, modality = "RNA",
                     p = rna_p),
  b_df %>% transmute(gene_md, log2fc = Protein, modality = "Protein",
                     p = prot_p)
) %>%
  mutate(modality = factor(modality, levels = c("RNA","Protein")),
         tag = case_when(is.na(p) ~ "NA",
                         p < 0.05 ~ "*",
                         p < 0.1  ~ ".",
                         TRUE ~ ""))

max_abs_b <- max(abs(b_long$log2fc), na.rm = TRUE) * 1.05

p_A2 <- ggplot(b_long, aes(log2fc, gene_md, fill = log2fc)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
  geom_col(colour = "black", linewidth = 0.2) +
  geom_text(data = b_tag %>% filter(!is.na(log2fc)),
            aes(x = log2fc + sign(log2fc) * max_abs_b * 0.05,
                y = gene_md, label = tag),
            inherit.aes = FALSE,
            size = (base_font_size + 0.5)/.pt) +
  facet_wrap(~ modality, nrow = 1) +
  scale_fill_gradient2(low = col_dn, mid = "#F5F5F5", high = col_up,
                       midpoint = 0, limits = c(-max_abs_b, max_abs_b),
                       oob = squish, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(x = "log2 FC (FAP-KO vs WT)", y = NULL,
       title = "Receiver receptors in mouse FAP-KO",
       subtitle = "★ ITGA2 = nominated lead α2β1 receiver",
       tag = "b") +
  theme_nc() +
  theme(axis.text.y = ggtext::element_markdown(size = base_font_size - 0.3),
        panel.spacing.x = unit(2, "mm"))

## =======================================================================
## Enhancement B — Fig 7 state signatures scored on Fig 8 samples (ULM)
## =======================================================================
state_df <- state_sets %>%
  select(state, genes_present) %>%
  separate_rows(genes_present, sep = ";") %>%
  transmute(state = gsub("_", " ", state),
            gene_human = genes_present,
            gene_mouse = to_mgi(genes_present))

write_tsv(state_df, file.path(fig_dir, "fig7_states_mouse_symbols.tsv"))

## Build RNA log2 FPKM matrix in mouse MGI symbols
rna_cols <- paste0("FPKM.", samp_order)
rna_mat <- fpkm_tab %>%
  select(gene_name, all_of(rna_cols)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>%
  slice_max(rowSums(across(all_of(rna_cols))), n = 1, with_ties = FALSE) %>%
  ungroup() %>% column_to_rownames("gene_name") %>% as.matrix()
colnames(rna_mat) <- samp_order
rna_log <- log2(rna_mat + 1)

## Protein matrix
prot_mat <- prot_tab %>%
  select(gene_name, all_of(samp_order)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  mutate(across(all_of(samp_order), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(gene_mgi = to_mgi(gene_name)) %>%
  group_by(gene_mgi) %>%
  summarise(across(all_of(samp_order), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  column_to_rownames("gene_mgi") %>% as.matrix()
prot_log <- log2(prot_mat + 1)
prot_log[!is.finite(prot_log)] <- NA
col_min_half <- apply(prot_log, 2, function(x) min(x, na.rm = TRUE) - 1)
for (j in seq_len(ncol(prot_log))) prot_log[is.na(prot_log[, j]), j] <- col_min_half[j]

net_state <- state_df %>% transmute(source = state,
                                     target = gene_mouse,
                                     mor = 1)

ulm_rna <- suppressMessages(decoupleR::run_ulm(
  mat = rna_log, net = net_state,
  .source = "source", .target = "target", .mor = "mor",
  minsize = 3))
ulm_prot <- suppressMessages(decoupleR::run_ulm(
  mat = prot_log, net = net_state,
  .source = "source", .target = "target", .mor = "mor",
  minsize = 3))

score_rna  <- ulm_rna  %>% transmute(state = source, sample = condition,
                                      score = score) %>%
  mutate(dataset = "RNA")
score_prot <- ulm_prot %>% transmute(state = source, sample = condition,
                                      score = score) %>%
  mutate(dataset = "Proteomics")
score_all <- bind_rows(score_rna, score_prot) %>%
  left_join(samp_grp, by = "sample") %>%
  mutate(state = factor(state,
                        levels = c("follicular lineage high",
                                   "lineage preserved epithelial",
                                   "plasticity high transition",
                                   "terminal dedifferentiated")),
         dataset = factor(dataset, levels = c("RNA", "Proteomics")))

## group stats
state_stats <- score_all %>% group_by(state, dataset) %>%
  summarise(wt_mean = mean(score[group == "WT host"]),
            fap_mean = mean(score[group == "FAP-deficient host"]),
            delta_fap_minus_wt = fap_mean - wt_mean,
            wilcox_p = tryCatch(
              wilcox.test(score ~ group)$p.value, error = function(e) NA_real_),
            .groups = "drop") %>%
  mutate(lbl = vapply(wilcox_p, p_label, character(1)))

write_tsv(state_stats, file.path(fig_dir, "fig7_states_in_fig8_group_stats.tsv"))

## Panel B — state signature box (RNA + Protein)
p_B <- ggplot(score_all, aes(group, score, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.9,
              shape = 21, stroke = 0.2, colour = "black") +
  geom_text(data = state_stats,
            aes(x = 1.5, y = Inf, label = lbl),
            inherit.aes = FALSE, vjust = 1.3, fontface = "italic",
            size = 6.5/.pt) +
  facet_grid(dataset ~ state, scales = "free_y") +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_x_discrete(labels = c("WT host" = "WT",
                              "FAP-deficient host" = "KO")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = NULL, y = "ULM activity (z)",
       title = "Epithelial-state signatures scored in FAP-KO",
       subtitle = "follicular / lineage-preserved up; plasticity / terminal-dediff down in FAP-KO — matches Fig 3",
       tag = "c") +
  theme_nc() +
  theme(strip.text.x = element_text(size = base_font_size - 1),
        strip.text.y = element_text(size = base_font_size, face = "bold"),
        panel.spacing.x = unit(3, "mm"),
        panel.spacing.y = unit(3, "mm"))

## =======================================================================
## assemble
## =======================================================================
design <- "
AABB
AABB
CCCC
CCCC
"

## In-figure plot title stripped to match scripts 122/123 convention
## (Supplementary figure title lives in the SI docx caption, not on the figure).
fig_s12 <- p_A1 + p_A2 + p_B +
  plot_layout(design = design) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(fig_dir, "FigS22_fig7_fig8_cross_echo_v2_itga5_highlight.pdf")
png_path <- file.path(fig_dir, "FigS22_fig7_fig8_cross_echo_v2_itga5_highlight.png")

## --- vector export at TRUE SI print size (180mm); height reduced after dropping Panel C
dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
if (capabilities("cairo")) grDevices::cairo_pdf(pdf_path, width = 180/25.4, height = 140/25.4) else if (capabilities("aqua")) grDevices::quartz(file = pdf_path, type = "pdf", width = 180/25.4, height = 140/25.4) else grDevices::pdf(pdf_path, width = 180/25.4, height = 140/25.4)
print(fig_s12); dev.off()
ggsave(png_path, fig_s12, width = 180, height = 140, units = "mm", dpi = 600, device = ragg::agg_png)

cat("Fig S12 cross-echo written:\n  ", pdf_path, "\n  ", png_path, "\n")

## Print summary
cat("\n==== A: Fig 7 L-R pairs FC in Fig 8 ====\n")
print(lr_detail %>% dplyr::select(priority_rank, pair_label,
                                  ligand_rna_fc, ligand_prot_fc,
                                  receptor_rna_fc, receptor_prot_fc) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))

cat("\n==== B: Fig 7 state signatures in Fig 8 ====\n")
print(state_stats %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))
