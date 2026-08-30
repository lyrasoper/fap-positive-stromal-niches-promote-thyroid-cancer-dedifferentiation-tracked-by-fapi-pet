#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S15 - NicheNet ligand activity on the bulk contrast.
## NicheNet stroma -> tumor ligand-target analysis for the
## WT vs FAP-deficient host BPC xenograft (bulk RNA + DIA proteomics).
##
## Biological framing:
##   * WT host has an intact FAP+ CAF niche that releases stromal ligands.
##   * FAP-deficient host loses part of that ligand pool.
##   * NicheNet ranks which stromal ligands lost in FAP-KO best explain
##     the tumor gene programs collapsed in FAP-KO (stroma-dependent program)
##     and separately, the tumor programs released in FAP-KO
##     (e.g., thyroid lineage re-expression).
##   * Every NicheNet-prioritized ligand / target is cross-checked against
##     the orthogonal proteomics delta.


suppressPackageStartupMessages({
  library(nichenetr); library(dplyr); library(tidyr); library(readr)
  library(ggplot2);   library(patchwork); library(ggrepel); library(scales)
  library(circlize);  library(ComplexHeatmap); library(grid)
})

## ---- paths ---------------------------------------------------------------
proj     <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
prior_dir<- file.path(proj, "reference/nichenet_mouse")
out_root <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407/nichenet")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

## ---- palette -------------------------------------------------------------
col_wt  <- "#3B6B9C"; col_fap <- "#C9493A"; col_grey <- "#B7B7B7"
col_up  <- col_fap;   col_dn  <- col_wt

base_font_size <- 7
theme_nc <- function(base = base_font_size){
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = base + 1, hjust = 0),
      plot.subtitle = element_text(size = base - 0.5, colour = "grey30"),
      plot.tag = element_text(face = "bold", size = base + 3, hjust = 0),
      axis.title = element_text(size = base),
      axis.text  = element_text(size = base - 0.5, colour = "black"),
      axis.line  = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      strip.background = element_blank(),
      strip.text = element_text(size = base),
      legend.key.size = unit(3.2, "mm"),
      legend.text = element_text(size = base - 0.5),
      legend.title = element_text(size = base - 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(3, 4, 3, 4)
    )
}

## ---- load priors & bulk data --------------------------------------------
message("Loading NicheNet mouse priors ...")
ligand_target_matrix <- readRDS(file.path(prior_dir, "ligand_target_matrix_nsga2r_final_mouse.rds"))
lr_network           <- readRDS(file.path(prior_dir, "lr_network_mouse_21122021.rds"))
weighted_networks    <- readRDS(file.path(prior_dir, "weighted_networks_nsga2r_final_mouse.rds"))
weighted_networks_lr <- weighted_networks$lr_sig %>%
  dplyr::inner_join(lr_network %>% distinct(from, to), by = c("from","to"))

bulk_out <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
deg_rna  <- read_tsv(file.path(bulk_out, "rna/rna_differential_expression.tsv"), show_col_types = FALSE)
deg_prot <- read_tsv(file.path(bulk_out, "proteomics/proteomics_differential_expression.tsv"), show_col_types = FALSE)

fpkm_tab <- read_tsv(file.path(proj, "001_mRNA_Summary/1.GeneExpression/1_genes_fpkm_expression.txt"),
                     show_col_types = FALSE, progress = FALSE)

## ---- normalize gene symbols: convert bulk to NicheNet mouse-symbol space ----
## NicheNet mouse prior uses MGI-style first-letter-capital symbols (Fap, Col1a1).
## Our RNA table already uses MGI symbols (e.g. Fap, Pax8). Proteomics uses
## UPPER symbols (FAP, COL1A1); convert to MGI style by capitalizing first letter.
to_mgi <- function(x){
  x <- toupper(x)
  paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
}
deg_prot <- deg_prot %>% mutate(gene_mgi = to_mgi(gene))

## ---- define background & genesets ---------------------------------------
## background = all genes with testable RNA data AND present in the prior rows
bg_genes <- deg_rna$gene %>% unique()
bg_genes <- bg_genes[bg_genes %in% rownames(ligand_target_matrix)]

## ligand pool = NicheNet ligands (rows of the matrix) also expressed in bulk
ligand_pool <- colnames(ligand_target_matrix)

## Strong RNA DEGs (FDR<0.05), signed
deg_strong <- deg_rna %>% filter(significant == "yes") %>%
  mutate(gene = as.character(gene))
dn_genes <- deg_strong %>% filter(regulation == "down") %>% pull(gene)  # lost in FAP-KO
up_genes <- deg_strong %>% filter(regulation == "up")   %>% pull(gene)  # gained in FAP-KO

## Keep only genes present in prior matrix background for the geneset of interest
gs_dn <- intersect(dn_genes, rownames(ligand_target_matrix))
gs_up <- intersect(up_genes, rownames(ligand_target_matrix))

## Potential ligands:
##  - Set A (stroma-lost ligands): ligands DOWN in FAP-KO (FDR<0.05 RNA).
##    Used to explain which genes go DOWN.
##  - Set B (stroma-released ligands): ligands UP in FAP-KO (inhibitors
##    whose loss frees the tumor). Used to explain which genes go UP.
potential_ligands_down <- intersect(ligand_pool,
                                    deg_strong %>% filter(regulation == "down") %>% pull(gene))
potential_ligands_up   <- intersect(ligand_pool,
                                    deg_strong %>% filter(regulation == "up")   %>% pull(gene))

message(sprintf("Background genes: %d | DN geneset: %d | UP geneset: %d",
                length(bg_genes), length(gs_dn), length(gs_up)))
message(sprintf("Potential DOWN ligands (stromal-lost): %d | UP ligands: %d",
                length(potential_ligands_down), length(potential_ligands_up)))

## ---- core: predict ligand activities ------------------------------------
run_activity <- function(geneset, potential_ligands, label){
  if (length(geneset) < 5 || length(potential_ligands) < 1){
    message(sprintf("[%s] insufficient input, skipped", label))
    return(NULL)
  }
  res <- predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = bg_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  ) %>% arrange(desc(aupr_corrected)) %>%
    mutate(rank = row_number(), analysis = label)
  res
}

act_dn_by_dnL <- run_activity(gs_dn, potential_ligands_down,
                              "stroma-lost ligands -> collapsed tumor program")
act_up_by_upL <- run_activity(gs_up, potential_ligands_up,
                              "stroma-up ligands -> released tumor program")

## Cross-alternative: use ALL ligands in prior that are also detected in bulk
## (regardless of DE sign) as the potential_ligands set.  This is the classic
## NicheNet bulk setting; it lets us see which ligand wins when we do not
## pre-filter by DE direction.
potential_all <- intersect(ligand_pool, bg_genes)
act_dn_any <- run_activity(gs_dn, potential_all,
                           "any detected ligand -> collapsed tumor program")
act_up_any <- run_activity(gs_up, potential_all,
                           "any detected ligand -> released tumor program")

all_activities <- bind_rows(act_dn_by_dnL, act_up_by_upL, act_dn_any, act_up_any)
write_tsv(all_activities, file.path(out_root, "ligand_activities_all.tsv"))

## ---- pick the primary analysis for figure panels ------------------------
## Primary = stroma-lost ligand explains collapsed program
primary <- act_dn_by_dnL
stopifnot(!is.null(primary))

## Cross-check each top ligand with proteomics log2FC
prot_fc_lookup <- deg_prot %>%
  transmute(gene_mgi, prot_log2fc = log2_fc, prot_p = p_value) %>%
  distinct(gene_mgi, .keep_all = TRUE)

rna_fc_lookup <- deg_rna %>%
  transmute(gene, rna_log2fc = log2_fc, rna_q = q_value) %>%
  distinct(gene, .keep_all = TRUE)

primary <- primary %>%
  left_join(rna_fc_lookup,  by = c("test_ligand" = "gene")) %>%
  left_join(prot_fc_lookup, by = c("test_ligand" = "gene_mgi"))

write_tsv(primary, file.path(out_root, "primary_stromalost_ligand_ranking.tsv"))

## ---- panel A: ligand-activity bar chart (top 25) ------------------------
top_n <- 25
top_ligands <- primary %>% slice_head(n = top_n) %>%
  mutate(test_ligand = factor(test_ligand, levels = rev(test_ligand)))

p_A <- ggplot(top_ligands, aes(aupr_corrected, test_ligand, fill = rna_log2fc)) +
  geom_col(colour = "black", linewidth = 0.25) +
  scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                       midpoint = 0, limits = c(-3, 3), oob = squish,
                       name = "RNA log2FC\n(FAP-KO vs WT)") +
  labs(x = "Ligand activity (AUPR corrected)",
       y = NULL,
       title = sprintf("Top %d stroma-lost ligands explaining\nthe collapsed tumor program in FAP-KO", top_n),
       subtitle = sprintf("Geneset: %d DOWN genes (FDR<0.05) | Candidate ligands: %d",
                          length(gs_dn), length(potential_ligands_down)),
       tag = "A") +
  theme_nc()

## ---- panel B: ligand x target regulatory-potential heatmap --------------
top_ligand_ids <- as.character(top_ligands$test_ligand)
ligand_target_df <- purrr::map_dfr(top_ligand_ids, function(lg){
  get_weighted_ligand_target_links(
    ligand = lg, geneset = gs_dn,
    ligand_target_matrix = ligand_target_matrix, n = 200
  )
}) %>% dplyr::filter(!is.na(target))

## pick the N most-frequently-targeted genes across top ligands
top_targets <- ligand_target_df %>%
  group_by(target) %>% summarise(n = n(), sum_w = sum(weight)) %>%
  arrange(desc(sum_w)) %>% slice_head(n = 50) %>% pull(target)

lt_mat <- ligand_target_df %>%
  filter(target %in% top_targets, !is.na(target), !is.na(ligand)) %>%
  pivot_wider(id_cols = target, names_from = ligand,
              values_from = weight, values_fill = 0,
              values_fn = max) %>%
  tibble::column_to_rownames("target") %>% as.matrix()
cat("lt_mat raw dim:", dim(lt_mat), "\n")
## order ligands by activity, targets by row-sum
lt_mat <- lt_mat[, intersect(top_ligand_ids, colnames(lt_mat)), drop = FALSE]
cat("lt_mat after ligand filter dim:", dim(lt_mat), "\n")
if (nrow(lt_mat) > 0 && ncol(lt_mat) > 0){
  lt_mat <- lt_mat[order(rowSums(lt_mat), decreasing = TRUE), , drop = FALSE]
  keep_n <- min(40, nrow(lt_mat))
  if (keep_n >= 1) lt_mat <- lt_mat[seq_len(keep_n), , drop = FALSE]
}

lt_long <- as.data.frame(lt_mat) %>% tibble::rownames_to_column("target") %>%
  pivot_longer(-target, names_to = "ligand", values_to = "weight") %>%
  mutate(target = factor(target, levels = rev(rownames(lt_mat))),
         ligand = factor(ligand, levels = colnames(lt_mat)))

p_B <- ggplot(lt_long, aes(ligand, target, fill = weight)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "white", high = "#4B1E70",
                      name = "Regulatory\npotential") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Ligand -> target regulatory potential (stroma-lost axis)",
       tag = "B") +
  theme_nc() +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0,
                                       size = base_font_size - 0.5),
        axis.text.y = element_text(size = base_font_size - 1.5),
        axis.ticks = element_blank(), axis.line = element_blank(),
        legend.position = "right",
        legend.key.height = unit(8, "mm"),
        legend.key.width  = unit(2.6, "mm"))

## ---- panel C: ligand RNA vs protein concordance ------------------------
lig_cross <- top_ligands %>%
  transmute(ligand = as.character(test_ligand),
            activity = aupr_corrected,
            rna_log2fc, prot_log2fc) %>%
  mutate(concordant = !is.na(prot_log2fc) & sign(rna_log2fc) == sign(prot_log2fc))

p_C <- ggplot(lig_cross, aes(rna_log2fc, prot_log2fc)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              linewidth = 0.25, colour = "grey70") +
  geom_point(aes(size = activity, fill = concordant),
             shape = 21, colour = "black", stroke = 0.3) +
  geom_text_repel(aes(label = ligand), size = (base_font_size - 1.5)/.pt,
                  min.segment.length = 0, max.overlaps = Inf,
                  segment.size = 0.25, box.padding = 0.3) +
  scale_fill_manual(values = c(`TRUE` = "#4B1E70", `FALSE` = "grey85"),
                    labels = c(`TRUE` = "RNA/Protein same sign",
                               `FALSE` = "protein NA or opposite"),
                    name = NULL) +
  scale_size_continuous(range = c(1.6, 5), name = "Ligand AUPR") +
  labs(x = "Ligand RNA log2FC (FAP-KO vs WT)",
       y = "Ligand protein log2FC (FAP-KO vs WT)",
       title = "RNA <-> protein concordance of top stroma-lost ligands",
       tag = "C") +
  theme_nc() +
  theme(legend.position = "right")

## ---- panel D: target RNA vs protein concordance ------------------------
tgt_cross <- tibble(target = rownames(lt_mat)) %>%
  left_join(rna_fc_lookup,  by = c("target" = "gene")) %>%
  left_join(prot_fc_lookup, by = c("target" = "gene_mgi")) %>%
  filter(!is.na(prot_log2fc))

if (nrow(tgt_cross) >= 2){
  cor_s <- suppressWarnings(cor(tgt_cross$rna_log2fc,
                                tgt_cross$prot_log2fc,
                                method = "spearman",
                                use = "pairwise.complete.obs"))
  p_D <- ggplot(tgt_cross, aes(rna_log2fc, prot_log2fc)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                linewidth = 0.25, colour = "grey70") +
    geom_point(size = 1.5, fill = col_dn, shape = 21, stroke = 0.2, colour = "black") +
    geom_text_repel(aes(label = target), size = (base_font_size - 1.5)/.pt,
                    max.overlaps = Inf, segment.size = 0.25,
                    min.segment.length = 0) +
    labs(x = "Target RNA log2FC", y = "Target protein log2FC",
         title = "Top NicheNet targets: RNA <-> protein",
         subtitle = sprintf("Spearman rho = %.2f (n = %d)", cor_s, nrow(tgt_cross)),
         tag = "D") +
    theme_nc()
} else {
  p_D <- ggplot() + theme_void() +
    annotate("text", x = 1, y = 1, label = "Insufficient target-protein overlap")
}

## ---- panel E: ligand -> receptor weighted links for best 8 ligands ------
best_ligs <- head(top_ligand_ids, 12)
rec_pool <- intersect(lr_network$to, bg_genes)
lr_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands = best_ligs,
  expressed_receptors   = rec_pool,
  lr_network            = lr_network,
  weighted_networks_lr_sig = weighted_networks_lr
) %>% filter(weight > 0.08)

if (nrow(lr_df) > 0){
  lr_mat <- lr_df %>%
    filter(!is.na(to), !is.na(from)) %>%
    pivot_wider(id_cols = to, names_from = from,
                values_from = weight, values_fill = 0,
                values_fn = max) %>%
    tibble::column_to_rownames("to") %>% as.matrix()
  lr_mat <- lr_mat[, intersect(best_ligs, colnames(lr_mat)), drop = FALSE]
  if (nrow(lr_mat) > 0 && ncol(lr_mat) > 0){
    lr_mat <- lr_mat[order(rowSums(lr_mat), decreasing = TRUE), , drop = FALSE]
    keep_r <- min(30, nrow(lr_mat))
    if (keep_r >= 1) lr_mat <- lr_mat[seq_len(keep_r), , drop = FALSE]
  }

  lr_long <- as.data.frame(lr_mat) %>% tibble::rownames_to_column("receptor") %>%
    pivot_longer(-receptor, names_to = "ligand", values_to = "weight") %>%
    mutate(receptor = factor(receptor, levels = rev(rownames(lr_mat))),
           ligand   = factor(ligand, levels = colnames(lr_mat)))

  p_E <- ggplot(lr_long, aes(ligand, receptor, fill = weight)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient(low = "white", high = "#165A42",
                        name = "LR weight") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL,
         title = "Ligand -> receptor prior weight (top 12 stroma-lost ligands)",
         tag = "E") +
    theme_nc() +
    theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0,
                                         size = base_font_size - 0.5),
          axis.text.y = element_text(size = base_font_size - 1.5),
          axis.line = element_blank(), axis.ticks = element_blank(),
          legend.position = "right",
          legend.key.height = unit(8, "mm"),
          legend.key.width  = unit(2.6, "mm"))
} else {
  p_E <- ggplot() + theme_void()
}

## ---- panel F: SECONDARY analysis — what explains the redifferentiation? ------
act_up_any_df <- act_up_any
if (!is.null(act_up_any_df)){
  top_up <- act_up_any_df %>% slice_head(n = 15) %>%
    left_join(rna_fc_lookup,  by = c("test_ligand" = "gene")) %>%
    left_join(prot_fc_lookup, by = c("test_ligand" = "gene_mgi")) %>%
    mutate(test_ligand = factor(test_ligand, levels = rev(test_ligand)))

  p_F <- ggplot(top_up, aes(aupr_corrected, test_ligand, fill = rna_log2fc)) +
    geom_col(colour = "black", linewidth = 0.25) +
    scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                         midpoint = 0, limits = c(-3, 3), oob = squish,
                         name = "RNA log2FC") +
    labs(x = "Ligand activity (AUPR corrected)", y = NULL,
         title = "Ligands explaining the released (redifferentiation) program",
         subtitle = sprintf("Geneset: %d UP genes | Candidate ligands: %d",
                            length(gs_up), length(potential_all)),
         tag = "F") +
    theme_nc()
} else {
  p_F <- ggplot() + theme_void()
}

## ---- assemble figure ----------------------------------------------------
design <- "
AABB
AABB
CCBB
DDEE
DDEE
FFEE
"
fig_nn <- p_A + p_B + p_C + p_D + p_E + p_F +
  plot_layout(design = design) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(out_root, "nichenet_stroma_tumor_overview.pdf")
png_path <- file.path(out_root, "nichenet_stroma_tumor_overview.png")
ggsave(pdf_path, fig_nn, width = 220, height = 260, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(png_path, fig_nn, width = 220, height = 260, units = "mm",
       dpi = 600, device = ragg::agg_png)

## ---- also dump tidy ligand-target table for sharing ---------------------
write_tsv(ligand_target_df,
          file.path(out_root, "top_ligand_target_weighted_links.tsv"))
write_tsv(lig_cross,
          file.path(out_root, "top_ligand_rna_protein_concordance.tsv"))
if (exists("tgt_cross")) write_tsv(tgt_cross,
                                   file.path(out_root, "top_target_rna_protein_concordance.tsv"))
if (exists("lr_df"))     write_tsv(lr_df,
                                   file.path(out_root, "top_ligand_receptor_weighted_links.tsv"))

cat("\nNicheNet outputs written to:\n  ", out_root, "\n")
