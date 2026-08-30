#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S15 - MOFA2 multi-omic factor extraction.
## Step 3 of multi-omics integration:
## MOFA2 joint factor analysis across RNA + DIA proteomics
## for WT vs FAP-deficient host BPC xenograft (4v4).
##
## Deliverables:
##   * Trained MOFA model (.hdf5)
##   * Factor x Sample score matrix, Factor x Gene weights per modality
##   * Variance-explained plot per modality & factor
##   * Factor vs group (WT/FAP-KO) separation test (Wilcox / t-test)
##   * Factor vs external phenotype scores (TCGA-BRS, TDS, dediff axes)
##   * Top-loading gene enrichment for each significant factor
##     (Hallmark via fgsea + simple ranked-list sanity)
##   * Publication-style 6-panel overview figure


suppressPackageStartupMessages({
  library(MOFA2)
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(fgsea)
})

proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
bulk_out <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
out_dir  <- file.path(bulk_out, "mofa2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_wt  <- "#3B6B9C"; col_fap <- "#C9493A"; col_up <- col_fap; col_dn <- col_wt
col_purple <- "#4B1E70"
grp_levels <- c("WT host", "FAP-deficient host")
grp_cols <- setNames(c(col_wt, col_fap), grp_levels)

base_font_size <- 7
theme_nc <- function(base = base_font_size){
  theme_classic(base_size = base) +
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

to_mgi <- function(x){
  x <- toupper(x); paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
}

## ---- build expression matrices ------------------------------------------
fpkm_tab <- read_tsv(file.path(proj, "001_mRNA_Summary/1.GeneExpression/1_genes_fpkm_expression.txt"),
                     show_col_types = FALSE, progress = FALSE)
prot_tab <- read_tsv(file.path(proj, "002_DIA_Summary/02.ProteinExp/protein_annotation_profile.txt"),
                     show_col_types = FALSE, progress = FALSE)

samp_order <- c("WT_BPC1","WT_BPC2","WT_BPC3","WT_BPC4",
                "FAP_BPC1","FAP_BPC2","FAP_BPC3","FAP_BPC4")
samp_grp <- tibble(sample = samp_order,
                   group = factor(ifelse(grepl("^WT", samp_order),
                                         "WT host", "FAP-deficient host"),
                                  levels = grp_levels))

## RNA: log2(FPKM+1) in MGI symbol space, filter low expression
rna_cols <- paste0("FPKM.", samp_order)
rna_mat <- fpkm_tab %>%
  select(gene_name, all_of(rna_cols)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>%
  slice_max(rowSums(across(all_of(rna_cols))), n = 1, with_ties = FALSE) %>%
  ungroup() %>% column_to_rownames("gene_name") %>% as.matrix()
colnames(rna_mat) <- samp_order
rna_log <- log2(rna_mat + 1)

# Drop near-zero RNA features
rna_keep <- rowSums(rna_mat >= 1) >= 2
rna_log <- rna_log[rna_keep, , drop = FALSE]

# Keep top variable
rna_var <- apply(rna_log, 1, var)
rna_top <- names(sort(rna_var, decreasing = TRUE))[1:min(5000, length(rna_var))]
rna_log_top <- rna_log[rna_top, , drop = FALSE]
rna_log_top <- rna_log_top - rowMeans(rna_log_top) # mean-center by feature

## Protein: log2 intensities; impute NA by per-column min/2
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
# Drop rows with SD 0
prot_log <- prot_log[apply(prot_log, 1, sd) > 0, , drop = FALSE]
# Top variable
prot_var <- apply(prot_log, 1, var)
prot_top <- names(sort(prot_var, decreasing = TRUE))[1:min(3000, length(prot_var))]
prot_log_top <- prot_log[prot_top, , drop = FALSE]
prot_log_top <- prot_log_top - rowMeans(prot_log_top)

cat(sprintf("RNA view: %d features x %d samples\n", nrow(rna_log_top), ncol(rna_log_top)))
cat(sprintf("Prot view: %d features x %d samples\n", nrow(prot_log_top), ncol(prot_log_top)))

## ---- build MOFA object -------------------------------------------------
data_list <- list(RNA = rna_log_top, Proteomics = prot_log_top)
MOFA <- create_mofa(data_list)

# MOFA2 reserves the "group" column; use "host_group" for the WT/FAP label
# and keep a matching "group" column pointing to the single MOFA group
mofa_grp_name <- groups_names(MOFA)
samples_metadata(MOFA) <- data.frame(
  sample = samp_order,
  group  = mofa_grp_name,
  host_group = samp_grp$group
)

## ---- set options --------------------------------------------------------
data_opts  <- get_default_data_options(MOFA)
data_opts$scale_views <- TRUE

model_opts <- get_default_model_options(MOFA)
model_opts$num_factors <- 5  # 8 samples -> keep modest
model_opts$likelihoods <- c(RNA = "gaussian", Proteomics = "gaussian")

train_opts <- get_default_training_options(MOFA)
train_opts$convergence_mode <- "medium"
train_opts$maxiter <- 1500
train_opts$seed <- 42
train_opts$drop_factor_threshold <- 0.01  # drop factors <1% variance

MOFA <- prepare_mofa(MOFA,
                     data_options     = data_opts,
                     model_options    = model_opts,
                     training_options = train_opts)

## ---- train --------------------------------------------------------------
model_hdf5 <- file.path(out_dir, "model.hdf5")
set.seed(42)
MOFA <- run_mofa(MOFA, outfile = model_hdf5, use_basilisk = FALSE)

## ---- extract outputs ----------------------------------------------------
var_expl <- MOFA2::get_variance_explained(MOFA)$r2_per_factor[[1]] %>%
  as.data.frame() %>% rownames_to_column("factor") %>%
  pivot_longer(-factor, names_to = "view", values_to = "pct")

factors_mat <- MOFA2::get_factors(MOFA, factors = "all", as.data.frame = FALSE)[[1]]
factor_df <- as.data.frame(factors_mat) %>% rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "factor", values_to = "value") %>%
  left_join(samp_grp %>% rename(host_group = group), by = "sample") %>%
  mutate(group = host_group)

weights_list <- MOFA2::get_weights(MOFA, views = "all", factors = "all",
                                   as.data.frame = TRUE) %>%
  as_tibble()

write_tsv(var_expl,   file.path(out_dir, "variance_explained_per_factor.tsv"))
write_tsv(factor_df,  file.path(out_dir, "factor_scores_per_sample.tsv"))
write_tsv(weights_list, file.path(out_dir, "factor_weights_all.tsv"))

## ---- factor vs group test ----------------------------------------------
factor_stats <- factor_df %>% group_by(factor) %>%
  summarise(wilcox_p = tryCatch(wilcox.test(value ~ group)$p.value, error = function(e) NA_real_),
            t_p      = tryCatch(t.test(value ~ group)$p.value,      error = function(e) NA_real_),
            delta_fap_minus_wt =
              mean(value[group == "FAP-deficient host"], na.rm = TRUE) -
              mean(value[group == "WT host"], na.rm = TRUE),
            .groups = "drop") %>%
  arrange(wilcox_p)
write_tsv(factor_stats, file.path(out_dir, "factor_group_stats.tsv"))
cat("\n==== factor vs group =====\n")
print(factor_stats)

## ---- link factors to external phenotype scores -------------------------
dediff_scores <- read_tsv(file.path(bulk_out, "dediff_axes/integrated/dediff_axes_scores_per_sample.tsv"),
                          show_col_types = FALSE)
brs_scores    <- read_tsv(file.path(bulk_out, "tcga_brs_official/bulk_rna/bulk_official_tcga_brs_scores.tsv"),
                          show_col_types = FALSE)

phen_wide <- dediff_scores %>%
  filter(dataset == "RNA", !is.na(score_z)) %>%
  select(sample, signature, score_z) %>%
  bind_rows(brs_scores %>%
              transmute(sample, signature = "TCGA-BRS (raw)",
                        score_z = scale(brs_raw)[, 1])) %>%
  pivot_wider(names_from = signature, values_from = score_z)

factor_wide <- factor_df %>%
  select(sample, factor, value) %>%
  pivot_wider(names_from = factor, values_from = value)

joined_phen <- factor_wide %>% inner_join(phen_wide, by = "sample")

factor_cols <- grep("^Factor", colnames(joined_phen), value = TRUE)
phen_cols   <- setdiff(colnames(joined_phen), c("sample", factor_cols))

cor_tab <- expand.grid(factor = factor_cols, phen = phen_cols,
                       stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(rho = suppressWarnings(cor(joined_phen[[factor]],
                                    joined_phen[[phen]],
                                    method = "spearman", use = "pairwise.complete.obs")),
         p   = suppressWarnings(cor.test(joined_phen[[factor]],
                                         joined_phen[[phen]],
                                         method = "spearman",
                                         use = "pairwise.complete.obs",
                                         exact = FALSE)$p.value)) %>%
  ungroup() %>%
  mutate(q = p.adjust(p, method = "BH"))
write_tsv(cor_tab, file.path(out_dir, "factor_phenotype_correlations.tsv"))

## ---- enrichment of Hallmark on top weights -----------------------------
gmt_path <- file.path(proj, "reference/mh.all.v2026.1.Mm.symbols.gmt")
hallmark <- if (file.exists(gmt_path)){
  lines <- readLines(gmt_path)
  lapply(strsplit(lines, "\t"), function(x){
    if (length(x) >= 3) setNames(list(x[-(1:2)]), x[1]) else NULL
  }) %>% purrr::list_flatten()
} else list()

## strip MOFA-added view suffixes from feature names
strip_view_suffix <- function(x){
  sub("_(RNA|Proteomics)$", "", x)
}

run_fgsea_factor <- function(wt_df, view_sel){
  w <- wt_df %>% filter(view == view_sel) %>%
    mutate(gene = strip_view_suffix(feature),
           gene = if (view_sel == "Proteomics") to_mgi(gene) else gene)
  ranks <- w$value; names(ranks) <- w$gene
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- ranks[is.finite(ranks)]
  if (length(hallmark) == 0 || length(ranks) < 20) return(NULL)
  set.seed(1)
  fgsea::fgseaMultilevel(pathways = hallmark, stats = ranks,
                         minSize = 10, maxSize = 500) %>%
    as_tibble() %>% mutate(view = view_sel)
}

enrich_all <- list()
for (fac in factor_cols){
  wt_f <- weights_list %>% filter(factor == fac)
  if (nrow(wt_f) == 0) next
  e_rna  <- run_fgsea_factor(wt_f, "RNA")
  e_prot <- run_fgsea_factor(wt_f, "Proteomics")
  enrich_all[[fac]] <- bind_rows(e_rna, e_prot) %>% mutate(factor = fac)
}
enrich_tbl <- bind_rows(enrich_all)
if (nrow(enrich_tbl) > 0) write_tsv(enrich_tbl %>% select(-leadingEdge),
                                    file.path(out_dir, "factor_hallmark_fgsea.tsv"))

## ---- visualization -----------------------------------------------------
## Panel A: variance explained heatmap
p_A <- ggplot(var_expl, aes(view, factor, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            size = (base_font_size - 1.5)/.pt) +
  scale_fill_gradient(low = "white", high = col_purple,
                      name = "Var. expl.\n(%)") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "MOFA2 variance explained",
       subtitle = "Per factor x modality",
       tag = "A") +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key.height = unit(6, "mm"),
        legend.key.width  = unit(2.6, "mm"))

## Panel B: factor scores x sample (boxplot/ strip)
order_f <- factor_stats$factor
factor_df <- factor_df %>% mutate(factor = factor(factor, levels = order_f))

p_lab <- factor_stats %>%
  mutate(label = ifelse(wilcox_p < 0.05,
                        sprintf("P = %.3f*", wilcox_p),
                        sprintf("P = %.2f", wilcox_p)),
         factor = factor(factor, levels = order_f))

p_B <- ggplot(factor_df, aes(group, value, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.9,
              shape = 21, stroke = 0.2, colour = "black") +
  geom_text(data = p_lab, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.3,
            size = (base_font_size - 1.5)/.pt) +
  facet_wrap(~ factor, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_x_discrete(labels = c("WT host" = "WT", "FAP-deficient host" = "KO")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "Factor value",
       title = "Factor scores by host group",
       tag = "B") +
  theme_nc()

## Panel C: factor x phenotype correlation heatmap
cor_plot <- cor_tab %>%
  mutate(factor = factor(factor, levels = rev(factor_cols)),
         phen   = factor(phen, levels = phen_cols),
         rho_capped = pmax(pmin(rho, 1), -1),
         sig_tag = case_when(q < 0.05 ~ "*", p < 0.05 ~ ".", TRUE ~ ""))

p_C <- ggplot(cor_plot, aes(phen, factor, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sig_tag),
            size = (base_font_size + 0.5)/.pt, vjust = 0.7) +
  scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman\nrho") +
  labs(x = NULL, y = NULL,
       title = "Factor vs external phenotype scores",
       subtitle = "* FDR<0.05, . P<0.05",
       tag = "C") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.key.height = unit(7, "mm"),
        legend.key.width  = unit(2.6, "mm"))

## Panel D: top weights for factor with strongest group-separation
best_factor <- factor_stats$factor[1]
w_best <- weights_list %>% filter(factor == best_factor) %>%
  mutate(feature_clean = strip_view_suffix(feature)) %>%
  group_by(view) %>%
  arrange(desc(abs(value))) %>%
  slice_head(n = 15) %>% ungroup() %>%
  mutate(feature_clean = factor(feature_clean,
                                levels = unique(feature_clean[order(value)])),
         view = factor(view, levels = c("RNA", "Proteomics")))

p_D <- ggplot(w_best, aes(value, feature_clean, fill = value)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_col(colour = "black", linewidth = 0.2) +
  facet_wrap(~ view, nrow = 1, scales = "free_y") +
  scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                       midpoint = 0, guide = "none",
                       limits = c(-max(abs(w_best$value)),
                                   max(abs(w_best$value)))) +
  labs(x = "Weight", y = NULL,
       title = sprintf("Top 15 features loading %s (per view)", best_factor),
       subtitle = sprintf("Group Wilcox P = %.3f",
                          factor_stats$wilcox_p[factor_stats$factor == best_factor]),
       tag = "D") +
  theme_nc() +
  theme(axis.text.y = element_text(size = base_font_size - 1.5,
                                   face = "italic"),
        panel.spacing.x = unit(2, "mm"))

## Panel E: Hallmark enrichment bubble for the best factor
if (nrow(enrich_tbl) > 0){
  e_top <- enrich_tbl %>% filter(factor == best_factor) %>%
    arrange(padj) %>% group_by(view) %>% slice_head(n = 10) %>% ungroup()
  if (nrow(e_top) > 0){
    e_top <- e_top %>%
      mutate(pathway = sub("^HALLMARK_", "", pathway),
             pathway_ordered = factor(pathway,
                                      levels = unique(pathway[order(NES)])))
    p_E <- ggplot(e_top, aes(NES, pathway_ordered,
                             size = -log10(pmax(padj, 1e-10)),
                             fill = NES)) +
      geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
      geom_point(shape = 21, colour = "black", stroke = 0.25) +
      facet_wrap(~ view, nrow = 1, scales = "free") +
      scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                           midpoint = 0, guide = "none") +
      scale_size_continuous(range = c(1.5, 5),
                            name = "-log10(padj)") +
      labs(x = "NES", y = NULL,
           title = sprintf("Hallmark enrichment of %s weights",
                           best_factor), tag = "E") +
      theme_nc() +
      theme(axis.text.y = element_text(size = base_font_size - 1.5),
            panel.spacing.x = unit(2, "mm"))
  } else {
    p_E <- ggplot() + theme_void() +
      annotate("text", x = 1, y = 1, label = "No Hallmark hits")
  }
} else {
  p_E <- ggplot() + theme_void() +
    annotate("text", x = 1, y = 1, label = "Hallmark GMT not found")
}

## Panel F: RNA vs Protein weight concordance on the best factor
w_match <- weights_list %>%
  filter(factor == best_factor) %>%
  mutate(gene = strip_view_suffix(feature),
         gene = ifelse(view == "Proteomics", to_mgi(gene), gene)) %>%
  select(view, gene, value) %>%
  pivot_wider(names_from = view, values_from = value,
              values_fn = mean) %>%
  filter(!is.na(RNA), !is.na(Proteomics))

if (nrow(w_match) >= 5){
  rho_w <- suppressWarnings(cor(w_match$RNA, w_match$Proteomics, method = "spearman"))
  top_lab <- w_match %>%
    mutate(score = abs(RNA) + abs(Proteomics)) %>%
    arrange(desc(score)) %>% slice_head(n = 15) %>% pull(gene)
  p_F <- ggplot(w_match, aes(RNA, Proteomics)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                linewidth = 0.25, colour = "grey70") +
    geom_point(colour = col_purple, size = 0.8, alpha = 0.6) +
    geom_text_repel(data = subset(w_match, gene %in% top_lab),
                    aes(label = gene), size = (base_font_size - 1.5)/.pt,
                    min.segment.length = 0, segment.size = 0.2,
                    max.overlaps = Inf, box.padding = 0.3,
                    fontface = "italic") +
    labs(x = "Weight (RNA view)",
         y = "Weight (Proteomics view)",
         title = sprintf("Cross-view weight alignment on %s", best_factor),
         subtitle = sprintf("Spearman %.2f (n = %d shared features)",
                            rho_w, nrow(w_match)),
         tag = "F") +
    theme_nc()
} else {
  p_F <- ggplot() + theme_void() +
    annotate("text", x = 1, y = 1, label = "Insufficient shared features")
}

design <- "
AABB
CCBB
DDEE
DDEE
FFEE
"
fig <- p_A + p_B + p_C + p_D + p_E + p_F +
  plot_layout(design = design) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(out_dir, "mofa2_overview.pdf")
png_path <- file.path(out_dir, "mofa2_overview.png")
ggsave(pdf_path, fig, width = 240, height = 260, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(png_path, fig, width = 240, height = 260, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("\nMOFA2 outputs written to:\n  ", out_dir, "\n")
