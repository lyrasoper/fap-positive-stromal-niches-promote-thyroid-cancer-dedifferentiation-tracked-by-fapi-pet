#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S14 - decoupleR ULM cross-omics scoring.
## decoupleR-based cross-modality pathway & TF activity analysis
## WT vs FAP-deficient host BPC xenograft (bulk RNA + DIA proteomics).
##
## Produces, on BOTH modalities using consistent scoring:
##   * PROGENy 14-pathway activities (mouse)
##   * MSigDB Hallmark activities (mouse, local .gmt)
##   * CollecTRI TF activities (mouse) if available
## Then integrates via delta-of-delta (RNA vs Protein shift agreement).


suppressPackageStartupMessages({
  library(decoupleR); library(progeny)
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
})

## ---- paths ---------------------------------------------------------------
proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
out_root <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407/decoupler_crossomics")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

col_wt  <- "#3B6B9C"; col_fap <- "#C9493A"; col_up <- col_fap; col_dn <- col_wt
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

## ---- utility -------------------------------------------------------------
to_mgi <- function(x){
  x <- toupper(x)
  paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
}

## ---- build expression matrices ------------------------------------------
fpkm_tab <- read_tsv(file.path(proj, "001_mRNA_Summary/1.GeneExpression/1_genes_fpkm_expression.txt"),
                     show_col_types = FALSE, progress = FALSE)
prot_tab <- read_tsv(file.path(proj, "002_DIA_Summary/02.ProteinExp/protein_annotation_profile.txt"),
                     show_col_types = FALSE, progress = FALSE)

samp_order <- c("WT_BPC1","WT_BPC2","WT_BPC3","WT_BPC4",
                "FAP_BPC1","FAP_BPC2","FAP_BPC3","FAP_BPC4")

## RNA matrix: FPKM -> log2(FPKM+1), aggregate duplicates by max mean, keep MGI symbol
rna_cols <- paste0("FPKM.", samp_order)
rna_mat <- fpkm_tab %>%
  select(gene_name, all_of(rna_cols)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>% slice_max(rowSums(across(all_of(rna_cols))), n = 1, with_ties = FALSE) %>%
  ungroup() %>% column_to_rownames("gene_name") %>% as.matrix()
colnames(rna_mat) <- samp_order
rna_log <- log2(rna_mat + 1)

## Protein matrix: keep as-is (intensities), log2, then impute NA by column min/2
prot_mat <- prot_tab %>%
  select(gene_name, all_of(samp_order)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  mutate(across(all_of(samp_order), ~ suppressWarnings(as.numeric(.)))) %>%
  group_by(gene_name) %>% summarise(across(all_of(samp_order), ~ mean(.x, na.rm = TRUE)),
                                    .groups = "drop") %>%
  mutate(gene_mgi = to_mgi(gene_name)) %>%
  distinct(gene_mgi, .keep_all = TRUE) %>% select(-gene_name) %>%
  column_to_rownames("gene_mgi") %>% as.matrix()
prot_log <- log2(prot_mat + 1)
prot_log[!is.finite(prot_log)] <- NA

col_min_half <- apply(prot_log, 2, function(x) min(x, na.rm = TRUE) - 1)
for (j in seq_len(ncol(prot_log))){
  prot_log[is.na(prot_log[, j]), j] <- col_min_half[j]
}

cat(sprintf("RNA matrix %d genes x %d samples; Protein matrix %d genes x %d samples\n",
            nrow(rna_log), ncol(rna_log), nrow(prot_log), ncol(prot_log)))

## group metadata
samp_grp <- tibble(sample = samp_order,
                   group = factor(ifelse(grepl("^WT", samp_order), "WT host", "FAP-deficient host"),
                                  levels = grp_levels))

## ---- PROGENy 14 pathways -----------------------------------------------
progeny_net <- progeny::getModel(organism = "Mouse", top = 100) %>%
  as.data.frame() %>% rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "source", values_to = "weight") %>%
  filter(weight != 0) %>% mutate(mor = sign(weight), likelihood = 1)

run_progeny <- function(mat, tag){
  # decoupleR wmean on PROGENy
  res <- run_wmean(mat = mat, net = progeny_net, .source = "source",
                   .target = "gene", .mor = "weight",
                   times = 500, minsize = 3)
  res %>% filter(statistic == "norm_wmean") %>%
    transmute(pathway = source, sample = condition, score,
              dataset = tag)
}

progeny_rna  <- run_progeny(rna_log,  "RNA")
progeny_prot <- run_progeny(prot_log, "Proteomics")
progeny_all  <- bind_rows(progeny_rna, progeny_prot) %>%
  left_join(samp_grp, by = "sample")

progeny_delta <- progeny_all %>%
  group_by(dataset, pathway, group) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_score) %>%
  mutate(delta_fap_minus_wt = `FAP-deficient host` - `WT host`)

# Wilcoxon per modality
progeny_stats <- progeny_all %>%
  group_by(dataset, pathway) %>%
  summarise(wilcox_p = wilcox.test(score ~ group)$p.value, .groups = "drop")
progeny_delta <- progeny_delta %>% left_join(progeny_stats,
                                              by = c("dataset", "pathway"))

write_tsv(progeny_all,   file.path(out_root, "progeny_scores_per_sample.tsv"))
write_tsv(progeny_delta, file.path(out_root, "progeny_delta_with_pvalues.tsv"))

## ---- CollecTRI TF activity (optional) ----------------------------------
collectri_net <- try(decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE), silent = TRUE)
if (!inherits(collectri_net, "try-error") && nrow(collectri_net) > 0){
  run_collectri <- function(mat, tag){
    res <- run_ulm(mat = mat, net = collectri_net,
                   .source = "source", .target = "target", .mor = "mor",
                   minsize = 5)
    res %>% transmute(tf = source, sample = condition, score, p_value,
                      dataset = tag)
  }
  tf_rna  <- run_collectri(rna_log,  "RNA")
  tf_prot <- run_collectri(prot_log, "Proteomics")
  tf_all  <- bind_rows(tf_rna, tf_prot) %>% left_join(samp_grp, by = "sample")
  tf_delta <- tf_all %>% group_by(dataset, tf, group) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = mean_score) %>%
    mutate(delta_fap_minus_wt = `FAP-deficient host` - `WT host`)
  tf_stats <- tf_all %>% group_by(dataset, tf) %>%
    summarise(wilcox_p = tryCatch(wilcox.test(score ~ group)$p.value, error = function(e) NA_real_),
              .groups = "drop")
  tf_delta <- tf_delta %>% left_join(tf_stats, by = c("dataset", "tf"))
  write_tsv(tf_all,   file.path(out_root, "collectri_tf_scores_per_sample.tsv"))
  write_tsv(tf_delta, file.path(out_root, "collectri_tf_delta_with_pvalues.tsv"))
} else {
  tf_delta <- NULL
  message("CollecTRI unavailable (likely no network). Skipping TF activity.")
}

## ---- MSigDB Hallmark (via local gmt) -----------------------------------
gmt_path <- file.path(proj, "reference/mh.all.v2026.1.Mm.symbols.gmt")
hallmark_net <- NULL
if (file.exists(gmt_path)){
  lines <- readLines(gmt_path)
  hm_rows <- lapply(strsplit(lines, "\t"), function(x){
    if (length(x) < 3) return(NULL)
    data.frame(source = x[1], target = x[-(1:2)],
               mor = 1, stringsAsFactors = FALSE)
  })
  hallmark_net <- do.call(rbind, hm_rows)
}

if (!is.null(hallmark_net)){
  run_hallmark <- function(mat, tag){
    res <- run_ulm(mat = mat, net = hallmark_net,
                   .source = "source", .target = "target", .mor = "mor",
                   minsize = 5)
    res %>% transmute(hallmark = source, sample = condition, score, p_value,
                      dataset = tag)
  }
  hm_rna  <- run_hallmark(rna_log,  "RNA")
  hm_prot <- run_hallmark(prot_log, "Proteomics")
  hm_all  <- bind_rows(hm_rna, hm_prot) %>% left_join(samp_grp, by = "sample")
  hm_delta <- hm_all %>% group_by(dataset, hallmark, group) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = mean_score) %>%
    mutate(delta_fap_minus_wt = `FAP-deficient host` - `WT host`)
  hm_stats <- hm_all %>% group_by(dataset, hallmark) %>%
    summarise(wilcox_p = tryCatch(wilcox.test(score ~ group)$p.value, error = function(e) NA_real_),
              .groups = "drop")
  hm_delta <- hm_delta %>% left_join(hm_stats, by = c("dataset", "hallmark"))
  write_tsv(hm_all,   file.path(out_root, "hallmark_scores_per_sample.tsv"))
  write_tsv(hm_delta, file.path(out_root, "hallmark_delta_with_pvalues.tsv"))
} else {
  hm_delta <- NULL
  message("MSigDB Hallmark .gmt not found; skipping hallmark activity")
}

## ---- visualize: cross-modality delta concordance for PROGENy ------------
progeny_wide <- progeny_delta %>%
  select(dataset, pathway, delta_fap_minus_wt) %>%
  pivot_wider(names_from = dataset, values_from = delta_fap_minus_wt) %>%
  filter(!is.na(RNA) & !is.na(Proteomics))

plot_delta_concordance <- function(wide_df, label_col = "pathway", title_str){
  spear <- suppressWarnings(cor(wide_df$RNA, wide_df$Proteomics,
                                method = "spearman", use = "pairwise.complete.obs"))
  pear  <- suppressWarnings(cor(wide_df$RNA, wide_df$Proteomics,
                                method = "pearson",  use = "pairwise.complete.obs"))
  ggplot(wide_df, aes(RNA, Proteomics)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                linewidth = 0.25, colour = "grey70") +
    geom_point(aes(fill = RNA), shape = 21, size = 2.4, stroke = 0.25, colour = "black") +
    scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                         midpoint = 0, name = "RNA delta", guide = "none") +
    geom_text_repel(aes(label = .data[[label_col]]),
                    size = (base_font_size - 1.5)/.pt,
                    segment.size = 0.2, max.overlaps = Inf,
                    min.segment.length = 0, box.padding = 0.3) +
    labs(x = "RNA delta (FAP-KO minus WT)",
         y = "Proteomics delta",
         title = title_str,
         subtitle = sprintf("Spearman %.2f | Pearson %.2f", spear, pear)) +
    theme_nc()
}

p_prog_delta <- plot_delta_concordance(progeny_wide, "pathway",
                                       "PROGENy 14 pathways delta concordance") +
  labs(tag = "A")

## forest-style ordered bar for each modality
plot_forest <- function(delta_df, name_col, title_str, top_k = 14){
  d <- delta_df %>% filter(!is.na(delta_fap_minus_wt))
  d <- d %>%
    mutate(name = .data[[name_col]],
           sig_tag = case_when(wilcox_p < 0.05 ~ "*",
                               wilcox_p < 0.1  ~ ".",
                               TRUE ~ "")) %>%
    arrange(dataset, delta_fap_minus_wt) %>%
    group_by(dataset) %>%
    mutate(name_ordered = factor(name, levels = name))

  ggplot(d, aes(delta_fap_minus_wt, name_ordered,
                fill = delta_fap_minus_wt)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_col(colour = "black", linewidth = 0.2) +
    geom_text(aes(label = sig_tag), hjust = -0.3,
              size = (base_font_size + 0.5)/.pt) +
    facet_wrap(~ dataset, nrow = 1, scales = "free_y") +
    scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                         midpoint = 0, guide = "none",
                         limits = c(-max(abs(d$delta_fap_minus_wt), na.rm = TRUE),
                                     max(abs(d$delta_fap_minus_wt), na.rm = TRUE))) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(x = "Activity delta (FAP-KO minus WT)", y = NULL,
         title = title_str,
         subtitle = "* P<0.05, . P<0.1 (Wilcoxon)") +
    theme_nc()
}

p_prog_forest <- plot_forest(progeny_delta %>% filter(!is.na(delta_fap_minus_wt)),
                             "pathway", "PROGENy activities by modality") +
  labs(tag = "B")

## TF top
if (!is.null(tf_delta)){
  tf_top <- tf_delta %>% filter(dataset == "RNA", !is.na(wilcox_p)) %>%
    arrange(wilcox_p) %>% slice_head(n = 20) %>% pull(tf)
  tf_focus <- tf_delta %>% filter(tf %in% tf_top)
  p_tf_forest <- plot_forest(tf_focus, "tf",
                             "Top TF activities (by RNA Wilcox p)") +
    labs(tag = "C")
} else {
  p_tf_forest <- ggplot() + theme_void() +
    annotate("text", x=1, y=1, label = "CollecTRI unavailable")
}

## Hallmark
if (!is.null(hm_delta)){
  hm_top <- hm_delta %>% filter(dataset == "RNA", !is.na(wilcox_p)) %>%
    arrange(wilcox_p) %>% slice_head(n = 20) %>% pull(hallmark)
  hm_focus <- hm_delta %>% filter(hallmark %in% hm_top) %>%
    mutate(hallmark = sub("^HALLMARK_", "", hallmark))
  p_hm_forest <- plot_forest(hm_focus, "hallmark",
                             "Top Hallmark activities (by RNA Wilcox p)") +
    labs(tag = "D")
} else {
  p_hm_forest <- ggplot() + theme_void() +
    annotate("text", x=1, y=1, label = "Hallmark gmt not provided")
}

fig_integr <- (p_prog_delta | p_prog_forest) /
              (p_tf_forest  | p_hm_forest) +
  plot_layout(heights = c(1, 1.2)) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(out_root, "decoupler_crossomics_overview.pdf")
png_path <- file.path(out_root, "decoupler_crossomics_overview.png")
ggsave(pdf_path, fig_integr, width = 240, height = 230, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(png_path, fig_integr, width = 240, height = 230, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("\ndecoupleR cross-omics outputs written to:\n  ", out_root, "\n")
