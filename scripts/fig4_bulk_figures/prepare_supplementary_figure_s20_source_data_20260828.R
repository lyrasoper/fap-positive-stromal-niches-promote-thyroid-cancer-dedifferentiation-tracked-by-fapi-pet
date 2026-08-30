#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S20 - panel-ready source-data bundle.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

# Prepare a panel-ready, immutable source-data bundle for Supplementary Fig. S20.
# The legacy analysis builder is evaluated in an isolated environment with its
# export functions disabled, so the upstream project figures are not overwritten.

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

bulk_project <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
builder <- file.path(project_root, "scripts", "fig4_bulk_omics",
                     "build_figS17_progeny_cross_cohort.R")
if (!file.exists(builder)) {
  stop("Cannot find the cross-cohort PROGENy builder at ", builder)
}

required <- c("dplyr", "tidyr", "readr", "tibble")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

output_dir <- file.path(project_root, "outputs", "suppfig20_beautified_R_20260828")
source_dir <- file.path(output_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

# Prevent the sourced legacy builder from rewriting any upstream TSV or figure.
analysis_env <- new.env(parent = globalenv())
analysis_env$write_tsv <- function(...) invisible(NULL)
analysis_env$ggsave <- function(...) invisible(NULL)

message("Recomputing the established PROGENy objects without legacy exports ...")
builder_log <- capture.output(
  suppressWarnings(suppressMessages(source(builder, local = analysis_env))),
  type = "output"
)
writeLines(builder_log, file.path(output_dir, "source_data_preparation.log"))

needed_objects <- c("cross_tbl", "factor1_cor", "ladder_long", "kt_stats")
missing_objects <- needed_objects[!vapply(needed_objects, exists, logical(1), envir = analysis_env)]
if (length(missing_objects)) {
  stop("The legacy analysis did not create: ", paste(missing_objects, collapse = ", "))
}

pathway_levels <- c(
  "MAPK", "PI3K", "EGFR", "VEGF", "JAK-STAT", "NFkB", "TNFa",
  "TGFb", "WNT", "Hypoxia", "p53", "Trail", "Estrogen", "Androgen"
)

# Panel a: raw pathway-score differences and the explicit within-cohort scaling
# used only for the heatmap colour. Cell text remains on the raw-score scale.
panel_a <- analysis_env$cross_tbl %>%
  mutate(pathway = as.character(pathway)) %>%
  pivot_longer(
    cols = -pathway,
    names_to = c("cohort", ".value"),
    names_pattern = "(mouse|tcga|landa)_(delta|p)"
  ) %>%
  mutate(
    cohort = recode(
      cohort,
      mouse = "Mouse: Fap-deficient minus WT",
      tcga = "TCGA: BRAF V600E minus RAS driver",
      landa = "GSE76039: ATC minus PDTC"
    )
  ) %>%
  group_by(cohort) %>%
  mutate(scaled_delta = delta / max(abs(delta), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    pathway = factor(pathway, levels = pathway_levels),
    n_group_1 = case_when(
      cohort == "Mouse: Fap-deficient minus WT" ~ 4L,
      cohort == "TCGA: BRAF V600E minus RAS driver" ~ 213L,
      cohort == "GSE76039: ATC minus PDTC" ~ 20L
    ),
    n_group_2 = case_when(
      cohort == "Mouse: Fap-deficient minus WT" ~ 4L,
      cohort == "TCGA: BRAF V600E minus RAS driver" ~ 48L,
      cohort == "GSE76039: ATC minus PDTC" ~ 17L
    ),
    significance_symbol = case_when(p < 0.05 ~ "*", p < 0.10 ~ "dot", TRUE ~ "")
  ) %>%
  arrange(pathway, cohort)

# Panel b: effect sizes and both nominal and BH-adjusted P values.
panel_b <- analysis_env$factor1_cor %>%
  transmute(
    pathway = factor(as.character(pathway), levels = pathway_levels),
    rho,
    p_nominal = p,
    p_bh_14_pathways = p.adjust(p, method = "BH")
  ) %>%
  arrange(pathway)

# Panel c: per-tumour scores after the established within-cohort z transform.
panel_c_scores <- analysis_env$ladder_long %>%
  transmute(
    sample,
    ordered_group = as.character(ladder),
    cohort = if_else(grepl("^TCGA", ordered_group), "TCGA", "GSE76039"),
    pathway = as.character(pathway),
    progeny_activity_within_cohort_z = z
  ) %>%
  mutate(
    ordered_group = factor(
      ordered_group,
      levels = c("TCGA RAS", "TCGA other", "TCGA BRAF", "Landa PDTC", "Landa ATC")
    ),
    pathway = factor(pathway, levels = c("MAPK", "WNT", "TGFb", "JAK-STAT"))
  ) %>%
  arrange(pathway, ordered_group, sample)

panel_c_trends <- analysis_env$kt_stats %>%
  transmute(
    pathway = as.character(pathway),
    kendall_tau = as.numeric(tau),
    p_nominal = p
  ) %>%
  arrange(match(pathway, c("MAPK", "WNT", "TGFb", "JAK-STAT")))

pair_specs <- tribble(
  ~comparison, ~group_1, ~group_2,
  "TCGA RAS vs BRAF V600E", "TCGA RAS", "TCGA BRAF",
  "GSE76039 PDTC vs ATC", "Landa PDTC", "Landa ATC"
)

panel_c_pairwise <- tidyr::crossing(
  pathway = c("MAPK", "WNT", "TGFb", "JAK-STAT"),
  pair_specs
) %>%
  rowwise() %>%
  mutate(
    n_group_1 = sum(
      panel_c_scores$pathway == pathway &
        panel_c_scores$ordered_group == group_1 &
        !is.na(panel_c_scores$progeny_activity_within_cohort_z)
    ),
    n_group_2 = sum(
      panel_c_scores$pathway == pathway &
        panel_c_scores$ordered_group == group_2 &
        !is.na(panel_c_scores$progeny_activity_within_cohort_z)
    ),
    p_nominal = {
      x <- panel_c_scores$progeny_activity_within_cohort_z[
        panel_c_scores$pathway == pathway & panel_c_scores$ordered_group == group_1
      ]
      y <- panel_c_scores$progeny_activity_within_cohort_z[
        panel_c_scores$pathway == pathway & panel_c_scores$ordered_group == group_2
      ]
      stats::wilcox.test(x, y, alternative = "two.sided", exact = FALSE)$p.value
    }
  ) %>%
  ungroup() %>%
  group_by(comparison) %>%
  mutate(p_bh_4_pathways_within_comparison = p.adjust(p_nominal, method = "BH")) %>%
  ungroup()

# Panel d: the ten prespecified displayed sets, retaining layer-specific FDR.
gsea_rna_path <- file.path(
  bulk_project, "outputs", "wt_vs_fapko_bulk_omics_20260407",
  "rna", "rna_all_gsea_pathways.tsv"
)
gsea_protein_path <- file.path(
  bulk_project, "outputs", "wt_vs_fapko_bulk_omics_20260407",
  "proteomics", "proteomics_all_gsea_pathways.tsv"
)
if (!file.exists(gsea_rna_path) || !file.exists(gsea_protein_path)) {
  stop("Cannot find the RNA/protein GSEA result tables.")
}

curated <- tribble(
  ~pathway_key, ~display_label, ~pathway_class,
  "ecm receptor interaction", "ECM-receptor interaction", "ECM / adhesion",
  "focal adhesion", "Focal adhesion", "ECM / adhesion",
  "collagen biosynthesis and modifying enzymes", "Collagen biosynthesis/modifying", "ECM / adhesion",
  "collagen formation", "Collagen formation", "ECM / adhesion",
  "extracellular matrix organization", "ECM organization", "ECM / adhesion",
  "integrin cell surface interactions", "Integrin cell-surface interactions", "ECM / adhesion",
  "thyroid cancer", "Thyroid cancer", "Thyroid / lineage",
  "thyroid hormone synthesis", "Thyroid hormone synthesis", "Thyroid / lineage",
  "autoimmune thyroid disease", "Autoimmune thyroid disease", "Thyroid / lineage",
  "mapk signaling pathway", "MAPK signaling", "MAPK"
)

rna_gsea <- read_tsv(gsea_rna_path, show_col_types = FALSE) %>%
  transmute(
    pathway_key = tolower(pathway), layer = "RNA", nes = NES,
    p_nominal = p_value, fdr_bh = fdr, direction
  ) %>%
  distinct(pathway_key, .keep_all = TRUE)

protein_gsea <- read_tsv(gsea_protein_path, show_col_types = FALSE) %>%
  transmute(
    pathway_key = tolower(pathway), layer = "Protein", nes = NES,
    p_nominal = p_value, fdr_bh = fdr, direction
  ) %>%
  distinct(pathway_key, .keep_all = TRUE)

panel_d <- curated %>%
  left_join(bind_rows(rna_gsea, protein_gsea), by = "pathway_key") %>%
  filter(!is.na(nes)) %>%
  mutate(fdr_lt_0_05 = fdr_bh < 0.05) %>%
  arrange(
    factor(pathway_class, levels = c("ECM / adhesion", "Thyroid / lineage", "MAPK")),
    display_label,
    factor(layer, levels = c("RNA", "Protein"))
  )

write_tsv(panel_a, file.path(source_dir, "Supplementary_Figure_S20_panel_A_cross_cohort.tsv"))
write_tsv(panel_b, file.path(source_dir, "Supplementary_Figure_S20_panel_B_factor1_correlations.tsv"))
write_tsv(panel_c_scores, file.path(source_dir, "Supplementary_Figure_S20_panel_C_tumour_scores.tsv"))
write_tsv(panel_c_trends, file.path(source_dir, "Supplementary_Figure_S20_panel_C_kendall_trends.tsv"))
write_tsv(panel_c_pairwise, file.path(source_dir, "Supplementary_Figure_S20_panel_C_wilcoxon.tsv"))
write_tsv(panel_d, file.path(source_dir, "Supplementary_Figure_S20_panel_D_GSEA.tsv"))

manifest <- tribble(
  ~panel, ~file, ~content,
  "A", "Supplementary_Figure_S20_panel_A_cross_cohort.tsv", "Raw pathway-score deltas, nominal Wilcoxon P values, and explicit within-cohort colour scaling",
  "B", "Supplementary_Figure_S20_panel_B_factor1_correlations.tsv", "Spearman effect sizes with nominal and 14-pathway BH-adjusted P values",
  "C", "Supplementary_Figure_S20_panel_C_tumour_scores.tsv", "Per-tumour, within-cohort z-standardized PROGENy activities",
  "C", "Supplementary_Figure_S20_panel_C_kendall_trends.tsv", "Kendall trend statistics across the displayed group order",
  "C", "Supplementary_Figure_S20_panel_C_wilcoxon.tsv", "Displayed within-cohort two-sided Wilcoxon comparisons",
  "D", "Supplementary_Figure_S20_panel_D_GSEA.tsv", "Displayed RNA/protein GSEA results and layer-specific BH FDR values"
)
write_tsv(manifest, file.path(source_dir, "SOURCE_DATA_MANIFEST.tsv"))

message("Panel-ready source data written to: ", source_dir)
