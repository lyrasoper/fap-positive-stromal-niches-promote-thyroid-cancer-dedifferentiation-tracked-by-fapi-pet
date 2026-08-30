#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S16 - Factor 1 survival analysis in TCGA-THCA.
## Step 3: Kaplan-Meier / Cox regression of Factor1 signature score
## on TCGA-THCA with Xena-curated survival (OS, DSS, PFI).
## Uses Factor1 weighted score pre-computed in analyze_factor1_tcga_projection.R.


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(survival); library(survminer)
})

proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
bulk_out <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
score_path <- file.path(bulk_out, "factor1_tcga_projection/factor1_tcga_merged_annot.tsv")
surv_path  <- file.path(proj, "reference/tcga_thca/THCA_survival.txt")
cli_path   <- file.path(proj, "reference/tcga_thca/THCA_clinicalMatrix.tsv")
out_dir    <- file.path(bulk_out, "factor1_tcga_survival")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_wt  <- "#3B6B9C"; col_fap <- "#C9493A"; col_purple <- "#4B1E70"
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

## ---- load data -----------------------------------------------------------
scored <- read_tsv(score_path, show_col_types = FALSE)
surv <- read_tsv(surv_path, show_col_types = FALSE)
cli  <- read_tsv(cli_path, show_col_types = FALSE)

# keep only primary-tumor survival rows & short sample IDs
surv <- surv %>%
  mutate(sample = substr(sample, 1, 15))
cli  <- cli %>%
  mutate(sample = substr(sampleID, 1, 15))

cli_sel <- cli %>%
  transmute(sample,
            age = as.numeric(age_at_initial_pathologic_diagnosis),
            gender,
            stage = pathologic_stage,
            histology = histological_type,
            ajcc_T = pathologic_T,
            ajcc_N = pathologic_N,
            ajcc_M = pathologic_M) %>%
  distinct(sample, .keep_all = TRUE)

# assemble master
master <- scored %>%
  dplyr::select(sample, score_weighted, score_ulm,
                brs_group, brs_label_f, histological_type,
                brs_raw, tds, erk_score) %>%
  left_join(surv %>% dplyr::select(sample, OS, OS.time, DSS, DSS.time,
                                    PFI, PFI.time),
            by = "sample") %>%
  left_join(cli_sel %>% dplyr::select(-histology), by = "sample") %>%
  filter(!is.na(OS.time))
cat(sprintf("Samples with survival & score: %d\n", nrow(master)))
cat(sprintf("OS events: %d | PFI events: %d | DSS events: %d\n",
            sum(master$OS, na.rm = TRUE),
            sum(master$PFI, na.rm = TRUE),
            sum(master$DSS, na.rm = TRUE)))

# stratification: tertile by score_weighted within cohort
master <- master %>%
  mutate(score_tertile = cut(score_weighted,
                             breaks = quantile(score_weighted, c(0, 1/3, 2/3, 1),
                                               na.rm = TRUE),
                             include.lowest = TRUE,
                             labels = c("low", "mid", "high")),
         score_high = score_weighted >= quantile(score_weighted, 2/3, na.rm = TRUE),
         score_binary = factor(ifelse(score_high, "top third", "bottom 2/3"),
                               levels = c("bottom 2/3", "top third")))

write_tsv(master, file.path(out_dir, "factor1_tcga_survival_master.tsv"))

## ---- KM & Cox ------------------------------------------------------------
run_surv <- function(event_col, time_col, label){
  df <- master %>%
    filter(!is.na(.data[[event_col]]), !is.na(.data[[time_col]])) %>%
    mutate(.time = .data[[time_col]], .event = .data[[event_col]])
  if (sum(df$.event) < 5){
    cat(sprintf("%s: too few events (%d), skip\n",
                label, sum(df$.event)))
    return(NULL)
  }
  # KM by tertile
  sfit <- survfit(Surv(.time, .event) ~ score_tertile, data = df)
  diff <- survdiff(Surv(.time, .event) ~ score_tertile, data = df)
  lr_p <- 1 - pchisq(diff$chisq, df = length(diff$n) - 1)

  # Cox (continuous score)
  cox_uni <- coxph(Surv(.time, .event) ~ score_weighted, data = df)
  cox_uni_s <- summary(cox_uni)

  # Multivariable (age + stage)
  df2 <- df %>% mutate(stage_simple =
                         ifelse(grepl("IV", stage), "late", "early"))
  cox_mv <- try(coxph(Surv(.time, .event) ~
                        score_weighted + age + stage_simple, data = df2),
                silent = TRUE)
  cox_mv_s <- if (!inherits(cox_mv, "try-error")) summary(cox_mv) else NULL

  list(label = label, fit = sfit, diff = diff, logrank_p = lr_p,
       cox_uni = cox_uni_s, cox_mv = cox_mv_s, df = df)
}

surv_os  <- run_surv("OS",  "OS.time",  "OS")
surv_dss <- run_surv("DSS", "DSS.time", "DSS")
surv_pfi <- run_surv("PFI", "PFI.time", "PFI")

## ---- KM plots using survminer::ggsurvplot ------------------------------
make_km <- function(res, outcome_label){
  if (is.null(res)) return(NULL)
  cox_uni <- res$cox_uni
  hr  <- cox_uni$conf.int["score_weighted", "exp(coef)"]
  lci <- cox_uni$conf.int["score_weighted", "lower .95"]
  uci <- cox_uni$conf.int["score_weighted", "upper .95"]
  p_cox <- cox_uni$coefficients["score_weighted", "Pr(>|z|)"]

  g <- ggsurvplot(
    res$fit, data = res$df,
    risk.table = TRUE, conf.int = FALSE,
    palette = c(low = col_wt, mid = "grey60", high = col_fap),
    legend.labs = paste(levels(res$df$score_tertile), "Factor1"),
    legend.title = "",
    pval = sprintf("Log-rank %s\nCox HR %.2f (%.2f-%.2f), %s",
                   p_label(res$logrank_p), hr, lci, uci, p_label(p_cox)),
    pval.size = 2.5,
    risk.table.height = 0.28,
    ggtheme = theme_nc() + theme(legend.position = "top"),
    font.main = 7, font.x = 7, font.y = 7,
    font.tickslab = 6, font.legend = 6,
    tables.theme = theme_nc() + theme(axis.line = element_blank(),
                                       axis.ticks = element_blank())
  )
  g$plot <- g$plot +
    labs(x = sprintf("%s time (days)", outcome_label),
         y = sprintf("%s probability", outcome_label),
         title = sprintf("TCGA-THCA %s by Factor1 tertile", outcome_label))
  g
}

km_os  <- make_km(surv_os,  "OS")
km_dss <- make_km(surv_dss, "DSS")
km_pfi <- make_km(surv_pfi, "PFI")

## save individual KM plots (survminer objects include separate pieces)
save_km <- function(km, fname){
  if (is.null(km)) return(invisible(NULL))
  g <- cowplot::plot_grid(km$plot, km$table, ncol = 1,
                          rel_heights = c(3, 1.2))
  ggsave(file.path(out_dir, fname), g,
         width = 100, height = 90, units = "mm",
         device = grDevices::pdf, useDingbats = FALSE)
  ggsave(sub("\\.pdf$", ".png", file.path(out_dir, fname)), g,
         width = 100, height = 90, units = "mm",
         dpi = 600, device = ragg::agg_png)
}
save_km(km_os,  "KM_OS.pdf")
save_km(km_dss, "KM_DSS.pdf")
save_km(km_pfi, "KM_PFI.pdf")

## ---- summary stats table -----------------------------------------------
make_row <- function(res, outcome){
  if (is.null(res)) return(NULL)
  c <- res$cox_uni
  cm <- res$cox_mv
  out <- tibble(
    outcome = outcome,
    n = nrow(res$df),
    events = sum(res$df$.event, na.rm = TRUE),
    logrank_p = res$logrank_p,
    cox_HR = c$conf.int["score_weighted", "exp(coef)"],
    cox_HR_lo = c$conf.int["score_weighted", "lower .95"],
    cox_HR_hi = c$conf.int["score_weighted", "upper .95"],
    cox_p = c$coefficients["score_weighted", "Pr(>|z|)"]
  )
  if (!is.null(cm)){
    out <- out %>% mutate(
      cox_mv_HR_score = cm$conf.int["score_weighted", "exp(coef)"],
      cox_mv_p_score  = cm$coefficients["score_weighted", "Pr(>|z|)"]
    )
  }
  out
}
stats_tbl <- bind_rows(
  make_row(surv_os,  "OS"),
  make_row(surv_dss, "DSS"),
  make_row(surv_pfi, "PFI")
)
write_tsv(stats_tbl, file.path(out_dir, "factor1_tcga_survival_stats.tsv"))
cat("\n==== survival summary ====\n")
print(stats_tbl)

## ---- Forest plot of Cox HR (univariate + mv) ----------------------------
build_forest_df <- function(res, outcome){
  if (is.null(res)) return(NULL)
  bind_rows(
    tibble(outcome = outcome, model = "univariate",
           HR = res$cox_uni$conf.int["score_weighted", "exp(coef)"],
           lo = res$cox_uni$conf.int["score_weighted", "lower .95"],
           hi = res$cox_uni$conf.int["score_weighted", "upper .95"],
           p  = res$cox_uni$coefficients["score_weighted", "Pr(>|z|)"]),
    if (!is.null(res$cox_mv))
      tibble(outcome = outcome, model = "age + stage adj",
             HR = res$cox_mv$conf.int["score_weighted", "exp(coef)"],
             lo = res$cox_mv$conf.int["score_weighted", "lower .95"],
             hi = res$cox_mv$conf.int["score_weighted", "upper .95"],
             p  = res$cox_mv$coefficients["score_weighted", "Pr(>|z|)"])
    else NULL
  )
}
forest_df <- bind_rows(
  build_forest_df(surv_os,  "OS"),
  build_forest_df(surv_dss, "DSS"),
  build_forest_df(surv_pfi, "PFI")
) %>% mutate(row_lbl = paste0(outcome, " | ", model),
             p_str = vapply(p, p_label, character(1)),
             panel_lbl = sprintf("HR %.2f (%.2f-%.2f), %s",
                                 HR, lo, hi, p_str))

p_forest <- ggplot(forest_df, aes(HR, row_lbl)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             linewidth = 0.3, colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4,
                 colour = "grey25") +
  geom_point(aes(fill = log2(HR)), shape = 21, size = 2.2,
             stroke = 0.3, colour = "black") +
  geom_text(aes(x = hi, label = panel_lbl),
            hjust = -0.1, size = (base_font_size - 1.5)/.pt) +
  scale_fill_gradient2(low = col_wt, mid = "white", high = col_fap,
                       midpoint = 0, guide = "none") +
  scale_x_log10() +
  labs(x = "Cox HR per 1 unit Factor1 score (log scale)",
       y = NULL,
       title = "Factor1 score - survival in TCGA-THCA",
       subtitle = "PFI most informative (more events than OS/DSS)") +
  theme_nc()

ggsave(file.path(out_dir, "cox_forest.pdf"), p_forest,
       width = 140, height = 60, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(file.path(out_dir, "cox_forest.png"), p_forest,
       width = 140, height = 60, units = "mm",
       dpi = 600, device = ragg::agg_png)

## ---- combined panel: forest + PFI KM ----------------------------------
if (!is.null(km_pfi)){
  g_pfi <- cowplot::plot_grid(km_pfi$plot, km_pfi$table, ncol = 1,
                              rel_heights = c(3, 1.2))
  fig_comb <- cowplot::plot_grid(p_forest, g_pfi, nrow = 1,
                                 rel_widths = c(1.25, 1))
  ggsave(file.path(out_dir, "factor1_tcga_survival_overview.pdf"),
         fig_comb, width = 220, height = 95, units = "mm",
         device = grDevices::pdf, useDingbats = FALSE)
  ggsave(file.path(out_dir, "factor1_tcga_survival_overview.png"),
         fig_comb, width = 220, height = 95, units = "mm",
         dpi = 600, device = ragg::agg_png)
}

cat("\nSurvival outputs written to:\n  ", out_dir, "\n")
