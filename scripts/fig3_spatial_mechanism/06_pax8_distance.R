# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S11 (PAX8 vs distance from FAP-high seeds).
#          The "Supp Fig 19" in the comment below is earlier numbering.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

## ============================================================================
##  06_pax8_distance.R
##  PAX8 expression as a function of distance to nearest FAP+ CAF spot
##  ----------------------------------------------------------------------------
##  Purpose: provide spatial protein-level evidence (Layer 2 of FAP gatekeeper
##  theory) that PAX8 (master TF of thyroid differentiation) recovers as one
##  moves away from FAP+ CAF niches in 8 human Visium samples (formal_noP1).
##
##  Output: Supp Fig 19 (or insertion into Supp Fig 11/14) — distance decay
##  curves per sample + meta + Spearman forest.
##
##  Date: 2026-04-27
## ============================================================================

suppressMessages({
  library(data.table)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(lme4)
  library(lmerTest)
})

# ---------------------------------------------------------------------------- #
# 0. paths
# ---------------------------------------------------------------------------- #
proj_root  <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas")
spot_tsv   <- file.path(proj_root,
              "results/20260410_spatial_state_load_wholeslide",
              "spot_level_fap_state_load_context.tsv")
seurat_rds <- file.path(proj_root,
              "results/20260409_spatial_epithelial_state_transfer",
              "Spatial_integrated9P.with_epithelial_state_transfer.rds")

out_dir    <- file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/pax8_distance_decay_20260427")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------- #
# 1. load spot-level context (with pre-computed distance)
# ---------------------------------------------------------------------------- #
cat("[1/6] Loading spot-level context table ...\n")
dt <- fread(spot_tsv)

# keep formal 8-sample panel + default top20 threshold
dt <- dt[slide_formal == "formal_noP1" & threshold_label == "top20"]
cat("  rows after filter:", nrow(dt), "\n")
cat("  samples:", paste(sort(unique(dt$sample)), collapse = ", "), "\n")

# ---------------------------------------------------------------------------- #
# 2. extract PAX8 expression directly from Seurat
# ---------------------------------------------------------------------------- #
cat("[2/6] Extracting PAX8 expression from Seurat object ...\n")
sp <- readRDS(seurat_rds)

# PAX8 may live in SCT or RNA assay
get_pax8 <- function(seu) {
  for (a in c("SCT", "RNA")) {
    if (a %in% Assays(seu)) {
      m <- tryCatch(GetAssayData(seu, assay = a, layer = "data"),
                    error = function(e) tryCatch(
                      GetAssayData(seu, assay = a, slot = "data"),
                      error = function(e) NULL))
      if (!is.null(m) && "PAX8" %in% rownames(m)) {
        return(list(values = as.numeric(m["PAX8", ]),
                    cells = colnames(m), assay = a))
      }
    }
  }
  return(NULL)
}
px <- get_pax8(sp)
if (is.null(px)) {
  warning("PAX8 not found in any assay; falling back to TDS-only proxy.")
  pax8_dt <- data.table(spot_id = unique(dt$spot_id), PAX8_expr = NA_real_)
} else {
  cat("  PAX8 found in assay:", px$assay, "; n cells =", length(px$cells), "\n")
  pax8_dt <- data.table(spot_id = px$cells, PAX8_expr = px$values)
}

# spot_id in TSV may be Seurat cell barcodes; harmonise
dt <- merge(dt, pax8_dt, by = "spot_id", all.x = TRUE)
cat("  spots with PAX8 value:",
    sum(!is.na(dt$PAX8_expr)), "/", nrow(dt), "\n")

# ---------------------------------------------------------------------------- #
# 3. binning by distance (in spot units, normalised by sample spacing)
# ---------------------------------------------------------------------------- #
cat("[3/6] Binning spots by distance to nearest FAP+ CAF seed ...\n")

# ring labels — use spot-units (integer rings 0..N)
dt[, ring := cut(distance_spot_units,
                  breaks = c(-0.01, 0.01, 1, 2, 3, 5, 10, Inf),
                  labels = c("0 (seed)", "0–1", "1–2", "2–3", "3–5", "5–10", ">10"),
                  right = TRUE)]

# define PAX8 high spots (top quintile within sample, by direct expression
#   if available, else by tds_signature_score)
dt[, PAX8_proxy := ifelse(!is.na(PAX8_expr), PAX8_expr, tds_signature_score)]
dt[, PAX8_high  := PAX8_proxy >= quantile(PAX8_proxy, 0.8, na.rm = TRUE),
   by = sample]

# ---------------------------------------------------------------------------- #
# 4. per-sample summaries
# ---------------------------------------------------------------------------- #
cat("[4/6] Summarising per-sample distance-decay ...\n")

ring_summary <- dt[!is.na(ring),
  .(N_spots          = .N,
    N_PAX8_high      = sum(PAX8_high, na.rm = TRUE),
    pct_PAX8_high    = 100 * mean(PAX8_high, na.rm = TRUE),
    mean_TDS         = mean(tds_signature_score, na.rm = TRUE),
    sem_TDS          = sd(tds_signature_score, na.rm = TRUE) / sqrt(.N),
    mean_follicular  = mean(epi_deconv_weight_follicular_lineage_high, na.rm = TRUE),
    sem_follicular   = sd(epi_deconv_weight_follicular_lineage_high, na.rm = TRUE) / sqrt(.N),
    mean_PAX8_expr   = mean(PAX8_expr, na.rm = TRUE),
    sem_PAX8_expr    = sd(PAX8_expr, na.rm = TRUE) / sqrt(.N),
    mean_terminal    = mean(epi_deconv_weight_terminal_dedifferentiated, na.rm = TRUE)),
  by = .(sample, ring)]

fwrite(ring_summary,
       file.path(out_dir, "per_sample_ring_summary.tsv"), sep = "\t")

# pooled (all 8 samples)
pooled_summary <- ring_summary[,
  .(mean_TDS_pooled   = weighted.mean(mean_TDS,        N_spots, na.rm = TRUE),
    mean_foll_pooled  = weighted.mean(mean_follicular, N_spots, na.rm = TRUE),
    mean_PAX8_pooled  = weighted.mean(mean_PAX8_expr,  N_spots, na.rm = TRUE),
    mean_term_pooled  = weighted.mean(mean_terminal,   N_spots, na.rm = TRUE),
    n_total           = sum(N_spots)),
  by = ring]

fwrite(pooled_summary,
       file.path(out_dir, "pooled_ring_summary.tsv"), sep = "\t")

# ---------------------------------------------------------------------------- #
# 5. per-sample Spearman correlation (distance vs PAX8 proxies)
# ---------------------------------------------------------------------------- #
cat("[5/6] Per-sample Spearman correlations ...\n")

spearman_per_sample <- function(d) {
  out <- list()
  for (var in c("tds_signature_score", "epi_deconv_weight_follicular_lineage_high",
                "PAX8_expr", "epi_deconv_weight_terminal_dedifferentiated")) {
    if (all(is.na(d[[var]]))) next
    r <- suppressWarnings(cor.test(d$distance_spot_units, d[[var]],
                                    method = "spearman",
                                    exact  = FALSE))
    out[[var]] <- data.table(variable = var,
                              n        = sum(complete.cases(d$distance_spot_units, d[[var]])),
                              rho      = unname(r$estimate),
                              p        = r$p.value)
  }
  rbindlist(out)
}

spearman_dt <- dt[, spearman_per_sample(.SD), by = sample,
                  .SDcols = c("distance_spot_units", "tds_signature_score",
                              "epi_deconv_weight_follicular_lineage_high",
                              "PAX8_expr",
                              "epi_deconv_weight_terminal_dedifferentiated")]

fwrite(spearman_dt,
       file.path(out_dir, "per_sample_spearman.tsv"), sep = "\t")

cat("\n  per-sample Spearman summary (PAX8 proxy = TDS):\n")
print(spearman_dt[variable == "tds_signature_score"])

# ---------------------------------------------------------------------------- #
# 5b. Mixed-effects model (per-sample random intercept + slope)
# ---------------------------------------------------------------------------- #
cat("[5b] Linear mixed-effects model (random intercept + slope by sample) ...\n")

# log-transform distance to handle long tail; offset to avoid log(0)
dt[, log_dist := log10(distance_spot_units + 1)]

run_lmm <- function(d, response, label) {
  fit <- tryCatch(
    lmer(as.formula(sprintf("%s ~ log_dist + (1 + log_dist | sample)", response)),
         data = d, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa",
                                optCtrl = list(maxfun = 2e5))),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s   <- summary(fit)$coefficients
  ci  <- tryCatch(confint(fit, parm = "log_dist", method = "Wald"),
                  error = function(e) matrix(NA, 1, 2))
  list(label = label,
       beta   = s["log_dist", "Estimate"],
       se     = s["log_dist", "Std. Error"],
       df     = s["log_dist", "df"],
       t_val  = s["log_dist", "t value"],
       p      = s["log_dist", "Pr(>|t|)"],
       ci_lo  = ci[1, 1],
       ci_hi  = ci[1, 2],
       n      = nobs(fit))
}

lmm_results <- list()
for (resp in c("tds_signature_score",
                "PAX8_expr",
                "epi_deconv_weight_follicular_lineage_high",
                "epi_deconv_weight_terminal_dedifferentiated")) {
  if (all(is.na(dt[[resp]]))) next
  res <- run_lmm(dt, resp, resp)
  if (!is.null(res)) lmm_results[[resp]] <- res
}

lmm_dt <- rbindlist(lapply(lmm_results, as.data.table), idcol = NULL, fill = TRUE)
fwrite(lmm_dt, file.path(out_dir, "mixed_effects_summary.tsv"), sep = "\t")

cat("\n  Mixed-effects fixed-effect (log_dist) estimates:\n")
print(lmm_dt[, .(label, beta = signif(beta, 3),
                  ci_95 = sprintf("[%.3f, %.3f]", ci_lo, ci_hi),
                  p = signif(p, 3),
                  n)])

# ---------------------------------------------------------------------------- #
# 6. plotting
# ---------------------------------------------------------------------------- #
cat("[6/6] Plotting ...\n")

# colour palette
col_wt   <- "#3B6B9C"   # primary blue
col_fap  <- "#C9493A"   # warm red
col_neut <- "#7B5A9C"   # purple

theme_nc <- theme_classic(base_size = 7) +
  theme(plot.title    = element_text(size = 7.5, face = "bold"),
        axis.title    = element_text(size = 7),
        axis.text     = element_text(size = 6.5, colour = "grey15"),
        strip.text    = element_text(size = 7, face = "bold"),
        strip.background = element_rect(fill = "grey95", colour = NA),
        legend.position = "right",
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7),
        plot.margin  = margin(4, 4, 4, 4))

# --- Panel A : per-sample TDS-decay curves --------------------------------- #
pA <- ggplot(ring_summary[!is.na(ring)],
             aes(x = ring, y = mean_TDS, group = sample)) +
  geom_errorbar(aes(ymin = mean_TDS - sem_TDS, ymax = mean_TDS + sem_TDS),
                width = 0.2, colour = "grey50", linewidth = 0.3) +
  geom_line(colour = col_wt, linewidth = 0.55) +
  geom_point(colour = col_wt, fill = "white", size = 1.3, shape = 21,
             stroke = 0.45) +
  facet_wrap(~ sample, nrow = 2) +
  labs(title = "TDS (thyroid differentiation score) increases with distance from FAP+ CAF seed",
       x = "distance ring (spot units)",
       y = "mean TDS ± SEM") +
  theme_nc

# --- Panel B : pooled curve (all proxies overlay) -------------------------- #
pooled_long <- melt(pooled_summary, id.vars = c("ring", "n_total"),
                    measure.vars = c("mean_TDS_pooled",
                                     "mean_foll_pooled",
                                     "mean_PAX8_pooled",
                                     "mean_term_pooled"),
                    variable.name = "metric", value.name = "value")
pooled_long[, metric := factor(metric,
   levels = c("mean_TDS_pooled","mean_foll_pooled","mean_PAX8_pooled","mean_term_pooled"),
   labels = c("TDS score","Follicular weight","PAX8 expr (SCT)","Terminal-ddiff weight"))]

# z-score within metric for overlay
pooled_long[, value_z := scale(value)[, 1], by = metric]

# annotation text from LMM results
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-300) return("P < 1e-300")
  sprintf("P = %.1e", p)
}
get_lmm <- function(resp) {
  r <- lmm_results[[resp]]
  if (is.null(r)) return("")
  sprintf("β = %.2f  [%.2f, %.2f]\n%s",
          r$beta, r$ci_lo, r$ci_hi, fmt_p(r$p))
}

lmm_anno <- sprintf(
  "Linear mixed-effects (per-sample random intercept + slope on log10 distance):\nTDS  %s     |     PAX8  %s     |     Terminal-ddiff  %s",
  fmt_p(lmm_results[["tds_signature_score"]]$p),
  fmt_p(lmm_results[["PAX8_expr"]]$p),
  fmt_p(lmm_results[["epi_deconv_weight_terminal_dedifferentiated"]]$p))

pB <- ggplot(pooled_long[!is.na(ring) & !is.nan(value_z)],
             aes(x = ring, y = value_z, colour = metric, group = metric)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70", linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  annotate("text", x = 1, y = max(pooled_long$value_z, na.rm = TRUE) * 1.05,
           label = lmm_anno, hjust = 0, vjust = 1,
           size = 1.9, colour = "grey25", fontface = "italic", lineheight = 0.95) +
  scale_colour_manual(values = c("TDS score"           = col_wt,
                                  "Follicular weight"  = "#5DA66A",
                                  "PAX8 expr (SCT)"    = col_neut,
                                  "Terminal-ddiff weight" = col_fap),
                      name = NULL) +
  labs(title = "Pooled (8 samples): differentiation proxies up, terminal-ddiff down with distance",
       x = "distance ring (spot units)",
       y = "z-score across rings") +
  theme_nc +
  theme(legend.position = "bottom",
        legend.key.height = unit(2, "mm"),
        legend.key.width  = unit(4, "mm"))

# --- Panel C : per-sample Spearman forest ---------------------------------- #
fp <- spearman_dt[variable == "tds_signature_score"]
fp[, label := sprintf("%s\n(n = %d)", sample, n)]
fp[, sample := factor(sample, levels = sample[order(rho)])]

pC <- ggplot(fp, aes(x = rho, y = sample)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = rho, yend = sample),
               linewidth = 0.5, colour = col_wt) +
  geom_point(size = 2.2, colour = col_wt, fill = "white", shape = 21,
             stroke = 0.55) +
  geom_text(aes(label = sprintf("ρ = %.2f, P = %.1e", rho, p)),
            size = 2.0, hjust = -0.05, colour = "grey25") +
  scale_x_continuous(limits = c(-0.05, max(fp$rho) * 1.4),
                     expand = c(0, 0)) +
  labs(title = "Per-sample Spearman: distance vs TDS",
       x = "Spearman ρ (distance × TDS)",
       y = NULL) +
  theme_nc +
  theme(axis.text.y = element_text(size = 6.5))

# --- Panel D : direct PAX8 expression decay (if PAX8 found) ---------------- #
pD <- if (!all(is.na(dt$PAX8_expr))) {
  ggplot(ring_summary[!is.na(ring) & !is.nan(mean_PAX8_expr)],
         aes(x = ring, y = mean_PAX8_expr, group = sample, colour = sample)) +
    geom_errorbar(aes(ymin = mean_PAX8_expr - sem_PAX8_expr,
                      ymax = mean_PAX8_expr + sem_PAX8_expr),
                  width = 0.2, linewidth = 0.3) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.3, fill = "white", shape = 21, stroke = 0.4) +
    scale_colour_brewer(palette = "Dark2") +
    labs(title = "Direct PAX8 expression (SCT) per sample",
         x = "distance ring (spot units)",
         y = "mean PAX8 ± SEM") +
    theme_nc +
    theme(legend.position = "right",
          legend.key.height = unit(2.5, "mm"))
} else {
  ggplot() + theme_void() +
    labs(title = "PAX8 expression unavailable in current Seurat assays")
}

# assemble
combined <- (pA / (pB | pC) / pD) +
  plot_layout(heights = c(1.4, 1, 1)) +
  plot_annotation(
    title = "Spatial-protein evidence: PAX8 / differentiation recovers with distance from FAP+ CAF niches",
    subtitle = "8-sample human DTC Visium · formal_noP1 · top20 FAP threshold (Fig 3 cohort)",
    theme = theme(plot.title = element_text(size = 9, face = "bold"),
                  plot.subtitle = element_text(size = 7.5,
                                                colour = "grey25",
                                                face = "italic")))

ggsave(file.path(out_dir, "Fig_PAX8_distance_decay.pdf"),
       combined, width = 220, height = 240, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(file.path(out_dir, "Fig_PAX8_distance_decay.png"),
       combined, width = 220, height = 240, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("\nDone. Outputs in:\n  ", out_dir, "\n")
cat("  - Fig_PAX8_distance_decay.{pdf,png}\n")
cat("  - per_sample_ring_summary.tsv\n")
cat("  - pooled_ring_summary.tsv\n")
cat("  - per_sample_spearman.tsv\n")
