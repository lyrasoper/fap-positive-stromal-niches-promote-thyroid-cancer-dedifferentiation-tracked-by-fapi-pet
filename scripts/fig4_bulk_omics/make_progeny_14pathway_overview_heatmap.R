#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S20 - 14-pathway PROGENy overview heatmap.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

script_path <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
}

outdir <- file.path(root, "outputs", "wt_vs_fapko_bulk_omics_20260407", "progeny", "integrated")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pathway_groups <- list(
  "hormone-like signaling" = c("Androgen", "Estrogen"),
  "growth/RTK signaling" = c("EGFR", "MAPK", "PI3K", "WNT"),
  "stromal/angiogenic signaling" = c("TGFb", "VEGF", "Hypoxia"),
  "inflammatory/stress signaling" = c("JAK-STAT", "NFkB", "TNFa", "Trail", "p53")
)

group_levels <- names(pathway_groups)
pathway_map <- bind_rows(lapply(seq_along(pathway_groups), function(i) {
  category <- group_levels[i]
  pathways <- unname(pathway_groups[[i]])
  tibble(
    category = category,
    pathway = pathways,
    pathway_order = seq_along(pathways)
  )
}))

read_stats <- function(path, dataset_name) {
  read.delim(path, check.names = FALSE) |>
    transmute(
      pathway,
      dataset = dataset_name,
      delta_z = as.numeric(delta_fap_minus_wt_z),
      p_value = as.numeric(p_value),
      q_value = as.numeric(q_value)
    )
}

stats_df <- bind_rows(
  read_stats(file.path(outdir, "..", "rna", "rna_progeny_pathway_stats.tsv"), "RNA"),
  read_stats(file.path(outdir, "..", "proteomics", "proteomics_progeny_pathway_stats.tsv"), "Proteomics")
) 

plot_df <- pathway_map |>
  tidyr::crossing(dataset = c("RNA", "Proteomics")) |>
  left_join(stats_df, by = c("pathway", "dataset")) |>
  mutate(
    category = factor(category, levels = group_levels),
    dataset = factor(dataset, levels = c("RNA", "Proteomics")),
    pathway_facet = paste(category, pathway, sep = "___"),
    sig_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.05 ~ "**",
      p_value < 0.10 ~ "*",
      TRUE ~ ""
    )
  ) |>
  mutate(
    pathway_facet = factor(
      pathway_facet,
      levels = rev(unique(paste(pathway_map$category, pathway_map$pathway, sep = "___")))
    )
  )

wide_df <- plot_df |>
  transmute(category, pathway = sub("^.*___", "", as.character(pathway_facet)), dataset, delta_z, p_value, q_value) |>
  pivot_wider(
    names_from = dataset,
    values_from = c(delta_z, p_value, q_value),
    names_vary = "slowest"
  ) |>
  arrange(category, desc(pathway))

write.table(
  wide_df,
  file.path(outdir, "07_progeny_14pathway_overview.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  pathway_map,
  file.path(outdir, "07_progeny_14pathway_categories.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p <- ggplot(plot_df, aes(x = dataset, y = pathway_facet, fill = delta_z)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sig_label), size = 3.4) +
  scale_fill_gradient2(
    low = "#2563EB",
    mid = "white",
    high = "#DC2626",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    title = "PROGENy 14-pathway overview",
    subtitle = "* p < 0.10, ** p < 0.05 (raw p); red = higher in FAP-deficient host",
    x = NULL,
    y = NULL,
    fill = "delta z\n(FAP-WT)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey40"),
    strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 0),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9),
    panel.spacing.y = unit(0.3, "lines")
  )

p_compact <- ggplot(plot_df, aes(x = dataset, y = pathway_facet, fill = delta_z)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = sig_label), size = 2.9) +
  scale_fill_gradient2(
    low = "#5E81AC",
    mid = "#FBFBFA",
    high = "#D08770",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    x = NULL,
    y = NULL,
    fill = "delta z"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey96", colour = "grey55", linewidth = 0.5),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 6.9),
    axis.text.x = element_text(face = "bold", size = 8.5),
    axis.text.y = element_text(face = "bold", size = 7.8),
    legend.title = element_text(size = 8.2, face = "bold"),
    legend.text = element_text(size = 7.6),
    legend.key.height = unit(0.42, "in"),
    legend.key.width = unit(0.16, "in"),
    panel.spacing.y = unit(0.16, "lines"),
    plot.margin = margin(3, 4, 3, 4)
  )

ggsave(
  file.path(outdir, "07_progeny_14pathway_overview_heatmap.png"),
  p,
  width = 6.2,
  height = 8.8,
  dpi = 300
)
ggsave(
  file.path(outdir, "07_progeny_14pathway_overview_heatmap.pdf"),
  p,
  width = 6.2,
  height = 8.8
)

ggsave(
  file.path(outdir, "07_progeny_14pathway_overview_heatmap_compact.png"),
  p_compact,
  width = 5.6,
  height = 7.4,
  dpi = 300
)
ggsave(
  file.path(outdir, "07_progeny_14pathway_overview_heatmap_compact.pdf"),
  p_compact,
  width = 5.6,
  height = 7.4
)

message("PROGENy 14-pathway overview heatmap saved to: ", outdir)
