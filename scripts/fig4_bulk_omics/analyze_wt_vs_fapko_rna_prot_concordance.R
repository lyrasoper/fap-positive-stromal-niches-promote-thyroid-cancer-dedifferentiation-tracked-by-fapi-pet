#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S14 - RNA x protein concordance, Fap-deficient vs wild-type host.
## Step 1 of multi-omics integration:
## RNA <-> Protein sign-concordance quadrant + rank aggregation for
## WT vs FAP-deficient host BPC xenograft.
##
## Deliverables:
##   * 4-quadrant scatter log2FC(RNA) vs log2FC(Protein), labeled hits
##   * Quadrant count table + Fisher enrichment for sign concordance
##   * Robust rank aggregation (RobustRankAggreg + Fisher's combined p)
##   * Ranked "core FAP-dependent" gene table (union & intersection)
##   * Category-stratified concordance summary
##   * Composite 4-panel figure, publication-ready


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(RobustRankAggreg)
})

proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
bulk_out <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
out_dir  <- file.path(bulk_out, "rna_protein_concordance")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_wt  <- "#3B6B9C"; col_fap <- "#C9493A"; col_up <- col_fap; col_dn <- col_wt
col_grey <- "#B7B7B7"; col_purple <- "#4B1E70"

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

## ---- load & harmonize ---------------------------------------------------
deg_rna  <- read_tsv(file.path(bulk_out, "rna/rna_differential_expression.tsv"),
                     show_col_types = FALSE)
deg_prot <- read_tsv(file.path(bulk_out, "proteomics/proteomics_differential_expression.tsv"),
                     show_col_types = FALSE)

deg_rna  <- deg_rna  %>% mutate(gene_mgi = to_mgi(gene))
deg_prot <- deg_prot %>% mutate(gene_mgi = to_mgi(gene))

## inner-join on MGI symbols
joined <- deg_rna %>%
  transmute(gene_rna = gene, gene_mgi, rna_log2fc = log2_fc,
            rna_p = p_value, rna_q = q_value,
            rna_sig = significant == "yes") %>%
  inner_join(
    deg_prot %>% transmute(gene_prot = gene, gene_mgi,
                           prot_log2fc = log2_fc, prot_p = p_value,
                           prot_q = q_value,
                           prot_sig = significant == "yes"),
    by = "gene_mgi"
  ) %>% filter(!is.na(rna_log2fc), !is.na(prot_log2fc))

cat(sprintf("Shared (by MGI symbol) RNA-Protein features: %d\n", nrow(joined)))

## ---- quadrant analysis --------------------------------------------------
joined <- joined %>%
  mutate(rna_dir  = ifelse(rna_log2fc  >= 0, "RNA up",  "RNA down"),
         prot_dir = ifelse(prot_log2fc >= 0, "Prot up", "Prot down"),
         quadrant = paste(rna_dir, prot_dir, sep = " / "),
         sign_concordant = sign(rna_log2fc) == sign(prot_log2fc),
         both_sig = rna_sig & prot_sig,
         either_sig = rna_sig | prot_sig)

qd <- joined %>% count(quadrant) %>% arrange(desc(n))

## Fisher's exact test: is sign concordance enriched over chance?
sig_genes <- joined %>% filter(either_sig)
fe_tab <- with(sig_genes,
               table(RNA = rna_log2fc >= 0, Prot = prot_log2fc >= 0))
fe <- fisher.test(fe_tab)

## Spearman & Pearson correlation of log2FCs
rho_s <- suppressWarnings(cor(joined$rna_log2fc, joined$prot_log2fc, method = "spearman"))
rho_p <- suppressWarnings(cor(joined$rna_log2fc, joined$prot_log2fc, method = "pearson"))

cat(sprintf("Spearman = %.3f | Pearson = %.3f | Fisher OR (sig only) = %.2f | P = %.1e\n",
            rho_s, rho_p, fe$estimate, fe$p.value))

## ---- rank aggregation (RRA) --------------------------------------------
## Build two ranked gene lists sorted by evidence strength; since we have
## per-modality p-values, rank by p_value ascending, break ties by |log2fc|.
rank_list_rna <- joined %>% arrange(rna_p, desc(abs(rna_log2fc))) %>% pull(gene_mgi)
rank_list_prot <- joined %>% arrange(prot_p, desc(abs(prot_log2fc))) %>% pull(gene_mgi)

rra_res <- RobustRankAggreg::aggregateRanks(
  glist = list(RNA = rank_list_rna, Prot = rank_list_prot),
  N = length(rank_list_rna)
) %>% as_tibble() %>%
  transmute(gene_mgi = Name, rra_score = Score,
            rra_p = p.adjust(Score, method = "BH"))

## Fisher's combined p
fisher_combine <- function(p1, p2){
  p1 <- pmax(p1, 1e-300); p2 <- pmax(p2, 1e-300)
  stat <- -2 * (log(p1) + log(p2))
  pchisq(stat, df = 4, lower.tail = FALSE)
}

joined <- joined %>% left_join(rra_res, by = "gene_mgi") %>%
  mutate(fisher_p = fisher_combine(rna_p, prot_p),
         fisher_q = p.adjust(fisher_p, method = "BH"))

core_concordant <- joined %>%
  filter(sign_concordant) %>%
  arrange(fisher_p) %>%
  mutate(combined_effect = sign(rna_log2fc) * (abs(rna_log2fc) + abs(prot_log2fc)) / 2,
         rank_combined = rank(-abs(combined_effect)))

write_tsv(joined,         file.path(out_dir, "rna_prot_joined_all.tsv"))
write_tsv(core_concordant, file.path(out_dir, "core_concordant_ranked.tsv"))
write_tsv(qd,             file.path(out_dir, "quadrant_counts.tsv"))

## ---- category annotation for top genes -----------------------------------
stromal_over <- read_tsv(file.path(bulk_out, "integrated/19_stromal_ecm_marker_overview.tsv"),
                         show_col_types = FALSE)
cat_lut <- stromal_over %>%
  transmute(gene_mgi = to_mgi(gene), category)

## curated lineage / tumor genes
lineage_genes <- to_mgi(c("PAX8","NKX2-1","TG","TPO","DIO1","DIO2","FOXE1",
                          "SLC5A5","SLC5A8","TSHR","THRA","THRB","DUOX1","DUOX2"))
mapk_genes   <- to_mgi(c("DUSP4","DUSP5","DUSP6","ETV4","ETV5","SPRY2","SPRY4"))

annot_tbl <- bind_rows(
  tibble(gene_mgi = lineage_genes, category = "Thyroid lineage"),
  tibble(gene_mgi = mapk_genes,    category = "MAPK output"),
  cat_lut
) %>% distinct(gene_mgi, .keep_all = TRUE)

joined <- joined %>% left_join(annot_tbl, by = "gene_mgi") %>%
  mutate(category = ifelse(is.na(category), "other", category))

## ---- Panel A: 4-quadrant scatter ---------------------------------------
lbl_genes <- joined %>%
  filter(sign_concordant, rna_sig) %>%
  filter(!(category == "other" & abs(rna_log2fc) < 1.5)) %>%
  arrange(fisher_p) %>% slice_head(n = 30) %>%
  pull(gene_mgi)

lbl_genes <- union(lbl_genes,
                   to_mgi(c("Fap","Pax8","Tpo","Tg","Postn","Col1a1","Col3a1",
                            "Lrrc15","Lefty1","Gdnf","Mmp2","Aebp1","Cxcl12")))

plot_df <- joined %>% mutate(is_label = gene_mgi %in% lbl_genes,
                             flag = case_when(
                               both_sig & sign_concordant ~ "both sig, concordant",
                               rna_sig & sign_concordant ~ "RNA sig, concordant",
                               either_sig & !sign_concordant ~ "discordant (either sig)",
                               TRUE ~ "ns"))
plot_df$flag <- factor(plot_df$flag, levels = c(
  "both sig, concordant", "RNA sig, concordant",
  "discordant (either sig)", "ns"
))

count_box <- joined %>% count(rna_dir, prot_dir, sign_concordant) %>%
  mutate(x = ifelse(rna_dir == "RNA up",  Inf, -Inf),
         y = ifelse(prot_dir == "Prot up", Inf, -Inf),
         hj = ifelse(rna_dir == "RNA up",  1.1, -0.1),
         vj = ifelse(prot_dir == "Prot up", 1.4, -0.4),
         lbl = sprintf("n = %d", n))

max_fc <- max(abs(c(joined$rna_log2fc, joined$prot_log2fc)), na.rm = TRUE) * 0.95

p_A <- ggplot(plot_df, aes(rna_log2fc, prot_log2fc)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey45") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              linewidth = 0.3, colour = "grey60") +
  geom_point(data = subset(plot_df, flag == "ns"),
             colour = col_grey, size = 0.4, alpha = 0.35) +
  geom_point(data = subset(plot_df, flag == "discordant (either sig)"),
             colour = "#E2B200", size = 0.8, alpha = 0.75) +
  geom_point(data = subset(plot_df, flag == "RNA sig, concordant"),
             colour = col_purple, size = 0.8, alpha = 0.8) +
  geom_point(data = subset(plot_df, flag == "both sig, concordant"),
             colour = "black", fill = col_up, shape = 21,
             size = 1.5, stroke = 0.3) +
  geom_text(data = count_box,
            aes(x = x, y = y, label = lbl,
                hjust = hj, vjust = vj, colour = sign_concordant),
            size = (base_font_size - 0.5)/.pt, fontface = "plain",
            show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = col_purple, `FALSE` = "#E2B200")) +
  geom_text_repel(data = subset(plot_df, is_label),
                  aes(label = gene_mgi),
                  size = (base_font_size - 1.5)/.pt,
                  min.segment.length = 0, segment.size = 0.2,
                  max.overlaps = Inf, box.padding = 0.3,
                  fontface = "italic", colour = "black") +
  coord_cartesian(xlim = c(-max_fc, max_fc),
                  ylim = c(-max_fc, max_fc)) +
  labs(x = "RNA log2FC (FAP-KO vs WT)",
       y = "Protein log2FC (FAP-KO vs WT)",
       title = sprintf("RNA <-> Protein concordance (%d shared features)",
                       nrow(joined)),
       subtitle = sprintf("Spearman %.2f | Pearson %.2f | Fisher sign-concordance OR %.2f (P = %.1e)",
                          rho_s, rho_p, fe$estimate, fe$p.value),
       tag = "A") +
  theme_nc() +
  theme(legend.position = "none")

## ---- Panel B: quadrant bar ---------------------------------------------
bar_df <- count_box %>%
  mutate(label = paste0(rna_dir, "\n", prot_dir),
         quadrant = case_when(
           rna_dir == "RNA up"   & prot_dir == "Prot up"   ~ "Q1 (both up)",
           rna_dir == "RNA down" & prot_dir == "Prot down" ~ "Q3 (both down)",
           rna_dir == "RNA up"   & prot_dir == "Prot down" ~ "Q4 (RNA up / Prot down)",
           rna_dir == "RNA down" & prot_dir == "Prot up"   ~ "Q2 (RNA down / Prot up)",
           TRUE ~ "other"),
         quadrant = factor(quadrant,
                           levels = c("Q2 (RNA down / Prot up)",
                                      "Q3 (both down)",
                                      "Q4 (RNA up / Prot down)",
                                      "Q1 (both up)")))
p_B <- ggplot(bar_df, aes(n, quadrant, fill = sign_concordant)) +
  geom_col(colour = "black", linewidth = 0.25) +
  geom_text(aes(label = n), hjust = -0.2, size = (base_font_size - 1)/.pt) +
  scale_fill_manual(values = c(`TRUE` = col_purple, `FALSE` = "#E2B200"),
                    labels = c(`TRUE` = "concordant", `FALSE` = "discordant"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(x = "Shared genes (n)", y = NULL,
       title = "Quadrant composition", tag = "B") +
  theme_nc() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"))

## ---- Panel C: top core concordant genes bar (log2FC both modalities) ---
top_core <- core_concordant %>% slice_head(n = 30) %>%
  transmute(gene_mgi, rna_log2fc, prot_log2fc,
            direction = ifelse(combined_effect > 0, "up in FAP-KO", "down in FAP-KO")) %>%
  pivot_longer(c(rna_log2fc, prot_log2fc),
               names_to = "modality", values_to = "log2fc") %>%
  mutate(modality = recode(modality,
                           rna_log2fc = "RNA",
                           prot_log2fc = "Protein"),
         modality = factor(modality, levels = c("RNA", "Protein")),
         gene_mgi = factor(gene_mgi,
                           levels = rev(core_concordant$gene_mgi[1:30])))

p_C <- ggplot(top_core, aes(log2fc, gene_mgi, fill = log2fc)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_col(colour = "black", linewidth = 0.2) +
  facet_wrap(~ modality, nrow = 1) +
  scale_fill_gradient2(low = col_dn, mid = "white", high = col_up,
                       midpoint = 0, guide = "none",
                       limits = c(-3, 3), oob = squish) +
  labs(x = "log2FC (FAP-KO vs WT)", y = NULL,
       title = "Top 30 core concordant genes by Fisher combined p",
       tag = "C") +
  theme_nc() +
  theme(axis.text.y = element_text(size = base_font_size - 1.5),
        panel.spacing.x = unit(2, "mm"))

## ---- Panel D: category-stratified concordance --------------------------
cat_stats <- joined %>%
  group_by(category) %>%
  summarise(n = n(),
            concordant = sum(sign_concordant),
            prop_concordant = mean(sign_concordant),
            .groups = "drop") %>%
  arrange(desc(n)) %>% filter(n >= 3) %>%
  mutate(category = factor(category, levels = rev(category)),
         label = sprintf("%d / %d", concordant, n))

p_D <- ggplot(cat_stats, aes(prop_concordant, category, fill = prop_concordant)) +
  geom_col(colour = "black", linewidth = 0.25) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             linewidth = 0.3, colour = "grey45") +
  geom_text(aes(label = label), hjust = -0.15,
            size = (base_font_size - 1)/.pt) +
  scale_fill_gradient(low = "#F5F5F5", high = col_purple,
                      limits = c(0, 1), guide = "none") +
  scale_x_continuous(limits = c(0, 1.15),
                     breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.02, 0.1)),
                     labels = label_percent()) +
  labs(x = "Proportion sign-concordant",
       y = NULL,
       title = "Concordance by annotated category",
       subtitle = "dashed = 50% (chance)",
       tag = "D") +
  theme_nc() +
  theme(axis.text.y = element_text(size = base_font_size - 1.5))

## ---- assemble & save ---------------------------------------------------
design <- "
AAABBB
AAACCC
AAACCC
DDDCCC
DDDCCC
"
fig <- p_A + p_B + p_C + p_D +
  plot_layout(design = design) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(out_dir, "rna_protein_concordance_overview.pdf")
png_path <- file.path(out_dir, "rna_protein_concordance_overview.png")
ggsave(pdf_path, fig, width = 220, height = 220, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(png_path, fig, width = 220, height = 220, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("\nOutputs written to:\n  ", out_dir, "\n")

## numerical summary block for stdout
cat("\n==================== Summary ====================\n")
cat(sprintf("Shared features: %d | RNA-sig: %d | Prot-sig: %d | both-sig: %d\n",
            nrow(joined), sum(joined$rna_sig), sum(joined$prot_sig),
            sum(joined$both_sig)))
cat(sprintf("Sign-concordant: %d (%.1f%%)\n",
            sum(joined$sign_concordant),
            100 * mean(joined$sign_concordant)))
cat(sprintf("Spearman rho: %.3f | Pearson: %.3f\n", rho_s, rho_p))
cat(sprintf("Fisher OR (sig only): %.2f | P: %.2e\n", fe$estimate, fe$p.value))
cat("\nTop 10 core concordant genes (Fisher combined p):\n")
print(core_concordant %>% slice_head(n = 10) %>%
        select(gene_mgi, rna_log2fc, prot_log2fc, rna_q, prot_p,
               fisher_p, fisher_q, rra_score) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))
