#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Figure 4D (ECM/adhesion rows) - NABA matrisome and CAF signatures.
## Standardized stromal program scoring for WT vs FAP-deficient host xenograft.
## Replaces the hand-curated 7-panel "stromal ECM categories" with
## peer-reviewed, citable gene sets:
##
##   NABA matrisome (Naba et al., Mol Cell Proteomics 2012; MSigDB C2:CGP)
##     * NABA_CORE_MATRISOME
##     * NABA_MATRISOME_ASSOCIATED
##     * NABA_ECM_GLYCOPROTEINS
##     * NABA_ECM_REGULATORS
##     * NABA_SECRETED_FACTORS
##     * NABA_COLLAGENS (subset, hardcoded from human set + ortholog)
##     * NABA_PROTEOGLYCANS (subset, hardcoded)
##
##   CAF subtype signatures (peer-reviewed):
##     * Elyada 2019 (Cancer Discov) PDAC  : myCAF, iCAF, apCAF
##     * Dominguez 2020 (Cancer Discov)    : LRRC15+ myCAF, universal CAF
##     * Bartoschek 2018 (Nat Commun)      : vCAF, mCAF, dCAF, cCAF
##     * Kieffer 2020 (Cancer Discov) pan  : CAF-S1 (ecm-myCAF/wound-myCAF)
##
## Each signature is scored per sample on:
##   - bulk RNA (log2 FPKM+1)
##   - DIA proteomics (log2 intensity, NA-imputed)
## using decoupleR::run_ulm (univariate linear model), then compared
## between WT and FAP-KO.
##
## Output: 4-panel NC-style figure + tidy tables suitable as Fig. 8e
## (replacement) or Fig. S8 supplementary.


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(msigdbr); library(decoupleR)
})

proj     <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
bulk_out <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
out_dir  <- file.path(bulk_out, "naba_caf_signatures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_wt  <- "#3B6B9C"; col_fap <- "#C9493A"; col_up <- col_fap; col_dn <- col_wt
col_purple <- "#4B1E70"
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
      axis.line     = element_line(linewidth = 0.6),
      axis.ticks    = element_line(linewidth = 0.6),
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

to_mgi <- function(x){
  x <- toupper(x); paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))))
}

## ============================================================
## 1. Build the signature library
## ============================================================

## ---- NABA from msigdbr (ortholog-mapped to mouse) ----
naba_msig <- suppressMessages(msigdbr(species = "Mus musculus", collection = "C2"))
naba_sets <- c("NABA_CORE_MATRISOME",
               "NABA_MATRISOME_ASSOCIATED",
               "NABA_ECM_GLYCOPROTEINS",
               "NABA_ECM_REGULATORS",
               "NABA_SECRETED_FACTORS")
naba_df <- naba_msig %>%
  filter(gs_name %in% naba_sets) %>%
  transmute(signature = sub("^NABA_", "Naba_", gs_name),
            gene = gene_symbol,
            reference = "Naba 2012 MSigDB") %>%
  distinct()

## ---- Extra NABA subsets (collagens, proteoglycans, ECM-affiliated) ----
## Naba 2012 original classifications, using canonical members:
## (These specific subsets are not in the mouse ortholog C2:CGP list.
##  We hardcode the full human lists and lowercase for MGI symbol
##  convention.)
naba_collagens_h <- c("COL1A1","COL1A2","COL2A1","COL3A1","COL4A1","COL4A2","COL4A3","COL4A4","COL4A5","COL4A6",
                      "COL5A1","COL5A2","COL5A3","COL6A1","COL6A2","COL6A3","COL6A5","COL6A6","COL7A1","COL8A1",
                      "COL8A2","COL9A1","COL9A2","COL9A3","COL10A1","COL11A1","COL11A2","COL12A1","COL13A1",
                      "COL14A1","COL15A1","COL16A1","COL17A1","COL18A1","COL19A1","COL20A1","COL21A1","COL22A1",
                      "COL23A1","COL24A1","COL25A1","COL26A1","COL27A1","COL28A1")
naba_proteogl_h  <- c("ACAN","ASPN","BCAN","BGN","CHAD","CHADL","DCN","EPYC","ESM1","FMOD","HAPLN1","HAPLN2",
                      "HAPLN3","HAPLN4","HSPG2","IMPG1","IMPG2","KERA","LUM","NYX","OGN","OMD","OPTC","PODN",
                      "PODNL1","PRELP","PRG2","PRG3","PRG4","SPOCK1","SPOCK2","SPOCK3","SRGN","VCAN")
naba_ecm_aff_h   <- c("ANXA1","ANXA2","ANXA3","ANXA4","ANXA5","ANXA6","ANXA7","ANXA8","ANXA9","ANXA10",
                      "ANXA11","ANXA13","C1QA","C1QB","C1QC","C1QTNF1","C1QTNF3","CLEC3B","CTSA","CTSB",
                      "CTSC","CTSD","CTSE","CTSF","CTSG","CTSH","CTSK","CTSL","CTSS","CTSW","CTSZ",
                      "GPC1","GPC2","GPC3","GPC4","GPC5","GPC6","MUC1","MUC2","MUC4","MUC5AC","MUC5B",
                      "MUC6","MUC7","MUC12","MUC13","MUC15","MUC16","MUC17","MUC20","MUC21","SLIT1",
                      "SLIT2","SLIT3","SEMA3A","SEMA3B","SEMA3C","SEMA3D","SEMA3E","SEMA3F","SEMA3G")

naba_extra <- bind_rows(
  tibble(signature = "Naba_COLLAGENS",    gene = to_mgi(naba_collagens_h),
         reference = "Naba 2012 (hardcoded)"),
  tibble(signature = "Naba_PROTEOGLYCANS", gene = to_mgi(naba_proteogl_h),
         reference = "Naba 2012 (hardcoded)"),
  tibble(signature = "Naba_ECM_AFFILIATED", gene = to_mgi(naba_ecm_aff_h),
         reference = "Naba 2012 (hardcoded)")
)

## ---- Published CAF subtype signatures ----
## Elyada et al., Cancer Discov 2019 (10.1158/2159-8290.CD-19-0094)
## PDAC single-cell CAF subtypes
elyada_myCAF  <- c("Acta2","Tagln","Myl9","Tpm1","Tpm2","Mmp11","Postn","Hopx",
                   "Thy1","Gja4","Ctgf","Thbs1","Col1a1","Col3a1","Lrrc15","Fap")
elyada_iCAF   <- c("Il6","Cxcl12","Cxcl1","Cxcl2","Cxcl3","Ccl2","Ccl7","Ccl8",
                   "Has1","Has2","Pdgfra","Lmna","Dpt","Cfd","Clec3b","Ptgs2","Nr4a1")
elyada_apCAF  <- c("H2-Ab1","H2-Aa","H2-Eb1","Cd74","Saa3","Slpi")

## Dominguez et al., Cancer Discov 2020 (10.1158/2159-8290.CD-19-1384)
## LRRC15+ myCAF pan-cancer / TGFb-driven CAF
dominguez_LRRC15 <- c("Lrrc15","Col10a1","Col11a1","Cthrc1","Inhba","Mmp11","Mfap5",
                      "Postn","Sulf1","Comp","Fap","Ltbp2","Tagln","Acta2","Tnc",
                      "Thbs2","Cilp","Fbn1")
## Dominguez "universal CAF" markers
dominguez_univ   <- c("Fap","Thy1","Col1a1","Col1a2","Col3a1","Sparc","Fn1",
                      "Vim","Pdpn","Pdgfra","Pdgfrb","S100a4")

## Bartoschek et al., Nat Commun 2018 (10.1038/s41467-018-07582-3) — breast
bart_vCAF <- c("Nr2f2","Mcam","Rgs5","Mylk","Esam","Higd1b","Cox4i2","Pdgfrb","Thy1","Sept4")
bart_mCAF <- c("Pdgfra","Dcn","Col14a1","Lum","Gsn","Mfap5","Smoc2","Serpinf1","Abca8a","Cxcl14")
bart_dCAF <- c("Scrg1","Cryab","Fabp7","Sox10","Plp1","Mpz","Ngfr")
bart_cCAF <- c("Mki67","Top2a","Stmn1","Birc5","Cenpa","Cenpf")

## Kieffer 2020 (Cancer Discov) — CAF-S1 ecm-myCAF (breast/HNSCC)
kieffer_ecm_myCAF <- c("Fap","Lrrc15","Col1a1","Col1a2","Col3a1","Mfap5","Postn",
                       "Mmp11","Comp","Tnc","Pdpn","Thy1","S100a4")

caf_sets <- bind_rows(
  tibble(signature = "Elyada_myCAF",  gene = elyada_myCAF, reference = "Elyada 2019 Cancer Discov"),
  tibble(signature = "Elyada_iCAF",   gene = elyada_iCAF,  reference = "Elyada 2019 Cancer Discov"),
  tibble(signature = "Elyada_apCAF",  gene = elyada_apCAF, reference = "Elyada 2019 Cancer Discov"),
  tibble(signature = "Dominguez_LRRC15_myCAF", gene = dominguez_LRRC15,
         reference = "Dominguez 2020 Cancer Discov"),
  tibble(signature = "Dominguez_universal_CAF", gene = dominguez_univ,
         reference = "Dominguez 2020 Cancer Discov"),
  tibble(signature = "Bartoschek_vCAF",  gene = bart_vCAF, reference = "Bartoschek 2018 Nat Commun"),
  tibble(signature = "Bartoschek_mCAF",  gene = bart_mCAF, reference = "Bartoschek 2018 Nat Commun"),
  tibble(signature = "Bartoschek_dCAF",  gene = bart_dCAF, reference = "Bartoschek 2018 Nat Commun"),
  tibble(signature = "Bartoschek_cCAF",  gene = bart_cCAF, reference = "Bartoschek 2018 Nat Commun"),
  tibble(signature = "Kieffer_ecm_myCAF", gene = kieffer_ecm_myCAF,
         reference = "Kieffer 2020 Cancer Discov")
)

## ---- Hallmark fibrosis/EMT as context ----
gmt <- readLines(file.path(proj, "reference/mh.all.v2026.1.Mm.symbols.gmt"))
hm_list <- lapply(strsplit(gmt, "\t"), function(x){
  if (length(x) >= 3) setNames(list(x[-(1:2)]), x[1]) else NULL
}) %>% purrr::list_flatten()
hm_focus <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
              "HALLMARK_TGF_BETA_SIGNALING",
              "HALLMARK_ANGIOGENESIS",
              "HALLMARK_MYOGENESIS",
              "HALLMARK_INFLAMMATORY_RESPONSE",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE")
hallmark_df <- imap_dfr <- bind_rows(
  lapply(hm_focus, function(n){
    tibble(signature = sub("^HALLMARK_", "Hallmark_", n),
           gene = hm_list[[n]],
           reference = "MSigDB Hallmark (Liberzon 2015)")
  })
)

## ---- consolidate ----
all_sigs <- bind_rows(naba_df, naba_extra, caf_sets, hallmark_df) %>%
  mutate(collection = case_when(
    grepl("^Naba_",  signature) ~ "NABA matrisome",
    grepl("^Hallmark_", signature) ~ "MSigDB Hallmark",
    TRUE ~ "CAF subtype (published)"
  ))
write_tsv(all_sigs, file.path(out_dir, "signature_library.tsv"))
cat(sprintf("Signature library: %d signatures, %d gene-signature rows\n",
            length(unique(all_sigs$signature)), nrow(all_sigs)))

## ============================================================
## 2. Load expression matrices
## ============================================================
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

rna_cols <- paste0("FPKM.", samp_order)
rna_mat <- fpkm_tab %>%
  select(gene_name, all_of(rna_cols)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>%
  slice_max(rowSums(across(all_of(rna_cols))), n = 1, with_ties = FALSE) %>%
  ungroup() %>% column_to_rownames("gene_name") %>% as.matrix()
colnames(rna_mat) <- samp_order
rna_log <- log2(rna_mat + 1)

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

## ============================================================
## 3. Score signatures using decoupleR::run_ulm
## ============================================================
net_df <- all_sigs %>%
  transmute(source = signature, target = gene, mor = 1)

score_mod <- function(mat, tag){
  res <- suppressMessages(decoupleR::run_ulm(
    mat = mat, net = net_df,
    .source = "source", .target = "target", .mor = "mor",
    minsize = 5))
  res %>% transmute(signature = source, sample = condition,
                    score = score, p_value = p_value, dataset = tag)
}

ulm_rna  <- score_mod(rna_log,  "RNA")
ulm_prot <- score_mod(prot_log, "Proteomics")
ulm_all  <- bind_rows(ulm_rna, ulm_prot) %>% left_join(samp_grp, by = "sample")

write_tsv(ulm_all, file.path(out_dir, "signature_scores_per_sample.tsv"))

## group-level stats
stats <- ulm_all %>% group_by(dataset, signature) %>%
  summarise(
    wt_mean   = mean(score[group == "WT host"], na.rm = TRUE),
    fap_mean  = mean(score[group == "FAP-deficient host"], na.rm = TRUE),
    delta_fap_minus_wt = fap_mean - wt_mean,
    wilcox_p  = tryCatch(wilcox.test(score ~ group)$p.value,
                         error = function(e) NA_real_),
    .groups = "drop"
  ) %>% mutate(
    bh_q = p.adjust(wilcox_p, method = "BH"),
    collection = case_when(
      grepl("^Naba_",  signature) ~ "NABA matrisome",
      grepl("^Hallmark_", signature) ~ "MSigDB Hallmark",
      TRUE ~ "CAF subtype (published)"
    )
  ) %>%
  arrange(collection, wilcox_p)
write_tsv(stats, file.path(out_dir, "signature_stats_by_group.tsv"))

cat("\n==== top 8 most-significant signatures per collection ====\n")
print(stats %>% group_by(collection) %>% slice_head(n = 8))

## ============================================================
## 4. Panels
## ============================================================
## (A) delta x modality heatmap — NABA
naba_stats <- stats %>% filter(collection == "NABA matrisome")
naba_wide <- naba_stats %>%
  transmute(signature, dataset, delta_fap_minus_wt) %>%
  pivot_wider(names_from = dataset, values_from = delta_fap_minus_wt)

heatmap_delta <- function(df, title_str, tag_str){
  dset <- df %>% transmute(signature, dataset, delta_fap_minus_wt,
                           wilcox_p) %>%
    mutate(sig_tag = case_when(wilcox_p < 0.05 ~ "*",
                               wilcox_p < 0.1  ~ ".",
                               TRUE ~ ""),
           dataset = factor(dataset, levels = c("RNA", "Proteomics")))
  ord <- df %>% filter(dataset == "RNA") %>%
    arrange(delta_fap_minus_wt) %>% pull(signature)
  dset <- dset %>% mutate(signature = factor(signature, levels = rev(ord)))
  max_abs <- max(abs(dset$delta_fap_minus_wt), na.rm = TRUE)
  ggplot(dset, aes(dataset, signature, fill = delta_fap_minus_wt)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%s%.2f", sig_tag, delta_fap_minus_wt),
                  colour = abs(delta_fap_minus_wt) > max_abs * 0.6),
              size = (base_font_size - 1.5)/.pt, show.legend = FALSE) +
    scale_fill_gradient2(low = col_dn, mid = "#F5F5F5", high = col_up,
                         midpoint = 0, limits = c(-max_abs, max_abs),
                         oob = squish, name = "delta\nFAP-KO\nminus WT") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, title = title_str, tag = tag_str,
         subtitle = "* P<0.05, . P<0.1 (Wilcoxon)") +
    theme_nc() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_text(size = base_font_size - 0.5),
          legend.key.height = unit(8, "mm"),
          legend.key.width  = unit(2.6, "mm"))
}

p_A <- heatmap_delta(naba_stats, "NABA matrisome activities", "a")
p_B <- heatmap_delta(stats %>% filter(collection == "CAF subtype (published)"),
                     "CAF-subtype signature activities", "b")
p_C <- heatmap_delta(stats %>% filter(collection == "MSigDB Hallmark"),
                     "MSigDB Hallmark (fibrosis/EMT context)", "c")

## (D) per-sample score dot plot for the ONE most-decisive NABA & CAF signature
spotlight_sigs <- c(
  naba_stats %>% filter(dataset == "RNA") %>% arrange(wilcox_p) %>%
    slice_head(n = 2) %>% pull(signature),
  stats %>% filter(collection == "CAF subtype (published)",
                    dataset == "RNA") %>%
    arrange(wilcox_p) %>% slice_head(n = 2) %>% pull(signature)
)
spot_df <- ulm_all %>% filter(signature %in% spotlight_sigs) %>%
  mutate(signature = factor(signature, levels = spotlight_sigs),
         dataset = factor(dataset, levels = c("RNA", "Proteomics")))

p_lab <- stats %>% filter(signature %in% spotlight_sigs) %>%
  mutate(signature = factor(signature, levels = spotlight_sigs),
         dataset = factor(dataset, levels = c("RNA", "Proteomics")),
         label = vapply(wilcox_p, p_label, character(1)))

p_D <- ggplot(spot_df, aes(group, score, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 0.8,
              shape = 21, stroke = 0.2, colour = "black") +
  geom_text(data = p_lab, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.3,
            size = (base_font_size - 1.5)/.pt) +
  facet_grid(signature ~ dataset, scales = "free_y", switch = "y") +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_x_discrete(labels = c("WT host" = "WT", "FAP-deficient host" = "KO")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
  labs(x = NULL, y = "ULM activity score",
       title = "Top standardized signatures per sample",
       tag = "d") +
  theme_nc() +
  theme(strip.text.y.left = element_text(angle = 0, size = base_font_size - 1.5,
                                         face = "plain"),
        strip.placement = "outside",
        panel.spacing.y = unit(1.5, "mm"))

## assemble
design <- "
AAABB
AAABB
AAABB
CCCDD
CCCDD
CCCDD
"
fig <- p_A + p_B + p_C + p_D +
  plot_layout(design = design) &
  theme(plot.tag.position = c(0.01, 1.02))

pdf_path <- file.path(out_dir, "standardized_stromal_signatures.pdf")
png_path <- file.path(out_dir, "standardized_stromal_signatures.png")

ggsave(pdf_path, fig, width = 240, height = 240, units = "mm",
       device = grDevices::pdf, useDingbats = FALSE)
ggsave(png_path, fig, width = 240, height = 240, units = "mm",
       dpi = 600, device = ragg::agg_png)

cat("\nStandardized-signature outputs written to:\n  ", out_dir, "\n")
