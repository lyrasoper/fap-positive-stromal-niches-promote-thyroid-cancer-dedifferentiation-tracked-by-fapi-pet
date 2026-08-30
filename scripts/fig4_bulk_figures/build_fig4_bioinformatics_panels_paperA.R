#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Figure 4C-E - the bioinformatics panels of Figure 4 (paper_A render).
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.
## Fig 4 bioinformatics panels — paper_A standalone re-render (2026-06-13)
## Reuses the EXACT computations from build_fig8_v14_itga5_highlight.R (proven
## against the real WT-vs-FAPKO BPC bulk omics) but renders each 生信 candidate
## panel as an INDIVIDUAL beautified figure so the user can pick 2–5 for main.
## Beautification vs original: heatmap/diverging midpoint pure-white -> off-white
## (#F5F5F5); cleaner standalone titles; italic stats; lowercase bold tags kept;
## ITGA5 ★ kept. NO data/value changes. Panels j/k (B-direction cartoon) and
## i (external phospho) are intentionally excluded from the main-candidate set.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(tibble)
  library(patchwork); library(ggrepel); library(ggtext); library(scales); library(grid)
})

proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
out_root <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
OUT <- file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/_fig4_bioinformatics_panels_20260613")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- palette / theme (paper_A) ----------------------------------------
col_wt   <- "#3B6B9C"; col_fap  <- "#C9493A"; col_grey <- "#B7B7B7"
col_up   <- col_fap;   col_dn   <- col_wt
col_purple <- "#4B1E70"; col_gold <- "#E2B200"
col_mid  <- "#F5F5F5"                       # ← off-white midpoint (was pure white)
col_lineage <- "#C9493A"; col_stromal <- "#3B6B9C"
grp_levels <- c("WT host", "FAP-deficient host")
grp_cols <- setNames(c(col_wt, col_fap), grp_levels)

base_font_size <- 7
theme_nc <- function(base = base_font_size){
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(face="plain", size = base + 1, hjust = 0),
      plot.subtitle = element_text(size = base - 1.2, colour = "grey35"),
      plot.tag      = element_text(face = "bold", size = base + 3, hjust = 0),
      axis.title    = element_text(size = base),
      axis.text     = element_text(size = base - 0.5, colour = "black"),
      axis.line     = element_line(linewidth = 0.3),
      axis.ticks    = element_line(linewidth = 0.3),
      strip.background = element_blank(),
      strip.text    = element_text(size = base, face = "plain"),
      legend.key.size = unit(3.2, "mm"),
      legend.text   = element_text(size = base - 0.5),
      legend.title  = element_text(size = base - 0.5, face = "plain"),
      legend.background = element_blank(), legend.box.background = element_blank(),
      panel.grid    = element_blank(),
      plot.margin   = margin(3, 4, 3, 4)
    )
}
tag_group <- function(df) df %>% mutate(group = factor(ifelse(grepl("^WT", sample), "WT host", "FAP-deficient host"), levels = grp_levels))
p_label <- function(p){ if (is.na(p)) return("ns"); if (p < 0.001) sprintf("italic(P) == %.1e", p) else sprintf("italic(P) == %.3f", p) }
safe_z <- function(x){ sx <- sd(x, na.rm=TRUE); mx <- mean(x, na.rm=TRUE); if (!is.finite(sx) || sx==0) return(rep(0,length(x))); (x-mx)/sx }
to_mgi <- function(x){ x <- toupper(x); paste0(substr(x,1,1), tolower(substr(x,2,nchar(x)))) }

save_panel <- function(p, name, w, h){
  pdf_p <- file.path(OUT, paste0(name, ".pdf")); png_p <- file.path(OUT, paste0(name, ".png"))
  ## quartz pdf = macOS native, UTF-8 safe (renders ★ and en-dash; base pdf does not)
  grDevices::quartz(file = pdf_p, type = "pdf", width = w/25.4, height = h/25.4, family = "Helvetica")
  print(p); dev.off()
  ggsave(png_p, p, width = w, height = h, units = "mm", dpi = 600, device = ragg::agg_png)
  cat("  wrote", name, sprintf("(%gx%g mm)\n", w, h))
}

## ---- data --------------------------------------------------------------
deg_rna  <- read_tsv(file.path(out_root, "rna/rna_differential_expression.tsv"), show_col_types = FALSE)
deg_prot <- read_tsv(file.path(out_root, "proteomics/proteomics_differential_expression.tsv"), show_col_types = FALSE)
joined_all <- read_tsv(file.path(out_root, "rna_protein_concordance/rna_prot_joined_all.tsv"), show_col_types = FALSE)
fpkm_tab <- read_tsv(file.path(proj, "001_mRNA_Summary/1.GeneExpression/1_genes_fpkm_expression.txt"), show_col_types = FALSE, progress = FALSE)
prot_tab <- read_tsv(file.path(proj, "002_DIA_Summary/02.ProteinExp/protein_annotation_profile.txt"), show_col_types = FALSE, progress = FALSE)
naba_scores   <- read_tsv(file.path(out_root, "naba_caf_signatures/signature_scores_per_sample.tsv"), show_col_types = FALSE)
progeny_scores <- read_tsv(file.path(out_root, "decoupler_crossomics/progeny_scores_per_sample.tsv"), show_col_types = FALSE)
dediff_scores <- read_tsv(file.path(out_root, "dediff_axes/integrated/dediff_axes_scores_per_sample.tsv"), show_col_types = FALSE)
nnls_df  <- read_tsv(file.path(out_root, "caf_deconvolution/nnls_cell_fractions.tsv"), show_col_types = FALSE)
ulm_df   <- read_tsv(file.path(out_root, "caf_deconvolution/ulm_activity_scores.tsv"), show_col_types = FALSE)
decon_stats <- read_tsv(file.path(out_root, "caf_deconvolution/deconvolution_group_stats.tsv"), show_col_types = FALSE)
core_tbl <- read_tsv(file.path(out_root, "ecm_mycaf_core_genes/core_gene_stats.tsv"), show_col_types = FALSE)

samp_order <- c("WT_BPC1","WT_BPC2","WT_BPC3","WT_BPC4","FAP_BPC1","FAP_BPC2","FAP_BPC3","FAP_BPC4")
samp_labs <- setNames(c("WT1","WT2","WT3","WT4","KO1","KO2","KO3","KO4"), samp_order)

fpkm_long <- function(gene_symbol){
  hit <- fpkm_tab %>% filter(gene_name == gene_symbol); if (nrow(hit)==0) return(tibble()); hit <- hit[1,]
  fc <- grep("^FPKM\\.", colnames(hit), value=TRUE)
  tibble(sample = sub("^FPKM\\.","",fc), value = as.numeric(hit[1,fc])) %>% tag_group() %>% mutate(gene=gene_symbol)
}
prot_long <- function(gene_symbol){
  hit <- prot_tab %>% filter(toupper(gene_name)==toupper(gene_symbol)); if (nrow(hit)==0) return(tibble()); hit <- hit[1,]
  sc <- samp_order[samp_order %in% colnames(hit)]; if(!length(sc)) return(tibble())
  tibble(sample=sc, value=as.numeric(hit[1,sc])) %>% tag_group() %>% mutate(gene=gene_symbol)
}

stromal_bundle <- c("Fap","Postn","Col1a1","Col3a1","Lrrc15","Mmp2","Mmp11","Thy1","Acta2","Sparc","Aebp1","Cthrc1","Mfap5","Comp","Tnc")
lineage_bundle <- c("Pax8","Tpo","Tg","Foxe1","Nkx2-1","Slc5a5","Dio1","Dio2","Duox2")
all_bundle <- c(stromal_bundle, lineage_bundle)
main_volcano_labels <- c("Fap","Col1a1","Postn","Cthrc1","Tnc","Mmp2","Acta2","Pax8","Tpo","Slc5a5","Duox2")

## =======================================================================
## b — on-target FAP loss
## =======================================================================
deg_fap_rna  <- deg_rna  %>% filter(gene=="Fap")  %>% slice_head(n=1)
deg_fap_prot <- deg_prot %>% filter(gene=="FAP")  %>% slice_head(n=1)
on_target_df <- tibble(
  readout = factor(c("Fap mRNA","FAP protein"), levels=c("FAP protein","Fap mRNA")),
  log2fc = c(deg_fap_rna$log2_fc[1], deg_fap_prot$log2_fc[1]),
  stat = c(sprintf("italic(q) == %.1e", deg_fap_rna$q_value[1]), sprintf("italic(P) == %.2g", deg_fap_prot$p_value[1])))
p_b <- ggplot(on_target_df, aes(log2fc, readout)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey55") +
  geom_segment(aes(x=0, xend=log2fc, yend=readout), linewidth=1.1, colour=col_wt, lineend="round") +
  geom_point(shape=21, size=3.2, stroke=0.35, fill=col_wt, colour="black") +
  geom_text(aes(label=sprintf("%.2f", log2fc)), position=position_nudge(y=0.30), size=(base_font_size-0.8)/.pt, fontface="bold", colour="grey15") +
  geom_text(aes(x=1.05, label=stat), hjust=1, size=(base_font_size-1.7)/.pt, colour="grey30", parse=TRUE) +
  coord_cartesian(xlim=c(floor(min(on_target_df$log2fc))-0.45, 1.1), clip="off") +
  labs(x="log2FC (FAP-KO vs WT)", y=NULL, title="On-target stromal FAP loss",
       subtitle="n=4/grp · transcript loss robust; protein concordant, n.s. (near DIA floor)", tag="b") +
  theme_nc() + theme(axis.text.y=element_text(face="bold", size=base_font_size-0.6), panel.grid.major.y=element_blank())
save_panel(p_b, "panel_b_on_target_FAP_loss", 58, 40)

## =======================================================================
## c — RNA volcano (ECM-lineage axis)
## =======================================================================
volcano_prep <- function(deg, gene_col="gene", bundle_mgi=TRUE){
  gv <- deg[[gene_col]]; lk <- if (bundle_mgi) to_mgi(gv) else gv
  deg %>% mutate(gene_mgi=lk, log10p=-log10(pmax(p_value,.Machine$double.xmin)),
    bundle=case_when(gene_mgi %in% stromal_bundle ~ "stromal / ECM",
                     gene_mgi %in% lineage_bundle ~ "thyroid lineage",
                     significant=="yes" & log2_fc>=0 ~ "other up (FDR)",
                     significant=="yes" & log2_fc<0  ~ "other down (FDR)", TRUE ~ "ns"),
    bundle=factor(bundle, levels=c("stromal / ECM","thyroid lineage","other up (FDR)","other down (FDR)","ns")),
    lab=ifelse(gene_mgi %in% main_volcano_labels, gene_mgi, NA_character_))
}
vol_rna <- volcano_prep(deg_rna, "gene", FALSE) %>% mutate(log10p=pmin(log10p,18))
bundle_cols <- c("stromal / ECM"=col_stromal,"thyroid lineage"=col_lineage,"other up (FDR)"="#F4C7B8","other down (FDR)"="#B8CBE0","ns"=col_grey)
p_c_data <- vol_rna %>% mutate(log2_fc_plot=pmax(pmin(log2_fc,4),-4))
p_c <- ggplot(p_c_data, aes(log2_fc_plot, log10p, colour=bundle)) +
  annotate("rect", xmin=-4, xmax=-1, ymin=-Inf, ymax=Inf, fill=alpha(col_stromal,0.055), colour=NA) +
  annotate("rect", xmin=1, xmax=4, ymin=-Inf, ymax=Inf, fill=alpha(col_lineage,0.05), colour=NA) +
  annotate("text", x=-3.25, y=11.85, label="ECM/FAP niche\ncollapse", size=(base_font_size-1.6)/.pt, colour=col_stromal, fontface="bold", lineheight=0.9) +
  annotate("text", x=3.20, y=11.85, label="lineage\nrebound", size=(base_font_size-1.6)/.pt, colour=col_lineage, fontface="bold", lineheight=0.9) +
  geom_point(data=subset(p_c_data, bundle=="ns"), size=0.4, alpha=0.45, stroke=0) +
  geom_point(data=subset(p_c_data, bundle=="other down (FDR)"), size=0.55, alpha=0.65, stroke=0) +
  geom_point(data=subset(p_c_data, bundle=="other up (FDR)"), size=0.55, alpha=0.65, stroke=0) +
  geom_point(data=subset(p_c_data, bundle %in% c("stromal / ECM","thyroid lineage")), aes(fill=bundle), shape=21, size=1.6, stroke=0.35, colour="black", alpha=0.95) +
  geom_text_repel(data=subset(p_c_data, !is.na(lab)), aes(label=lab), size=(base_font_size-1.5)/.pt, fontface="italic", segment.size=0.2, segment.colour="grey35", min.segment.length=0, max.overlaps=Inf, box.padding=0.25, show.legend=FALSE) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", linewidth=0.22, colour="grey50") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth=0.22, colour="grey50") +
  scale_colour_manual(values=bundle_cols, name=NULL) + scale_fill_manual(values=bundle_cols, guide="none") +
  scale_x_continuous(breaks=c(-4,-2,0,2,4), labels=c("<-4","-2","0","2",">4")) +
  scale_y_continuous(expand=expansion(mult=c(0.01,0.06))) +
  labs(x="log2 FC (FAP-KO vs WT host)", y=expression(-log[10]~italic(P)),
       title="FAP loss redirects tumours along an ECM–lineage axis",
       subtitle="Genome-wide RNA; mechanism zones highlighted", tag="c") +
  theme_nc() + theme(legend.position="top", legend.direction="horizontal", legend.box.spacing=unit(1,"mm"),
                     legend.text=element_text(size=base_font_size-1.6), legend.key.size=unit(2.4,"mm"),
                     legend.spacing.x=unit(0.5,"mm"))
save_panel(p_c, "panel_c_volcano_ecm_lineage", 100, 80)

## =======================================================================
## d — mechanism-gene module shift + focused gene heatmap
## =======================================================================
mechanism_gene_tbl <- tibble(
  gene_mgi=c("Fap","Cthrc1","Col1a1","Postn","Mfap5","Col3a1","Pax8","Tpo","Slc5a5","Duox2"),
  display=c("FAP","CTHRC1","COL1A1","POSTN","MFAP5","COL3A1","PAX8","TPO","SLC5A5/NIS","DUOX2"),
  class=c(rep("ECM/FAP niche",6), rep("thyroid lineage",4)))
d_mech <- mechanism_gene_tbl %>%
  left_join(deg_rna %>% transmute(gene_mgi=gene, RNA=log2_fc, rna_p=p_value, rna_q=q_value), by="gene_mgi") %>%
  left_join(deg_prot %>% mutate(gene_mgi=to_mgi(gene)) %>% transmute(gene_mgi, Protein=log2_fc, prot_p=p_value, prot_q=q_value), by="gene_mgi") %>%
  mutate(display=factor(display, levels=rev(display)), class=factor(class, levels=c("ECM/FAP niche","thyroid lineage")))
d_mech_long <- d_mech %>% pivot_longer(c(RNA,Protein), names_to="modality", values_to="log2fc") %>%
  mutate(modality=factor(modality, levels=c("RNA","Protein")),
         p_value=ifelse(modality=="RNA", rna_p, prot_p), q_value=ifelse(modality=="RNA", rna_q, prot_q),
         detected=!is.na(log2fc), log2fc_plot=pmax(pmin(log2fc,2.2),-2.2),
         tile_label=case_when(!detected~"ND", log2fc>2.2~">2.2", log2fc< -2.2~"<-2.2", TRUE~sprintf("%.1f",log2fc)),
         sig_label=case_when(!detected~"", !is.na(q_value)&q_value<0.1~"*", !is.na(p_value)&p_value<0.1~".", TRUE~""),
         signif=detected & ((!is.na(q_value)&q_value<0.1)|(!is.na(p_value)&p_value<0.1)))
score_rna_module <- function(genes, module_label){
  rc <- paste0("FPKM.", samp_order)
  d <- fpkm_tab %>% select(gene_name, all_of(rc)) %>% mutate(gu=toupper(gene_name)) %>% filter(gu %in% toupper(genes)) %>%
    group_by(gu) %>% summarise(across(all_of(rc), ~mean(as.numeric(.x),na.rm=TRUE)), .groups="drop")
  mat <- log2(as.matrix(d[,rc,drop=FALSE])+1); z <- t(scale(t(mat))); z[!is.finite(z)] <- 0
  tibble(sample=samp_order, group=factor(ifelse(grepl("^WT",samp_order),"WT host","FAP-deficient host"),levels=grp_levels),
         module=module_label, score_z=colMeans(z,na.rm=TRUE))
}
d_module_scores <- bind_rows(score_rna_module(c("Fap","Cthrc1","Col1a1","Postn","Mfap5","Col3a1"),"ECM/FAP niche"),
                             score_rna_module(c("Pax8","Tpo","Slc5a5","Duox2"),"thyroid lineage")) %>%
  mutate(module=factor(module, levels=c("thyroid lineage","ECM/FAP niche")))
d_module_stats <- d_module_scores %>% group_by(module) %>%
  group_modify(~tibble(delta=mean(.x$score_z[.x$group=="FAP-deficient host"])-mean(.x$score_z[.x$group=="WT host"]))) %>%
  ungroup() %>% mutate(colour_key=ifelse(delta<0,"ECM/FAP niche","thyroid lineage"))
p_d_module <- ggplot(d_module_stats, aes(delta, module)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey45") +
  geom_segment(aes(x=0, xend=delta, yend=module), linewidth=0.65, colour="grey35") +
  geom_point(aes(fill=colour_key), shape=21, size=2.6, stroke=0.25, colour="black") +
  scale_fill_manual(values=c("ECM/FAP niche"=col_stromal,"thyroid lineage"=col_lineage), guide="none") +
  scale_x_continuous(limits=c(-1,1), breaks=c(-1,0,1)) +
  labs(x="KO - WT module z", y=NULL, title="Module shift", subtitle="RNA score") +
  theme_nc() + theme(plot.title=element_text(size=base_font_size), axis.text.y=element_text(face="bold", size=base_font_size-1.2),
                     panel.grid.major.y=element_line(colour="grey92", linewidth=0.25))
p_d_heat <- ggplot(d_mech_long, aes(modality, display, fill=log2fc_plot)) +
  geom_tile(aes(alpha=signif), colour="white", linewidth=0.35) +
  scale_alpha_manual(values=c(`FALSE`=0.32,`TRUE`=1), guide="none") +
  geom_text(aes(label=tile_label), size=(base_font_size-2.0)/.pt, colour="black") +
  geom_text(aes(label=sig_label), nudge_x=0.34, nudge_y=0.27, size=(base_font_size-1.3)/.pt, colour="black") +
  facet_grid(class~., scales="free_y", space="free_y", switch="y") +
  scale_fill_gradient2(low=col_dn, mid=col_mid, high=col_up, midpoint=0, limits=c(-2.2,2.2), oob=squish, na.value="#EEEEEE", name="log2FC") +
  labs(x=NULL, y=NULL, title="Focused gene deltas", subtitle="RNA + DIA; ND, faded=n.s.") +
  theme_nc() + theme(axis.line=element_blank(), axis.ticks=element_blank(),
                     axis.text.x=element_text(size=base_font_size-0.7, face="bold"),
                     axis.text.y=element_text(size=base_font_size-0.9, face="italic"),
                     strip.placement="outside", strip.text.y.left=element_text(angle=0, face="bold", size=base_font_size-1.2),
                     panel.spacing.y=unit(1.3,"mm"), legend.key.width=unit(2.2,"mm"), legend.key.height=unit(5,"mm"))
p_d <- (p_d_module | p_d_heat) + plot_layout(widths=c(0.68,1.32))
p_d <- wrap_elements(full=p_d) + labs(tag="d")
save_panel(p_d, "panel_d_mechanism_gene_deltas", 110, 66)

## =======================================================================
## e1 — CAF composition stable (NNLS)
## =======================================================================
ct_cols <- c("myCAF"="#C9493A","iCAF"="#E08F3A","apCAF"="#5E4A7B","Thyroid_tumor"="#3B6B9C","Immune"="#88B0A7","Endothelial"="#B7B7B7")
ct_levels <- c("myCAF","iCAF","apCAF","Thyroid_tumor","Immune","Endothelial")
mycaf_p <- decon_stats %>% filter(celltype=="myCAF", method=="NNLS fraction") %>% pull(wilcox_p)
icaf_p  <- decon_stats %>% filter(celltype=="iCAF",  method=="NNLS fraction") %>% pull(wilcox_p)
e1_df <- decon_stats %>% filter(method=="NNLS fraction") %>%
  mutate(celltype=factor(celltype, levels=rev(ct_levels)), delta_pp=100*delta_fap_minus_wt,
         sig_tag=case_when(wilcox_p<0.05~"*", wilcox_p<0.1~".", TRUE~""))
p_e1 <- ggplot(e1_df, aes(delta_pp, celltype)) +
  geom_vline(xintercept=0, linewidth=0.35, colour="grey45") +
  geom_segment(aes(x=0, xend=delta_pp, yend=celltype), linewidth=0.55, colour="grey45") +
  geom_point(aes(fill=celltype), shape=21, size=2.4, stroke=0.25, colour="black") +
  geom_text(aes(label=sig_tag), nudge_x=0.2, size=(base_font_size+0.5)/.pt) +
  scale_fill_manual(values=ct_cols, guide="none") +
  scale_x_continuous(labels=function(x) paste0(x," pp"), limits=c(-1.6,1.6), expand=expansion(mult=c(0.08,0.18))) +
  labs(x="Delta fraction (FAP-KO - WT)", y=NULL, title="CAF composition stable after FAP loss",
       subtitle=sprintf("NNLS; myCAF italic(P)==%.2f, iCAF italic(P)==%.2f", mycaf_p, icaf_p), tag="e1") +
  theme_nc() + theme(panel.grid.major.y=element_line(colour="grey92", linewidth=0.25), axis.text.y=element_text(size=base_font_size-0.9))
## subtitle has plotmath via parse-free; render plainly:
p_e1 <- p_e1 + labs(subtitle=sprintf("NNLS; myCAF P=%.2f, iCAF P=%.2f", mycaf_p, icaf_p))
save_panel(p_e1, "panel_e1_caf_composition_stable", 64, 54)

## =======================================================================
## e2 — Fap-selective myCAF scaffold perturbation
## =======================================================================
core12 <- core_tbl %>% filter(n_sig>=2) %>% arrange(rna_log2fc) %>% mutate(gene_mgi=factor(gene_mgi, levels=gene_mgi))
e2_df <- core12 %>% select(gene_mgi, RNA=rna_log2fc, Protein=prot_log2fc) %>%
  pivot_longer(c(RNA,Protein), names_to="modality", values_to="log2fc") %>% mutate(modality=factor(modality, levels=c("RNA","Protein")))
tag_df <- bind_rows(
  core12 %>% select(gene_mgi, rna_q) %>% mutate(modality="RNA",  tag=case_when(rna_q<0.05~"*", rna_q<0.1~".", TRUE~"")),
  core12 %>% select(gene_mgi, prot_p) %>% mutate(modality="Protein", tag=case_when(prot_p<0.05~"*", prot_p<0.1~".", TRUE~""))) %>%
  left_join(e2_df %>% select(gene_mgi, modality, log2fc), by=c("gene_mgi","modality")) %>% mutate(modality=factor(modality, levels=c("RNA","Protein")))
p_e2 <- ggplot(e2_df, aes(log2fc, gene_mgi, fill=log2fc)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey45") + geom_col(colour="black", linewidth=0.2) +
  geom_text(data=tag_df, aes(x=log2fc+sign(log2fc)*0.12, y=gene_mgi, label=tag), inherit.aes=FALSE, size=(base_font_size+0.5)/.pt) +
  facet_wrap(~modality, nrow=1) +
  scale_fill_gradient2(low=col_dn, mid=col_mid, high=col_up, midpoint=0, limits=c(-2.2,2.2), oob=squish, guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0.08,0.1)), breaks=c(-2,-1,0,1)) +
  labs(x="log2 FC (FAP-KO vs WT)", y=NULL, title="Fap-selective myCAF scaffold perturbation",
       subtitle="12 recurrent scaffold genes; broad scaffold retained", tag="e2") +
  theme_nc() + theme(axis.text.y=element_text(face="italic", size=base_font_size-1.2), panel.spacing.x=unit(2,"mm"))
save_panel(p_e2, "panel_e2_mycaf_scaffold", 72, 68)

## =======================================================================
## e3 — tumour-lineage program rises; myCAF activity unchanged (ULM)
## =======================================================================
thy_df <- ulm_df %>% filter(celltype=="Thyroid_tumor") %>% mutate(group=factor(group, levels=grp_levels))
thy_p  <- decon_stats %>% filter(celltype=="Thyroid_tumor", method=="ULM activity") %>% pull(wilcox_p)
myc_df <- ulm_df %>% filter(celltype=="myCAF") %>% mutate(group=factor(group, levels=grp_levels))
myc_p  <- decon_stats %>% filter(celltype=="myCAF", method=="ULM activity") %>% pull(wilcox_p)
combo_df <- bind_rows(thy_df, myc_df) %>% mutate(celltype=factor(celltype, levels=c("myCAF","Thyroid_tumor")))
combo_p <- tibble(celltype=factor(c("myCAF","Thyroid_tumor"), levels=c("myCAF","Thyroid_tumor")),
                  label=c(sprintf("P=%.3f", myc_p), sprintf("P=%.3f", thy_p)))
p_e3 <- ggplot(combo_df, aes(group, score, fill=group)) +
  geom_boxplot(width=0.55, outlier.shape=NA, linewidth=0.3, alpha=0.85) +
  geom_jitter(width=0.12, height=0, size=0.9, shape=21, stroke=0.2, colour="black") +
  geom_text(data=combo_p, aes(x=1.5, y=Inf, label=label), inherit.aes=FALSE, vjust=1.3, size=(base_font_size-1.5)/.pt) +
  facet_wrap(~celltype, nrow=1, scales="free_y", labeller=as_labeller(c(myCAF="myCAF activity", Thyroid_tumor="Tumor-lineage activity"))) +
  scale_fill_manual(values=grp_cols, guide="none") + scale_x_discrete(labels=c("WT host"="WT","FAP-deficient host"="KO")) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.25))) +
  labs(x=NULL, y="ULM signature activity (z)", title="Tumour-lineage program rises", subtitle="myCAF activity unchanged", tag="e3") +
  theme_nc() + theme(strip.text=element_text(size=base_font_size-0.5))
save_panel(p_e3, "panel_e3_ulm_lineage_activity", 66, 54)

## =======================================================================
## f — mechanism axes shift (forest)  ★ TOP MAIN CANDIDATE
## =======================================================================
axis_from_scores <- function(df, value_col="score", name, label, axis_class, modality, source, min_n=3){
  d <- df %>% filter(!is.na(.data[[value_col]]), group %in% grp_levels) %>%
    mutate(group=factor(group, levels=grp_levels), score_for_axis=as.numeric(scale(.data[[value_col]])))
  if (nrow(d) < min_n*2 || any(table(d$group) < min_n)) return(tibble())
  g <- d %>% group_by(group) %>% summarise(mean_z=mean(score_for_axis), sd_z=sd(score_for_axis), n=n(), .groups="drop")
  wt <- g$mean_z[g$group=="WT host"]; ko <- g$mean_z[g$group=="FAP-deficient host"]
  p <- tryCatch(wilcox.test(score_for_axis~group, data=d, exact=FALSE)$p.value, error=function(e) NA_real_)
  tibble(signature=name, axis_label=label, axis_class=axis_class, delta=ko-wt,
         se=sqrt(sum(g$sd_z^2/g$n)), p_value=p)
}
score_custom_rna_axis <- function(genes, name, label, axis_class){
  rc <- paste0("FPKM.", samp_order)
  d <- fpkm_tab %>% select(gene_name, all_of(rc)) %>% mutate(gu=toupper(gene_name)) %>% filter(gu %in% toupper(genes)) %>%
    group_by(gu) %>% summarise(across(all_of(rc), ~mean(as.numeric(.x),na.rm=TRUE)), .groups="drop")
  if (nrow(d)<3) return(tibble())
  mat <- log2(as.matrix(d[,rc,drop=FALSE])+1); z <- t(scale(t(mat))); z[!is.finite(z)] <- 0
  tibble(sample=samp_order, score=colMeans(z,na.rm=TRUE),
         group=factor(ifelse(grepl("^WT",samp_order),"WT host","FAP-deficient host"),levels=grp_levels)) %>%
    axis_from_scores("score", name, label, axis_class, "RNA", "")
}
dediff_axis <- function(sig, label, axis_class) dediff_scores %>% filter(dataset=="RNA", signature==sig, !is.na(score_z)) %>% axis_from_scores("score_z", sig, label, axis_class, "RNA","")
naba_axis <- function(sig, modality, label, axis_class) naba_scores %>% filter(dataset==modality, signature==sig) %>% axis_from_scores("score", sig, label, axis_class, modality, "")
progeny_axis <- function(pid, modality, label, axis_class) progeny_scores %>% filter(dataset==modality, pathway==pid) %>% axis_from_scores("score", paste0("PROGENy_",pid), label, axis_class, modality, "")
all_axes <- bind_rows(
  naba_axis("Kieffer_ecm_myCAF","Proteomics","CAF-ECM fibrosis\n(Kieffer; DIA)","ECM / adhesion niche"),
  naba_axis("Naba_ECM_GLYCOPROTEINS","Proteomics","NABA ECM glycoproteins\n(DIA)","ECM / adhesion niche"),
  naba_axis("Hallmark_EPITHELIAL_MESENCHYMAL_TRANSITION","Proteomics","Hallmark EMT/fibrosis\n(DIA)","ECM / adhesion niche"),
  score_custom_rna_axis(c("PTK2","SRC","PXN","VCL","TLN1"),"FAK/SRC core","FAK/SRC adhesion\ncore (RNA)","ECM / adhesion niche"),
  progeny_axis("MAPK","RNA","PROGENy MAPK output\n(RNA)","MAPK / ERK output"),
  dediff_axis("TDS / lineage identity","TDS / lineage identity\n(TCGA 16-gene; RNA)","lineage / dediff"),
  dediff_axis("BRAF-RAS score","BRAF-RAS score\n(RNA)","lineage / dediff")) %>%
  mutate(q_value=p.adjust(p_value, method="BH"),
         ci_lo=delta-1.96*se, ci_hi=delta+1.96*se,
         sig_tag="",  # honesty: 2-sided Wilcoxon + BH across 7 axes -> no axis q<0.05 at n=4; show q inline, no significance stars
         axis_class=factor(axis_class, levels=c("ECM / adhesion niche","MAPK / ERK output","lineage / dediff")))
axis_levels <- c("BRAF-RAS score","TDS / lineage identity","PROGENy_MAPK","FAK/SRC core","Hallmark_EPITHELIAL_MESENCHYMAL_TRANSITION","Naba_ECM_GLYCOPROTEINS","Kieffer_ecm_myCAF")
all_axes <- all_axes %>% mutate(axis_label=factor(axis_label, levels=all_axes$axis_label[match(axis_levels, all_axes$signature)]))
## beautify: main name normal + assay tag (DIA/RNA) small grey, via ggtext
md_axis_lab <- function(x){
  vapply(strsplit(as.character(x),"\n",fixed=TRUE), function(p){
    main <- p[1]; sub <- if(length(p)>1) gsub("^\\(|\\)$","",p[2]) else ""
    if(nzchar(sub)) sprintf("%s<br><span style='font-size:5pt;color:#8f8f8f'>%s</span>", main, sub) else main
  }, character(1))
}
all_axes <- all_axes %>% mutate(pt_size = 2.0 + 0.8*pmin(abs(delta),1.3)/1.3,
                                val_lab = sprintf("%+.1f", delta),          # effect size Δz
                                p_lab   = sprintf("P=%.2g", p_value))       # exact NOMINAL p (uncorrected); BH-q all >0.05 (no axis FDR-sig)
## dump per-axis statistics table (for Fig 4D supplementary stats table)
write.csv(all_axes %>% dplyr::transmute(axis=as.character(signature), class=as.character(axis_class),
            delta_z=round(delta,3), ci_lo=round(ci_lo,3), ci_hi=round(ci_hi,3),
            nominal_p=signif(p_value,3), bh_q=signif(q_value,3)),
          file.path(OUT,"fig4D_axis_stats.csv"), row.names=FALSE)
p_f <- ggplot(all_axes, aes(delta, axis_label)) +
  annotate("rect", xmin=-Inf, xmax=0,   ymin=-Inf, ymax=Inf, fill=col_dn, alpha=0.06) +
  annotate("rect", xmin=0,    xmax=Inf, ymin=-Inf, ymax=Inf, fill=col_up, alpha=0.06) +
  geom_vline(xintercept=0, linewidth=0.35, colour="grey40") +
  geom_errorbar(aes(xmin=ci_lo, xmax=ci_hi), width=0, linewidth=0.45, colour="grey30", orientation="y") +
  geom_point(aes(fill=delta, size=pt_size), shape=21, stroke=0.28, colour="black") +
  geom_text(aes(label=val_lab), position=position_nudge(y=0.33), size=(base_font_size-1.7)/.pt, colour="grey25") +
  geom_text(aes(label=p_lab), position=position_nudge(y=-0.34), size=(base_font_size-2.7)/.pt, colour="grey55") +
  facet_grid(axis_class~., scales="free_y", space="free_y", switch="y") +
  scale_fill_gradient2(low=col_dn, mid=col_mid, high=col_up, midpoint=0, guide="none", limits=c(-1.5,1.5), oob=squish) +
  scale_size_identity() +
  scale_x_continuous(expand=expansion(mult=c(0.18,0.20))) +
  scale_y_discrete(labels=md_axis_lab) +
  labs(x="Δ group-mean z-score (FAP-KO − WT)", y=NULL, title="Mechanism axes align directionally with FAP loss",
       subtitle="<span style='color:#3b6ea5'>← reduced in FAP-KO</span>   ·   <span style='color:#b0463b'>increased →</span>    (n=4/grp; two-sided Wilcoxon; concordant direction — no axis survives BH-FDR, all q>0.05)", tag="f") +
  theme_nc() + theme(strip.text.y.left=element_text(angle=0, face="bold", size=base_font_size-1.2), strip.placement="outside",
                     strip.background=element_rect(fill="grey95", colour=NA),
                     axis.text.y.left=ggtext::element_markdown(lineheight=1.0, halign=1),
                     plot.subtitle=ggtext::element_markdown(size=base_font_size-1.5),
                     panel.spacing.y=unit(1.4,"mm"), panel.grid.major.y=element_line(colour="grey92", linewidth=0.3),
                     plot.title=element_text(margin=margin(l=2)))
save_panel(p_f, "panel_f_mechanism_axes_forest", 106, 76)

## =======================================================================
## g — ECM/focal-adhesion protein heatmap (ITGA5 ★)  ★ TOP MAIN CANDIDATE
## =======================================================================
ecm_markers <- c("Fap","Fn1","Postn","Tnc","Col1a1","Col3a1","Col6a1","Col12a1","Cthrc1","Thbs2","Mmp2","Mmp14","Itga2","Itgb1","Itga5","Lox")
ecm_prot <- bind_rows(lapply(ecm_markers, prot_long)) %>%
  mutate(value_log=log2(value+1), gene_display=toupper(gene), sample=factor(sample, levels=samp_order),
         sample_label=factor(samp_labs[as.character(sample)], levels=unname(samp_labs))) %>% filter(!is.na(value_log))
ecm_heat <- ecm_prot %>% group_by(gene_display) %>% mutate(z=safe_z(value_log)) %>% ungroup()
ecm_stats <- ecm_heat %>% group_by(gene_display) %>%
  summarise(delta_z=mean(z[group=="FAP-deficient host"],na.rm=TRUE)-mean(z[group=="WT host"],na.rm=TRUE), .groups="drop") %>% arrange(delta_z)
ecm_gene_levels <- ecm_stats$gene_display
ecm_gene_md <- dplyr::case_when(
  ecm_gene_levels=="ITGA2" ~ "*ITGA2*<span style='font-size:5pt;color:#8f8f8f'> (Fig 3)</span>",
  ecm_gene_levels=="ITGB1" ~ "*ITGB1*<span style='font-size:5pt;color:#8f8f8f'> (Fig 3)</span>",
  TRUE ~ paste0("*", ecm_gene_levels, "*"))
ecm_gene_lookup <- setNames(ecm_gene_md, ecm_gene_levels)
ecm_heat <- ecm_heat %>% mutate(gene_display_md=factor(ecm_gene_lookup[as.character(gene_display)], levels=rev(ecm_gene_md)))
p_g <- ggplot(ecm_heat, aes(sample_label, gene_display_md, fill=z)) +
  geom_tile(colour="white", linewidth=0.28) +
  facet_grid(.~group, scales="free_x", space="free_x") +
  scale_fill_gradient2(low=col_dn, mid=col_mid, high=col_up, midpoint=0, limits=c(-2,2), oob=squish, name="protein\nz-score") +
  labs(x=NULL, y=NULL, title="ECM/focal-adhesion protein heatmap",
       subtitle="DIA protein, z per gene · collagen/matrix ligands trend ↓ (n.s., n=4); receptors unchanged. α2β1 nominated from human spatial (Fig 3), n.s. here", tag="g") +
  theme_nc() + theme(axis.line=element_blank(), axis.ticks=element_blank(),
                     axis.text.x=element_text(size=base_font_size-1.2), axis.text.y=ggtext::element_markdown(size=base_font_size-0.7),
                     strip.text=element_text(face="bold", size=base_font_size-0.8), panel.spacing.x=unit(1.5,"mm"),
                     legend.key.height=unit(6,"mm"), legend.key.width=unit(2.5,"mm"))
save_panel(p_g, "panel_g_ecm_focal_adhesion_heatmap", 120, 82)

## transposed wide-short version of g for the one-row layout (genes on x, samples on y)
ecm_heat_row <- ecm_heat %>% mutate(gene_x_md = factor(ecm_gene_lookup[as.character(gene_display)], levels = ecm_gene_md))
p_g_row <- ggplot(ecm_heat_row, aes(gene_x_md, sample_label, fill=z)) +
  geom_tile(colour="white", linewidth=0.28) +
  facet_grid(group ~ ., scales="free_y", space="free_y", switch="y") +
  scale_fill_gradient2(low=col_dn, mid=col_mid, high=col_up, midpoint=0, limits=c(-2,2), oob=squish, name="protein z") +
  labs(x=NULL, y=NULL, title="ECM/focal-adhesion proteins · WT vs FAP-KO host",
       subtitle="DIA protein, z per gene · ligands trend ↓ (n.s.), receptors unchanged; α2β1 nominated from Fig 3", tag="g") +
  theme_nc() + theme(axis.line=element_blank(), axis.ticks=element_blank(),
    axis.text.x=ggtext::element_markdown(angle=90, hjust=1, vjust=0.5, size=base_font_size-1.6),
    axis.text.y=element_text(size=base_font_size-1.6),
    strip.placement="outside", strip.text.y.left=element_text(angle=0, face="bold", size=base_font_size-1.0),
    panel.spacing.y=unit(1,"mm"), legend.position="right",
    legend.key.height=unit(4,"mm"), legend.key.width=unit(2,"mm"),
    plot.subtitle=element_text(size=base_font_size-2.0))
save_panel(p_g_row, "panel_g_ecm_heatmap_TRANSPOSED", 120, 46)

## =======================================================================
## h — thyroid-lineage RNA readouts rise  ★ MAIN CANDIDATE (RNA-level, honest)
## =======================================================================
lineage_targets <- tibble(gene_mgi=c("Pax8","Tpo","Slc5a5"), display=c("PAX8","TPO","SLC5A5/NIS"))
lineage_rna <- bind_rows(lapply(seq_len(nrow(lineage_targets)), function(i)
  fpkm_long(lineage_targets$gene_mgi[i]) %>% mutate(gene_display=lineage_targets$display[i]))) %>%
  mutate(value_log=log2(value+1), readout=paste0(gene_display," RNA"),
         sample=factor(sample, levels=samp_order), sample_label=factor(samp_labs[as.character(sample)], levels=unname(samp_labs))) %>% filter(!is.na(value_log))
readout_levels_h <- c("PAX8 RNA","TPO RNA","SLC5A5/NIS RNA")
lineage_heat <- lineage_rna %>% filter(readout %in% readout_levels_h) %>% group_by(readout) %>% mutate(z=safe_z(value_log)) %>% ungroup() %>%
  mutate(readout=factor(readout, levels=rev(readout_levels_h)))
lineage_stats <- lineage_heat %>% group_by(readout) %>%
  summarise(delta_z=mean(z[group=="FAP-deficient host"],na.rm=TRUE)-mean(z[group=="WT host"],na.rm=TRUE),
            p_value=tryCatch(wilcox.test(value_log[group=="WT host"], value_log[group=="FAP-deficient host"], exact=FALSE)$p.value, error=function(e) NA_real_), .groups="drop") %>%
  mutate(sig_tag=case_when(p_value<0.05~"*", p_value<0.1~".", TRUE~""), readout=factor(readout, levels=rev(readout_levels_h)))
h_heat <- ggplot(lineage_heat, aes(sample_label, readout, fill=z)) +
  geom_tile(colour="white", linewidth=0.28) + facet_grid(.~group, scales="free_x", space="free_x") +
  scale_fill_gradient2(low=col_dn, mid=col_mid, high=col_up, midpoint=0, limits=c(-2,2), oob=squish, name="z-score") +
  labs(x=NULL, y=NULL) + theme_nc() +
  theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text.x=element_text(size=base_font_size-1.4),
        axis.text.y=element_text(face="italic", size=base_font_size-0.8), strip.text=element_text(face="bold", size=base_font_size-0.8),
        panel.spacing.x=unit(1.3,"mm"), legend.key.height=unit(6,"mm"), legend.key.width=unit(2.5,"mm"))
h_delta <- ggplot(lineage_stats, aes(delta_z, readout)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey45") +
  geom_segment(aes(x=0, xend=delta_z, yend=readout), linewidth=0.5, colour="grey45") +
  geom_point(aes(fill=delta_z), shape=21, size=2.2, stroke=0.25, colour="black") +
  geom_text(aes(label=sig_tag), nudge_x=0.14, size=(base_font_size+0.5)/.pt) +
  scale_fill_gradient2(low=col_dn, mid=col_mid, high=col_up, midpoint=0, limits=c(-1.5,1.5), oob=squish, guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0.18,0.25))) +
  labs(x="KO - WT (delta z)", y=NULL) + theme_nc() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_line(colour="grey92", linewidth=0.25))
h_heat <- h_heat + theme(strip.text=element_text(face="bold", size=base_font_size-1.4))
p_h <- (h_heat | h_delta) + plot_layout(widths=c(1.35,0.75)) +
  plot_annotation(title="Thyroid-lineage RNA readouts rise after stromal FAP loss",
                  subtitle="RNA-level triad: PAX8, TPO, SLC5A5/NIS",
                  theme=theme(plot.title=element_text(size=base_font_size+1, face="plain", family="Helvetica", margin=margin(l=16, b=1)),
                              plot.subtitle=element_text(size=base_font_size-1, colour="grey35", family="Helvetica", margin=margin(l=16)),
                              plot.margin=margin(3,5,3,4)))
p_h <- wrap_elements(full=p_h) + labs(tag="h") +
  theme(plot.tag.position=c(0.01, 1.0))
save_panel(p_h, "panel_h_lineage_rna_readouts", 112, 56)

## =======================================================================
## COMBINED preview — user-selected main-Fig-4 bioinformatics block: b f g h
## Causal reading: b on-target -> f axes shift // g ECM-protein(ITGA5) -> h lineage
## =======================================================================
## ---- beautified ONE-ROW main block: b | f | g(transposed) | h(delta-only) ----
## row-friendly lineage h: delta summary only (per-sample heatmap -> standalone/supp)
hr <- lineage_stats %>% mutate(lbl = factor(sub(" RNA","",as.character(readout)),
                               levels = sub(" RNA","",levels(readout))))
p_h_row <- ggplot(hr, aes(delta_z, lbl)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey45") +
  geom_segment(aes(x=0, xend=delta_z, yend=lbl), linewidth=0.6, colour="grey45") +
  geom_point(aes(fill=delta_z), shape=21, size=2.8, stroke=0.25, colour="black") +
  geom_text(aes(label=sig_tag), nudge_x=0.10, size=(base_font_size+0.5)/.pt) +
  scale_fill_gradient2(low=col_dn, mid=col_mid, high=col_up, midpoint=0, limits=c(-1.5,1.5), oob=squish, guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0.12,0.32))) +
  labs(x="KO − WT (Δz)", y=NULL, title="Lineage rebound", subtitle="PAX8 / TPO / NIS RNA", tag="h") +
  theme_nc() + theme(axis.text.y=element_text(face="italic", size=base_font_size-0.4),
                     panel.grid.major.y=element_line(colour="grey92", linewidth=0.25))

## trimmed single-line titles for a clean row
b_r <- p_b     + labs(title="On-target FAP loss", subtitle=NULL, tag="a")
f_r <- p_f     + labs(title="Mechanism axes align directionally", subtitle="n=4 · nominal P (uncorrected) · no axis FDR-sig", tag="b")
g_r <- p_g_row + labs(title="ECM/FA proteins (DIA)", subtitle="ligands trend ↓ (n.s.); receptors unchanged", tag="c")

## standalone compact lineage panel (delta-only) — pairs with the real NIS/TTF1 IHC, NOT in the bio row
grDevices::quartz(file=file.path(OUT,"panel_h_lineage_delta_compact.pdf"), type="pdf",
                  width=64/25.4, height=42/25.4, family="Helvetica")
print(p_h_row + labs(tag="h")); dev.off()
ggsave(file.path(OUT,"panel_h_lineage_delta_compact.png"), p_h_row + labs(tag="h"),
       width=64, height=42, units="mm", dpi=400, device=ragg::agg_png)

## MAIN bioinformatics row = b + f + g  (3 panels, comfortable at the true 180 mm column)
fig4_bio <- (b_r | f_r | g_r) + plot_layout(widths = c(0.72, 1.26, 1.72)) +
  plot_annotation(theme = theme(plot.margin = margin(3, 9, 3, 3)))
Wc <- 180; Hc <- 60
grDevices::quartz(file = file.path(OUT, "COMBINED_b_f_g_MAINROW_180mm.pdf"),
                  type = "pdf", width = Wc/25.4, height = Hc/25.4, family = "Helvetica")
print(fig4_bio); dev.off()
ggsave(file.path(OUT, "COMBINED_b_f_g_MAINROW_180mm.png"), fig4_bio,
       width = Wc, height = Hc, units = "mm", dpi = 500, device = ragg::agg_png)
cat("  wrote COMBINED_b_f_g_MAINROW_180mm (", Wc, "x", Hc, "mm) + panel_h_lineage_delta_compact\n")

## ---- USER-PREFERRED layout: a on-target | b honest-forest | c TALL ECM heatmap (C/D/E style) ----
a_t <- b_r + labs(title=NULL, subtitle=NULL, tag="C") + theme(plot.margin = margin(2, 11, 2, 2))
b_t <- f_r + labs(title=NULL, subtitle=NULL, tag="D") + theme(plot.margin = margin(2, 11, 2, 2))
g_tall <- p_g + labs(title=NULL,
                     subtitle=NULL,
                     tag="E") +
  facet_grid(. ~ group, scales="free_x", space="free_x",
             labeller=as_labeller(c("WT host"="WT","FAP-deficient host"="FAP-KO"))) +
  theme(plot.subtitle=element_text(size=base_font_size-1.8, colour="grey35"),
        strip.text=element_text(face="bold", size=base_font_size-0.8),
        plot.margin=margin(2,2,2,4))
fig4_bio_tallE <- (a_t | b_t | g_tall) + plot_layout(widths = c(0.56, 1.02, 1.42))
Wc2 <- 230; Hc2 <- 64
grDevices::quartz(file = file.path(OUT, "COMBINED_abc_tallE_230mm.pdf"),
                  type = "pdf", width = Wc2/25.4, height = Hc2/25.4, family = "Helvetica")
print(fig4_bio_tallE); dev.off()
ggsave(file.path(OUT, "COMBINED_abc_tallE_230mm.png"), fig4_bio_tallE,
       width = Wc2, height = Hc2, units = "mm", dpi = 500, device = ragg::agg_png)
cat("  wrote COMBINED_abc_tallE_230mm (", Wc2, "x", Hc2, "mm) — user-preferred tall-E layout\n")

cat("\nDONE. Panels in:\n  ", OUT, "\n")
