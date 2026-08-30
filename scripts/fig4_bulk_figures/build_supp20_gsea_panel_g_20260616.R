#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S20 panel g - two-layer GSEA.
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.
## Supp 20 panel g — two-layer GSEA (FAP-KO), rendered at 240 mm width to match
## FigS17 (image51, 5669 px @ 600 dpi) for clean vertical compositing. paper_A.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr); library(ggtext); library(tidyr); library(tibble); library(ggh4x) })
out_root <- file.path(Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>"), "outputs/wt_vs_fapko_bulk_omics_20260407")
OUT <- file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/_fig4_bioinformatics_panels_20260613")
base_font_size <- 8
theme_nc <- function(base=base_font_size) theme_classic(base_size=base, base_family="Helvetica") +
  theme(plot.title=element_text(face="plain", size=base+1, hjust=0), plot.subtitle=element_text(size=base-1.2, colour="grey35"),
        plot.tag=element_text(face="bold", size=base+4, hjust=0), axis.title=element_text(size=base),
        axis.text=element_text(size=base-0.5, colour="black"), axis.line=element_line(linewidth=0.3), axis.ticks=element_line(linewidth=0.3),
        legend.text=element_text(size=base-0.5), legend.title=element_text(size=base-0.5), panel.grid=element_blank(), plot.margin=margin(3,6,3,4))

gsea_rna  <- read_tsv(file.path(out_root,"rna/rna_all_gsea_pathways.tsv"), show_col_types=FALSE)
gsea_prot <- read_tsv(file.path(out_root,"proteomics/proteomics_all_gsea_pathways.tsv"), show_col_types=FALSE)
rna_l  <- gsea_rna  %>% transmute(pw=tolower(pathway), rna_nes=NES,  rna_fdr=fdr)  %>% distinct(pw,.keep_all=TRUE)
prot_l <- gsea_prot %>% transmute(pw=tolower(pathway), prot_nes=NES, prot_fdr=fdr) %>% distinct(pw,.keep_all=TRUE)
curated <- tribble(~pw,~label,~class,
  "ecm receptor interaction","ECM–receptor interaction","ECM / adhesion",
  "focal adhesion","Focal adhesion","ECM / adhesion",
  "collagen biosynthesis and modifying enzymes","Collagen biosynthesis/modifying","ECM / adhesion",
  "collagen formation","Collagen formation","ECM / adhesion",
  "extracellular matrix organization","ECM organization","ECM / adhesion",
  "integrin cell surface interactions","Integrin cell-surface interactions","ECM / adhesion",
  "thyroid cancer","Thyroid cancer","Thyroid / lineage",
  "thyroid hormone synthesis","Thyroid hormone synthesis","Thyroid / lineage",
  "autoimmune thyroid disease","Autoimmune thyroid disease","Thyroid / lineage",
  "mapk signaling pathway","MAPK signaling","MAPK")
gd <- curated %>% left_join(rna_l,by="pw") %>% left_join(prot_l,by="pw")
ord <- gd %>% arrange(class, dplyr::coalesce(rna_nes,prot_nes)) %>% pull(label)
gd  <- gd %>% mutate(label=factor(label,levels=ord), class=factor(class,levels=c("ECM / adhesion","Thyroid / lineage","MAPK")))
gl  <- gd %>% pivot_longer(c(rna_nes,prot_nes),names_to="layer",values_to="nes") %>%
  mutate(layer=ifelse(layer=="rna_nes","RNA","Protein"), fdr=ifelse(layer=="RNA",rna_fdr,prot_fdr), fdr_sig=!is.na(fdr)&fdr<0.05) %>% filter(!is.na(nes))
lay_cols <- c("RNA"="#0F7B6C","Protein"="#9C6114")
## alternating row-band background (per class facet) — zebra striping for readability
zebra <- gd %>% group_by(class) %>% arrange(label) %>% mutate(.row = dplyr::row_number()) %>%
  ungroup() %>% filter(.row %% 2L == 0L)
## per-group colour-coded facet strips (soft tints, distinct from teal/gold dots; matching dark bold text)
strip_y <- strip_themed(
  background_y = list(element_rect(fill="#D3DEEC", colour=NA),   # ECM / adhesion  -> soft blue
                      element_rect(fill="#E1D8EA", colour=NA),   # Thyroid / lineage -> soft violet
                      element_rect(fill="#F0DBD3", colour=NA)),  # MAPK            -> soft rose
  text_y = list(element_text(colour="#2F5980", face="bold", angle=0, size=base_font_size-0.5),
                element_text(colour="#574272", face="bold", angle=0, size=base_font_size-0.5),
                element_text(colour="#94483C", face="bold", angle=0, size=base_font_size-0.5)))
p_g <- ggplot(gl, aes(nes,label)) +
  geom_rect(data=zebra, aes(ymin=.row-0.5, ymax=.row+0.5), xmin=-Inf, xmax=Inf,
            fill="grey95", colour=NA, inherit.aes=FALSE) +
  geom_vline(xintercept=0, linewidth=0.35, colour="grey45") +
  geom_segment(data=gd %>% filter(!is.na(rna_nes)&!is.na(prot_nes)), aes(x=rna_nes,xend=prot_nes,y=label,yend=label),
               inherit.aes=FALSE, colour="grey75", linewidth=0.7, lineend="round") +
  geom_point(aes(fill=layer,shape=fdr_sig), size=3.2, stroke=0.35, colour="black") +
  scale_shape_manual(values=c(`TRUE`=23,`FALSE`=21), labels=c(`TRUE`="FDR<0.05",`FALSE`="n.s."), name=NULL) +
  scale_fill_manual(values=lay_cols, name="Layer", guide=guide_legend(override.aes=list(shape=21))) +
  facet_grid2(class~., scales="free_y", space="free_y", switch="y", strip=strip_y) +
  scale_x_continuous(expand=expansion(mult=c(0.04,0.06))) +
  labs(x="GSEA NES (FAP-KO vs WT)   ·   ← lower in FAP-KO        higher in FAP-KO →", y=NULL,
       title="Two-layer gene-set enrichment — compensatory ECM transcription versus protein-level matrix reduction",
       subtitle="RNA induces an ECM/adhesion + thyroid-lineage programme (only ECM–receptor reaches FDR<0.05); at the protein layer the structural collagen/ECM modules are *reduced*, not accumulated", tag="g") +
  theme_nc() + theme(strip.placement="outside",
                     plot.subtitle=ggtext::element_markdown(size=base_font_size-1.4, colour="grey35"), plot.title=element_text(size=base_font_size+0.5),
                     panel.spacing.y=unit(2.2,"mm"), panel.grid.major.y=element_blank(),
                     panel.grid.major.x=element_line(colour="grey95", linewidth=0.25), legend.position="right")
ggsave(file.path(OUT,"supp20_panel_g_gsea_240mm.png"), p_g, width=240, height=82, units="mm", dpi=600, device=ragg::agg_png)
cat("wrote supp20_panel_g_gsea_240mm.png\n")
