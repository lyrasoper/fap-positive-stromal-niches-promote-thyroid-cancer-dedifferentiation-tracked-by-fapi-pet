#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Figure 4C-E - causal-chain layout variant of the same panels.
## Fig 4 MINIMAL — causal chain version
## FAP loss → ECM/integrin ↓ → MAPK/ERK ↓ → thyroid-lineage ↑
## 6 panels: A design+PCA · B FAP↓ + CAF stable · C volcano · D focused heatmap
##           · E module shifts · F concise model. paper_A house style.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(ragg)
})
proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
out  <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
FIGD <- file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/fig4_minimal")
dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)

## ---- palette (semantic) ----
col_wt <- "#3B6B9C"; col_ko <- "#C9493A"; grey <- "#9AA0A6"
up_c <- "#C0392B"; dn_c <- "#2C6FA6"; ecm_c <- "#C9493A"; int_c <- "#E08214"; lin_c <- "#2E8B57"
## causal-chain stage palette — used consistently in a, e, f (read colour once, applies figure-wide)
st_caf<-"#A4303F"; st_ecm<-"#D06B52"; st_int<-"#E2A24A"; st_mapk<-"#CE7B3C"; st_lin<-"#4E8C72"
card <- function(x,y,w,h,label,fill,tcol="white",sz=2.15,r=0.22) list(
  annotation_custom(grid::roundrectGrob(r=grid::unit(r,"snpc"), gp=grid::gpar(fill=fill,col=NA)),
                    xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2),
  annotate("text",x=x,y=y,label=label,colour=tcol,size=sz,fontface="bold",lineheight=0.82))
flow <- function(x1,x2,y,lab=NULL,col="grey45"){
  g<-list(annotate("segment",x=x1,xend=x2,y=y,yend=y,
                   arrow=arrow(length=grid::unit(1.2,"mm"),type="closed"),linewidth=0.45,colour=col))
  if(!is.null(lab)) g<-c(g,list(annotate("text",x=(x1+x2)/2,y=y+0.2,label=lab,size=1.55,colour="grey45",fontface="italic")))
  g}
pa_theme <- function(base=7){
  theme_classic(base_size=base, base_family="Helvetica") +
  theme(axis.text=element_text(size=base-1, colour="black"),
        axis.title=element_text(size=base),
        plot.title=element_text(size=base+1, face="bold", hjust=0),
        plot.subtitle=element_text(size=base-1, colour="grey35"),
        legend.text=element_text(size=base-1.5), legend.title=element_text(size=base-1),
        axis.line=element_line(linewidth=0.4), axis.ticks=element_line(linewidth=0.4),
        legend.key.size=unit(3,"mm"), plot.tag=element_text(size=base+3, face="bold"))
}
tag_lab <- function(p, tag) p + labs(tag=tag)

## ============ DATA ============
deg_rna  <- read_tsv(file.path(out,"rna/rna_differential_expression.tsv"), show_col_types=FALSE)
deg_prot <- read_tsv(file.path(out,"proteomics/proteomics_differential_expression.tsv"), show_col_types=FALSE)
joined   <- read_tsv(file.path(out,"rna_protein_concordance/rna_prot_joined_all.tsv"), show_col_types=FALSE)
nnls     <- read_tsv(file.path(out,"caf_deconvolution/nnls_cell_fractions.tsv"), show_col_types=FALSE)
dstat    <- read_tsv(file.path(out,"caf_deconvolution/deconvolution_group_stats.tsv"), show_col_types=FALSE)
naba     <- read_tsv(file.path(out,"naba_caf_signatures/signature_scores_per_sample.tsv"), show_col_types=FALSE)
prog     <- read_tsv(file.path(out,"decoupler_crossomics/progeny_scores_per_sample.tsv"), show_col_types=FALSE)
dediff   <- read_tsv(file.path(out,"dediff_axes/integrated/dediff_axes_scores_per_sample.tsv"), show_col_types=FALSE)
fpkm <- read_tsv(file.path(proj,"001_mRNA_Summary/1.GeneExpression/1_genes_fpkm_expression.txt"), show_col_types=FALSE)
prot <- read_tsv(file.path(proj,"002_DIA_Summary/02.ProteinExp/protein_annotation_profile.txt"), show_col_types=FALSE)

## ============ PANEL A — design + PCA ============
# schematic (minimal)
p_schem <- ggplot() +
  card(1.0,2.0,1.8,0.74,"WT host\n(FAP+/+)",col_wt,sz=2.05) +
  card(1.0,0.98,1.8,0.74,"FAP-KO host\n(FAP+/-)",col_ko,sz=2.05) +
  flow(1.97,2.75,2.0) + flow(1.97,2.75,0.98) +
  card(3.55,1.5,1.45,0.72,"BPC\nxenograft","#5B636E",sz=2.0) +
  flow(4.35,4.95,1.5) +
  card(5.95,1.5,1.75,0.72,"RNA-seq +\nDIA proteomics","#343A40",sz=2.0) +
  coord_cartesian(xlim=c(0.0,6.9), ylim=c(0.55,2.5), clip="off") +
  labs(title="a  Stromal FAP perturbation → paired omics") + theme_void(base_family="Helvetica") +
  theme(plot.title=element_text(size=7.6, face="bold", hjust=0), plot.margin=margin(2,2,2,2))

# PCA from FPKM
fcols <- grep("^FPKM\\.", names(fpkm), value=TRUE)
mat <- as.matrix(fpkm[, fcols]); rownames(mat) <- make.unique(fpkm$gene_name)
mat <- log2(mat+1); v <- apply(mat,1,var); mat <- mat[order(-v)[1:2000],]
pc <- prcomp(t(mat), center=TRUE, scale.=TRUE)
ve <- round(100*pc$sdev^2/sum(pc$sdev^2),1)
pcd <- as_tibble(pc$x[,1:2], rownames="sample") %>%
  mutate(group=ifelse(grepl("FAP",sample),"FAP-KO","WT"))
p_pca <- ggplot(pcd, aes(PC1,PC2,fill=group)) +
  geom_point(size=2.2, shape=21, colour="white", stroke=0.4) +
  scale_fill_manual(values=c(WT=col_wt,`FAP-KO`=col_ko), name=NULL) +
  labs(x=sprintf("PC1 (%.0f%%)",ve[1]), y=sprintf("PC2 (%.0f%%)",ve[2]),
       title="PCA · top-2,000 variable genes") +
  pa_theme() + theme(legend.position=c(0.82,0.18), plot.title=element_text(size=7.5,face="plain"))

## ============ PANEL B — FAP on-target + CAF composition stable ============
# Fap per-sample RNA (protein all-NA in raw DIA → show as log2FC annotation)
fap_rna <- tibble(sample=fcols, val=log2(as.numeric(fpkm[fpkm$gene_name=="Fap", fcols][1,])+1)) %>%
  mutate(group=ifelse(grepl("FAP",sample),"FAP-KO","WT"))
fap_rna_lfc  <- deg_rna  %>% filter(toupper(gene)=="FAP") %>% pull(log2_fc) %>% .[1]
fap_prot_lfc <- deg_prot %>% filter(toupper(gene)=="FAP") %>% pull(log2_fc) %>% .[1]
ymx <- max(fap_rna$val)
p_fap <- ggplot(fap_rna, aes(group,val)) +
  geom_boxplot(aes(fill=group), width=0.46, outlier.shape=NA, linewidth=0.35, alpha=0.9) +
  geom_jitter(aes(fill=group), width=0.08, height=0, size=1.7, shape=21, colour="white", stroke=0.45) +
  annotate("segment", x=1, xend=2, y=ymx+0.17, yend=ymx+0.17, linewidth=0.4, colour="grey30") +
  annotate("segment", x=1, xend=1, y=ymx+0.10, yend=ymx+0.17, linewidth=0.4, colour="grey30") +
  annotate("segment", x=2, xend=2, y=ymx+0.10, yend=ymx+0.17, linewidth=0.4, colour="grey30") +
  annotate("text", x=1.5, y=ymx+0.30, label="***", size=2.9, colour="grey20") +
  scale_fill_manual(values=c(WT=col_wt,`FAP-KO`=col_ko), guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0.04,0.16))) +
  labs(x=NULL, y="Fap (log2 FPKM)", title="b  FAP on-target loss",
       subtitle=sprintf("RNA log2FC %.1f (q=2e-7) · protein %.1f", fap_rna_lfc, fap_prot_lfc)) +
  pa_theme() + theme(axis.text.x=element_text(size=6.5), plot.subtitle=element_text(size=5.4))

# NNLS CAF fractions stable (clean labels + per-pair n.s.)
ct_lab <- c(myCAF="myCAF", iCAF="iCAF", apCAF="apCAF", Thyroid_tumor="Tumour", Immune="Immune", Endothelial="Endoth.")
nn <- nnls %>% mutate(group=ifelse(grepl("FAP",sample),"FAP-KO","WT"),
                      celltype=factor(celltype, levels=unique(celltype)))
ns_lab <- nn %>% group_by(celltype) %>% summarise(top=max(fraction), .groups="drop") %>%
  left_join(dstat %>% filter(method=="NNLS fraction") %>% transmute(celltype,p=wilcox_p), by="celltype") %>%
  mutate(lab=ifelse(p<0.05, sprintf("%.3f",p), "n.s."), y=top+0.016)
p_nnls <- ggplot(nn, aes(celltype, fraction, fill=group)) +
  geom_boxplot(width=0.62, outlier.shape=NA, linewidth=0.3, position=position_dodge(0.7)) +
  geom_text(data=ns_lab, aes(celltype, y, label=lab), inherit.aes=FALSE, size=1.5, colour="grey55") +
  scale_fill_manual(values=c(WT=col_wt,`FAP-KO`=col_ko), name=NULL) +
  scale_x_discrete(labels=ct_lab) +
  scale_y_continuous(expand=expansion(mult=c(0.04,0.1))) +
  labs(x=NULL, y="NNLS fraction", title="CAF composition preserved",
       subtitle="all cell-type fractions n.s. (WT vs FAP-KO)") +
  pa_theme() + theme(axis.text.x=element_text(angle=30,hjust=1,size=6),
                     legend.position=c(0.87,0.9), legend.key.size=unit(2.6,"mm"),
                     plot.title=element_text(size=7.5,face="plain"), plot.subtitle=element_text(size=5.4))

## ============ PANEL C — RNA volcano ============
lab_genes <- c("Fap","Cthrc1","Col1a1","Col3a1","Fn1","Postn","Thbs2","Tnc","Mmp2",
               "Itga2","Itgb1","Itga5","Pax8","Tg","Tpo","Slc5a5","Duox2","Nkx2-1")
vol <- deg_rna %>% mutate(nlp=-log10(p_value),
              cls=case_when(gene %in% c("Fap","Cthrc1","Col1a1","Col3a1","Fn1","Postn","Thbs2","Tnc","Mmp2","Itga2","Itgb1","Itga5")~"ECM/integrin",
                            gene %in% c("Pax8","Tg","Tpo","Slc5a5","Duox2","Nkx2-1")~"Lineage",
                            regulation=="up"~"up", regulation=="down"~"down", TRUE~"ns")) %>%
  mutate(nlp=pmin(nlp,30))
p_vol <- ggplot(vol, aes(log2_fc,nlp)) +
  geom_point(data=filter(vol,cls %in% c("ns","up","down")), aes(colour=cls), size=0.35, alpha=0.45) +
  geom_point(data=filter(vol,cls %in% c("ECM/integrin","Lineage")), aes(fill=cls), size=1.6, shape=21, colour="black", stroke=0.3) +
  geom_text_repel(data=filter(vol, gene %in% lab_genes),
                  aes(label=gene, colour=cls), size=2, fontface="italic", max.overlaps=Inf,
                  box.padding=0.2, segment.size=0.2, min.segment.length=0) +
  scale_colour_manual(values=c(ns=grey,up="#E8B4AD",down="#A9C4DB",`ECM/integrin`=ecm_c,Lineage=lin_c), guide="none") +
  scale_fill_manual(values=c(`ECM/integrin`=ecm_c,Lineage=lin_c), name=NULL) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey60", linetype=2) +
  labs(x="log2 FC (FAP-KO vs WT)", y="-log10 P",
       title="c  ECM/stromal ↓ · lineage genes ↑ (FAP-KO)") +
  pa_theme() + theme(legend.position=c(0.85,0.88), plot.title=element_text(size=7.3))

## ============ PANEL D — focused mechanism heatmap ============
groups <- tibble(
  gene=c("Fap","Cthrc1","Mfap5",  "Col1a1","Col3a1","Fn1","Postn","Thbs2",
         "Itga2","Itgb1","Itga5",  "Tg","Tpo","Pax8","Slc5a5","Duox2"),
  grp=c(rep("FAP/CAF",3), rep("ECM ligands",5), rep("Integrin receivers",3), rep("Lineage",5)))
rna_v  <- deg_rna  %>% transmute(gene, RNA=log2_fc)
prot_v <- deg_prot %>% transmute(gU=toupper(gene), Protein=log2_fc) %>% distinct(gU, .keep_all=TRUE)
hd <- groups %>% mutate(gU=toupper(gene)) %>%
  left_join(rna_v, by="gene") %>% left_join(prot_v, by="gU") %>% select(-gU) %>%
  pivot_longer(c(RNA,Protein), names_to="modality", values_to="lfc") %>%
  mutate(grp=factor(grp, levels=c("FAP/CAF","ECM ligands","Integrin receivers","Lineage")),
         gene=factor(gene, levels=rev(groups$gene)),
         modality=factor(modality, levels=c("RNA","Protein")),
         lfc_cl=pmax(pmin(lfc,3),-3), nd=is.na(lfc),
         txt=ifelse(!is.na(lfc) & abs(lfc)>=0.6, sprintf("%+.1f",lfc), NA_character_),
         txtcol=ifelse(abs(pmax(pmin(lfc,3),-3))>1.6,"white","grey25"))
p_heat <- ggplot(hd, aes(modality, gene, fill=lfc_cl)) +
  geom_tile(colour="white", linewidth=0.5) +
  geom_text(data=filter(hd,nd), aes(label="ND"), size=1.8, colour="grey55") +
  geom_text(data=filter(hd, !is.na(txt) & gene!="Itga2"), aes(label=txt, colour=txtcol), size=1.4) +
  scale_colour_identity() +
  geom_text(data=filter(hd,!nd & gene=="Itga2"), aes(label="★"), size=2.1, colour="black") +
  facet_grid(grp~., scales="free_y", space="free_y", switch="y") +
  scale_fill_gradient2(low=dn_c, mid="#F5F5F5", high=up_c, midpoint=0, na.value="white",
                       name="log2 FC\n(KO-WT)", limits=c(-3,3)) +
  labs(x=NULL, y=NULL, title="d  Focused mechanism panel") +
  pa_theme() + theme(axis.text.y=element_text(size=5.6, face="italic"),
                     strip.placement="outside", strip.background=element_blank(),
                     strip.text.y.left=element_text(size=5.6, face="bold", angle=0),
                     panel.spacing=unit(0.6,"mm"), legend.key.height=unit(4,"mm"),
                     legend.key.width=unit(2.5,"mm"))

## ============ PANEL E — module shifts (causal chain) ============
# each module read at its informative layer (ECM: protein-rich; signalling/lineage: RNA)
ddelta <- function(df, key, vcol, ds) df %>% filter(dataset==ds) %>% group_by(k=.data[[key]]) %>%
  summarise(d=mean(.data[[vcol]][grepl("FAP",group)])-mean(.data[[vcol]][!grepl("FAP",group)]), .groups="drop")
ded_p <- ddelta(dediff,"signature","score_z","Proteomics"); ded_r <- ddelta(dediff,"signature","score_z","RNA")
nab_p <- ddelta(naba,"signature","score","Proteomics");     pro_r <- ddelta(prog,"pathway","score","RNA")
gv <- function(t,k) t$d[t$k==k]
mods <- tribble(~layer,~item,~d,
  "ECM / integrin","ECM-integrin-FAK", gv(ded_p,"ECM-integrin-FAK"),
  "ECM / integrin","Core matrisome",   gv(nab_p,"Naba_CORE_MATRISOME"),
  "ECM / integrin","ecm-myCAF",        gv(nab_p,"Kieffer_ecm_myCAF"),
  "Signalling","MAPK (PROGENy)",       gv(pro_r,"MAPK"),
  "Signalling","ERK score",            gv(ded_r,"ERK score"),
  "Signalling","YAP/TAZ",              gv(ded_r,"YAP/TAZ"),
  "Lineage","TDS / lineage",           gv(ded_r,"TDS / lineage identity")) %>%
  mutate(layer=factor(layer,levels=c("ECM / integrin","Signalling","Lineage")),
         item=factor(item,levels=rev(item)), dir=ifelse(d>0,"up","down"))
p_mod <- ggplot(mods, aes(d, item, fill=layer)) +
  geom_col(width=0.66) + geom_vline(xintercept=0, linewidth=0.3, colour="grey50") +
  facet_grid(layer~., scales="free_y", space="free_y", switch="y") +
  scale_fill_manual(values=c(`ECM / integrin`=st_ecm, Signalling=st_mapk, Lineage=st_lin), guide="none") +
  labs(x="Module shift (FAP-KO − WT)", y=NULL,
       title="e  ECM/integrin ↓ · MAPK ↓ · lineage ↑",
       subtitle="ECM at protein layer; signalling & lineage at RNA") +
  pa_theme() + theme(strip.placement="outside", strip.background=element_blank(),
                     strip.text.y.left=element_text(size=5.6, face="bold", angle=0),
                     axis.text.y=element_text(size=6), panel.spacing=unit(0.8,"mm"),
                     plot.title=element_text(size=7.2), plot.subtitle=element_text(size=5.3))

## ============ PANEL F — concise model ============
yy<-1.12; cw<-1.82; ch<-0.94
p_mod2 <- ggplot() +
  card(1.0,yy,cw,ch,"FAP⁺ CAF\nniche",st_caf) + flow(1.97,2.53,yy,"secrete") +
  card(3.5,yy,cw,ch,"ECM ligands\nCOL/FN1/POSTN",st_ecm) + flow(4.47,5.03,yy,"engage") +
  card(6.0,yy,cw,ch,"α2β1\nITGA2/ITGB1",st_int) + flow(6.97,7.53,yy,"activate") +
  card(8.5,yy,cw,ch,"ERK/MAPK\noutput",st_mapk) + flow(9.47,10.03,yy,"suppress") +
  card(11.0,yy,cw,ch,"thyroid lineage\nTG/TPO/PAX8",st_lin) +
  annotate("text", x=6.0, y=2.04, label="FAP-KO breaks the chain → thyroid-lineage identity restored",
           size=2.3, fontface="bold.italic", colour=st_lin) +
  annotate("segment", x=11.0, xend=1.0, y=0.4, yend=0.4, linewidth=0.55, colour=st_lin, linetype="dashed",
           arrow=arrow(length=unit(1.4,"mm"), type="closed")) +
  annotate("text", x=6.0, y=0.22, label="stromal FAP loss", size=1.7, colour=st_lin, fontface="italic") +
  coord_cartesian(xlim=c(0.05,12.0), ylim=c(0.05,2.3), clip="off") +
  labs(title="f  Working model — collagen → α2β1 → MAPK → lineage") + theme_void(base_family="Helvetica") +
  theme(plot.title=element_text(size=8, face="bold", hjust=0))

## ============ ASSEMBLE ============
rowA  <- wrap_plots(p_schem, p_pca, widths=c(1,1.05))
rowBC <- wrap_plots(p_fap, p_nnls, p_vol, widths=c(0.62,1.0,1.35))
rowDE <- wrap_plots(p_heat, p_mod, widths=c(1,1))
final <- wrap_plots(rowA, rowBC, rowDE, p_mod2, ncol=1, heights=c(0.82,1.05,1.28,0.40))

agg_png(file.path(FIGD,"Fig4_minimal_causal.png"), width=180, height=214, units="mm", res=360)
print(final); dev.off()
cat("DONE ->", file.path(FIGD,"Fig4_minimal_causal.png"), "\n")
