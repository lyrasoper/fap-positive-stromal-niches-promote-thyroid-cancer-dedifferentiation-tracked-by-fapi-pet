#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S18 - layered RNA/protein multi-omic support.
## Supplementary figure supporting the MINIMAL Fig 4 (causal-chain version)
## Transparent multi-omics backing: a RNA×protein concordance · b protein volcano
## c layer-resolved module shifts (RNA vs protein — justifies main panel e layer choice)
## d extended mechanism-gene heatmap. paper_A house style.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(patchwork); library(ggrepel); library(ragg)
})
proj <- Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>")
out  <- file.path(proj, "outputs/wt_vs_fapko_bulk_omics_20260407")
FIGD <- file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/fig4_minimal")
dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)
col_wt<-"#3B6B9C"; col_ko<-"#C9493A"; grey<-"#9AA0A6"; up_c<-"#C0392B"; dn_c<-"#2C6FA6"
rna_col<-"#4E79A7"; pro_col<-"#E1812C"; lin_c<-"#2E8B57"
pa_theme <- function(base=7) theme_classic(base_size=base, base_family="Helvetica") +
  theme(axis.text=element_text(size=base-1,colour="black"), axis.title=element_text(size=base),
        plot.title=element_text(size=base+1,face="bold",hjust=0), plot.subtitle=element_text(size=base-1.3,colour="grey35"),
        legend.text=element_text(size=base-1.5), legend.title=element_text(size=base-1),
        axis.line=element_line(linewidth=0.4), axis.ticks=element_line(linewidth=0.4),
        legend.key.size=unit(3,"mm"))

deg_rna <- read_tsv(file.path(out,"rna/rna_differential_expression.tsv"), show_col_types=FALSE)
deg_prot<- read_tsv(file.path(out,"proteomics/proteomics_differential_expression.tsv"), show_col_types=FALSE)
joined  <- read_tsv(file.path(out,"rna_protein_concordance/rna_prot_joined_all.tsv"), show_col_types=FALSE)
naba    <- read_tsv(file.path(out,"naba_caf_signatures/signature_scores_per_sample.tsv"), show_col_types=FALSE)
prog    <- read_tsv(file.path(out,"decoupler_crossomics/progeny_scores_per_sample.tsv"), show_col_types=FALSE)
dediff  <- read_tsv(file.path(out,"dediff_axes/integrated/dediff_axes_scores_per_sample.tsv"), show_col_types=FALSE)

mech <- c("Fap","Cthrc1","Col1a1","Col3a1","Fn1","Postn","Thbs2","Tnc",
          "Itga2","Itgb1","Itga5","Pax8","Tg","Tpo","Slc5a5","Duox2")

## ---- a: RNA × protein DIRECTIONAL concordance (honest: sign-level, not Pearson) ----
cc <- joined %>% filter(is.finite(rna_log2fc), is.finite(prot_log2fc))
ft <- fisher.test(table(cc$rna_log2fc>0, cc$prot_log2fc>0))
conc_all  <- mean(sign(cc$rna_log2fc)==sign(cc$prot_log2fc))
resp <- cc %>% filter(either_sig==TRUE); conc_resp <- mean(sign(resp$rna_log2fc)==sign(resp$prot_log2fc))
cc <- cc %>% mutate(grp=ifelse(either_sig,"responsive","unchanged"),
                    conc=ifelse(sign(rna_log2fc)==sign(prot_log2fc),"same sign","opposite"),
                    mech=gene_mgi %in% mech)
p_a <- ggplot(cc, aes(rna_log2fc, prot_log2fc)) +
  geom_hline(yintercept=0,linewidth=0.3,colour="grey70")+geom_vline(xintercept=0,linewidth=0.3,colour="grey70")+
  geom_point(data=filter(cc,grp=="unchanged"), size=0.2, colour="grey85", alpha=0.4) +
  geom_point(data=filter(cc,grp=="responsive"), aes(colour=conc), size=0.7, alpha=0.85) +
  geom_point(data=filter(cc,mech), shape=21, fill=col_ko, colour="black", stroke=0.3, size=1.5) +
  geom_text_repel(data=filter(cc,mech), aes(label=gene_mgi), size=1.9, fontface="italic",
                  max.overlaps=Inf, segment.size=0.2, box.padding=0.2) +
  scale_colour_manual(values=c(`same sign`="#2E8B57", opposite="#C46A1F"), name=NULL) +
  labs(x="RNA log2 FC (KO−WT)", y="Protein log2 FC (KO−WT)",
       title="a  RNA × protein directional concordance",
       subtitle=sprintf("Fisher OR = %.1f, P = %.0e · %.0f%% same-sign (responsive %.0f%%)",
                        ft$estimate, ft$p.value, 100*conc_all, 100*conc_resp)) +
  coord_cartesian(xlim=c(-4,4), ylim=c(-4,4)) + pa_theme() + theme(legend.position=c(0.80,0.13))

## ---- b: protein DIA volcano ----
vp <- deg_prot %>% mutate(nlp=pmin(-log10(p_value),12),
        cls=case_when(toupper(gene) %in% toupper(c("Fap","Cthrc1","Col1a1","Col3a1","Fn1","Postn","Thbs2","Tnc","Itga2","Itgb1","Itga5"))~"ECM/integrin",
                      toupper(gene) %in% toupper(c("Pax8","Tg","Tpo","Slc5a5","Duox2"))~"Lineage",
                      regulation=="up"~"up", regulation=="down"~"down", TRUE~"ns"))
p_b <- ggplot(vp, aes(log2_fc,nlp)) +
  geom_point(data=filter(vp,cls %in% c("ns","up","down")), aes(colour=cls), size=0.35, alpha=0.45) +
  geom_point(data=filter(vp,cls %in% c("ECM/integrin","Lineage")), aes(fill=cls), size=1.5, shape=21, colour="black", stroke=0.3) +
  geom_text_repel(data=filter(vp, toupper(gene) %in% toupper(mech)), aes(label=gene, colour=cls),
                  size=1.9, fontface="italic", max.overlaps=Inf, segment.size=0.2, box.padding=0.2) +
  scale_colour_manual(values=c(ns=grey,up="#E8B4AD",down="#A9C4DB",`ECM/integrin`=col_ko,Lineage=lin_c), guide="none") +
  scale_fill_manual(values=c(`ECM/integrin`=col_ko,Lineage=lin_c), name=NULL) +
  geom_vline(xintercept=0,linewidth=0.3,colour="grey60",linetype=2) +
  labs(x="protein log2 FC (KO−WT)", y="-log10 P", title="b  DIA proteomics perturbation landscape") +
  pa_theme() + theme(legend.position=c(0.84,0.92))

## ---- c: layer-resolved module shifts (RNA vs protein) ----
dl <- function(df,key,vcol) df %>% group_by(ds=dataset, k=.data[[key]]) %>%
  summarise(d=mean(.data[[vcol]][grepl("FAP",group)])-mean(.data[[vcol]][!grepl("FAP",group)]), .groups="drop")
mod <- bind_rows(
  dl(dediff,"signature","score_z") %>% filter(k %in% c("ECM-integrin-FAK","ERK score","YAP/TAZ","BRAF-RAS score","TDS / lineage identity")),
  dl(naba,"signature","score") %>% filter(k %in% c("Naba_CORE_MATRISOME","Kieffer_ecm_myCAF")),
  dl(prog,"pathway","score") %>% filter(k %in% c("MAPK","EGFR","NFkB"))) %>%
  mutate(k=recode(k, "Naba_CORE_MATRISOME"="Core matrisome","Kieffer_ecm_myCAF"="ecm-myCAF",
                  "TDS / lineage identity"="TDS / lineage","MAPK"="MAPK (PROGENy)"),
         layer=case_when(k %in% c("ECM-integrin-FAK","Core matrisome","ecm-myCAF")~"ECM / integrin",
                         k %in% c("MAPK (PROGENy)","ERK score","YAP/TAZ","EGFR","NFkB")~"Signalling",
                         TRUE~"Lineage / state"),
         layer=factor(layer,levels=c("ECM / integrin","Signalling","Lineage / state")),
         ds=factor(ds,levels=c("RNA","Proteomics")))
ord <- mod %>% group_by(k) %>% summarise(m=mean(d)) %>% arrange(m) %>% pull(k)
mod <- mod %>% mutate(k=factor(k,levels=ord))
p_c <- ggplot(mod, aes(d,k,fill=ds)) +
  geom_col(position=position_dodge(0.7), width=0.62) +
  geom_vline(xintercept=0,linewidth=0.3,colour="grey50") +
  facet_grid(layer~., scales="free_y", space="free_y", switch="y") +
  scale_fill_manual(values=c(RNA=rna_col,Proteomics=pro_col), name=NULL) +
  labs(x="Module shift (FAP-KO − WT)", y=NULL, title="c  Layer-resolved module shifts",
       subtitle="ECM modules drop at protein layer; signalling & lineage shift at RNA (basis for main Fig. 4e)") +
  pa_theme() + theme(strip.placement="outside", strip.background=element_blank(),
                     strip.text.y.left=element_text(size=6,face="bold",angle=0),
                     axis.text.y=element_text(size=6), panel.spacing=unit(0.8,"mm"),
                     legend.position=c(0.85,0.93), plot.subtitle=element_text(size=5.6))

## ---- d: extended mechanism heatmap ----
groups <- tibble(gene=c("Fap","Cthrc1","Mfap5","Lox","Mmp2","Mmp14",
                        "Col1a1","Col3a1","Col6a1","Fn1","Postn","Thbs2","Tnc","Vcan",
                        "Itga2","Itgb1","Itga5","Cd44","Ddr1",
                        "Tg","Tpo","Pax8","Slc5a5","Duox2","Nkx2-1","Foxe1"),
  grp=c(rep("FAP/CAF & remodeling",6), rep("ECM ligands",8), rep("Integrin/receptor",5), rep("Lineage",7)))
rna_v <- deg_rna %>% transmute(gene, RNA=log2_fc)
prot_v<- deg_prot %>% transmute(gU=toupper(gene), Protein=log2_fc) %>% distinct(gU,.keep_all=TRUE)
hd <- groups %>% mutate(gU=toupper(gene)) %>% left_join(rna_v,by="gene") %>% left_join(prot_v,by="gU") %>% select(-gU) %>%
  pivot_longer(c(RNA,Protein),names_to="modality",values_to="lfc") %>%
  mutate(grp=factor(grp,levels=c("FAP/CAF & remodeling","ECM ligands","Integrin/receptor","Lineage")),
         gene=factor(gene,levels=rev(groups$gene)), modality=factor(modality,levels=c("RNA","Protein")),
         lfc_cl=pmax(pmin(lfc,3),-3), nd=is.na(lfc))
p_d <- ggplot(hd, aes(modality,gene,fill=lfc_cl)) +
  geom_tile(colour="white",linewidth=0.5) +
  geom_text(data=filter(hd,nd), aes(label="ND"), size=1.7, colour="grey55") +
  geom_text(data=filter(hd,!nd & gene=="Itga2"), aes(label="★"), size=2, colour="black") +
  facet_grid(grp~., scales="free_y", space="free_y", switch="y") +
  scale_fill_gradient2(low=dn_c, mid="#F5F5F5", high=up_c, midpoint=0, na.value="white", name="log2 FC", limits=c(-3,3)) +
  labs(x=NULL, y=NULL, title="d  Extended mechanism panel", subtitle="RNA + DIA protein layers") +
  pa_theme() + theme(axis.text.y=element_text(size=5.4,face="italic"), strip.placement="outside",
                     strip.background=element_blank(), strip.text.y.left=element_text(size=5.4,face="bold",angle=0),
                     panel.spacing=unit(0.6,"mm"), legend.key.height=unit(4,"mm"), legend.key.width=unit(2.5,"mm"))

final <- (p_a | p_b) / (p_c | p_d) + plot_layout(heights=c(1,1.15))
agg_png(file.path(FIGD,"Fig4_supp_layered_omics.png"), width=180, height=170, units="mm", res=360)
print(final); dev.off()
cat("DONE Fisher OR=", round(ft$estimate,2), " conc_all=", round(conc_all,3), "\n")
