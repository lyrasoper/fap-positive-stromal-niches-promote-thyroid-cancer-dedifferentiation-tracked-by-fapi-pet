#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Figure 3 (ITGA2/ITGB1 collagen-receptor axis).
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

## Supplementary — the α2β1 (ITGA2/ITGB1) collagen-receptor axis (Fig 3 mechanism)
## a COL→ITGA2/ITGB1 top LR pairs · b receiver Δ (terminal vs lineage)
## c ITGA2 receiver activity spatial maps · d ITGA2 distance gradient from FAP-CAF niche
## ITGA5 (α5β1) NOT featured here — secondary, lives in the public-proteome breadth.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(ggplot2);library(patchwork);library(ragg);library(stringr)})
r14 <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/results")
FIGD<- file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/itga2_axis"); dir.create(FIGD,showWarnings=FALSE,recursive=TRUE)
ITGA2_c<-"#C0392B"; ITGB1_c<-"#E08214"; other_c<-"#9AA0A6"; term_c<-"#C0392B"; lin_c<-"#2E8B57"
pa <- function(base=7) theme_classic(base_size=base,base_family="Helvetica")+
  theme(axis.text=element_text(size=base-1,colour="black"),axis.title=element_text(size=base),
        plot.title=element_text(size=base+1,face="bold",hjust=0),plot.subtitle=element_text(size=base-1.3,colour="grey35"),
        legend.text=element_text(size=base-1.5),legend.title=element_text(size=base-1),
        axis.line=element_line(linewidth=0.4),axis.ticks=element_line(linewidth=0.4),legend.key.size=unit(3,"mm"))
rclass <- function(r) ifelse(grepl("ITGA2",r),"ITGA2",ifelse(grepl("ITGB1",r),"ITGB1",ifelse(grepl("ITGA5",r),"ITGA5","other")))

lr <- read_tsv(file.path(r14,"20260410_fap_fibro_epi_mechanism_prioritization/ranked_candidate_interactions.tsv"),show_col_types=FALSE)
ctx<- read_tsv(file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/fig7_redesign/augmented_context_v3.tsv"),show_col_types=FALSE)

## ---- a: top LR pairs (receptor-coloured) ----
la <- lr %>% slice_min(priority_rank, n=12) %>%
  mutate(pair=paste0(ligand," → ",receptor), rc=rclass(receptor),
         pair=factor(pair, levels=rev(pair)))
p_a <- ggplot(la, aes(overall_priority_score, pair, fill=rc)) +
  geom_col(width=0.72) +
  scale_fill_manual(values=c(ITGA2=ITGA2_c,ITGB1=ITGB1_c,ITGA5=other_c,other=other_c), name=NULL) +
  labs(x="LIANA priority score", y=NULL, title="a  Collagen → α2β1 leads the LR axis",
       subtitle="ranks 1–4 = COL/THBS2 → ITGA2; ITGB1 next; ITGA5 lower-ranked") +
  pa() + theme(axis.text.y=element_text(size=5.6), legend.position=c(0.84,0.32), plot.title=element_text(size=7.3),
               plot.subtitle=element_text(size=5.3))

## ---- b: receiver Δ (terminal vs lineage) per receptor ----
rb <- lr %>% group_by(receptor) %>% summarise(cell=first(cell_delta_terminal_vs_lineage),
                spat=first(spatial_delta_terminal_vs_lineage), .groups="drop") %>%
  mutate(rc=rclass(receptor)) %>% arrange(cell) %>% mutate(receptor=factor(receptor,levels=receptor)) %>%
  pivot_longer(c(cell,spat), names_to="layer", values_to="delta") %>%
  mutate(layer=recode(layer, cell="single-cell", spat="spatial"))
p_b <- ggplot(rb, aes(delta, receptor, fill=layer)) +
  geom_col(position=position_dodge(0.7), width=0.62) + geom_vline(xintercept=0,linewidth=0.3,colour="grey50")+
  scale_fill_manual(values=c(`single-cell`="#4E79A7", spatial="#F2A900"), name=NULL) +
  labs(x="receiver Δ (terminal − lineage)", y=NULL, title="b  Receptor enrichment (terminal-dediff)",
       subtitle="ITGA2 / ITGB1 most terminal-enriched") +
  pa() + theme(axis.text.y=element_text(size=6), legend.position=c(0.8,0.2), plot.subtitle=element_text(size=5.3))

## ---- c: ITGA2 receiver activity spatial maps ----
reps <- c("P12","P17","P26","P44","P83")
cc <- ctx %>% filter(sample %in% reps) %>% group_by(sample) %>%
  mutate(a2=pmin(itga2_receiver_program_activity, quantile(itga2_receiver_program_activity,0.99,na.rm=TRUE))) %>% ungroup() %>%
  mutate(sample=factor(sample, levels=reps))
p_c <- ggplot(cc, aes(x_plot, -y_plot, colour=a2)) +
  geom_point(size=0.22) + facet_wrap(~sample, nrow=1) +
  scale_colour_gradientn(colours=c("#F0F0F0","#FBD9C6","#E8845C","#C0392B","#7B241C"), name="ITGA2\nactivity") +
  coord_equal() + labs(title="c  ITGA2 receiver-program activity (representative Visium samples)") +
  theme_void(base_family="Helvetica") +
  theme(plot.title=element_text(size=7.5,face="bold",hjust=0), strip.text=element_text(size=6.5,face="bold"),
        legend.title=element_text(size=5.6), legend.text=element_text(size=5), legend.key.size=unit(3,"mm"),
        legend.position="right")

## ---- d: distance gradient (ITGA2 tracks the niche; opposite to lineage) ----
brks<-c(-Inf,1,2,3,5,Inf); labs5<-c("0–1","1–2","2–3","3–5",">5")
# z-score each variable PER SAMPLE (across all spots), then average by ring across samples
dd <- ctx %>% filter(is.finite(distance_spot_units)) %>% group_by(sample) %>%
  mutate(ITGA2=scale(itga2_receiver_program_activity)[,1],
         `terminal load`=scale(epi_deconv_load_terminal_dedifferentiated)[,1],
         `lineage load`=scale(epi_deconv_load_lineage_preserved_epithelial)[,1]) %>% ungroup() %>%
  mutate(ring=cut(distance_spot_units, breaks=brks, labels=labs5)) %>%
  pivot_longer(c(ITGA2,`terminal load`,`lineage load`), names_to="var", values_to="z") %>%
  group_by(sample,ring,var) %>% summarise(z=mean(z,na.rm=TRUE), .groups="drop") %>%
  group_by(ring,var) %>% summarise(m=mean(z,na.rm=TRUE), se=sd(z,na.rm=TRUE)/sqrt(n()), .groups="drop")
p_d <- ggplot(dd, aes(ring,m,colour=var,group=var)) +
  geom_hline(yintercept=0,linewidth=0.3,colour="grey75")+
  geom_ribbon(aes(ymin=m-1.96*se,ymax=m+1.96*se,fill=var),alpha=0.15,colour=NA)+
  geom_line(linewidth=0.7)+geom_point(size=1.3)+
  scale_colour_manual(values=c(ITGA2=ITGA2_c,`terminal load`="#7B241C",`lineage load`=lin_c),name=NULL)+
  scale_fill_manual(values=c(ITGA2=ITGA2_c,`terminal load`="#7B241C",`lineage load`=lin_c),guide="none")+
  labs(x="distance from nearest FAP-CAF seed (spot rings)", y="mean z-score",
       title="d  ITGA2 tracks the niche + dediff state",
       subtitle="ITGA2 + terminal high at niche, fall with distance; lineage-preserved rises")+
  pa()+theme(legend.position=c(0.8,0.85), plot.subtitle=element_text(size=5.3))

## ---- e: ITGA5 as the SECONDARY (α5β1, FN1-driven) receiver — also niche-tracking ----
de <- ctx %>% filter(is.finite(distance_spot_units)) %>% group_by(sample) %>%
  mutate(ITGA2=scale(itga2_receiver_program_activity)[,1],
         ITGA5=scale(itga5_receiver_program_activity)[,1]) %>% ungroup() %>%
  mutate(ring=cut(distance_spot_units, breaks=brks, labels=labs5)) %>%
  pivot_longer(c(ITGA2,ITGA5), names_to="var", values_to="z") %>%
  group_by(sample,ring,var) %>% summarise(z=mean(z,na.rm=TRUE), .groups="drop") %>%
  group_by(ring,var) %>% summarise(m=mean(z,na.rm=TRUE), se=sd(z,na.rm=TRUE)/sqrt(n()), .groups="drop")
p_e <- ggplot(de, aes(ring,m,colour=var,group=var)) +
  geom_hline(yintercept=0,linewidth=0.3,colour="grey75")+
  geom_ribbon(aes(ymin=m-1.96*se,ymax=m+1.96*se,fill=var),alpha=0.13,colour=NA)+
  geom_line(linewidth=0.7)+geom_point(size=1.3)+
  scale_colour_manual(values=c(ITGA2=ITGA2_c,ITGA5="#6A51A3"),name=NULL)+
  scale_fill_manual(values=c(ITGA2=ITGA2_c,ITGA5="#6A51A3"),guide="none")+
  labs(x="distance from FAP-CAF seed (rings)", y="mean z-score",
       title="e  α5β1 (ITGA5) — secondary FN1-driven receiver",
       subtitle="ITGA5 also tracks the niche but is rank-7 (a) — the collagen→ITGA2 axis leads")+
  pa()+theme(legend.position=c(0.82,0.86), plot.subtitle=element_text(size=5.2))

final <- (p_a | p_b) / p_c / (p_d | p_e) + plot_layout(heights=c(1.0,0.78,0.95))
agg_png(file.path(FIGD,"Supp_ITGA2_axis.png"), width=180, height=178, units="mm", res=360)
print(final); dev.off()
cat("DONE -> ", file.path(FIGD,"Supp_ITGA2_axis.png"), "\n")
