#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S13 (ITGA5/fibronectin comparator axis).
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

## Supplementary — the α5β1 (ITGA5/ITGB1) fibronectin-receptor axis (parallel to the ITGA2 supp)
## Built for an honest side-by-side: ITGA5 is the FN1/POSTN receiver, weaker & lower-ranked than ITGA2.
## a LR pairs (ITGA5 highlighted) · b receiver Δ · c ITGA5 spatial maps · d ITGA5 distance gradient
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(ggplot2);library(patchwork);library(ragg);library(stringr)})
r14 <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/results")
FIGD<- file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/itga2_axis"); dir.create(FIGD,showWarnings=FALSE,recursive=TRUE)
A5_c<-"#6A51A3"; other_c<-"#C9CDD2"; lin_c<-"#2E8B57"
pa <- function(base=7) theme_classic(base_size=base,base_family="Helvetica")+
  theme(axis.text=element_text(size=base-1,colour="black"),axis.title=element_text(size=base),
        plot.title=element_text(size=base+1,face="bold",hjust=0),plot.subtitle=element_text(size=base-1.3,colour="grey35"),
        legend.text=element_text(size=base-1.5),legend.title=element_text(size=base-1),
        axis.line=element_line(linewidth=0.4),axis.ticks=element_line(linewidth=0.4),legend.key.size=unit(3,"mm"))

lr <- read_tsv(file.path(r14,"20260410_fap_fibro_epi_mechanism_prioritization/ranked_candidate_interactions.tsv"),show_col_types=FALSE)
ctx<- read_tsv(file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/fig7_redesign/augmented_context_v3.tsv"),show_col_types=FALSE)

## ---- a: full LR ranking, ITGA5 pairs highlighted ----
la <- lr %>% slice_min(priority_rank, n=12) %>%
  mutate(pair=paste0(ligand," → ",receptor), is5=ifelse(grepl("ITGA5",receptor),"ITGA5 pair","other"),
         pair=factor(pair, levels=rev(pair)))
p_a <- ggplot(la, aes(overall_priority_score, pair, fill=is5)) +
  geom_col(width=0.72) +
  scale_fill_manual(values=c(`ITGA5 pair`=A5_c, other=other_c), name=NULL) +
  labs(x="LIANA priority score", y=NULL, title="a  ITGA5 (α5β1) pairs in the ranked LR axis",
       subtitle="ITGA5 = fibronectin/POSTN receiver (POSTN→ITGA5 rank 7, FN1→ITGA5 lower) — below the collagen→ITGA2 pairs") +
  pa() + theme(axis.text.y=element_text(size=5.6), legend.position=c(0.82,0.32), plot.title=element_text(size=7.3),
               plot.subtitle=element_text(size=5.1))

## ---- b: receiver Δ, ITGA5 highlighted ----
rb <- lr %>% group_by(receptor) %>% summarise(cell=first(cell_delta_terminal_vs_lineage),
                spat=first(spatial_delta_terminal_vs_lineage), .groups="drop") %>%
  arrange(cell) %>% mutate(receptor=factor(receptor,levels=receptor)) %>%
  pivot_longer(c(cell,spat), names_to="layer", values_to="delta") %>%
  mutate(layer=recode(layer, cell="single-cell", spat="spatial"),
         hi=ifelse(grepl("ITGA5",receptor),"ITGA5","other"))
p_b <- ggplot(rb, aes(delta, receptor, fill=layer, alpha=hi)) +
  geom_col(position=position_dodge(0.7), width=0.62) + geom_vline(xintercept=0,linewidth=0.3,colour="grey50")+
  scale_fill_manual(values=c(`single-cell`="#4E79A7", spatial="#F2A900"), name=NULL) +
  scale_alpha_manual(values=c(ITGA5=1, other=0.35), guide="none") +
  labs(x="receiver Δ (terminal − lineage)", y=NULL, title="b  Receptor enrichment (ITGA5 highlighted)",
       subtitle="ITGA5 is terminal-enriched but modest vs ITGA2/ITGB1") +
  pa() + theme(axis.text.y=element_text(size=6), legend.position=c(0.8,0.2), plot.subtitle=element_text(size=5.3))

## ---- c: ITGA5 receiver activity spatial maps ----
reps <- c("P12","P17","P26","P44","P83")
cc <- ctx %>% filter(sample %in% reps) %>% group_by(sample) %>%
  mutate(a5=pmin(itga5_receiver_program_activity, quantile(itga5_receiver_program_activity,0.99,na.rm=TRUE))) %>% ungroup() %>%
  mutate(sample=factor(sample, levels=reps))
p_c <- ggplot(cc, aes(x_plot, -y_plot, colour=a5)) +
  geom_point(size=0.22) + facet_wrap(~sample, nrow=1) +
  scale_colour_gradientn(colours=c("#F0F0F0","#DCD0EC","#9E7FCB","#6A51A3","#3F007D"), name="ITGA5\nactivity") +
  coord_equal() + labs(title="c  ITGA5 receiver-program activity (representative Visium samples)") +
  theme_void(base_family="Helvetica") +
  theme(plot.title=element_text(size=7.5,face="bold",hjust=0), strip.text=element_text(size=6.5,face="bold"),
        legend.title=element_text(size=5.6), legend.text=element_text(size=5), legend.key.size=unit(3,"mm"), legend.position="right")

## ---- d: ITGA5 distance gradient ----
brks<-c(-Inf,1,2,3,5,Inf); labs5<-c("0–1","1–2","2–3","3–5",">5")
dd <- ctx %>% filter(is.finite(distance_spot_units)) %>% group_by(sample) %>%
  mutate(ITGA5=scale(itga5_receiver_program_activity)[,1],
         `terminal load`=scale(epi_deconv_load_terminal_dedifferentiated)[,1],
         `lineage load`=scale(epi_deconv_load_lineage_preserved_epithelial)[,1]) %>% ungroup() %>%
  mutate(ring=cut(distance_spot_units, breaks=brks, labels=labs5)) %>%
  pivot_longer(c(ITGA5,`terminal load`,`lineage load`), names_to="var", values_to="z") %>%
  group_by(sample,ring,var) %>% summarise(z=mean(z,na.rm=TRUE), .groups="drop") %>%
  group_by(ring,var) %>% summarise(m=mean(z,na.rm=TRUE), se=sd(z,na.rm=TRUE)/sqrt(n()), .groups="drop")
p_d <- ggplot(dd, aes(ring,m,colour=var,group=var)) +
  geom_hline(yintercept=0,linewidth=0.3,colour="grey75")+
  geom_ribbon(aes(ymin=m-1.96*se,ymax=m+1.96*se,fill=var),alpha=0.15,colour=NA)+
  geom_line(linewidth=0.7)+geom_point(size=1.3)+
  scale_colour_manual(values=c(ITGA5=A5_c,`terminal load`="#7B241C",`lineage load`=lin_c),name=NULL)+
  scale_fill_manual(values=c(ITGA5=A5_c,`terminal load`="#7B241C",`lineage load`=lin_c),guide="none")+
  labs(x="distance from nearest FAP-CAF seed (spot rings)", y="mean z-score",
       title="d  ITGA5 activity vs distance from the FAP-CAF niche",
       subtitle="ITGA5 also tracks the niche spatially — its weakness is rank/enrichment (a, b), not the gradient")+
  pa()+theme(legend.position=c(0.8,0.85), plot.subtitle=element_text(size=5.3))

final <- (p_a | p_b) / p_c / p_d + plot_layout(heights=c(1.0,0.78,0.95))
agg_png(file.path(FIGD,"Supp_ITGA5_axis.png"), width=180, height=176, units="mm", res=360)
print(final); dev.off()
cat("DONE\n")
