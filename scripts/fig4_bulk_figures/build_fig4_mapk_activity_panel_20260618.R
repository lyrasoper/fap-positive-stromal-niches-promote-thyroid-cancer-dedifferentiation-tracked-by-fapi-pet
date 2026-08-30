# Cancer Research submission - figure code release
# Builds: Figure 4D - the PROGENy MAPK-output activity row.
## MAPK-output activity panel (PROGENy footprint) — honest downstream support for ECM-integrin/FAK -> MAPK
suppressMessages({library(ggplot2); library(dplyr); library(readr); library(ragg)})
OUT  <- file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/_fig4_bioinformatics_panels_20260613")
DATA <- file.path(Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>"), "outputs/wt_vs_fapko_bulk_omics_20260407/decoupler_crossomics")
base <- 7
col_wt <- "#3B6B9C"; col_fap <- "#C9493A"
grp_levels <- c("WT host","FAP-deficient host"); grp_cols <- setNames(c(col_wt,col_fap), grp_levels)

theme_nc <- function() theme_minimal(base_size=base, base_family="Helvetica") +
  theme(panel.grid=element_blank(),
        axis.line=element_line(linewidth=0.3, colour="grey20"),
        axis.ticks=element_line(linewidth=0.3, colour="grey20"),
        axis.text=element_text(colour="grey20"),
        plot.title=element_text(face="bold", size=base+1),
        plot.subtitle=element_text(size=base-2, colour="grey35", lineheight=1.05),
        plot.tag=element_text(face="bold", size=base+3, hjust=0),
        strip.text=element_text(face="bold", size=base-0.4),
        panel.spacing.x=unit(3,"mm"),
        legend.position="none",
        plot.margin=margin(3,5,3,3))

relabel <- function(d) factor(ifelse(d=="Proteomics","Protein (DIA)","RNA"), levels=c("RNA","Protein (DIA)"))

ps <- read_tsv(file.path(DATA,"progeny_scores_per_sample.tsv"), show_col_types=FALSE) %>%
  filter(pathway=="MAPK") %>%
  mutate(group=factor(group, levels=grp_levels), dataset=relabel(dataset))

st <- read_tsv(file.path(DATA,"progeny_delta_with_pvalues.tsv"), show_col_types=FALSE) %>%
  filter(pathway=="MAPK") %>%
  mutate(dataset=relabel(dataset),
         lab=sprintf("Δ = %+.2f · P = %.2f", delta_fap_minus_wt, wilcox_p))

ann <- ps %>% group_by(dataset) %>%
  summarise(y=max(score)+0.16*(max(score)-min(score)), .groups="drop") %>%
  left_join(st %>% select(dataset,lab), by="dataset")

p <- ggplot(ps, aes(group, score)) +
  geom_jitter(aes(fill=group), shape=21, size=2.1, stroke=0.3, width=0.11, height=0, colour="black", alpha=0.95) +
  stat_summary(fun=mean, geom="crossbar", width=0.46, linewidth=0.45, aes(colour=group)) +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.18, linewidth=0.4, colour="grey35") +
  geom_text(data=ann, aes(x=1.5, y=y, label=lab), size=(base-2)/.pt, colour="grey25", inherit.aes=FALSE) +
  facet_wrap(~dataset, scales="free_y") +
  scale_fill_manual(values=grp_cols) + scale_colour_manual(values=grp_cols) +
  scale_x_discrete(labels=c("WT host"="WT","FAP-deficient host"="FAP-KO")) +
  scale_y_continuous(expand=expansion(mult=c(0.10,0.20))) +
  labs(x=NULL, y="MAPK pathway activity (PROGENy)",
       title="MAPK output declines with stromal FAP loss",
       subtitle=paste0("PROGENy footprint · FAP-KO < WT at RNA and protein (n=4/grp; two-sided Wilcoxon, n.s.)\n",
                       "Downstream readout of the ECM-integrin/FAK step — direct pFAK-Y397 / pERK requires WB/IHC"),
       tag="") +
  theme_nc()

ggsave(file.path(OUT,"panel_mapk_activity.png"), p, width=92, height=66, units="mm", dpi=500, device=ragg::agg_png)
ggsave(file.path(OUT,"panel_mapk_activity.pdf"), p, width=92, height=66, units="mm")
cat("wrote panel_mapk_activity (92x66mm)\n")
