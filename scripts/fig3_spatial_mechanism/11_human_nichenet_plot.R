# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S13A and S13B.

# paper_A figure: HUMAN NicheNet — FAP+ CAF ligands -> epithelial integrin receptors,
# highlighting that CD49b/ITGA2 IS an expressed receptor in human (collagen-alpha2beta1 axis).
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
source("scripts/_shared/paper_A_style.R")
OUT <- "outputs/fig3_ligand_activity/human/nichenet_out"; FAM <- "Helvetica"
la  <- read.csv(file.path(OUT,"human_ligand_activities.csv"))
lri <- read.csv(file.path(OUT,"human_lr_links_integrin.csv"))

is_ecm <- function(g) grepl("^COL|^FBN|^LUM$|^HSPG2$|^MMP|^SDC|^VCAM|^ICAM|^JAM|^CD44$|^LGALS|^THBS|^TNC$|^FN1$|^COMP$|^POSTN$|^ENG$", g)
topL <- la %>% slice_max(aupr_corrected, n=15) %>%
  mutate(test_ligand=factor(test_ligand, levels=rev(test_ligand)),
         cls=ifelse(is_ecm(as.character(test_ligand)),"ECM / adhesion ligand","other"))
th <- pa_theme(base_size=7, base_family=FAM, border=FALSE) +
  theme(plot.tag=element_text(face="bold", size=9, family=FAM), axis.text=element_text(size=6),
        axis.title=element_text(size=6.8), plot.title=element_text(size=7, hjust=0))
pa <- ggplot(topL, aes(aupr_corrected, test_ligand, fill=cls)) +
  geom_col(width=.68, colour="grey30", linewidth=.2) +
  scale_fill_manual(values=c("ECM / adhesion ligand"="#C03840","other"="#9CB3C9")) +
  scale_x_continuous(expand=expansion(mult=c(0,.06))) +
  labs(x="Ligand activity (AUPR corr.)", y=NULL, title="Human FAP+ CAF ligands driving dedifferentiation", tag="a") +
  th + theme(legend.position=c(.66,.16), legend.title=element_blank(),
             legend.text=element_text(size=5), legend.key.size=unit(2.4,"mm"))

recs <- c("ITGB1","ITGA2","ITGAV","ITGA3","ITGB8")
lri2 <- lri %>% filter(from %in% as.character(topL$test_ligand) & to %in% recs)
lri2$from <- factor(lri2$from, levels=rev(levels(topL$test_ligand)))
lri2$to   <- factor(lri2$to, levels=intersect(recs, unique(lri2$to)))
# highlight CD49b/ITGA2 column
pb <- ggplot(lri2, aes(to, from)) +
  geom_tile(aes(fill=to=="ITGA2"), colour="white", width=.9, height=.9) +
  scale_fill_manual(values=c(`TRUE`="#C03840",`FALSE`="#2C5E9D"), guide="none") +
  labs(x="Integrin receptor (human epithelium)", y=NULL, title="Ligand -> integrin receptor (CD49b/ITGA2 in red)", tag="b") +
  pa_theme(base_size=7, base_family=FAM, border=TRUE) +
  theme(plot.tag=element_text(face="bold", size=9, family=FAM), axis.text=element_text(size=6),
        axis.text.x=element_text(angle=40, hjust=1), axis.title=element_text(size=6.8),
        panel.grid=element_blank(), plot.title=element_text(size=7, hjust=0))

g <- pa + pb + plot_layout(widths=c(1.2,1))
ggsave(file.path(OUT,"Human_nichenet_ligand_integrin.pdf"), g, width=175, height=72, units="mm", useDingbats=FALSE)
ggsave(file.path(OUT,"Human_nichenet_ligand_integrin.png"), g, width=175, height=72, units="mm", dpi=400)
cat("saved Human_nichenet_ligand_integrin.{pdf,png}\n")
cat("ITGA2/CD49b receptor links (human):", sum(lri2$to=="ITGA2"), "ligands\n")
