# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S13 (human NicheNet ligand activity).
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

# HUMAN NicheNet: FIBRO_FAP_CAF_core (sender) -> epithelium (receiver) ligands
# explaining the human dedifferentiation signature; integrin-receptor focus.
suppressPackageStartupMessages({library(nichenetr); library(Seurat); library(dplyr); library(tidyr); library(Matrix)})
CACHE <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/results/20260421_fig7_liana_nichenet/nichenet_cache")
SRC   <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/results/20260421_fig7_liana_nichenet/subset_ccc_labelled.rds")
OUT   <- "outputs/fig3_ligand_activity/human/nichenet_out"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

ltm <- readRDS(file.path(CACHE,"ligand_target_matrix.rds"))
lr  <- readRDS(file.path(CACHE,"lr_network.rds")) %>% distinct(from, to)
o   <- readRDS(SRC); DefaultAssay(o) <- "RNA"; o <- NormalizeData(o, verbose=FALSE)

SEND <- "FIBRO_FAP_CAF_core"
RECV <- c("EPI_lineage_preserved_epithelial","EPI_terminal_dedifferentiated")
Idents(o) <- "ccc_label"
deg <- FindMarkers(o, ident.1="EPI_terminal_dedifferentiated", ident.2="EPI_lineage_preserved_epithelial",
                   min.pct=0.1, logfc.threshold=0.25, verbose=FALSE)
geneset_oi <- intersect(rownames(deg)[deg$p_val_adj<0.05 & deg$avg_log2FC>0.25], rownames(ltm))
cat("geneset_oi (human dediff up):", length(geneset_oi), "\n")

cnt <- GetAssayData(o, assay="RNA", layer="counts")
expr <- function(cells, thr=0.10){ m<-cnt[,cells,drop=FALSE]; rownames(m)[Matrix::rowMeans(m>0)>=thr] }
recv_cells <- colnames(o)[o$ccc_label %in% RECV]
send_cells <- colnames(o)[o$ccc_label == SEND]
background    <- intersect(expr(recv_cells,0.05), rownames(ltm))
expr_receiver <- expr(recv_cells,0.10)
expr_sender   <- expr(send_cells,0.10)
cat("sender cells:", length(send_cells), "| receiver cells:", length(recv_cells), "\n")

ligands<-unique(lr$from); receptors<-unique(lr$to)
potential_ligands <- lr %>% filter(from %in% intersect(ligands,expr_sender) & to %in% intersect(receptors,expr_receiver)) %>%
  pull(from) %>% unique() %>% intersect(colnames(ltm))
integ <- c("ITGB1","ITGA2","ITGA5","ITGAV","ITGB3","ITGB5","ITGA1","ITGA3","ITGA6","ITGA11","ITGB4","ITGB6","ITGB8")
cat("** CD49b/ITGA2 an EXPRESSED epithelial receptor in human?:", "ITGA2" %in% expr_receiver, "**\n")
cat("integrin receptors expressed (human epi):", paste(intersect(integ, expr_receiver), collapse=", "), "\n")

la <- predict_ligand_activities(geneset=geneset_oi, background_expressed_genes=background,
        ligand_target_matrix=ltm, potential_ligands=potential_ligands) %>%
      arrange(desc(aupr_corrected)) %>% mutate(rank=row_number())
write.csv(la, file.path(OUT,"human_ligand_activities.csv"), row.names=FALSE)
cat("\n== top 20 human ligands ==\n"); print(as.data.frame(head(la,20)))

best <- la %>% slice_max(aupr_corrected, n=25) %>% pull(test_ligand)
lr_int <- lr %>% filter(from %in% best & to %in% intersect(expr_receiver, integ)) %>% arrange(from,to)
write.csv(lr_int, file.path(OUT,"human_lr_links_integrin.csv"), row.names=FALSE)
cat("\n== top human ligands -> integrin receptors (incl. CD49b/ITGA2) ==\n"); print(as.data.frame(lr_int))
saveRDS(list(deg=deg, geneset_oi=geneset_oi, expr_receiver=expr_receiver, expr_sender=expr_sender), file.path(OUT,"human_nichenet_inputs.rds"))
cat("\nHUMAN_NICHENET_DONE | best:", paste(head(best,10),collapse=", "), "\n")
