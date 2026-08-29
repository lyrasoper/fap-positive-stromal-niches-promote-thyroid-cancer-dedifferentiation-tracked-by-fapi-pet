# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S10 (NicheNet input preparation).

# NicheNet input prep: sender = FAP+ inflaCAF, receiver = malignant epithelium.
# geneset of interest = dedifferentiation signature (genes UP in terminal/EMT/ECM vs Differentiated).
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
PROJECT_ROOT  <- Sys.getenv("PROJECT_ROOT",  unset = "<PROJECT_ROOT>")
EXTERNAL_DATA <- Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>")
OUT <- file.path(PROJECT_ROOT, "outputs/fig3_ligand_activity/nichenet_out")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

obj_path <- file.path(EXTERNAL_DATA, "data/fig3/epi_fib_imm_endo_seurat_object.rds")
if(!file.exists(obj_path)) stop(sprintf("Missing %s — produced by an upstream step; see README.", obj_path))
obj <- readRDS(obj_path)
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose=FALSE)

states <- c("Differentiated TCs","Proliferative TCs","Inflammatory TCs",
            "EMT-Like TCs","ECM-Remodeling TCs","Terminally TCs")
epi <- subset(obj, cells = rownames(obj@meta.data)[obj@meta.data$celltype=="Epithelial" &
                                                   obj@meta.data$epi_sub %in% states])

# receiver geneset_oi: dedifferentiation signature (UP in dedifferentiated vs differentiated)
Idents(epi) <- epi$epi_sub
set.seed(1)
deg <- FindMarkers(epi, ident.1=c("EMT-Like TCs","ECM-Remodeling TCs","Terminally TCs"),
                   ident.2="Differentiated TCs", min.pct=0.1, logfc.threshold=0.25, verbose=FALSE)
geneset_oi <- rownames(deg)[deg$p_val_adj < 0.05 & deg$avg_log2FC > 0.25]
cat("dedifferentiation geneset_oi (up):", length(geneset_oi), "genes\n")

expr_genes <- function(o, thr=0.10){
  m <- GetAssayData(o, assay="RNA", layer="counts"); rownames(m)[Matrix::rowMeans(m > 0) >= thr]
}
background      <- expr_genes(epi, 0.05)                 # receiver background
expr_receiver  <- expr_genes(epi, 0.10)                  # receptors live here
caf <- subset(obj, cells = rownames(obj@meta.data)[obj@meta.data$celltype_detail=="FAP+ inflaCAF"])
expr_sender    <- expr_genes(caf, 0.10)                  # ligands live here
cat("sender FAP+ inflaCAF cells:", ncol(caf), "| receiver epi cells:", ncol(epi), "\n")
cat("expr sender:", length(expr_sender), "| expr receiver:", length(expr_receiver),
    "| background:", length(background), "\n")

saveRDS(list(geneset_oi=geneset_oi, background=background,
             expr_sender=expr_sender, expr_receiver=expr_receiver, deg=deg),
        file.path(OUT, "nichenet_inputs.rds"))
write.csv(deg, file.path(OUT, "dedifferentiation_DEG.csv"))
cat("saved nichenet_inputs.rds\n")
