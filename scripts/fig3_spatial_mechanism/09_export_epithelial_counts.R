# Export the human epithelial count matrix (NT/PTC/ATC atlas) used by
# 12_human_cd49b_vs_tds.R.
#
# Genes kept = 3,000 variable features union the integrin / thyroid-lineage
# targets, so the exported matrix stays small enough to load without the full
# Seurat object.
#
# Writes outputs/fig3_ligand_activity/human/epi_counts_export/
#   counts.mtx    genes x cells, raw counts (MatrixMarket)
#   genes.txt     row names
#   barcodes.txt  column names
#   metadata.csv  per-cell group label and thyroid differentiation score

suppressPackageStartupMessages({library(Seurat); library(Matrix)})

PROJECT_ROOT  <- Sys.getenv("PROJECT_ROOT",  unset = "<PROJECT_ROOT>")
EXTERNAL_DATA <- Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>")

SRC <- file.path(EXTERNAL_DATA, "scRNA_atlas/results/epithelial_foundation/epithelial_object.rds")
if (!file.exists(SRC)) stop(sprintf("Missing %s — set EXTERNAL_DATA; see README.", SRC))

OUT <- file.path(PROJECT_ROOT, "outputs/fig3_ligand_activity/human/epi_counts_export")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(SRC)
if ("RNA" %in% Assays(obj)) DefaultAssay(obj) <- "RNA"
cat("loaded:", paste(dim(obj), collapse = " x "), "| assay:", DefaultAssay(obj), "\n")

grp <- if ("sample_group" %in% colnames(obj@meta.data)) "sample_group" else "Tissue_label_detail8"
obj@meta.data$grp <- as.character(obj@meta.data[[grp]])

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)

targets <- c("ITGB1", "ITGA2", "ITGA5", "ITGAV", "PTK2", "PAX8", "FOXE1", "NKX2-1",
             "TG", "TPO", "SLC5A5", "TSHR", "DIO1", "GLIS3")
keep <- intersect(union(VariableFeatures(obj), targets), rownames(obj))
cts  <- GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts")[keep, ]

Matrix::writeMM(cts, file.path(OUT, "counts.mtx"))
writeLines(rownames(cts), file.path(OUT, "genes.txt"))
writeLines(colnames(cts), file.path(OUT, "barcodes.txt"))

mcols <- intersect(c("grp", "TDS_Score1", "tds", "emt_core_score", "epi_origin_label"),
                   colnames(obj@meta.data))
write.csv(obj@meta.data[, mcols, drop = FALSE], file.path(OUT, "metadata.csv"))

cat("wrote", file.path(OUT, "counts.mtx"), "(", nrow(cts), "x", ncol(cts), ")\n")
