# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S13C and S13D (ITGA2/CD49b vs TDS).

# VERIFY (low-mem, cached matrix): does CD49b/ITGA2 (and the alpha2beta1 collagen-receptor
# signature) anti-correlate with the thyroid differentiation score (TDS) in human epithelium?
# If yes -> positive evidence "receptor up in dedifferentiated cells" (replaces the defensive null panel).
suppressPackageStartupMessages({library(Matrix)})
# epi_counts_export/{counts.mtx,genes.txt,barcodes.txt,metadata.csv} are written by
# 09_export_epithelial_counts.R in this directory.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>")
IN <- file.path(PROJECT_ROOT, "outputs/fig3_ligand_activity/human/epi_counts_export")
for(f in c("counts.mtx","genes.txt","barcodes.txt","metadata.csv")){
  p <- file.path(IN,f)
  if(!file.exists(p)) stop(sprintf("Missing %s — run 09_export_epithelial_counts.R first.", p))
}
X  <- as(Matrix::readMM(file.path(IN,"counts.mtx")), "CsparseMatrix")   # genes x cells
rownames(X) <- readLines(file.path(IN,"genes.txt")); colnames(X) <- readLines(file.path(IN,"barcodes.txt"))
meta <- read.csv(file.path(IN,"metadata.csv"), row.names=1)
cells <- intersect(colnames(X), rownames(meta)); X <- X[,cells]; meta <- meta[cells,]
Xn <- log1p(sweep(X, 2, pmax(Matrix::colSums(X),1), "/")*1e4)
TDS <- meta$TDS_Score1
sc <- function(g){ g<-intersect(g,rownames(Xn)); if(!length(g)) return(rep(NA,ncol(Xn)))
  if(length(g)==1) as.numeric(Xn[g,]) else rowMeans(scale(t(as.matrix(Xn[g,,drop=FALSE])))) }

cat("== per-cell Spearman vs TDS (negative = up in dedifferentiated) ==\n")
feats <- list(`ITGA2 (CD49b)`="ITGA2", `ITGB1 (CD29)`="ITGB1",
              `alpha2beta1 (ITGA2+ITGB1)`=c("ITGA2","ITGB1"),
              `collagen-receptor (ITGA2/ITGB1/ITGA1/ITGA11)`=c("ITGA2","ITGB1","ITGA1","ITGA11"),
              `focal-adhesion/FAK (PTK2/PXN/TLN1/VCL/ZYX)`=c("PTK2","PXN","TLN1","VCL","ZYX"))
for(nm in names(feats)){ v<-sc(feats[[nm]]); ok<-is.finite(v)&is.finite(TDS)
  cat(sprintf("  %-46s rho=%+.3f  p=%.1e  (genes: %s)\n", nm,
      cor(v[ok],TDS[ok],method="spearman"),
      cor.test(v[ok],TDS[ok],method="spearman",exact=FALSE)$p.value,
      paste(intersect(feats[[nm]],rownames(Xn)),collapse="/"))) }

cat("\n== alpha2beta1 signature mean by group (NT/PTC/ATC) ==\n")
a2b1 <- sc(c("ITGA2","ITGB1"))
print(round(tapply(a2b1, meta$grp, mean, na.rm=TRUE),3))
cat("TDS mean by group:\n"); print(round(tapply(TDS, meta$grp, mean, na.rm=TRUE),3))

## save tidy per-cell table for plotting (ITGA2 vs TDS panel)
out <- data.frame(ITGA2=as.numeric(Xn["ITGA2",]), ITGB1=as.numeric(Xn["ITGB1",]),
                  TDS=TDS, grp=meta$grp)
out_csv <- file.path(PROJECT_ROOT, "outputs/fig3_ligand_activity/human/cd49b_vs_tds.csv")
dir.create(dirname(out_csv), recursive=TRUE, showWarnings=FALSE)
write.csv(out, out_csv, row.names=FALSE)
cat("\nsaved cd49b_vs_tds.csv  rows:", nrow(out), "\n")
