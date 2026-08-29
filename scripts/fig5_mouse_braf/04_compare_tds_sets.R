# Cancer Research submission - figure code release
# Builds: Figure 5I gene-set determination: TDS variants. Prints to stdout.

suppressPackageStartupMessages(library(GSVA))
DD <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "mouse_bulk/mTC")
GMT <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "mouse_bulk/reference/mh.all.v2026.1.Mm.symbols.gmt")
mat <- as.matrix(read.csv(file.path(DD,"GSE118022_groupinf.csv"), row.names=1, check.names=FALSE))
storage.mode(mat) <- "numeric"
meta <- read.csv(file.path(DD,"GSE118022_expr.csv"), check.names=FALSE)
mat <- mat[, meta$sample, drop=FALSE]
g <- factor(meta$group, levels=c("mNT","mPTC","mATC","remATC"))
pg <- function(l){x<-strsplit(readLines(l,warn=FALSE),"\t");s<-lapply(x,function(f)f[-(1:2)][nzchar(f[-(1:2)])]);names(s)<-sapply(x,`[`,1);s}
h <- pg(GMT)
TDS_A <- c("Dio1","Dio2","Duox1","Duox2","Foxe1","Glis3","Nkx2-1","Pax8","Slc26a4","Slc5a5","Slc5a8","Tg","Thra","Thrb","Tpo","Tshr")   # integrate_thyroid_markers
TDS_S5<- c("Tg","Tpo","Slc5a5","Slc26a4","Tshr","Pax8","Nkx2-1","Dio1","Dio2","Duox1","Duox2","Iyd","Tff3","Foxe1","Glis3","Slc5a8")  # Supp Table S5
ECM14 <- c("Fap","Col1a1","Col1a2","Col3a1","Fn1","Sparc","Postn","Dcn","Lum","Lox","Acta2","Tagln","Ctgf","Serpine1")
obs1 <- c(0.218,0.221,0.228,0.230,0.238,0.267,0.024,0.028,0.046,0.091,
  -0.359,-0.343,-0.309,-0.275,-0.070,-0.341,-0.340,-0.328,-0.295,-0.288,-0.266,-0.188,-0.129,0.042)
srt <- function(v) unlist(lapply(levels(g), function(k) sort(v[g==k])), use.names=FALSE)
for (nm in c("TDS_A","TDS_S5")) {
  tds <- get(nm)
  cat(sprintf("\n== %s : 定义 %d，检出 %d (%s)\n", nm, length(tds),
      length(intersect(tds,rownames(mat))), paste(setdiff(tds,rownames(mat)),collapse="/")))
  sets <- list(Differentiation_TDS=tds,
               EMT_program=h[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
               TGFB_program=h[["HALLMARK_TGF_BETA_SIGNALING"]],
               FAP_CAF_program=ECM14)
  sets <- lapply(sets, function(x) intersect(x, rownames(mat)))
  sc <- gsva(ssgseaParam(exprData=mat, geneSets=sets, minSize=4, maxSize=Inf, alpha=0.25, normalize=TRUE), verbose=FALSE)
  p <- srt(sc["Differentiation_TDS",])
  cat(sprintf("   vs 图面: rmse=%.4f maxdev=%.4f r=%.5f\n", sqrt(mean((p-obs1)^2)), max(abs(p-obs1)), cor(p,obs1)))
  cmps <- list(c("mNT","mPTC"),c("mPTC","mATC"),c("mATC","remATC"),c("mPTC","remATC"))
  for (m in rownames(sc)) {
    raw <- sapply(cmps, function(cp) suppressWarnings(wilcox.test(sc[m,g==cp[1]], sc[m,g==cp[2]], exact=TRUE)$p.value))
    ad <- p.adjust(raw,"holm")
    sy <- ifelse(ad<0.001,"***",ifelse(ad<0.01,"**",ifelse(ad<0.05,"*","ns")))
    cat(sprintf("   %-20s %s\n", m, paste(sprintf("%s=%s(%.3g)", sapply(cmps,paste,collapse="-"), sy, ad), collapse="  ")))
  }
}
