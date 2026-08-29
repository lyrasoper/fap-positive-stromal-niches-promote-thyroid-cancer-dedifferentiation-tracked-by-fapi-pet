# Cancer Research submission - figure code release
# Builds: Figure 5I gene-set determination (which archived gene sets reproduce the
#          original panel). Prints to stdout; writes nothing.

suppressPackageStartupMessages(library(GSVA))
DD <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "mouse_bulk/mTC")
GMT <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "mouse_bulk/reference/mh.all.v2026.1.Mm.symbols.gmt")
mat <- as.matrix(read.csv(file.path(DD,"GSE118022_groupinf.csv"), row.names=1, check.names=FALSE))
storage.mode(mat) <- "numeric"
meta <- read.csv(file.path(DD,"GSE118022_expr.csv"), check.names=FALSE)
mat <- mat[, meta$sample, drop=FALSE]
g <- factor(meta$group, levels=c("mNT","mPTC","mATC","remATC"))
pg <- function(l){x<-strsplit(readLines(l,warn=FALSE),"\t");n<-sapply(x,`[`,1);s<-lapply(x,function(f)f[-(1:2)][nzchar(f[-(1:2)])]);names(s)<-n;s}
h <- pg(GMT)

cand <- list(
 TDS16       = c("Dio1","Dio2","Duox1","Duox2","Foxe1","Glis3","Nkx2-1","Pax8","Slc26a4","Slc5a5","Slc5a8","Tg","Thra","Thrb","Tpo","Tshr"),
 TDS_mdiff   = c("Tg","Tpo","Slc5a5","Nkx2-1","Pax8","Foxe1","Dio1","Dio2","Duox2","Iyd","Tshr","Slc26a4"),
 ECM14       = c("Fap","Col1a1","Col1a2","Col3a1","Postn","Fn1","Sparc","Acta2","Lox","Dcn","Lum","Tagln","Ctgf","Serpine1"),
 ECM10       = c("Fap","Col1a1","Col1a2","Col3a1","Col5a1","Fn1","Dcn","Lum","Postn","Acta2"),
 EMT_hall    = h[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
 EMT_small   = c("Vim","Zeb1","Zeb2","Snai1","Snai2","Twist1","Fn1","Cdh2","Mmp2","Mmp9","Col1a1","Col1a2"),
 TGFB_hall   = h[["HALLMARK_TGF_BETA_SIGNALING"]],
 TGFB_small  = c("Tgfb1","Tgfbr1","Tgfbr2","Smad2","Smad3","Serpine1","Ctgf","Tagln","Acta2","Col1a1","Col1a2")
)
cand <- lapply(cand, function(x) intersect(x, rownames(mat)))
sc <- gsva(ssgseaParam(exprData=mat, geneSets=cand, minSize=4, maxSize=Inf, alpha=0.25, normalize=TRUE), verbose=FALSE)

plotted <- list(
 P1=list(mNT=c(0.218,0.221,0.228,0.230,0.238,0.267), mPTC=c(0.024,0.028,0.046,0.091),
         mATC=c(-0.359,-0.343,-0.309,-0.275,-0.070),
         remATC=c(-0.341,-0.340,-0.328,-0.295,-0.288,-0.266,-0.188,-0.129,0.042)),
 P2=list(mNT=c(0.216,0.226,0.230,0.242,0.242,0.266), mPTC=c(0.273,0.301,0.340,0.372),
         mATC=c(0.420,0.429,0.456,0.470,0.474),
         remATC=c(0.235,0.367,0.393,0.396,0.398,0.409,0.413,0.427,0.446)),
 P3=list(mNT=c(0.225,0.231,0.237,0.240,0.241,0.260), mPTC=c(0.296,0.304,0.332,0.369),
         mATC=c(0.418,0.429,0.438,0.461,0.478),
         remATC=c(0.276,0.369,0.393,0.398,0.400,0.400,0.432,0.445,0.464)),
 P4=list(mNT=c(0.313,0.326,0.331,0.357,0.362,0.365), mPTC=c(0.422,0.461,0.480,0.550),
         mATC=c(0.581,0.587,0.610,0.620,0.640),
         remATC=c(0.398,0.415,0.560,0.561,0.573,0.579,0.581,0.581,0.601))
)
flat <- function(L) unlist(L[c("mNT","mPTC","mATC","remATC")], use.names=FALSE)
srt  <- function(v) unlist(lapply(levels(g), function(k) sort(v[g==k])), use.names=FALSE)

cat(sprintf("%-6s %-12s %8s %8s %8s\n","panel","geneset","rmse","maxdev","cor"))
for (pn in names(plotted)) {
  obs <- flat(plotted[[pn]])
  for (cn in rownames(sc)) {
    pred <- srt(sc[cn,])
    cat(sprintf("%-6s %-12s %8.4f %8.4f %8.4f\n", pn, cn,
        sqrt(mean((pred-obs)^2)), max(abs(pred-obs)), cor(pred,obs)))
  }
  cat("\n")
}
