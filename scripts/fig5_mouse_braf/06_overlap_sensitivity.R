# Cancer Research submission - figure code release
# Builds: Figure 5I sensitivity analysis: leave-one-out for the 13/14-gene overlap
#          between the Hallmark EMT program and the FAP-CAF program.
#          Prints to stdout.

# EMT_program（Hallmark）与 FAP_CAF_program 共享 13/14 基因。
# 留一验证：把 FAP-CAF 的 14 个基因从 EMT / TGFβ 中剔除后重算，看判读是否改变。
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
TDS <- c("Dio1","Dio2","Duox1","Duox2","Foxe1","Glis3","Nkx2-1","Pax8","Slc26a4","Slc5a5","Slc5a8","Tg","Thra","Thrb","Tpo","Tshr")
ECM <- c("Fap","Col1a1","Col1a2","Col3a1","Fn1","Sparc","Postn","Dcn","Lum","Lox","Acta2","Tagln","Ctgf","Serpine1")
EMT <- h[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]]; TGF <- h[["HALLMARK_TGF_BETA_SIGNALING"]]
CMP <- list(c("mNT","mPTC"),c("mPTC","mATC"),c("mATC","remATC"),c("mPTC","remATC"))
sym <- function(p) ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","ns")))
run <- function(sets){
  sets <- lapply(sets, function(x) intersect(x, rownames(mat)))
  sc <- gsva(ssgseaParam(exprData=mat, geneSets=sets, minSize=4, maxSize=Inf, alpha=0.25, normalize=TRUE), verbose=FALSE)
  out <- list()
  for (m in rownames(sc)) {
    raw <- sapply(CMP, function(cp) suppressWarnings(wilcox.test(sc[m,g==cp[1]], sc[m,g==cp[2]], exact=TRUE)$p.value))
    out[[m]] <- list(q=p.adjust(raw,"holm"), med=tapply(sc[m,], g, median))
  }
  list(sc=sc, st=out)
}
A <- run(list(Differentiation_TDS=TDS, EMT_program=EMT, TGFB_program=TGF, FAP_CAF_program=ECM))
B <- run(list(Differentiation_TDS=TDS, EMT_program=setdiff(EMT,ECM),
              TGFB_program=setdiff(TGF,ECM), FAP_CAF_program=ECM))
cat(sprintf("EMT: %d → %d 基因；TGFβ: %d → %d（剔除 FAP-CAF 的 14 个基因后）\n\n",
    length(intersect(EMT,rownames(mat))), length(intersect(setdiff(EMT,ECM),rownames(mat))),
    length(intersect(TGF,rownames(mat))), length(intersect(setdiff(TGF,ECM),rownames(mat)))))
nm <- sapply(CMP, paste, collapse="-")
same <- 0; tot <- 0
for (m in names(A$st)) {
  a <- sym(A$st[[m]]$q); b <- sym(B$st[[m]]$q)
  cat(sprintf("%-20s primary  %s\n%-20s leave-out %s   %s\n", m, paste(sprintf("%s=%s(%.3g)", nm, a, A$st[[m]]$q), collapse="  "),
      "", paste(sprintf("%s=%s(%.3g)", nm, b, B$st[[m]]$q), collapse="  "),
      if (all(a==b)) "一致" else "★ 有变化"))
  same <- same + sum(a==b); tot <- tot + length(a)
}
cat(sprintf("\n判读一致 %d / %d\n", same, tot))
cat("\nSpearman(primary vs leave-out) 每个模块：\n")
for (m in rownames(A$sc)) cat(sprintf("  %-20s %.4f\n", m, cor(A$sc[m,], B$sc[m,], method="spearman")))
