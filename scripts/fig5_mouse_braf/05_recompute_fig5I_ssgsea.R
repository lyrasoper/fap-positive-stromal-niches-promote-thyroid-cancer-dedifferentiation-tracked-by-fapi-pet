# Cancer Research submission - figure code release
# Builds: Figure 5I ssGSEA recomputation (GSVA) and pairwise statistics.
#          Writes outputs/fig5I/, which 07_build_panel_I.py reads.

# Figure 5I 重算：真 ssGSEA（GSVA）+ 预设两两检验（2026-08-24）
#
# 原图注标 "ssGSEA score"，但归档中唯一能复现的打分脚本
# (111_bulk_mTC/scripts/plot_rank_signature_per_dataset.py) 用的是秩百分位均值，
# 且其数值与图面不符；EMT / TGFβ 两个模块的逐样本分数从未归档。
# 因此本脚本以公开可复现的方法整体重算，并输出完整逐样本源数据与统计。
#
# 方法：GSVA::gsva(ssgseaParam(...))，kcdf = "Gaussian"（输入为 log2 表达），
#       ssgsea.norm = TRUE（GSVA 默认，全样本极差归一）。
#
# 基因集来源（全部可引用、可复现）：
#   Differentiation_TDS  TDS16（TCGA THCA, Cell 2014）小鼠同源符号
#   EMT_program          MSigDB Hallmark v2026.1.Mm  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
#   TGFB_program         MSigDB Hallmark v2026.1.Mm  HALLMARK_TGF_BETA_SIGNALING
#   FAP_CAF_program      本文 FAP/ECM 模块 14 基因（与 Fig 1 / Supp 一致）
#
# 敏感性：同时用原小规模自定义 EMT / TGFβ 基因列表复算一遍。

suppressPackageStartupMessages({
  library(GSVA)
})

DATA_DIR <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "mouse_bulk/mTC")
GMT      <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "mouse_bulk/reference/mh.all.v2026.1.Mm.symbols.gmt")
OUT      <- file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/fig5I")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

GROUP_ORDER <- c("mNT", "mPTC", "mATC", "remATC")
# 与归档脚本 plot_rank_signature_per_dataset.py 的 COMPARISONS["GSE118022"] 一致
COMPARISONS <- list(c("mNT", "mPTC"), c("mPTC", "mATC"),
                    c("mATC", "remATC"), c("mPTC", "remATC"))

parse_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  sets <- lapply(lines, function(l) {
    f <- strsplit(l, "\t", fixed = TRUE)[[1]]
    f[-(1:2)][nzchar(f[-(1:2)])]
  })
  names(sets) <- vapply(strsplit(lines, "\t", fixed = TRUE), `[`, "", 1)
  sets
}

# ── 数据 ────────────────────────────────────────────────────────────────────
# 注意：这两个文件的名字与内容是反的（GSE118022_expr.csv 才是分组表）
mat <- read.csv(file.path(DATA_DIR, "GSE118022_groupinf.csv"),
                row.names = 1, check.names = FALSE)
mat <- as.matrix(mat); storage.mode(mat) <- "numeric"
meta <- read.csv(file.path(DATA_DIR, "GSE118022_expr.csv"), check.names = FALSE)
stopifnot(all(meta$sample %in% colnames(mat)))
mat <- mat[, meta$sample, drop = FALSE]
meta$group <- factor(meta$group, levels = GROUP_ORDER)
cat(sprintf("矩阵 %d 基因 × %d 样本；分组 %s\n", nrow(mat), ncol(mat),
            paste(sprintf("%s=%d", levels(meta$group), table(meta$group)),
                  collapse = ", ")))

hall <- parse_gmt(GMT)

tds16 <- c("Dio1","Dio2","Duox1","Duox2","Foxe1","Glis3","Nkx2-1","Pax8",
           "Slc26a4","Slc5a5","Slc5a8","Tg","Thra","Thrb","Tpo","Tshr")
fap_caf <- c("Fap","Col1a1","Col1a2","Col3a1","Postn","Fn1","Sparc",
             "Acta2","Lox","Dcn","Lum","Tagln","Ctgf","Serpine1")
emt_small  <- c("Vim","Zeb1","Zeb2","Snai1","Snai2","Twist1","Fn1","Cdh2",
                "Mmp2","Mmp9","Col1a1","Col1a2")
tgfb_small <- c("Tgfb1","Tgfbr1","Tgfbr2","Smad2","Smad3","Serpine1","Ctgf",
                "Tagln","Acta2","Col1a1","Col1a2")

sets_primary <- list(
  Differentiation_TDS = tds16,
  EMT_program         = hall[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
  TGFB_program        = hall[["HALLMARK_TGF_BETA_SIGNALING"]],
  FAP_CAF_program     = fap_caf
)
sets_sens <- list(
  Differentiation_TDS = tds16,
  EMT_program         = emt_small,
  TGFB_program        = tgfb_small,
  FAP_CAF_program     = fap_caf
)

coverage <- do.call(rbind, lapply(names(sets_primary), function(n) {
  g <- sets_primary[[n]]; k <- intersect(g, rownames(mat))
  data.frame(module = n, genes_defined = length(g), genes_detected = length(k),
             genes_missing = paste(setdiff(g, rownames(mat)), collapse = ";"),
             stringsAsFactors = FALSE)
}))
write.csv(coverage, file.path(OUT, "fig5I_gene_coverage.csv"), row.names = FALSE)
print(coverage[, 1:3])

run_ssgsea <- function(sets) {
  sets <- lapply(sets, function(g) intersect(g, rownames(mat)))
  par <- ssgseaParam(exprData = mat, geneSets = sets,
                     minSize = 4, maxSize = Inf,
                     alpha = 0.25, normalize = TRUE)
  gsva(par, verbose = FALSE)
}

scores_primary <- run_ssgsea(sets_primary)
scores_sens    <- run_ssgsea(sets_sens)

MODULES <- c("Differentiation_TDS", "EMT_program", "TGFB_program", "FAP_CAF_program")

tidy <- function(sc, tag) {
  do.call(rbind, lapply(MODULES, function(m)
    data.frame(sample_id = colnames(sc), group = meta$group,
               tissue = meta$tissue, module = m,
               score = as.numeric(sc[m, ]), gene_set = tag,
               stringsAsFactors = FALSE)))
}
long <- rbind(tidy(scores_primary, "primary"), tidy(scores_sens, "sensitivity"))
write.csv(long, file.path(OUT, "fig5I_ssGSEA_scores_recomputed.csv"), row.names = FALSE)

# ── 统计：预设 4 组两两比较 × 4 模块，精确 Wilcoxon 秩和 + Holm（模块内） ──
stat_rows <- list()
for (tag in c("primary", "sensitivity")) {
  d0 <- long[long$gene_set == tag, ]
  for (m in MODULES) {
    d <- d0[d0$module == m, ]
    kw <- kruskal.test(score ~ factor(group, levels = GROUP_ORDER), data = d)
    raw <- vapply(COMPARISONS, function(cp) {
      x <- d$score[d$group == cp[1]]; y <- d$score[d$group == cp[2]]
      suppressWarnings(wilcox.test(x, y, exact = TRUE)$p.value)
    }, numeric(1))
    adj <- p.adjust(raw, method = "holm")
    for (i in seq_along(COMPARISONS)) {
      stat_rows[[length(stat_rows) + 1]] <- data.frame(
        gene_set = tag, module = m,
        kruskal_p = kw$p.value,
        comparison = paste(COMPARISONS[[i]], collapse = " vs "),
        n1 = sum(d$group == COMPARISONS[[i]][1]),
        n2 = sum(d$group == COMPARISONS[[i]][2]),
        median1 = median(d$score[d$group == COMPARISONS[[i]][1]]),
        median2 = median(d$score[d$group == COMPARISONS[[i]][2]]),
        p_raw = raw[i], p_holm = adj[i], stringsAsFactors = FALSE)
    }
  }
}
stats <- do.call(rbind, stat_rows)
sym <- function(p) ifelse(is.na(p), "NA",
        ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", "ns"))))
stats$symbol_holm <- sym(stats$p_holm)
stats$symbol_raw  <- sym(stats$p_raw)
write.csv(stats, file.path(OUT, "fig5I_pairwise_stats.csv"), row.names = FALSE)

cat("\n── primary（Hallmark EMT/TGFβ）──\n")
p <- stats[stats$gene_set == "primary", ]
for (m in MODULES) {
  s <- p[p$module == m, ]
  cat(sprintf("%-20s Kruskal P=%.3g\n", m, s$kruskal_p[1]))
  for (i in seq_len(nrow(s)))
    cat(sprintf("   %-18s n=%d/%d  median %+.3f→%+.3f  raw=%.4g  holm=%.4g  %s(holm) %s(raw)\n",
        s$comparison[i], s$n1[i], s$n2[i], s$median1[i], s$median2[i],
        s$p_raw[i], s$p_holm[i], s$symbol_holm[i], s$symbol_raw[i]))
}
cat("\n── 各组中位数（primary）──\n")
agg <- aggregate(score ~ module + group, data = long[long$gene_set == "primary", ], FUN = median)
print(reshape(agg, idvar = "module", timevar = "group", direction = "wide"))

cat("\n── 敏感性（原小基因集）符号是否一致 ──\n")
s2 <- stats[stats$gene_set == "sensitivity", ]
cmp <- merge(p[, c("module","comparison","symbol_holm")],
             s2[, c("module","comparison","symbol_holm")],
             by = c("module","comparison"), suffixes = c("_primary","_sens"))
print(cmp)
cat(sprintf("\n一致 %d / %d\n", sum(cmp$symbol_holm_primary == cmp$symbol_holm_sens), nrow(cmp)))

sink(file.path(OUT, "sessionInfo.txt")); print(sessionInfo()); sink()
cat("\n完成，输出于 ", OUT, "\n")
