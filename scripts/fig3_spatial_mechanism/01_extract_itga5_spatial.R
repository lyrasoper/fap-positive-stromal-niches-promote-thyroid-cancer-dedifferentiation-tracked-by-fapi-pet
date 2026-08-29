# Cancer Research submission - figure code release
# Builds: Figure 3 / Supplementary Fig. S13 (ITGA5 comparator extraction).
# Note: some output paths and inline comments below still use the figure
#       numbering of an earlier manuscript version (fig7_* = current Fig. 3,
#       fig8_* = current Fig. 4). See README.md.

# 01_extract_itga5_spatial.R
# Purpose: For Fig 7 ITGA5 emphasis (Option C), extract:
#   1. ITGA5 per-spot expression from Spatial_integrated9P → add column to
#      outputs/fig7_redesign/mechanism_gene_expression_per_spot.tsv
#   2. ITGA5 single-cell receiver delta (terminal_dedif vs lineage_preserved)
#      from subset_ccc_labelled → add row to
#      .../20260421_fig7_liana_nichenet/tables/A_receiver_receptor_deltas.tsv

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

PROJ <- Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>")
ATLAS_BASE <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "scRNA_atlas/results")
SPATIAL_RDS <- file.path(ATLAS_BASE, "20260409_spatial_epithelial_state_transfer",
                          "Spatial_integrated9P.with_epithelial_state_transfer.rds")
SCSUB_RDS <- file.path(ATLAS_BASE, "20260421_fig7_liana_nichenet", "subset_ccc_labelled.rds")
PER_SPOT_TSV <- file.path(PROJ, "outputs/fig7_redesign/mechanism_gene_expression_per_spot.tsv")
AUG_CTX_TSV <- file.path(PROJ, "outputs/fig7_redesign/augmented_context_v3.tsv")
RECEIVER_TSV <- file.path(ATLAS_BASE, "20260421_fig7_liana_nichenet/tables/A_receiver_receptor_deltas.tsv")

# ----- 1. SPATIAL EXTRACTION -----
cat("[1/2] Loading spatial integrated atlas (2.4 GB) ...\n")
if (!file.exists(SPATIAL_RDS)) stop(sprintf("Missing %s — produced by an upstream step; see README.", SPATIAL_RDS))
spat <- readRDS(SPATIAL_RDS)
cat("  Spatial assays:", paste(names(spat@assays), collapse=", "), "\n")
cat("  Default assay:", DefaultAssay(spat), "\n")

# Try Spatial first, fallback to SCT or RNA
if ("Spatial" %in% names(spat@assays)) DefaultAssay(spat) <- "Spatial"

# 检查 ITGA5 是否存在
genes <- rownames(spat)
cat("  ITGA5 in genes?", "ITGA5" %in% genes, "\n")
if (!("ITGA5" %in% genes)) {
  cat("  → ITGA5 not found, trying SCT assay ...\n")
  if ("SCT" %in% names(spat@assays)) {
    DefaultAssay(spat) <- "SCT"
    genes <- rownames(spat)
    cat("  ITGA5 in SCT?", "ITGA5" %in% genes, "\n")
  }
}

# Get ITGA5 normalized expression (log-normalized, same scale as existing genes)
itga5_expr <- FetchData(spat, vars = "ITGA5", layer = "data")
cat("  ITGA5 expr extracted:", nrow(itga5_expr), "spots\n")
cat("  Range:", round(range(itga5_expr[,1]), 3), "  Mean:", round(mean(itga5_expr[,1]), 3), "\n")

# Match by barcode to existing per-spot file
if (!file.exists(PER_SPOT_TSV)) stop(sprintf("Missing %s — produced by an upstream step; see README.", PER_SPOT_TSV))
per_spot <- read_tsv(PER_SPOT_TSV, show_col_types = FALSE)
cat("\n  Existing per_spot rows:", nrow(per_spot), "\n")
cat("  Existing samples:", paste(unique(per_spot$sample), collapse=","), "\n")

# 提取 spot barcodes 对应
spat_barcodes <- rownames(itga5_expr)
cat("  Spatial barcode example:", head(spat_barcodes, 3), "\n")
cat("  per_spot example:", head(per_spot$spot, 3), "\n")

itga5_df <- data.frame(spot = spat_barcodes, ITGA5 = itga5_expr[,1], stringsAsFactors = FALSE)

# Join to per_spot
merged <- per_spot %>% left_join(itga5_df, by = "spot")
matched <- sum(!is.na(merged$ITGA5))
cat("\n  Matched ITGA5 to", matched, "/", nrow(merged), "spots\n")

if (matched < 0.8 * nrow(merged)) {
  # 可能 barcode 格式不匹配，尝试加 sample suffix
  cat("  Low match. Trying alternative barcode format ...\n")
  # Check if spatial uses underscore sample suffix differently
  sample_suffix <- sub(".*_(\\d+)$", "\\1", spat_barcodes[1])
  cat("  Spatial sample suffix pattern:", sample_suffix, "\n")
}

# Reorder columns to put ITGA5 after ITGA2
existing_cols <- names(per_spot)
itga2_idx <- which(existing_cols == "ITGA2")
new_order <- c(existing_cols[1:itga2_idx], "ITGA5", existing_cols[(itga2_idx+1):length(existing_cols)])
merged <- merged[, new_order]

dir.create(dirname(PER_SPOT_TSV), recursive = TRUE, showWarnings = FALSE)
write_tsv(merged, PER_SPOT_TSV)
cat("\n✓ Updated:", PER_SPOT_TSV, "with ITGA5 column\n")

# 计算 per-sample ITGA5 receiver program activity 
# (与现有 itga2_receiver_program_activity 同方法：z-score normalized within sample)
cat("\n  Computing per-sample ITGA5 z-score for augmented context ...\n")
itga5_act <- merged %>%
  group_by(sample) %>%
  mutate(itga5_receiver_program_activity = scale(ITGA5)[,1]) %>%
  ungroup() %>%
  select(spot, itga5_receiver_program_activity)

if (!file.exists(AUG_CTX_TSV)) stop(sprintf("Missing %s — produced by an upstream step; see README.", AUG_CTX_TSV))
aug <- read_tsv(AUG_CTX_TSV, show_col_types = FALSE)
aug2 <- aug %>% left_join(itga5_act, by = c("spot_id" = "spot"))
cat("  Aug context: matched", sum(!is.na(aug2$itga5_receiver_program_activity)),
    "/", nrow(aug2), "\n")
dir.create(dirname(AUG_CTX_TSV), recursive = TRUE, showWarnings = FALSE)
write_tsv(aug2, AUG_CTX_TSV)
cat("✓ Updated:", AUG_CTX_TSV, "with itga5_receiver_program_activity\n")

# 释放内存
rm(spat); gc(verbose = FALSE)

# ----- 2. SINGLE-CELL RECEIVER DELTA -----
cat("\n[2/2] Loading scRNA subset for receiver delta (0.7 GB) ...\n")
if (!file.exists(SCSUB_RDS)) stop(sprintf("Missing %s — produced by an upstream step; see README.", SCSUB_RDS))
sc <- readRDS(SCSUB_RDS)
cat("  Cells:", ncol(sc), "  Genes:", nrow(sc), "\n")
cat("  Metadata cols:", paste(head(names(sc@meta.data), 20), collapse=", "), "\n")

# 找 epithelial state column
state_col <- NULL
for (c in c("epi_state","epithelial_state","state","cell_state","celltype_state")) {
  if (c %in% names(sc@meta.data)) { state_col <- c; break }
}
cat("  State column:", state_col, "\n")
if (is.null(state_col)) {
  cat("  Searching for column with 'terminal' values ...\n")
  for (c in names(sc@meta.data)) {
    vals <- unique(as.character(sc@meta.data[[c]]))
    if (any(grepl("terminal|lineage_preserved|dedif", vals))) {
      state_col <- c
      cat("  Found:", c, "  values:", paste(head(vals), collapse=", "), "\n")
      break
    }
  }
}

if (!is.null(state_col)) {
  states <- as.character(sc@meta.data[[state_col]])
  term_cells <- colnames(sc)[grepl("terminal_dedif", states)]
  lin_cells <- colnames(sc)[grepl("lineage_preserved", states)]
  cat("\n  Terminal dedif cells:", length(term_cells), "\n")
  cat("  Lineage preserved cells:", length(lin_cells), "\n")
  
  # Get ITGA5 expression
  if ("RNA" %in% names(sc@assays)) DefaultAssay(sc) <- "RNA"
  if ("ITGA5" %in% rownames(sc)) {
    itga5_sc <- FetchData(sc, vars = "ITGA5", layer = "data")[,1]
    e_term <- itga5_sc[term_cells]
    e_lin <- itga5_sc[lin_cells]
    
    # Compute delta + Wilcoxon
    delta_log <- mean(e_term, na.rm=TRUE) - mean(e_lin, na.rm=TRUE)
    pct_term <- mean(e_term > 0, na.rm=TRUE)
    pct_lin <- mean(e_lin > 0, na.rm=TRUE)
    wilcox <- wilcox.test(e_term, e_lin, alternative = "two.sided")
    
    cat("\n  ITGA5 receiver delta:\n")
    cat("    Δ log norm expression:", round(delta_log, 4), "\n")
    cat("    % expressing (terminal):", round(100*pct_term, 1), "%\n")
    cat("    % expressing (lineage):", round(100*pct_lin, 1), "%\n")
    cat("    Wilcoxon p:", format.pval(wilcox$p.value, digits=3), "\n")
    
    # 加到 receiver deltas 表
    if (file.exists(RECEIVER_TSV)) {
      rd <- read_tsv(RECEIVER_TSV, show_col_types = FALSE)
      cat("\n  Existing receiver deltas:", nrow(rd), "rows\n")
      cat("  Columns:", paste(names(rd), collapse=","), "\n")
      
      # 构造新行（匹配现有列结构）
      new_row <- list(
        receptor = "ITGA5",
        delta_log_norm = delta_log,
        wilcox_p = wilcox$p.value,
        pct_terminal = pct_term,
        pct_lineage = pct_lin
      )
      # 添加缺失列为 NA
      for (c in names(rd)) if (!c %in% names(new_row)) new_row[[c]] <- NA
      new_row <- new_row[names(rd)]
      
      # 检查是否已存在 ITGA5
      receptor_col <- intersect(c("receptor","gene","symbol"), names(rd))[1]
      if (!is.null(receptor_col) && "ITGA5" %in% rd[[receptor_col]]) {
        cat("  ITGA5 already in table, updating ...\n")
        rd <- rd %>% filter(!!sym(receptor_col) != "ITGA5")
      }
      rd2 <- bind_rows(rd, as_tibble(new_row))
      
      # BH-adjust p
      if ("wilcox_p" %in% names(rd2)) {
        rd2$wilcox_p_bh <- p.adjust(rd2$wilcox_p, method = "BH")
      }
      
      write_tsv(rd2, RECEIVER_TSV)
      cat("✓ Updated:", RECEIVER_TSV, "with ITGA5 row\n")
    }
  } else {
    cat("  ERROR: ITGA5 not in scRNA RNA assay\n")
  }
} else {
  cat("  ERROR: no state column found\n")
}

cat("\n=== ALL EXTRACTIONS DONE ===\n")
