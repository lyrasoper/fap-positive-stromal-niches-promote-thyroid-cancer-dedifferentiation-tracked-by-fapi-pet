# Cancer Research submission - figure code release
# Builds: Figure 2 (Visium spot coordinates bridge).

# Bridge: extract spot x/y tissue coordinates from the integrated Visium
# Seurat object so that scripts/fig2_spatial/03_panels_paperA.py can render
# Fig 2 panels B (TDS heatmap on tissue) and C (FAP heatmap on tissue)
# in matplotlib paper_A style.
#
# Output: outputs/figure2_rebuilt/source_data/Fig2_spots_with_coords.csv
#   columns: spot, slice, x, y, FAP_logCP10K, TDS16, TDS16_z

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>")
RDS_PATH    <- file.path(Sys.getenv("EXTERNAL_DATA", unset = "<EXTERNAL_DATA>"), "spatial_visium/Spatial_integrated9P.rds")
OUT_CSV     <- file.path(PROJECT_ROOT,
                         "outputs/figure2_rebuilt/source_data/Fig2_spots_with_coords.csv")
EXISTING    <- file.path(PROJECT_ROOT,
                         "outputs/figure2_rebuilt/source_data/Fig2_spatial_spot_level.csv")

message("[1/4] Loading Seurat RDS (~1.4 GB, ~30 s) ...")
obj <- readRDS(RDS_PATH)
message(sprintf("       loaded: %d cells × %d genes, %d image slot(s)",
                ncol(obj), nrow(obj), length(obj@images)))

message("[2/4] Extracting tissue coordinates per slice ...")
coord_list <- list()
for (img_name in names(obj@images)) {
  img <- obj@images[[img_name]]
  coords <- GetTissueCoordinates(img)
  # GetTissueCoordinates returns a data.frame with rownames = spot barcodes,
  # cols typically `imagerow`, `imagecol` (or x, y depending on Seurat version)
  if (!"imagerow" %in% colnames(coords) && "x" %in% colnames(coords)) {
    coords$imagerow <- coords$y
    coords$imagecol <- coords$x
  }
  coord_list[[img_name]] <- data.frame(
    spot  = rownames(coords),
    slice = img_name,
    x     = coords$imagecol,
    y     = -coords$imagerow,    # flip y so tissue is "right-side-up" in matplotlib
    stringsAsFactors = FALSE
  )
}
coords_all <- do.call(rbind, coord_list)
message(sprintf("       extracted %d spot coordinates across %d slices",
                nrow(coords_all), length(coord_list)))
message("       slice spot counts:")
print(table(coords_all$slice))

message("[3/4] Merging with existing TDS/FAP CSV ...")
existing <- read.csv(EXISTING, stringsAsFactors = FALSE)
message(sprintf("       existing CSV: %d rows × %d cols",
                nrow(existing), ncol(existing)))
# The `spot` column in the existing CSV has format like "AACAC...-1_2" where
# the suffix encodes slice. The img_name in Seurat may or may not match.
# Try direct join first; if not, fall back to (slice, spot_prefix) match.

# Detect format
message("       existing CSV spot examples: ", paste(head(existing$spot, 3), collapse=", "))
message("       coords spot examples:       ", paste(head(coords_all$spot, 3), collapse=", "))

merged <- merge(coords_all, existing, by = c("spot", "slice"),
                all.x = FALSE, all.y = FALSE)
message(sprintf("       merged: %d rows (intersection)", nrow(merged)))

if (nrow(merged) == 0) {
  message("       direct merge failed; trying barcode-only merge ...")
  # Drop slice column from existing for fallback, match by spot only
  merged <- merge(coords_all, existing[, !(colnames(existing) %in% "slice")],
                  by = "spot", all.x = FALSE, all.y = FALSE)
  message(sprintf("       merged (barcode-only): %d rows", nrow(merged)))
}

message("[4/4] Writing CSV ...")
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
write.csv(merged, OUT_CSV, row.names = FALSE)
message(sprintf("       wrote: %s (%.1f KB)",
                OUT_CSV, file.info(OUT_CSV)$size / 1024))
message("DONE.")
