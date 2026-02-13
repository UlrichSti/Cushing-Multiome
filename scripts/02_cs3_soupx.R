# scripts/02_cs3_soupx.R
# Run SoupX on CS3 (RNA only) using raw + filtered matrices
# Inputs:
#   - cache/CS3_precluster_seurat.rds  (from scripts/01_cs3_precluster.R)
#   - CS3/raw_feature_bc_matrix.h5
#   - CS3/filtered_feature_bc_matrix.h5
# Outputs:
#   - cache/CS3_soupx_seurat.rds          (NOT for git)
#   - reports/CS3_soupx_report.csv        (small, OK for git)

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
})

sample <- "CS3"
raw_h5 <- file.path(sample, "raw_feature_bc_matrix.h5")
filtered_h5 <- file.path(sample, "filtered_feature_bc_matrix.h5")
precluster_rds <- file.path("cache", paste0(sample, "_precluster_seurat.rds"))

if (!file.exists(precluster_rds)) stop("Missing precluster object: ", precluster_rds)
if (!file.exists(raw_h5)) stop("Missing: ", raw_h5)
if (!file.exists(filtered_h5)) stop("Missing: ", filtered_h5)

# --- Load precluster Seurat (for clusters + UMAP) ---
so <- readRDS(precluster_rds)

# --- Read raw + filtered; pick RNA modality ("Gene Expression") ---
raw_x <- Read10X_h5(raw_h5)
filt_x <- Read10X_h5(filtered_h5)

raw_rna <- if (is.list(raw_x)) {
  if ("Gene Expression" %in% names(raw_x)) raw_x[["Gene Expression"]] else raw_x[[1]]
} else raw_x

filt_rna <- if (is.list(filt_x)) {
  if ("Gene Expression" %in% names(filt_x)) filt_x[["Gene Expression"]] else filt_x[[1]]
} else filt_x

# --- Align cells between Seurat + filtered matrix ---
common <- intersect(colnames(filt_rna), colnames(so))
if (length(common) < 100) stop("Too few shared barcodes between Seurat and filtered matrix. Check sample/paths.")

# Subset everything to the same cells (keeps names consistent)
filt_rna <- filt_rna[, common, drop = FALSE]
raw_rna  <- raw_rna[,  common, drop = FALSE]
so <- subset(so, cells = common)

# Build SoupChannel AFTER subsetting
sc <- SoupChannel(tod = raw_rna, toc = filt_rna)

cat("SoupX toc cells:", ncol(sc$toc), " Seurat cells:", ncol(so), "\n")

# --- Provide clusters + 2D embedding to SoupX (improves estimation) ---
cl <- Idents(so)
sc <- setClusters(sc, setNames(as.character(cl), names(cl)))

umap <- Embeddings(so, "umap")
umap <- umap[common, , drop = FALSE]
sc <- setDR(sc, umap)

# --- Estimate contamination and correct counts ---
sc <- autoEstCont(sc, doPlot = FALSE)
out <- adjustCounts(sc)   # sparse dgCMatrix; do NOT coerce to dense

# --- Put corrected counts back into Seurat ---
so_soupx <- so
# Seurat v5+: use layer instead of slot
so_soupx[["RNA"]] <- SetAssayData(so_soupx[["RNA"]], layer = "counts", new.data = out)

# --- Save outputs ---
dir.create("cache", showWarnings = FALSE)
dir.create("reports", showWarnings = FALSE)

saveRDS(so_soupx, file = file.path("cache", paste0(sample, "_soupx_seurat.rds")))

report <- data.frame(
  sample = sample,
  n_cells = ncol(so_soupx),
  est_contamination_mean = mean(sc$metaData$rho, na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(report, file = file.path("reports", paste0(sample, "_soupx_report.csv")), row.names = FALSE)

cat("Saved:", file.path("cache", paste0(sample, "_soupx_seurat.rds")), "\n")
cat("Report:", file.path("reports", paste0(sample, "_soupx_report.csv")), "\n")
