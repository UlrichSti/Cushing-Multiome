# scripts/02_soupx_all.R
# Run SoupX per sample (RNA only) using:
#   tod = RAW RNA matrix (all droplets, do NOT subset to filtered cells)
#   toc = FILTERED RNA matrix (cells)
#   clusters + UMAP from precluster Seurat object
# Outputs:
#   cache/<sample>_soupx_seurat.rds   (local; cache/ should be gitignored)
#   reports/<sample>_soupx_report.csv (small; tracked)

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
})

samples <- read.delim("config/samples.tsv", sep = "\t", stringsAsFactors = FALSE)

pick_rna <- function(x, sample_id) {
  # 10x Multiome H5 returns a list of modalities; we REQUIRE RNA explicitly.
  if (!is.list(x)) return(x)
  if ("Gene Expression" %in% names(x)) return(x[["Gene Expression"]])
  stop(sample_id, ": 'Gene Expression' not found. Available modalities: ",
       paste(names(x), collapse = ", "))
}

run_one <- function(sample_id, path) {
  message("== SoupX: ", sample_id)
  
  precluster_rds <- file.path("cache", paste0(sample_id, "_precluster_seurat.rds"))
  raw_h5 <- file.path(path, "raw_feature_bc_matrix.h5")
  filtered_h5 <- file.path(path, "filtered_feature_bc_matrix.h5")
  
  if (!file.exists(precluster_rds)) stop("Missing precluster object: ", precluster_rds)
  if (!file.exists(raw_h5)) stop("Missing: ", raw_h5)
  if (!file.exists(filtered_h5)) stop("Missing: ", filtered_h5)
  
  # Load preclustered Seurat (for clusters + UMAP)
  so <- readRDS(precluster_rds)
  
  # Read raw + filtered matrices and pick RNA (strict)
  raw_x  <- Read10X_h5(raw_h5)
  filt_x <- Read10X_h5(filtered_h5)
  
  raw_rna  <- pick_rna(raw_x,  sample_id)   # tod: keep ALL droplets
  filt_rna <- pick_rna(filt_x, sample_id)   # toc: filtered cells
  
  # Align Seurat object to filtered matrix cells (toc)
  common <- intersect(colnames(filt_rna), colnames(so))
  if (length(common) < 100) stop("Too few shared barcodes for ", sample_id, ". Check paths/sample match.")
  
  filt_rna <- filt_rna[, common, drop = FALSE]
  so <- subset(so, cells = common)
  
  # Build SoupChannel: IMPORTANT - raw_rna is NOT subsetted
  sc <- SoupChannel(tod = raw_rna, toc = filt_rna)
  
  cat("SoupX toc cells:", ncol(sc$toc), " | tod droplets:", ncol(sc$tod),
      " | Seurat cells:", ncol(so), "\n")
  
  # Provide clusters + 2D embedding
  cl <- Idents(so)
  sc <- setClusters(sc, setNames(as.character(cl), names(cl)))
  
  umap <- Embeddings(so, "umap")[common, , drop = FALSE]
  sc <- setDR(sc, umap)
  
  # Estimate contamination + adjust counts
  sc <- autoEstCont(sc, doPlot = FALSE)
  out <- adjustCounts(sc)
  
  # Put corrected counts back into Seurat (Seurat v5 uses 'layer')
  so_soupx <- so
  so_soupx[["RNA"]] <- SetAssayData(so_soupx[["RNA"]], layer = "counts", new.data = out)
  
  # Save outputs
  dir.create("cache", showWarnings = FALSE)
  dir.create("reports", showWarnings = FALSE)
  
  out_rds <- file.path("cache", paste0(sample_id, "_soupx_seurat.rds"))
  saveRDS(so_soupx, out_rds)
  
  report <- data.frame(
    sample = sample_id,
    n_cells = ncol(so_soupx),
    est_contamination_mean = mean(sc$metaData$rho, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  out_csv <- file.path("reports", paste0(sample_id, "_soupx_report.csv"))
  write.csv(report, out_csv, row.names = FALSE)
  
  message("Saved: ", out_rds)
  message("Report: ", out_csv)
}

for (i in seq_len(nrow(samples))) {
  run_one(samples$sample_id[i], samples$path[i])
}
