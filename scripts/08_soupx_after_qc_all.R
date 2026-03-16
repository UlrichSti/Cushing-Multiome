# scripts/08_soupx_after_qc_all.R
# Run SoupX on final QC-filtered cells for all samples
# Input:
#   cache/<sample>_final_qc_seurat.rds
#   <sample>/raw_feature_bc_matrix.h5
#   <sample>/filtered_feature_bc_matrix.h5
# Output:
#   cache/<sample>_final_qc_soupx_seurat.rds
#   reports/<sample>_final_soupx_report.csv
#   reports/final_soupx_summary.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
})

samples <- read.delim("config/samples.tsv", sep = "\t", stringsAsFactors = FALSE)

pick_rna <- function(x, sample_id) {
  if (!is.list(x)) return(x)
  if ("Gene Expression" %in% names(x)) return(x[["Gene Expression"]])
  stop(sample_id, ": 'Gene Expression' not found. Available modalities: ",
       paste(names(x), collapse = ", "))
}

run_one <- function(sample_id, path) {
  message("== Final SoupX: ", sample_id)
  
  in_rds <- file.path("cache", paste0(sample_id, "_final_qc_seurat.rds"))
  raw_h5 <- file.path(path, "raw_feature_bc_matrix.h5")
  filtered_h5 <- file.path(path, "filtered_feature_bc_matrix.h5")
  
  if (!file.exists(in_rds)) stop("Missing: ", in_rds)
  if (!file.exists(raw_h5)) stop("Missing: ", raw_h5)
  if (!file.exists(filtered_h5)) stop("Missing: ", filtered_h5)
  
  so <- readRDS(in_rds)
  
  # quick reclustering on final QC cells for SoupX
  so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, verbose = FALSE)
  so <- ScaleData(so, verbose = FALSE)
  so <- RunPCA(so, npcs = 30, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
  so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
  so <- RunUMAP(so, dims = 1:30, verbose = FALSE)
  
  raw_x <- Read10X_h5(raw_h5)
  filt_x <- Read10X_h5(filtered_h5)
  
  raw_rna <- pick_rna(raw_x, sample_id)   # tod = all droplets
  filt_rna <- pick_rna(filt_x, sample_id) # toc = filtered cells
  
  final_cells <- colnames(so)
  common <- intersect(final_cells, colnames(filt_rna))
  if (length(common) < 50) stop("Too few shared cells for ", sample_id)
  
  so <- subset(so, cells = common)
  filt_rna <- filt_rna[, common, drop = FALSE]
  
  sc <- SoupChannel(tod = raw_rna, toc = filt_rna)
  
  cl <- Idents(so)
  sc <- setClusters(sc, setNames(as.character(cl), names(cl)))
  
  umap <- Embeddings(so, "umap")[common, , drop = FALSE]
  sc <- setDR(sc, umap)
  
  sc <- autoEstCont(sc, doPlot = FALSE)
  out <- adjustCounts(sc)
  
  so_soupx <- so
  so_soupx[["RNA"]] <- SetAssayData(so_soupx[["RNA"]], layer = "counts", new.data = out)
  
  dir.create("cache", showWarnings = FALSE)
  dir.create("reports", showWarnings = FALSE)
  
  out_rds <- file.path("cache", paste0(sample_id, "_final_qc_soupx_seurat.rds"))
  saveRDS(so_soupx, out_rds)
  
  report <- data.frame(
    sample = sample_id,
    n_cells = ncol(so_soupx),
    est_contamination_mean = mean(sc$metaData$rho, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  out_csv <- file.path("reports", paste0(sample_id, "_final_soupx_report.csv"))
  write.csv(report, out_csv, row.names = FALSE)
  
  message("Saved: ", out_rds)
  message("Report: ", out_csv)
  
  report
}

all_reports <- do.call(
  rbind,
  lapply(seq_len(nrow(samples)), function(i) {
    run_one(samples$sample_id[i], samples$path[i])
  })
)

write.csv(all_reports, file = "reports/final_soupx_summary.csv", row.names = FALSE)

message("Saved: reports/final_soupx_summary.csv")
print(all_reports)