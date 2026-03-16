# scripts/04_qc_filter_raw_all.R
# Apply baseline QC filters to raw filtered RNA counts for all samples
# Gene filter at object creation: min.cells = 2
# Cell filters:
#   nCount_RNA >= 400
#   nFeature_RNA >= 100
#   nFeature_RNA <= 6000
#   percent.mt <= 15

suppressPackageStartupMessages({
  library(Seurat)
})

samples <- read.delim("config/samples.tsv", sep = "\t", stringsAsFactors = FALSE)

pick_rna <- function(x, sample_id) {
  if (!is.list(x)) return(x)
  if ("Gene Expression" %in% names(x)) return(x[["Gene Expression"]])
  stop(sample_id, ": 'Gene Expression' not found. Available modalities: ",
       paste(names(x), collapse = ", "))
}

run_one <- function(sample_id, path) {
  message("== QC filter: ", sample_id)
  
  filtered_h5 <- file.path(path, "filtered_feature_bc_matrix.h5")
  if (!file.exists(filtered_h5)) stop("Missing: ", filtered_h5)
  
  x <- Read10X_h5(filtered_h5)
  rna <- pick_rna(x, sample_id)
  
  # gene filter: keep genes present in at least 2 cells
  so <- CreateSeuratObject(
    counts = rna,
    project = sample_id,
    min.cells = 2
  )
  
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  
  n_before <- ncol(so)
  
  so_qc <- subset(
    so,
    subset =
      nCount_RNA >= 400 &
      nFeature_RNA >= 100 &
      nFeature_RNA <= 6000 &
      percent.mt <= 15
  )
  
  n_after <- ncol(so_qc)
  
  dir.create("cache", showWarnings = FALSE)
  dir.create("reports", showWarnings = FALSE)
  
  out_rds <- file.path("cache", paste0(sample_id, "_qc_raw_seurat.rds"))
  saveRDS(so_qc, out_rds)
  
  report <- data.frame(
    sample = sample_id,
    n_cells_before = n_before,
    n_cells_after = n_after,
    n_removed = n_before - n_after,
    pct_removed = round(100 * (n_before - n_after) / n_before, 2),
    median_nCount_RNA_after = median(so_qc$nCount_RNA),
    median_nFeature_RNA_after = median(so_qc$nFeature_RNA),
    median_percent_mt_after = median(so_qc$percent.mt),
    stringsAsFactors = FALSE
  )
  
  out_csv <- file.path("reports", paste0(sample_id, "_qc_filter_report.csv"))
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

write.csv(all_reports, file = "reports/qc_filter_summary.csv", row.names = FALSE)

message("Saved: reports/qc_filter_summary.csv")
print(all_reports)