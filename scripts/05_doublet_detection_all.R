# scripts/05_doublet_detection_all.R
# Doublet detection on QC-filtered raw RNA objects for all samples
# Input:
#   cache/<sample>_qc_raw_seurat.rds
# Output:
#   cache/<sample>_qc_singlets_seurat.rds
#   reports/<sample>_doublet_report.csv
#   reports/doublet_summary.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
})

samples <- read.delim("config/samples.tsv", sep = "\t", stringsAsFactors = FALSE)

run_one <- function(sample_id) {
  message("== Doublet detection: ", sample_id)
  
  in_rds <- file.path("cache", paste0(sample_id, "_qc_raw_seurat.rds"))
  if (!file.exists(in_rds)) stop("Missing: ", in_rds)
  
  so <- readRDS(in_rds)
  
  sce <- as.SingleCellExperiment(so)
  sce <- scDblFinder(sce)
  
  # add results back to Seurat metadata
  so$scDblFinder.score <- colData(sce)$scDblFinder.score
  so$scDblFinder.class <- colData(sce)$scDblFinder.class
  
  n_before <- ncol(so)
  
  so_singlet <- subset(so, subset = scDblFinder.class == "singlet")
  
  n_after <- ncol(so_singlet)
  n_doublets <- n_before - n_after
  
  dir.create("cache", showWarnings = FALSE)
  dir.create("reports", showWarnings = FALSE)
  
  out_rds <- file.path("cache", paste0(sample_id, "_qc_singlets_seurat.rds"))
  saveRDS(so_singlet, out_rds)
  
  report <- data.frame(
    sample = sample_id,
    n_cells_before = n_before,
    n_singlets = n_after,
    n_doublets = n_doublets,
    pct_doublets = round(100 * n_doublets / n_before, 2),
    mean_doublet_score = mean(so$scDblFinder.score, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  out_csv <- file.path("reports", paste0(sample_id, "_doublet_report.csv"))
  write.csv(report, out_csv, row.names = FALSE)
  
  message("Saved: ", out_rds)
  message("Report: ", out_csv)
  
  report
}

all_reports <- do.call(
  rbind,
  lapply(samples$sample_id, run_one)
)

write.csv(all_reports, file = "reports/doublet_summary.csv", row.names = FALSE)

message("Saved: reports/doublet_summary.csv")
print(all_reports)