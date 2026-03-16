# scripts/07_combined_qc_filter_all.R
# Combine RNA QC + singlet status + ATAC QC into final per-sample filtered objects
# Input:
#   cache/<sample>_qc_singlets_seurat.rds
#   reports/<sample>_atac_qc_cells.csv
# Output:
#   cache/<sample>_final_qc_seurat.rds
#   reports/<sample>_combined_qc_report.csv
#   reports/combined_qc_summary.csv

suppressPackageStartupMessages({
  library(Seurat)
})

samples <- read.delim("config/samples.tsv", sep = "\t", stringsAsFactors = FALSE)

run_one <- function(sample_id) {
  message("== Combined QC: ", sample_id)
  
  in_rds <- file.path("cache", paste0(sample_id, "_qc_singlets_seurat.rds"))
  atac_csv <- file.path("reports", paste0(sample_id, "_atac_qc_cells.csv"))
  
  if (!file.exists(in_rds)) stop("Missing: ", in_rds)
  if (!file.exists(atac_csv)) stop("Missing: ", atac_csv)
  
  so <- readRDS(in_rds)
  atac_qc <- read.csv(atac_csv, stringsAsFactors = FALSE)
  
  rownames(atac_qc) <- atac_qc$cell
  
  common <- intersect(colnames(so), rownames(atac_qc))
  if (length(common) < 100) stop("Too few shared cells for ", sample_id)
  
  so <- subset(so, cells = common)
  atac_qc <- atac_qc[common, , drop = FALSE]
  
  # add ATAC QC metrics into Seurat metadata
  so$nCount_ATAC <- atac_qc$nCount_ATAC
  so$nFeature_ATAC <- atac_qc$nFeature_ATAC
  so$nucleosome_signal <- atac_qc$nucleosome_signal
  so$TSS.enrichment <- atac_qc$TSS.enrichment
  
  n_before <- ncol(so)
  
  so_final <- subset(
    so,
    subset =
      nCount_RNA >= 400 &
      nFeature_RNA >= 100 &
      nFeature_RNA <= 6000 &
      percent.mt <= 15 &
      nCount_ATAC >= 100 &
      nFeature_ATAC >= 50 &
      TSS.enrichment >= 2 &
      nucleosome_signal <= 4
  )
  
  n_after <- ncol(so_final)
  
  dir.create("cache", showWarnings = FALSE)
  dir.create("reports", showWarnings = FALSE)
  
  out_rds <- file.path("cache", paste0(sample_id, "_final_qc_seurat.rds"))
  saveRDS(so_final, out_rds)
  
  report <- data.frame(
    sample = sample_id,
    n_cells_before = n_before,
    n_cells_after = n_after,
    n_removed = n_before - n_after,
    pct_removed = round(100 * (n_before - n_after) / n_before, 2),
    median_nCount_RNA_after = median(so_final$nCount_RNA),
    median_nFeature_RNA_after = median(so_final$nFeature_RNA),
    median_percent_mt_after = median(so_final$percent.mt),
    median_nCount_ATAC_after = median(so_final$nCount_ATAC),
    median_nFeature_ATAC_after = median(so_final$nFeature_ATAC),
    median_nucleosome_signal_after = median(so_final$nucleosome_signal),
    median_TSS_enrichment_after = median(so_final$TSS.enrichment),
    stringsAsFactors = FALSE
  )
  
  out_csv <- file.path("reports", paste0(sample_id, "_combined_qc_report.csv"))
  write.csv(report, out_csv, row.names = FALSE)
  
  message("Saved: ", out_rds)
  message("Report: ", out_csv)
  
  report
}

all_reports <- do.call(
  rbind,
  lapply(samples$sample_id, run_one)
)

write.csv(all_reports, file = "reports/combined_qc_summary.csv", row.names = FALSE)

message("Saved: reports/combined_qc_summary.csv")
print(all_reports)