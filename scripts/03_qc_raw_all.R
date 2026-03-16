# scripts/03_qc_raw_all.R
# QC on raw filtered RNA counts for all samples
# Input: <sample>/filtered_feature_bc_matrix.h5
# Output:
#   reports/<sample>_qc_cells.csv
#   reports/qc_summary.csv

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
  message("== QC: ", sample_id)
  
  filtered_h5 <- file.path(path, "filtered_feature_bc_matrix.h5")
  if (!file.exists(filtered_h5)) stop("Missing: ", filtered_h5)
  
  x <- Read10X_h5(filtered_h5)
  rna <- pick_rna(x, sample_id)
  
  so <- CreateSeuratObject(counts = rna, project = sample_id)
  
  # mitochondrial percentage
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  
  qc_df <- data.frame(
    cell = colnames(so),
    sample = sample_id,
    nCount_RNA = so$nCount_RNA,
    nFeature_RNA = so$nFeature_RNA,
    percent.mt = so$percent.mt,
    stringsAsFactors = FALSE
  )
  
  dir.create("reports", showWarnings = FALSE)
  
  out_cells <- file.path("reports", paste0(sample_id, "_qc_cells.csv"))
  write.csv(qc_df, out_cells, row.names = FALSE)
  
  summary_df <- data.frame(
    sample = sample_id,
    n_cells = ncol(so),
    median_nCount_RNA = median(so$nCount_RNA),
    median_nFeature_RNA = median(so$nFeature_RNA),
    median_percent_mt = median(so$percent.mt),
    mean_percent_mt = mean(so$percent.mt),
    stringsAsFactors = FALSE
  )
  
  message("Saved: ", out_cells)
  summary_df
}

all_summaries <- do.call(
  rbind,
  lapply(seq_len(nrow(samples)), function(i) {
    run_one(samples$sample_id[i], samples$path[i])
  })
)

write.csv(all_summaries, file = "reports/qc_summary.csv", row.names = FALSE)

message("Saved: reports/qc_summary.csv")
print(all_summaries)