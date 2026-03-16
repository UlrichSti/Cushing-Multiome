# scripts/09_merge_final_soupx_all.R
# Merge all final QC + SoupX corrected RNA Seurat objects
# Input:
#   cache/<sample>_final_qc_soupx_seurat.rds
# Output:
#   cache/merged_final_qc_soupx_seurat.rds
#   reports/merge_summary.csv

suppressPackageStartupMessages({
  library(Seurat)
})

samples <- read.delim("config/samples.tsv", sep = "\t", stringsAsFactors = FALSE)

objs <- lapply(samples$sample_id, function(sample_id) {
  in_rds <- file.path("cache", paste0(sample_id, "_final_qc_soupx_seurat.rds"))
  if (!file.exists(in_rds)) stop("Missing: ", in_rds)
  
  so <- readRDS(in_rds)
  so$sample <- sample_id
  so <- RenameCells(so, add.cell.id = sample_id)
  so
})

names(objs) <- samples$sample_id

merged <- objs[[1]]
if (length(objs) > 1) {
  for (i in 2:length(objs)) {
    merged <- merge(merged, y = objs[[i]])
  }
}

dir.create("cache", showWarnings = FALSE)
dir.create("reports", showWarnings = FALSE)

saveRDS(merged, file = "cache/merged_final_qc_soupx_seurat.rds")

merge_summary <- data.frame(
  sample = samples$sample_id,
  n_cells = sapply(objs, ncol),
  stringsAsFactors = FALSE
)

write.csv(merge_summary, "reports/merge_summary.csv", row.names = FALSE)

cat("Saved: cache/merged_final_qc_soupx_seurat.rds\n")
cat("Saved: reports/merge_summary.csv\n")
print(merge_summary)
cat("Merged object cells:", ncol(merged), "\n")
cat("Merged object genes:", nrow(merged), "\n")