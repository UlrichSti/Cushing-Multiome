# scripts/12_wnn_analysis.R
# WNN analysis on final multiome object

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
})

multiome <- readRDS("cache/merged_final_multiome_seurat.rds")

# --------------------
# RNA
# --------------------
DefaultAssay(multiome) <- "RNA"

multiome <- NormalizeData(multiome, verbose = FALSE)
multiome <- FindVariableFeatures(multiome, verbose = FALSE)
multiome <- ScaleData(multiome, verbose = FALSE)
multiome <- RunPCA(multiome, npcs = 30, verbose = FALSE)

# --------------------
# ATAC
# --------------------
DefaultAssay(multiome) <- "ATAC"

multiome <- RunTFIDF(multiome)
multiome <- FindTopFeatures(multiome, min.cutoff = "q0")
multiome <- RunSVD(multiome)

# --------------------
# WNN
# --------------------
multiome <- FindMultiModalNeighbors(
  multiome,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

multiome <- RunUMAP(
  multiome,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  verbose = TRUE
)

multiome <- FindClusters(
  multiome,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.5,
  verbose = TRUE
)

# --------------------
# Save
# --------------------
saveRDS(multiome, file = "cache/merged_final_multiome_wnn.rds")

dir.create("reports", showWarnings = FALSE)

report <- data.frame(
  n_cells = ncol(multiome),
  n_genes = nrow(multiome[["RNA"]]),
  n_peaks = nrow(multiome[["ATAC"]]),
  n_clusters = length(unique(Idents(multiome))),
  stringsAsFactors = FALSE
)

write.csv(report, "reports/wnn_summary.csv", row.names = FALSE)

cat("Saved: cache/merged_final_multiome_wnn.rds\n")
cat("Saved: reports/wnn_summary.csv\n")