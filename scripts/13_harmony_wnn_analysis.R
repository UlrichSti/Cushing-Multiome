# scripts/13_harmony_wnn_analysis.R
# Harmony integration + WNN on final multiome object

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(harmony)
})

multiome <- readRDS("cache/merged_final_multiome_seurat.rds")

# make sure sample info exists
if (!"orig.ident" %in% colnames(multiome@meta.data)) {
  stop("orig.ident not found in metadata.")
}

# --------------------
# RNA
# --------------------
DefaultAssay(multiome) <- "RNA"

multiome <- NormalizeData(multiome, verbose = FALSE)
multiome <- FindVariableFeatures(multiome, verbose = FALSE)
multiome <- ScaleData(multiome, verbose = FALSE)
multiome <- RunPCA(multiome, npcs = 30, verbose = FALSE)

multiome <- RunHarmony(
  object = multiome,
  group.by.vars = "orig.ident",
  reduction = "pca",
  assay.use = "RNA",
  project.dim = FALSE,
  reduction.save = "harmony_rna",
  verbose = TRUE
)

# --------------------
# ATAC
# --------------------
DefaultAssay(multiome) <- "ATAC"

multiome <- RunTFIDF(multiome)
multiome <- FindTopFeatures(multiome, min.cutoff = "q0")
multiome <- RunSVD(multiome)

multiome <- RunHarmony(
  object = multiome,
  group.by.vars = "orig.ident",
  reduction = "lsi",
  assay.use = "ATAC",
  project.dim = FALSE,
  reduction.save = "harmony_atac",
  verbose = TRUE
)

# --------------------
# WNN on harmony reductions
# --------------------
multiome <- FindMultiModalNeighbors(
  multiome,
  reduction.list = list("harmony_rna", "harmony_atac"),
  dims.list = list(1:30, 2:30),
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

saveRDS(multiome, file = "cache/merged_final_multiome_harmony_wnn.rds")

dir.create("reports", showWarnings = FALSE)
report <- data.frame(
  n_cells = ncol(multiome),
  n_genes = nrow(multiome[["RNA"]]),
  n_peaks = nrow(multiome[["ATAC"]]),
  n_clusters = length(unique(Idents(multiome))),
  stringsAsFactors = FALSE
)

write.csv(report, "reports/harmony_wnn_summary.csv", row.names = FALSE)

cat("Saved: cache/merged_final_multiome_harmony_wnn.rds\n")
cat("Saved: reports/harmony_wnn_summary.csv\n")