# scripts/01_cs3_precluster.R
# Goal: load CS3 RNA (filtered), do quick clustering for SoupX, save Seurat object to cache/

suppressPackageStartupMessages({
  library(Seurat)
})

sample <- "CS3"
filtered_h5 <- file.path(sample, "filtered_feature_bc_matrix.h5")

x <- Read10X_h5(filtered_h5)

# 10x Multiome H5 usually contains a list, pick Gene Expression (RNA)
rna <- if (is.list(x)) {
  if ("Gene Expression" %in% names(x)) x[["Gene Expression"]] else x[[1]]
} else {
  x
}

so <- CreateSeuratObject(counts = rna, project = sample)

# quick standard preprocessing (enough for SoupX)
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
so <- RunUMAP(so, dims = 1:30, verbose = FALSE)

dir.create("cache", showWarnings = FALSE)
saveRDS(so, file = file.path("cache", paste0(sample, "_precluster_seurat.rds")))

cat("Saved:", file.path("cache", paste0(sample, "_precluster_seurat.rds")), "\n")
