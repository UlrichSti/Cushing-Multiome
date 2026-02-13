# scripts/01_precluster_all.R
# Load filtered RNA for each sample, quick cluster/UMAP (for SoupX), save to cache/

suppressPackageStartupMessages({
  library(Seurat)
})

samples <- read.delim("config/samples.tsv", sep = "\t", stringsAsFactors = FALSE)

pick_rna <- function(x) {
  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) return(x[["Gene Expression"]])
    return(x[[1]])
  }
  x
}

run_one <- function(sample_id, path) {
  message("== Precluster: ", sample_id)
  
  filtered_h5 <- file.path(path, "filtered_feature_bc_matrix.h5")
  if (!file.exists(filtered_h5)) stop("Missing: ", filtered_h5)
  
  x <- Read10X_h5(filtered_h5)
  rna <- pick_rna(x)
  
  so <- CreateSeuratObject(counts = rna, project = sample_id)
  
  so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, verbose = FALSE)
  so <- ScaleData(so, verbose = FALSE)
  so <- RunPCA(so, npcs = 30, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
  so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
  so <- RunUMAP(so, dims = 1:30, verbose = FALSE)
  
  dir.create("cache", showWarnings = FALSE)
  out_rds <- file.path("cache", paste0(sample_id, "_precluster_seurat.rds"))
  saveRDS(so, out_rds)
  message("Saved: ", out_rds)
}

for (i in seq_len(nrow(samples))) {
  run_one(samples$sample_id[i], samples$path[i])
}
