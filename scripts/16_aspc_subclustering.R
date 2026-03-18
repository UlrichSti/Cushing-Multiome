# scripts/16_aspc_subclustering.R
# Subset ASPCs from the main object and perform Harmony + WNN subclustering

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(harmony)
  library(ggplot2)
  library(dplyr)
})

# --------------------
# load main object
# --------------------
obj <- readRDS("cache/merged_final_multiome_harmony_wnn_main_analysis.rds")

if (!"broad_celltype2" %in% colnames(obj@meta.data)) {
  stop("broad_celltype2 not found in metadata.")
}
if (!"orig.ident" %in% colnames(obj@meta.data)) {
  stop("orig.ident not found in metadata.")
}

obj$patient <- obj$orig.ident
obj$condition <- ifelse(grepl("^CS", obj$patient), "CS", "KT")

# --------------------
# subset ASPCs
# --------------------
aspc <- subset(obj, subset = broad_celltype2 == "ASPC")

cat("ASPC cells:", ncol(aspc), "\n")

# --------------------
# small count summaries
# --------------------
dir.create("reports", showWarnings = FALSE)

counts_by_patient <- aspc@meta.data %>%
  count(patient) %>%
  arrange(desc(n))

counts_by_condition <- aspc@meta.data %>%
  count(condition) %>%
  arrange(desc(n))

write.csv(counts_by_patient, "reports/aspc_counts_by_patient.csv", row.names = FALSE)
write.csv(counts_by_condition, "reports/aspc_counts_by_condition.csv", row.names = FALSE)

# --------------------
# RNA reprocessing
# --------------------
DefaultAssay(aspc) <- "RNA"

aspc <- NormalizeData(aspc, verbose = FALSE)
aspc <- FindVariableFeatures(aspc, verbose = FALSE)
aspc <- ScaleData(aspc, verbose = FALSE)
aspc <- RunPCA(aspc, npcs = 30, verbose = FALSE)

aspc <- RunHarmony(
  object = aspc,
  group.by.vars = "orig.ident",
  reduction = "pca",
  assay.use = "RNA",
  project.dim = FALSE,
  reduction.save = "harmony_rna",
  verbose = TRUE
)

# --------------------
# ATAC reprocessing
# --------------------
DefaultAssay(aspc) <- "ATAC"

aspc <- RunTFIDF(aspc)
aspc <- FindTopFeatures(aspc, min.cutoff = "q0")
aspc <- RunSVD(aspc)

aspc <- RunHarmony(
  object = aspc,
  group.by.vars = "orig.ident",
  reduction = "lsi",
  assay.use = "ATAC",
  project.dim = FALSE,
  reduction.save = "harmony_atac",
  verbose = TRUE
)

# --------------------
# WNN
# --------------------
aspc <- FindMultiModalNeighbors(
  aspc,
  reduction.list = list("harmony_rna", "harmony_atac"),
  dims.list = list(1:30, 2:30),
  verbose = TRUE
)

aspc <- RunUMAP(
  aspc,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  verbose = TRUE
)

aspc <- FindClusters(
  aspc,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.4,
  verbose = TRUE
)

# --------------------
# save object
# --------------------
dir.create("cache", showWarnings = FALSE)
saveRDS(aspc, file = "cache/aspc_harmony_wnn.rds")

# --------------------
# reports
# --------------------
cluster_sizes <- data.frame(
  cluster = names(table(Idents(aspc))),
  n_cells = as.integer(table(Idents(aspc))),
  stringsAsFactors = FALSE
)

write.csv(cluster_sizes, "reports/aspc_cluster_sizes.csv", row.names = FALSE)

summary_df <- data.frame(
  n_cells = ncol(aspc),
  n_genes = nrow(aspc[["RNA"]]),
  n_peaks = nrow(aspc[["ATAC"]]),
  n_clusters = length(unique(Idents(aspc))),
  stringsAsFactors = FALSE
)

write.csv(summary_df, "reports/aspc_summary.csv", row.names = FALSE)

# --------------------
# overview plots
# --------------------
p_clusters <- DimPlot(
  aspc,
  reduction = "wnn.umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + theme_classic() + ggtitle("ASPC subclusters")

p_patient <- DimPlot(
  aspc,
  reduction = "wnn.umap",
  group.by = "orig.ident"
) + theme_classic() + ggtitle("ASPC by patient")

p_condition <- DimPlot(
  aspc,
  reduction = "wnn.umap",
  group.by = "condition"
) + theme_classic() + ggtitle("ASPC by condition")

ggsave("reports/aspc_umap_clusters.png", p_clusters, width = 8, height = 6, dpi = 300)
ggsave("reports/aspc_umap_by_patient.png", p_patient, width = 8, height = 6, dpi = 300)
ggsave("reports/aspc_umap_by_condition.png", p_condition, width = 8, height = 6, dpi = 300)

cat("Saved: cache/aspc_harmony_wnn.rds\n")
cat("Saved: reports/aspc_summary.csv\n")
cat("Saved: reports/aspc_cluster_sizes.csv\n")