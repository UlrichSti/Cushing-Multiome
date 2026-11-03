# ============================================================
# FIGURE 1 - HUMAN ADIPOSE LANDSCAPE IN CUSHING SYNDROME
# Reproducible manuscript script
#
# Run from the repository root, or set:
#   CUSHING_MULTIOME_PROJECT=/path/to/Cushing-Multiome
#
# INPUT OBJECTS (local, not committed):
#   data/objects/AllCells_RNA_final.rds
#   data/objects/AllCells_WNN_final.rds
#
# CANONICAL MANUSCRIPT OBJECTS CREATED:
#   data/objects/AllCells_RNA_manuscript.rds
#   data/objects/AllCells_WNN_manuscript.rds
#
# THIS SCRIPT DOES NOT RECALCULATE:
#   PCA, Harmony, clustering, RNA UMAP, or WNN UMAP
#
# Locked manuscript annotation changes:
#   RNA cluster 13 -> Myeloid cells
#   RNA cluster 14 -> Adipocytes
# ============================================================

project_dir <- Sys.getenv("CUSHING_MULTIOME_PROJECT", unset = ".")
setwd(project_dir)

if (!file.exists("README.md")) {
  stop("Run from the Cushing-Multiome repository root or set CUSHING_MULTIOME_PROJECT.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
})

dir.create("Figure_1/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("Figure_1/Source_Data", recursive = TRUE, showWarnings = FALSE)
dir.create("Figure_1/Supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("data/objects", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. CREATE / LOCK CANONICAL MANUSCRIPT RNA OBJECT
# ============================================================

rna <- readRDS("data/objects/AllCells_RNA_final.rds")

rna$celltype_RNA_manuscript <- as.character(rna$celltype_RNA)
rna$celltype_RNA_manuscript[rna$seurat_clusters == "13"] <- "Myeloid cells"
rna$celltype_RNA_manuscript[rna$seurat_clusters == "14"] <- "Adipocytes"

celltype_order <- c(
  "Adipocytes",
  "ASPCs",
  "Endothelial cells",
  "Smooth muscle cells",
  "Pericytes",
  "Myeloid cells",
  "Lymphoid cells",
  "Lymphatic endothelial cells",
  "Schwann cells"
)

rna$celltype_RNA_manuscript <- factor(
  rna$celltype_RNA_manuscript,
  levels = celltype_order
)

cat("\n======================================\n")
cat("FINAL RNA MANUSCRIPT ANNOTATION\n")
cat("======================================\n")
print(table(rna$seurat_clusters, rna$celltype_RNA_manuscript))
print(table(rna$celltype_RNA_manuscript))

mapping_table <- data.frame(
  cluster = levels(factor(rna$seurat_clusters)),
  celltype = sapply(
    levels(factor(rna$seurat_clusters)),
    function(x) {
      unique(as.character(
        rna$celltype_RNA_manuscript[rna$seurat_clusters == x]
      ))[1]
    }
  )
)

write.csv(
  mapping_table,
  "Figure_1/Source_Data/Fig1C_RNA_cluster_annotation.csv",
  row.names = FALSE
)

saveRDS(rna, "data/objects/AllCells_RNA_manuscript.rds")

# ============================================================
# 2. FIGURE 1C - GLOBAL RNA UMAP
# ============================================================

celltype_cols <- c(
  "Adipocytes" = "#F2B6B6",
  "ASPCs" = "#F3D19C",
  "Endothelial cells" = "#A8D5BA",
  "Smooth muscle cells" = "#B7C9E2",
  "Pericytes" = "#C8B6E2",
  "Myeloid cells" = "#E7B5D6",
  "Lymphoid cells" = "#A9D7E8",
  "Lymphatic endothelial cells" = "#BFD8A8",
  "Schwann cells" = "#D4C5A9"
)

p_fig1c <- DimPlot(
  rna,
  reduction = "rna.umap",
  group.by = "celltype_RNA_manuscript",
  cols = celltype_cols,
  label = FALSE,
  raster = FALSE,
  pt.size = 0.25
) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

print(p_fig1c)

ggsave("Figure_1/Figures/Fig1C_Global_RNA_UMAP.pdf", p_fig1c, width = 8, height = 6)
ggsave("Figure_1/Figures/Fig1C_Global_RNA_UMAP.png", p_fig1c, width = 8, height = 6, dpi = 600)

# ============================================================
# 3. FIGURE 1D - GLOBAL RNA MARKER DOTPLOT
# ============================================================

DefaultAssay(rna) <- "RNA"

global_markers <- c(
  "ADIPOQ", "PLIN1",
  "PDGFRA", "DCN",
  "PECAM1", "CDH5",
  "MYH11", "MYOCD",
  "RGS5", "PDGFRB",
  "MRC1", "F13A1",
  "IL7R", "CD3E",
  "PROX1", "PTPRE",
  "XKR4"
)

p_fig1d <- DotPlot(
  rna,
  features = global_markers,
  group.by = "celltype_RNA_manuscript",
  cols = c("grey90", "black")
) +
  RotatedAxis() +
  labs(
    x = NULL,
    y = NULL,
    colour = "Average expression",
    size = "Percent expressed"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "plain")
  )

print(p_fig1d)

ggsave("Figure_1/Figures/Fig1D_Global_RNA_DotPlot.pdf", p_fig1d, width = 8, height = 7)
ggsave("Figure_1/Figures/Fig1D_Global_RNA_DotPlot.png", p_fig1d, width = 8, height = 7, dpi = 600)

# ============================================================
# 4. RNA QC SOURCE DATA
# ============================================================

qc_per_cell_rna <- data.frame(
  Cell = Cells(rna),
  Patient = rna$patient,
  Condition = rna$condition,
  nCount_RNA = rna$nCount_RNA,
  nFeature_RNA = rna$nFeature_RNA,
  percent_mt = rna$percent.mt,
  row.names = NULL
)

nuclei_per_patient <- as.data.frame(
  table(Patient = rna$patient, Condition = rna$condition)
)
nuclei_per_patient <- nuclei_per_patient[nuclei_per_patient$Freq > 0, ]
colnames(nuclei_per_patient)[3] <- "Nuclei"

qc_summary_patient_rna <- do.call(
  rbind,
  lapply(
    split(qc_per_cell_rna, qc_per_cell_rna$Patient),
    function(x) {
      data.frame(
        Patient = x$Patient[1],
        Condition = x$Condition[1],
        Nuclei = nrow(x),
        Median_nCount_RNA = median(x$nCount_RNA, na.rm = TRUE),
        Mean_nCount_RNA = mean(x$nCount_RNA, na.rm = TRUE),
        Median_nFeature_RNA = median(x$nFeature_RNA, na.rm = TRUE),
        Mean_nFeature_RNA = mean(x$nFeature_RNA, na.rm = TRUE),
        Median_percent_mt = median(x$percent_mt, na.rm = TRUE),
        Mean_percent_mt = mean(x$percent_mt, na.rm = TRUE)
      )
    }
  )
)
rownames(qc_summary_patient_rna) <- NULL

write.xlsx(
  list(
    Nuclei_per_patient = nuclei_per_patient,
    QC_summary_patient = qc_summary_patient_rna,
    QC_per_nucleus = qc_per_cell_rna
  ),
  file = "Figure_1/Source_Data/FigS1_RNA_QC_SourceData.xlsx",
  overwrite = TRUE
)

# ============================================================
# 5. CREATE / LOCK CANONICAL MANUSCRIPT WNN OBJECT
# ============================================================

wnn <- readRDS("data/objects/AllCells_WNN_final.rds")

stopifnot(identical(sort(Cells(rna)), sort(Cells(wnn))))

wnn$celltype_RNA_manuscript <- as.character(
  rna$celltype_RNA_manuscript[match(Cells(wnn), Cells(rna))]
)
wnn$celltype_RNA_manuscript <- factor(
  wnn$celltype_RNA_manuscript,
  levels = celltype_order
)

stopifnot(!any(is.na(wnn$celltype_RNA_manuscript)))

cat("\n======================================\n")
cat("FINAL WNN MANUSCRIPT ANNOTATION\n")
cat("======================================\n")
print(table(wnn$celltype_RNA_manuscript))

saveRDS(wnn, "data/objects/AllCells_WNN_manuscript.rds")

# ============================================================
# 6. SUPPLEMENTARY - GLOBAL WNN UMAP
# ============================================================

p_wnn_umap <- DimPlot(
  wnn,
  reduction = "wnn.umap",
  group.by = "celltype_RNA_manuscript",
  cols = celltype_cols,
  label = FALSE,
  raster = FALSE,
  pt.size = 0.25
) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

print(p_wnn_umap)

ggsave("Figure_1/Supplementary/FigS1_WNN_UMAP_RNAlabels.pdf", p_wnn_umap, width = 8, height = 6)
ggsave("Figure_1/Supplementary/FigS1_WNN_UMAP_RNAlabels.png", p_wnn_umap, width = 8, height = 6, dpi = 600)

# ============================================================
# 7. MULTIOME QC + WNN SOURCE DATA
# ============================================================

qc_per_nucleus <- data.frame(
  Cell = Cells(wnn),
  Patient = wnn$patient,
  Condition = wnn$condition,
  Celltype = as.character(wnn$celltype_RNA_manuscript),
  nCount_RNA = wnn$nCount_RNA,
  nFeature_RNA = wnn$nFeature_RNA,
  percent_mt = wnn$percent.mt,
  nCount_ATAC = wnn$nCount_ATAC,
  nFeature_ATAC = wnn$nFeature_ATAC,
  TSS_enrichment = wnn$TSS.enrichment,
  nucleosome_signal = wnn$nucleosome_signal,
  RNA_weight = wnn$RNA.weight,
  ATAC_weight = wnn$ATAC.weight,
  row.names = NULL
)

qc_summary_patient_multiome <- qc_per_nucleus %>%
  group_by(Patient, Condition) %>%
  summarise(
    Nuclei = n(),
    Median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    Median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    Median_percent_mt = median(percent_mt, na.rm = TRUE),
    Median_nCount_ATAC = median(nCount_ATAC, na.rm = TRUE),
    Median_nFeature_ATAC = median(nFeature_ATAC, na.rm = TRUE),
    Median_TSS_enrichment = median(TSS_enrichment, na.rm = TRUE),
    Median_nucleosome_signal = median(nucleosome_signal, na.rm = TRUE),
    Mean_RNA_weight = mean(RNA_weight, na.rm = TRUE),
    Mean_ATAC_weight = mean(ATAC_weight, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 8. SUPPLEMENTARY - WNN MODALITY WEIGHTS
# ============================================================

wnn_weight_summary <- qc_per_nucleus %>%
  group_by(Celltype) %>%
  summarise(
    Nuclei = n(),
    Mean_RNA_weight = mean(RNA_weight, na.rm = TRUE),
    Mean_ATAC_weight = mean(ATAC_weight, na.rm = TRUE),
    Median_RNA_weight = median(RNA_weight, na.rm = TRUE),
    Median_ATAC_weight = median(ATAC_weight, na.rm = TRUE),
    .groups = "drop"
  )

wnn_weight_summary$Celltype <- factor(
  wnn_weight_summary$Celltype,
  levels = celltype_order
)

wnn_weight_long <- wnn_weight_summary %>%
  select(Celltype, Mean_RNA_weight, Mean_ATAC_weight) %>%
  pivot_longer(
    cols = c(Mean_RNA_weight, Mean_ATAC_weight),
    names_to = "Modality",
    values_to = "Mean_weight"
  )

wnn_weight_long$Modality <- factor(
  wnn_weight_long$Modality,
  levels = c("Mean_RNA_weight", "Mean_ATAC_weight"),
  labels = c("RNA", "ATAC")
)

p_wnn_weight <- ggplot(
  wnn_weight_long,
  aes(x = Celltype, y = Mean_weight, fill = Modality)
) +
  geom_col(width = 0.75) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Mean WNN modality weight", fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_wnn_weight)

ggsave("Figure_1/Supplementary/FigS1_WNN_ModalityWeights.pdf", p_wnn_weight, width = 8, height = 5)
ggsave("Figure_1/Supplementary/FigS1_WNN_ModalityWeights.png", p_wnn_weight, width = 8, height = 5, dpi = 600)

# ============================================================
# 9. WNN UMAP SOURCE DATA
# ============================================================

wnn_embeddings <- Embeddings(wnn, "wnn.umap")

wnn_umap_source <- data.frame(
  Cell = rownames(wnn_embeddings),
  WNN_UMAP_1 = wnn_embeddings[, 1],
  WNN_UMAP_2 = wnn_embeddings[, 2],
  row.names = NULL
)

wnn_umap_source$Patient <- wnn$patient[
  match(wnn_umap_source$Cell, Cells(wnn))
]
wnn_umap_source$Condition <- wnn$condition[
  match(wnn_umap_source$Cell, Cells(wnn))
]
wnn_umap_source$Celltype <- as.character(
  wnn$celltype_RNA_manuscript[
    match(wnn_umap_source$Cell, Cells(wnn))
  ]
)

write.xlsx(
  list(
    QC_per_nucleus = qc_per_nucleus,
    QC_summary_patient = qc_summary_patient_multiome,
    WNN_weights_celltype = wnn_weight_summary,
    WNN_weights_per_nucleus = qc_per_nucleus %>%
      select(Cell, Patient, Condition, Celltype, RNA_weight, ATAC_weight),
    WNN_UMAP = wnn_umap_source
  ),
  file = "Figure_1/Source_Data/FigS1_Multiome_QC_SourceData.xlsx",
  overwrite = TRUE
)

# ============================================================
# 10. FIGURE 1E - BROAD CELL-TYPE PROPORTIONS
# Biological replicate = patient
# Control vs Cushing = Welch t-test, unadjusted
# ============================================================

cell_counts <- rna@meta.data %>%
  mutate(
    Patient = patient,
    Condition = condition,
    Celltype = as.character(celltype_RNA_manuscript)
  ) %>%
  count(Patient, Condition, Celltype, name = "N_cells")

all_patients <- unique(rna@meta.data[, c("patient", "condition")])
colnames(all_patients) <- c("Patient", "Condition")
all_celltypes <- levels(rna$celltype_RNA_manuscript)

cell_counts <- expand_grid(
  Patient = unique(all_patients$Patient),
  Celltype = all_celltypes
) %>%
  left_join(all_patients, by = "Patient") %>%
  left_join(cell_counts, by = c("Patient", "Condition", "Celltype")) %>%
  mutate(N_cells = replace_na(N_cells, 0))

cell_proportions <- cell_counts %>%
  group_by(Patient, Condition) %>%
  mutate(
    Total_cells = sum(N_cells),
    Proportion = N_cells / Total_cells,
    Proportion_percent = Proportion * 100
  ) %>%
  ungroup()

stats <- cell_proportions %>%
  group_by(Celltype) %>%
  summarise(
    p_value = tryCatch(
      t.test(Proportion_percent ~ Condition, data = cur_data(), var.equal = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

condition_summary <- cell_proportions %>%
  group_by(Celltype, Condition) %>%
  summarise(
    N_patients = n(),
    Mean_percent = mean(Proportion_percent),
    SD_percent = sd(Proportion_percent),
    Median_percent = median(Proportion_percent),
    .groups = "drop"
  )

summary_wide <- condition_summary %>%
  select(Celltype, Condition, N_patients, Mean_percent, SD_percent, Median_percent) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(N_patients, Mean_percent, SD_percent, Median_percent)
  )

final_table <- summary_wide %>%
  left_join(stats, by = "Celltype") %>%
  arrange(match(Celltype, levels(rna$celltype_RNA_manuscript)))

write.xlsx(
  list(
    Patient_proportions = cell_proportions,
    Condition_summary = condition_summary,
    Statistics = final_table
  ),
  file = "Figure_1/Source_Data/Fig1E_BroadCellType_Proportions_Statistics.xlsx",
  overwrite = TRUE
)

cat("\n======================================\n")
cat("FIGURE 1 ANALYSIS COMPLETE\n")
cat("======================================\n")
print(final_table)
