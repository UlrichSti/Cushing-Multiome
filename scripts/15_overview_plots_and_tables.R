# scripts/15_overview_plots_and_tables.R
# Generate overview dotplots, UMAPs, fraction tables, and barplots
# for the main analysis object (without adrenal tissue), plus comparison
# dotplot with/without adrenal tissue.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

multiome <- readRDS("cache/merged_final_multiome_harmony_wnn_annotated.rds")
multiome_main <- readRDS("cache/merged_final_multiome_harmony_wnn_main_analysis.rds")

DefaultAssay(multiome) <- "RNA"
DefaultAssay(multiome_main) <- "RNA"

# ---- color palette ----
celltype_colors2 <- c(
  "Adipocytes" = "#E69F00",
  "ASPC" = "#56B4E9",
  "Endothelial cells" = "#009E73",
  "Smooth muscle cells" = "#D55E00",
  "Pericytes" = "#CC79A7",
  "Lymphoid immune cells" = "#0072B2",
  "Myeloid immune cells" = "#B22222",
  "Schwann cells" = "#A6761D",
  "LEC" = "#00BFC4",
  "Adrenal tissue" = "#7F7F7F"
)

order_with_adrenal <- c(
  "Adipocytes",
  "ASPC",
  "Endothelial cells",
  "Smooth muscle cells",
  "Pericytes",
  "Lymphoid immune cells",
  "Myeloid immune cells",
  "Schwann cells",
  "LEC",
  "Adrenal tissue"
)

order_no_adrenal <- c(
  "Adipocytes",
  "ASPC",
  "Endothelial cells",
  "Smooth muscle cells",
  "Pericytes",
  "Lymphoid immune cells",
  "Myeloid immune cells",
  "Schwann cells",
  "LEC"
)

multiome$broad_celltype2 <- factor(multiome$broad_celltype2, levels = order_with_adrenal)
multiome_main$broad_celltype2 <- factor(multiome_main$broad_celltype2, levels = order_no_adrenal)

# ---- dotplots ----
markers_with_adrenal <- c(
  "ADIPOQ", "PLIN1",
  "PDGFRA", "DCN",
  "PECAM1", "CDH5",
  "MYH11", "MYOCD",
  "RGS5", "PDGFRB",
  "IL7R", "CD3E",
  "MRC1", "F13A1",
  "XKR4",
  "PROX1", "PTPRE",
  "CYP11B1", "FDX1"
)
markers_with_adrenal <- markers_with_adrenal[markers_with_adrenal %in% rownames(multiome[["RNA"]])]

p_dot_with_adrenal <- DotPlot(
  multiome,
  features = markers_with_adrenal,
  assay = "RNA",
  group.by = "broad_celltype2"
) + RotatedAxis() + ggtitle("DotPlot with adrenal tissue")

markers_no_adrenal <- c(
  "ADIPOQ", "PLIN1",
  "PDGFRA", "DCN",
  "PECAM1", "CDH5",
  "MYH11", "MYOCD",
  "RGS5", "PDGFRB",
  "IL7R", "CD3E",
  "MRC1", "F13A1",
  "XKR4",
  "PROX1", "PTPRE"
)
markers_no_adrenal <- markers_no_adrenal[markers_no_adrenal %in% rownames(multiome_main[["RNA"]])]

p_dot_no_adrenal <- DotPlot(
  multiome_main,
  features = markers_no_adrenal,
  assay = "RNA",
  group.by = "broad_celltype2"
) + RotatedAxis() + ggtitle("DotPlot without adrenal tissue")

# ---- UMAPs ----
p_umap_with_adrenal <- DimPlot(
  multiome,
  reduction = "wnn.umap",
  group.by = "broad_celltype2",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3
) +
  scale_color_manual(values = celltype_colors2) +
  theme_classic() +
  ggtitle("UMAP with adrenal tissue")

p_umap_no_adrenal <- DimPlot(
  multiome_main,
  reduction = "wnn.umap",
  group.by = "broad_celltype2",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3
) +
  scale_color_manual(values = celltype_colors2[names(celltype_colors2) != "Adrenal tissue"]) +
  theme_classic() +
  ggtitle("UMAP without adrenal tissue")

p_umap_condition <- DimPlot(
  multiome_main,
  reduction = "wnn.umap",
  group.by = "broad_celltype2",
  split.by = "condition",
  pt.size = 0.3,
  ncol = 2
) +
  scale_color_manual(values = celltype_colors2[names(celltype_colors2) != "Adrenal tissue"]) +
  theme_classic() +
  ggtitle("UMAP by condition (no adrenal)")

p_umap_patient <- DimPlot(
  multiome_main,
  reduction = "wnn.umap",
  group.by = "broad_celltype2",
  split.by = "patient",
  pt.size = 0.3,
  ncol = 4
) +
  scale_color_manual(values = celltype_colors2[names(celltype_colors2) != "Adrenal tissue"]) +
  theme_classic() +
  ggtitle("UMAP by patient (no adrenal)")

# ---- fraction tables ----
fractions_condition <- multiome_main@meta.data %>%
  count(condition, broad_celltype2) %>%
  group_by(condition) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()

fractions_patient <- multiome_main@meta.data %>%
  count(patient, broad_celltype2) %>%
  group_by(patient) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup()

# ---- barplots ----
p_bar_condition <- ggplot(
  fractions_condition,
  aes(x = condition, y = fraction, fill = broad_celltype2)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_colors2[names(celltype_colors2) != "Adrenal tissue"], drop = FALSE) +
  theme_classic() +
  ylab("Fraction") +
  xlab("Condition")

p_bar_patient <- ggplot(
  fractions_patient,
  aes(x = patient, y = fraction, fill = broad_celltype2)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_colors2[names(celltype_colors2) != "Adrenal tissue"], drop = FALSE) +
  theme_classic() +
  ylab("Fraction") +
  xlab("Patient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- save outputs ----
dir.create("reports", showWarnings = FALSE)

write.csv(fractions_condition, "reports/fractions_by_condition_no_adrenal.csv", row.names = FALSE)
write.csv(fractions_patient, "reports/fractions_by_patient_no_adrenal.csv", row.names = FALSE)

ggsave("reports/dotplot_broad_celltypes_with_adrenal.png", p_dot_with_adrenal, width = 10, height = 5, dpi = 300)
ggsave("reports/dotplot_broad_celltypes_without_adrenal.png", p_dot_no_adrenal, width = 9, height = 5, dpi = 300)

ggsave("reports/umap_broad_celltypes2_with_adrenal.png", p_umap_with_adrenal, width = 8, height = 6, dpi = 300)
ggsave("reports/umap_broad_celltypes2_without_adrenal.png", p_umap_no_adrenal, width = 8, height = 6, dpi = 300)
ggsave("reports/umap_by_condition_no_adrenal.png", p_umap_condition, width = 12, height = 6, dpi = 300)
ggsave("reports/umap_by_patient_no_adrenal.png", p_umap_patient, width = 16, height = 10, dpi = 300)

ggsave("reports/barplot_cell_fractions_by_condition_no_adrenal.png", p_bar_condition, width = 7, height = 5, dpi = 300)
ggsave("reports/barplot_cell_fractions_by_patient_no_adrenal.png", p_bar_patient, width = 10, height = 6, dpi = 300)

ggsave("reports/umap_by_condition_no_adrenal.pdf", p_umap_condition, width = 12, height = 6)
ggsave("reports/umap_by_patient_no_adrenal.pdf", p_umap_patient, width = 16, height = 10)
ggsave("reports/barplot_cell_fractions_by_condition_no_adrenal.pdf", p_bar_condition, width = 7, height = 5)
ggsave("reports/barplot_cell_fractions_by_patient_no_adrenal.pdf", p_bar_patient, width = 10, height = 6)

cat("Saved annotation overview tables and plots to reports/\n")