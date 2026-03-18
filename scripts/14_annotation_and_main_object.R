# scripts/14_annotation_and_main_object.R
# Add broad cell type annotations, define main analysis object without adrenal tissue,
# and save both annotated objects.

suppressPackageStartupMessages({
  library(Seurat)
})

multiome <- readRDS("cache/merged_final_multiome_harmony_wnn.rds")

# ---- cluster-to-celltype mapping ----
celltype_map <- c(
  "0"  = "Adipocytes",
  "1"  = "Adipocytes",
  "2"  = "Macrophages",
  "3"  = "ASPC",
  "4"  = "Adipocytes",
  "5"  = "Adipocytes",
  "6"  = "ASPC",
  "7"  = "Endothelial cells",
  "8"  = "Smooth muscle cells",
  "9"  = "Adrenal tissue",
  "10" = "Endothelial cells",
  "11" = "Pericytes",
  "12" = "T cells",
  "13" = "Adipocytes",
  "14" = "Endothelial cells",
  "15" = "Adipocytes",
  "16" = "ASPC",
  "17" = "B cells",
  "18" = "Adipocytes",
  "19" = "Adipocytes",
  "20" = "LEC",
  "21" = "Adipocytes",
  "22" = "Adipocytes",
  "23" = "Adipocytes",
  "24" = "Adipocytes",
  "25" = "Adipocytes",
  "26" = "Adipocytes",
  "27" = "Adipocytes",
  "28" = "Adipocytes",
  "29" = "Adipocytes",
  "30" = "Schwann cells",
  "31" = "Adipocytes",
  "32" = "Macrophages",
  "33" = "Macrophages"
)

multiome$broad_celltype <- unname(celltype_map[as.character(Idents(multiome))])

# ---- simplified broad cell types ----
multiome$broad_celltype2 <- multiome$broad_celltype
multiome$broad_celltype2[multiome$broad_celltype2 %in% c("Macrophages", "B cells")] <- "Myeloid immune cells"
multiome$broad_celltype2[multiome$broad_celltype2 %in% c("T cells")] <- "Lymphoid immune cells"

# ---- patient / condition ----
multiome$patient <- multiome$orig.ident
multiome$condition <- ifelse(grepl("^CS", multiome$patient), "CS", "KT")

# ---- flag adrenal contamination ----
multiome$use_for_main_analysis <- TRUE
multiome$use_for_main_analysis[multiome$broad_celltype == "Adrenal tissue"] <- FALSE

# ---- create main analysis object without adrenal tissue ----
multiome_main <- subset(multiome, subset = use_for_main_analysis)

# ---- save objects ----
dir.create("cache", showWarnings = FALSE)
saveRDS(multiome, "cache/merged_final_multiome_harmony_wnn_annotated.rds")
saveRDS(multiome_main, "cache/merged_final_multiome_harmony_wnn_main_analysis.rds")

# ---- small reports ----
dir.create("reports", showWarnings = FALSE)

write.csv(
  as.data.frame(table(multiome$broad_celltype2), stringsAsFactors = FALSE),
  "reports/broad_celltype2_counts_with_adrenal.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(table(multiome_main$broad_celltype2), stringsAsFactors = FALSE),
  "reports/broad_celltype2_counts_no_adrenal.csv",
  row.names = FALSE
)

cat("Saved: cache/merged_final_multiome_harmony_wnn_annotated.rds\n")
cat("Saved: cache/merged_final_multiome_harmony_wnn_main_analysis.rds\n")
cat("Saved: reports/broad_celltype2_counts_with_adrenal.csv\n")
cat("Saved: reports/broad_celltype2_counts_no_adrenal.csv\n")