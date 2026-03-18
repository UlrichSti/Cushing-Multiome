# scripts/17_aspc_DE_and_DA.R
# Differential expression genes and differential accessibility peaks
# for all ASPCs combined: CS vs KT

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
})

# --------------------
# load ASPC object
# --------------------
aspc <- readRDS("cache/aspc_harmony_wnn.rds")

if (!"condition" %in% colnames(aspc@meta.data)) {
  aspc$patient <- aspc$orig.ident
  aspc$condition <- ifelse(grepl("^CS", aspc$patient), "CS", "KT")
}

dir.create("reports", showWarnings = FALSE)

# --------------------
# RNA DE genes
# --------------------
DefaultAssay(aspc) <- "RNA"
Idents(aspc) <- "condition"

deg_aspc <- FindMarkers(
  aspc,
  ident.1 = "CS",
  ident.2 = "KT",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

deg_aspc$gene <- rownames(deg_aspc)

write.csv(
  deg_aspc,
  "reports/aspc_DEG_all_CS_vs_KT.csv",
  row.names = FALSE
)

# --------------------
# ATAC differential peaks
# --------------------
DefaultAssay(aspc) <- "ATAC"
Idents(aspc) <- "condition"

da_peaks <- FindMarkers(
  aspc,
  ident.1 = "CS",
  ident.2 = "KT",
  logfc.threshold = 0.25,
  min.pct = 0.05
)

da_peaks$peak <- rownames(da_peaks)

write.csv(
  da_peaks,
  "reports/aspc_DApeaks_all_CS_vs_KT.csv",
  row.names = FALSE
)

# --------------------
# annotate peaks to nearest gene
# --------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(aspc[["ATAC"]]) <- annotations

peak_annot <- ClosestFeature(aspc, regions = rownames(da_peaks))
peak_annot$peak <- peak_annot$query_region

da_peaks_annot <- merge(
  da_peaks,
  peak_annot,
  by = "peak",
  all.x = TRUE
)

write.csv(
  da_peaks_annot,
  "reports/aspc_DApeaks_all_CS_vs_KT_annotated.csv",
  row.names = FALSE
)

# --------------------
# simple summary
# --------------------
summary_df <- data.frame(
  comparison = "CS_vs_KT",
  n_deg_total = nrow(deg_aspc),
  n_deg_sig = sum(deg_aspc$p_val_adj < 0.05, na.rm = TRUE),
  n_dapeaks_total = nrow(da_peaks),
  n_dapeaks_sig = sum(da_peaks$p_val_adj < 0.05, na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(
  summary_df,
  "reports/aspc_DE_DA_summary.csv",
  row.names = FALSE
)

cat("Saved: reports/aspc_DEG_all_CS_vs_KT.csv\n")
cat("Saved: reports/aspc_DApeaks_all_CS_vs_KT.csv\n")
cat("Saved: reports/aspc_DApeaks_all_CS_vs_KT_annotated.csv\n")
cat("Saved: reports/aspc_DE_DA_summary.csv\n")