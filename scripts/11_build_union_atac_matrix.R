# scripts/11_build_union_atac_matrix.R
# Build a union-peak ATAC matrix for all final QC cells and add it to merged RNA object

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(Matrix)
})

samples <- read.delim("config/samples.tsv", sep = "\t", stringsAsFactors = FALSE)

# load merged RNA object
merged <- readRDS("cache/merged_final_qc_soupx_seurat.rds")

# load union peaks
union_peaks <- readRDS("data/union_peaks.rds")

# keep standard chromosomes only
standard_chr <- paste0("chr", c(1:22, "X", "Y"))
union_peaks <- union_peaks[as.character(seqnames(union_peaks)) %in% standard_chr]

cat("Union peaks used:", length(union_peaks), "\n")

atac_list <- list()
fragment_list <- list()

for (i in seq_len(nrow(samples))) {
  sample_id <- samples$sample_id[i]
  path <- samples$path[i]
  
  message("== Building ATAC matrix: ", sample_id)
  
  in_rds <- file.path("cache", paste0(sample_id, "_final_qc_soupx_seurat.rds"))
  fragments_file <- file.path(path, "atac_fragments.tsv.gz")
  
  if (!file.exists(in_rds)) stop("Missing: ", in_rds)
  if (!file.exists(fragments_file)) stop("Missing: ", fragments_file)
  
  so <- readRDS(in_rds)
  final_cells_original <- colnames(so)
  
  frags <- CreateFragmentObject(
    path = fragments_file,
    cells = final_cells_original
  )
  
  fragment_list[[sample_id]] <- frags
  
  atac_counts <- FeatureMatrix(
    fragments = frags,
    features = union_peaks,
    cells = final_cells_original
  )
  
  # rename cells to match merged object
  colnames(atac_counts) <- paste(sample_id, colnames(atac_counts), sep = "_")
  
  atac_list[[sample_id]] <- atac_counts
  
  message("ATAC matrix built: ", sample_id, " | cells: ", ncol(atac_counts))
}

# combine all sample matrices
combined_atac <- do.call(cbind, atac_list)

cat(
  "Combined ATAC matrix dimensions:",
  nrow(combined_atac), "peaks x", ncol(combined_atac), "cells\n"
)

# ensure same order as merged object
common_cells <- intersect(colnames(merged), colnames(combined_atac))
cat("Shared cells with merged RNA object:", length(common_cells), "\n")

merged <- subset(merged, cells = common_cells)
combined_atac <- combined_atac[, common_cells, drop = FALSE]

# create ATAC assay
chrom_assay <- CreateChromatinAssay(
  counts = combined_atac,
  ranges = union_peaks,
  genome = "hg38",
  fragments = fragment_list
)

merged[["ATAC"]] <- chrom_assay

saveRDS(merged, file = "cache/merged_final_multiome_seurat.rds")

dir.create("reports", showWarnings = FALSE)
report <- data.frame(
  n_cells = ncol(merged),
  n_genes = nrow(merged[["RNA"]]),
  n_peaks = nrow(merged[["ATAC"]]),
  stringsAsFactors = FALSE
)

write.csv(report, "reports/final_multiome_object_summary.csv", row.names = FALSE)

cat("Saved: cache/merged_final_multiome_seurat.rds\n")
cat("Saved: reports/final_multiome_object_summary.csv\n")