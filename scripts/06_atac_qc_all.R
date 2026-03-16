# scripts/06_atac_qc_all.R
# ATAC QC for all snMultiome samples
# Input per sample:
#   - atac_peaks.bed
#   - atac_fragments.tsv.gz
#   - filtered_feature_bc_matrix.h5
# Output:
#   - reports/<sample>_atac_qc_cells.csv
#   - reports/atac_qc_summary.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
})

samples <- read.delim("config/samples.tsv", sep = "\t", stringsAsFactors = FALSE)

pick_modality <- function(x, wanted, sample_id) {
  if (!is.list(x)) stop(sample_id, ": H5 does not contain multiple modalities.")
  if (wanted %in% names(x)) return(x[[wanted]])
  stop(sample_id, ": modality '", wanted, "' not found. Available: ",
       paste(names(x), collapse = ", "))
}

run_one <- function(sample_id, path) {
  message("== ATAC QC: ", sample_id)
  
  peaks_bed <- file.path(path, "atac_peaks.bed")
  fragments_file <- file.path(path, "atac_fragments.tsv.gz")
  filtered_h5 <- file.path(path, "filtered_feature_bc_matrix.h5")
  
  if (!file.exists(peaks_bed)) stop("Missing: ", peaks_bed)
  if (!file.exists(fragments_file)) stop("Missing: ", fragments_file)
  if (!file.exists(filtered_h5)) stop("Missing: ", filtered_h5)
  
  # read peak ranges
  peaks_df <- read.table(peaks_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(peaks_df) <- c("chr", "start", "end")
  peaks_gr <- makeGRangesFromDataFrame(peaks_df)
  
  # read filtered cells from H5 and extract ATAC counts
  x <- Read10X_h5(filtered_h5)
  atac <- pick_modality(x, "Peaks", sample_id)
  
  # create chromatin assay
  chrom_assay <- CreateChromatinAssay(
    counts = atac,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = fragments_file,
    ranges = peaks_gr,
    min.cells = 1,
    min.features = 1
  )
  
  so <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    project = sample_id
  )
  
  # add annotation required for TSS enrichment
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "hg38"
  Annotation(so[["ATAC"]]) <- annotation
  
  # ATAC QC metrics
  so <- NucleosomeSignal(so)
  so <- TSSEnrichment(so, fast = FALSE)
  
  qc_df <- data.frame(
    cell = colnames(so),
    sample = sample_id,
    nCount_ATAC = so$nCount_ATAC,
    nFeature_ATAC = so$nFeature_ATAC,
    nucleosome_signal = so$nucleosome_signal,
    TSS.enrichment = so$TSS.enrichment,
    stringsAsFactors = FALSE
  )
  
  dir.create("reports", showWarnings = FALSE)
  
  out_cells <- file.path("reports", paste0(sample_id, "_atac_qc_cells.csv"))
  write.csv(qc_df, out_cells, row.names = FALSE)
  
  summary_df <- data.frame(
    sample = sample_id,
    n_cells = ncol(so),
    median_nCount_ATAC = median(so$nCount_ATAC, na.rm = TRUE),
    median_nFeature_ATAC = median(so$nFeature_ATAC, na.rm = TRUE),
    median_nucleosome_signal = median(so$nucleosome_signal, na.rm = TRUE),
    median_TSS_enrichment = median(so$TSS.enrichment, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  message("Saved: ", out_cells)
  summary_df
}

all_summaries <- do.call(
  rbind,
  lapply(seq_len(nrow(samples)), function(i) {
    run_one(samples$sample_id[i], samples$path[i])
  })
)

dir.create("reports", showWarnings = FALSE)
write.csv(all_summaries, file = "reports/atac_qc_summary.csv", row.names = FALSE)

message("Saved: reports/atac_qc_summary.csv")
print(all_summaries)