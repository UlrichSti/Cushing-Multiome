# scripts/10_union_peaks.R
# Create a union peak set across all samples

suppressPackageStartupMessages({
  library(GenomicRanges)
})

samples <- read.delim("config/samples.tsv", sep = "\t", stringsAsFactors = FALSE)

all_peaks <- list()

for (i in seq_len(nrow(samples))) {
  sample_id <- samples$sample_id[i]
  path <- samples$path[i]
  
  message("== Reading peaks: ", sample_id)
  
  peaks_file <- file.path(path, "atac_peaks.bed")
  if (!file.exists(peaks_file)) stop("Missing: ", peaks_file)
  
  peaks_df <- read.table(peaks_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(peaks_df) <- c("chr", "start", "end")
  
  gr <- makeGRangesFromDataFrame(peaks_df)
  all_peaks[[sample_id]] <- gr
}

# combine all peaks into one GRanges object
combined <- Reduce(c, all_peaks)

cat("Total peaks before merge:", length(combined), "\n")

# merge overlapping peaks
union_peaks <- reduce(combined)

cat("Total peaks after reduce:", length(union_peaks), "\n")

# optional: filter extreme peak sizes (recommended)
peak_widths <- width(union_peaks)

union_peaks <- union_peaks[
  peak_widths > 20 & peak_widths < 10000
]

cat("Total peaks after filtering:", length(union_peaks), "\n")

# save
dir.create("data", showWarnings = FALSE)

saveRDS(union_peaks, file = "data/union_peaks.rds")

# also save as BED
bed_df <- data.frame(
  chr = seqnames(union_peaks),
  start = start(union_peaks),
  end = end(union_peaks)
)

write.table(
  bed_df,
  file = "data/union_peaks.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("Saved: data/union_peaks.rds\n")
cat("Saved: data/union_peaks.bed\n")