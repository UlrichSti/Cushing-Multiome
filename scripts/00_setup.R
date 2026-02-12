# scripts/00_setup.R
# Minimal: read sample sheet and check files exist

samples <- read.delim('config/samples.tsv', sep='\t', stringsAsFactors=FALSE)

required <- c(
  'atac_fragments.tsv.gz',
  'atac_fragments.tsv.gz.tbi',
  'atac_peaks.bed',
  'filtered_feature_bc_matrix.h5',
  'raw_feature_bc_matrix.h5'
)

check_sample <- function(s) {
  files <- file.path(s, required)
  ok <- file.exists(files)
  data.frame(sample_id=rep(s, length(required)), file=required, exists=ok)
}

out <- do.call(rbind, lapply(samples, check_sample))
print(out)

dir.create('results', showWarnings=FALSE)
write.csv(out, file='results/file_check.csv', row.names=FALSE)
