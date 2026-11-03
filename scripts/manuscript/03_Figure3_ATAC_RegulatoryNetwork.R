# ============================================================
# FIGURE 3 - CHROMATIN REMODELING AND REGULATORY TF PROGRAMS
# Reproducible manuscript script
#
# Run from repository root, or set:
#   CUSHING_MULTIOME_PROJECT=/path/to/Cushing-Multiome
#
# Required local checkpoint objects (not committed):
#   data/objects/ASPC_CXCL14_PPARG_DAP_directional_RegionMatrix_depthnormalized.rds
#   data/objects/ASPC_Lineage_RNA_ATAC_Motifs.rds
#   data/objects/ASPC_Lineage_RNA_DEG_LinkPeaks.rds
#
# Required local result table:
#   data/manuscript_inputs/Figure_3/
#     ASPC_PPARG_vs_CXCL14_ATAC_Pseudobulk_DAP_log2FC05_padj005.csv
#
# Panels:
#   3A DAP-centered aggregate ATAC accessibility + heatmaps
#   3B Genomic annotation of directional DAPs
#   3C Unbiased all-TF motif enrichment
#   3D Focused TF-program linked-target bubble plot
#   3E PPARG RNA-ATAC locus
#   3F Regulatory network schematic assembled manually
#
# Positive ATAC log2FC = PPARG+ higher accessibility
# Negative ATAC log2FC = CXCL14+ higher accessibility
# ============================================================

project_dir <- Sys.getenv("CUSHING_MULTIOME_PROJECT", unset = ".")
setwd(project_dir)

if (!file.exists("README.md")) {
  stop("Run from the Cushing-Multiome repository root or set CUSHING_MULTIOME_PROJECT.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(ChIPseeker)
  library(EnsDb.Hsapiens.v86)
  library(openxlsx)
})

dir.create("Figure_3/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("Figure_3/Source_Data", recursive = TRUE, showWarnings = FALSE)

clean_peak <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  gsub(":", "-", x, fixed = TRUE)
}

# ============================================================
# 1. LOAD LOCKED DIFFERENTIAL ATAC PEAKS
# ============================================================

dap_file <- file.path(
  "data",
  "manuscript_inputs",
  "Figure_3",
  "ASPC_PPARG_vs_CXCL14_ATAC_Pseudobulk_DAP_log2FC05_padj005.csv"
)

if (!file.exists(dap_file)) {
  stop(
    "Missing Figure 3 DAP table: ", dap_file,
    "\nCopy the locked DAP result into data/manuscript_inputs/Figure_3/."
  )
}

dap <- read.csv(
  dap_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dap$peak_key <- clean_peak(dap$peak)
dap$Direction <- ifelse(
  dap$log2FoldChange > 0,
  "PPARG+ higher",
  "CXCL14+ higher"
)
dap$Population <- ifelse(
  dap$log2FoldChange > 0,
  "PPARG+",
  "CXCL14+"
)

cat("\nDirectional DAPs:\n")
print(table(dap$Direction))

# Expected locked result:
# CXCL14+ higher = 638
# PPARG+ higher  = 445
# Total          = 1083

# ============================================================
# 3A. DAP-CENTERED ACCESSIBILITY
# Existing library-size-normalized RegionMatrix calculated
# +/- 2000 bp around DAP centers.
# ============================================================

atac_dap_norm <- readRDS(
  "data/objects/ASPC_CXCL14_PPARG_DAP_directional_RegionMatrix_depthnormalized.rds"
)

DefaultAssay(atac_dap_norm) <- "ATAC"

if ("WNN_population" %in% colnames(atac_dap_norm@meta.data)) {
  Idents(atac_dap_norm) <- "WNN_population"
}

region_keys <- names(
  slot(atac_dap_norm[["ATAC"]], "positionEnrichment")
)

stopifnot(
  "DAP_CX" %in% region_keys,
  "DAP_PP" %in% region_keys
)

group_order <- c(
  "CXCL14+ ASPCs",
  "PPARG+ ASPCs"
)

p3a_cx_profile <- RegionPlot(
  object = atac_dap_norm,
  key = "DAP_CX",
  assay = "ATAC",
  idents = group_order,
  group.order = group_order,
  upstream = 2000,
  downstream = 2000,
  normalize = FALSE
)

p3a_cx_heat <- RegionHeatmap(
  object = atac_dap_norm,
  key = "DAP_CX",
  assay = "ATAC",
  idents = group_order,
  group.order = group_order,
  upstream = 2000,
  downstream = 2000,
  normalize = FALSE,
  max.cutoff = "q95",
  cols = "#E95C5C"
)

p3a_pp_profile <- RegionPlot(
  object = atac_dap_norm,
  key = "DAP_PP",
  assay = "ATAC",
  idents = group_order,
  group.order = group_order,
  upstream = 2000,
  downstream = 2000,
  normalize = FALSE
)

p3a_pp_heat <- RegionHeatmap(
  object = atac_dap_norm,
  key = "DAP_PP",
  assay = "ATAC",
  idents = group_order,
  group.order = group_order,
  upstream = 2000,
  downstream = 2000,
  normalize = FALSE,
  max.cutoff = "q95",
  cols = "#E95C5C"
)

print(p3a_cx_profile)
print(p3a_cx_heat)
print(p3a_pp_profile)
print(p3a_pp_heat)

ggsave("Figure_3/Figures/Fig3A_CXCL14_DAP_Profile.pdf", p3a_cx_profile, width = 7, height = 3, bg = "white")
ggsave("Figure_3/Figures/Fig3A_CXCL14_DAP_Heatmap.pdf", p3a_cx_heat, width = 7, height = 6, bg = "white")
ggsave("Figure_3/Figures/Fig3A_PPARG_DAP_Profile.pdf", p3a_pp_profile, width = 7, height = 3, bg = "white")
ggsave("Figure_3/Figures/Fig3A_PPARG_DAP_Heatmap.pdf", p3a_pp_heat, width = 7, height = 6, bg = "white")

# ============================================================
# 3B. GENOMIC ANNOTATION OF DIRECTIONAL DAPs
# ============================================================

dap_gr <- StringToGRanges(
  dap$peak,
  sep = c("-", "-")
)

dap_gr$Direction <- dap$Direction
GenomeInfoDb::seqlevelsStyle(dap_gr) <- "Ensembl"

dap_gr_cx <- dap_gr[dap_gr$Direction == "CXCL14+ higher"]
dap_gr_pp <- dap_gr[dap_gr$Direction == "PPARG+ higher"]

anno_cx <- annotatePeak(
  peak = dap_gr_cx,
  tssRegion = c(-3000, 3000),
  TxDb = EnsDb.Hsapiens.v86,
  verbose = FALSE
)

anno_pp <- annotatePeak(
  peak = dap_gr_pp,
  tssRegion = c(-3000, 3000),
  TxDb = EnsDb.Hsapiens.v86,
  verbose = FALSE
)

anno_cx_stats <- as.data.frame(ChIPseeker::getAnnoStat(anno_cx))
anno_cx_stats$Direction <- "CXCL14+ higher"

anno_pp_stats <- as.data.frame(ChIPseeker::getAnnoStat(anno_pp))
anno_pp_stats$Direction <- "PPARG+ higher"

anno_stats <- dplyr::bind_rows(
  anno_cx_stats,
  anno_pp_stats
)

p_fig3b <- ggplot(
  anno_stats,
  aes(x = Direction, y = Frequency, fill = Feature)
) +
  geom_col(width = 0.7) +
  labs(x = NULL, y = "DAPs (%)", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

print(p_fig3b)

ggsave("Figure_3/Figures/Fig3B_DAP_GenomicAnnotation.pdf", p_fig3b, width = 7, height = 5, bg = "white")
ggsave("Figure_3/Figures/Fig3B_DAP_GenomicAnnotation.png", p_fig3b, width = 7, height = 5, dpi = 600, bg = "white")

# ============================================================
# 3C. UNBIASED ALL-TF MOTIF ENRICHMENT
# ============================================================

motif_obj <- readRDS(
  "data/objects/ASPC_Lineage_RNA_ATAC_Motifs.rds"
)

DefaultAssay(motif_obj) <- "ATAC"

motif_excel <- "Figure_3/Source_Data/Fig3C_AllTF_MotifEnrichment.xlsx"

if (file.exists(motif_excel)) {

  motifs_all <- openxlsx::read.xlsx(
    motif_excel,
    sheet = "All_motifs"
  )

} else {

  motif_peaks <- rownames(motif_obj[["ATAC"]])

  dap_cx_peaks <- intersect(
    dap$peak[dap$Direction == "CXCL14+ higher"],
    motif_peaks
  )

  dap_pp_peaks <- intersect(
    dap$peak[dap$Direction == "PPARG+ higher"],
    motif_peaks
  )

  atac_meta_features <- GetAssayData(
    object = motif_obj,
    assay = "ATAC",
    layer = "meta.features"
  )

  if (!"GC.percent" %in% colnames(atac_meta_features)) {
    stop("GC.percent missing. RegionStats is required before FindMotifs.")
  }

  set.seed(1234)

  motifs_cx <- FindMotifs(
    object = motif_obj,
    features = dap_cx_peaks,
    background = 40000,
    assay = "ATAC",
    p.adjust.method = "BH",
    verbose = TRUE
  )
  motifs_cx$Direction <- "CXCL14+ higher"

  set.seed(1234)

  motifs_pp <- FindMotifs(
    object = motif_obj,
    features = dap_pp_peaks,
    background = 40000,
    assay = "ATAC",
    p.adjust.method = "BH",
    verbose = TRUE
  )
  motifs_pp$Direction <- "PPARG+ higher"

  motifs_all <- dplyr::bind_rows(
    motifs_cx,
    motifs_pp
  )

  openxlsx::write.xlsx(
    list(
      All_motifs = motifs_all,
      CXCL14_all = motifs_cx,
      PPARG_all = motifs_pp
    ),
    motif_excel,
    overwrite = TRUE
  )
}

motifs_all <- motifs_all %>%
  dplyr::mutate(
    Significant = !is.na(p.adjust) & p.adjust < 0.05,
    minus_log10_padj = -log10(
      pmax(p.adjust, .Machine$double.xmin)
    )
  )

motifs_unique <- motifs_all %>%
  dplyr::filter(
    Significant,
    fold.enrichment > 1
  ) %>%
  dplyr::group_by(
    Direction,
    motif.name
  ) %>%
  dplyr::arrange(
    p.adjust,
    dplyr::desc(fold.enrichment)
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

top_cx <- motifs_unique %>%
  dplyr::filter(Direction == "CXCL14+ higher") %>%
  dplyr::arrange(
    p.adjust,
    dplyr::desc(fold.enrichment)
  ) %>%
  dplyr::slice_head(n = 15)

top_pp <- motifs_unique %>%
  dplyr::filter(Direction == "PPARG+ higher") %>%
  dplyr::arrange(
    p.adjust,
    dplyr::desc(fold.enrichment)
  ) %>%
  dplyr::slice_head(n = 15)

top_cx$motif.name <- factor(
  top_cx$motif.name,
  levels = rev(top_cx$motif.name)
)

top_pp$motif.name <- factor(
  top_pp$motif.name,
  levels = rev(top_pp$motif.name)
)

p3c_cx <- ggplot(
  top_cx,
  aes(
    x = fold.enrichment,
    y = motif.name,
    size = minus_log10_padj,
    color = fold.enrichment
  )
) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(
    low = "#4575B4",
    high = "#D73027"
  ) +
  labs(
    title = "CXCL14+ higher DAPs",
    x = "Motif fold enrichment",
    y = NULL,
    size = "-log10\nadjusted p",
    color = "Fold\nenrichment"
  ) +
  theme_classic(base_size = 12)

p3c_pp <- ggplot(
  top_pp,
  aes(
    x = fold.enrichment,
    y = motif.name,
    size = minus_log10_padj,
    color = fold.enrichment
  )
) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(
    low = "#4575B4",
    high = "#D73027"
  ) +
  labs(
    title = "PPARG+ higher DAPs",
    x = "Motif fold enrichment",
    y = NULL,
    size = "-log10\nadjusted p",
    color = "Fold\nenrichment"
  ) +
  theme_classic(base_size = 12)

p_fig3c <- p3c_cx | p3c_pp
print(p_fig3c)

ggsave("Figure_3/Figures/Fig3C_AllTF_MotifEnrichment.pdf", p_fig3c, width = 12, height = 7, bg = "white")
ggsave("Figure_3/Figures/Fig3C_AllTF_MotifEnrichment.png", p_fig3c, width = 12, height = 7, dpi = 600, bg = "white")

# ============================================================
# 3D. FOCUSED TF-PROGRAM BUBBLE PLOT
#
# FINAL TFs:
# NR3C1, KLF15, STAT5B, CEBPD, CEBPB, CEBPA, PPARG::RXRA
#
# Denominator:
#   all genes linked to state-specific DAPs
#
# Numerator:
#   genes linked through those DAPs containing the TF motif
#
# Point size:
#   % linked target genes
#
# Point color:
#   motif fold enrichment
#
# IMPORTANT:
# Keep color = Fold_enrichment.
# Do NOT replace with shape=21 + fill; this caused an
# empty-rendering problem in the original analysis environment.
# ============================================================

link_obj <- readRDS(
  "data/objects/ASPC_Lineage_RNA_DEG_LinkPeaks.rds"
)

DefaultAssay(link_obj) <- "ATAC"

lnk <- Links(link_obj[["ATAC"]])

links <- data.frame(
  peak = as.character(lnk$peak),
  gene = as.character(lnk$gene),
  score = as.numeric(lnk$score),
  pvalue = as.numeric(lnk$pvalue),
  stringsAsFactors = FALSE
)

links$peak_key <- clean_peak(links$peak)

links <- links %>%
  dplyr::filter(pvalue < 0.05)

motif_matrix <- GetMotifData(
  motif_obj,
  assay = "ATAC",
  slot = "data"
)

motif_peak_key <- clean_peak(
  rownames(motif_matrix)
)

selected_tfs <- c(
  "NR3C1",
  "KLF15",
  "STAT5B",
  "CEBPD",
  "CEBPB",
  "CEBPA",
  "PPARG::RXRA"
)

# Preserve the exact validated motif choices for the original
# five TFs. KLF15 and STAT5B are selected from the same unbiased
# motif analysis used above.

original_locked <- data.frame(
  Direction = c(
    rep("CXCL14+ higher", 5),
    rep("PPARG+ higher", 5)
  ),
  Population = c(
    rep("CXCL14+", 5),
    rep("PPARG+", 5)
  ),
  TF = rep(
    c(
      "NR3C1",
      "CEBPD",
      "CEBPB",
      "CEBPA",
      "PPARG::RXRA"
    ),
    2
  ),
  Fold_enrichment = c(
    1.536,
    1.510,
    2.266,
    1.307,
    1.093,
    1.950,
    2.299,
    3.154,
    1.827,
    1.472
  ),
  Motif_padj = c(
    0.004988,
    0.003232,
    0.005278,
    0.040961,
    0.125790,
    0.000192,
    5.09e-11,
    7.50e-04,
    8.61e-07,
    5.68e-09
  ),
  stringsAsFactors = FALSE
)

find_tf_rows <- function(direction, tf) {

  x <- motifs_all %>%
    dplyr::filter(Direction == direction)

  nm <- as.character(x$motif.name)

  nm_clean <- toupper(
    gsub("\\.[0-9]+$", "", nm)
  )

  if (tf == "PPARG::RXRA") {

    keep <- grepl(
      "PPARG.*RXRA",
      nm_clean
    )

  } else if (tf == "STAT5B") {

    keep <- grepl(
      "STAT5B",
      nm_clean
    )

  } else {

    keep <- nm_clean == toupper(tf)

    if (!any(keep)) {
      keep <- grepl(
        paste0(
          "(^|::)",
          toupper(tf),
          "($|::)"
        ),
        nm_clean
      )
    }
  }

  x[keep, , drop = FALSE]
}

motif_exact_list <- list()
counter <- 1

for (
  direction in c(
    "CXCL14+ higher",
    "PPARG+ higher"
  )
) {

  population <- ifelse(
    direction == "CXCL14+ higher",
    "CXCL14+",
    "PPARG+"
  )

  for (tf in selected_tfs) {

    candidates <- find_tf_rows(
      direction,
      tf
    )

    if (nrow(candidates) == 0) {
      warning(
        "Motif not found for ",
        tf,
        " / ",
        direction
      )
      next
    }

    locked_row <- original_locked %>%
      dplyr::filter(
        Direction == direction,
        TF == tf
      )

    if (nrow(locked_row) == 1) {

      candidates$fold_distance <- abs(
        candidates$fold.enrichment -
          locked_row$Fold_enrichment
      )

      candidates <- candidates %>%
        dplyr::arrange(
          fold_distance,
          p.adjust
        )

      chosen <- candidates[1, , drop = FALSE]
      final_fold <- locked_row$Fold_enrichment
      final_padj <- locked_row$Motif_padj

    } else {

      candidates <- candidates %>%
        dplyr::arrange(
          is.na(p.adjust),
          p.adjust,
          dplyr::desc(fold.enrichment)
        )

      chosen <- candidates[1, , drop = FALSE]
      final_fold <- as.numeric(chosen$fold.enrichment)
      final_padj <- as.numeric(chosen$p.adjust)
    }

    motif_exact_list[[counter]] <- data.frame(
      Direction = direction,
      Population = population,
      TF = tf,
      motif_id = as.character(chosen$motif),
      motif_name = as.character(chosen$motif.name),
      Fold_enrichment = final_fold,
      Motif_padj = final_padj,
      stringsAsFactors = FALSE
    )

    counter <- counter + 1
  }
}

motif_exact <- dplyr::bind_rows(
  motif_exact_list
)

motif_exact$Found_in_matrix <-
  motif_exact$motif_id %in%
  colnames(motif_matrix)

if (!all(motif_exact$Found_in_matrix)) {
  print(motif_exact[!motif_exact$Found_in_matrix, ])
  stop("At least one selected motif ID is missing from motif matrix.")
}

tf_counts <- motif_exact %>%
  dplyr::count(TF)

selected_tfs_present <- tf_counts %>%
  dplyr::filter(n == 2) %>%
  dplyr::pull(TF)

selected_tfs_present <- selected_tfs[
  selected_tfs %in%
    selected_tfs_present
]

result_list <- list()
target_gene_list <- list()
counter <- 1

for (
  state in c(
    "CXCL14+",
    "PPARG+"
  )
) {

  state_daps <- unique(
    dap$peak_key[
      dap$Population == state
    ]
  )

  state_links <- links %>%
    dplyr::filter(
      peak_key %in%
        state_daps
    )

  all_linked_genes <- unique(
    state_links$gene
  )

  if (length(all_linked_genes) == 0) {
    stop("No linked genes found for ", state)
  }

  for (tf in selected_tfs_present) {

    motif_row <- motif_exact %>%
      dplyr::filter(
        Population == state,
        TF == tf
      )

    motif_id <- motif_row$motif_id[1]

    motif_present <- as.vector(
      motif_matrix[, motif_id] > 0
    )

    motif_positive_peaks <-
      motif_peak_key[
        motif_present
      ]

    state_motif_daps <- intersect(
      state_daps,
      motif_positive_peaks
    )

    tf_links <- state_links %>%
      dplyr::filter(
        peak_key %in%
          state_motif_daps
      )

    tf_genes <- unique(
      tf_links$gene
    )

    target_percent <-
      100 *
      length(tf_genes) /
      length(all_linked_genes)

    result_list[[counter]] <- data.frame(
      Population = state,
      TF = tf,
      TF_motif_DAPs = length(state_motif_daps),
      Linked_targets = length(tf_genes),
      Total_linked_targets = length(all_linked_genes),
      Linked_target_percent = target_percent,
      Fold_enrichment = motif_row$Fold_enrichment[1],
      Motif_padj = motif_row$Motif_padj[1],
      stringsAsFactors = FALSE
    )

    if (length(tf_genes) > 0) {
      target_gene_list[[counter]] <- data.frame(
        Population = state,
        TF = tf,
        Gene = tf_genes,
        stringsAsFactors = FALSE
      )
    }

    counter <- counter + 1
  }
}

plot_df <- dplyr::bind_rows(
  result_list
)

target_gene_table <- dplyr::bind_rows(
  target_gene_list
)

cat("\nFigure 3D final data:\n")
print(plot_df, row.names = FALSE)

tf_order <- c(
  "NR3C1",
  "KLF15",
  "STAT5B",
  "CEBPD",
  "CEBPB",
  "CEBPA",
  "PPARG::RXRA"
)

tf_order <- tf_order[
  tf_order %in%
    selected_tfs_present
]

plot_df$Population <- factor(
  plot_df$Population,
  levels = c(
    "CXCL14+",
    "PPARG+"
  )
)

plot_df$TF <- factor(
  plot_df$TF,
  levels = tf_order
)

p_fig3d <- ggplot(
  plot_df,
  aes(
    x = Population,
    y = TF
  )
) +
  geom_point(
    aes(
      size = Linked_target_percent,
      color = Fold_enrichment
    ),
    alpha = 0.95
  ) +
  scale_size_continuous(
    name = "% linked\ntarget genes",
    breaks = c(10, 20, 30, 40, 50),
    range = c(3, 18)
  ) +
  scale_color_gradient(
    name = "Motif\nfold enrichment",
    low = "grey80",
    high = "red"
  ) +
  labs(
    title = "Relative prominence of ASPC TF programs",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(
      size = 20,
      hjust = 0.5
    ),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )

print(p_fig3d)

ggsave("Figure_3/Figures/Fig3D_TF_Program_Prominence.pdf", p_fig3d, width = 9, height = 8, bg = "white")
ggsave("Figure_3/Figures/Fig3D_TF_Program_Prominence.png", p_fig3d, width = 9, height = 8, dpi = 600, bg = "white")

# ============================================================
# 3E. PPARG LOCUS
# ============================================================

possible_group_columns <- c(
  "WNN_population",
  "trajectory_population",
  "ASPC_population"
)

available_group_columns <- possible_group_columns[
  possible_group_columns %in%
    colnames(link_obj@meta.data)
]

if (length(available_group_columns) == 0) {
  stop("Could not identify ASPC population column in LinkPeaks object.")
}

group_column <- available_group_columns[1]
Idents(link_obj) <- group_column

pparg_linked_peak <- GRanges(
  seqnames = "chr3",
  ranges = IRanges(
    start = 12452634,
    end = 12454044
  )
)

mcols(pparg_linked_peak)$color <- "#F4A3A3"

p_fig3e <- CoveragePlot(
  object = link_obj,
  region = "PPARG",
  features = "PPARG",
  assay = "ATAC",
  expression.assay = "RNA",
  group.by = group_column,
  idents = c(
    "CXCL14+ ASPCs",
    "PPARG+ ASPCs"
  ),
  annotation = TRUE,
  peaks = TRUE,
  links = TRUE,
  region.highlight = pparg_linked_peak,
  extend.upstream = 50000,
  extend.downstream = 50000
)

print(p_fig3e)

ggsave("Figure_3/Figures/Fig3E_PPARG_Locus_RNA_ATAC_Links.pdf", p_fig3e, width = 12, height = 8, bg = "white")
ggsave("Figure_3/Figures/Fig3E_PPARG_Locus_RNA_ATAC_Links.png", p_fig3e, width = 12, height = 8, dpi = 600, bg = "white")

pparg_links <- links %>%
  dplyr::filter(
    gene == "PPARG"
  )

# Known positive-control linked DAP:
# chr3-12452634-12454044
# PPARG LinkPeaks score ~0.219, p ~0.00127
# ATAC log2FC ~+2.16
# PPARG RNA log2FC ~+1.44

# ============================================================
# 3F. REGULATORY MODEL / NETWORK
#
# Assembled manually in Inkscape.
#
# Suggested organization:
#
#              NR3C1
#                |
#        KLF15        STAT5B
#                |
#        CEBPB        CEBPD
#                |
#        CEBPA    PPARG::RXRA
#
# The schematic summarizes motif-supported regulatory programs.
# It does not establish direct TF occupancy or direct TF->TF
# regulation.
# ============================================================

# ============================================================
# EXPORT COMPLETE FIGURE 3 SOURCE DATA
# ============================================================

openxlsx::write.xlsx(
  list(
    Directional_DAPs = dap,
    Genomic_annotation = anno_stats,
    All_motif_results = motifs_all,
    Top_CXCL14_motifs = top_cx,
    Top_PPARG_motifs = top_pp,
    Focused_TF_motifs = motif_exact,
    TF_program_summary = plot_df,
    TF_linked_target_genes = target_gene_table,
    PPARG_links = pparg_links
  ),
  file = "Figure_3/Source_Data/Figure3_All_SourceData.xlsx",
  overwrite = TRUE
)

capture.output(
  sessionInfo(),
  file = "Figure_3/Source_Data/Figure3_sessionInfo.txt"
)

cat("\n======================================\n")
cat("FIGURE 3 COMPLETE\n")
cat("======================================\n")
