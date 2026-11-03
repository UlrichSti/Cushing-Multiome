# ============================================================
# FIGURE 2 - ASPC STATES AND ADIPOGENIC PROGRESSION
# Reproducible manuscript script
#
# Run from repository root, or set:
#   CUSHING_MULTIOME_PROJECT=/path/to/Cushing-Multiome
#
# Required local checkpoint objects (not committed):
#   data/objects/ASPC_WNN40_res04_final.rds
#   data/objects/Adipogenesis_Trajectory_final.rds
#
# Required local result tables:
#   data/manuscript_inputs/Figure_2/ASPC_PPARG_vs_CXCL14_Pseudobulk_ALL.csv
#   data/manuscript_inputs/Figure_2/ASPC_PPARG_vs_CXCL14_Pseudobulk_DEGs_log2FC05_padj005.csv
#
# Panels:
#   2A ASPC WNN UMAP
#   2B ASPC marker dotplot
#   2C Patient-level ASPC proportions and CXCL14/PPARG balance
#   2D Saved adipogenic pseudotime on exact ASPC WNN UMAP
#   2E Patient-first pseudotime density and population means
#   2F CXCL14+ vs PPARG+ RNA pseudobulk volcano
#   2G GO Biological Process enrichment
#   2H TF expression dynamics across pseudotime
#
# Locked analysis principles:
#   - Patient = biological replicate.
#   - No recalculation of ASPC WNN UMAP or saved pseudotime.
#   - RNA pseudobulk: DESeq2 design ~ patient + population.
#   - Positive RNA log2FC = PPARG+ higher.
#   - DEG threshold: padj < 0.05 & |log2FC| > 0.5.
#   - TF pseudotime: patient x bin pseudobulk, 12 equal-cell bins,
#     minimum 20 nuclei/patient/bin, DESeq2 LRT.
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
  library(ggrepel)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(splines)
  library(pheatmap)
  library(openxlsx)
})

dir.create("Figure_2/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("Figure_2/Source_Data", recursive = TRUE, showWarnings = FALSE)
dir.create("Figure_2/Exploratory", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD LOCKED ASPC OBJECT
# ============================================================

aspc <- readRDS("data/objects/ASPC_WNN40_res04_final.rds")

population_col <- c(
  "WNN_population",
  "trajectory_population",
  "ASPC_population"
)

population_col <- population_col[
  population_col %in% colnames(aspc@meta.data)
][1]

if (is.na(population_col)) {
  stop("Could not identify ASPC population metadata column.")
}

patient_col <- c("patient", "Patient", "orig.ident")
patient_col <- patient_col[patient_col %in% colnames(aspc@meta.data)][1]

condition_col <- c("condition", "Condition")
condition_col <- condition_col[condition_col %in% colnames(aspc@meta.data)][1]

if (is.na(patient_col) || is.na(condition_col)) {
  stop("Patient or condition metadata not found.")
}

aspc$Manuscript_population <- as.character(aspc@meta.data[[population_col]])
aspc$Patient <- as.character(aspc@meta.data[[patient_col]])
aspc$Condition <- as.character(aspc@meta.data[[condition_col]])

aspc_order <- c(
  "DPP4+ ASPCs",
  "CXCL14+ ASPCs",
  "PPARG+ ASPCs"
)

aspc$Manuscript_population <- factor(
  aspc$Manuscript_population,
  levels = aspc_order
)

aspc_cols <- c(
  "DPP4+ ASPCs" = "#FFD29D",
  "CXCL14+ ASPCs" = "#E98324",
  "PPARG+ ASPCs" = "#8C2F00"
)

cat("\nASPC population counts:\n")
print(table(aspc$Manuscript_population))

# Expected locked counts:
# DPP4+   881
# CXCL14+ 2724
# PPARG+  1134
# total   4739

# ============================================================
# 2A. ASPC WNN UMAP
# ============================================================

p_fig2a <- DimPlot(
  aspc,
  reduction = "wnn.umap",
  group.by = "Manuscript_population",
  cols = aspc_cols,
  raster = FALSE,
  label = FALSE,
  pt.size = 0.35
) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

print(p_fig2a)

ggsave("Figure_2/Figures/Fig2A_ASPC_WNN_UMAP.pdf", p_fig2a, width = 7, height = 6)
ggsave("Figure_2/Figures/Fig2A_ASPC_WNN_UMAP.png", p_fig2a, width = 7, height = 6, dpi = 600)

wnn_emb <- Embeddings(aspc, "wnn.umap")
fig2a_source <- data.frame(
  Cell = rownames(wnn_emb),
  WNN_UMAP_1 = wnn_emb[, 1],
  WNN_UMAP_2 = wnn_emb[, 2],
  Patient = aspc$Patient[match(rownames(wnn_emb), Cells(aspc))],
  Condition = aspc$Condition[match(rownames(wnn_emb), Cells(aspc))],
  Population = as.character(aspc$Manuscript_population[
    match(rownames(wnn_emb), Cells(aspc))
  ]),
  row.names = NULL
)

# ============================================================
# 2B. MARKER DOTPLOT
# ============================================================

DefaultAssay(aspc) <- "RNA"

aspc_markers <- c(
  "DPP4", "PI16",
  "CXCL14", "COL15A1",
  "PPARG", "FABP4"
)

p_fig2b <- DotPlot(
  aspc,
  features = aspc_markers,
  group.by = "Manuscript_population",
  cols = c("grey90", "black")
) +
  RotatedAxis() +
  labs(
    x = NULL,
    y = NULL,
    colour = "Average expression",
    size = "Percent expressed"
  ) +
  theme_classic()

print(p_fig2b)

ggsave("Figure_2/Figures/Fig2B_ASPC_Marker_DotPlot.pdf", p_fig2b, width = 6.5, height = 4.5)
ggsave("Figure_2/Figures/Fig2B_ASPC_Marker_DotPlot.png", p_fig2b, width = 6.5, height = 4.5, dpi = 600)

# ============================================================
# 2C. PATIENT-LEVEL ASPC COMPOSITION
# ============================================================

patient_info <- aspc@meta.data %>%
  transmute(
    Patient = as.character(.data[[patient_col]]),
    Condition = as.character(.data[[condition_col]])
  ) %>%
  distinct()

aspc_counts <- aspc@meta.data %>%
  transmute(
    Patient = as.character(.data[[patient_col]]),
    Condition = as.character(.data[[condition_col]]),
    Population = as.character(Manuscript_population)
  ) %>%
  count(Patient, Condition, Population, name = "N_cells")

aspc_props <- tidyr::expand_grid(
  Patient = patient_info$Patient,
  Population = aspc_order
) %>%
  left_join(patient_info, by = "Patient") %>%
  left_join(aspc_counts, by = c("Patient", "Condition", "Population")) %>%
  mutate(N_cells = replace_na(N_cells, 0L)) %>%
  group_by(Patient, Condition) %>%
  mutate(
    Total_ASPC = sum(N_cells),
    Proportion = N_cells / Total_ASPC,
    Percent = 100 * Proportion
  ) %>%
  ungroup()

aspc_prop_stats <- aspc_props %>%
  group_by(Population) %>%
  summarise(
    p_value = tryCatch(
      t.test(Percent ~ Condition, data = cur_data(), var.equal = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

aspc_wide <- aspc_props %>%
  select(Patient, Condition, Population, Proportion) %>%
  pivot_wider(names_from = Population, values_from = Proportion)

aspc_balance <- aspc_wide %>%
  mutate(
    CXCL14_PPARG_ratio =
      `CXCL14+ ASPCs` / `PPARG+ ASPCs`,
    log2_CXCL14_PPARG_ratio =
      log2(CXCL14_PPARG_ratio)
  )

balance_p <- tryCatch(
  t.test(log2_CXCL14_PPARG_ratio ~ Condition, data = aspc_balance, var.equal = FALSE)$p.value,
  error = function(e) NA_real_
)

cat("\nASPC proportion Welch tests:\n")
print(aspc_prop_stats)
cat("\nCXCL14/PPARG balance Welch test p =", balance_p, "\n")

# Locked values from the accepted analysis include:
# CXCL14+ proportion p = 0.04215096
# PPARG+ proportion  p = 0.09953971

p_fig2c <- ggplot(
  aspc_props,
  aes(x = Condition, y = Percent, fill = Condition)
) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.35) +
  geom_jitter(width = 0.08, size = 2.5) +
  facet_wrap(~Population, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = NULL, y = "% of ASPCs")

print(p_fig2c)

ggsave("Figure_2/Figures/Fig2C_ASPC_Patient_Proportions.pdf", p_fig2c, width = 9, height = 4.5)
ggsave("Figure_2/Figures/Fig2C_ASPC_Patient_Proportions.png", p_fig2c, width = 9, height = 4.5, dpi = 600)

# ============================================================
# 3. LOAD SAVED ADIPOGENIC TRAJECTORY
# ============================================================

traj <- readRDS("data/objects/Adipogenesis_Trajectory_final.rds")

if (!"pseudotime" %in% colnames(traj@meta.data)) {
  stop("Saved trajectory object does not contain 'pseudotime'.")
}

# ============================================================
# 2D. PROJECT SAVED PSEUDOTIME ON EXACT FIG2A ASPC WNN UMAP
# ============================================================

aspc$pseudotime <- traj$pseudotime[
  match(Cells(aspc), Cells(traj))
]

cat("\nTransferred pseudotime values:", sum(!is.na(aspc$pseudotime)), "/", ncol(aspc), "\n")

if (sum(!is.na(aspc$pseudotime)) != ncol(aspc)) {
  warning("Not all ASPC barcodes received pseudotime.")
}

pt_df <- data.frame(
  Cell = Cells(aspc),
  WNN_UMAP_1 = Embeddings(aspc, "wnn.umap")[Cells(aspc), 1],
  WNN_UMAP_2 = Embeddings(aspc, "wnn.umap")[Cells(aspc), 2],
  pseudotime = aspc$pseudotime,
  Patient = aspc$Patient,
  Condition = aspc$Condition,
  Population = as.character(aspc$Manuscript_population)
)

p_fig2d <- ggplot(
  pt_df,
  aes(x = WNN_UMAP_1, y = WNN_UMAP_2, color = pseudotime)
) +
  geom_point(size = 0.35) +
  scale_color_viridis_c(option = "plasma") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(color = "Pseudotime")

print(p_fig2d)

ggsave("Figure_2/Figures/Fig2D_ASPC_Pseudotime_WNN_UMAP.pdf", p_fig2d, width = 7, height = 6)
ggsave("Figure_2/Figures/Fig2D_ASPC_Pseudotime_WNN_UMAP.png", p_fig2d, width = 7, height = 6, dpi = 600)

# ============================================================
# 2E. PATIENT-FIRST PSEUDOTIME DENSITY
# ============================================================

condition_cols <- c("KT" = "#A0A0A4FF", "CS" = "#0000FFFF")
condition_fills <- c("KT" = "#8F8F95FF", "CS" = "#90BFF9FF")

density_input <- pt_df %>%
  filter(!is.na(pseudotime))

grid_x <- seq(
  min(density_input$pseudotime),
  max(density_input$pseudotime),
  length.out = 256
)

patient_density_list <- lapply(
  split(density_input, density_input$Patient),
  function(x) {
    den <- density(
      x$pseudotime,
      from = min(grid_x),
      to = max(grid_x),
      n = length(grid_x)
    )
    data.frame(
      Patient = x$Patient[1],
      Condition = x$Condition[1],
      pseudotime = den$x,
      density = den$y
    )
  }
)

patient_density <- bind_rows(patient_density_list)

condition_density <- patient_density %>%
  group_by(Condition, pseudotime) %>%
  summarise(
    Mean_density = mean(density),
    SEM_density = sd(density) / sqrt(n()),
    .groups = "drop"
  )

patient_pt_mean <- density_input %>%
  group_by(Patient, Condition) %>%
  summarise(Mean_pseudotime = mean(pseudotime), .groups = "drop")

condition_pt_mean <- patient_pt_mean %>%
  group_by(Condition) %>%
  summarise(
    Mean_pseudotime = mean(Mean_pseudotime),
    SEM_pseudotime = sd(Mean_pseudotime) / sqrt(n()),
    .groups = "drop"
  )

p_fig2e <- ggplot(
  condition_density,
  aes(x = pseudotime, y = Mean_density, color = Condition, fill = Condition)
) +
  geom_ribbon(
    aes(
      ymin = pmax(0, Mean_density - SEM_density),
      ymax = Mean_density + SEM_density
    ),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(linewidth = 1.1) +
  geom_vline(
    data = condition_pt_mean,
    aes(xintercept = Mean_pseudotime, color = Condition),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_rect(
    data = condition_pt_mean,
    aes(
      xmin = Mean_pseudotime - SEM_pseudotime,
      xmax = Mean_pseudotime + SEM_pseudotime,
      ymin = -Inf,
      ymax = Inf,
      fill = Condition
    ),
    inherit.aes = FALSE,
    alpha = 0.08
  ) +
  scale_color_manual(values = condition_cols) +
  scale_fill_manual(values = condition_fills) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  labs(
    x = "Adipogenic pseudotime",
    y = "Density",
    color = NULL,
    fill = NULL
  )

print(p_fig2e)

ggsave("Figure_2/Figures/Fig2E_ASPC_Pseudotime_Density.pdf", p_fig2e, width = 7, height = 5)
ggsave("Figure_2/Figures/Fig2E_ASPC_Pseudotime_Density.png", p_fig2e, width = 7, height = 5, dpi = 600)

# ============================================================
# 2F. RNA PSEUDOBULK VOLCANO
# Locked DESeq2 result:
#   design ~ patient + population
#   positive log2FC = PPARG+ higher
# ============================================================

rna_all_file <- "data/manuscript_inputs/Figure_2/ASPC_PPARG_vs_CXCL14_Pseudobulk_ALL.csv"
rna_deg_file <- "data/manuscript_inputs/Figure_2/ASPC_PPARG_vs_CXCL14_Pseudobulk_DEGs_log2FC05_padj005.csv"

if (!file.exists(rna_all_file) || !file.exists(rna_deg_file)) {
  stop(
    "Figure 2 RNA pseudobulk result CSVs are missing. ",
    "Copy the locked results into data/manuscript_inputs/Figure_2/."
  )
}

rna_all <- read.csv(rna_all_file, stringsAsFactors = FALSE)
rna_deg <- read.csv(rna_deg_file, stringsAsFactors = FALSE)

gene_col <- if ("gene" %in% colnames(rna_all)) "gene" else if ("Gene" %in% colnames(rna_all)) "Gene" else NA_character_
if (is.na(gene_col)) stop("Gene column not found in RNA pseudobulk result.")

rna_all$Gene <- as.character(rna_all[[gene_col]])

rna_all$Significance <- "NS"
rna_all$Significance[
  !is.na(rna_all$padj) &
  rna_all$padj < 0.05 &
  rna_all$log2FoldChange > 0.5
] <- "PPARG+ higher"

rna_all$Significance[
  !is.na(rna_all$padj) &
  rna_all$padj < 0.05 &
  rna_all$log2FoldChange < -0.5
] <- "CXCL14+ higher"

rna_all$minus_log10_padj <- -log10(pmax(rna_all$padj, 1e-300))

label_genes <- c(
  "PPARG", "CEBPA", "CEBPB", "CEBPD",
  "KLF9", "NR3C1", "FKBP5", "TSC22D3"
)

label_df <- rna_all %>%
  filter(Gene %in% label_genes)

p_fig2f <- ggplot(
  rna_all,
  aes(x = log2FoldChange, y = minus_log10_padj, color = Significance)
) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(
    values = c(
      "CXCL14+ higher" = "#2166AC",
      "NS" = "grey75",
      "PPARG+ higher" = "#B2182B"
    )
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = Gene),
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  theme_classic() +
  labs(
    x = "log2 fold change (PPARG+ / CXCL14+)",
    y = "-log10 adjusted p-value",
    color = NULL
  )

print(p_fig2f)

ggsave("Figure_2/Figures/Fig2F_CXCL14_vs_PPARG_RNA_Volcano.pdf", p_fig2f, width = 7, height = 6)
ggsave("Figure_2/Figures/Fig2F_CXCL14_vs_PPARG_RNA_Volcano.png", p_fig2f, width = 7, height = 6, dpi = 600)

# ============================================================
# 2G. GO BIOLOGICAL PROCESS ENRICHMENT
# ============================================================

universe_symbols <- unique(rna_all$Gene[!is.na(rna_all$padj)])

cx_symbols <- unique(
  rna_all$Gene[
    !is.na(rna_all$padj) &
    rna_all$padj < 0.05 &
    rna_all$log2FoldChange < -0.5
  ]
)

pp_symbols <- unique(
  rna_all$Gene[
    !is.na(rna_all$padj) &
    rna_all$padj < 0.05 &
    rna_all$log2FoldChange > 0.5
  ]
)

map_symbols <- function(x) {
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = unique(x),
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
}

universe_entrez <- unique(na.omit(unname(map_symbols(universe_symbols))))
cx_entrez <- unique(na.omit(unname(map_symbols(cx_symbols))))
pp_entrez <- unique(na.omit(unname(map_symbols(pp_symbols))))

ego_cx <- enrichGO(
  gene = cx_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

ego_pp <- enrichGO(
  gene = pp_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

go_cx <- as.data.frame(ego_cx) %>%
  arrange(p.adjust) %>%
  slice_head(n = 8) %>%
  mutate(Direction = "CXCL14+ higher")

go_pp <- as.data.frame(ego_pp) %>%
  arrange(p.adjust) %>%
  slice_head(n = 8) %>%
  mutate(Direction = "PPARG+ higher")

go_plot <- bind_rows(go_cx, go_pp) %>%
  mutate(
    minus_log10_padj = -log10(pmax(p.adjust, 1e-300)),
    Label = factor(
      Description,
      levels = rev(unique(Description[order(Direction, p.adjust)]))
    )
  )

p_fig2g <- ggplot(
  go_plot,
  aes(x = Count, y = Label, size = Count, color = minus_log10_padj)
) +
  geom_point() +
  facet_wrap(~Direction, scales = "free_y") +
  scale_color_gradient(low = "#4575B4", high = "#D73027") +
  theme_classic() +
  labs(
    x = "Gene count",
    y = NULL,
    size = "Gene count",
    color = "-log10 adj. p"
  )

print(p_fig2g)

ggsave("Figure_2/Figures/Fig2G_GO_BP_Enrichment.pdf", p_fig2g, width = 11, height = 6)
ggsave("Figure_2/Figures/Fig2G_GO_BP_Enrichment.png", p_fig2g, width = 11, height = 6, dpi = 600)

# ============================================================
# 2H. TRANSCRIPTION FACTORS ALONG ASPC PSEUDOTIME
# 12 equal-cell global pseudotime bins.
# Patient x bin raw-count pseudobulk.
# Require >=20 nuclei per patient/bin.
# LRT: full ~ patient + ns(Mean_pseudotime, df=3)
#      reduced ~ patient
# ============================================================

aspc_pt <- aspc
DefaultAssay(aspc_pt) <- "RNA"

valid_cells <- Cells(aspc_pt)[!is.na(aspc_pt$pseudotime)]
aspc_pt <- subset(aspc_pt, cells = valid_cells)

pt_rank <- rank(aspc_pt$pseudotime, ties.method = "first")
aspc_pt$Pseudotime_bin <- cut(
  pt_rank,
  breaks = quantile(pt_rank, probs = seq(0, 1, length.out = 13)),
  include.lowest = TRUE,
  labels = paste0("Bin", sprintf("%02d", 1:12))
)

aspc_pt$Patient_bin <- paste(
  aspc_pt$Patient,
  aspc_pt$Pseudotime_bin,
  sep = "__"
)

bin_meta <- aspc_pt@meta.data %>%
  mutate(Cell = rownames(.)) %>%
  group_by(Patient, Pseudotime_bin, Patient_bin) %>%
  summarise(
    N_nuclei = n(),
    Mean_pseudotime = mean(pseudotime),
    .groups = "drop"
  ) %>%
  filter(N_nuclei >= 20)

keep_groups <- bin_meta$Patient_bin
keep_cells <- Cells(aspc_pt)[aspc_pt$Patient_bin %in% keep_groups]
aspc_pt <- subset(aspc_pt, cells = keep_cells)

# Join Seurat v5 RNA layers before retrieving raw counts.
if (inherits(aspc_pt[["RNA"]], "Assay5")) {
  aspc_pt[["RNA"]] <- JoinLayers(aspc_pt[["RNA"]])
}

rna_counts <- GetAssayData(aspc_pt, assay = "RNA", layer = "counts")

pb_groups <- factor(aspc_pt$Patient_bin)
group_levels <- levels(pb_groups)

pb_counts <- vapply(
  group_levels,
  function(g) {
    Matrix::rowSums(rna_counts[, pb_groups == g, drop = FALSE])
  },
  FUN.VALUE = numeric(nrow(rna_counts))
)

rownames(pb_counts) <- rownames(rna_counts)
colnames(pb_counts) <- group_levels

pb_meta <- bin_meta %>%
  filter(Patient_bin %in% group_levels) %>%
  arrange(match(Patient_bin, group_levels))

rownames(pb_meta) <- pb_meta$Patient_bin

stopifnot(identical(colnames(pb_counts), rownames(pb_meta)))

keep_genes <- Matrix::rowSums(pb_counts) >= 10
pb_counts <- pb_counts[keep_genes, , drop = FALSE]

dds_tf <- DESeqDataSetFromMatrix(
  countData = round(pb_counts),
  colData = pb_meta,
  design = ~ Patient + ns(Mean_pseudotime, df = 3)
)

dds_tf <- DESeq(
  dds_tf,
  test = "LRT",
  reduced = ~ Patient,
  quiet = TRUE
)

tf_pt <- as.data.frame(results(dds_tf))
tf_pt$gene <- rownames(tf_pt)

tf_symbols <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "GOALL"),
  keytype = "GOALL",
  columns = c("SYMBOL", "ONTOLOGYALL")
) %>%
  filter(
    GOALL == "GO:0003700",
    ONTOLOGYALL == "MF",
    !is.na(SYMBOL)
  ) %>%
  pull(SYMBOL) %>%
  unique()

tf_pt <- tf_pt %>%
  filter(gene %in% tf_symbols) %>%
  arrange(padj)

cat("\nTF pseudotime LRT results:\n")
print(tibble::as_tibble(tf_pt), n = Inf)

tf_sig <- tf_pt %>%
  filter(!is.na(padj), padj < 0.05)

if (nrow(tf_sig) == 0) {
  stop("No significant pseudotime-associated TFs at padj < 0.05.")
}

vsd <- vst(dds_tf, blind = FALSE)
vst_mat <- assay(vsd)

# Descriptive patient-balanced bin means.
bin_patient_mat <- vst_mat[tf_sig$gene, , drop = FALSE]
bin_levels <- levels(aspc_pt$Pseudotime_bin)

bin_mean_list <- lapply(
  bin_levels,
  function(b) {
    sample_ids <- pb_meta$Patient_bin[pb_meta$Pseudotime_bin == b]
    if (length(sample_ids) == 0) {
      rep(NA_real_, nrow(bin_patient_mat))
    } else {
      rowMeans(bin_patient_mat[, sample_ids, drop = FALSE])
    }
  }
)

bin_mean <- do.call(cbind, bin_mean_list)
rownames(bin_mean) <- rownames(bin_patient_mat)
colnames(bin_mean) <- bin_levels

# Descriptive Spearman rho using mean bin values.
bin_mid_pt <- pb_meta %>%
  group_by(Pseudotime_bin) %>%
  summarise(Mean_pseudotime = mean(Mean_pseudotime), .groups = "drop")

rho <- apply(
  bin_mean,
  1,
  function(x) {
    ok <- !is.na(x)
    if (sum(ok) < 3) return(NA_real_)
    cor(
      x[ok],
      bin_mid_pt$Mean_pseudotime[
        match(colnames(bin_mean)[ok], bin_mid_pt$Pseudotime_bin)
      ],
      method = "spearman"
    )
  }
)

tf_sig$Spearman_rho <- rho[tf_sig$gene]

# Top up to 40 significant TFs, ordered by peak-expression bin.
top_tf <- tf_sig %>%
  arrange(padj) %>%
  slice_head(n = min(40, n())) %>%
  pull(gene)

heat_mat <- bin_mean[top_tf, , drop = FALSE]
heat_z <- t(scale(t(heat_mat)))

peak_bin <- apply(
  heat_z,
  1,
  function(x) {
    if (all(is.na(x))) return(Inf)
    which.max(x)
  }
)

heat_z <- heat_z[
  order(peak_bin, decreasing = FALSE),
  ,
  drop = FALSE
]

pdf(
  "Figure_2/Figures/Fig2H_TF_Pseudotime_Heatmap.pdf",
  width = 8,
  height = 10
)

pheatmap(
  heat_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  border_color = NA,
  main = "TF expression along adipogenic pseudotime"
)

dev.off()

png(
  "Figure_2/Figures/Fig2H_TF_Pseudotime_Heatmap.png",
  width = 2400,
  height = 3000,
  res = 300
)

pheatmap(
  heat_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  border_color = NA,
  main = "TF expression along adipogenic pseudotime"
)

dev.off()

# ============================================================
# EXPORT FIGURE 2 SOURCE DATA
# ============================================================

openxlsx::write.xlsx(
  list(
    Fig2A_WNN_UMAP = fig2a_source,
    Fig2C_ASPC_proportions = aspc_props,
    Fig2C_Proportion_stats = aspc_prop_stats,
    Fig2C_CXCL14_PPARG_balance = aspc_balance,
    Fig2D_Pseudotime_UMAP = pt_df,
    Fig2E_Patient_density = patient_density,
    Fig2E_Condition_density = condition_density,
    Fig2E_Patient_mean_pseudotime = patient_pt_mean,
    Fig2E_Condition_mean_pseudotime = condition_pt_mean,
    Fig2F_RNA_all = rna_all,
    Fig2F_RNA_DEGs = rna_deg,
    Fig2G_CXCL14_GO = go_cx,
    Fig2G_PPARG_GO = go_pp,
    Fig2H_TF_LRT = tf_pt,
    Fig2H_TF_significant = tf_sig,
    Fig2H_Pseudobulk_metadata = pb_meta
  ),
  "Figure_2/Source_Data/Figure2_All_SourceData.xlsx",
  overwrite = TRUE
)

capture.output(
  sessionInfo(),
  file = "Figure_2/Source_Data/Figure2_sessionInfo.txt"
)

cat("\n======================================\n")
cat("FIGURE 2 COMPLETE\n")
cat("======================================\n")
