# Manuscript analysis inputs

The final manuscript scripts in `scripts/manuscript/` intentionally use repository-relative paths. Large checkpoint objects and locked result tables are kept local and are ignored by Git.

Create the following local directories from the repository root:

```text
data/
├── objects/
├── manuscript_inputs/
│   ├── Figure_2/
│   └── Figure_3/
```

## Figure 1

Required objects:

```text
data/objects/AllCells_RNA_final.rds
data/objects/AllCells_WNN_final.rds
```

The script creates:

```text
data/objects/AllCells_RNA_manuscript.rds
data/objects/AllCells_WNN_manuscript.rds
```

## Figure 2

Required objects:

```text
data/objects/ASPC_WNN40_res04_final.rds
data/objects/Adipogenesis_Trajectory_final.rds
```

Required locked RNA pseudobulk results:

```text
data/manuscript_inputs/Figure_2/
├── ASPC_PPARG_vs_CXCL14_Pseudobulk_ALL.csv
└── ASPC_PPARG_vs_CXCL14_Pseudobulk_DEGs_log2FC05_padj005.csv
```

These tables correspond to the patient-level DESeq2 model:

```text
~ patient + population
```

with PPARG+ ASPCs compared with CXCL14+ ASPCs.

## Figure 3

Required objects:

```text
data/objects/ASPC_CXCL14_PPARG_DAP_directional_RegionMatrix_depthnormalized.rds
data/objects/ASPC_Lineage_RNA_ATAC_Motifs.rds
data/objects/ASPC_Lineage_RNA_DEG_LinkPeaks.rds
```

Required locked differential ATAC table:

```text
data/manuscript_inputs/Figure_3/
└── ASPC_PPARG_vs_CXCL14_ATAC_Pseudobulk_DAP_log2FC05_padj005.csv
```

The final DAP set contains 1,083 peaks:

- 638 with greater accessibility in CXCL14+ ASPCs
- 445 with greater accessibility in PPARG+ ASPCs

## Running the scripts

From the repository root:

```r
source("scripts/manuscript/01_Figure1_Global_RNA_Multiome.R")
source("scripts/manuscript/02_Figure2_ASPC_Adipogenesis.R")
source("scripts/manuscript/03_Figure3_ATAC_RegulatoryNetwork.R")
```

Alternatively, set the environment variable `CUSHING_MULTIOME_PROJECT` to the repository root before running a script from another working directory.

The manuscript scripts reuse saved dimensional reductions, cluster assignments and trajectory checkpoints where explicitly stated. They should not be modified to silently recompute PCA, Harmony, WNN clustering, UMAP or pseudotime unless a new analysis version is intentionally created.
