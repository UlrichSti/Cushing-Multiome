# Manuscript analysis scripts

This directory is reserved for the finalized manuscript-facing analysis scripts.

Planned canonical scripts:

- `01_Figure1_Global_RNA_Multiome.R`
- `02_Figure2_ASPC_Adipogenesis.R`
- `03_Figure3_ATAC_RegulatoryNetwork.R`

These scripts should start from finalized checkpoint objects rather than rerunning upstream PCA, Harmony, WNN clustering or UMAP unless explicitly required by the analysis.

Generated figures and large RDS checkpoint objects should not be committed to the repository. Source-data tables that are sufficiently small may be committed when appropriate.
