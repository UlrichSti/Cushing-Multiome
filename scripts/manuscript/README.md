# Manuscript analysis scripts

These scripts reproduce the finalized computational panels for Figures 1–3 of the Cushing ASPC manuscript.

## Canonical scripts

- `01_Figure1_Global_RNA_Multiome.R`  
  Global adipose cell landscape, broad marker expression, patient-level cell-type proportions, RNA/ATAC QC summaries and WNN modality weights.

- `02_Figure2_ASPC_Adipogenesis.R`  
  ASPC states, patient-level ASPC composition, saved adipogenic pseudotime, patient-first pseudotime distributions, CXCL14+ versus PPARG+ RNA pseudobulk differential expression, GO enrichment and transcription-factor dynamics along pseudotime.

- `03_Figure3_ATAC_RegulatoryNetwork.R`  
  Directional differential accessibility, genomic annotation, unbiased motif enrichment, focused motif-supported TF target programs and the PPARG regulatory locus.

## Reproducibility rules

The manuscript scripts start from finalized checkpoint objects. In particular, they do not silently rerun PCA, Harmony, WNN clustering, UMAP or saved pseudotime when the accepted manuscript figure depends on a frozen embedding or trajectory.

Patients are treated as biological replicates in group-level abundance, pseudotime and pseudobulk analyses.

Large RDS checkpoint objects and locked intermediate result tables are not committed to GitHub. Their expected repository-relative locations are documented in `docs/Manuscript_inputs.md`.

Generated `Figure_1/`, `Figure_2/` and `Figure_3/` directories are ignored by Git.
