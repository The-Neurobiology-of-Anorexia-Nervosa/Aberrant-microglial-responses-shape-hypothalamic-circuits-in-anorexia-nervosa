# Aberrant microglial responses shape hypothalamic circuits in anorexia nervosa

This repository contains analysis code supporting the manuscript
" aberrant microglial responses shape hypothalamic circuits in anorexia nervosa."

## Repository structure

### `01_rnaseq`

RNA-seq quality control, differential-expression analysis, GSEA and
figure generation.

### `02_spatial_transcriptomics`

Spatial transcriptomic quality control, differential-expression analysis,
GSEA and figure generation.

### `03_microglia_states`

Fisher's exact tests and FGSEA comparing study-derived results with
published human and mouse microglial-state markers.

## Data availability

Large sequencing files are not stored in GitHub. Generated and public datasets,
their accession numbers and their use in the manuscript are documented in
`metadata/dataset_manifest.csv`.

Computer-specific data paths are defined using a local configuration file
that is not committed to GitHub.

## Reproducibility

The analysis scripts use project-relative paths through the R package `here`.
Package versions will be recorded using `renv`.

Detailed instructions are available in the README file of each analysis module.

## License

The code in this repository is available under the MIT License.