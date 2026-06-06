# GTEx Expression

Scripts and selected data files for evaluating whether ASD-PPI prey genes are more highly expressed than expected in GTEx brain tissues.

## Scripts

- `scripts/01_downloadGtexData.R`: Downloads and formats GTEx median TPM data.
- `scripts/02_GTEx_asdPPI100_analysis.R`: Runs expression and permutation analyses for ASD-PPI prey genes in GTEx.

## Data

- `data/gtex_medianTPM.RData`: Formatted GTEx median TPM object.
- `data/rawdata/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`: GTEx v8 median gene-level TPM file.

## Outputs

- `output/asdPPI_gtexBrain_expression_boxplot.pdf`: GTEx brain expression comparison plot.
- `output/gtexBrain_expression_asdPPIprey_vs_permuted.pdf`: Permutation comparison plot for ASD-PPI prey expression in GTEx brain samples.
