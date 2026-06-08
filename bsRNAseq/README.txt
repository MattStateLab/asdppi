# BrainSpan RNA-seq

Scripts and selected data files for evaluating ASD-PPI and hcASD102 gene expression in BrainSpan developmental RNA-seq data.

## scripts

- `scripts/01_download_format_bsRNAseq.R`: Downloads and formats BrainSpan RNA-seq data.
- `scripts/02_asdPPI_bsRNAseq_prenatal_expression.R`: Compares prenatal BrainSpan expression for ASD-PPI bait, prey, and other HEK293T-expressed genes.
- `scripts/03_asdPPI_bsRNAseq_permutations.R`: Runs permutation analyses for ASD-PPI prey and hcASD102 genes.
- `scripts/04_make_bsRNAseq_asdPPI_hcASD102_permutation_plots.R`: Creates comparison plots from the permutation outputs.

## data

- `data/bsRNAseq.RData`: Formatted BrainSpan RNA-seq object with gene metadata, sample metadata, and expression matrix.
- `data/rawdata/`: BrainSpan metadata and documentation files used by the formatting script.

## output

- `output/asdPPI_bsRNASeq_prenatal_expression_boxplot.pdf`: Prenatal expression comparison plot.
- `output/bsRNAseq_prey100_ads102_comparison_permutationPlots_byPeriod.pdf`: Period-level comparison of ASD-PPI prey and hcASD102 enrichment patterns.
- `output/bsRNAseq_prey100_ads102_comparison_permutationPlots_byNatal.pdf`: Prenatal versus postnatal comparison.
- `output/bsRNAseq_prey100_ads102_comparison_permutationPlots_allSamples.pdf`: Sample-level comparison of ASD-PPI prey and hcASD102 enrichment patterns.
