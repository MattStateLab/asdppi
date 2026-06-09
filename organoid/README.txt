# FOXP1 Organoid Scripts

Scripts for FOXP1 organoid PIP-seq, CUT&Tag, and IHC analyses.

## scripts

- `FOXP1-PIPseq-script_R513H.R`: Processes FOXP1 R513H/R514H PIP-seq data with Seurat, including QC, clustering, cell-type annotation, differential expression, and coexpression summaries.
- `FOXP1-PIPseq-script_L327P_v2.R`: Processes FOXP1 L327P PIP-seq data and compares L327P and R513H/R514H differential-expression patterns.
- `FOXP1_fig8_stats_LME_v3.R`: Runs linear model and mixed-effects model statistics for figure 8 marker-quantification datasets.
- `FOXP1_fig8_stats_LME_v4_permutationBW.R`: Runs updated figure 8 model statistics with permutation testing.
- `Diffpeaks_anno_FOXP1_R514H.R`: Performs CUT&Tag peak QC, annotation, GREAT enrichment, counting, and differential peak analysis.
- `cuttag_FOXP1.sh`: SLURM shell pipeline for FOXP1 CUT&Tag trimming, alignment, duplicate marking, and peak calling.
- `organoid_MEA_analysis`: Runs kilosort2 followed by mean firing rate analysis and raster plotting on sorted electrode data, as well as burst, ISI, and STTC analysis.
