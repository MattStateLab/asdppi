# SSC Enrichment

Analyses of ASD-PPI prey enrichment for de novo variant burden and constraint metrics using Satterstrom et al. 2020 SSC data.

## Scripts

- `scripts/01_format_satterstrom2020_SSC_variants.R`: Downloads and formats Satterstrom et al. 2020 supplementary variant and gene tables.
- `scripts/02_eval_asdPPI_geneConstraintMetrics.R`: Compares pLI, missense Z, synonymous Z, and s_het metrics for ASD-PPI bait, prey, and other HEK293T proteins.
- `scripts/03_asdPPI_geneticRisk_Satterstrom2020_SSC.R`: Tests whether ASD-PPI prey are enriched for damaging variants in ASD probands versus control siblings.
- `scripts/04_asdPPI_prey_ASDrisk_downsamplingAnalysis.R`: Downsampling analysis for bait count and ASD de novo variant enrichment.

## Data

- `data/Satterstrom_2020_Table_S1_ASD_denovo_genes.xlsx`: Satterstrom et al. 2020 de novo variant table.
- `data/Satterstrom_2020_Table_S2_ASD_genes.xlsx`: Satterstrom et al. 2020 TADA gene table.
- `data/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz`: ExAC constraint metrics.
- `data/Cassa2017_Table_S1_shet.xlsx`: Cassa et al. 2017 s_het metrics.

## Outputs

- `output/Satterstrom2020_SSC_DNM.csv`: Formatted SSC de novo mutation annotations.
- `output/Satterstrom2020_geneAnnotations_SSC_caseControl_asdPPI100.csv`: Gene-level ASD-PPI and SSC variant annotations.
- `output/asdPPI_constraintMetrics_boxplots.pdf`: Constraint metric comparison plot.
- `output/FET_SSCsubset_exome.pdf`: Enrichment tests using exome-wide gene sets.
- `output/FET_SSCsubset_HEK293T_preyOnly.pdf`: Enrichment tests for HEK293T and ASD-PPI prey categories.
- `output/FET_SSCsubset_downsampling_asdPPI.pdf`: Downsampling enrichment plot.
