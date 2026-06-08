# Risk Gene Enrichment

Fisher's exact test analyses for enrichment of ASD, developmental delay, schizophrenia, and SFARI risk genes among ASD-PPI prey.

## scripts

- `script/01_asdPPI_prey_riskGeneEnrichment.R`: Tests enrichment of published risk-gene sets among ASD-PPI prey.
- `script/01_asdPPI_preyOnly_riskGeneEnrichment.R`: Repeats risk-gene enrichment analyses with hcASD genes removed from the prey set.
- `script/02_asdPPI_downsampling_hcASD255discovery.R`: Downsampling analysis for ASD-PPI network size and hcASD255 discovery.
- `script/03_hcASD2022_eulerDiagram.R`: Plots overlap among hcASD gene sets from Satterstrom 2020, Fu 2022, Trost 2022, and Zhou 2022.

## data

- `data/geneAnnotation_hek293T_ASD_DD_SCZ.csv`: Gene annotations used for HEK293T expression, ASD-PPI prey status, and risk-gene membership.

## output

- `output/FET_asd255_scz34_asdPPIprey.pdf`: ASD and schizophrenia risk-gene enrichment among ASD-PPI prey.
- `output/asdPPI_downsampling_hcASD255_FET.pdf`: Downsampling summary for hcASD255 enrichment.
- `output/2022_hcASD_geneset_overlap.pdf`: hcASD gene-set overlap plot.
- `output/preyOnly/`: Enrichment outputs for prey analyses excluding hcASD genes.
