# Shared Data

Shared ASD-PPI, hcASD102, and HEK293T proteome data used across analysis modules.

## Files

- `ppi100.RData`: Flat ASD-PPI interaction table.
- `hek293tProteome/01_define_hek293tproteome.R`: Formats the HEK293T proteome from Bekker-Jensen et al. 2017 supplementary tables.
- `hek293tProteome/BekkerJensen2017_293tProteome.RData`: HEK293T proteome object with gene annotations and protein-abundance metrics.
- `asd102/download_satterstrom_hcASD102.R`: Downloads Satterstrom et al. 2020 Table S2 and creates the hcASD102 object.
- `asd102/asd102.Rdata`: hcASD102 gene list and metadata.
- `asd102/Satterstrom_2020_Table_S2_ASD_genes.xlsx`: Satterstrom et al. 2020 Table S2.
