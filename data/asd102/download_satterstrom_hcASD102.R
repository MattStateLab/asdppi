#Purpose: download Satterstrom 2020 Supplemntary Table 2 and extract 102 hcASD genes

# Run from the repository root.
dd <- "data/asd102/"

#set up workspace
library(tidyverse)
library(readxl)

# download Table S2 from Satterstrom et al 2020 (PMID: 31981491)
if(!file.exists(paste0(dd, "Satterstrom_2020_Table_S2_ASD_genes.xlsx"))){
  download.file(url = "https://www.cell.com/cms/10.1016/j.cell.2019.12.036/attachment/cc634b34-b459-4aff-acf0-cb58a3013bd1/mmc2.xlsx",
                destfile = paste0(dd, "Satterstrom_2020_Table_S2_ASD_genes.xlsx"))
}


# Create list of 102 hcASD genes 
satterstrom_asd <- read_excel(path = paste0(dd, "Satterstrom_2020_Table_S2_ASD_genes.xlsx"), sheet = 2)
asd102 = list()
asd102[["hgnc"]] = satterstrom_asd %>% filter(qval_dnccPTV < 0.1) %>% pull(hugoGene)
asd102[["genesMetadata"]] = satterstrom_asd %>% filter(qval_dnccPTV < 0.1) %>% dplyr::select(gene, hugoGene, hgnc_id, entrez_id, ensembl_gene_id, refseq_accession,	uniprot_ids, location, chr)

#Save as Rdata object
save(asd102, file=paste0(dd, "asd102.Rdata"))
