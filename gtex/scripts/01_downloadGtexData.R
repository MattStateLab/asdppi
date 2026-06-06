################################################################################
# Download and format GTex data
################################################################################
# Run from the repository root.
wd <- "gtex/"
dd <- "gtex/data/"
od <- "gtex/output/"


##### ----- set up workspace ----- #####
#rm(list=ls())
# Set up workspace
library(tidyverse)
library(CePa) #to import GTEx .gct files
library(Rmisc)
library(rstatix)
library(ggpubr)
library(ggrepel)
library(biomaRt)


##### ----- Download GTEx data ----- #####
# Download gene TMP expression files from GTEx website (if not previously done)
### Dataset imported here is from GTEx Analysis V8 (dbGaP Accession phs000424.v8.p2)

if(!dir.exists(paste0(dd, "rawdata/"))){
  dir.create(paste0(dd, "rawdata/"))

  # Download GTEx v8 RNAseq data (Median gene-level TPM by tissue. Median expression was calculated from the file GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz.)
  url <- "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
  file = basename (url)
  download.file(url, paste0(dd, "rawdata/", file))
}

##### ----- Import and reformat GTEX RNAseq data (median TPM) ----- #####
## Import GTEx RNAseq data (for now, focus on median TPM)
gtexMedian <- read.gct(gzfile(paste0(dd, "rawdata/","GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")))

## Annotate GTEx data - convert Ensenbl IDs to gene names
mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
geneIDs <- rownames(gtexMedian)
geneIDs <- str_replace(geneIDs, 
                       pattern = ".[0-9]+$",
                       replacement = "")
geneIDtoName <- getBM(filters = "ensembl_gene_id",
                      attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                      values = geneIDs,
                      mart = mart)

## Trim GTEx data - keep only genes that are: 1) protein coding; 2) have hgnc_symbol that is not ""
geneName <- geneIDtoName %>% 
  filter(gene_biotype == "protein_coding") %>%
  filter(hgnc_symbol != "") %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol)
ensembl_gene_id <- geneIDs
gtex0 <- cbind(ensembl_gene_id, as.data.frame(gtexMedian))
gtex <- left_join(geneName, gtex0, by = "ensembl_gene_id", values = T) 
## update HGNC symbols, remove duplicaets (MPHOSPH6)
nameChanges = data.frame(original = gtex$hgnc_symbol, updated = alias2SymbolTable(gtex$hgnc_symbol, species = "Hs")) %>% mutate_all(as.character)
nameChanges %>% filter(original != updated | is.na(updated))
gtex$hgnc_symbol = alias2SymbolTable(gtex$hgnc_symbol, species = "Hs")
gtex = gtex %>% filter(!duplicated(hgnc_symbol)) %>% filter(!is.na(hgnc_symbol))


## create gene metadata, sample metadata, and exprdata tables
genes = gtex[, 1:2] %>%
  mutate(rowNum = c(1:dim(gtex)[1])) %>%
  dplyr::select(rowNum, ensembl_gene_id, hgnc_symbol)

exprData = gtex[, -c(1:2)]
rownames(exprData) = genes$hgnc_symbol

samples = tibble(colNum = 1:dim(exprData)[2],
                sample = colnames(exprData)) %>%
  mutate(tissue = str_split(sample, "\\.", simplify = T)[,1])  %>%
  mutate(tissue = replace(tissue, tissue == "Minor", "Salivary")) %>%
  mutate(tissue = replace(tissue, tissue == "Whole", "Blood")) %>%
  mutate(tissue = replace(tissue, tissue %in% c("Esophagus", "Small", "Colon", "Pancreas", "Stomach"), "GI")) %>%
  mutate(brain = ifelse(tissue == "Brain", "brain", "notBrain"))

gtex = list("genes" = genes,
            "samples" = samples,
            "exprData" = exprData)

#save formatted GTEx data
# save data
save(gtex, file=paste0(dd, "gtex_medianTPM.RData"))
