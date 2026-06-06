################################################################################
# Download and format Brainspan Developmental RNAseq data
###############################################################################

# ------------------------------------------------------------------------------
#set up workspace 
# ------------------------------------------------------------------------------
# Run from the repository root.
wd <- "bsRNAseq/"
dd <- "bsRNAseq/data/"
od <- "bsRNAseq/output/"


#rm(list=ls())
# Set up workspace
library(tidyverse)
library(data.table)
library(readxl)
library(limma)
library(ggpubr)
library(biomaRt)

# ------------------------------------------------------------------------------
# 1) Download BrainSpan developmental RNAseq data (from Brainspan Website)
# ------------------------------------------------------------------------------

## Download RNASEQ data
## source: http://www.brainspan.org/static/download.html
## This data set contains normalized data from original RPKM (reads per kilobase per million; see the whitepaper at www.brainspan.org). Values averaged to genes.
## The data folder contains 3 separate files
# 1. expression_matrix.csv -- the rows are genes and the columns samples; the first column is the row number 
# 2. rows_metadata.csv -- the genes are listed in the same order as the rows in expression_matrix.csv
# 3. columns_metadata.csv -- the samples are listed in the same order as the columns in expression_matrix
if(!dir.exists(paste0(dd, "/rawdata"))){
  dir.create(paste0(dd, "/rawdata"), recursive = TRUE)
  #Download RNAseq summarized to genes
  download.file("https://www.brainspan.org/api/v2/well_known_file_download/267666525",
                paste0(dd, "/rawdata/brainspan_rnaseq_genes.zip"))
  unzip(paste0(dd, "/rawdata/brainspan_rnaseq_genes.zip"), exdir = paste0(dd, "/rawdata"))
}

# ------------------------------------------------------------------------------
# 2) FORMAT BRAINSPAN RNAseq DATA
# ------------------------------------------------------------------------------
# IMPORT FILES
exprData_o = read.csv(paste0(dd, "rawdata/expression_matrix.csv"), 
                      header = F)
genes_o = read.csv(paste0(dd, "rawdata/rows_metadata.csv"), stringsAsFactors = F)
samples_o = read.csv(paste0(dd, "rawdata/columns_metadata.csv"), stringsAsFactors = F)


# FORMAT GENES METADATA (goes from 52376 genes to 20430 genes)
### update gene metadata with latest hgnc_symbol
mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
geneIDtoName <- getBM(filters = "entrezgene_id",
                      attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "gene_biotype"),
                      values = genes_o$entrez_id,
                      mart = mart) 

### keep only genes that meet the following criteria
# 1. protein coding (removed #34040)
# 2. have matching hgnc_symbol (removed #0)
# 3. hgnc_symbol is not "" (removed #87)
# 4. duplicated entrez_id (mostly pseudogenes, removed #2505)
# 5. duplicated ensembl_gene_id (removed #31)
# 6. duplicated hcnc_symbol (removed #1)

genes <- left_join(genes_o, geneIDtoName, by = c("entrez_id" = "entrezgene_id")) %>%
  # keep only protein-coding genes
  filter(gene_biotype == "protein_coding") %>%  #21176
  filter(!is.na(hgnc_symbol)) %>% #21176
  filter(hgnc_symbol != "") %>% #21089
  filter(!duplicated(entrez_id)) %>% #18584
  filter(!is.na(ensembl_gene_id.y)) %>% #18584
  filter(!duplicated(ensembl_gene_id.y)) %>% #18553
  filter(!duplicated(hgnc_symbol)) %>% #18552 
  dplyr::select(row_num, gene_id, ensembl_gene_id_original = ensembl_gene_id.x, ensembl_gene_id_updated = ensembl_gene_id.y,
                gene_symbol_original = gene_symbol, entrez_id, hgnc_symbol, gene_biotype)

# > dim(genes_o)
# [1] 52376     5
# > dim(genes)
# [1] 18552     8

# Check whether gene symbols used in PPI should be udpated
### yes - there are 1 genes with more uptodate HGNC gene symbols
nameChanges = data.frame(original = genes$hgnc_symbol, updated = alias2SymbolTable(genes$hgnc_symbol, species = "Hs")) %>% mutate_all(as.character)
nameChanges %>% filter(original != updated | is.na(updated))
genes$hgnc_symbol = alias2SymbolTable(genes$hgnc_symbol, species = "Hs")


# FORMAT SAMPLES METADATA
### annotate Brainspan microarray data by age/period
samples = samples_o %>% mutate(period = as.character(age)) %>% 
  mutate(period = replace(period, period == "8 pcw", 2),
         period = replace(period, period == "9 pcw", 2),
         period = replace(period, period == "12 pcw", 3),
         period = replace(period, period == "13 pcw", 4),
         period = replace(period, period == "16 pcw", 5),
         period = replace(period, period == "17 pcw", 5),
         period = replace(period, period == "19 pcw", 6),
         period = replace(period, period == "21 pcw", 6),
         period = replace(period, period == "24 pcw", 7),
         period = replace(period, period == "25 pcw", 7),
         period = replace(period, period == "26 pcw", 7),
         period = replace(period, period == "35 pcw", 7),
         period = replace(period, period == "37 pcw", 7),
         period = replace(period, period == "4 mos", 8),
         period = replace(period, period == "10 mos", 9),
         period = replace(period, period == "1 yrs", 10),
         period = replace(period, period == "2 yrs", 10),
         period = replace(period, period == "3 yrs", 10),
         period = replace(period, period == "4 yrs", 10),
         period = replace(period, period == "8 yrs", 11),
         period = replace(period, period == "11 yrs", 11),
         period = replace(period, period == "13 yrs", 12),
         period = replace(period, period == "15 yrs", 12),
         period = replace(period, period == "18 yrs", 12),
         period = replace(period, period == "19 yrs", 12),
         period = replace(period, period == "21 yrs", 13),
         period = replace(period, period == "23 yrs", 13),
         period = replace(period, period == "30 yrs", 13),
         period = replace(period, period == "36 yrs", 13),
         period = replace(period, period == "37 yrs", 13),
         period = replace(period, period == "40 yrs", 14)) %>%
  mutate(period = as.numeric(period)) %>%
  mutate(brain = donor_id) %>%
  mutate(region = structure_acronym) %>%
  mutate("natal" = ifelse(period <= 7, "prenatal", "postnatal")) %>%
  mutate("natal" = factor(natal, levels = c("prenatal", "postnatal"))) %>%
  mutate(regionList = region) %>%
  mutate(regionList = replace(regionList, regionList == "DFC", "dfcList")) %>%
  mutate(regionList = replace(regionList, regionList == "VFC", "vfcList")) %>%
  mutate(regionList = replace(regionList, regionList == "MFC", "mfcList")) %>%
  mutate(regionList = replace(regionList, regionList == "OFC", "ofcList")) %>%
  mutate(regionList = replace(regionList, regionList %in% c("M1C", "M1C-S1C"), "m1cList")) %>%
  mutate(regionList = replace(regionList, regionList == "S1C", "s1cList")) %>%
  mutate(regionList = replace(regionList, regionList %in% c("IPC", "PCx"), "ipcList")) %>%
  mutate(regionList = replace(regionList, regionList == "A1C", "a1cList")) %>%
  mutate(regionList = replace(regionList, regionList %in% c("STC", "TCx"), "stcList")) %>%
  mutate(regionList = replace(regionList, regionList == "ITC", "itcList")) %>%
  mutate(regionList = replace(regionList, regionList %in% c("V1C", "Ocx"), "v1cList")) %>%
  mutate(regionList = replace(regionList, regionList == "HIP", "hipList")) %>%
  mutate(regionList = replace(regionList, regionList == "AMY", "amyList")) %>%
  mutate(regionList = replace(regionList, regionList %in% c("STR", "MGE", "LGE", "CGE"), "strList")) %>%
  mutate(regionList = replace(regionList, regionList %in% c("DTH", "MD"), "thmList")) %>%
  mutate(regionList = replace(regionList, regionList %in% c("CBC", "CB", "URL"), "cblList"))%>%
  mutate(broadRegion = regionList) %>%
  mutate(broadRegion = replace(broadRegion, broadRegion %in% c("dfcList", "vfcList", "mfcList", "ofcList","m1cList", "s1cList"), "PFC-MSC")) %>%
  mutate(broadRegion = replace(broadRegion, broadRegion %in% c("ipcList", "a1cList", "stcList", "itcList","v1cList"), "IPC-V1C")) %>%
  mutate(broadRegion = replace(broadRegion, broadRegion %in% c("hipList", "amyList", "strList"), "HIP-STR")) %>%
  mutate(broadRegion = replace(broadRegion, broadRegion %in% c("thmList", "cblList"), "MD-CBC")) %>%
  mutate(ncxRegion = regionList) %>%
  mutate(ncxRegion = replace(ncxRegion, ncxRegion %in% c("dfcList", "vfcList", "mfcList", "ofcList","m1cList", "s1cList", "ipcList", "a1cList","stcList", "itcList", "v1cList"), "NCX")) %>%
  mutate(ncxRegion = replace(ncxRegion, ncxRegion %in% c("hipList", "amyList", "strList", "thmList","cblList"), "NON-NCX")) %>%
  mutate(donor = donor_name) %>%
  dplyr::select(brain, donor, gender, age, period, natal, region, regionList, broadRegion, ncxRegion, structure_id, structure_name)

# FORMAT EXPRESSION DATA
### trim expression data to selected genes, remove first column (which just contains row numbers)
exprData = exprData_o[genes$row_num, -1]
rownames(exprData) = genes$hgnc_symbol
colnames(exprData) = 1:ncol(exprData)

bsRNAseq = list("samples" = samples,
                "genes" = genes,
                "exprData" = exprData)

#save data
save(bsRNAseq,
     file = paste0(dd, "bsRNAseq.RData"))
