################################################################################
# Download HEK293T proteome from Bekker Jensen et al 2017  
################################################################################

# ------------------------------------------------------------------------------
#set up workspace 
# ------------------------------------------------------------------------------

#rm(list=ls())
# Set up workspace
library(tidyverse)
library(readxl)
library(biomaRt)
library(limma)

# Run from the repository root.
dd <- "data/hek293tProteome/"

# ------------------------------------------------------------------------------
#Download 293T proteome data
# ------------------------------------------------------------------------------
# Download supplementary tables from Bekker-Jensen et al 2017 (PMID: 28601559)

if(!file.exists(paste0(dd, "BekkerJensen2017_TableS7.xlsx"))){
  # Download Table S7, which contains quantitative proteomics of protein groups detected in various cell lines (including 293T)
  url <- "https://www.cell.com/cms/10.1016/j.cels.2017.05.009/attachment/aa34044f-55ec-4dc5-bcc2-9e52e4945e5f/mmc8.xlsx" 
  file = basename (url)
  download.file(url, paste0(dd, "BekkerJensen2017_TableS7.xlsx"))
  
  # Download Table S8, which contains summary statistics of proteins detected in cell lines (including 293T)
  url <- "https://www.cell.com/cms/10.1016/j.cels.2017.05.009/attachment/a4390bfc-aefa-4727-aba6-5e902248828d/mmc9.xlsx"
  file = basename (url)
  download.file(url, paste0(dd, "BekkerJensen2017_TableS8.xlsx"))
}

# ------------------------------------------------------------------------------
#Extract 293T proteome data 
# ------------------------------------------------------------------------------

# load Table S7, which contains quantitative proteomics of protein groups detected in various cell lines (including 293T)
allProteins <- read_excel(paste0(dd, "BekkerJensen2017_TableS7.xlsx"), 
                          skip = 2)

# load Table S8, which contains summary statistics of proteins detected in cell lines (including 293T)
allProteinsSummary <- read_excel(paste0(dd, "BekkerJensen2017_TableS8.xlsx"), 
                                 skip = 2)
# BW NOTE: from empirical testing, it appears that the reported number of proteins detected per cell line (11699 proteins for 293T) appears to have been calculated by averaging the number of of "Identication Type" = "By MS/MS" for E1 and E2 replicates. This was not explicitly stated in the methods section, but I have verified using data from Table S7 that this calculation results in matching numbers with Table S8.

# Data from Table S8:
# `Cell line /tissue` `MS/MS`    PSMs Peptides Proteins
# <chr>                 <dbl>   <dbl>    <dbl>    <dbl>
#   1 All HeLa            9800974 3515750   584291    14237
# 2 A549                1010861  392980   151297    11400
# 3 SH-SY5Y             1036038  466300   174432    11997
# 4 HEK293               998737  454790   169493    11699 <- want to match this number

# Checking data from Table S7 for 293T cell line
# test <- allProteins %>% select(1:5, grep("293", colnames(allProteins)))
# sum(test$`Identification type HEK293-46fracs-E1` == "By MS/MS", na.rm = TRUE)
# [1] 11700
# > sum(test$`Identification type HEK293-46fracs-E2` == "By MS/MS", na.rm = TRUE)
# > mean(c(11700, 11698))
# [1] 11699 <- MATCH

# Extract names of unique genes with protein products identified from 293T cells in one or both replicates
### Keep only subset of table pertaining to HEK293T cells
hek293t <- allProteins %>% dplyr::select(1:5, grep("293", colnames(allProteins)))

### rename columns for easier use with dplyr (remove spaces, simplify names)
colnames(hek293t) <- c("mostAbundantProteinCodingGene", "ensembl_gene_id", "proteinIDs", "proteinNames", "geneNames", 
                       "razorUniquePeptides_E1", "razorUniquePeptides_E2", 
                       "iBAQ_E1", "iBAQ_E2", 
                       "LFQ_E1", "LFQ_E2", 
                       "IDtype_E1", "IDtype_E2")
# Format data
hek293t =  hek293t %>%
  # Keep only genes that were identified "By MS/MS" in at least one replicate
  filter(IDtype_E1 == "By MS/MS" | IDtype_E2 == "By MS/MS")  %>%   
  mutate(ensembl_gene_id = strsplit(ensembl_gene_id, ";")) %>% 
  unnest(ensembl_gene_id)

### calculate average iBAQ and LFQ values between replicates
      # LFQ (Label-free quantification) intensities are based on the (raw) intensities and normalized on multiple levels to make sure that profiles of LFQ intensities across samples accurately reflect the relative amounts of the proteins.
      # iBAQ (Intensity Based Absolute Quantification) values calculated by MaxQuant are the (raw) intensities divided by the number of theoretical peptides. Thus, iBAQ values are proportional to the molar quantities of the proteins. The iBAQ algorithm can roughly estimate the relative abundance of the proteins within each sample.””
hek293t <- hek293t %>%
  mutate(iBAQ_avg = apply(hek293t[ , c("iBAQ_E1", "iBAQ_E2")], 1, mean)) %>%
  mutate(LFQ_avg = apply(hek293t[ , c("LFQ_E1", "LFQ_E2")], 1, mean))

### match Ensembl gene IDs to HGNC symbol
mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
geneIDtoName <- getBM(filters = "ensembl_gene_id",
                      attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      values = hek293t$ensembl_gene_id,
                      mart = mart,
                      useCache=F)
### annotate data with gene names
hek293tProteome = left_join(hek293t, geneIDtoName, by= "ensembl_gene_id")


### keep rows with unique ensembl_gene_id and hgnc_symbol
hek293tProteome <- hek293tProteome %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) %>%
  filter(!(is.na(hgnc_symbol))) %>%   #remove rows (#1) without hgnc_symbol
  dplyr::select(ensembl_gene_id, hgnc_symbol, iBAQ_avg, LFQ_avg)  #trim columns

### update HGNC symbols
hek293tProteome = hek293tProteome %>% mutate(hgnc_symbol = alias2SymbolTable(hgnc_symbol)) %>% filter(!is.na(hgnc_symbol))


#save data
save(hek293tProteome,
     file=paste0(dd, "BekkerJensen2017_293tProteome.RData"))
