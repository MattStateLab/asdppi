################################################################################
# Assess constraint metrics for ASD PPI genes
################################################################################
# ------------------------------------------------------------------------------
#set up workspace 
# ------------------------------------------------------------------------------
#rm(list=ls())
# Set up workspace
library(tidyverse)
library(data.table)
library(readxl)
library(limma)
library(ggpubr)

# Run from the repository root.
wd <- "ssc_enrichment/"
dd <- "ssc_enrichment/data/"
od <- "ssc_enrichment/output/"

#import genesAnnotation table that contains Satterstrom 2020 genes HGNC symbols annotated with HEKproteome status, ASD-PPI bait and prey status, and number of individuals in case-control or SSC simplex data with different types of genetic variants
geneAnnotations = read.csv(paste0(od, "Satterstrom2020_geneAnnotations_SSC_caseControl_asdPPI100.csv"))


#annotate genes with ExAC pLI, misZ, synZ
###Download ExAC gene constraint scores TSV: "The ExAC data set contains data from 60,706 exomes, all mapped to the GRCh37/hg19 reference sequence."
if(!file.exists(paste0(dd, "forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz"))){
  url = "https://storage.googleapis.com/gcp-public-data--gnomad/legacy/exac_browser/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz"
  file = basename(url)
  download.file(url, paste0(dd, file))
}
exac = fread(paste0(dd, "forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz"))
### update gene names
exac = exac %>% mutate(hgnc_symbol = alias2SymbolTable(gene, species = "Hs" ))
### add ExAC data to geneAnnotations
geneAnnotations_exac = left_join(geneAnnotations, dplyr::distinct(exac, hgnc_symbol, .keep_all = T), by = "hgnc_symbol")


#annotate genes with s_het scores
###download and import s_het data (Cassa et al 2017, PMID 28369035, Table S1).  
if(!file.exists(paste0(dd, "Cassa2017_Table_S1_shet.xlsx"))){
  url = "https://static-content.springer.com/esm/art%3A10.1038%2Fng.3831/MediaObjects/41588_2017_BFng3831_MOESM71_ESM.xlsx"
  download.file(url, paste0(dd, "Cassa2017_Table_S1_shet.xlsx"))
}
shet = read_excel(paste0(dd, "Cassa2017_Table_S1_shet.xlsx"))
### update gene names
shet = shet %>% mutate(hgnc_symbol = alias2SymbolTable(gene_symbol, species = "Hs" ))
###Add s_het data to gene annotations
geneAnnotations_exac_shet = shet %>% dplyr::select(s_het, hgnc_symbol) %>% left_join(geneAnnotations_exac, ., by="hgnc_symbol")


#make plot
###Compare pLI, misZ, synZ and shet of ASD-PPI bait, prey, and 293T proteome  
toPlot = geneAnnotations_exac_shet %>% dplyr::select(hgnc_symbol, hek293TproteinExpr, bait, prey, ASD102, pLI, mis_z, s_het, syn_z) %>%
  mutate(geneset = "other") %>%
  mutate(geneset = replace(geneset, prey ==1, "prey")) %>%
  mutate(geneset = replace(geneset, ASD102==1, "ASD102")) %>%
  mutate(asdPPI = "NA") %>%
  mutate(asdPPI = replace(asdPPI, hek293TproteinExpr==1, "other")) %>%
  mutate(asdPPI = replace(asdPPI, prey ==1, "prey")) %>%
  mutate(asdPPI = replace(asdPPI, bait ==1, "bait")) %>%
  mutate(asdPPI = factor(asdPPI, levels = c("bait", "prey", "other"))) %>%
  mutate(geneset = factor(geneset, levels = c("prey", "ASD102", "other"))) %>%
  filter(!is.na(asdPPI))

###Calculate median s_het
toPlot %>% group_by(asdPPI) %>% dplyr::summarize(median.shet = median(s_het, na.rm = T))

### plot significance
### note: need to adjust for multiple hypothesis testing afterwards
fig = toPlot %>% pivot_longer(pLI:syn_z, values_to = "Score", names_to="Metric") %>%
  mutate(asdPPI = recode_factor(asdPPI, 'bait'='Bait', 'prey'='Prey (- hcASD102)', 'other'="Other")) %>%
  mutate(Metric = factor(Metric, levels = c("pLI", "mis_z", "s_het","syn_z"))) %>%
  ggplot(aes(x = asdPPI, y = Score, fill = asdPPI)) + facet_wrap(~Metric, ncol=4, scale="free") +
  geom_boxplot(show.legend = F) +
  #geom_jitter(size=0.5, alpha=0.05) +
  stat_compare_means(comparisons=list(c("Bait", "Prey (- hcASD102)"), c("Prey (- hcASD102)", "Other"), c("Bait", "Other")), method = "t.test") +
# stat_compare_means(comparisons=list(c("Bait", "Prey (- hcASD102)"), c("Prey (- hcASD102)", "Other"), c("Bait", "Other")), method = "t.test", label = "p.signif") +
  labs(y = "pLI", x = "Gene category",
       title = "ASD-PPI: pLI, s_het, mis_z, and syn_z scores") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
fig
ggsave(paste0(od, "asdPPI_constraintMetrics_boxplots.pdf"), width = 7, height = 7)
# ggsave(paste0(od, "asdPPI_constraintMetrics_boxplots.pdf"), width = 8, height = 7)

## get data for supplemental tables
# toPlot %>% filter(!is.na(syn_z)) %>% pull(asdPPI) %>% table
# b = toPlot %>% filter(asdPPI=="bait") %>% pull(syn_z)
# p = toPlot %>% filter(asdPPI=="prey") %>% pull(syn_z)
# o= toPlot %>% filter(asdPPI=="other") %>% pull(syn_z)
# 
# t.test(b,p)
# t.test(b,o)
# t.test(p,o)

# **mean pLI for bait > prey > other 293T proteome genes** The pLI scores indicate the tolerance level of a given gene to loss of function (LoF) based on the number of protein-truncating variants (Lek et al 2016). Prior studies have indicated that ASD-associated genes tend to have higher pLI (citations). We compared the pLI scores of ASD-PPI bait genes and prey genes with the pLI score of the remainder of the 293T proteome (Bekker-Jensen 2019). We confirmed that the pLI score of ASD-PPI bait genes is higher than expected, and found that the mean pLI score of prey genes is also significantly higher than that of other 293T-expressed proteins. 
# 
# **mean misZ for bait > prey > other 293T proteome genes** the missenze Z scores (misZ) captures the number of observed missense variants in a gene compared to the expected number of missense variants in the general population. Higher (more positive) Z scores indicate that the transcript is more intolerant of variation (more constrained).  
# 
# **mean s_het for bait > prey > other 293T proteome genes** s_het estimates the selection against heterozygous loss of gene function. Higher s_het reflects greater intolerance heterozygous PTVs. 

