################################################################################
# Assess median expression percentile of ASD-PPI genes in bsRNAseq prenatal brain tissue
################################################################################

# ------------------------------------------------------------------------------
#set up workspace 
# ------------------------------------------------------------------------------
# Run from the repository root.
wd <- "bsRNAseq/"
dd <- "bsRNAseq/data/"
od <- "bsRNAseq/output/"
rd <- "data/"

#rm(list=ls())
# Set up workspace
library(tidyverse)
library(limma) 
library(parallel)
library(ggpubr)


# import Brainspan Data
### downloaded and formatted in BWASD02_downloadBrainspanDevelopmenalData.R
### bsRNAseq = brainspan RNAseq with 524 samples from Brainspan website
load(file = paste0(dd, "bsRNAseq.RData"))

# import 293T proteome data        
### data from from Bekker-Jensen et al 2017 (PMID: 28601559)
### 293T-specific data extracted and formatted in 01_define_hek293tproteome.R
load(paste0(rd, "hek293tProteome/BekkerJensen2017_293tProteome.RData"))
hek293tProteome_o = hek293tProteome

# Import PPI data 
load(paste0(rd, "ppi100.RData"))
ppiFull = ppi100

# ------------------------------------------------------------------------------
# Format and annotate data 
# ------------------------------------------------------------------------------

# Extract bait and prey genes (defining prey genes to be exclusive of bait genes)
bait <- unique(ppiFull$bait)
prey <- setdiff(unique(ppiFull$prey), bait)

# Reformat BrainSpan data by converting expression level to rank within each sample (where higher rank reflects higher relative gene expression)
bsRNAseq[["exprData_ranked"]] = apply(bsRNAseq$exprData, 2, function(x) rank(x, ties.method = "average"))

# Annotate 293T proteome by:
# 1) ASD-PPI bait
# 2) ASD-PPI prey
# 3) In BrainSpan microarray data
# 4) in BrainSpan RNAseq data
hek293tProteome <- hek293tProteome_o %>%
  mutate(asdPPI = hgnc_symbol) %>%
  mutate(asdPPI = replace(asdPPI, hgnc_symbol %in% prey, "prey")) %>%
  mutate(asdPPI = replace(asdPPI, hgnc_symbol %in% bait, "bait")) %>%
  mutate(asdPPI = replace(asdPPI, !asdPPI %in% c("prey", "bait"), "other")) %>%
  mutate(bsRNAseqStatus = hgnc_symbol) %>%
  mutate(bsRNAseqStatus = ifelse(bsRNAseqStatus %in% bsRNAseq[["genes"]][, "hgnc_symbol"], "present", "absent")) %>%
  filter(hgnc_symbol != "") # remove genes that do not have hugo gene name (ex:ENSG00000187186)

bait293T = hek293tProteome %>% filter(asdPPI == "bait") %>% .$hgnc_symbol
prey293T = hek293tProteome %>% filter(asdPPI == "prey") %>% .$hgnc_symbol 
other293T = hek293tProteome %>% filter(asdPPI == "other") %>% .$hgnc_symbol


# ------------------------------------------------------------------------------
# DEFINE FUNCTION
# ------------------------------------------------------------------------------

# Function to calculate enrichment of ASD-PPI in brainspan expression data
### Inputs:
# 1. bsData = expression data list (with gene info, sample info, and exprData)
# 2. proteomeData: hek293T proteome genes annotated by hgnc_symbol (matching BrainSpan annotations), asdPPI gene category, and whether measured in BrainSpan data
### Output: summary table of sample (rows) and median geneset rank + median geneset percentile expression for bait geneset, prey geneset, bait-prey geneset
hek293tASDppiRank <- function(bsData = bsRNAseq,
                                  proteomeData = hek293tProteome){
  bsExprData = bsData$exprData
  bsExprData_ranked =  apply(bsExprData, 2, function(x) rank(x, ties.method = "average")); rownames(bsExprData_ranked) = rownames(bsExprData)
  bsExprData_percentile = apply(bsExprData, 2, dplyr::percent_rank); rownames(bsExprData_percentile) = rownames(bsExprData)
  bsGenes = bsData$genes %>%
    filter(!is.na(hgnc_symbol)) %>% .$hgnc_symbol
  
  #define genesets
  hek293T_bait = proteomeData$hgnc_symbol[which(proteomeData$asdPPI == "bait")]
  hek293T_prey = proteomeData$hgnc_symbol[which(proteomeData$asdPPI == "prey")]
  hek293T_baitprey = c(hek293T_bait, hek293T_prey)
  hek293T_other = proteomeData$hgnc_symbol[which(proteomeData$asdPPI == "other")]
  notHEK293T = bsGenes[!(bsGenes %in% proteomeData$hgnc_symbol)]
  allHEK293T = bsGenes[(bsGenes %in% proteomeData$hgnc_symbol)]
  allBrainspan = bsGenes
  
  #define numbers for wStat normalization
  nbait = length(proteomeData$hgnc_symbol[which(proteomeData$asdPPI == "bait")])
  nprey = length(proteomeData$hgnc_symbol[which(proteomeData$asdPPI == "prey")])
  nbaitprey = nbait + nprey
  nother = length(proteomeData$hgnc_symbol[which(proteomeData$asdPPI == "other")])
  
  # calculate median geneset rank
  medRank_bait = apply(bsExprData_ranked[which(rownames(bsExprData_ranked) %in% hek293T_bait),], 2, median)
  medRank_prey = apply(bsExprData_ranked[which(rownames(bsExprData_ranked) %in% hek293T_prey),], 2, median)
  medRank_baitprey= apply(bsExprData_ranked[which(rownames(bsExprData_ranked) %in% hek293T_baitprey),], 2, median)
  medRank_other= apply(bsExprData_ranked[which(rownames(bsExprData_ranked) %in% hek293T_other),], 2, median)
  medRank_notHEK293T = apply(bsExprData_ranked[which(rownames(bsExprData_ranked) %in% notHEK293T),], 2, median)
  medRank_allHEK293T = apply(bsExprData_ranked[which(rownames(bsExprData_ranked) %in% allHEK293T),], 2, median)
  medRank_allBrainspan = apply(bsExprData_ranked[which(rownames(bsExprData_ranked) %in% allBrainspan),], 2, median)
  
  # calculate median geneset percentile expression
  medPercentile_bait = apply(bsExprData_percentile[which(rownames(bsExprData_percentile) %in% hek293T_bait),], 2, median)
  medPercentile_prey = apply(bsExprData_percentile[which(rownames(bsExprData_percentile) %in% hek293T_prey),], 2, median)
  medPercentile_baitprey= apply(bsExprData_percentile[which(rownames(bsExprData_percentile) %in% hek293T_baitprey),], 2, median)
  medPercentile_other= apply(bsExprData_percentile[which(rownames(bsExprData_percentile) %in% hek293T_other),], 2, median)
  medPercentile_notHEK293T = apply(bsExprData_percentile[which(rownames(bsExprData_percentile) %in% notHEK293T),], 2, median)
  medPercentile_allHEK293T = apply(bsExprData_percentile[which(rownames(bsExprData_percentile) %in% allHEK293T),], 2, median)
  medPercentile_allBrainspan = apply(bsExprData_percentile[which(rownames(bsExprData_percentile) %in% allBrainspan),], 2, median)
  
  
  summary = tibble(samples = colnames(bsExprData),
                   regionList = bsData$samples$regionList,
                   broadRegion = bsData$samples$broadRegion,
                   ncxRegion = bsData$samples$ncxRegion,
                   period = bsData$samples$period,
                   natal = bsData$samples$natal,
                   medRank_bait = medRank_bait,
                   medPercentile_bait = medPercentile_bait,
                   medRank_prey = medRank_prey,
                   medPercentile_prey = medPercentile_prey, 
                   medRank_baitprey = medRank_baitprey,
                   medPercentile_baitprey = medPercentile_baitprey,
                   medRank_other = medRank_other,
                   medPercentile_other = medPercentile_other,
                   medRank_notHEK293T = medRank_notHEK293T,
                   medRank_allHEK293T = medRank_allHEK293T,
                   medRank_allBrainspan = medRank_allBrainspan, 
                   medPercentile_notHEK293T = medPercentile_notHEK293T,
                   medPercentile_allHEK293T = medPercentile_allHEK293T,
                   medPercentile_allBrainspan = medPercentile_allBrainspan)
  
  return(summary)
}


# ------------------------------------------------------------------------------
#  Assess relative expression of 293T bait, prey, and other genes in BrainSpan prenatal samples
# ------------------------------------------------------------------------------

# Calculate Wilcoxon rank sum statistic and associated p values and pVal values for 293T bait, prey, and baitprey geneset compared to 293T background ("other") genes in BrainSpan prenatal samples
bsRNAseq_WRS = hek293tASDppiRank(bsData = bsRNAseq,
                                     proteomeData = hek293tProteome)

# Format data for plotting (pivot longer)
bsRNAseq_WRS_toPlot = bsRNAseq_WRS %>% pivot_longer(
  cols = medRank_bait:medPercentile_allBrainspan,
  names_to = c(".value", "geneset"),
  names_sep = "_")
bsRNAseq_WRS_toPlot = bsRNAseq_WRS_toPlot %>%
  mutate(geneset = factor(geneset, levels = c("bait", "prey", "baitprey", "other", "allHEK293T", "notHEK293T", "allBrainspan"))) %>%
  mutate(natal = factor(natal, levels = c("prenatal", "postnatal"))) 

# make plot
datasetName="bsRNAseq"
bsRNAseq_WRS_toPlot %>% filter(natal == "prenatal",
                         geneset %in% c("bait", "prey", "other")) %>%
  mutate(geneset = recode_factor(geneset, 'bait'='Bait', 'prey'='Prey (- hcASD102)', 'other'="Other")) %>%
  ggplot(aes(x = geneset, y = medPercentile * 100, fill = geneset)) +
  geom_boxplot(show.legend = F) +
  #geom_jitter(size = 1, alpha = 0.15) +
  stat_compare_means(comparisons=list(c("Bait", "Prey (- hcASD102)"), c("Prey (- hcASD102)", "Other"), c("Bait", "Other")), method = "t.test", aes(label=..p.adj..)) +
  labs(title = paste0("Median geneset expression (percentile)\n", datasetName, ", prenatal samples"), 
       x="Geneset", y="Median expression (percentile)")+
  theme_bw()
# ggsave(paste0(od, "asdPPI_bsRNASeq_prenatal_expression_boxplot.pdf"), width = 5, height=6)
ggsave(paste0(od, "asdPPI_bsRNASeq_prenatal_expression_boxplot.pdf"), width = 3.5, height=6)


#t.test findings for supplementary table
group1 = bsRNAseq_WRS_toPlot %>% filter(natal=="prenatal" & geneset=="bait") %>% pull(medPercentile)
group2 = bsRNAseq_WRS_toPlot %>% filter(natal=="prenatal" & geneset=="prey") %>% pull(medPercentile)
group3 = bsRNAseq_WRS_toPlot %>% filter(natal=="prenatal" & geneset=="other") %>% pull(medPercentile) 
t.test(group1, group2)
t.test(group1, group3)
t.test(group2, group3)
