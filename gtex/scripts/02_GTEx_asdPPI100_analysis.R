################################################################################
# Assess whether ASD PPI prey genes are more enriched than expected in the brain (GTEx median data)

################################################################################
# ------------------------------------------------------------------------------
#set up workspace 
# ------------------------------------------------------------------------------

# Run from the repository root.
wd <- ""
dd <- "gtex/data/"
od <- "gtex/output/"

#rm(list=ls())
# Set up workspace
library(tidyverse)
library(Rmisc)
library(rstatix)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(parallel)

# import GTEx and PPI data
### Import GTEx data (downloaded and formatted in 01_downloadGtexData.R
load(paste0(dd, "/gtex_medianTPM.RData"))

# Import PPI data 
load(paste0(wd, "data/ppi100.RData"))
ppiFull = ppi100

# import 293T proteome data        
### data from from Bekker-Jensen et al 2017 (PMID: 28601559)
### HEK293T-specific data extracted and formatted in 01_define_hek293tproteome.R)
load(paste0(wd, "data/hek293tProteome/BekkerJensen2017_293tProteome.RData"))
hek293tProteome_o = hek293tProteome

# ------------------------------------------------------------------------------
# Format and annotate GTEx and PPI data 
# ------------------------------------------------------------------------------

# identify bait and prey genes
bait <- unique(ppiFull$bait)
prey <- setdiff(unique(ppiFull$prey), bait)

# Reformat GTEx data by converting expression level to rank within each sample (where higher rank reflects higher relative gene expression)
gtex[["exprData_ranked"]] = apply(gtex$exprData, 2, function(x) rank(x, ties.method = "average"))

# Annotate 293T proteome by:
# 1) ASD-PPI bait
# 2) ASD-PPI prey
# 3) In GTEx RNAseq data
hek293tProteome <- hek293tProteome_o %>%
  mutate(asdPPI = hgnc_symbol) %>%
  mutate(asdPPI = replace(asdPPI, hgnc_symbol %in% prey, "prey")) %>%
  mutate(asdPPI = replace(asdPPI, hgnc_symbol %in% bait, "bait")) %>%
  mutate(asdPPI = replace(asdPPI, !asdPPI %in% c("prey", "bait"), "other")) %>%
  mutate(gtexStatus = hgnc_symbol) %>%
  mutate(gtexStatus = ifelse(gtexStatus %in% gtex[["genes"]][, "hgnc_symbol"], "present", "absent")) %>%
  filter(hgnc_symbol != "") # remove genes that do not have hugo gene name (ex:ENSG00000233757)

bait293T = hek293tProteome %>% filter(asdPPI == "bait") %>% .$hgnc_symbol
prey293T = hek293tProteome %>% filter(asdPPI == "prey") %>% .$hgnc_symbol 
other293T = hek293tProteome %>% filter(asdPPI == "other") %>% .$hgnc_symbol

# ------------------------------------------------------------------------------
# DEFINE FUNCTIONS
# ------------------------------------------------------------------------------

# Function to perform Wilcoxon ranksum to calculate enrichment of ASD-PPI in GTEx expression data
### Inputs:
# 1. inputData = expression data list (with gene info, sample info, and exprData)
# 2. proteomeData: hek293T proteome genes annotated by hgnc_symbol (matching GTEx annotations), asdPPI gene category, and whether measured in GTEx RNAseq
### Output: summary table of sample (rows) and median rank and median percentiles for bait geneset, prey geneset, bait-prey geneset
hek293tASDppiStats <- function(inputData = gtex,
                                  proteomeData = hek293tProteome){
  exprData = inputData$exprData
  exprData_ranked = inputData$exprData_ranked
  exprData_percentile = apply(exprData, 2, dplyr::percent_rank)
  
  inputGenes = inputData$genes %>%
    filter(!is.na(hgnc_symbol)) %>% .$hgnc_symbol
  
  #define genesets
  hek293T_bait = proteomeData$hgnc_symbol[which(proteomeData$asdPPI == "bait")]
  hek293T_prey = proteomeData$hgnc_symbol[which(proteomeData$asdPPI == "prey")]
  hek293T_baitprey = c(hek293T_bait, hek293T_prey)
  hek293T_other = proteomeData$hgnc_symbol[which(proteomeData$asdPPI == "other")]
  notHEK293T = inputGenes[!(inputGenes %in% proteomeData$hgnc_symbol)]
  allHEK293T = inputGenes[(inputGenes %in% proteomeData$hgnc_symbol)]
  alLGtex = inputGenes
  
  # calculate median geneset rank
  medRank_bait = apply(exprData_ranked[which(rownames(exprData_ranked) %in% hek293T_bait),], 2, median)
  medRank_prey = apply(exprData_ranked[which(rownames(exprData_ranked) %in% hek293T_prey),], 2, median)
  medRank_baitprey= apply(exprData_ranked[which(rownames(exprData_ranked) %in% hek293T_baitprey),], 2, median)
  medRank_other= apply(exprData_ranked[which(rownames(exprData_ranked) %in% hek293T_other),], 2, median)
  medRank_notHEK293T = apply(exprData_ranked[which(rownames(exprData_ranked) %in% notHEK293T),], 2, median)
  medRank_allHEK293T = apply(exprData_ranked[which(rownames(exprData_ranked) %in% allHEK293T),], 2, median)
  medRank_allGtex = apply(exprData_ranked[which(rownames(exprData_ranked) %in% inputGenes),], 2, median)
  
  # calculate median geneset percentile expression
  medPercentile_bait = apply(exprData_percentile[which(rownames(exprData_percentile) %in% hek293T_bait),], 2, median)
  medPercentile_prey = apply(exprData_percentile[which(rownames(exprData_percentile) %in% hek293T_prey),], 2, median)
  medPercentile_baitprey= apply(exprData_percentile[which(rownames(exprData_percentile) %in% hek293T_baitprey),], 2, median)
  medPercentile_other= apply(exprData_percentile[which(rownames(exprData_percentile) %in% hek293T_other),], 2, median)
  medPercentile_notHEK293T = apply(exprData_percentile[which(rownames(exprData_percentile) %in% notHEK293T),], 2, median)
  medPercentile_allHEK293T = apply(exprData_percentile[which(rownames(exprData_percentile) %in% allHEK293T),], 2, median)
  medPercentile_allGtex = apply(exprData_percentile[which(rownames(exprData_percentile) %in% inputGenes),], 2, median)
  
  summary = tibble(samples = colnames(exprData),
                   tissue = inputData$samples$tissue,
                   brain = inputData$samples$brain,
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
                   medRank_allGtex = medRank_allGtex,
                   medPercentile_notHEK293T = medPercentile_notHEK293T,
                   medPercentile_allHEK293T = medPercentile_allHEK293T,
                   medPercentile_allGtex = medPercentile_allGtex)
  
  return(summary)
}


# ------------------------------------------------------------------------------
#  ANALYSIS 1- Assess for enrichment of 293T bait, prey, and other genes in gtex samples
# ------------------------------------------------------------------------------

# Calculate median percentile expression for 293T bait, prey, and baitprey geneset compared to 293T background ("other") genes in gtex tissue samples
gtex_stats = hek293tASDppiStats(inputData = gtex,
                                 proteomeData = hek293tProteome)

# Format data for plotting (pivot longer)
gtex_stats_toPlot = gtex_stats %>% pivot_longer(
  cols = medRank_bait:medRank_allGtex,
  names_to = c(".value", "geneset"),
  names_sep = "_")
gtex_stats_toPlot = gtex_stats_toPlot %>%
  mutate(geneset = factor(geneset, levels = c("bait", "prey", "baitprey", "other", "allHEK293T", "notHEK293T", "allGtex"))) 

# plot median geneset expression for bait, prey, and other in GTEx brain samples
### paired t-tests to compare
p1 =  gtex_stats_toPlot %>%  
  filter(geneset %in% c("bait", "prey", "other") & tissue == "Brain") %>%
  mutate(geneset = recode_factor(geneset, 'bait'='Bait', 'prey'='Prey (- hcASD102)', 'other'="Other")) %>%
  ggplot(aes(x = geneset, y = medPercentile*100, fill = geneset)) +
  geom_boxplot(show.legend=F) +
  #geom_jitter(size = 1, alpha = 0.3) +
  stat_compare_means(comparisons=list(c("Bait", "Prey (- hcASD102)"), c("Prey (- hcASD102)", "Other"), c("Bait", "Other")), method = "t.test", aes(label=..p.adj..)) +
  labs(title = "Median geneset expression (percentile)\nGTEx brain samples", y = "Median expression (percentile)", x = "Geneset") +
  theme_bw()
p1
ggsave(paste0(od, "asdPPI_gtexBrain_expression_boxplot.pdf"), width=5, height=6)

## check stats for supplementary table
# b = gtex_stats_toPlot %>% filter(geneset=="bait" & tissue == "Brain") %>% mutate(medPercentile = medPercentile*100) %>% pull(medPercentile)
# o = gtex_stats_toPlot %>% filter(geneset=="other" & tissue == "Brain") %>% mutate(medPercentile = medPercentile*100) %>% pull(medPercentile)
# p = gtex_stats_toPlot %>% filter(geneset=="prey" & tissue == "Brain") %>% mutate(medPercentile = medPercentile*100) %>% pull(medPercentile)
# t.test(b,p)
# t.test(b,o)
# t.test(p,o)


# ------------------------------------------------------------------------------
# ANALYSIS 2 - Calculate whether prey genes are expressed at higher levels than expected in GTEx samples using permutation testing
# ------------------------------------------------------------------------------

# We want to assess whether ASD-PPI prey genes are more highly expressed than expected in GTEx brain samples compared to 293T background genes, when matched for proteomic expression level
# Specifically, we want to calculate whether prey genes have higher than expected median geneset rank in GTEx samples compared to a null distribution of randomly permuted "other" genes. Each permuted "other" geneset must:
# 1) have same number of genes as prey genes
# 2) have median expr rank within 1 stDev that of prey genes


# FUNCTION to create permuted genesets from 293T proteome-expressed genes (excluding bait) that fit the following criteria:
# 1) have same number of genes as prey genes (that are expressed in 293T proteome AND measured in GTEx)
# 2) have 293T proteome median expr rank within 1 stDev that of prey genes. Proteome expression level is determiend by iBAQ score
# output is a list, where
# 1) permutedGenesets: list in which each entry is a vector containing gene names that fit criteria detemrined above
#2) permutedGenesets_medRank293TiBAQ: vector in which each entry is median 293T iBAQ rank of permuted geneset
# 3) preyInfo: vector with information on prey genes(median 293T iBAQ rank, accepted lower and upper range in selecting permuted genesets, numPrey)
createPermutedGenesets <- function(proteomeSummaryData = genesets,
                                   inputData = "gtex",
                                   proteomeQuantification = "iBAQ",
                                   nPermutedSets = 10000){
  #calculate range of acceptable median geneset ranks (median rank of "prey" geneset +/- 1 stDev)
  range = genesets[[paste0("293T_", inputData, "_prey")]] %>% dplyr::select(paste0(proteomeQuantification, "_ranked")) %>% pull(.) %>% quantile
  range = range[c(2, 4)]
  
  # calculate number of prey genes
  numPrey = dim(genesets[[paste0("293T_", inputData, "_prey")]])[1]
  
  # identify set of genes that are expressed in 293T proteome and measured in inputData (excluding bait genes); and associated proteomic quantification rank
  toPermute = genesets[[paste0("293T_", inputData, "_notBait")]] %>% dplyr::select(hgnc_symbol, paste0(proteomeQuantification, "_ranked"))
  numNotBait = dim(toPermute)[1]
  # create vector to use to weigh probability of gene selection by proteome rank
  probNotBait = toPermute[,2] %>% pull(.)
  
  permutedGenesets = list()
  permutedGenesets_medRank293T = list()
  count = 0
  while(count < nPermutedSets){
    genes = sample(1:numNotBait, numPrey, replace = FALSE, prob = probNotBait)
    while((median(toPermute[genes, 2] %>% pull(.)) < range[1]) |
          (median(toPermute[genes,2] %>% pull(.)) > range[2]) ){
      genes = sample(1:numNotBait, numPrey, replace = FALSE)
    }
    count = count + 1
    permutedGenesets[[count]] = toPermute$hgnc_symbol[genes]
    permutedGenesets_medRank293T[[count]] = median(toPermute[genes, 2] %>% pull(.))
  }
  preyInfo = c("293T_iBAQmedRank" = median(genesets[[paste0("293T_", inputData, "_prey")]]$iBAQ_ranked),
               "lower_range" = range[1],
               "upper_range" = range[2],
               "numPrey" = numPrey)
  summary = list("permutedGenesets" = permutedGenesets,
                 "permutedGenesets_medRank293TiBAQ" = unlist(permutedGenesets_medRank293T),
                 "preyInfo" = preyInfo)
  
  return(summary)
}


# Format 293T proteome dataset for plotting (annotate by bait, prey, other gene status; rank proteome expression level)
hek293tProteome_toPlot = hek293tProteome %>%
  mutate(geneset = hgnc_symbol) %>%
  mutate(geneset = replace(geneset, geneset %in% bait293T, "bait")) %>%
  mutate(geneset = replace(geneset, geneset %in% prey293T, "prey")) %>%
  mutate(geneset = replace(geneset, geneset %in% other293T, "other")) %>%
  mutate(geneset = factor(geneset, levels = c("bait", "prey", "other"))) %>%
  mutate(iBAQ_ranked = rank(iBAQ_avg, ties.method = "average")) %>%
  mutate(LFQ_ranked = rank(LFQ_avg, ties.method = "average"))

# Create list of containing subsets of information from hek293tProteome table (genes with data in 293T proteome and GTEx, and associated 293T proteome quantification)
genesets = list()
genesets[["293T_gtex_bait"]] = hek293tProteome_toPlot %>% filter(asdPPI == "bait", gtexStatus == "present") %>% dplyr::select(hgnc_symbol, iBAQ_ranked, LFQ_ranked)
genesets[["293T_gtex_prey"]] = hek293tProteome_toPlot %>% filter(asdPPI == "prey", gtexStatus == "present") %>% dplyr::select(hgnc_symbol, iBAQ_ranked, LFQ_ranked)
genesets[["293T_gtex_other"]] = hek293tProteome_toPlot %>% filter(asdPPI == "other", gtexStatus == "present") %>% dplyr::select(hgnc_symbol, iBAQ_ranked, LFQ_ranked)
genesets[["293T_gtex_notBait"]] = hek293tProteome_toPlot %>% filter(asdPPI != "bait", gtexStatus == "present") %>% dplyr::select(hgnc_symbol, iBAQ_ranked, LFQ_ranked)


## Create permuted genesets from notBait 293T proteome genes that match the size and median 293T proteome iBAQ rank of ASD-PPI prey genes
## if permuted genesets have been previously generated, skip this step and load directly 


if(!file.exists(paste0(od, "/gtex_asdPPIprey_permutationData.RData"))){
  RNGkind("L'Ecuyer-CMRG")
  set.seed(28)
  permuted_by293TiBAQ_gtex = createPermutedGenesets(proteomeSummaryData = genesets,
                                                    inputData = "gtex",
                                                    proteomeQuantification = "iBAQ",
                                                    nPermutedSets = 100000)
  
  
  ## calculate median gene expression rank of each permuted genesets in GTEx data
  permuted_by293TiBAQ_gtex_medianExprRank = mclapply(permuted_by293TiBAQ_gtex$permutedGenesets, 
                                                     function(x) apply(gtex$exprData_ranked[x,], 2, median), mc.cores = (detectCores()- 8))
  
  
  # SAVE DATA
  save(permuted_by293TiBAQ_gtex,
       permuted_by293TiBAQ_gtex_medianExprRank,
       file = paste0(od, "gtex_asdPPIprey_permutationData.RData"))
}else(
  load(file = paste0(od, "gtex_asdPPIprey_permutationData.RData"))
)


##calculate median gene expression rank of prey geneset in GTEx data
prey_gtex_medianExprRank = apply(gtex$exprData_ranked[genesets$`293T_gtex_prey`$hgnc_symbol, ], 2, median)

## reformat data into dataframe
### add prey/observed data to first row
permuted_gtex_medianExprRank =  matrix(unlist(permuted_by293TiBAQ_gtex_medianExprRank), nrow = length(permuted_by293TiBAQ_gtex_medianExprRank), byrow = TRUE)
permuted_gtex_medianExprRank = rbind(prey_gtex_medianExprRank, permuted_gtex_medianExprRank)
colnames(permuted_gtex_medianExprRank) = colnames(gtex$exprData_ranked)
rownames(permuted_gtex_medianExprRank) = c("observed.prey", paste0("permuted.", c(1:(nrow(permuted_gtex_medianExprRank)-1))))

## calculate significance for each sample
signif_gtex_iBAQ = apply(permuted_gtex_medianExprRank, 2, function(x) sum(x[-1] >= x[1]) / length(x))

# format data for plotting
gtex_iBAQ_Permuted_toPlot = gtex$samples %>% 
  mutate(medRank_prey = prey_gtex_medianExprRank,
         medRank_median.permuted.notBait = apply(permuted_gtex_medianExprRank[-1,], 2, median),
         pVal_permuted = signif_gtex_iBAQ) %>%
  pivot_longer(cols = medRank_prey:medRank_median.permuted.notBait,
               names_to = c(".value", "geneset"),
               names_sep = "_") %>%
  mutate(geneset = factor(geneset, levels = c("prey", "median.permuted.notBait")))

gtex_iBAQ_Permuted_medianRanks_toPlot = t(permuted_gtex_medianExprRank) %>% as.data.frame %>% cbind(gtex$samples, .)

permuted =  gtex_iBAQ_Permuted_medianRanks_toPlot %>% dplyr::select(starts_with("permuted")) %>% 
  apply(., 1, median)

toPlot = gtex_iBAQ_Permuted_toPlot %>% filter(geneset == "prey") %>%
  mutate(log10.pVal_permuted = -log10(pVal_permuted+0.00001)) %>%
  mutate(permuted.medRank = permuted,
         rankdiff = medRank - permuted.medRank) %>% 
  filter(brain == "brain") %>% mutate(sample = gsub("Brain...", "", sample)) %>%
  mutate(pval.text = ifelse(pVal_permuted<0.05, sample, "")) %>%
  mutate(pointColor = ifelse(pVal_permuted<0.05, "red", "black"))

## Make plots
p2 = ggplot(toPlot, aes(x = rankdiff, y = log10.pVal_permuted)) +
  geom_point(alpha = 0.3, size = 4) +
  geom_point(data = subset(toPlot, pVal_permuted<0.05), aes(x = rankdiff, y=log10.pVal_permuted), color = "red", size = 4) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dotted", size = 1) +
  theme_bw() +
  labs(title = "Median prey expression rank, observed versus permutations\nGTEx, 13 brain sites",
       x = "Median geneset rank, observed - median(permuted)",
       y = "-log10(permuted pvalue)") +
  geom_hline(yintercept = -log10(0.05/13), color = "orange", linetype = "dotted", size = 1) +
  geom_text(aes(y = -log10(0.05), x = 50, label = "p = 0.05"), colour = "red", vjust = 1.5, hjust = 0, size = 4) +
  geom_text(aes(y = -log10(0.05/13), x = 50, label = "p.adj = 0.05"), colour = "orange", vjust = 1.5, hjust = 0, size = 4) +
  geom_text_repel(aes(x = rankdiff, y = log10.pVal_permuted, label = pval.text), colour = "black", size = 5)

p2
ggsave(paste0(od, "gtexBrain_expression_asdPPIprey_vs_permuted.pdf"), width=6, height=6)
