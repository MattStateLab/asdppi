################################################################################
# Assess whether ASD-PPI prey and hcASD102 genes are more highly expressed in BrainSpan samples compared to permuted genesets from the 293T proteome.
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

#Import ASD102 genes
load(paste0(rd, "asd102/asd102.Rdata"))
asd102genes = asd102$hgnc; rm(asd102)

# ------------------------------------------------------------------------------
# Format and annotate data 
# ------------------------------------------------------------------------------

# Extract bait and prey genes (defining prey genes to be exclusive of bait genes)
bait <- unique(ppiFull$bait)
prey <- setdiff(unique(ppiFull$prey), bait)

# Reformat BrainSpan data by converting expression level to rank within each sample (where higher rank reflects higher relative gene expression)
bsRNAseq[["exprData_ranked"]] = apply(bsRNAseq$exprData, 2, function(x) rank(x, ties.method = "average"))

# ------------------------------------------------------------------------------
# DEFINE FUNCTIONS
# ------------------------------------------------------------------------------


# FUNCTION to create permuted genesets from 293T proteome-expressed genes (excluding bait) that fit the following criteria:
# 1) have same number of genes as prey genes (that are expressed in 293T proteome AND measured in brainspan)
# 2) have 293T proteome median expr rank within 1 quantile  of mean median expr rank of prey genes. Proteome expression level is determiend by iBAQ score
# output is a list, where
# 1) permutedGenesets: list in which each entry is a vector containing gene names that fit criteria detemrined above
#2) permutedGenesets_medRank293TiBAQ: vector in which each entry is median 293T iBAQ rank of permuted geneset
# 3) preyInfo: vector with information on prey genes(median 293T iBAQ rank, accepted lower and upper range in selecting permuted genesets, numPrey)
createPermutedGenesets <- function(proteomeSummaryData = genesets,
                                   bsData = "bsRNAseq",
                                   proteomeQuantification = "iBAQ",
                                   nPermutedSets = 10000){
  #calculate range of acceptable median geneset ranks (median rank of "prey" geneset +/- 1 quantile)
  range = genesets[[paste0("293T_", bsData, "_prey")]] %>% dplyr::select(paste0(proteomeQuantification, "_ranked")) %>% pull(.) %>% quantile
  range = range[c(2, 4)]
  
  # calculate number of prey genes
  numPrey = dim(genesets[[paste0("293T_", bsData, "_prey")]])[1]
  
  # identify set of genes that are expressed in 293T proteome and measured in bsData (excluding bait genes); and associated proteomic quantification rank
  toPermute = genesets[[paste0("293T_", bsData, "_notBait")]] %>% dplyr::select(hgnc_symbol, paste0(proteomeQuantification, "_ranked"))
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
  preyInfo = c("293T_iBAQmedRank" = median(genesets[[paste0("293T_", bsData, "_prey")]]$iBAQ_ranked),
               "lower_range" = range[1],
               "upper_range" = range[2],
               "numPrey" = numPrey)
  summary = list("permutedGenesets" = permutedGenesets,
                 "permutedGenesets_medRank293TiBAQ" = unlist(permutedGenesets_medRank293T),
                 "preyInfo" = preyInfo)
  
  return(summary)
}

# FUNCTION to create permuted genesets from 293T proteome-expressed genes (excluding bait) that fit the following criteria:
# 1) have same number of genes as asd102 genes (that are expressed in 293T proteome AND measured in brainspan)
# 2) have 293T proteome median expr rank within 1 quantile  of mean median expr rank of asd102 genes. Proteome expression level is determiend by iBAQ score
# output is a list, where
# 1) permutedGenesets: list in which each entry is a vector containing gene names that fit criteria detemrined above
#2) permutedGenesets_medRank293TiBAQ: vector in which each entry is median 293T iBAQ rank of permuted geneset
# 3) asd102Info: vector with information on asd102 genes(median 293T iBAQ rank, accepted lower and upper range in selecting permuted genesets, numasd102)
createPermutedGenesets.asd102 <- function(proteomeSummaryData = genesets,
                                   bsData = "bsRNAseq",
                                   proteomeQuantification = "iBAQ",
                                   nPermutedSets = 100){
  #calculate range of acceptable median geneset ranks (median rank of "asd102" geneset +/- 1 quantile)
  range = genesets[[paste0("293T_", bsData, "_asd102")]] %>% dplyr::select(paste0(proteomeQuantification, "_ranked")) %>% pull(.) %>% quantile
  range = range[c(2, 4)]
  
  # calculate number of asd102 genes
  numasd102 = dim(genesets[[paste0("293T_", bsData, "_asd102")]])[1]
  
  # identify set of genes that are expressed in 293T proteome and measured in bsData (excluding bait genes); and associated proteomic quantification rank
  toPermute = genesets[[paste0("293T_", bsData, "_all")]] %>% dplyr::select(hgnc_symbol, paste0(proteomeQuantification, "_ranked"))
  numRNAseq = dim(toPermute)[1]
  # create vector to use to weigh probability of gene selection by proteome rank
  probExpr = toPermute[,2] %>% pull(.)
  
  permutedGenesets = list()
  permutedGenesets_medRank293T = list()
  count = 0
  while(count < nPermutedSets){
    genes = sample(1:numRNAseq, numasd102, replace = FALSE, prob = probExpr)
    while((median(toPermute[genes, 2] %>% pull(.)) < range[1]) |
          (median(toPermute[genes,2] %>% pull(.)) > range[2]) ){
      genes = sample(1:numRNAseq, numasd102, replace = FALSE)
    }
    count = count + 1
    permutedGenesets[[count]] = toPermute$hgnc_symbol[genes]
    permutedGenesets_medRank293T[[count]] = median(toPermute[genes, 2] %>% pull(.))
  }
  asd102Info = c("293T_iBAQmedRank" = median(genesets[[paste0("293T_", bsData, "_asd102")]]$iBAQ_ranked),
                 "lower_range" = range[1],
                 "upper_range" = range[2],
                 "numasd102" = numasd102)
  summary = list("permutedGenesets" = permutedGenesets,
                 "permutedGenesets_medRank293TiBAQ" = unlist(permutedGenesets_medRank293T),
                 "asd102Info" = asd102Info)
  
  return(summary)
}
# ------------------------------------------------------------------------------
# PPI PREY PERMUTATION ANALYSIS: Calculate whether prey genes are expressed at higher levels than expected in BrainSpan samples using permutation testing
# ------------------------------------------------------------------------------
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

# Format 293T proteome
hek293tProteome_toPlot = hek293tProteome %>%
  mutate(geneset = hgnc_symbol) %>%
  mutate(geneset = replace(geneset, geneset %in% bait293T, "bait")) %>%
  mutate(geneset = replace(geneset, geneset %in% prey293T, "prey")) %>%
  mutate(geneset = replace(geneset, geneset %in% other293T, "other")) %>%
  mutate(geneset = factor(geneset, levels = c("bait", "prey", "other"))) %>%
  mutate(iBAQ_ranked = rank(iBAQ_avg, ties.method = "average")) %>%
  mutate(LFQ_ranked = rank(LFQ_avg, ties.method = "average"))

#Looking at the mean and median proteome rank of different genesets - it is clear that prey genes tend to have higher proteomic expression

# hek293tProteome_toPlot %>%
#   group_by(geneset) %>%
#   dplyr::summarise(mean_iBAQrank = mean(iBAQ_ranked),
#                    med_iBAQrank = median(iBAQ_ranked),
#                    mean_LFQrank = mean(LFQ_ranked),
#                    med_LFQrank = median(LFQ_ranked))
# geneset mean_iBAQrank med_iBAQrank mean_LFQrank med_LFQrank
# <fct>           <dbl>        <dbl>        <dbl>       <dbl>
#   1 bait            5232.        5350         6057.       6312.
# 2 prey            7247.        7885         7592.       8313 
# 3 other           5412.        5280.        5370.       5236.

# We want to assess whether ASD-PPI prey genes are more highly expressed than expected in BrainSpan fetal samples compared to 293T background genes, when matched for proteomic expression level
# Specifically, we want to calculate whether prey genes have higher than expected median geneset rank in BrainSpan samples compared to a null distribution of randomly permuted "other" genes. Each permuted "other" geneset must:
# 1) have same number of genes as prey genes
# 2) have median expr rank within 1 quantile of mean median expr rank of prey genes

# Create list of containing subsets of information from hek293tProteome table (genes with data in 293T proteome and BrainSpan, and associated 293T proteome quantification)
genesets = list()
genesets[["293T_bsRNAseq_bait"]] = hek293tProteome_toPlot %>% filter(asdPPI == "bait", bsRNAseqStatus == "present") %>% dplyr::select(hgnc_symbol, iBAQ_ranked, LFQ_ranked)
genesets[["293T_bsRNAseq_prey"]] = hek293tProteome_toPlot %>% filter(asdPPI == "prey", bsRNAseqStatus == "present") %>% dplyr::select(hgnc_symbol, iBAQ_ranked, LFQ_ranked)
genesets[["293T_bsRNAseq_other"]] = hek293tProteome_toPlot %>% filter(asdPPI == "other", bsRNAseqStatus == "present") %>% dplyr::select(hgnc_symbol, iBAQ_ranked, LFQ_ranked)
genesets[["293T_bsRNAseq_notBait"]] = hek293tProteome_toPlot %>% filter(asdPPI != "bait", bsRNAseqStatus == "present") %>% dplyr::select(hgnc_symbol, iBAQ_ranked, LFQ_ranked)

## Create permuted genesets from  HEK293T proteome genes (excluding ASD-PPI bait) that match the size and median 293T proteome iBAQ rank of ASD-PPI prey genes, and format for plotting
## if permuted genesets have been previously generated, skip this step adn load directly 
if(!file.exists(paste0(od, "/bsRNAseq_asdPPIprey_permutationData_formattedToPlot.RData"))){
  RNGkind("L'Ecuyer-CMRG")
  set.seed(28)
  permuted_by293TiBAQ_bsRNAseq = createPermutedGenesets(proteomeSummaryData = genesets,
                                                        bsData = "bsRNAseq",
                                                        proteomeQuantification = "iBAQ",
                                                        nPermutedSets = 100000)
  
  ## calculate median gene expression rank of each permuted genesets in brainspan data
  permuted_by293TiBAQ_bsRNAseq_medianExprRank = mclapply(permuted_by293TiBAQ_bsRNAseq$permutedGenesets, 
                                                         function(x) apply(bsRNAseq$exprData_ranked[x,], 2, median), mc.cores = (detectCores()- 8))
  
  
  ##calculate median gene expression rank of prey geneset in brainspan data
  prey_bsRNAseq_medianExprRank = apply(bsRNAseq$exprData_ranked[genesets$`293T_bsRNAseq_prey`$hgnc_symbol, ], 2, median)
  
  ## reformat data into dataframe
  ### add prey/observed data to first row
  permuted_bsRNAseq_medianExprRank =  matrix(unlist(permuted_by293TiBAQ_bsRNAseq_medianExprRank), nrow = length(permuted_by293TiBAQ_bsRNAseq_medianExprRank), byrow = TRUE)
  permuted_bsRNAseq_medianExprRank = rbind(prey_bsRNAseq_medianExprRank, permuted_bsRNAseq_medianExprRank)
  colnames(permuted_bsRNAseq_medianExprRank) = colnames(bsRNAseq$exprData_ranked)
  rownames(permuted_bsRNAseq_medianExprRank) = c("observed.prey", paste0("permuted.", c(1:(nrow(permuted_bsRNAseq_medianExprRank)-1))))
  
  ## calculate significance for each sample
  signif_bsRNAseq_iBAQ = apply(permuted_bsRNAseq_medianExprRank, 2, function(x) sum(x[-1] >= x[1]) / length(x))
  
  # format data for plotting
  bsRNAseq_iBAQ_Permuted_toPlot = bsRNAseq$samples %>% 
    mutate(medRank_prey = prey_bsRNAseq_medianExprRank,
           medRank_median.permuted.notBait = apply(permuted_bsRNAseq_medianExprRank[-1,], 2, median),
           pVal_permuted = signif_bsRNAseq_iBAQ) %>%
    pivot_longer(cols = medRank_prey:medRank_median.permuted.notBait,
                 names_to = c(".value", "geneset"),
                 names_sep = "_") %>%
    mutate(geneset = factor(geneset, levels = c("prey", "median.permuted.notBait")),
           natal = factor(natal, levels = c("prenatal", "postnatal")))
  
  ## for plotting median of medExprRank of observed prey vs permuted sets
  bsRNAseq_iBAQ_Permuted_medianRanks_toPlot = t(permuted_bsRNAseq_medianExprRank) %>% as.data.frame %>% cbind(bsRNAseq$samples, .)
  
  # SAVE DATA
  save(permuted_by293TiBAQ_bsRNAseq,
       permuted_by293TiBAQ_bsRNAseq_medianExprRank,
       bsRNAseq_iBAQ_Permuted_toPlot, 
       bsRNAseq_iBAQ_Permuted_medianRanks_toPlot, 
       file = paste0(od, "bsRNAseq_asdPPIprey_permutationData_formattedToPlot.RData"))
}else{
  load(file = paste0(od, "bsRNAseq_asdPPIprey_permutationData_formattedToPlot.RData"))
}


# ------------------------------------------------------------------------------
# hcASD102 PERMUTATION ANALYSIS: Calculate whether hcASD102 genes are expressed at higher levels than expected in BrainSpan samples using permutation testing
# ------------------------------------------------------------------------------

# Annotate 293T proteome by:
# 1) ASD102 status
# 2) in BrainSpan RNAseq data
hek293tProteome <- hek293tProteome_o %>%
  mutate(asd102 = hgnc_symbol) %>%
  mutate(asd102 = replace(asd102, asd102 %in% asd102genes, "asd102")) %>%
  mutate(asd102 = replace(asd102, !asd102 == "asd102", "other")) %>%
  mutate(bsRNAseqStatus = hgnc_symbol) %>%
  mutate(bsRNAseqStatus = ifelse(bsRNAseqStatus %in% bsRNAseq[["genes"]][, "hgnc_symbol"], "present", "absent")) %>%
  filter(hgnc_symbol != "") # remove genes that do not have hugo gene name (ex:ENSG00000187186)

asd102.293T = hek293tProteome %>% filter(asd102 == "asd102") %>% .$hgnc_symbol
other.293T = hek293tProteome %>% filter(asd102 == "other") %>% .$hgnc_symbol

# Format 293T proteome 
hek293tProteome_toPlot = hek293tProteome %>%
  mutate(geneset = hgnc_symbol) %>%
  mutate(geneset = replace(geneset, geneset %in% asd102.293T, "asd102")) %>%
  mutate(geneset = replace(geneset, geneset %in% other.293T, "other")) %>%
  mutate(geneset = factor(geneset, levels = c("asd102", "other"))) %>%
  mutate(iBAQ_ranked = rank(iBAQ_avg, ties.method = "average")) %>%
  mutate(LFQ_ranked = rank(LFQ_avg, ties.method = "average"))



# We want to assess whether ASD102 genes are more highly expressed than expected in BrainSpan fetal samples compared to 293T background genes, when matched for proteomic expression level
# Specifically, we want to calculate whether prey genes have higher than expected median geneset rank in BrainSpan samples compared to a null distribution of randomly permuted "other" genes. Each permuted "other" geneset must:
# 1) have same number of genes as asd102 genes
# 2) have median expr rank within 1 quantile of mean median expr rank of asd102 genes

# Create list of containing subsets of information from hek293tProteome table (genes with data in 293T proteome and BrainSpan, and associated 293T proteome quantification)
genesets = list()
genesets[["293T_bsRNAseq_asd102"]] = hek293tProteome_toPlot %>% filter(asd102 == "asd102", bsRNAseqStatus == "present") %>% dplyr::select(hgnc_symbol, iBAQ_ranked, LFQ_ranked)
genesets[["293T_bsRNAseq_other"]] = hek293tProteome_toPlot %>% filter(asd102 == "other", bsRNAseqStatus == "present") %>% dplyr::select(hgnc_symbol, iBAQ_ranked, LFQ_ranked)
genesets[["293T_bsRNAseq_all"]] = hek293tProteome_toPlot %>% filter(bsRNAseqStatus == "present") %>% dplyr::select(hgnc_symbol, iBAQ_ranked, LFQ_ranked)


## Create permuted genesets from  HEK293T proteome genes (excluding ASD-PPI bait) that match the size and median 293T proteome iBAQ rank of hcASD102 genes, and format for plotting
## if permuted genesets have been previously generated, skip this step adn load directly 
if(!file.exists(paste0(od, "/bsRNAseq_asd102_permutationData_formattedToPlot.RData"))){
  RNGkind("L'Ecuyer-CMRG")
  set.seed(28)
  permuted_by293TiBAQ_bsRNAseq = createPermutedGenesets.asd102(proteomeSummaryData = genesets,
                                                               bsData = "bsRNAseq",
                                                               proteomeQuantification = "iBAQ",
                                                               nPermutedSets = 100000)
  
  
  ## calculate median gene expression rank of each permuted genesets in brainspan data
  permuted_by293TiBAQ_bsRNAseq_medianExprRank = mclapply(permuted_by293TiBAQ_bsRNAseq$permutedGenesets, 
                                                         function(x) apply(bsRNAseq$exprData_ranked[x,], 2, median), mc.cores = (detectCores()- 6))
  
  
  ##calculate median gene expression rank of prey geneset in brainspan data
  asd102_bsRNAseq_medianExprRank = apply(bsRNAseq$exprData_ranked[genesets$`293T_bsRNAseq_asd102`$hgnc_symbol, ], 2, median)
  
  ## reformat data into dataframe
  ### add asd102/observed data to first row
  permuted_bsRNAseq_medianExprRank =  matrix(unlist(permuted_by293TiBAQ_bsRNAseq_medianExprRank), nrow = length(permuted_by293TiBAQ_bsRNAseq_medianExprRank), byrow = TRUE)
  permuted_bsRNAseq_medianExprRank = rbind(asd102_bsRNAseq_medianExprRank, permuted_bsRNAseq_medianExprRank)
  colnames(permuted_bsRNAseq_medianExprRank) = colnames(bsRNAseq$exprData_ranked)
  rownames(permuted_bsRNAseq_medianExprRank) = c("observed.asd102", paste0("permuted.", c(1:(nrow(permuted_bsRNAseq_medianExprRank)-1))))
  
  ## calculate significance for each sample
  signif_bsRNAseq_iBAQ = apply(permuted_bsRNAseq_medianExprRank, 2, function(x) sum(x[-1] >= x[1]) / length(x))
  
  # format data for plotting
  bsRNAseq_iBAQ_Permuted_toPlot = bsRNAseq$samples %>% 
    mutate(medRank_asd102 = asd102_bsRNAseq_medianExprRank,
           medRank_median.permuted.all = apply(permuted_bsRNAseq_medianExprRank[-1,], 2, median),
           pVal_permuted = signif_bsRNAseq_iBAQ) %>%
    pivot_longer(cols = medRank_asd102:medRank_median.permuted.all,
                 names_to = c(".value", "geneset"),
                 names_sep = "_") %>%
    mutate(geneset = factor(geneset, levels = c("asd102", "median.permuted.all")),
           natal = factor(natal, levels = c("prenatal", "postnatal")))
  
  ## for plotting median of medExprRank of observed prey vs permuted sets
  bsRNAseq_iBAQ_Permuted_medianRanks_toPlot = t(permuted_bsRNAseq_medianExprRank) %>% as.data.frame %>% cbind(bsRNAseq$samples, .)
  
  # SAVE DATA
  save(permuted_by293TiBAQ_bsRNAseq,
       permuted_by293TiBAQ_bsRNAseq_medianExprRank,
       bsRNAseq_iBAQ_Permuted_toPlot, 
       bsRNAseq_iBAQ_Permuted_medianRanks_toPlot, 
       file = paste0(od, "bsRNAseq_asd102_permutationData_formattedToPlot.RData"))
}else{
  load(file = paste0(od, "bsRNAseq_asd102_permutationData_formattedToPlot.RData"))
}

