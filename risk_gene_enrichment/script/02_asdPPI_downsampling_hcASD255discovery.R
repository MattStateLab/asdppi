# Goal: downsampling analysis to assess whether enlarging ASD-PPI network increases ability to detect new hcASD genes (Fu et al 2022 TADA FDR <0.1)
# Authors: Belinda Wang

################################################################################
# setup workspace
################################################################################
# [ ] replace genes metadata with supplementary table?

library(tidyverse)
library(readxl)
library(limma)

# Run from the repository root.
wd <- "risk_gene_enrichment/"
dd <- "risk_gene_enrichment/data/"
od <- "risk_gene_enrichment/output/"
rd <- "data/"

# Import PPI data 
load(paste0(rd, "ppi100.RData"))

# importa gees metadata
# - hgnc_symbol: HGNC gene symbol  
# - gene_gencodeV33: 
# - hek293Texpr: indicator for whether gene is one of the n=11167 proteins expressed in HEK293T cells (union of Bekker Jensen 2017 and prey100; 1=T; 0=F)  
# - bait100: indicator for whether gene is one of the n=100 hcASD genes used as bait in ASD-PPI_100 (1=T; 0=F)  
# - prey100: indicator for whether gene is one of the n=1074 ASD-PPI_100 prey (1=T; 0=F). Note: bait that were identified as prey are included here.
# - asd102_satterstrom2020: indicator for whether gene is one of the n=102 hcASD genes identified in Satterstrom et al 2020 (1=T; 0=F)  
# - asd255_fu2022: indicator for whether gene is one of the n=255 hcASD genes identified in Fu et al 2022 (1=T; 0=F)  
# - asd72_zhou2022: indicator for whether gene is one of the n=72 genes with study-wide significance (based on 5,754 constraint genes, p<8.69E-06) in  analysis of SPARK data in Zhou et al Nature 2022 (PMID 35982159)
# - dd285_kaplanis2020:  indicator for whether gene is one of the n=285 unique genes that are significantly associated with DD (after one-sided Bonferroni correction) in Kaplanis et al Nature 2020 (PMID 33057194) (1=T; 0=F)  
# - scz34_singh2022: indicator for whether gene is one of the n=34 genes with FDR<0.1 from SCHEMA WES  Singh et al Nature 2022 (PMID 35396579) (1=T; 0=F)  
# - sfari_all: indicator for whether gene is one of the n=1020 genes in SFARI gene list (gene.sfari.org, 9/2021 release) (1=T; 0=F)  
# - sfari_s: indicator for whether gene is one of the n=230 genes in SFARI listed as syndromic (gene.sfari.org, 9/2021 release) (1=T; 0=F)    
# - sfari_tier1: indicator for whether gene is one of the n=206 genes in SFARI tier 1 (highest confidence) (gene.sfari.org, 9/2021 release) (1=T; 0=F)   
# - sfari_tier1s: indicator for whether gene is one of the n=92 genes in SFARI tier 1 (highest confidence) and syndromic (gene.sfari.org, 9/2021 release) (1=T; 0=F)  
# - sfari_tier2: indicator for whether gene is one of the n=219 genes in SFARI tier 2 (medium confidence) (gene.sfari.org, 9/2021 release)  
# - sfari_tier3: indicator for whether gene is one of the n=514 genes in SFARI tier 3 (lowest confidence)  (gene.sfari.org, 9/2021 release)  
# - satterstrom2020gene: indicator for whether gene is one of the n=17332 genes assessed in Satterstrom et al 2020 (1=T; 0=F)  
# - fu2022gene: indicator for whether gene is one of the n=17821 genes assessed in Fu et al 2022 (1=T; 0=F)  
# - zhou2022gene: indicator for whether gene is one of the n=5746 highly constrained genes considered in Zhou et al 2022 (1=T, 0=F)
# - kaplanis2020gene: indicator for whether gene is one of the n=19051 genes assessed in Kaplanis et al 2020 (1=T; 0=F) 
# - singh2022gene: indicator for whether gene is one of the n=18098 genes assessed in Singh et al 2022 (1=T; 0=F)  
# - fuSatterstrom: indicator for whether gene is one of the n=16982 genes assessed in Fu et al 2022 and Satterstrom et al 2020 (1=T; 0=F)  
# - zhouSatterstrom: indicator for whether gene is one of the n=5219 genes assessed in Zhou et al 2022 and Satterstrom et al 2020 (1=T; 0=F)  
# - kaplanisSatterstrom: indicator for whether gene is one of the n=17304 genes assessed in Kaplanis et al 2020 and Satterstrom et al 2020 (1=T; 0=F)  
# - singhSatterstrom:  indicator for whether gene is one of the n=16897 genes assessed in Singh et al 2022 and Satterstrom et al 2020 (1=T; 0=F)  

genesMetadata = read.csv(paste0(dd, "geneAnnotation_hek293T_ASD_DD_SCZ.csv"))



################################################################################
# DEFINE FUNCTIONS
################################################################################

# create function to subset ppi to a smaller set of bait/prey for a given number of iterations, and return a list with the given bait and prey
bw.subsetPPI <- function(f_ppi, nBait = 5, nIter = 10){
  f_potentialBait = unique(f_ppi$bait)
  f_subsetBait = lapply(1:nIter, function(x) sample(f_potentialBait, nBait))
  f_subsetPrey = lapply(f_subsetBait, function(x) f_ppi %>% filter(bait %in% x) %>% pull(prey) %>% unique)
  names(f_subsetBait) = paste0(1:length(f_subsetBait))
  names(f_subsetPrey) = paste0(1:length(f_subsetPrey))
  result = list("subsetBait" = f_subsetBait,
                "subsetPrey" = f_subsetPrey,
                "setInfo" = c("nBait" = nBait, "nIter" = nIter))
  return(result)
}


# function that takes list of bait and prey, and calculates OR and p value for prey and number of hcASD255 genes among prey
### Contingency table (gene counts): hcASD255 gene, not hcASD255 gene, prey, other HEK293T gene
### gebe universe:  1) in HEK293T proteome (Bekker Jensen 2017 + ASD-PPI prey); 2) assessed in Satterstrom 2020; and 3) assessed in Fu et al 2022
bw.subsetPreyOR_hcASD255 <- function(f_preyVector, f_genesMetadata){
  ### update genesMetadata table
  f_genesMetadata = f_genesMetadata %>% mutate(prey = ifelse(hgnc_symbol %in% f_preyVector, 1, 0))
  
  # make contingency table to perform fisher's exact test
  dat = data.frame(c(f_genesMetadata %>% filter(hek293Texpr == 1 & fuSatterstrom == 1 & prey == 1 & asd255_fu2022 == 1) %>% nrow,
                     f_genesMetadata %>% filter(hek293Texpr == 1 & fuSatterstrom == 1 & prey == 0 & asd255_fu2022 == 1) %>% nrow),
                   c(f_genesMetadata %>% filter(hek293Texpr == 1 & fuSatterstrom == 1 & prey == 1 & asd255_fu2022 == 0) %>% nrow,
                     f_genesMetadata %>% filter(hek293Texpr == 1 & fuSatterstrom == 1 & prey == 0 & asd255_fu2022 == 0) %>% nrow),
                   row.names = c("Prey", "293Tother"))
  colnames(dat) =  c("hcASD255", "not.hcASD255")
  
  #perform fisher's exact test (pvalue from on one-sided, greater. Also do one-sided, less to grab upper CI for plotting)
  stats1 = fisher.test(dat, alternative = "greater")
  stats2 = fisher.test(dat, alternative = "less")
  #save data
  results = c(stats1$p.value, stats1$conf.int[1], stats2$conf.int[2], stats1$estimate, unlist(dat))
  names(results) = c("pval", "CIlower", "CIupper", "OR", "a", "c", "b", "d")
  return(results)
}

# wrapper function that creates subsets of prey for a given number of bait, and returns pval/OR for prey and exome-ASD102-prey
bw.subsetPreyORList <- function(f_ppi, nBait=5, nIter=10, f_genesMetadata){
  f_subset = bw.subsetPPI(f_ppi, nBait, nIter)
  f_subsetPrey = f_subset$subsetPrey
  temp1 = lapply(f_subsetPrey, function(x) bw.subsetPreyOR_hcASD255(x, f_genesMetadata))
  #annotate with iteration number, number of bait, number of prey
  temp2 = do.call(rbind, temp1) %>% as.data.frame %>% mutate(nBait=nBait, 
                                                             iteration = names(temp1),
                                                             numASD255 = a)
  return(temp2)
}

################################################################################
# assess relationship between number of bait and ability to identify prey that are hcASD255
################################################################################
#To determine whether increasing the number of baits used to create ASD-PPI is associated with increased ability to identify prey associated with ASD genetic risk, we will downsample the ASD-PPI network. Specifically, we randomly selected sets of bait (ranging from size 1-100 genes, 1000 iterations for each set size), and trim the ASD-PPI network to include only prey associated with the downsampled bait. For each downsampled network, we will calculate whether associated prey are enriched for Fu et al hcASD255 genes.  For each bait set size, we calculated the median and standard deviation of ORs across the 1000 iterations.

#DOWNSAMPLING ANALYSIS, INCLUDING hcASD102
###Downsample the number of baits (randomly selecting sets of size 1 to 100, 1000 iterations for esach set size). For each PPI subset, calculate OR for enrichment of hcASD255 genes in prey
if(!file.exists(paste0(od, "downsampledASDppi100_preyOR_hcASD255.Rdata"))){
  set.seed(28)
  subsetSummary = lapply(1:100, function(x) bw.subsetPreyORList(ppi100, nBait = x, nIter = 1000, genesMetadata))
  # save data
  save(subsetSummary, file=paste0(od, "downsampledASDppi100_preyOR_hcASD255.Rdata"))
}
load(file=paste0(od, "downsampledASDppi100_preyOR_hcASD255.Rdata"))

#Format data for plotting
summary = do.call(rbind, subsetSummary)
toPlot = summary %>% dplyr::group_by(nBait) %>% 
  dplyr::summarize(mean.pval = mean(pval, na.rm = T),
                   median.pval = median(pval, na.rm = T),
                   median.CIlower = median(CIlower[which(CIlower != 0)]), #I think we need to remove CIlower = 0 samples
                   median.CIupper = median(CIupper[which(CIupper !=Inf)]), # I think we need to remove CIupper = Inf samples
                   nIter = n(),
                   num = sum(!is.na(OR)),
                   sd.pval = sd(pval, na.rm = T),
                   se.pval = sd.pval/sqrt(num),
                   mean.OR = mean(OR, na.rm = T),
                   median.OR = median(OR, na.rm = T),
                   min.OR = min(OR, na.rm = T),
                   sd.OR = sd(OR, na.rm = T),
                   se.OR = sd.OR/sqrt(num),
                   median.numASD255 = median(numASD255),
                   sd.numASD255 = sd(numASD255),
                   se.numASD255 = sd(numASD255)) %>%
  mutate(SElower.OR = median.OR - sd.OR,
         SEupper.OR = median.OR + sd.OR)


#DOWNSAMPLING ANALYSIS, EXCLUDING hcASD102
### downsample the number of baits (randomly selecting sets of size 1 to 100, 1000 iterations for esach set size). For each PPI subset, calculate OR for prey (excluding ASD102) and exome-prey-asd102
genesMetadata_noASD102 = genesMetadata %>% filter(asd102_satterstrom2020==0)
if(!file.exists(paste0(od, "downsampledASDppi100_preyOR_hcASD255_nothcASD102.Rdata"))){
  set.seed(28)
  subsetSummary.noASD102 = lapply(1:100, function(x) bw.subsetPreyORList(ppi100, nBait = x, nIter = 1000, genesMetadata_noASD102))
  # save data
  save(subsetSummary.noASD102, file=paste0(od, "downsampledASDppi100_preyOR_hcASD255_nothcASD102.Rdata"))
}
load(file=paste0(od, "downsampledASDppi100_preyOR_hcASD255_nothcASD102.Rdata"))

#Format data for plotting
summary.noASD102 = do.call(rbind, subsetSummary.noASD102)


################################################################################
# MAKE FIGURE
################################################################################

# format hcASD255 downsampling results for plotting
temp = summary %>% dplyr::group_by(nBait) %>% 
  dplyr::summarize(mean.pval = mean(pval, na.rm = T),
                   median.pval = median(pval, na.rm = T),
                   median.CIlower = median(CIlower[which(CIlower != 0)]), #I think we need to remove CIlower = 0 samples
                   median.CIupper = median(CIupper[which(CIupper !=Inf)]), # I think we need to remove CIupper = Inf samples
                   nIter = n(),
                   num = sum(!is.na(OR)),
                   sd.pval = sd(pval, na.rm = T),
                   se.pval = sd.pval/sqrt(num),
                   mean.OR = mean(OR, na.rm = T),
                   median.OR = median(OR, na.rm = T),
                   min.OR = min(OR, na.rm = T),
                   sd.OR = sd(OR, na.rm = T),
                   se.OR = sd.OR/sqrt(num),
                   median.numASD255 = median(numASD255),
                   sd.numASD255 = sd(numASD255),
                   se.numASD255 = sd(numASD255)) %>%
  mutate(SElower.OR = median.OR - sd.OR,
         SEupper.OR = median.OR + sd.OR) %>% 
  mutate(hcASD_geneset="hcASD255")

# format hcASD255 (-hcASD102) downsampling results for plotting
temp.noASD102 = summary.noASD102 %>% dplyr::group_by(nBait) %>% 
  dplyr::summarize(mean.pval = mean(pval, na.rm = T),
                   median.pval = median(pval, na.rm = T),
                   median.CIlower = median(CIlower[which(CIlower != 0)]), #I think we need to remove CIlower = 0 samples
                   median.CIupper = median(CIupper[which(CIupper !=Inf)]), # I think we need to remove CIupper = Inf samples
                   nIter = n(),
                   num = sum(!is.na(OR)),
                   sd.pval = sd(pval, na.rm = T),
                   se.pval = sd.pval/sqrt(num),
                   mean.OR = mean(OR, na.rm = T),
                   median.OR = median(OR, na.rm = T),
                   min.OR = min(OR, na.rm = T),
                   sd.OR = sd(OR, na.rm = T),
                   se.OR = sd.OR/sqrt(num),
                   median.numASD255 = median(numASD255),
                   sd.numASD255 = sd(numASD255),
                   se.numASD255 = sd(numASD255)) %>%
  mutate(SElower.OR = median.OR - sd.OR,
         SEupper.OR = median.OR + sd.OR)%>% 
  mutate(hcASD_geneset = "hcASD255 (-hcASD102)")


# combine ddata from downsampling analysis for hcASD255 and hcASD255 (-hcASD102)
toPlot = rbind(temp, temp.noASD102)

# make plot
### calculate threshold for number of ASD-PPI bait for FET to have median p value that i
sig.nBait = temp %>% filter(median.pval<0.05) %>% pull(nBait) %>% min
sig.nBait.noASD102 = temp.noASD102 %>% filter(median.pval<0.05) %>% pull(nBait) %>% min
ggplot(toPlot, aes(x=nBait, y=median.numASD255, color=hcASD_geneset, fill = hcASD_geneset)) + geom_line() +
  geom_ribbon(aes(ymin=median.numASD255-sd.numASD255, ymax = median.numASD255+sd.numASD255), alpha=0.1, linetype=0) +
  geom_vline(xintercept = sig.nBait, linetype = "dashed", color = "#F8766D") +
  geom_vline(xintercept = sig.nBait.noASD102, linetype = "dashed", color = "#00BFC4") +
  labs(x="Number of bait", y="Median(number of hcASD255 among prey)",
       title="Effect of enlarging ASD-PPI dataset on identifying hcASD genes",
       caption = paste0("nBait for which hcASD255 median.pval<0.05: ", sig.nBait, "\nnBait for which hcASD255 (-hcASD102) median.pval<0.05: ", sig.nBait.noASD102)) +
  theme_bw() + theme(legend.position='bottom')
#ggsave(paste0(od, "asdPPI_downsampling_hcASD255_FET.pdf"), width=6, height=6)
ggsave(paste0(od, "asdPPI_downsampling_hcASD255_FET.pdf"), width=5, height=6)
