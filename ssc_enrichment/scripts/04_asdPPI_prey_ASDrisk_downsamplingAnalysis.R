################################################################################
# Downsampling analysis to assess relationship between number of bait and OR of prey enrichment for ASD de novo genetic risk
###############################################################################

# ------------------------------------------------------------------------------
#set up workspace 
# ------------------------------------------------------------------------------
# Run from the repository root.
wd <- "ssc_enrichment/"
dd <- "ssc_enrichment/data/"
od <- "ssc_enrichment/output/"
rd <- "data/"

#rm(list=ls())
library(tidyverse)
library(readxl)
library(limma)
library(ggplot2)
library(ggpubr)
library(data.table)

#IMPORT DATA
#### Import asdPPI_100
load(paste0(rd, "ppi100.RData"))
ppi = ppi100
ppi_bait <- unique(ppi$bait)
ppi_prey <- setdiff(unique(ppi$prey), ppi_bait)

#### Import gene annotations generated in 01_format_satterstrom2020_SSC_variants.R
# Definitions:  
# 1. damagingMissense: Polyphen Mis3 (damaging) or MPC MisB (MPC >=2)  
# 2. PTV: frameshift, stop gained, or canonical splice site disruption  
# 3. damaging: PTV + damagingMissense  
# 4. missense: missense coding mutations  
# 5. synonymous: synonymous coding mutations  

### geneAnnotations column guide (Satterstrom2020_geneAnnotations_SSC_caseControl_asdPPI100.csv)
# ennsembl_gene_id
# hgnc_symbol: HGNC gene symbol
# - ASD102: indicator for whether gene is hcASD gene (Satterstrom 2020); 1=yes, 2=no  
# - hek293TproteinExpr: : indicator for whether protein is expressed in HEK293T cells (Bekker-Jensen 2017 or prey); 1=yes, 2=no  
# - associationGene: indicator for whether gene is one of 17,484 autosomal genes assessed for genetic risk in Satterstrom 2020 Table S2; 1=yes, 2=no
# - bait: indicator for whether gene is ASD-PPI bait gene; 1=yes, 2=no  
# - prey: indicator for whether gene is ASD-PPI prey gene (excludes ASD-PPI bait); 1=yes, 2=no 
# - preyAll: indicator for whether gene is ASD-PPI prey gene (allows for ASD-PPI bait that were identified as prey); 1=yes, 2=no 
# - case.ptv: indicator for whether gene had PTV  identified in case (ASD) individual of DBS and Swedish case samples (Satterstrom 2020 Table S2); 1=yes, 2=no
# - control.ptv: number of de novo PTV variants identified in control individual of DBS and Swedish case samples (Satterstrom 2020 Table S2)
# - proband.dnPTV: number of de novo PTV variants identified in proband/ASD individual of famly study (Satterstrom 2020 Table S2)
# - proband.dnMisA: number of de novo MisA variants identified in proband/ASD individual of famly study (Satterstrom 2020 Table S2)
# - proband.dnMisB: number of de novo MisB variants identified in proband/ASD individual of famly study (Satterstrom 2020 Table S2)
# *  proband.dnm.damagingMissense = indicator for whether gene had damaging missense de novo mutation identified in ASD proband (Satterstrom 2020 Table S1); 1=yes, 2=no  
# *  sibling.dnm.damagingMissense = indicator for whether gene had damaging missense de novo mutation identified in unaffected sibling (Satterstrom 2020 Table S1); 1=yes, 2=no  
# *  proband.dnm.PTV = indicator for whether gene had de novo PTV identified in ASD proband (Satterstrom 2020 Table S1); 1=yes, 2=no  
# *  sibling.dnm.PTV = indicator for whether gene had de novo PTV identified in unaffected sibling (Satterstrom 2020 Table S1); 1=yes, 2=no  
# *  proband.dnm.damaging = indicator for whether gene had damaging  de novo mutation (PTV or damaging missense) identified in ASD proband (Satterstrom 2020 Table S1); 1=yes, 2=no  
# *  sibling.dnm.damaging = indicator for whether gene had damaging  de novo mutation (PTV or damaging missense) identified in unaffected sibling (Satterstrom 2020 Table S1); 1=yes, 2=no  
# *  proband.dnm.missense = indicator for whether gene had de novo missense mutation identified in ASD proband (Satterstrom 2020 Table S1); 1=yes, 2=no  
# *  sibling.dnm.missense = indicator for whether gene had de novo misssense mutation identified in unaffected sibling (Satterstrom 2020 Table S1); 1=yes, 2=no  
# *  proband.dnm.synonymous = indicator for whether gene had de novo synonymous mutation identified in ASD proband (Satterstrom 2020 Table S1); 1=yes, 2=no  
# *  sibling.dnm.synonymous = indicator for whether gene had de novo synonymous mutation identified in unaffected sibling (Satterstrom 2020 Table S1); 1=yes, 2=no  
# *  numProband.damaging = number of probands in SSC subset with a damaging variant in gene
# *  numSibling.damaging = number of siblings in SSC subset with a damaging variant in gene
geneAnnotations = read.csv(file=paste0(od, "Satterstrom2020_geneAnnotations_SSC_caseControl_asdPPI100.csv"))

### individualAnnotations column guide (Satterstrom2020_SSC_DNM.csv)
###individual de novo mutatiosn are annotated by functional class and bait/prey status to enable downstream individual-level analysis
# - hgnc_symbol: HGNC gene symbol
# - Child_ID: SSC child ID
# - Affected Status: indicator for whether child has ASD or not; 1=no, 2=yes
# - hek293TproteinExpr: : indicator for whether protein is expressed in HEK293T cells (Bekker-Jensen 2017 or prey); 1=yes, 2=no  
# - bait: indicator for whether gene is ASD-PPI bait gene; 1=yes, 0=no  
# - prey: indicator for whether gene is ASD-PPI prey gene (excludes ASD-PPI bait); 1=yes, 0=no 
# - preyAll: indicator for whether gene is ASD-PPI prey gene (allows for ASD-PPI bait that were identified as prey); 1=yes, 0=no 
# - ASD102: indicator for whether gene is hcASD gene (Satterstrom 2020); 1=yes, 0=no  
# - dnm.damagingMissense: indicator for whether child has damaging mutation in this gene; 1=yes, 0=no
# - dnm.PTV: indicator for whether child has a PTV in this gene; 1=yes, 0=no
# - dnm.damaging: indicator for whether child has a damaging variant (damaginMissense or PTV) in this gene; 1=yes, 0=no
# - dnm.missense: indicator for whether child has a missense variant in this gene; 1=yes, 0=no
# - dnm.synonymous: indicator for whether child has a synonymous variant in this gene; 1=yes, 0=no

individualAnnotations = read.csv(file=paste0(od, "Satterstrom2020_SSC_DNM.csv"))

#Load results from FET (summarDF_d, generated in 03_asdPPI_geneticRisk_Satterstrom2020_SSC.R) 
###FET: damaging mutation in geneset, no damaging mutation in genesent, ASD proband, ASD sibling.
### genesets: prey, all exome genes, exomeNotASD102, exomeNotASDPPI, HEK293Tother
summaryDF_d = read.csv(file=paste0(od, "FET_SSCsubset_results.csv"))



# ------------------------------------------------------------------------------
#DEFINE FUNCTIONS
# ------------------------------------------------------------------------------
#To determine whether increasing the number of baits used to create ASD-PPI is associated with increased ability to identify prey associated with ASD genetic risk, we will downsample the ASD-PPI network. Specifically, we randomly selected sets of bait (ranging from size 1-100 genes, 1000 iterations for each set size), and trim the ASD-PPI network to include only prey associated with the downsampled bait. For each downsampled network, we will calculate the association between having a damaging variant in different genesets of interest with ASD status. We assed 2 genesets: 1) downsampled ASD-PPI prey (excluding 102 hcASD risk genes); and 2) the human exome minus downsampleed ASD-PPI prey and hcASD genes. For each bait set size, we calculated the median and standard deviation of ORs across the 1000 iterations.


# DEFINE FUNCTIONS
# function to subset ppi to a smaller set of bait/prey for a given number of iterations, and return a list with the given bait and prey
bw.subsetPPI <- function(ppi_f, nBait = 5, nIter = 10){
  potentialBait_f = unique(ppi_f$bait)
  subsetBait_f = lapply(1:nIter, function(x) sample(potentialBait_f, nBait))
  subsetPrey_f = lapply(subsetBait_f, function(x) ppi_f %>% filter(bait %in% x) %>% pull(prey) %>% unique)
  names(subsetBait_f) = paste0(1:length(subsetBait_f))
  names(subsetPrey_f) = paste0(1:length(subsetPrey_f))
  result = list("subsetBait" = subsetBait_f,
                "subsetPrey" = subsetPrey_f,
                "setInfo" = c("nBait" = nBait, "nIter" = nIter))
  return(result)
}

# function that takes list of bait and prey, and calculates OR and p value for prey and calculates OR and p value for prey and exome-ASD102-prey
### Contignency table (individual counts): damaging mutation in geneset, no damaging mutation in genesent, ASD proband, ASD sibling.
bw.subsetPreyOR <- function(preyVector_f, geneAnnotations, numProbands, numSiblings){
  ### update annotation table
  annotation_f = geneAnnotations %>% mutate(prey = ifelse(hgnc_symbol %in% preyVector_f, 1, 0))
  
  ### calculate prey (excluding asd102) OR
  numProbands_variant = annotation_f %>% filter(prey == 1 & ASD102 ==0 & numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
  numProbands_noVariant = numProbands - numProbands_variant
  numSiblings_variant = annotation_f %>% filter(prey == 1 &  ASD102 ==0 & numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
  numSiblings_noVariant = numSiblings - numSiblings_variant
  dat_f = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
  stats1 = fisher.test(dat_f, alternative = "greater")
  stats2 = fisher.test(dat_f, alternative = "less")
  
  ### calculate exome-prey-asd102 OR
  numProbands_variant = annotation_f %>% filter(prey ==0 & ASD102 == 0 & numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
  numProbands_noVariant = numProbands - numProbands_variant
  numSiblings_variant = annotation_f %>% filter(prey == 0 & ASD102 == 0 & numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
  numSiblings_noVariant = numSiblings - numSiblings_variant
  dat2_f = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
  stats3 = fisher.test(dat2_f, alternative = "greater")
  stats4 = fisher.test(dat2_f, alternative = "less")
  
  ### calculate prey (including asd102 that happen to be prey) OR
  numProbands_variant = annotation_f %>% filter(prey == 1 & numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
  numProbands_noVariant = numProbands - numProbands_variant
  numSiblings_variant = annotation_f %>% filter(prey == 1  & numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
  numSiblings_noVariant = numSiblings - numSiblings_variant
  dat_f3 = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
  stats5 = fisher.test(dat_f3, alternative = "greater")
  stats6 = fisher.test(dat_f3, alternative = "less")
  
  ### calculate exome-prey(including asd102 prey) OR
  numProbands_variant = annotation_f %>% filter(prey ==0 & numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
  numProbands_noVariant = numProbands - numProbands_variant
  numSiblings_variant = annotation_f %>% filter(prey == 0  & numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
  numSiblings_noVariant = numSiblings - numSiblings_variant
  dat2_f = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
  stats7 = fisher.test(dat2_f, alternative = "greater")
  stats8 = fisher.test(dat2_f, alternative = "less")
  
  results = data.frame("geneset" = c("prey", "exomeMinusPreyASD102", "allPrey", "exomeMinusAllPrey"),
                       "pval" = c(stats1$p.value, stats3$p.value, stats5$p.value, stats7$p.value),
                       "CIlower" = c(stats1$conf.int[1], stats3$conf.int[1], stats5$conf.int[1], stats7$conf.int[1]),
                       "CIupper" = c(stats2$conf.int[2], stats4$conf.int[2], stats6$conf.int[2], stats8$conf.int[2]),
                       "OR" = c(stats1$estimate, stats3$estimate, stats5$estimate, stats7$estimate))
  rownames(results) = c("prey", "exomeMinusPreyASD102", "allPrey", "exomeMinusAllPrey")
  return(results)
}

# wrapper function that creates subsets of prey for a given number of bait, and returns pval/OR for prey and exome-ASD102-prey
bw.subsetPreyORList <- function(ppi_f, nBait=5, nIter=10, geneAnnotations, numProbands, numSiblings){
  subset_f = bw.subsetPPI(ppi_f, nBait, nIter)
  subsetPrey_f = subset_f$subsetPrey
  temp1 = lapply(subsetPrey_f, function(x) bw.subsetPreyOR(x, geneAnnotations, numProbands, numSiblings))
  #annotate with iteration number, number of bait, number of prey
  temp2 = lapply(1:length(temp1), function(x) temp1[[x]] %>% mutate("nBait" = nBait, iteration = names(temp1)[x]))
  temp3 = do.call(rbind, temp2)
  return(temp3)
}


# ------------------------------------------------------------------------------
#Create downsampled sets of ASD-PPI data
# ------------------------------------------------------------------------------
#Downsample the number of baits (randomly selecting sets of size 1 to 100 baits, 1000 iterations for each set size). For each PPI subset, calculate OR for prey (-hcASDASD102) and exome (-prey, -hcASD102)
if(!file.exists(paste0(od, "downsampledPPI100_preyOR_SSCsubset.Rdata"))){
  set.seed(28)
  numProbands = individualAnnotations %>% filter(Affected_Status ==2) %>% pull(Child_ID) %>% unique %>% length 
  numSiblings = individualAnnotations %>% filter(Affected_Status ==1) %>% pull(Child_ID) %>% unique %>% length
  subsetSummary = lapply(1:100, function(x) bw.subsetPreyORList(ppi, nBait = x, nIter = 1000, geneAnnotations, numProbands, numSiblings))
  # save data
  save(subsetSummary, file=paste0(od, "downsampledPPI100_preyOR_SSCsubset.Rdata"))
}
load(file=paste0(od, "downsampledPPI100_preyOR_SSCsubset.Rdata"))

# ------------------------------------------------------------------------------
#Make plots
# ------------------------------------------------------------------------------
#Format data for plotting
###Geneset annotations are:  
  #  prey: ASD-PPI prey, excluding ASD102 genes  
  # allPrey:ASD-PPI prey, includes 31 ASD102 genes  
  # exomeMinusPreyASD102: all genes in autosomal exome, excluding ASD-PPI prey and ASD102 genes  
  # exomeMinusAllPrey:all genes in autosomal exome, excluding allPrey  
summary = do.call(rbind, subsetSummary) %>%
  mutate(OR = replace(OR, which(OR == Inf), NA))

toPlotAll = summary %>% dplyr::group_by(geneset, nBait) %>% 
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
                   se.OR = sd.OR/sqrt(num)) %>%
  mutate(SElower.OR = median.OR - sd.OR,
         SEupper.OR = median.OR + sd.OR)


#plot on -log10 scale
### Note: for allPrey, there are some instances where median.pval - se.pval <0 (which makes no sense). For these instances, I set SElower.pval to be 0, and when we converted to log scale, I did not show these infinity values
toPlotAll.pval = toPlotAll %>% mutate(SElower.pval = median.pval - se.pval,
                                      SEupper.pval = median.pval + se.pval,
                                      SEupper.OR = median.OR - se.OR,
                                      SElower.OR = median.OR - se.OR,
                                      SDupper.OR = median.OR - se.OR,
                                      SDlower.OR = median.OR - se.OR) %>%
  mutate(SElower.pval = ifelse(SElower.pval < 0, 0, SElower.pval)) %>% ## if pval - 1se is less than 0, replace with 0
  mutate(log10.median.pval = -log10(median.pval),
         log10.SElower.pval = -log10(SElower.pval),
         log10.SEupper.pval = -log10(SEupper.pval)) %>%
  #mutate(log10.SElower.pval = ifelse(log10.SElower.pval == Inf, log10.median.pval, log10.SElower.pval)) %>% # if -log10(pval - 1se) is infinity, replace with log10.median.pval 
  mutate(geneset = factor(geneset, levels = c("exomeMinusPreyASD102", "prey", "exomeMinusAllPrey", "allPrey")))


### define variables
startPreyOR = toPlotAll.pval %>% filter(geneset == "prey") %>% filter(nBait==1) %>% pull(median.OR) %>% signif(., 3)
endPreyOR = toPlotAll.pval %>% filter(geneset == "prey") %>% filter(nBait==100) %>% pull(median.OR) %>% signif(., 3)
startExomeMinusPreyASD102OR = toPlotAll.pval %>% filter(geneset == "exomeMinusPreyASD102") %>% filter(nBait==1) %>% pull(median.OR) %>% signif(., 3)
endExomeMinusPreyASD102OR = toPlotAll.pval %>% filter(geneset == "exomeMinusPreyASD102") %>% filter(nBait==100) %>% pull(median.OR) %>% signif(., 3)
startPreyAllOR = toPlotAll.pval %>% filter(geneset == "allPrey") %>% filter(nBait==1) %>% pull(median.OR) %>% signif(., 3)
endPreyAllOR = toPlotAll.pval %>% filter(geneset == "allPrey") %>% filter(nBait==100) %>% pull(median.OR) %>% signif(., 3)
startExomeMinusPreyAllOR = toPlotAll.pval %>% filter(geneset == "exomeMinusAllPrey") %>% filter(nBait==1) %>% pull(median.OR) %>% signif(., 3)
endExomeMinusPreyAllOR = toPlotAll.pval %>% filter(geneset == "exomeMinusAllPrey") %>% filter(nBait==100) %>% pull(median.OR) %>% signif(., 3)

exomeAll_OR = summaryDF_d %>% filter(group=="exome") %>% pull(OR) %>%signif(.,3)
exomeAll_pval = summaryDF_d %>% filter(group=="exome") %>% pull(pval) %>%signif(.,3)
prey_pval = summaryDF_d  %>% filter(group=="prey") %>% pull(pval) %>%signif(.,3)
exomeMinusPreyASD102_pval = summaryDF_d  %>% filter(group=="exomeMinusPreyASD102") %>% pull(pval) %>%signif(.,3)
preyAll_pval = toPlotAll.pval %>% filter(nBait == 100, geneset == "allPrey") %>% pull(median.pval) %>% signif(., 3)
exomeMinusPreyAll_pval = toPlotAll.pval %>% filter(nBait == 100, geneset == "exomeMinusAllPrey") %>% pull(median.pval) %>% signif(., 3)

#### plot prey (no bait) and exomeMinusPreyASD102
pf = toPlotAll.pval %>% filter(geneset %in% c("prey", "exomeMinusPreyASD102")) %>%
  mutate(geneset = recode_factor(geneset, "exomeMinusPreyASD102"="Exome (- Prey, - hcASD102)",  "prey"="Prey (- hcASD102)")) %>%
  ggplot(aes(x=nBait, y=log10.median.pval, col=geneset)) + geom_line() + 
  geom_ribbon(aes(ymin = log10.SElower.pval, ymax = log10.SEupper.pval, fill = geneset), alpha = 0.1, linetype = "dotted") + 
  labs(x = "Number of bait", y = "-log10(median pval)", title = "Effect of enlarging ASD-PPI dataset on ability to identify ASD risk genes",
       subtitle = "Odds of SSC damaging mutations in prey vs exome\ncomparing ASD cases and sibling controls, individual count",
       caption=paste0("prey = ASD-PPI prey, excluding ASD102 genes\nexomeMinusPreyASD102 = exome genes, excluding prey and ASD102 genes
                        OR of entire exome = ", exomeAll_OR, " (p = ", exomeAll_pval, "); OR of exome - prey - ASD102 = ", endExomeMinusPreyASD102OR, " (p = ", exomeMinusPreyASD102_pval, "); OR of prey = ", endPreyOR,  " (p = ", prey_pval, ")")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + geom_text(aes(2, -log10(0.05), label = "p = 0.05", vjust = 1.5), col = "black") + 
  geom_vline(xintercept = 52, linetype = "dashed", color = "darkgray") + geom_text(aes(52, 3.25, label = "nBait = 52", hjust = -0.1), col = "darkgray")+
  geom_vline(xintercept = 16, linetype = "dashed", color = "darkgray") + geom_text(aes(16, 3.25, label = "nBait = 16", hjust = -0.1), col = "darkgray")+
  geom_text(aes(4, 0.17, label = paste0("OR = ", startPreyOR)), color = "#00BFC4") + 
  geom_text(aes(96, 2.75, label = paste0("OR = ", endPreyOR)), color = "#00BFC4") +
  geom_text(aes(4, 2.93, label = paste0("OR = ", startExomeMinusPreyASD102OR)), color = "#F8766D") + 
  geom_text(aes(96, 1.42, label = paste0("OR = ", endExomeMinusPreyASD102OR)), color = "#F8766D") +
  theme_bw() +  theme(legend.position='bottom')

pf
#ggsave(paste0(od, "FET_SSCsubset_downsampling_asdPPI.pdf"), width = 6.5, height = 7.25)
ggsave(paste0(od, "FET_SSCsubset_downsampling_asdPPI.pdf"), width = 5, height = 6.5)

# save data used to make plots
temp = toPlotAll.pval %>% filter(geneset %in% c("prey", "exomeMinusPreyASD102"))
write.csv(temp, file=paste0(od, "FET_subset_downsampling_asdPPI.csv"), row.names=F)
