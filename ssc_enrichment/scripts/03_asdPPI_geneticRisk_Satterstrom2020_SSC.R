################################################################################
# Assess for evidence of ASD genetic risk in ASD-PPI prey
###############################################################################

# ------------------------------------------------------------------------------
#set up workspace 
# ------------------------------------------------------------------------------
#*** [ ] remove other figures


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


# ------------------------------------------------------------------------------
#FET: damaging mutation in geneset, no damaging mutation in genesent, ASD proband, ASD sibling.
### genesets: prey, all exome genes, exomeNotASD102, exomeNotASDPPI, HEK293Tother
# ------------------------------------------------------------------------------

#FET, individual level, + de novo damaging variant or not, proband vs sibling
#Contignency table (individual counts): damaging mutation in geneset, no damaging mutation in genesent, ASD proband, ASD sibling.  
#Genesets of interest: prey, preyNotASD102, ASD102, all exome genes, exomeNotASD102, exome genes minus ASD-PPI  

# create list to store results
results_list_d = list()
numProbands = individualAnnotations %>% filter(Affected_Status ==2) %>% pull(Child_ID) %>% unique %>% length 
numSiblings = individualAnnotations %>% filter(Affected_Status ==1) %>% pull(Child_ID) %>% unique %>% length 

# Prey genes (allowing hcASD102)
numProbands_variant = geneAnnotations %>% filter(preyAll == 1 & numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
numProbands_noVariant = numProbands - numProbands_variant
numSiblings_variant = geneAnnotations %>% filter(preyAll == 1 & numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
numSiblings_noVariant = numSiblings - numSiblings_variant
dat = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
row.names(dat) = c("DamagingVar_prey", "No_damagingVar_prey")
colnames(dat) =  c("ASD", "sibling")
dat
stats1 = fisher.test(dat, alternative = "greater")
stats2 = fisher.test(dat, alternative = "less")
results_list_d[["prey"]] = c(stats1$p.value, stats1$conf.int[1], stats2$conf.int[2], stats1$estimate, unlist(dat))

# Prey, excluding ASD102
numProbands_variant = geneAnnotations %>% filter(prey == 1 & ASD102 != 1 & numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
numProbands_noVariant = numProbands - numProbands_variant
numSiblings_variant = geneAnnotations %>% filter(prey == 1 & ASD102 != 1 & numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
numSiblings_noVariant = numSiblings - numSiblings_variant
dat = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
row.names(dat) = c("DamagingVar_preyMinusASD102", "No_damagingVar_preyMinusASD102")
colnames(dat) =  c("ASD", "sibling")
dat
stats1 = fisher.test(dat, alternative = "greater")
stats2 = fisher.test(dat, alternative = "less")
results_list_d[["preyMinusASD102"]] = c(stats1$p.value, stats1$conf.int[1], stats2$conf.int[2], stats1$estimate, unlist(dat))

# ASD102 genes
numProbands_variant = geneAnnotations %>% filter(ASD102 == 1 & numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
numProbands_noVariant = numProbands - numProbands_variant
numSiblings_variant = geneAnnotations %>% filter(ASD102 == 1 & numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
numSiblings_noVariant = numSiblings - numSiblings_variant
dat = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
row.names(dat) = c("DamagingVar_ASD102", "No_damagingVar_ASD102")
colnames(dat) =  c("ASD", "sibling")
dat
stats1 = fisher.test(dat, alternative = "greater")
stats2 = fisher.test(dat, alternative = "less")
results_list_d[["asd102"]] = c(stats1$p.value, stats1$conf.int[1], stats2$conf.int[2], stats1$estimate, unlist(dat))

# all Satterstrom 2020 autosome genes
numProbands_variant = geneAnnotations %>% filter(numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
numProbands_noVariant = numProbands - numProbands_variant
numSiblings_variant = geneAnnotations %>% filter(numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
numSiblings_noVariant = numSiblings - numSiblings_variant
dat = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
row.names(dat) = c("DamagingVar_exome", "No_damagingVar_exome")
colnames(dat) =  c("ASD", "sibling")
dat
stats1 = fisher.test(dat, alternative = "greater")
stats2 = fisher.test(dat, alternative = "less")
results_list_d[["exome"]] = c(stats1$p.value, stats1$conf.int[1], stats2$conf.int[2], stats1$estimate, unlist(dat))

# all Satterstrom 2020 autosome genes, excluding ASD102
numProbands_variant = geneAnnotations %>% filter(ASD102 != 1 & numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
numProbands_noVariant = numProbands - numProbands_variant
numSiblings_variant = geneAnnotations %>% filter(ASD102 != 1 & numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
numSiblings_noVariant = numSiblings - numSiblings_variant
dat = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
row.names(dat) = c("DamagingVar_exomeMinusPreyASD102", "No_damagingVar_exomeMinusASD102")
colnames(dat) =  c("ASD", "sibling")
dat
stats1 = fisher.test(dat, alternative = "greater")
stats2 = fisher.test(dat, alternative = "less")
results_list_d[["exomeMinusASD102"]] = c(stats1$p.value, stats1$conf.int[1], stats2$conf.int[2], stats1$estimate, unlist(dat))

# all Satterstrom 2020 autosome genes, excluding ASD-PPI & ASD102
numProbands_variant = geneAnnotations %>% filter(prey != 1 & ASD102 != 1 & numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
numProbands_noVariant = numProbands - numProbands_variant
numSiblings_variant = geneAnnotations %>% filter(prey != 1 & ASD102 != 1 & numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
numSiblings_noVariant = numSiblings - numSiblings_variant
dat = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
row.names(dat) = c("DamagingVar_exomeMinusPreyASD102", "No_damagingVar_exomeMinusPreyASD102")
colnames(dat) =  c("ASD", "sibling")
dat
stats1 = fisher.test(dat, alternative = "greater")
stats2 = fisher.test(dat, alternative = "less")
results_list_d[["exomeMinusPreyASD102"]] = c(stats1$p.value, stats1$conf.int[1], stats2$conf.int[2], stats1$estimate, unlist(dat))

# hek293t (hek293t proteome)
numProbands_variant = geneAnnotations %>% filter(hek293TproteinExpr == 1 & numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
numProbands_noVariant = numProbands - numProbands_variant
numSiblings_variant = geneAnnotations %>% filter(hek293TproteinExpr == 1 & numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
numSiblings_noVariant = numSiblings - numSiblings_variant
dat = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
row.names(dat) = c("DamagingVar_hek293tMinusPrey", "No_damagingVar_hek293tMinusPrey")
colnames(dat) =  c("ASD", "sibling")
dat
stats1 = fisher.test(dat, alternative = "greater")
stats2 = fisher.test(dat, alternative = "less")
results_list_d[["hek293t"]] = c(stats1$p.value, stats1$conf.int[1], stats2$conf.int[2], stats1$estimate, unlist(dat))

# hek293tMinusPrey (hek293t proteome minus preyAll)
numProbands_variant = geneAnnotations %>% filter(preyAll == 0 & hek293TproteinExpr == 1 & numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
numProbands_noVariant = numProbands - numProbands_variant
numSiblings_variant = geneAnnotations %>% filter(preyAll == 0 & hek293TproteinExpr == 1 & numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
numSiblings_noVariant = numSiblings - numSiblings_variant
dat = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
row.names(dat) = c("DamagingVar_hek293tMinusPrey", "No_damagingVar_hek293tMinusPrey")
colnames(dat) =  c("ASD", "sibling")
dat
stats1 = fisher.test(dat, alternative = "greater")
stats2 = fisher.test(dat, alternative = "less")
results_list_d[["hek293tMinusPrey"]] = c(stats1$p.value, stats1$conf.int[1], stats2$conf.int[2], stats1$estimate, unlist(dat))

# hek293tMinusPreyASD102 (hek293t proteome minus prey and hcASD102)
numProbands_variant = geneAnnotations %>% filter(prey == 0 & ASD102 != 1 & hek293TproteinExpr == 1 & numProband.damaging > 0) %>% pull(numProband.damaging) %>% sum
numProbands_noVariant = numProbands - numProbands_variant
numSiblings_variant = geneAnnotations %>% filter(prey == 0 &  ASD102 != 1 & hek293TproteinExpr == 1 & numSibling.damaging > 0) %>% pull(numSibling.damaging) %>% sum
numSiblings_noVariant = numSiblings - numSiblings_variant
dat = data.frame(c(numProbands_variant, numProbands_noVariant), c(numSiblings_variant, numSiblings_noVariant))
row.names(dat) = c("DamagingVar_hek293tMinusPreyASD102", "No_damagingVar_hek293tMinusPreyASD102")
colnames(dat) =  c("ASD", "sibling")
dat
stats1 = fisher.test(dat, alternative = "greater")
stats2 = fisher.test(dat, alternative = "less")
results_list_d[["hek293tMinusPreyASD102"]] = c(stats1$p.value, stats1$conf.int[1], stats2$conf.int[2], stats1$estimate, unlist(dat))


# Plot data
summaryDF_d = do.call("rbind", results_list_d)
colnames(summaryDF_d) = c("pval", "CIlower", "CIupper", "OR", "a", "c", "b", "d")

summaryDF_d = summaryDF_d %>% as.data.frame %>% rownames_to_column("group") %>% 
  #mutate(adj.pval = pval*4) %>%
  #mutate(adj.pval = replace(adj.pval, adj.pval >1, 1)) %>%
  #mutate(adj.pval.text = paste0("p.adj = ", signif(adj.pval, 2))) %>%
  mutate(pval.text = paste0("p = ", signif(pval,2))) %>%
  # make factor to enable ordered plotting
  mutate(group = factor(group, levels = c( "hek293tMinusPreyASD102", "hek293tMinusPrey", "hek293t", "exomeMinusPreyASD102", "exomeMinusASD102", "preyMinusASD102", "prey", "asd102", "exome")))
summaryDF_d

# export summaryDF_d, which contains results from FET 
#FET: damaging mutation in geneset, no damaging mutation in genesent, ASD proband, ASD sibling.
### genesets: prey, all exome genes, exomeNotASD102, exomeNotASDPPI, HEK293Tother
write.csv(summaryDF_d, file=paste0(od, "FET_SSCsubset_results.csv"), row.names=F)

fig1a =  subset(summaryDF_d, group %in% c("exome", "exomeMinusASD102", "exomeMinusPreyASD102")) %>%
  mutate(group = recode_factor(group, "exomeMinusPreyASD102"="Exome (- Prey, hcASD102)", "exomeMinusASD102"="Exome (- hcASD102)", "exome"="Exome")) %>%
  mutate(p.adj = pval * nrow(.)) %>%
  mutate(padj.text = paste0("p.adj = ", signif(p.adj,2))) %>%
  ggplot( aes(x = group, y = OR, ymin = CIlower, ymax = CIupper)) +
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 1, color="gray") +
  geom_errorbar(width = 0.2) + 
  geom_point(size = 3) +
  labs(y = "Odds ratio", x = "",
       title = "Odds of SSC dn damaging variant vs not in geneset\ncomparing ASD cases and sibling controls",
       subtitle = "Fisher's exact test, one sided, greater\nIndividual count") +
  coord_flip() + theme_bw() + theme(axis.text = element_text(size = 12))+
  theme(legend.position="none") +
  geom_text(aes(label = padj.text), size = 4, vjust = 2, hjust = 0.5)
fig1a
ggsave(paste0(od, "FET_SSCsubset_exome.pdf"), width = 7, height = 3)

# fig1b = subset(summaryDF_d, group %in% c("prey",  "hek293tMinusPrey", "preyMinusASD102",  "hek293tMinusPreyASD102")) %>%
#   mutate(group = recode_factor(group,  "hek293tMinusPreyASD102"="HEK293T (- Prey, hcASD102)", "hek293tMinusPrey"="HEK293T (- Prey)","preyMinusASD102"="Prey (- hcASD102)", "prey"="Prey")) %>%
#   mutate(p.adj = pval * nrow(.)) %>%
#   mutate(padj.text = paste0("p.adj = ", signif(p.adj,2))) %>%
#   ggplot( aes(x = group, y = OR, ymin = CIlower, ymax = CIupper)) +
#   geom_hline(yintercept = 1.0, linetype = "dotted", size = 1, color="gray") +
#   geom_errorbar(width = 0.2) + 
#   geom_point(size = 3) +
#   labs(y = "Odds ratio", x = "",
#        title = "Odds of SSC dn damaging variant vs not in geneset\ncomparing ASD cases and sibling controls",
#        subtitle = "Fisher's exact test, one sided, greater\nIndivdual count") +
#   coord_flip() + theme_bw() + theme(axis.text = element_text(size = 12))+
#   theme(legend.position="none") +
#   geom_text(aes(label = padj.text), size = 4, vjust = 2, hjust = 0.5)

# fig1b = subset(summaryDF_d, group %in% c("prey", "hek293tMinusPreyASD102")) %>%
#   mutate(group = recode_factor(group,  "hek293tMinusPreyASD102"="HEK293T (- Prey, hcASD102)", "prey"="Prey")) %>%
#   mutate(p.adj = pval * nrow(.)) %>%
#   mutate(padj.text = paste0("p.adj = ", signif(p.adj,2))) %>%
#   ggplot( aes(x = group, y = OR, ymin = CIlower, ymax = CIupper)) +
#   geom_hline(yintercept = 1.0, linetype = "dotted", size = 1, color="gray") +
#   geom_errorbar(width = 0.2) + 
#   geom_point(size = 3) +
#   labs(y = "Odds ratio", x = "",
#        title = "Odds of SSC dn damaging variant vs not in geneset\ncomparing ASD cases and sibling controls",
#        subtitle = "Fisher's exact test, one sided, greater\nIndivdual count") +
#   coord_flip() + theme_bw() + theme(axis.text = element_text(size = 12))+
#   theme(legend.position="none") +
#   geom_text(aes(label = padj.text), size = 4, vjust = 2, hjust = 0.5)
# fig1b
# ggsave(paste0(od, "FET_SSCsubset_HEK293T_prey.pdf"), width = 7, height = 2.5)


fig1c = subset(summaryDF_d, group %in% c("prey", "preyMinusASD102", "hek293tMinusPreyASD102")) %>%
  mutate(group = recode_factor(group,  "hek293tMinusPreyASD102"="HEK293T (- Prey, hcASD102)", "preyMinusASD102" = "Prey (-hcASD102)", "prey"="Prey")) %>%
  mutate(p.adj = pval * nrow(.)) %>%
  mutate(padj.text = paste0("p.adj = ", signif(p.adj,2))) %>%
  ggplot( aes(x = group, y = OR, ymin = CIlower, ymax = CIupper)) +
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 1, color="gray") +
  geom_errorbar(width = 0.2) + 
  geom_point(size = 3) +
  labs(y = "Odds ratio", x = "",
       title = "Odds of SSC dn damaging variant vs not in geneset\ncomparing ASD cases and sibling controls",
       subtitle = "Fisher's exact test, one sided, greater\nIndivdual count") +
  coord_flip() + theme_bw() + theme(axis.text = element_text(size = 12))+
  theme(legend.position="none") +
  geom_text(aes(label = padj.text), size = 4, vjust = 2, hjust = 0.5)
fig1c
ggsave(paste0(od, "FET_SSCsubset_HEK293T_preyOnly.pdf"), width = 6, height = 3)
