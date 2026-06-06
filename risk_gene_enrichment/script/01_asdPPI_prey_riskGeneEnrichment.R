# Goal: evaluate whether ASD-PPI prey are enriched for risk genes that have been implicated in ASD, DD, or schizophrenia. 
# Authors: Belinda Wang

################################################################################
# setup workspace
################################################################################
library(tidyverse)
library(readxl)
library(limma)

# Run from the repository root.
wd <- "risk_gene_enrichment/"
dd <- "risk_gene_enrichment/data/"
od <- "risk_gene_enrichment/output/"

# import data
# - hgnc_symbol: HGNC gene symbol  
# - gene_gencodeV33: 
# - hek293Texpr: indicator for whether gene is one of the n=11167 proteins expressed in HEK293T cells (union of Bekker Jensen 2017 and prey100; 1=T; 0=F)  
# -  bait100: indicator for whether gene is one of the n=100 hcASD genes used as bait in ASD-PPI_100 (1=T; 0=F)  
# - prey100: indicator for whether gene is one of the n=1074 ASD-PPI_100 prey (1=T; 0=F). Note: bait that were identified as prey are included here.
# - asd102_satterstrom2020: indicator for whether gene is one of the n=102 hcASD genes identified in Satterstrom et al 2020 (1=T; 0=F)  
# - asd255_fu2022: indicator for whether gene is one of the n=255 hcASD genes identified in Fu et al 2022 (1=T; 0=F)  
# - asd134_trost2022: indicator for whether gene is one of the 134 genes with TADA FDR<0.1 in WES from MSSNG 2022 release (Trost 2022 PMID 36368308)
# - asd72_zhou2022: indicator for whether gene is one of the n=72 genes with study-wide significance (based on 5,754 constraint genes, p<8.69E-06) in  analysis of SPARK data in Zhou et al Nature 2022 (PMID 35982159)
# - dd285_kaplanis2020:  indicator for whether gene is one of the n=285 unique genes that are significantly associated with DD (after one-sided Bonferroni correction) in Kaplanis et al Nature 2020 (PMID 33057194) (1=T; 0=F)  
# - scz34_singh2022: indicator for whether gene is one of the n=34 genes with FDR<0.1 from SCHEMA WES  Singh et al Nature 2022 (PMID 35396579) (1=T; 0=F) 

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
# FET set 1
### Assess whether ASD-PPI_100 prey are enriched for risk genes implicated in other disorders (WES, with known assessed gene lists)
### Gene universe: HEK293T proteome + Satterstrom gene universe + disease geneset gene universe
################################################################################                  
#Here, we conduct FET in which we restrict the universe of genes to those that are 1) in HEK293T proteome (Bekker Jensen 2017 + ASD-PPI prey); 2) assessed in Satterstrom 2020; and 3) assessed in the disorder geneset WES study. 
# Genesets that we will focus on include:  
# - Fu et al 2022 (ASD WES)  
# - Singh et al 2022 (SCZ WES)  
# - Kaplanis et al 2020 (DD WES)  
# - Zhou et al 2022 (ASD WES, SPARK)
# - Trost et al 2022 (ASD WES, MSSNG)
# Contingency table: prey, not prey, in disease geneset, not in disease geneset

#Function to perform Fisher's exact test using different disease genesets
bw.FET.HEK293TwesUniverse <- function(f_genesMetadata = genesMetadata, geneUniverseColName = "fuSatterstrom", group1Colname = "prey100", group2Colname = "asd255_fu2022", results = results_list_a, compName){
  # make contingency table to perform fisher's exact test
  dat = data.frame(c(f_genesMetadata %>% filter(hek293Texpr == 1 & !!rlang::sym(geneUniverseColName) == 1 & !!rlang::sym(group1Colname) == 1 & !!rlang::sym(group2Colname) == 1) %>% nrow,
                     f_genesMetadata %>% filter(hek293Texpr == 1 & !!rlang::sym(geneUniverseColName) == 1 & !!rlang::sym(group1Colname) == 0 & !!rlang::sym(group2Colname) == 1) %>% nrow),
                   c(f_genesMetadata %>% filter(hek293Texpr == 1 & !!rlang::sym(geneUniverseColName) == 1 & !!rlang::sym(group1Colname) == 1 & !!rlang::sym(group2Colname) == 0) %>% nrow,
                     f_genesMetadata %>% filter(hek293Texpr == 1 & !!rlang::sym(geneUniverseColName) == 1 & !!rlang::sym(group1Colname) == 0 & !!rlang::sym(group2Colname) == 0) %>% nrow),
                   row.names = c("Prey", "293Tother"))
  colnames(dat) =  c(group2Colname, paste0("not.", group2Colname))
  
  #perform fisher's exact test (pvalue from on one-sided, greater. Also do one-sided, less to grab upper CI for plotting)
  stats1 = fisher.test(dat, alternative = "greater")
  stats2 = fisher.test(dat, alternative = "less")
  #save data
  results[[compName]] = c(stats1$p.value, stats1$conf.int[1], stats2$conf.int[2], stats1$estimate, unlist(dat))
  print(dat)
  return(results)
}

#Conduct FET
results_list_293Twes = list() # list to save reuslts from Fisher's exact test
results_list_293Twes  = bw.FET.HEK293TwesUniverse(genesMetadata, "fuSatterstrom", "prey100", "asd255_fu2022", results_list_293Twes, "asd255_fu2022")
results_list_293Twes  = bw.FET.HEK293TwesUniverse(genesMetadata, "singhSatterstrom", "prey100", "scz34_singh2022", results_list_293Twes, "scz34_singh2022")
results_list_293Twes  = bw.FET.HEK293TwesUniverse(genesMetadata, "kaplanisSatterstrom", "prey100", "dd285_kaplanis2020", results_list_293Twes, "dd285_kaplanis2020")
results_list_293Twes  = bw.FET.HEK293TwesUniverse(genesMetadata, "zhouSatterstrom", "prey100", "asd72_zhou2022", results_list_293Twes, "asd72_zhou2022")
results_list_293Twes  = bw.FET.HEK293TwesUniverse(genesMetadata, "fuSatterstrom", "prey100", "asd134_trost2022", results_list_293Twes, "asd134_trost2022")

summaryDF_293Twes = do.call("rbind", results_list_293Twes)
colnames(summaryDF_293Twes) = c("pval", "CIlower", "CIupper", "OR", "a", "c", "b", "d")
summaryDF_293Twes = summaryDF_293Twes %>% as.data.frame %>% rownames_to_column("group")
summaryDF_293Twes


# Make plots
# MAIN TEXT FIGURE: FET for prey enrichment of risk genes associated with ASD or SCZ from recent WES studies. 
###Gene universe = HEK293T proteme & Satterstrom 2020 & assessed in disease geneset WES, includes ASD102
toPlot = summaryDF_293Twes %>% filter(!group=="dd285_kaplanis2020") %>% 
  mutate(group = factor(group, levels = c("scz34_singh2022", "asd72_zhou2022", "asd134_trost2022", "asd255_fu2022"), labels = c("SCZ, Singh 2022 (n=34)", "ASD, Zhou 2022 (n=72)", "ASD, Trost 2022 (n=134)", "ASD, Fu 2022 (n=255)"))) %>%
  mutate(p.adj = pval*nrow(.)) %>%
  mutate(p.adj = ifelse(p.adj>1, 1, p.adj)) %>%
  mutate(padj.text = paste0("p.adj = ", signif(p.adj, 2)))

ggplot(toPlot, aes(x = group, y = OR, ymin = CIlower, ymax = CIupper)) +
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 1, color="gray") +
  geom_errorbar(width = 0.2) + 
  geom_point(size = 3) +
  labs(y = "Odds ratio", x = "",
       #title = "Enrichment of ASD or SCZ risk genes in ASD-PPI prey",
       subtitle = "Enrichment of ASD or SCZ risk genes in ASD-PPI prey\nFisher's exact test, one sided, greater") +
  coord_flip() + theme_bw() + 
  geom_text(aes(label = padj.text), size = 3, vjust = 2)
ggsave(paste0(od, "FET_asd255_scz34_asdPPIprey.pdf"), width=5.5, height=3)

################################################################################
# FET set 2
### Assess whether ASD-PPI_100 prey are enriched for risk genes implicated in ASD, DD, or SCZ
### Gene universe: HEK293T proteome
################################################################################ 

# Here, we will restrict the gene universe to the HEK293T proteome (union of proteins detected in Bekker-Jensen 2017 and ASD-PPI_100 prey). We have not further trimmed to genes that were assessed in the various individual studies because these data are not readily available for many genesets (ex: SFARI)
#Contingency table: prey, not prey, in disease geneset, not in disease geneset


#Function to perform Fisher's exact test using different disease genesets
bw.FET.HEK293Tuniverse <- function(f_genesMetadata = genesMetadata, group1Colname = "prey100", group2Colname = "asd255_fu2022", results = results_list_a, compName){
  # make contingency table to perform fisher's exact test
  dat = data.frame(c(f_genesMetadata %>% filter(hek293Texpr == 1 & !!rlang::sym(group1Colname) == 1 & !!rlang::sym(group2Colname) == 1) %>% nrow,
                     f_genesMetadata %>% filter(hek293Texpr == 1 & !!rlang::sym(group1Colname) == 0 & !!rlang::sym(group2Colname) == 1) %>% nrow),
                   c(f_genesMetadata %>% filter(hek293Texpr == 1 & !!rlang::sym(group1Colname) == 1 & !!rlang::sym(group2Colname) == 0) %>% nrow,
                     f_genesMetadata %>% filter(hek293Texpr == 1 & !!rlang::sym(group1Colname) == 0 & !!rlang::sym(group2Colname) == 0) %>% nrow),
                   row.names = c("Prey", "293Tother"))
  colnames(dat) =  c(group2Colname, paste0("not.", group2Colname))
  
  #perform fisher's exact test (pvalue from on one-sided, greater. Also do one-sided, less to grab upper CI for plotting)
  stats1 = fisher.test(dat, alternative = "greater")
  stats2 = fisher.test(dat, alternative = "less")
  #save data
  results[[compName]] = c(stats1$p.value, stats1$conf.int[1], stats2$conf.int[2], stats1$estimate, unlist(dat))
  print(dat)
  return(results)
}

#Conduct FET
results_list_293Texpr = list() # list to save reuslts from Fisher's exact test
results_list_293Texpr  = bw.FET.HEK293Tuniverse(genesMetadata, "prey100", "asd255_fu2022", results_list_293Texpr, "asd255_fu2022")
results_list_293Texpr  = bw.FET.HEK293Tuniverse(genesMetadata, "prey100", "sfari_all", results_list_293Texpr, "sfari_all")
results_list_293Texpr  = bw.FET.HEK293Tuniverse(genesMetadata, "prey100", "sfari_s", results_list_293Texpr, "sfari_s")
results_list_293Texpr  = bw.FET.HEK293Tuniverse(genesMetadata, "prey100", "sfari_tier1", results_list_293Texpr, "sfari_tier1")
results_list_293Texpr  = bw.FET.HEK293Tuniverse(genesMetadata, "prey100", "sfari_tier1s", results_list_293Texpr, "sfari_tier1s")
results_list_293Texpr  = bw.FET.HEK293Tuniverse(genesMetadata, "prey100", "sfari_tier2", results_list_293Texpr, "sfari_tier2")
results_list_293Texpr  = bw.FET.HEK293Tuniverse(genesMetadata, "prey100", "sfari_tier3", results_list_293Texpr, "sfari_tier3")
results_list_293Texpr  = bw.FET.HEK293Tuniverse(genesMetadata, "prey100", "asd134_trost2022", results_list_293Texpr, "asd134_trost2022")
results_list_293Texpr  = bw.FET.HEK293Tuniverse(genesMetadata, "prey100", "asd72_zhou2022", results_list_293Texpr, "asd72_zhou2022")
results_list_293Texpr  = bw.FET.HEK293Tuniverse(genesMetadata, "prey100", "scz34_singh2022", results_list_293Texpr, "scz34_singh2022")
results_list_293Texpr  = bw.FET.HEK293Tuniverse(genesMetadata, "prey100", "dd285_kaplanis2020", results_list_293Texpr, "dd285_kaplanis2020")

summaryDF_293Texpr = do.call("rbind", results_list_293Texpr)
colnames(summaryDF_293Texpr) = c("pval", "CIlower", "CIupper", "OR", "a", "c", "b", "d")
summaryDF_293Texpr = summaryDF_293Texpr %>% as.data.frame %>% rownames_to_column("group") %>%
  mutate(pval.text = paste0("p = ", signif(pval,2)))
summaryDF_293Texpr

#Make plots
# SUPPLEMENT FIGURE: FET for prey enrichment of risk genes associated with various genesets. 
### Gene universe = HEK293T proteome
toPlot2 = summaryDF_293Texpr %>% filter(group %in% c("asd255_fu2022", "sfari_all", "sfari_tier3", "sfari_tier2", "sfari_s", "sfari_tier1", "sfari_tier1s", "asd134_trost2022", "asd72_zhou2022", "scz34_singh2022", "dd285_kaplanis2020")) %>%
  mutate(group = factor(group, levels = c("scz34_singh2022", "sfari_all",  "sfari_tier3", "sfari_tier2", "sfari_tier1", "sfari_tier1s", "sfari_s", "dd285_kaplanis2020", "asd72_zhou2022", "asd134_trost2022", "asd255_fu2022"), 
                        labels = c("SCZ, Singh 2022 (n=34)", "SFARI all (n=1020)",  "SFARI category 3 (n=514)", "SFARI category 2 (n=219)", "SFARI category 1 (n=206)", "SFARI category 1 & syndromic (n=92)", "SFARI syndromic (n=230)", "DD, Kaplanis 2020 (n=285)", "ASD, Zhou 2022 (n=72)", "ASD, Trost 2022 (n=134)", "ASD, Fu 2022 (n=255)"))) %>%
  mutate(p.adj = pval*nrow(.)) %>%
  mutate(p.adj = ifelse(p.adj>1, 1, p.adj)) %>%
  mutate(padj.text = paste0("p.adj = ", signif(p.adj, 3)))


ggplot(toPlot2, aes(x = group, y = OR, ymin = CIlower, ymax = CIupper)) +
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 1, color="gray") +
  geom_errorbar(width = 0.2) + 
  geom_point(size = 3) +
  labs(y = "Enrichment ratio", x = "",
       #title = "Enrichment of varios disease genesets in ASD-PPI prey",
       subtitle = "Enrichment of various risk genes in ASD-PPI prey\nFisher's exact test, one sided, greater") +
  coord_flip() + theme_bw() + 
  geom_text(aes(label = padj.text), size = 3, vjust = 2)
ggsave(paste0(od, "FET_dzGenesets_asdPPIprey.pdf"), width=6, height=5)
       
