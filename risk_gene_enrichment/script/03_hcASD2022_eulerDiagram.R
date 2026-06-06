# Goal: evaluate overlap between hcASD risk gene sets
# Authors: Belinda Wang

################################################################################
# setup workspace
################################################################################
library(tidyverse)
library(eulerr)

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
# make euler diagram
################################################################################

hcASDlist = list("hcASD" = genesMetadata %>% filter(asd102_satterstrom2020==1) %>% pull(hgnc_symbol),
                 "Fu 2022" = genesMetadata %>% filter(asd255_fu2022==1) %>% pull(hgnc_symbol),
                 "Trost 2022" = genesMetadata %>% filter(asd134_trost2022==1) %>% pull(hgnc_symbol),
                 "Zhou 2022" = genesMetadata %>% filter(asd72_zhou2022==1) %>% pull(hgnc_symbol))

plot(euler(hcASDlist, shape = "ellipse"), quantities = TRUE)
#plot(euler(hcASDlist[2:4], shape = "ellipse"), quantities = TRUE)

ggsave(paste0(od, "2022_hcASD_geneset_overlap.pdf"), width=6, height=5)
