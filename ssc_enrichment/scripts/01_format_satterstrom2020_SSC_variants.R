################################################################################
# Download and format Satterstrom et al 2020 variant data
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

#### Import 293T background genes  
#Import proteins detected in HEK293T from Bekker Jensen et al 2017 (PMID: 28601559) - 11,161  genes
load(paste0(rd, "hek293tProteome/BekkerJensen2017_293tProteome.RData"))

####2c. Import hcASD102 genes
load(paste0(rd, "asd102/asd102.Rdata"))
asd102 = asd102$hgnc

### Download Satterstrom 2020 (PMID: 31981491) mutation data 
#####Download Satterstrom 2020 Tables S1 and S2.  
#####Table S1: final list of high confidence de novo variants included 15,789 de novo variants from 6,430 probands and 2,179 unaffected children.  
#####Table S2: This file displays the output of the TADA gene discovery model as well as the variant counts and mutation rates that were inputs for it and the novelty of the findings. 
if(!file.exists(paste0(dd, "Satterstrom_2020_Table_S2_ASD_genes.xlsx"))){
  download.file(url = "https://www.cell.com/cms/10.1016/j.cell.2019.12.036/attachment/44aca411-6be3-4158-b1d6-6b339d60136d/mmc1.xlsx",
                destfile = paste0(dd, "Satterstrom_2020_Table_S1_ASD_denovo_genes.xlsx"))
  download.file(url = "https://www.cell.com/cms/10.1016/j.cell.2019.12.036/attachment/cc634b34-b459-4aff-acf0-cb58a3013bd1/mmc2.xlsx",
                destfile = paste0(dd, "Satterstrom_2020_Table_S2_ASD_genes.xlsx"))
}


### Import Satterstrom 2020 Table S1: The final list of high confidence de novo variants included 15,789 de novo variants from 6,430 probands and 2,179 unaffected children 
denovo_all <- read_excel(path = paste0(dd, "Satterstrom_2020_Table_S1_ASD_denovo_genes.xlsx"), sheet = "De novo variants", na = c(".", ""),
                         col_types = c("guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess",	"guess",	"guess",	"numeric",	"guess", "guess",	"guess",	"text",	"text",	"numeric",	"guess",	"guess",	"numeric",	"numeric",	"numeric",	"numeric",	"numeric",	"numeric",	"numeric",	"numeric",	"numeric",	"guess",	"guess",	"guess",	"guess"))
denovo_all <- denovo_all %>% mutate(hgnc_symbol = alias2SymbolTable(GENE_NAME, species = "Hs")) %>% filter(!is.na(hgnc_symbol))

### Import Satterstrom 2020 Table S2: results from TADA, including input variant counts and mutation rates
asdMutations <- read_excel(path = paste0(dd, "Satterstrom_2020_Table_S2_ASD_genes.xlsx"), sheet = 2, na = c(".", ""))
### update HGNC gene symbols
asdMutations <- asdMutations %>% mutate(hgnc_symbol = alias2SymbolTable(hugoGene, species = "Hs")) %>% filter(!is.na(hgnc_symbol))

# ------------------------------------------------------------------------------
# Create gene annotations dataframe
# ------------------------------------------------------------------------------

#Annotate HEK293T proteome with ppi_bait and prey status
hek293tProteome_o = hek293tProteome 
hek293tProteome = hek293tProteome_o  %>% dplyr::select(ensembl_gene_id, hgnc_symbol) %>%  
  mutate(hek293TproteinExpr = 1,
         bait = ifelse(hgnc_symbol %in% ppi_bait, 1, 0),
         prey = ifelse(hgnc_symbol %in% ppi_prey, 1, 0),
         preyAll = ifelse(hgnc_symbol %in% unique(ppi$prey), 1, 0))

#Create gene annotation table including all HEK293T proteome genes AND ASD-PPI genes (23 ASD-PPI genes were missing from HEK293T proteome)
genesToAdd = c(ppi_bait[!(ppi_bait %in% hek293tProteome$hgnc_symbol)], ppi_prey[!ppi_prey %in% hek293tProteome$hgnc_symbol]) %>% unique
genesToAdd = data.frame("hgnc_symbol" = genesToAdd)

geneAnnotations = bind_rows(hek293tProteome, genesToAdd) %>%
  mutate(bait = ifelse(hgnc_symbol %in% ppi_bait, 1, 0),
         prey = ifelse(hgnc_symbol %in% ppi_prey, 1, 0),
         preyAll = ifelse(hgnc_symbol %in% ppi$prey, 1, 0)) %>%
  filter(!is.na(hgnc_symbol)) %>%
  mutate(hek293TproteinExpr = ifelse(hek293TproteinExpr ==1 | prey ==1, 1, 0) ) #make sure prey genes are annotated as being expressed in HEK293T cells


#ANNOTATE GENES WITH DATA FROM TABLE S2
toCombine = asdMutations %>% dplyr::select(hgnc_symbol, case.ptv, control.ptv, dn.ptv, dn.ptv, dn.misa, dn.misb) %>%
  dplyr::rename(proband.dnPTV=dn.ptv, proband.dnMisA = dn.misa, proband.dnMisB=dn.misb) %>%
  filter(!duplicated(hgnc_symbol)) %>%
  mutate(associationGene = 1)
# update gene annotations
geneAnnotations = full_join(geneAnnotations, toCombine, by="hgnc_symbol") %>% filter(!duplicated(hgnc_symbol)) %>% relocate(associationGene, .after = hek293TproteinExpr) %>%
  mutate(bait = replace(bait, is.na(bait), 0),
         prey = replace(prey, is.na(prey), 0),
         ASD102 = ifelse(hgnc_symbol %in% asd102, 1, 0),
         hek293TproteinExpr = replace(hek293TproteinExpr, is.na(hek293TproteinExpr), 0),
         associationGene = replace(associationGene, is.na(associationGene), 0))

#ANNOTATE GENES WITH DE NOVO VARIANTS FROM SSC
####Use data from Satterstrom 2020 Table S1, which is the final list of high confidence de novo variants included 15,789 de novo variants from 6,430 probands and 2,179 unaffected children
#### Filter for only variants from Simon's Simplex Collection
# - Do this by extracting samples for which Phenotype_ID is in format of 1****.p1 (proband) or 1****.s1 (sibling)  
# - Note: there is another ID column (Child_ID) which is annotated differently
    # total samples: 3438
    # probands: 1971
    # siblings: 1467 siblings
ssc_samples = denovo_all %>% pull(Phenotype_ID) %>% str_extract(., "\\d{5}\\.[:lower:]\\d{1}") %>% na.omit %>% unique
denovo = denovo_all %>% filter(Phenotype_ID %in% ssc_samples)

# Create list genes with de novo mutations (dnm) identified in probands versus siblings, classified by ASD status (proband or sibling) and mutation effect (PTV, damagingMissense, damaging, missense, synonymous). Definitions for mutation effects:   
# 1. damagingMissense: Polyphen Mis3 (damaging) or MPC MisB (MPC >=2)  
# 2. PTV: frameshift, stop gained, or canonical splice site disruption
# 3. damaging: PTV + damagingMissense  
# 4. missense: missense coding mutations (note: use VEP_functional_class_canonical (the VEP_functional_class_canonical_simplified includes some noncoding missense variants under 'missense_variant' label)
# 5. synonymous: synonymous coding mutations (note: use VEP_functional_class_canonical (the VEP_functional_class_canonical_simplified includes some splice site variatns in "synonymous' label")
damagingMissense =  denovo %>% 
  dplyr::select(hgnc_symbol, Child_ID, Affected_Status, VEP_functional_class_canonical, VEP_functional_class_canonical_simplified, Coding, Polyphen_prediction, MPC, Variant_type) %>%
  filter(MPC >=2 | Polyphen_prediction == "probably_damaging")

PTV =  denovo %>% 
  dplyr::select(hgnc_symbol, Child_ID, Affected_Status, VEP_functional_class_canonical, VEP_functional_class_canonical_simplified, Coding, Polyphen_prediction, MPC, Variant_type) %>%
  filter(VEP_functional_class_canonical_simplified %in% c("frameshift_variant", "stop_gained", "splice_acceptor_variant", "splice_donor_variant"))

damaging = denovo %>% 
  dplyr::select(hgnc_symbol, Child_ID, Affected_Status, VEP_functional_class_canonical, VEP_functional_class_canonical_simplified, Coding, Polyphen_prediction, MPC, Variant_type) %>%
  filter(VEP_functional_class_canonical_simplified %in% c("frameshift_variant", "stop_gained", "splice_acceptor_variant", "splice_donor_variant") |
           MPC >=2 | Polyphen_prediction == "probably_damaging")

missense = denovo %>% dplyr::select(hgnc_symbol, Child_ID, Affected_Status, VEP_functional_class_canonical) %>%
  filter(VEP_functional_class_canonical == "missense_variant")

synonymous  = denovo %>% dplyr::select(hgnc_symbol, Child_ID, Affected_Status, VEP_functional_class_canonical) %>%
  filter(VEP_functional_class_canonical == "synonymous_variant")

# list with number of dnm
dnm_counts = list(
  proband.dnm.damagingMissense = damagingMissense %>% filter(Affected_Status ==2) %>% pull(hgnc_symbol) %>% table %>% as.data.frame,
  sibling.dnm.damagingMissense = damagingMissense %>% filter(Affected_Status ==1) %>% pull(hgnc_symbol) %>% table %>% as.data.frame,
  proband.dnm.PTV = PTV %>% filter(Affected_Status ==2) %>% pull(hgnc_symbol) %>% table %>% as.data.frame,
  sibling.dnm.PTV = PTV %>% filter(Affected_Status ==1) %>% pull(hgnc_symbol) %>% table %>% as.data.frame,
  proband.dnm.damaging = damaging %>% filter(Affected_Status ==2) %>% pull(hgnc_symbol) %>% table %>% as.data.frame,
  sibling.dnm.damaging = damaging %>% filter(Affected_Status ==1) %>% pull(hgnc_symbol) %>% table %>% as.data.frame,
  proband.dnm.missense = missense %>% filter(Affected_Status ==2) %>% pull(hgnc_symbol) %>% table %>% as.data.frame,
  sibling.dnm.missense = missense %>% filter(Affected_Status ==1) %>% pull(hgnc_symbol) %>% table %>% as.data.frame,
  proband.dnm.synonymous = synonymous %>% filter(Affected_Status ==2) %>% pull(hgnc_symbol) %>% table %>% as.data.frame,
  sibling.dnm.synonymous = synonymous %>% filter(Affected_Status ==1) %>% pull(hgnc_symbol) %>% table %>% as.data.frame
)
dnm_counts = lapply(dnm_counts, setNames, c("hgnc_symbol", "counts"))

# Update geneAnnotations dataframe of genes with de novo counts from SSC
#Annotate 293T proteome genes
geneAnnotations = geneAnnotations %>%
  left_join(., dnm_counts$proband.dnm.damagingMissense, by = "hgnc_symbol") %>% mutate(counts = replace(counts, is.na(counts), 0)) %>% dplyr::rename(proband.dnm.damagingMissense=counts) %>%
  left_join(., dnm_counts$sibling.dnm.damagingMissense, by = "hgnc_symbol") %>% mutate(counts = replace(counts, is.na(counts), 0)) %>% dplyr::rename(sibling.dnm.damagingMissense=counts) %>%
  left_join(., dnm_counts$proband.dnm.PTV, by = "hgnc_symbol") %>% mutate(counts = replace(counts, is.na(counts), 0)) %>% dplyr::rename(proband.dnm.PTV=counts) %>%
  left_join(., dnm_counts$sibling.dnm.PTV, by = "hgnc_symbol") %>% mutate(counts = replace(counts, is.na(counts), 0)) %>% dplyr::rename(sibling.dnm.PTV=counts) %>%
  left_join(., dnm_counts$proband.dnm.damaging, by = "hgnc_symbol") %>% mutate(counts = replace(counts, is.na(counts), 0)) %>% dplyr::rename(proband.dnm.damaging=counts) %>%
  left_join(., dnm_counts$sibling.dnm.damaging, by = "hgnc_symbol") %>% mutate(counts = replace(counts, is.na(counts), 0)) %>% dplyr::rename(sibling.dnm.damaging=counts) %>%
  left_join(., dnm_counts$proband.dnm.missense, by = "hgnc_symbol") %>% mutate(counts = replace(counts, is.na(counts), 0)) %>% dplyr::rename(proband.dnm.missense=counts) %>%
  left_join(., dnm_counts$sibling.dnm.missense, by = "hgnc_symbol") %>% mutate(counts = replace(counts, is.na(counts), 0)) %>% dplyr::rename(sibling.dnm.missense=counts) %>%
  left_join(., dnm_counts$proband.dnm.synonymous, by = "hgnc_symbol") %>% mutate(counts = replace(counts, is.na(counts), 0)) %>% dplyr::rename(proband.dnm.synonymous=counts) %>%
  left_join(., dnm_counts$sibling.dnm.synonymous, by = "hgnc_symbol") %>% mutate(counts = replace(counts, is.na(counts), 0)) %>% dplyr::rename(sibling.dnm.synonymous=counts)

geneAnnotations %>% dplyr::select(!c(1:2)) %>% filter(associationGene == 1) %>% apply(., 2, sum)

#Create individualAnnotation data frame in which individual de novo mutatiosn are annotated by functional class and bait/prey status to enable downstream individual-level analysis
individualAnnotations = denovo %>%
  dplyr::select(hgnc_symbol, Child_ID, Affected_Status, VEP_functional_class_canonical, VEP_functional_class_canonical_simplified, Coding, Polyphen_prediction, MPC, Variant_type) %>%
  mutate(bait = ifelse(hgnc_symbol %in% ppi_bait, 1, 0),
         prey = ifelse(hgnc_symbol %in% ppi_prey, 1, 0),
         preyAll = ifelse(hgnc_symbol %in% ppi$prey, 1, 0),
         ASD102 = ifelse(hgnc_symbol %in% asd102, 1, 0),
         hek293TproteinExpr = ifelse(hgnc_symbol %in% ppi_prey | hgnc_symbol %in% hek293tProteome$hgnc_symbol, 1, 0)) %>%
  mutate(dnm.damagingMissense = ifelse(MPC >=2 | Polyphen_prediction == "probably_damaging", 1, 0),
  dnm.PTV = ifelse(VEP_functional_class_canonical_simplified %in% c("frameshift_variant", "stop_gained", "splice_acceptor_variant", "splice_donor_variant"), 1, 0),
  dnm.damaging = ifelse(VEP_functional_class_canonical_simplified %in% c("frameshift_variant", "stop_gained", "splice_acceptor_variant", "splice_donor_variant") | MPC >=2 | Polyphen_prediction == "probably_damaging", 1, 0),
  dnm.missense = ifelse(VEP_functional_class_canonical == "missense_variant", 1, 0),
  dnm.synonymous = ifelse(VEP_functional_class_canonical == "synonymous_variant", 1, 0)) %>%
  mutate(dnm.damagingMissense = replace(dnm.damagingMissense, is.na(dnm.damagingMissense), 0),
         dnm.damaging = replace(dnm.damaging, is.na(dnm.damaging), 0)) %>%
  dplyr::select(hgnc_symbol, Child_ID, Affected_Status, hek293TproteinExpr, bait, prey, preyAll, ASD102, dnm.damagingMissense, dnm.PTV, dnm.damaging, dnm.missense, dnm.synonymous)

individualAnnotations[,-c(1:3)] %>% apply(.,2,sum)


# update geneAnnotations data frame with 2 columns containining information re: number of probands with a damaging variant and number of siblings with a damaging variant
numProbands_damagingVar = individualAnnotations %>% filter(Affected_Status ==2 & dnm.damaging ==1) %>% dplyr::select(hgnc_symbol, Child_ID) %>% dplyr::group_by(hgnc_symbol) %>% dplyr::summarize(numProband.damaging = n_distinct(Child_ID))
numSiblings_damagingVar = individualAnnotations %>% filter(Affected_Status ==1 & dnm.damaging ==1) %>% dplyr::select(hgnc_symbol, Child_ID) %>% dplyr::group_by(hgnc_symbol) %>% dplyr::summarize(numSibling.damaging = n_distinct(Child_ID))

geneAnnotations = left_join(geneAnnotations, numProbands_damagingVar, by = "hgnc_symbol") %>%
  mutate(numProband.damaging = replace(numProband.damaging, is.na(numProband.damaging), 0)) %>%
  left_join(. , numSiblings_damagingVar, by = "hgnc_symbol") %>%
  mutate(numSibling.damaging = replace(numSibling.damaging, is.na(numSibling.damaging), 0)) %>%
  relocate(ASD102, .after=hgnc_symbol)

# ------------------------------------------------------------------------------
# Export geneAnnotations and individualAnnotations to csv file
# ------------------------------------------------------------------------------
write.csv(geneAnnotations, file=paste0(od, "Satterstrom2020_geneAnnotations_SSC_caseControl_asdPPI100.csv"), row.names=F)
write.csv(individualAnnotations, file=paste0(od, "Satterstrom2020_SSC_DNM.csv"), row.names=F)

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
