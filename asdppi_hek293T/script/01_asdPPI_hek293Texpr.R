################################################################################
# Assess relative protein expression of ASD-PPI bait and prey in HEK293T
################################################################################
# Run from the repository root.
wd <- "asdppi_hek293T/"
od <- "asdppi_hek293T/output/"
rd <- "data/"

#rm(list=ls())
# Set up workspace
library(tidyverse)
library(limma) #using alias2SymbolTable() function to update BrainSpan hgnc gene symbols
library(parallel)
library(ggpubr) ##using stat_compare_means() in plots

# import 293T proteome data        
### data from from Bekker-Jensen et al 2017 (PMID: 28601559)
### 293T-specific data extracted and formatted in 01_define_hek293tproteome.R
load(paste0(rd, "hek293tProteome/BekkerJensen2017_293tProteome.RData"))
hek293tProteome_o = hek293tProteome

# Import PPI data 
load(paste0(rd, "ppi100.RData"))
ppiFull = ppi100


# Extract bait and prey genes (defining prey genes to be exclusive of bait genes)
bait <- unique(ppiFull$bait)
prey <- setdiff(unique(ppiFull$prey), bait)


# Format 293T proteome dataset for plotting (annotate by bait, prey, other gene status; rank proteome expression level)
hek293tProteome_toPlot = hek293tProteome %>%
  mutate(geneset = hgnc_symbol) %>%
  mutate(geneset = replace(geneset, geneset %in% bait, "bait")) %>%
  mutate(geneset = replace(geneset, geneset %in% prey, "prey")) %>%
  mutate(geneset = replace(geneset, !geneset %in% c("bait", "prey"), "other")) %>%
  mutate(geneset = factor(geneset, levels = c("bait", "prey", "other"))) %>%
  mutate(iBAQ_ranked = rank(iBAQ_avg, ties.method = "average")) %>%
  mutate(LFQ_ranked = rank(LFQ_avg, ties.method = "average"))

# make plot
### boxplots comparing HEK293T expression of bait, prey, and hek293Tother
my_comparisons <- list(c("Bait", "Prey (-hcASD102)"), c("Prey (-hcASD102)", "Other"), c("Bait", "Other"))
fig = hek293tProteome_toPlot %>%
  mutate(geneset = recode_factor(geneset, 'bait'='Bait', 'prey'='Prey (-hcASD102)', 'other'="Other")) %>%
  mutate(log10_iBAQ_avg = log10(iBAQ_avg)) %>%
  mutate(log10_iBAQ_avg = ifelse(log10_iBAQ_avg == -Inf, 0, log10_iBAQ_avg)) %>%
  ggplot(aes(x = geneset, y = log10_iBAQ_avg, fill = geneset)) +
  geom_boxplot(show.legend = F) +
  #geom_jitter(size = 1, alpha = 0.15) +
  stat_compare_means(method="t.test", comparisons = my_comparisons, aes(label=..p.adj..))+
  labs(title = "HEK293T global proteome (Bekker-Jensen et al. 2017)\nAvg iBAQ expression", x = "Geneset") +
  theme_bw()

plot_list = list(fig)
pdf(paste0(od, "asdPPI_HEK293T_expression_boxplot.pdf"), width=5, height=6)
for (i in 1:length(plot_list)) {
  print(plot_list[[i]])
}
dev.off()

# # get numbers for supplementary table
# temp = hek293tProteome_toPlot %>%
#   mutate(log10_iBAQ_avg = log10(iBAQ_avg)) %>%
#   mutate(log10_iBAQ_avg = ifelse(log10_iBAQ_avg == -Inf, 0, log10_iBAQ_avg))
# 
# b = temp %>% filter(geneset=="bait") %>% pull(log10_iBAQ_avg)
# o =  temp %>% filter(geneset=="other") %>% pull(log10_iBAQ_avg)
# p = temp %>% filter(geneset=="prey") %>% pull(log10_iBAQ_avg)
# t.test(b, o)
# t.test(b,p)
# t.test(p,o)
