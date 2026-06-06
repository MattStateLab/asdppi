################################################################################
# Make bsRNAseq permutation plots
################################################################################

# ------------------------------------------------------------------------------
#set up workspace 
# ------------------------------------------------------------------------------
# Run from the repository root.
wd <- "bsRNAseq/"
dd <- "bsRNAseq/data/"
od <- "bsRNAseq/output/"
rd <- "data/"


###################################################################################
#Import and format data
###################################################################################
#import data for bsRNAseq,ASD-PPI_100 prey
load(file=paste0(od, "bsRNAseq_asdPPIprey_permutationData_formattedToPlot.RData"))
# format data (bsRNAseq, ASD-PPI_100 prey)
permuted =  bsRNAseq_iBAQ_Permuted_medianRanks_toPlot %>% select(starts_with("permuted")) %>% 
  apply(., 1, median)
toPlot_bsRNAseq_prey100 = bsRNAseq_iBAQ_Permuted_toPlot %>% filter(geneset == "prey") %>%
  mutate(log10.pVal_permuted = -log10(pVal_permuted+0.00001)) %>%
  mutate(permuted.medRank = permuted,
         rankdiff = medRank - permuted.medRank) %>%
  mutate(regionList = str_replace(regionList, "List", "")) %>%
  mutate(regionList = str_to_upper(regionList)) %>%
  mutate(regionList = factor(regionList, levels = c("OFC", "DFC", "VFC", "MFC", "M1C", "IPC", "S1C", "A1C", "ITC", "STC", "V1C", "HIP", "AMY", "THM", "STR", "CBL"))) %>%
  mutate(geneset = "Prey")

#import data for bsRNAseq, ASD102
load(file=paste0(od, "bsRNAseq_asd102_permutationData_formattedToPlot.RData"))
# format data (bsRNAseq, ASD102)
permuted =  bsRNAseq_iBAQ_Permuted_medianRanks_toPlot %>% select(starts_with("permuted")) %>% 
  apply(., 1, median)
toPlot_bsRNAseq_asd102 = bsRNAseq_iBAQ_Permuted_toPlot %>% filter(geneset == "asd102") %>%
  mutate(log10.pVal_permuted = -log10(pVal_permuted+0.00001)) %>%
  mutate(permuted.medRank = permuted,
         rankdiff = medRank - permuted.medRank) %>%
  mutate(regionList = str_replace(regionList, "List", "")) %>%
  mutate(regionList = str_to_upper(regionList)) %>%
  mutate(regionList = factor(regionList, levels = c("OFC", "DFC", "VFC", "MFC", "M1C", "IPC", "S1C", "A1C", "ITC", "STC", "V1C", "HIP", "AMY", "THM", "STR", "CBL"))) %>%
  mutate(geneset = "hcASD102")

toPlot = bind_rows(toPlot_bsRNAseq_prey100, toPlot_bsRNAseq_asd102)

###################################################################################
#Make plots and calculate Spearman's correlation coefficient
###################################################################################
# calculate correlation coefficients of median geneset rank differences for all samples in periods 2-14 for hcASD102 and ASD-PPI_100 prey
pearson_bySamples = cor(toPlot_bsRNAseq_prey100$rankdiff,toPlot_bsRNAseq_asd102$rankdiff, method="pearson")
spearman_bySamples = cor(toPlot_bsRNAseq_prey100$rankdiff,toPlot_bsRNAseq_asd102$rankdiff, method="spearman")

# calculate correlation coefficients of median rank differences across periods 2-14 for hcASD102 and ASD-PPI_100 prey
period_rankdiff_summary = toPlot %>% dplyr::group_by(geneset, period) %>% dplyr::summarize(median.rankdiff = median(rankdiff)) %>%
  pivot_wider(names_from=geneset, values_from=median.rankdiff)
pearson = cor(period_rankdiff_summary$hcASD102, period_rankdiff_summary$Prey, method="pearson")
spearman = cor(period_rankdiff_summary$hcASD102, period_rankdiff_summary$Prey, method="spearman")

median_toPlot = period_rankdiff_summary %>% pivot_longer(!period, values_to="medRankdiff", names_to="geneset")

p1 = toPlot %>% mutate(geneset = recode_factor(geneset, "hcASD102" = "hcASD102", "Prey" = "Prey (- hcASD102)")) %>%
  ggplot(aes(x = factor(period), y = rankdiff, fill = factor(period))) + geom_boxplot(show.legend = F) + facet_wrap(~geneset, scale="free") +
  theme_bw() + theme(legend.position = "none") +
  labs(title = "bsRNAseq, median geneset rank\n#Permutations = 100000",
       y = "Median geneset rank, observed - median(permuted)", 
       x = "Period") +
       # caption = paste0("Correlation of rank difference of hcASD102 and Prey across samples: Spearman = ", signif(spearman_bySamples,3), ", Pearson = ", signif(pearson_bySamples, 3),
       #                  "\nCorrelation of median(rank difference) of hcASD102 and Prey across periods: Spearman = ", signif(spearman, 3), ", Pearson = ", signif(pearson, 3))) +
  geom_smooth(method = "loess", aes(group=1), span = 1, color="black")
p1
ggsave(paste0(od, "bsRNAseq_prey100_ads102_comparison_permutationPlots_byPeriod.pdf"),
       width=7, height = 5)

# 
# #summary data for supplementary table
# toPlot %>% mutate(geneset = recode_factor(geneset, "hcASD102" = "hcASD102", "Prey" = "Prey (- hcASD102)")) %>%
#   mutate(period = factor(period)) %>%
#   group_by(period, geneset) %>% dplyr::summarize(meanRankdiff= mean(rankdiff),
#                                                  percentile_75 = quantile(rankdiff, 0.75),
#                                                  percentile_25 = quantile(rankdiff, 0.25),
#                                                  nsamples= n()) %>% 
#   arrange(geneset)


p2 = toPlot %>%
  mutate(natal = recode_factor(natal, "prenatal"="Prenatal", "postnatal"="Postnatal")) %>%
  mutate(geneset = recode_factor(geneset, "hcASD102" = "hcASD102", "Prey" = "Prey (- hcASD102)")) %>%
  ggplot(aes(x = natal, y = rankdiff, fill = natal)) + geom_boxplot(show.legend = F) + facet_wrap(~geneset, scale="free") +
  #geom_jitter(size = 0.5, alpha = 0.2) + 
  theme_bw() + theme(legend.position = "none") +
  stat_compare_means(method = "t.test") +
  labs(title = "bsRNAseq, median geneset rank\n#Permutations = 100000",
       y = "Median geneset rank, observed - median(permuted)", 
       x = "Period")
p2
ggsave(paste0(od, "bsRNAseq_prey100_ads102_comparison_permutationPlots_byNatal.pdf"),
       width=5, height = 5)

prey_rankdiff_by_sample = toPlot_bsRNAseq_prey100 %>% dplyr::select(brain:ncxRegion, rankdiff) %>% rename(prey_rankdiff = rankdiff) 
asd102_rankdiff_by_sample = toPlot_bsRNAseq_asd102 %>% dplyr::select(brain:ncxRegion, rankdiff) %>% rename(asd102_rankdiff = rankdiff)
toPlot2 = full_join(prey_rankdiff_by_sample, asd102_rankdiff_by_sample) %>% mutate(period = factor(period)) %>%
  rename(`Sample age` = natal)

p3 = ggplot(toPlot2, aes(x=asd102_rankdiff, y=prey_rankdiff)) + geom_point(alpha=0.3) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_smooth(method = "lm",color="black", formula = y ~ x) +
  theme_bw() +
  labs(x="hcASD102", y="Prey (- hcASD102)",
       title = "bsRNAseq\nMedian geneset rank, observed - median(permuted)\n#Permutations = 100000")
p3
ggsave(paste0(od, "bsRNAseq_prey100_ads102_comparison_permutationPlots_allSamples.pdf"),
       width=5, height = 5)
