##### ASD-PPI figure 8 statistics LMM ######
#version 3
#20260505
#modified code from BW


##### notes #####
## summary() is used instead of anova(), generated t statistics instead of F statistic as in anova function
## only summary() is applicable to data with multiple comparisons (like L327P) so used this for all
## Benjamini Hochberg correction for all
## lmer() is linear mixed effects model for all comparisons except those that do not have repeated measures
## this includes rosette levels quants for each organoid
## lm() is standard linear regression, used instead with no random effects (such as batch and rosettes), applies to FOXP4-KO and L327P DLX2 comparisons


## notes from meeting with ben, hannes, rasika, and belinda
#unique rosette ID
#consider collapsing rosettes, take avg
#outlier analysis
#remove intercept
#dont mention normality test
#can try permutation with lme

library(lme4)
library(lmerTest)

setwd("/media/chang/HDD-11/kelsey/FOXP1/R513H/fig8_stats/final/rosette_avg/")

## Fig 8B FOXP1-R513H IHC
# D39
df <- read.csv("R513H_IHC_D39.csv")
df$ID <- paste(df$differentiation, df$genotype, df$organoid, sep = "-")

markers <- c("PAX6", "BCL11B", "TBR1")
results <- list()

df$genotype <- relevel(factor(df$genotype), ref = "FOXP1 WT/WT")
df$differentiation <- factor(df$differentiation)
df$ID <- factor(df$ID)

for (m in markers) {
  f <- as.formula(paste(m, "~ genotype + (1 | differentiation)"))
  fit <- lmer(f, data = df)
  results[[m]] <- summary(fit)$coefficients
}

## unlist anova results, rbind, and save
df1 <- as.data.frame(results$PAX6)
df1$comparison <- rownames(df1)
df1$var <- "PAX6"
df2 <- as.data.frame(results$BCL11B)
df2$comparison <- rownames(df2)
df2$var <- "BCL11B"
df3 <- as.data.frame(results$TBR1)
df3$comparison <- rownames(df3)
df3$var <- "TBR1"

df <- rbind(df1, df2, df3)
rownames(df) <- NULL

df <- df[df$comparison != "(Intercept)", ]

# Calculate adjusted p-values
#toss intercept
#consider pooling experiments and BH on this
BH_padj <- p.adjust(df$`Pr(>|t|)`, method = "BH")
BH_padj <- as.data.frame(BH_padj)
df <- cbind(df, BH_padj)
df <- df[,c(6:7, 1:5, 8)]

write.csv(df, "outs/R513H_IHC_D39.csv")


## R513H D101
df <- read.csv("R513H_IHC_D101.csv")
df$ID <- paste(df$genotype, df$organoid, sep = "-")

markers <- c("PAX6", "BCL11B", "TBR1")
results <- list()

df$genotype <- relevel(factor(df$genotype), ref = "WT/WT")
#df$differentiation <- factor(df$differentiation)
df$ID <- factor(df$ID)

for (m in markers) {
  f <- as.formula(paste(m, "~ genotype"))
  fit <- lm(f, data = df)
  results[[m]] <- summary(fit)$coefficients
}

## unlist anova results, rbind, and save
df1 <- as.data.frame(results$PAX6)
df1$comparison <- rownames(df1)
df1$var <- "PAX6"
df2 <- as.data.frame(results$BCL11B)
df2$comparison <- rownames(df2)
df2$var <- "BCL11B"
df3 <- as.data.frame(results$TBR1)
df3$comparison <- rownames(df3)
df3$var <- "TBR1"

df <- rbind(df1, df2, df3)
rownames(df) <- NULL

df <- df[df$comparison != "(Intercept)", ]

# Calculate adjusted p-values
BH_padj <- p.adjust(df$`Pr(>|t|)`, method = "BH")
BH_padj <- as.data.frame(BH_padj)
df <- cbind(df, BH_padj)
df <- df[,c(5:6, 1:4, 7)]

write.csv(df, "outs/R513H_IHC_D101.csv")


## FOXP4 KO IHC
df <- read.csv("R513H_FOXP4-KO_IHC.csv")
markers <- c("BCL11B", "TBR1")
results <- list()

df$genotype <- relevel(factor(df$genotype), ref = "FOXP1-R513H_FOXP4-WT")

for (m in markers) {
  f <- as.formula(paste(m, "~ genotype"))
  fit <- lm(f, data = df)
  results[[m]] <- list(anova = summary(fit))
}

## unlist anova results, rbind, and save
df2 <- as.data.frame(results$BCL11B$anova$coefficients)
df2$comparison <- rownames(df2)
df2$var <- "BCL11B"
df3 <- as.data.frame(results$TBR1$anova$coefficients)
df3$comparison <- rownames(df3)
df3$var <- "TBR1"

df <- rbind(df2, df3)
rownames(df) <- NULL

df <- df[df$comparison != "(Intercept)", ]

# Calculate adjusted p-values
BH_padj <- p.adjust(df$`Pr(>|t|)`, method = "BH")
BH_padj <- as.data.frame(BH_padj)
df <- cbind(df, BH_padj)
df <- df[,c(5:6, 1:4, 7)]

write.csv(df, "outs/R513H_FOXP4-KO_IHC.csv")


## L327P D39
df <- read.csv("L327P_IHC_D39.csv")
df$ID <- paste(df$differentiation, df$genotype, df$organoid, sep = "-")

markers <- c("PAX6", "BCL11B", "TBR1")
results <- list()

df$genotype <- relevel(factor(df$genotype), ref = "FOXP1 WT/WT")
df$differentiation <- factor(df$differentiation)
df$ID <- factor(df$ID)

for (m in markers) {
  f <- as.formula(paste(m, "~ genotype + (1 | differentiation)"))
  fit <- lmer(f, data = df)
  results[[m]] <- summary(fit)
}

## unlist anova results, rbind, and save
df1 <- as.data.frame(results[["PAX6"]]$coefficients)
df1$comparison <- rownames(df1)
df1$var <- "PAX6"
df2 <- as.data.frame(results[["BCL11B"]]$coefficients)
df2$comparison <- rownames(df2)
df2$var <- "BCL11B"
df3 <- as.data.frame(results[["TBR1"]]$coefficients)
df3$comparison <- rownames(df3)
df3$var <- "TBR1"

df <- rbind(df1, df2, df3)
rownames(df) <- NULL

df <- df[df$comparison != "(Intercept)", ]

# Calculate adjusted p-values
BH_padj <- p.adjust(df$`Pr(>|t|)`, method = "BH")
BH_padj <- as.data.frame(BH_padj)
df <- cbind(df, BH_padj)
df <- df[,c(6:7, 1:5, 8)]

write.csv(df, "outs/L327P_IHC_D39.csv")


## L327P D101
df <- read.csv("L327P_IHC_D101.csv")
df$ID <- paste(df$genotype, df$organoid, sep = "-")

markers <- c("PAX6", "BCL11B", "TBR1")
results <- list()

df$genotype <- relevel(factor(df$genotype), ref = "FOXP1 WT/WT")
#df$differentiation <- factor(df$differentiation)
df$ID <- factor(df$ID)

for (m in markers) {
  f <- as.formula(paste(m, "~ genotype"))
  fit <- lm(f, data = df)
  results[[m]] <- summary(fit)
}

## unlist anova results, rbind, and save
df1 <- as.data.frame(results[["PAX6"]]$coefficients)
df1$comparison <- rownames(df1)
df1$var <- "PAX6"
df2 <- as.data.frame(results[["BCL11B"]]$coefficients)
df2$comparison <- rownames(df2)
df2$var <- "BCL11B"
df3 <- as.data.frame(results[["TBR1"]]$coefficients)
df3$comparison <- rownames(df3)
df3$var <- "TBR1"

df <- rbind(df1, df2, df3)
rownames(df) <- NULL

df <- df[df$comparison != "(Intercept)", ]

# Calculate adjusted p-values
BH_padj <- p.adjust(df$`Pr(>|t|)`, method = "BH")
BH_padj <- as.data.frame(BH_padj)
df <- cbind(df, BH_padj)
df <- df[,c(5:6, 1:4, 7)]

write.csv(df, "outs/L327P_IHC_D101.csv")


##### supplement #####
## DLX2
df <- read.csv("R513H_DLX2.csv")
df$ID <- paste(df$genotype, df$organoid, sep = "-")

markers <- c("DLX2")
results <- list()

df$genotype <- relevel(factor(df$genotype), ref = "FOXP1 WT/WT")
#df$differentiation <- factor(df$differentiation)
df$ID <- factor(df$ID)

for (m in markers) {
  f <- as.formula(paste(m, "~ genotype"))
  fit <- lm(f, data = df)
  results[[m]] <- summary(fit)
}

## unlist anova results, rbind, and save
df <- as.data.frame(results[["DLX2"]]$coefficients)
df$comparison <- rownames(df)
df$var <- "DLX2"

rownames(df) <- NULL

df <- df[df$comparison != "(Intercept)", ]

# Calculate adjusted p-values
BH_padj <- p.adjust(df$`Pr(>|t|)`, method = "BH")
BH_padj <- as.data.frame(BH_padj)
df <- cbind(df, BH_padj)
df <- df[,c(5:6, 1:4, 7)]

write.csv(df, "outs/R513H_DLX2.csv")


## L327P DLX2
df <- read.csv("L327P_DLX2.csv")

model <- lm(DLX2 ~ genotype, data = df)
results <- summary(model)

## unlist anova results, rbind, and save
df <- as.data.frame(results$coefficients)
df$comparison <- rownames(df)
df$var <- "DLX2"

rownames(df) <- NULL

df <- df[df$comparison != "(Intercept)", ]

# Calculate adjusted p-values
BH_padj <- p.adjust(df$`Pr(>|t|)`, method = "BH")
BH_padj <- as.data.frame(BH_padj)
df <- cbind(df, BH_padj)
df <- df[,c(5:6, 1:4, 7)]

write.csv(df, "outs/L327P_DLX2.csv")


## PPP2R5D
df <- read.csv("PPP2R5D_BCL11B.csv")
df$ID <- paste(df$genotype, df$organoid, sep = "-")

markers <- c("BCL11B")
results <- list()

df$genotype <- relevel(factor(df$genotype), ref = "PPP2R5D WT/WT")
#df$differentiation <- factor(df$differentiation)
df$ID <- factor(df$ID)

for (m in markers) {
  f <- as.formula(paste(m, "~ genotype"))
  fit <- lm(f, data = df)
  results[[m]] <- summary(fit)
}

## unlist anova results, rbind, and save
df <- as.data.frame(results[["BCL11B"]]$coefficients)
df$comparison <- rownames(df)
df$var <- "BCL11B"

rownames(df) <- NULL

df <- df[df$comparison != "(Intercept)", ]

# Calculate adjusted p-values
BH_padj <- p.adjust(df$`Pr(>|t|)`, method = "BH")
BH_padj <- as.data.frame(BH_padj)
df <- cbind(df, BH_padj)
df <- df[,c(5:6, 1:4, 7)]

write.csv(df, "outs/PPP2R5D_BCL11B.csv")


## PPP2R5D TBR1 and PAX6
df <- read.csv("PPP2R5D_TBR1_PAX6.csv")
df$ID <- paste(df$genotype, df$organoid, sep = "-")

markers <- c("TBR1", "PAX6")
results <- list()

df$genotype <- relevel(factor(df$genotype), ref = "PPP2R5D WT/WT")
#df$differentiation <- factor(df$differentiation)
df$ID <- factor(df$ID)

for (m in markers) {
  f <- as.formula(paste(m, "~ genotype"))
  fit <- lm(f, data = df)
  results[[m]] <- summary(fit)
}

## unlist anova results, rbind, and save
df1 <- as.data.frame(results[["PAX6"]]$coefficients)
df1$comparison <- rownames(df1)
df1$var <- "PAX6"
df3 <- as.data.frame(results[["TBR1"]]$coefficients)
df3$comparison <- rownames(df3)
df3$var <- "TBR1"

df <- rbind(df1, df3)
rownames(df) <- NULL

df <- df[df$comparison != "(Intercept)", ]

# Calculate adjusted p-values
BH_padj <- p.adjust(df$`Pr(>|t|)`, method = "BH")
BH_padj <- as.data.frame(BH_padj)
df <- cbind(df, BH_padj)
df <- df[,c(5:6, 1:4, 7)]

write.csv(df, "outs/PPP2R5D_TBR1_PAX6.csv")


##### rebuttal SATB2 quant #####
## R513H
df <- read.csv("R513H_SATB2.csv")
df$ID <- paste(df$genotype, df$organoid, sep = "-")

markers <- c("SATB2")
results <- list()

df$genotype <- relevel(factor(df$genotype), ref = "FOXP1 WT/WT")
#df$differentiation <- factor(df$differentiation)
df$ID <- factor(df$ID)

for (m in markers) {
  f <- as.formula(paste(m, "~ genotype"))
  fit <- lm(f, data = df)
  results[[m]] <- summary(fit)
}

## unlist anova results, rbind, and save
df <- as.data.frame(results[["SATB2"]]$coefficients)
df$comparison <- rownames(df)
df$var <- "SATB2"

rownames(df) <- NULL

df <- df[df$comparison != "(Intercept)", ]

# Calculate adjusted p-values
BH_padj <- p.adjust(df$`Pr(>|t|)`, method = "BH")
BH_padj <- as.data.frame(BH_padj)
df <- cbind(df, BH_padj)
df <- df[,c(5:6, 1:4, 7)]

write.csv(df, "outs/R513H_SATB2.csv")


## L327P SATB2
df <- read.csv("L327P_SATB2.csv")

df$ID <- paste(df$genotype, df$organoid, sep = "-")

markers <- c("SATB2")
results <- list()

df$genotype <- relevel(factor(df$genotype), ref = "FOXP1 WT/WT")
#df$differentiation <- factor(df$differentiation)
df$ID <- factor(df$ID)

for (m in markers) {
  f <- as.formula(paste(m, "~ genotype"))
  fit <- lm(f, data = df)
  results[[m]] <- summary(fit)
}

## unlist anova results, rbind, and save
df <- as.data.frame(results[["SATB2"]]$coefficients)
df$comparison <- rownames(df)
df$var <- "SATB2"

rownames(df) <- NULL
df <- df[df$comparison != "(Intercept)", ]

# Calculate adjusted p-values
BH_padj <- p.adjust(df$`Pr(>|t|)`, method = "BH")
BH_padj <- as.data.frame(BH_padj)
df <- cbind(df, BH_padj)
df <- df[,c(5:6, 1:4, 7)]

write.csv(df, "outs/L327P_SATB2.csv")


