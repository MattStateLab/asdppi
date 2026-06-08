##### ASD-PPI figure 8 statistics LMM + permutation tests ######
#version 4
#20260506
#modified code from BW


##### notes #####
## This script uses organoid-level data.
## Rosettes were averaged first, so each row is one organoid.
## For D39 datasets, genotype is a fixed effect and differentiation/batch is a random effect.
## For datasets without a batch column, genotype is tested with standard linear regression.
## Intercepts are removed before Benjamini Hochberg correction.
## Permutation tests shuffle genotype labels across organoids.
## For D39 datasets, genotype labels are shuffled within each differentiation/batch.

## notes
# boundary (singular) fit: see help('isSingular')
# for R513H and L327P BCL11B D39


library(lme4)
library(lmerTest)

## Set this to the folder with the organoid-averaged CSV files.
## Example:
setwd("/media/chang/HDD-11/kelsey/FOXP1/R513H/fig8_stats/final/rosette_avg/")

n_perm <- 10000
set.seed(20260506)


##### helper functions #####

make_clean_output <- function(df, reference_genotype) {
  ## Keep only genotype effects. The intercept is not a genotype comparison.
  df <- df[df$comparison != "(Intercept)", ]

  ## Add explicit comparison columns so the output table is easier to read.
  df$test_genotype <- sub("^genotype", "", df$comparison)
  df$reference_genotype <- reference_genotype
  df$contrast <- paste(df$test_genotype, "vs", df$reference_genotype)

  df
}

get_t_values <- function(df) {
  ## lmer output uses "t value"; lm output can use "t value" too.
  ## This helper keeps the permutation code below simple.
  if ("t value" %in% names(df)) {
    return(df$`t value`)
  }
  if ("t-value" %in% names(df)) {
    return(df$`t-value`)
  }
  stop("Could not find t value column.")
}

shuffle_genotypes <- function(df, batch_col = NULL) {
  ## This is the key permutation step.
  ## We do not shuffle marker values.
  ## We do not shuffle individual rosettes.
  ## We only shuffle genotype labels across organoids.
  ##
  ## If batch_col is provided, genotypes are shuffled within each batch.
  ## This preserves the number of WT and mutant organoids in each batch.

  df$genotype <- as.character(df$genotype)

  if (!is.null(batch_col)) {
    batch_levels <- unique(df[[batch_col]])
    for (b in batch_levels) {
      batch_rows <- which(df[[batch_col]] == b)
      df$genotype[batch_rows] <- sample(df$genotype[batch_rows])
    }
  } else {
    df$genotype <- sample(df$genotype)
  }

  df
}

fit_one_marker <- function(df,
                           marker,
                           reference_genotype,
                           model_type,
                           batch_col = NULL) {
  ## Drop organoids missing this marker.
  keep_cols <- c(marker, "genotype")
  if (!is.null(batch_col)) {
    keep_cols <- c(keep_cols, batch_col)
  }
  df <- df[complete.cases(df[, keep_cols]), ]

  df$genotype <- relevel(factor(df$genotype), ref = reference_genotype)
  if (!is.null(batch_col)) {
    df[[batch_col]] <- factor(df[[batch_col]])
  }

  if (model_type == "lmer") {
    f <- as.formula(paste(marker, "~ genotype + (1 |", batch_col, ")"))
    fit <- lmer(f, data = df)
    out <- as.data.frame(summary(fit)$coefficients)
  } else {
    f <- as.formula(paste(marker, "~ genotype"))
    fit <- lm(f, data = df)
    out <- as.data.frame(summary(fit)$coefficients)
    out$df <- df.residual(fit)
  }

  out$comparison <- rownames(out)
  out$var <- marker
  rownames(out) <- NULL

  out <- make_clean_output(out, reference_genotype)

  list(
    data = df,
    fit = fit,
    output = out
  )
}

run_permutation_test <- function(df,
                                 marker,
                                 reference_genotype,
                                 model_type,
                                 batch_col = NULL,
                                 n_perm = 10000) {
  ## First fit the real observed model.
  observed <- fit_one_marker(
    df = df,
    marker = marker,
    reference_genotype = reference_genotype,
    model_type = model_type,
    batch_col = batch_col
  )

  observed_output <- observed$output
  observed_terms <- observed_output$comparison
  observed_t <- get_t_values(observed_output)

  ## This matrix stores the genotype t statistic from each shuffled dataset.
  perm_t <- matrix(NA_real_, nrow = n_perm, ncol = length(observed_terms))
  colnames(perm_t) <- observed_terms

  for (i in seq_len(n_perm)) {
    ## One permutation:
    ## 1. Copy the observed organoid-level data.
    ## 2. Shuffle genotype labels across organoids.
    ## 3. Refit the same model.
    ## 4. Store the genotype t statistic.
    perm_df <- shuffle_genotypes(observed$data, batch_col = batch_col)

    perm_fit <- try(
      fit_one_marker(
        df = perm_df,
        marker = marker,
        reference_genotype = reference_genotype,
        model_type = model_type,
        batch_col = batch_col
      ),
      silent = TRUE
    )

    if (inherits(perm_fit, "try-error")) {
      next
    }

    perm_output <- perm_fit$output
    perm_t_values <- get_t_values(perm_output)
    names(perm_t_values) <- perm_output$comparison

    for (term in observed_terms) {
      if (term %in% names(perm_t_values)) {
        perm_t[i, term] <- perm_t_values[[term]]
      }
    }
  }

  ## Permutation p-value:
  ## how often was the shuffled absolute t statistic at least as large as observed?
  observed_output$permutation_p <- vapply(seq_along(observed_terms), function(j) {
    null_t <- perm_t[, j]
    null_t <- null_t[is.finite(null_t)]
    if (length(null_t) == 0) {
      return(NA_real_)
    }
    (sum(abs(null_t) >= abs(observed_t[[j]])) + 1) / (length(null_t) + 1)
  }, numeric(1))

  observed_output$n_successful_permutations <- colSums(is.finite(perm_t))
  observed_output
}

run_marker_panel <- function(input_csv,
                             output_csv,
                             markers,
                             reference_genotype,
                             model_type,
                             batch_col = NULL,
                             n_perm = 10000) {
  df <- read.csv(input_csv, check.names = FALSE)
  names(df) <- sub("^\ufeff", "", names(df))
  names(df) <- sub("^\357\273\277", "", names(df))

  results <- list()

  for (m in markers) {
    message("Running ", input_csv, " / ", m)
    results[[m]] <- run_permutation_test(
      df = df,
      marker = m,
      reference_genotype = reference_genotype,
      model_type = model_type,
      batch_col = batch_col,
      n_perm = n_perm
    )
  }

  df_out <- do.call(rbind, results)
  rownames(df_out) <- NULL

  ## BH corrections across all genotype tests in this panel.
  df_out$BH_padj <- p.adjust(df_out$`Pr(>|t|)`, method = "BH")
  df_out$permutation_BH <- p.adjust(df_out$permutation_p, method = "BH")
  df_out$model_significant_FDR_0.05 <- df_out$BH_padj <= 0.05
  df_out$permutation_significant_FDR_0.05 <- df_out$permutation_BH <= 0.05

  ## Put the human-readable columns first.
  first_cols <- c(
    "contrast",
    "test_genotype",
    "reference_genotype",
    "var",
    "Estimate",
    "Std. Error",
    "df",
    "t value",
    "Pr(>|t|)",
    "BH_padj",
    "model_significant_FDR_0.05",
    "permutation_p",
    "permutation_BH",
    "permutation_significant_FDR_0.05",
    "n_successful_permutations",
    "comparison"
  )
  first_cols <- first_cols[first_cols %in% names(df_out)]
  df_out <- df_out[, c(first_cols, setdiff(names(df_out), first_cols))]

  if (!dir.exists(dirname(output_csv))) {
    dir.create(dirname(output_csv), recursive = TRUE)
  }
  write.csv(df_out, output_csv, row.names = FALSE)
  df_out
}


##### Fig 8B FOXP1-R513H IHC #####

## D39
## Model: marker ~ genotype + (1 | differentiation)
## Permutation: shuffle genotype labels within differentiation.
R513H_D39 <- run_marker_panel(
  input_csv = "R513H_IHC_D39.csv",
  output_csv = "outs/R513H_IHC_D39_model_permutation.csv",
  markers = c("PAX6", "BCL11B", "TBR1"),
  reference_genotype = "FOXP1 WT/WT",
  model_type = "lmer",
  batch_col = "differentiation",
  n_perm = n_perm
)

## D101
## Model: marker ~ genotype
## Permutation: shuffle genotype labels across organoids.
R513H_D101 <- run_marker_panel(
  input_csv = "R513H_IHC_D101.csv",
  output_csv = "outs/R513H_IHC_D101_model_permutation.csv",
  markers = c("PAX6", "BCL11B", "TBR1"),
  reference_genotype = "WT/WT",
  model_type = "lm",
  n_perm = n_perm
)


##### FOXP4 KO IHC #####

FOXP4_KO <- run_marker_panel(
  input_csv = "R513H_FOXP4-KO_IHC.csv",
  output_csv = "outs/R513H_FOXP4-KO_IHC_model_permutation.csv",
  markers = c("BCL11B", "TBR1"),
  reference_genotype = "FOXP1-R513H_FOXP4-WT",
  model_type = "lm",
  n_perm = n_perm
)


##### L327P IHC #####

## D39
## Model: marker ~ genotype + (1 | differentiation)
## Permutation: shuffle genotype labels within differentiation.
L327P_D39 <- run_marker_panel(
  input_csv = "L327P_IHC_D39.csv",
  output_csv = "outs/L327P_IHC_D39_model_permutation.csv",
  markers = c("PAX6", "BCL11B", "TBR1"),
  reference_genotype = "FOXP1 WT/WT",
  model_type = "lmer",
  batch_col = "differentiation",
  n_perm = n_perm
)

## error for TBR1:
## Warning messages:
#1: Model failed to converge with 1 negative eigenvalue: -1.3e+00 
#2: Model failed to converge with 1 negative eigenvalue: -7.9e-02 
#3: Model failed to converge with 1 negative eigenvalue: -4.0e-01 
#4: Model failed to converge with 1 negative eigenvalue: -4.6e-01 
#5: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model failed to converge with max|grad| = 0.00540803 (tol = 0.002, component 1)
#                6: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                  Model failed to converge with max|grad| = 0.00540836 (tol = 0.002, component 1)

## D101
L327P_D101 <- run_marker_panel(
  input_csv = "L327P_IHC_D101.csv",
  output_csv = "outs/L327P_IHC_D101_model_permutation.csv",
  markers = c("PAX6", "BCL11B", "TBR1"),
  reference_genotype = "FOXP1 WT/WT",
  model_type = "lm",
  n_perm = n_perm
)


##### supplement #####

## DLX2
R513H_DLX2 <- run_marker_panel(
  input_csv = "R513H_DLX2.csv",
  output_csv = "outs/R513H_DLX2_model_permutation.csv",
  markers = c("DLX2"),
  reference_genotype = "FOXP1 WT/WT",
  model_type = "lm",
  n_perm = n_perm
)

L327P_DLX2 <- run_marker_panel(
  input_csv = "L327P_DLX2.csv",
  output_csv = "outs/L327P_DLX2_model_permutation.csv",
  markers = c("DLX2"),
  reference_genotype = "FOXP1 WT/WT",
  model_type = "lm",
  n_perm = n_perm
)

## PPP2R5D
PPP2R5D_BCL11B <- run_marker_panel(
  input_csv = "PPP2R5D_BCL11B.csv",
  output_csv = "outs/PPP2R5D_BCL11B_model_permutation.csv",
  markers = c("BCL11B"),
  reference_genotype = "PPP2R5D WT/WT",
  model_type = "lm",
  n_perm = n_perm
)

PPP2R5D_TBR1_PAX6 <- run_marker_panel(
  input_csv = "PPP2R5D_TBR1_PAX6.csv",
  output_csv = "outs/PPP2R5D_TBR1_PAX6_model_permutation.csv",
  markers = c("TBR1", "PAX6"),
  reference_genotype = "PPP2R5D WT/WT",
  model_type = "lm",
  n_perm = n_perm
)


##### rebuttal SATB2 quant #####

R513H_SATB2 <- run_marker_panel(
  input_csv = "R513H_SATB2.csv",
  output_csv = "outs/R513H_SATB2_model_permutation.csv",
  markers = c("SATB2"),
  reference_genotype = "FOXP1 WT/WT",
  model_type = "lm",
  n_perm = n_perm
)

L327P_SATB2 <- run_marker_panel(
  input_csv = "L327P_SATB2.csv",
  output_csv = "outs/L327P_SATB2_model_permutation.csv",
  markers = c("SATB2"),
  reference_genotype = "FOXP1 WT/WT",
  model_type = "lm",
  n_perm = n_perm
)
