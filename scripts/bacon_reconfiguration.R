# Script for running BACON

# Import libraries
library(argparse)
library(BiocParallel)
library(dplyr)
library(data.table)
library(bacon)
library(QCEWAS)
library(qqman)
library(ggplot2)
library(reshape2)
library(tibble)
library(cowplot)

# Import modified bacon functions
source("bacon_rng_fix.R")
source("bacon_init_fix.R")
source("modified_bacon_plots.R")

# Read in EWAS summary statistics
ewas <- fread("~/Downloads/jhs_Female_bmi_ewas.csv")
ewas <- ewas[, 1:6]

# Run bacon on tstatistics, effect-sizes, and standard errors
bc <- bacon(teststatistics = ewas$statistic,
             effectsizes = ewas$estimate,
             standarderrors = ewas$std.error,
             verbose = T,
             trim = 0.999,
             globalSeed = 42)

# Extract bacon-adjusted p-value
ewas$bacon.pval <- bacon::pval(bc)
# Extract bacon-adjusted t-statistic
ewas$bacon.statistic <- bacon::tstat(bc)
# Extract bacon-adjusted effect-size
ewas$bacon.es <- bacon::es(bc)
# Extract bacon-adjusted standard error
ewas$bacon.se <- bacon::se(bc)
# Estimate the original p-value lambda
ewas$lambda <- QCEWAS::P_lambda(ewas$p.value)
# Estimate the bacon-adjusted p-value lambda
ewas$b.lambda <- QCEWAS::P_lambda(ewas$bacon.pval)

# Run performance tests and export plots
ggtraces(bc) + labs(title = "filename")
ggsave("filename_bacon_traces.jpg", width = 16, height = 9.8, units = "cm")

posteriors(bc) 
fit(bc)
#dev.off()

posteriors(bc)
fit(bc)
qqman::qq(subset$p.value, main = "Before Bacon")
qqman::qq(subset$bacon.pval, main = "After Bacon")
bacon_bmi_ewas <- list()
set.seed(2974)
for (i in 1:length(bmi_ewas)){
  name <- names(bmi_ewas)[i]
  subset <- bmi_ewas[[name]]
  # Run bacon
  bacon_res <- bacon(teststatistics = subset$statistic,
                     verbose = T,
                     trim = 0.999)
  bacon_effects <- bacon(teststatistics = NULL,
                         effectsizes = subset$estimate,
                         standarderrors = subset$std.error,
                         verbose = T,
                         trim = 0.999)
  # Add bacon results to subset dataframe
  subset$bacon.pval <- bacon::pval(bacon_res)
  subset$bacon.statistic <- bacon::tstat(bacon_res)
  subset$bacon.es <- bacon::es(bacon_effects)
  subset$bacon.se <- bacon::se(bacon_effects)
  b.lambda <- QCEWAS::P_lambda(subset$bacon.pval)
  lambda <- QCEWAS::P_lambda(subset$p.value)
  subset$lambda <-lambda
  subset$b.lambda <- b.lambda
  bacon_bmi_ewas[[name]] <- subset
  cat("lambda before bacon: ", lambda, "\n", "lambda after bacon: ", b.lambda, "\n")
  # Run performance tests
  traces(bacon_res)
  posteriors(bacon_res)
  fit(bacon_effects)
  qqman::qq(subset$p.value, main = "Before Bacon")
  qqman::qq(subset$bacon.pval, main = "After Bacon")
}
