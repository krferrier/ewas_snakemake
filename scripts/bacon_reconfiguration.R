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
ewas <- fread("~/Lange_Lab/Data/JHS/Methylation/bmi_ewas/jhs_f_bmi_ewas.csv")
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

ggposteriors(bc) + labs(title = "filename")
ggsave("filename_bacon_posteriors.jpg", width = 16, height = 9.8, units = "cm")

ggfit(bc) + labs(title = "filename")
ggsave("filename_bacon_fit.jpg", width = 16, height = 9.8, units = "cm")

qqman::qq(subset$p.value, main = "Before Bacon")
qqman::qq(subset$bacon.pval, main = "After Bacon")
