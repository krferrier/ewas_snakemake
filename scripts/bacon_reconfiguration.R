# Script for running BACON

# Import libraries
library(argparse)
library(BiocParallel)
library(dplyr)
library(data.table)
library(bacon)

source("Claw Lab/ewas_snakemake/scripts/bacon_rng_fix.R")

set.seed(12345)
biases <- runif(6, -0.2, 0.2)
inflations <- runif(6, 1, 1.3)
es <- matrix(nrow=5000, ncol=6)
for(i in 1:6)
  es[,i] <- rnormmix(5000, c(0.9, biases[i], inflations[i], 0, 4, 1), shuffle=FALSE)
se <- replicate(6, 0.8*sqrt(4/rchisq(5000,df=4)))
colnames(es) <- colnames(se) <- LETTERS[1:ncol(se)]
rownames(es) <- rownames(se) <- 1:5000
head(rownames(es))

# run bacon 
bc1 <- bacon(teststatistics = NULL,
             effectsizes = es,
             standarderrors = se,
             verbose = T,
             trim = 0.999,
             rng = 432)
# run bacon again on the same exact data
bc2 <- bacon(teststatistics = NULL,
             effectsizes = es,
             standarderrors = se,
             verbose = T,
             trim = 0.999,
             rng = 432)
# check if outputs are the same
test.1 <- list(inflation(bc1), bias(bc1), pval(bc1))
test.2 <- list(inflation(bc1), bias(bc2), pval(bc2))
identical(test.1[[1]], test.2[[1]])
identical(test.1[[2]], test.2[[2]])
identical(test.1[[3]], test.2[[3]])


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
