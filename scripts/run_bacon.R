# Script for running BACON bias and inflation adjustment
# Import libraries
suppressPackageStartupMessages({
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
})

# Import modified bacon functions
source("scripts/updated_bacon/bacon_rng_fix.R")
source("scripts/updated_bacon/bacon_init_fix.R")
source("scripts/updated_bacon/modified_bacon_plots.R")

# Define command line arguments
parser <- argparse::ArgumentParser(description="Script for running BACON")
parser$add_argument('--input-file', '-i',
                     required=TRUE,
                     help="Path to ewas results data")
parser$add_argument('--out-dir',
                     required=TRUE,
                     help="Path to output directory")
parser$add_argument('--out-prefix',
                     type="character",
                     nargs="?",
                     const="all",
                     default="all",
                     help="Prefix for output files")
parser$add_argument('--out-type', 
                     type="character",
                     choices=c(".csv", ".csv.gz"), 
                     nargs="?",                    
                     const=".csv",
                     default=".csv",  
                     help="Output file type: CSV or CSV.GZ")

# parse arguments
args <- parser$parse_args()
ewas_results <- args$input_file
out_dir <- args$out_dir
filename <- args$out_prefix
plotname <- gsub("_", " ", filename)
out_type <- args$out_type

# Read in EWAS summary statistics
ewas <- fread(ewas_results)
assoc <- unique(ewas$term)
if(filename=="all"){
       filename <- assoc
       plotname <- assoc
} else{
       filename <- paste0(filename, "_", assoc)
       plotname <- paste0(plotname, " ", assoc)
       cat(paste0("Running BACON for subset: ", filename, "\n"))
}


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

# Export bacon-adjusted results
fwrite(ewas, file=paste0(out_dir, filename, "_ewas_bacon_results", out_type))

# Run performance tests and export plots
traces_plot <- ggtraces(bc) + labs(title = paste0(plotname, " traces"))
ggsave(paste0(out_dir, "bacon_plots/", filename, "_traces.jpg"),
       plot = traces_plot,
       width = 16, height = 9.8, units = "cm")

posteriors_plot <- ggposteriors(bc) + labs(title = paste0(plotname, " posteriors"))
ggsave(paste0(out_dir, "bacon_plots/", filename, "_posteriors.jpg"),
       plot = posteriors_plot,
       width = 10, height = 9.8, units = "cm")

fit_plot <- ggfit(bc) + labs(title = paste0(plotname, " fit"))
ggsave(paste0(out_dir, "bacon_plots/", filename, "_fit.jpg"),
       plot = fit_plot,
       width = 12, height = 9.8, units = "cm")

qq_plot <- bacon::plot(bc, type = c("qq")) + 
       labs(title = paste0(plotname, " qq plots")) +
       theme(legend.position = "none")
ggsave(paste0(out_dir, "bacon_plots/", filename, "_qqs.jpg"),
       plot = qq_plot,
       width = 12, height = 9.8, units = "cm")

