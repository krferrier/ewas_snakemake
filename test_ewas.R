# Script for performing Epigenome-wide association analysis that takes in command line arguments
# TO DO:
#   > Add argument for which phenotype to perform the ewas on
#   > Add argument for what type of file to output results as

# Import libraries
library(R.utils)
library(argparse)
library(doFuture)
library(progressr)
library(progress)
library(foreach)
library(vars)
library(dplyr)
library(data.table)
library(tibble)
library(tictoc)

# Source functions from outside scripts
source("stratify.R")
source("chunk.R")

################################################################################################################################
#                                   DEFINE AND PARSE COMMAND LINE ARGUMENTS                                                    # 
################################################################################################################################
# Define command line arguments
parser <- argparse::ArgumentParser(description = "Script for running EWAS")
parser$add_argument("--pheno", required=TRUE, help = "Path to phenotype data [samples, phenotypes] in CSV or CSV.GZ format with the first column being sample identifiers.")
parser$add_argument("--methyl", required=TRUE, help = "Path to methylation data [samples, CpGs] in CSV or CSV.GZ format with the first column being sample identifiers.")
parser$add_argument("--stratify", type = "character", nargs="*", help = "Variables to stratify")
parser$add_argument("--chunk-size", type = "integer", nargs="?", const = 1000, help = "number of CpGs per chunk")
parser$add_argument('--processing-type', '-pt',
                    type = "character",
                    nargs="?",
                    const="sequential",
                    default="sequential",
                    choices = c("sequential", "multisession", "multicore", "cluster"),
                    help="Parallelization type: sequential (default), multisession, multicore, or cluster")
parser$add_argument('--workers', type = "integer", nargs="?", const=1, default=1, help="Number of processes to run in parallel")
parser$add_argument('--out-dir', type = "character", required=TRUE, help="Path to output directory")

# parse arguments
args <- parser$parse_args()
pheno <- args$pheno
mvals <- args$methyl
stratify_vars <- args$stratify
chunk_size <- args$chunk_size
pt <- args$processing_type
n_workers <- args$workers
out_dir <- args$out_dir
print(out_dir)
# ################################################################################################################################
# #                                   CHECK THAT INPUT FILES ARE CSV OR CSV.GZ                                                   # 
# ################################################################################################################################
# Function to check for the file extensions '.csv' and '.csv.gz'
check_input <- function(f){
  if(endsWith(f, c(".csv")) | endsWith(f, c(".csv.gz"))){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# Throw an error if the phenotype data is not the correct format
if(check_input(pheno) == FALSE){
  stop("Phenotype data is not in CSV or CSV.GZ format. \n")
}
# Throw an error if the methylation data is not the correct format
if(check_input(mvals) == FALSE){
  stop("Methylation data is not in CSV or CSV.GZ format. \n")
}

####################################################################################################
#                                   READ IN DATA                                                   #
####################################################################################################

# Read in phenotype data
pheno <- fread(pheno) %>% column_to_rownames(var=colnames(.)[1])

# Read in methylation data
mvals <- fread(mvals) %>% column_to_rownames(var=colnames(.)[1])


####################################################################################################
#                                   STRATIFY DATA                                                  #
####################################################################################################

# Make a key for subsetting the data based on the --stratify variables given. If none were given, the subset.key is a list of all samples. 
subset.key <- make_subset.key(pheno, stratify_vars)

# Stratify the phenotype and methylation data with the subset.key. If no stratification variables were given, the data is not stratified. 
pheno <- stratify.pheno(pheno, stratify_vars)
mvals <- stratify.mvals(mvals)

####################################################################################################
#                                      CHUNK DATA                                                  #
####################################################################################################

# Create a list of cpg.chunks from a vector of all the CpGs to be tested, where each cpg.chunk is a list with n=chunk_size CpGs.
# If the total number of CpGs is not perfectly divisible by chunk_size, then the last cpg.chunk will have a length equal to the remainder. 
cpgs <- cpg.chunks(chunk_size, colnames(mvals[[1]]))


################################################################################################################################
#                                       LINEAR REGRESSION ANALYSIS                                                             # 
################################################################################################################################
# Set paramaters for processing the data sequentially (default) or in parallel.
registerDoFuture()
n.workers = n_workers
if(pt=="sequential"){
  plan(strategy = pt)
  cat("Processing run sequentially. \n")
}else{
  plan(strategy = pt, workers = n.workers)
  cat("Asynchronous parallel processing using", pt, "with ", n.workers, " worker(s). \n")
}

# Create a list for the results
results <- list()

# Loop through each strata
for(j in names(subset.key)){
  # Loop through each strata
  cat("Starting Subset: ", j, "\n")
  tic()
  # Create progress bar
  p <- progress::progress_bar$new(format = "[:bar] :percent ELAPSED::elapsed, ETA::eta",
                                  total = length(cpgs))
  p.sub <- pheno[[j]]
  m.sub <- chunk.df(mvals[[j]],cpgs)
  res.chunks <- list()
  # Loop through chunks
  for(ii in 1:length(cpgs)){
    m.chunk <- m.sub[[ii]]
    # Run linear regression
    ewas <- foreach(i = cpgs[[ii]], .packages = c("vars"), .verbose = F) %dopar% {
      .GlobalEnv$m.chunk <- m.chunk
      .GlobalEnv$p.sub <- p.sub
      model.base <- paste(colnames(p.sub), collapse = " + ")
      #string.formula <- paste0(i, " ~ ", model.base, " + n.haps")
      string.formula <- paste0("m.chunk$", i, " ~ ", model.base)
      fit <- lm(formula = string.formula, data = p.sub) %>%
        broom::tidy(conf.int = T) %>%
        dplyr::filter(term %in% c('BMI')) %>%
        dplyr::mutate(cpgid = i)
      fit
    }
    ewas <- rbindlist(ewas)
    res.chunks[[ii]] <- ewas
    p$tick()
  }
  res.chunks <- rbindlist(res.chunks)
  results[[j]] <- res.chunks
  toc()
}
save(results, file="test_ewas.RData")

# Export results by strata (if analysis was stratified)
for (i in names(subset.key)){
  sub.results <- results[[i]]
  file.name <- paste0(out_dir, i, "_ewas_results.csv")
  fwrite(sub.results, file = file.name)
}