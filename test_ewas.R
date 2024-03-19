# Script for performing Epigenome-Wide association analysis that takes in command line arguments

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
parser <- argparse::ArgumentParser(description="Script for running EWAS")
parser$add_argument("--pheno",
                    required=TRUE,
                    help="Path or URL to phenotype data [samples, phenotypes] with the first column being sample identifiers. \n
                    Acceptable formates include CSV, Jay, XLSX, and plain text. \n
                    Data can also be inside an archive such as .tar, .gz, .zip, .gz2, or .tgz .")
parser$add_argument("--methyl",
                    required=TRUE, 
                    help="Path or URL to methylation data [samples, CpGs] with the first column being sample identifiers. \n
                    Acceptable formates include CSV, Jay, XLSX, and plain text. \n
                    Data can also be inside an archive such as .tar, .gz, .zip, .gz2, or .tgz .")
parser$add_argument("--assoc", 
                    required=TRUE,
                    type="character", 
                    nargs=1, 
                    help="Variable to perform association with.")
parser$add_argument("--stratify", 
                    type="character", 
                    nargs="*", 
                    help="Variable(s) to stratify.")
parser$add_argument("--chunk-size", 
                    type="integer", 
                    nargs="?", 
                    const=1000, 
                    default=1000, 
                    help="number of CpGs per chunk.")
parser$add_argument('--processing-type', '-pt',
                    type="character",
                    nargs="?",
                    const="sequential",
                    default="sequential",
                    choices=c("sequential", "multisession", "multicore", "cluster"),
                    help="Parallelization type: sequential (default), multisession, multicore, or cluster")
parser$add_argument('--workers', 
                    type="integer", 
                    nargs="?", 
                    const=1, 
                    default=1, 
                    help="Number of processes to run in parallel")
parser$add_argument('--out-dir', 
                    type="character",
                    nargs="?",                    
                    const="~/", 
                    default="~/",  
                    help="Path to output directory")
parser$add_argument('--out-type', 
                    type="character",
                    choices=c(".csv", ".csv.gz"), 
                    nargs="?",                    
                    const=".csv",
                    default=".csv",  
                    help="Output file type: CSV or CSV.GZ")                    

# parse arguments
args <- parser$parse_args()
pheno <- args$pheno
mvals <- args$methyl
assoc_var <- args$assoc
stratify_vars <- args$stratify
chunk_size <- args$chunk_size
pt <- args$processing_type
n_workers <- args$workers
out_dir <- args$out_dir
out_type <- args$out_type

####################################################################################################
#                                   READ IN DATA                                                   #
####################################################################################################

# Read in phenotype data
pheno <- fread(pheno) %>% 
  column_to_rownames(var=colnames(.)[1]) # Move the sample IDs to the rownames

# Read in methylation data
mvals <- fread(mvals) %>% 
column_to_rownames(var=colnames(.)[1]) # Move the sample IDs to the rownames


####################################################################################################
#                                     CHECK DATA                                                   #
####################################################################################################

# Check that the association variable exists in the phenotype data
missing_vars <- setdiff(assoc_var, colnames(pheno))
if (length(missing_vars) > 0) {
  stop(paste("The following association variable(s) do not exist in the phenotype data:", paste(missing_vars, collapse = ", ")))
  # If the variable does exist in the phenotype data, set the variable as the first column
}else{
  pheno <- pheno %>% relocate(all_of(assoc_var))
}

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
  cat("Asynchronous parallel processing using", pt, "with", n.workers, "worker(s). \n")
}

# Create a list for the results
results <- list()

# Loop through each strata. If data was not stratified, the loop will run once with all samples.
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
      string.formula <- paste0("m.chunk$", i, " ~ ", model.base)
      fit <- lm(formula = string.formula, data = p.sub) %>%
        broom::tidy(conf.int = T) %>%
        dplyr::filter(term %in% c(assoc_var)) %>%
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

# Export results by strata (if analysis was stratified)
for (i in names(subset.key)){
  sub.results <- results[[i]]
  file.name <- paste0(out_dir, i, "_", assoc_var, "_ewas_results", out_type)
  fwrite(sub.results, file = file.name)
}