# Script for performing Epigenome-Wide association analysis that takes in command line arguments

# Import libraries
suppressPackageStartupMessages({
    library(R.utils)
    library(argparse)
    library(doFuture)
    library(progressr)
    library(progress)
    library(foreach)
    library(vars)
    library(dplyr)
    library(data.table)
    library(fst)
    library(tibble)
    library(tictoc)
    library(cli)
})

# Source functions from outside scripts
source("scripts/fxns/chunk_fxns.R")

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
parser$add_argument('--stratified',
                    choices=c("yes", "no"), 
                    default="no",
                    help="Stratified analysis: yes or no")
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
parser$add_argument('--out-prefix',
                    type="character",
                    nargs="?",
                    const="all",
                    default="all",
                    help="Prefix for output files")                       

# parse arguments
args <- parser$parse_args()
pheno <- args$pheno
mvals <- args$methyl
assoc_var <- args$assoc
stratified <- args$stratified
chunk_size <- args$chunk_size
pt <- args$processing_type
n_workers <- args$workers
out_dir <- args$out_dir
out_type <- args$out_type
out_prefix <- args$out_prefix

# Set number of threads to use if using fst to read/write data
threads_fst(nr_of_threads = n_workers)

####################################################################################################
#                                   READ IN DATA                                                   #
####################################################################################################
# Read in the phenotype data
if(endsWith(pheno, '.fst')){
  pheno <- read_fst(pheno)  %>% 
    column_to_rownames(var=colnames(.)[1]) # Move the sample IDs to the rownames
} else {
  pheno <- fread(pheno) %>% 
    column_to_rownames(var=colnames(.)[1]) # Move the sample IDs to the rownames
}

# Read in methylation data
if(endsWith(mvals, '.fst')){
  mvals <- read_fst(mvals)  %>% 
    column_to_rownames(var=colnames(.)[1]) # Move the sample IDs to the rownames
} else {
  mvals <- fread(mvals) %>% 
    column_to_rownames(var=colnames(.)[1]) # Move the sample IDs to the rownames
}

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
#                                      CHUNK DATA                                                  #
####################################################################################################

# Create a list of cpg.chunks from a vector of all the CpGs to be tested, where each cpg.chunk is a list with n=chunk_size CpGs.
# If the total number of CpGs is not perfectly divisible by chunk_size, then the last cpg.chunk will have a length equal to the remainder. 
cpgs <- cpg.chunks(chunk_size, colnames(mvals))
mvals <- chunk.df(mvals, cpgs)


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
    options(future.globals.maxSize= +Inf)
    cat("Asynchronous parallel processing using", pt, "with", n.workers, "worker(s). \n")
}

cat("Starting EWAS: ", out_prefix, "\n")

# Start timer 
tic()
# Create progress bar
handlers(global = TRUE)
if (interactive() == TRUE) {
  handlers(list(
    handler_progress(format = "[:bar] :percent ELAPSED::elapsed, ETA::eta")
  ))
} else {
  handlers("cli")
  options(cli.progress_handlers = "progressr")
}

ewas <- function(mvals, pheno){
  tic()
  p <- progressor(along = 1:length(mvals))
  # Loop through chunks
  results <- foreach(ii = 1:length(mvals), .verbose = F, .combine = "rbind", .packages = c("vars"), .inorder = T) %dopar% {
    m.chunk <- mvals[[ii]]    
    p()
    # Loop through cpgs in a chunk
    foreach(i = colnames(m.chunk), .verbose = F, .combine= "rbind", .packages = c("vars"), .inorder = F) %dopar% {
      .GlobalEnv$m.chunk <- m.chunk
      model.base <- paste(colnames(pheno), collapse = " + ")
      string.formula <- paste0("m.chunk$", i, " ~ ", model.base)
      fit <- lm(formula = string.formula, data = pheno) %>%
        broom::tidy(conf.int = F) %>%
        dplyr::filter(term %in% c(assoc_var)) %>%
        dplyr::mutate(cpgid = i)
      fit
    }
  }
  toc()
  return(results)
}

results <- ewas(mvals, pheno)

# Export results 
results$n <- nrow(pheno)
if (stratified == "yes"){
    filename <- paste0(out_dir, out_prefix, "_", assoc_var, "_ewas_results", out_type)
} else{
    filename <- paste0(out_dir, assoc_var, "_ewas_results", out_type)
}

fwrite(results, file = filename)

