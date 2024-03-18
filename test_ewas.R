# update script for performing linear regressions to take in command line arguments

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

################################################################################################################################
#                                 DUMMY DATA FOR TESTING                                                                       # 
################################################################################################################################
# dummy methylation data for testing
num_rows <- 10
num_cols <- 10000
column_names <- paste0('cg', 1:num_cols)
row_names <- paste0('sample_', 1:num_rows)

mvals <- data.frame(matrix(nrow = num_rows, ncol = num_cols))
rownames(mvals) <- row_names
colnames(mvals) <- column_names

for (i in 1:num_rows) {
  for (j in 1:num_cols) {
    mvals[i, j] <- runif(1, min = -5, max = 5)
  }
}

# dummy phenotype data for testing
pheno_all <- data.frame("sampleID" = c(paste0("sample_", 1:10)),
                        "sex" = c("M", "F", "M", "M", "F", "F", "F", "M", "F", "F"),
                        "re" = c(1,2,3,1,2,3,1,2,3,1),
                        "BMI" = c(21, 23, 26, 19, 25, 29, 22, 31, 24, 32))
pheno_all <- pheno_all %>% column_to_rownames("sampleID")

################################################################################################################################
#                                               STRATIFY                                                                       # 
################################################################################################################################
# Define command line arguments
parser <- argparse::ArgumentParser(description = "Stratify data")
parser$add_argument("--stratify", nargs="*", help = "Variables to stratify")
args <- parser$parse_args()

# Convert the stratify argument to character vector
stratify_vars <- args$stratify
stratify_vars <- "sex"

# Function to create a subset.key for stratifying analysis. If no stratify variables are given, the subset.key will be a list of 1 dataframe with all samples.
# Check if variables were provided for stratification
make_subset.key <- function(data, stratify_vars) {
  if (length(stratify_vars) != 0) {
    # Check if stratify_vars exist as column names in the dataframe
    missing_vars <- setdiff(stratify_vars, colnames(data))
    if (length(missing_vars) > 0) {
      stop(paste("The following stratify variable(s) do not exist in the phenotype data:", paste(missing_vars, collapse = ", ")))
    }
    # Perform splitting based on the provided variables
    data <- data %>%
      dplyr::select(all_of(stratify_vars))
    
    split_function <- function(data) {
      group_vars <- do.call(paste, c(data[stratify_vars], sep = "_"))
      split(data, group_vars)
    }
    data <- split_function(data)
  } else {
    # If no stratify variables provided, store the data as a single element list
    data <- list(data)
    names(data) <- c("all")
  }
  # Extract row names for each subset
  data <- lapply(data, function(i) {
    ids <- rownames(i)
    return(ids)
  })
  return(data)
}

subset.key <- make_subset.key(pheno_all, stratify_vars)
subset.key

# Stratify other datasets with the subset.key
stratify.pheno <- function(pheno.df, stratify_vars){
  pheno <- lapply(subset.key, function(i){
    df.sub = pheno.df[rownames(pheno.df) %in% i,]
    df.sub <- df.sub[order(rownames(df.sub)),]
    if(length(stratify_vars) != 0){
      df.sub <- df.sub %>% dplyr::select(-all_of(stratify_vars))
    }
    return(df.sub)
  })
  names(pheno) <- names(subset.key)
  return(pheno)
}

pheno_strat <- stratify.pheno(pheno_all, stratify_vars)
pheno_strat

stratify.mvals <- function(mvals.df){
  mvals <- lapply(subset.key, function(i){
    df.sub <- mvals.df[rownames(mvals.df) %in% i,]
    df.sub <- df.sub[order(rownames(df.sub)), order(colnames(df.sub))]
    return(df.sub)
  })
  names(mvals) <- names(subset.key)
  return(mvals)
}

mvals_strat <- stratify.mvals(mvals)
mvals_strat

################################################################################################################################
#                                               CHUNK                                                                          # 
################################################################################################################################
# Define command line arguments
parser <- argparse::ArgumentParser(description = "Chunk methylation data")
parser$add_argument("--chunk-size", type = "integer", nargs="?", const = 1000, help = "number of CpGs per chunk")
args <- parser$parse_args()

# parse the stratify argument to integer
chunk_size <- args$chunk_size
chunk_size <- 100

# vector of the cpgs 
cpgs <- colnames(mvals)

cpg.chunks <- function(n, cpgs){
  if(n<=0){
    stop("The number of CpGs specified in --chunk-size must be an integer greater than 0")
  }
  if(n>length(cpgs)){
    warning("The number of CpGs specified in --chunk-size (default=1000 if no value was given) is larger than the total number of CpGs being tested. No chunking will be performed.")
  }
  # Set parameters for creating chunks of n cpgs. This will create a list comprised of n-sized subsets and a subset with the remainder that was not divisible by n.
  n <- n # number of cpgs per subset
  nr <- length(cpgs) # total number of cpgs 
  # function to create a list of equal-sized chunks + one chunk with the remainder
  chunk <- rep(seq_len(ceiling(nr/n)),each = n,length.out = nr) 
  
  # Chunk the cpgs into n-sized lists
  chunks <- split(cpgs, f = chunk)
  return(chunks)
}

# chunk vector of cpgs to use for splitting methylation data
cpgs <- cpg.chunks(chunk_size, cpgs)

chunk.df <- function(df, cpgs){
  df <- lapply(cpgs, function(x){
    df.sub <- df[,colnames(df) %in% x]
    return(df.sub)
  })
  return(df)
}

chunked_df <- chunk.df(mvals, cpgs)
chunked_df

################################################################################################################################
#                                               LINEAR REGRESSIONS                                                             # 
################################################################################################################################
# Define command line arguments
parser <- argparse::ArgumentParser(description = "Perform linear regression analysis")
parser$add_argument('--parallel-type', '-pt',
                    type = "character",
                    nargs="?",
                    const="sequential",
                    default="sequential",
                    choices = c("sequential", "multisession", "multicore", "cluster"),
                    help="Parallelization type: sequential (default), multisession, multicore, or cluster")
parser$add_argument('--workers', type = "integer", nargs="?", const=1, default=1, help="Number of processes to run in parallel")
args <- parser$parse_args()

# parse the stratify argument to integer
pt <- args$parallel_type
n_workers <- args$workers
pt <- "multisession"
n_workers <- 2

pheno_all <- stratify.pheno(pheno_all, stratify_vars)
mvals <- stratify.mvals(mvals)
results <- list()

# Set paramaters for parallelizing
registerDoFuture()
n.workers = n_workers
if(pt=="sequential"){
  plan(strategy = pt)
  cat("Processing run sequentially. \n")
}else{
  plan(strategy = pt, workers = n.workers)
  cat("Asynchronous parallel processing using", pt, "with ", n.workers, " worker(s). \n")
}


# Loop through each strata
for(j in names(subset.key)){
  # Loop through each strata
  cat("Starting Subset: ", j, "\n")
  tic()
  # Create progress bar
  p <- progress::progress_bar$new(format = "[:bar] :percent ELAPSED::elapsed, ETA::eta",
                                  total = length(cpgs))
  p.sub <- pheno_all[[j]]
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