# Test argparsing for stratification of a dataframe

library(argparse)
library(dplyr)
library(tibble)

# Define command line arguments
parser <- argparse::ArgumentParser(description = "Stratify data")
parser$add_argument("--stratify", nargs="*", help = "Variables to stratify")
args <- parser$parse_args()

# Convert the stratify argument to character vector
stratify_vars <- args$stratify
stratify_vars 

# Create phenotype data for testing
pheno_all <- data.frame("sampleID" = c(paste0("sample_", 1:10)),
                        "sex" = c("M", "F", "M", "M", "F", "F", "F", "M", "F", "F"),
                        "re" = c(1,2,3,1,2,3,1,2,3,1),
                        "BMI" = c(21, 23, 26, 19, 25, 29, 22, 31, 24, 32))
pheno_all <- pheno_all %>% column_to_rownames("sampleID")

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


