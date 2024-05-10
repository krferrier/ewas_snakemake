# Script for stratifying phenotype and methylation data.
# If no stratify variables are given, the subset.key will be a list of 1 dataframe with all samples.
suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(dplyr)
  library(tibble)
})

# Import functions for stratifying data
source("scripts/fxns/stratify_fxns.R")

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
parser$add_argument("--stratify", 
                    type="character", 
                    nargs="*", 
                    help="Variable(s) to stratify.")
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
stratify_vars <- args$stratify
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
#                                   STRATIFY DATA                                                  #
####################################################################################################

# Make a key for subsetting the data based on the --stratify variables given. 
# If none were given, the subset.key is a list of all samples. 
subset.key <- make_subset.key(pheno, stratify_vars)

# Stratify the phenotype and methylation data with the subset.key. 
# If no stratification variables were given, the data is not stratified. 
pheno <- stratify.pheno(pheno, stratify_vars)
mvals <- stratify.mvals(mvals)

####################################################################################################
#                                EXPORT STRATIFIED DATA                                            #
####################################################################################################
for (i in names(subset.key)){
  p.sub <- pheno[[i]]
  filename <- paste0(out_dir, i, "/", i, "_pheno", out_type)
  fwrite(p.sub, file = filename, row.names=T)
}
for (i in names(subset.key)){
  m.sub <- mvals[[i]]
  filename <- paste0(out_dir, i, "/", i, "_mvals", out_type)
  fwrite(m.sub, file = filename, row.names=T) 
}