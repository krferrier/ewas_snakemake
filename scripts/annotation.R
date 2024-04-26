# Script for running adding annotation data to ewas results
# Import libraries
suppressPackageStartupMessages({
    library(argparse)
    library(dplyr)
    library(data.table)
})
# Define command line arguments
parser <- argparse::ArgumentParser(description="Script for adding annotation data to ewas results")
parser$add_argument('--input-file', '-i',
                    required=TRUE,
                    help="Path to ewas results data")
parser$add_argument('--out-dir',
                    required=TRUE,
                    help="Path to output directory")
parser$add_argument('--stratified',
                    choices=c("yes", "no"), 
                    default="no",
                    help="Results from a stratified analysis: yes or no")
parser$add_argument("--assoc", 
                    required=TRUE,
                    type="character", 
                    nargs=1, 
                    help="Association variable EWAS was performed with.")
parser$add_argument('--out-type', 
                    type="character",
                    choices=c(".csv", ".csv.gz"), 
                    nargs="?",                    
                    const=".csv",
                    default=".csv",  
                    help="Output file type: CSV or CSV.GZ")

# parse arguments
args <- parser$parse_args()
results <- args$input_file
out_dir <- args$out_dir
stratified <- args$stratified
assoc <- args$assoc
out_type <- args$out_type

# Read in EWAS summary statistics
ewas <- fread(results)

# Load annotation data
hg38 <- fread("annotation_files/EPIC_hg38.tsv.gz")
hg38 <- hg38 %>% dplyr::rename(cpgid = "probeID") 
cpg.to.rs <- fread("annotation_files/EPIC_snp_key.tsv.gz")
cpg.to.rs <- cpg.to.rs %>%
    dplyr::select(probeID, snpID, snpChrm, snpEnd, snpRef,
                snpAlt, GAN, GAC, GAF,distance) %>%
    distinct(probeID, .keep_all = T) %>% 
    dplyr::rename("cpgid" = "probeID")
annotation <- left_join(hg38, cpg.to.rs, by = "cpgid")
rm(hg38, cpg.to.rs)

if(stratified=="no"){
    ewas <- left_join(ewas, annotation, by = "cpgid")
    ewas <- ewas[order(ewas$bacon.pval),]
} else{
        ewas <- left_join(ewas, annotation, by = c("MarkerName"="cpgid"))
        ewas <- ewas %>% dplyr::select(-Allele1, -Allele2)
        ewas$"P-value" <- as.numeric(ewas$"P-value")
        ewas <- ewas[order(ewas$"P-value"),]
}
file_name <- paste0(out_dir, assoc, "_ewas_annotated_results", out_type)
fwrite(ewas, file = file_name)
