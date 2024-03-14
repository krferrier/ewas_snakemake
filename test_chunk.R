# update chunking function to take in command-line arguments 

library(argparse)

# Define command line arguments
parser <- argparse::ArgumentParser(description = "Chunk methylation data")
parser$add_argument("--chunk-size", type = "integer", nargs="?", const = 1000, help = "number of CpGs per chunk")
args <- parser$parse_args()

# parse the stratify argument to integer
chunk_size <- args$chunk_size
chunk_size 

# dummy methylation data for testing
num_rows <- 10
num_cols <- 100
column_names <- paste0('cg', 1:num_cols)
row_names <- paste0('sample_', 1:num_rows)

df <- data.frame(matrix(nrow = num_rows, ncol = num_cols))
rownames(df) <- row_names
colnames(df) <- column_names

for (i in 1:num_rows) {
  for (j in 1:num_cols) {
    df[i, j] <- runif(1, min = -5, max = 5)
  }
}

# vector of the cpgs 
cpgs <- colnames(df)

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

chunked_df <- chunk.df(df, cpgs)
chunked_df