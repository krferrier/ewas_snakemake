####################################################################################################
#                                FUNCTIONS FOR CHUNKING                                            #
####################################################################################################

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

chunk.df <- function(df, cpgs){
  df <- lapply(cpgs, function(x){
    df.sub <- df[,colnames(df) %in% x]
    return(df.sub)
  })
  return(df)
}