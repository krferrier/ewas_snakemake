####################################################################################################
#                              FUNCTIONS FOR STRATIFYING                                           #
####################################################################################################


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

# Stratify other datasets with the subset.key
stratify.pheno <- function(pheno.df, stratify_vars){
    pheno <- lapply(subset.key, function(i){
                df.sub <- pheno.df[rownames(pheno.df) %in% i,]
                df.sub <- df.sub[order(rownames(df.sub)),] 
                if(length(stratify_vars) != 0){
                    df.sub <- df.sub %>% 
                    dplyr::select(-all_of(stratify_vars)) 
                }
                return(df.sub)
            })
    names(pheno) <- names(subset.key)
    return(pheno)
}

stratify.mvals <- function(mvals.df){
    mvals <- lapply(subset.key, function(i){
                df.sub <- mvals.df[rownames(mvals.df) %in% i,]
                df.sub <- df.sub[order(rownames(df.sub)), order(colnames(df.sub))]
                return(df.sub)
            })
    names(mvals) <- names(subset.key)
    return(mvals)
}
