import pandas as pd

def generate_observed_combinations(df, stratify_cols):
    # Convert integer columns to strings
    df[stratify_cols] = df[stratify_cols].astype(str)
    
    # Group the DataFrame by 'stratify_cols' and aggregate the values
    observed_combinations = df.groupby(stratify_cols).size().reset_index()
    
    # Join the values of each combination with "_"
    observed_combinations['combination'] = observed_combinations[stratify_cols].apply('_'.join, axis=1)
    
    return observed_combinations['combination'].tolist()

pheno = pd.read_csv("pheno.csv")
stratify_cols = ["sex"]
generate_observed_combinations(pheno, stratify_cols)