import pandas as pd

def generate_observed_combinations(df, stratify_cols):
    if len(stratify_cols) > 0:
        # Convert integer columns to strings
        df[stratify_cols] = df[stratify_cols].astype(str)
        
        # Group the DataFrame by 'stratify_cols' and aggregate the values
        observed_combinations = df.groupby(stratify_cols).size().reset_index()
        
        # Join the values of each combination with "_"
        observed_combinations['combination'] = observed_combinations[stratify_cols].apply('_'.join, axis=1)
        
        return observed_combinations['combination'].tolist()
    else:
        # If no stratify_cols provided, return a list with only "all"
        return ['all']

