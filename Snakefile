import pandas as pd
from helper_fxns import generate_observed_combinations
configfile: "config.yml"

#----SET VARIABLES----#
STRATIFIED = config["stratified_ewas"]
GROUPS = generate_observed_combinations(df=pd.read_csv(config["pheno"]), stratify_cols=config["stratify_variables"])
PLOTS = ["traces", "posteriors", "fit", "qqs"]

def rule_all_combined():
    files1 = [config["out_directory"] + config["association_variable"] + "_ewas_results" + config["out_type"],
        config["out_directory"] + config["association_variable"] + "_ewas_bacon_results" + config["out_type"],
        config["out_directory"] + config["association_variable"] + "_ewas_annotated_results" + config["out_type"]]
    files2 = expand(config["out_directory"] + "bacon_plots/" + config["association_variable"] + "_{plot}.jpg", plot=PLOTS)
    files = files1 + files2
    return(files)

def rule_all_stratified(): 
    files1 = expand(config["out_directory"] + "{group}/{group}_" + config["association_variable"] + "_ewas_results" + config["out_type"], group=GROUPS)
    files2 = expand(config["out_directory"] + "{group}/{group}_" + config["association_variable"] + "_ewas_bacon_results" + config["out_type"], group=GROUPS)
    files3 = expand(config["out_directory"] + "{group}/bacon_plots/{group}_" + config["association_variable"] + "_{plot}.jpg", group=GROUPS, plot=PLOTS)
    files4 = ["scripts/meta_analysis_script.sh", 
        config["out_directory"] + config["association_variable"] + "_ewas_meta_analysis_results_1.txt",
        config["out_directory"] + config["association_variable"] + "_ewas_meta_annotated_results" + config["out_type"]]
    files = files1 + files2 + files3 + files4
    return(files)

def target_files(wildcards):
    if STRATIFIED == "yes":
        files = rule_all_stratified()
    else:
        files = rule_all_combined()
    return(files)

rule all:
    input:
        target_files

include: "rules/combined_ewas.smk"
include: "rules/stratified_ewas.smk"