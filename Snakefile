import pandas as pd
from helper_fxns import generate_observed_combinations
configfile: "config.yml"

#----SET VARIABLES----#
PHENO = config["pheno"]
MVALS = config["mvals"]
ASSOC = config["association_variable"]
STRATIFIED = config["stratified_ewas"]
STRAT_VARS = config["stratify_variables"]
CHUNK_SIZE = config["chunk_size"]
PROCESSING_TYPE = config["processing_type"]
N_WORKERS = config["workers"]
OUT_DIR = config["out_directory"]
OUT_TYPE = config["out_type"]
PLOTS = ["traces", "posteriors", "fit", "qqs"]

# Stratified EWAS
GROUPS = generate_observed_combinations(df=pd.read_csv(config["pheno"]), stratify_cols=config["stratify_variables"])

#---- INPUT & OUTPUT FILES ----#
# Final output results, stratified or not
annotated_results = OUT_DIR + ASSOC + "_ewas_annotated_results" + OUT_TYPE
manhattan_qq_plot = OUT_DIR + ASSOC + "_ewas_manhattan_qq_plots.jpg"

# Combined (not stratified) EWAS outputs
raw_results = OUT_DIR + ASSOC + "_ewas_results" + OUT_TYPE
bacon_results = OUT_DIR + ASSOC + "_ewas_bacon_results" + OUT_TYPE
bacon_plots = expand(OUT_DIR + "bacon_plots/" + ASSOC + "_{plot}.jpg", plot=PLOTS)

# Stratified EWAS outputs
strat_raw_results = expand(OUT_DIR + "{group}/{group}_" + ASSOC + "_ewas_results" + OUT_TYPE, group=GROUPS)
strat_bacon_results = expand(OUT_DIR + "{group}/{group}_" + ASSOC + "_ewas_bacon_results" + OUT_TYPE, group=GROUPS)
strat_bacon_plots = expand(OUT_DIR + "{group}/bacon_plots/{group}_" + ASSOC + "_{plot}.jpg", group=GROUPS, plot=PLOTS)
meta_analysis_results = OUT_DIR + ASSOC + "_ewas_meta_analysis_results_1.txt"

#---- DETERMINE INPUT FILES FOR RULE ALL ----#
if STRATIFIED == "yes":
    in_files = [PHENO, MVALS, strat_raw_results, strat_bacon_results,
                strat_bacon_plots, meta_analysis_results,
                annotated_results, manhattan_qq_plot]
else:
    in_files = [PHENO, MVALS, raw_results, bacon_results, 
                bacon_plots, annotated_results, 
                manhattan_qq_plot]

#---- BEGIN WORKFLOW ----#
rule all:
    input:
        in_files


include: "rules/combined_ewas.smk"
include: "rules/stratified_ewas.smk"
include: "rules/annotate.smk"
include: "rules/plots.smk"