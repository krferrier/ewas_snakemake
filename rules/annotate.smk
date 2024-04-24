def get_file(wildcards):
    if config["stratified_ewas"] == "yes":
        in_file = config["out_directory"] + config["association_variable"] + "_ewas_meta_analysis_results_1.txt"
    else:
        in_file = config["out_directory"] + config["association_variable"] + "_ewas_bacon_results" + config["out_type"]
    return(in_file)

rule add_annotation:
    input: 
        in_file = get_file,
        script = "scripts/annotation.R"
    params:
        o_dir = config["out_directory"],
        strat = config["stratified_ewas"],
        assoc = config["association_variable"],
        o_type = config["out_type"]
    output: 
        config["out_directory"] + config["association_variable"] + "_ewas_annotated_results" + config["out_type"]
    conda:
        "../envs/ewas.yaml"
    shell:
        f"""
        Rscript {{input.script}} \
        --input-file {{input.in_file}} \
        --out-dir {{params.o_dir}} \
        --stratified {{params.strat}} \
        --assoc {{params.assoc}} \
        --out-type {{params.o_type}}
        """