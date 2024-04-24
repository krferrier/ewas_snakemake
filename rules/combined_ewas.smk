rule run_combined_ewas:
    input:
        script = "scripts/ewas.R",
        pheno_file = config["pheno"],
        methyl_file = config["mvals"]
    params:
        assoc_var = config["association_variable"],
        cs = config["chunk_size"],
        pt = config["processing_type"],
        n_workers = config["workers"],
        o_dir = config["out_directory"],
        o_type = config["out_type"]
    output: 
        config["out_directory"] + config["association_variable"] + "_ewas_results" + config["out_type"]
    conda:
        "../envs/ewas.yaml"
    shell:
        f"""
        Rscript {{input.script}} \
        --pheno {{input.pheno_file}} \
        --methyl {{input.methyl_file}} \
        --assoc {{params.assoc_var}} \
        --chunk-size {{params.cs}} \
        --processing-type {{params.pt}} \
        --workers {{params.n_workers}} \
        --out-dir {{params.o_dir}} \
        --out-type {{params.o_type}}
        """


rule run_bacon:
    input:
        in_file = config["out_directory"] + config["association_variable"] + "_ewas_results" + config["out_type"],
        script = "scripts/run_bacon.R"
    params:
        o_dir = config["out_directory"],
        o_type = config["out_type"],
        o_prefix = ""
    output: 
        config["out_directory"] + config["association_variable"] + "_ewas_bacon_results" + config["out_type"],
        expand(config["out_directory"] + "bacon_plots/" + config["association_variable"] + "_{plot}.jpg", plot = PLOTS)
    conda:
        "../envs/ewas.yaml"
    shell:
        f"""
        Rscript {{input.script}} \
        --input-file {{input.in_file}} \
        --out-dir {{params.o_dir}} \
        --out-prefix {{params.o_prefix}} \
        --out-type {{params.o_type}} \
        """
    
    