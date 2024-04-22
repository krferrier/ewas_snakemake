rule run_stratified_ewas:
    input:
        script = "scripts/ewas.R",
        pheno_file = config["pheno"],
        methyl_file = config["mvals"]
    params:
        assoc_var = config["association_variable"],
        strat_vars = ' '.join(config["stratify_variables"]),
        cs = config["chunk_size"],
        pt = config["processing_type"],
        n_workers = config["workers"],
        o_dir = config["out_directory"],
        o_type = config["out_type"]
    output: 
        expand(config["out_directory"] + "{group}/{group}_" + config["association_variable"] + "_ewas_results" + config["out_type"], group=GROUPS)
    conda:
        "../envs/ewas.yml"
    shell:
        f"""
        Rscript {{input.script}} \
        --pheno {{input.pheno_file}} \
        --methyl {{input.methyl_file}} \
        --assoc {{params.assoc_var}} \
        --stratify {{params.strat_vars}} \
        --chunk-size {{params.cs}} \
        --processing-type {{params.pt}} \
        --workers {{params.n_workers}} \
        --out-dir {{params.o_dir}} \
        --out-type {{params.o_type}}
        """

for group in GROUPS:
    rule:
        name:
            f"run_bacon_{group}"
        input:
            in_file = config["out_directory"] + f"{group}/{group}_" + config["association_variable"] + "_ewas_results" + config["out_type"],
            script = "scripts/run_bacon.R"
        params:
            o_dir = config["out_directory"] + f"{group}/",
            plots_dir = config["out_directory"] + "bacon_plots/",
            o_type = config["out_type"],
            o_prefix = f"{group}"
        output: 
            config["out_directory"] + f"{group}/{group}_" + config["association_variable"] + "_ewas_bacon_results" + config["out_type"],
            expand(config["out_directory"] + f"{group}/bacon_plots/{group}_" + config["association_variable"] + "_{plot}.jpg", plot = PLOTS)
        shell:
            f"""
            Rscript {{input.script}} \
            --input-file {{input.in_file}} \
            --out-dir {{params.o_dir}} \
            --out-prefix {{params.o_prefix}} \
            --out-type {{params.o_type}} \
            """

rule make_metal_script:
    input:
        script = "scripts/metal_cmd.sh",
        in_files = expand(config["out_directory"] + "{group}/{group}_" + config["association_variable"] + "_ewas_bacon_results" + config["out_type"], group=GROUPS)
    params:
        out_prefix = config["out_directory"] + config["association_variable"] + "_ewas_meta_analysis_results_"
    output:
        "scripts/meta_analysis_script.sh"
    shell:
        f"sh {{input.script}} {{input.in_files}} {{params.out_prefix}}"
    

rule run_metal:
    input: 
        script = "scripts/meta_analysis_script.sh"
    output:
        config["out_directory"] + config["association_variable"] + "_ewas_meta_analysis_results_1.txt"
    singularity:
        "library://krferrier/metal/meta_analysis:metal"
    shell: 
        f"metal SOURCE {{input.script}}"