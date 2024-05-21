rule stratify_data:
    input:
        script = "scripts/stratify.R",
        pheno_file = PHENO,
        methyl_file = MVALS
    params:
        strat_vars = ' '.join(STRAT_VARS),
        o_dir = OUT_DIR,
        n_threads = N_WORKERS
    output: 
        temp(expand(OUT_DIR + "{group}/{group}_pheno.fst", group = GROUPS)),
        temp(expand(OUT_DIR + "{group}/{group}_mvals.fst", group = GROUPS))
    conda:
        "../envs/ewas.yaml"    
    shell:
        f"""
        Rscript {{input.script}} \
        --pheno {{input.pheno_file}} \
        --methyl {{input.methyl_file}} \
        --stratify {{params.strat_vars}} \
        --out-dir {{params.o_dir}} \
        --threads {{params.n_threads}}
        """


for group in GROUPS:
    rule:
        name:
            f"run_ewas_{group}"
        input:
            script = "scripts/ewas.R",
            pheno_file = OUT_DIR + f"{group}/{group}_pheno.fst",
            methyl_file= OUT_DIR + f"{group}/{group}_mvals.fst"
        params:
            assoc_var = ASSOC,
            stratified = STRATIFIED,
            cs = config["chunk_size"],
            pt = config["processing_type"],
            n_workers = N_WORKERS,
            o_dir = OUT_DIR + f"{group}/",
            o_type = OUT_TYPE,
            o_prefix = f"{group}"
        output: 
            OUT_DIR + f"{group}/{group}_" + ASSOC + "_ewas_results" + OUT_TYPE
        log:
            f"log/{group}_ewas.log"
        conda:
            "../envs/ewas.yaml"
        shell:
            f"""
            export R_PROGRESSR_ENABLE=TRUE 
            Rscript {{input.script}} \
            --pheno {{input.pheno_file}} \
            --methyl {{input.methyl_file}} \
            --assoc {{params.assoc_var}} \
            --stratified {{params.stratified}} \
            --chunk-size {{params.cs}} \
            --processing-type {{params.pt}} \
            --workers {{params.n_workers}} \
            --out-dir {{params.o_dir}} \
            --out-type {{params.o_type}} \
            --out-prefix {{params.o_prefix}} \
            > {{log}} 2>&1
            """
    rule:
        name:
            f"run_bacon_{group}"
        input:
            in_file = OUT_DIR + f"{group}/{group}_" + ASSOC + "_ewas_results" + OUT_TYPE,
            script = "scripts/run_bacon.R"
        params:
            o_dir = OUT_DIR + f"{group}/",
            o_type = OUT_TYPE,
            o_prefix = f"{group}"
        output: 
            OUT_DIR + f"{group}/{group}_" + ASSOC + "_ewas_bacon_results" + OUT_TYPE,
            expand(OUT_DIR + f"{group}/bacon_plots/{group}_" + ASSOC + "_{plot}.jpg", plot = PLOTS)
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

rule make_metal_script:
    input:
        script = "scripts/metal_cmd.sh",
        in_files = expand(OUT_DIR + "{group}/{group}_" + ASSOC + "_ewas_bacon_results" + OUT_TYPE, group=GROUPS)
    params:
        out_prefix = OUT_DIR + ASSOC + "_ewas_meta_analysis_results_"
    output:
        "scripts/meta_analysis_script.sh"
    shell:
        f"sh {{input.script}} {{input.in_files}} {{params.out_prefix}}"
    

rule run_metal:
    input: 
        script = "scripts/meta_analysis_script.sh"
    output:
        meta_analysis_results
    singularity:
        "library://krferrier/metal/meta_analysis:metal"
    shell: 
        f"metal {{input.script}}"
