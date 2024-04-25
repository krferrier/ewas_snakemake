rule plot_results:
    input: 
        in_file = annotated_results,
        script = "scripts/plots.R"
    params:
        o_dir = OUT_DIR,
        strat = STRATIFIED,
        assoc = ASSOC
    output: 
        manhattan_qq_plot
    conda:
        "../envs/ewas.yaml"
    shell:
        f"""
        Rscript {{input.script}} \
        --input-file {{input.in_file}} \
        --out-dir {{params.o_dir}} \
        --stratified {{params.strat}} \
        --assoc {{params.assoc}} 
        """