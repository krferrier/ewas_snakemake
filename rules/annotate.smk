def get_file(wildcards):
    if config["stratified_ewas"] == "yes":
        in_file = meta_analysis_results
    else:
        in_file = bacon_results
    return(in_file)

rule get_annotation_data:
    output:
        protected("annotation_files/EPIC_hg38.tsv.gz"),
        protected("annotation_files/EPIC_snp_key.tsv.gz")
    shell:
        f"""
        wget https://zhouserver.research.chop.edu/InfiniumAnnotation/20210615/EPIC/EPIC.hg38.manifest.gencode.v36.tsv.gz \
        -O annotation_files/EPIC_hg38.tsv.gz
        wget https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/EPIC/EPIC.hg38.commonsnp.tsv.gz \
        -O annotation_files/EPIC_snp_key.tsv.gz
        """

rule add_annotation:
    input: 
        rules.get_annotation_data.output,
        in_file = get_file,
        script = "scripts/annotation.R"
    params:
        o_dir = OUT_DIR,
        strat = STRATIFIED,
        assoc = ASSOC,
        o_type = OUT_TYPE
    output: 
        annotated_results
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