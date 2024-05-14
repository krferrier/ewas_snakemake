# Overview
This repository is a snakemake workflow for performing Epigenome-Wide Association Studies (EWAS) of methylation measured by the Illumina EPIC (850k) methylation array. This workflow can perform a standard EWAS or a stratified EWAS based on the variables provided for stratification. In either the standard or stratified EWAS, a linear regression model is used where the outcome is the methylation value (Beta or M-values) and the trait/phenotype you want to perform association testing with is the main predictor. All other variables included in the phenotype dataframe will be added as covariates to the model.

Results from the linear regression analyses will be adjusted for bias and inflation using a Bayesian approach implemented by the [BACON](https://www.bioconductor.org/packages/release/bioc/html/bacon.html) R package. It is strongly recommended to assess the performance of the bias and inflation adjustment by viewing the traces, posteriors, fit, and qq plots output at this step. Further details on how to assess the performance plots are provided below. If the EWAS was stratified, after adjustment for bias and inflation, the results from all the strata will be combined using an inverse-variance weighted meta-analysis approach with the command line tool [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation).

The final results will be annotated with hg38/GRCh38 human genome build information collated by [Wanding Zhou](https://zwdzwd.github.io/InfiniumAnnotation). A manhattan and qq plot of the final results will also be output.

# Dependencies
* Conda/Mamba
* Snakemake
* Singularity (only for stratified EWAS)

# Using this Workflow

## Modifying the Configuration File
This workflow uses a configuration file, `config.yml`, to specify the paths for input and output files, what kind of EWAS to perform (standard or stratified), and parallelization parameters. 
### Input Files
This workflow is intended to be used with phenotype and methylation data that has already been cleaned and has no missing/NA data. The first column of the phenotype data and the methylation data must be the sample IDs. Examples of what the phenotype and methylation data should look like can be found in the `data/` directory. 
### Output Files
For a standard EWAS the output files will include the raw linear regression results, the BACON bias- and inflation-adjusted results and plots, and the final annotated results in either .csv or .csv.gz format, and a manhattan and qq plot .jpg. These files will all be output to the 'out_directory' specified in the config file. A stratified EWAS will output the raw and BACON-adjusted results and plots in separate subdirectories for each strata, the results from the meta-analysis of all the strata, as well as the final annotated results and manhattan and qq plot. 
### Parallelization Parameters
The rule 'ewas' first chunks the methylation dataset into sets of CpGs where the length of each set is specified by the parameter `chunks` in the config file (default of 1000 CpGs if a number is not provided). Then, linear regressions are run for each chunk and the results combined back into one dataframe once all chunks have been processed. Run sequentially, this step could take several hours. However, each chunk can be processed in parallel to reduce the total computation time. The main step of the EWAS utilizes the [BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) R package to perform the linear regressions in parallel. You can specify the type of parallelization to be performed at this step in the config file. The options include: sequential, multisession (threads), multicore, or cluster. You must also specify the number of workers that you want to be used for the parallelization in the 'workers' parameter of the config file. For more details on the different types of parallelization you can reference the BiocParallel documentation.

In you are performing a stratified analysis, you can have each strata processed in parallel using the `-j` parameter in the snakemake command. For resource allocation, keep in mind that each EWAS performed will use the number of workers specified in the config file, so you will need n=(workers x jobs) resources available. For example, if you use multicore parallelization with 2 workers and `-j 2` in the snakemake command, you will need to have 4 cores available to run the analysis. 

## Run EWAS

Once you have modified the `config.yml` file to match your projects specifications and computational resources you are ready to run the EWAS. 

For a standard (non-stratified) EWAS:
```shell
snakemake -j 1 --use-conda
```

For a stratified EWAS:
```shell
snakemake -j 1 --use-conda --use-singularity
```
