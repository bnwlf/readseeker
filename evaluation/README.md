# Model Evaluation

## Requirements
To reproduce the evaluation you need:
- a working singularity and conda installation
- a conda environment with snakemake
- an environment to run jupyter notebooks

## Step 1 - Data Download
First download the readseeker model from [huggingface](https://huggingface.co/bnwlf/ReadSeeker)

## Step 2 - Download the required Reference and Sample Data
1. (Optional)To fully reproduce the process create a new conda environment using `conda create env -f conda/download_env.yaml` (It's only required to generate the also already given bedfiles.)
2. Download the sample and reference data by using the bash respective commands from [download_prepare_reference_sequences.sh](download_prepare_reference_sequences.sh) or execute the whole script to generate the bed files. (**Attention:** Executing the whole script will take a long period of time(>1d))

## Step 3 - Build the singularity image
`singularity build --fakeroot singularity/dnabert.sif singularity/dnabert.def`

## Step 4 - Execute the benchmark workflow
Execute the snakemake pipeline "[realreads_newtrained_singularity2.snakemake](realreads_newtrained_singularity2.snakemake)". Use the provided `evaluation_config.yaml` as config file and activate the `--use-conda` and 
`--use-singularity` flags.

## Step 5 - Plots and Statistics
Finally rerun the jupyter notebook to generate the plots and statistics.
