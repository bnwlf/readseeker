# Model Evaluation

## Requirements
To reproduce the evaluation you need:
- a working singularity and conda installation
- a conda environment with snakemake and jupyter notebook
- an environment to run jupyter notebooks
- Apptainer / Singularity engine
## Step 1 - Data Download
First download the readseeker model from [huggingface](https://huggingface.co/bnwlf/ReadSeeker)

## Step 2 - Download the required Reference and Sample Data
1. (Optional)To fully reproduce the process create a new conda environment using `conda create env -f conda/download_env.yaml` (It's only required to generate the also already given bedfiles.)
2. Download the sample and reference data by using the bash respective commands from [download_prepare_reference_sequences.sh](download_prepare_reference_sequences.sh) or execute the whole script to generate the bed files. (**Attention:** Executing the whole script will take a long period of time(>1d))

## Step 3 - Build the singularity image
`apptainer build  singularity/dnabert.sif singularity/dnabert.def`

## Step 4 - Execute the benchmark workflow
Execute the snakemake pipeline "[realreads_newtrained_singularity2.snakemake](evaluation.snakemake)". Use the provided `evaluation_config.yaml` as config file and activate the `--use-conda` and 
`--use-singularity` flags. In case you want to use the model with CUDA acceleration you need apply the `--singularity-args "--nv -B .:/dum"` flags to snakemake.
## Step 5 - Plots and Statistics
Finally rerun the jupyter notebook "[ReadSeeker_Evaluation_and_Statistics_stable_test_and_extension.ipynb](ReadSeeker_Evaluation_and_Statistics_stable_test_and_extension.ipynb)" to generate the plots and statistic tables from the paper.


## Known issues
In some cases conda will install a new libc version as apptainer/singularity was linked. This will result in errors like:
```
INFO:    Running post scriptlet
FATAL:   exec /.singularity.d/libs/fakeroot failed: a shared library is likely missing in the image
FATAL:   While performing build: while running engine: exit status 255
```
In this case try do deactivate conda or install the version of libc  used for apptainer!

See also: https://github.com/apptainer/apptainer/issues/783
