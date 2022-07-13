# MAPACHE

***Mapache*** is a lightweight mapping pipeline for ancient DNA using the workflow manager *snakemake*

You are welcome to try out mapache. Please beware that the pipeline is under development, and we are happy to get any feedback/propositions.

You can see an extensive documentation on how to use mapache in the [Wiki](https://github.com/sneuensc/mapache/wiki). 
If you already have some experience with DNA mapping and/or Snakemake, you can follow the quick guide below.


The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=sneuensc/mapache).


# Quick guide
To install mapache, you need to clone the repository and create a conda environment following the instructions below. This will automatically install the software needed to map FASTQ files to a reference genome and create a report with the mapping statistics.

To execute mapache, you can move to the cloned directory (`mapache`) or symlink its content (the directories `config/`,  `results/`, `workflow/`, and `test_data/` if you want to run the test) to your working directory.

You will mostly interact with mapache through its **configuration file** (`config/config.yaml`), where you can tweak the parameteres of the pipeline, and the **samples file** (`config/samples.tsv`), in which you can list all the input FASTQ files.


## Installation

```
#--------------------------------------------------------------#
## create mamba environment
conda create -n base-mamba -c conda-forge mamba

## activate mamba environment
conda activate base-mamba

#--------------------------------------------------------------#
## clone mapache repository
git clone https://github.com/sneuensc/mapache.git
cd mapache

## create conda environment for mapache
mamba env create -n mapache --file config/mapache-env.yaml

#--------------------------------------------------------------#
## now you can activate the mapache environment
conda activate mapache

```
<p align="center">

<img src="https://github.com/sneuensc/mapache/wiki/images/mapache_install_loop.gif" width="600" height="400"/>

</p>

## Prepare samples file

A samples file for the test dataset is provided (`config/samples.tsv`). If you want to run it on your own datasets, you can prepare a tab-separated file containing the same columns as the example.

Example of a sample file for single-end libraries:

```
SM      LB    ID       Data                            MAPQ   PL
ind1    L1    L1_1     reads/ind1.L1_R1_001.fastq.gz   30     ILLUMINA
ind1    L1    L1_2     reads/ind1.L1_R1_002.fastq.gz   30     ILLUMINA
ind1    L2    L2_1     reads/ind1.L2_R1_001.fastq.gz   30     ILLUMINA
ind1    L2    L2_2     reads/ind1.L2_R1_002.fastq.gz   30     ILLUMINA
```

Example of a sample file for paired-end and single-end libraries:

```
SM      LB    ID       Data1                            Data2                           MAPQ   PL
ind2    L1    L1_1     reads/ind2.L1_R1_001.fastq.gz    reads/ind2.L1_R2_001.fastq.gz   30     ILLUMINA
ind2    L1    L1_2     reads/ind2.L1_R1_002.fastq.gz    reads/ind2.L1_R2_001.fastq.gz   30     ILLUMINA
ind2    L2    L2_1     reads/ind2.L2_R1_002.fastq.gz    NULL                            30     ILLUMINA
```

In the first example, four fastq files will be mapped. They were generated from two different libraries (here, labelled as `L1` and `L2`) from a single sample (`ind1`). The reads will be mapped and retained if the mapping quality is above 30 (`MAPQ` column).

In the second example, there is still only one sample (`ind2`), and two libraries, sequenced in paired-end (`L1`) and single-end (`L2`) mode.

The columns `SM`, `LB`, `ID` and `PL` will be used to annotate the header of the BAM files produced (SM, LB, RG and PL tags, respectively). 

You can add more samples/libraries/fastq files in the input dataset by adding a new row including the 6 or 7 fields indicated above.


The columns of the sample file are:
- **SM:** Sample name. Libraries are merged according to this name.
- **LB:** Library name. Fastq files are merged according to this name.
- **ID:** An ID for the fastq library (examples: id1, fq_1, ind1_lib1_fq2, etc.)
- **Data** (single-end format)**:** Path to the fastq file. The file may be gzipped or not. Path may be absolute or relative to the working directory.
- **Data1** (paired-end format)**:** Path to the forward fastq file (R1) for paired-end data or the fastq file for single-end data. The file may be gzipped or not. Path may be absolute or relative to the working directory.
- **Data2** (paired-end format)**:** Path to the reverse fastq file (R2) for paired-end data or `NULL` for single-end data. The file may be gzipped or not. Path may be absolute or relative to the working directory.
- **MAPQ:** Fastq-file specific mapping quality filtering threshold.
- **PL:** Sequencing platform or any other text (single word) to annotate in the read group line `@RG` in the the bam file.

## Adapt config file
By default, mapache is configured to map ancient data to a human reference genome. You need to specify the path to the reference genome in FASTA format in the configuration file (`config/config.yaml`) provided by mapache.

### Specify reference genome

Mapache will map the input FASTQ files to the reference genome(s) indicated in the config file under the keyword `genome`.
The path to the genome in FASTA format should be indicated with the keyword `fasta`.

In the example below, all the samples will be mapped to two versions of the human genome. The final BAM files will be named with the format `{sample_ID}.{genome_name}.bam` (e.g., ind1.hg19.bam, ind1.GRCh38.bam).

```
genome: 
    hg19: 
        fasta: /path/to/hg19/genome/hs.build37.1/hs.build37.1.fa
    GRCh38:
        fasta: /path/to/GRCh38/genome/GCA_000001405.15_GRCh38.fa
```

### Enable/disable steps and specify memory and runtime
Mapache will perform different steps (see below) in order to produce a BAM file. Most of the steps are optional, but run by default.

The config file contains entries corresponding to each of the steps that can be run and customized. For example, the entry with the options for AdapterRemoval2 looks like this:

```
# adapter_removal (optional)
adapterremoval:
    run: True
    threads: 4 
    mem: 4 ## in GB
    time: 2
    params: '--minlength 30 --trimns --trimqualities'

```

If you need to modify the pipeline (for instance, to ommit one step), you can set its option to `run: True` or `run: False` in the config file. Additionally, if you intend to run mapache on an HPC system, you can specify the number of `threads`, memory (`mem` in GB), and runtime (`time` in hours) to be allocated to each step. Finally, you can pass additional parameters to the tool to be executed (e.g. `--minlength 30` to AdapterRemoval2) with the keyword `params`.

The following main steps are included in ***mapache***:

```
for each fastq file:
    1. adapter removal with AdapterRemoval (optional; run by default)
    2. mapping with bwa aln (default), bwa mem or bowtie2
for each library:
    3. marking duplicates with MarkDuplicates (optional; duplicates are removed by default)
    4. mapDamage is run, but by default not used for the next steps (optional)
for each sample:
    5. re-alignment of reads with GATK (optional; run by default)
    6. re-computation of the md flag with samtools (optional; run by default)
```



# Running mapache

To run mapache on your working directory, you will need to copy or create symbolic links to the following directories:
- workflow (mandatory)
- config (mandatory)
- test_data (optional, only if you want to run the tutorial)
- slurm (optional, useful to submit jobs in HPC systems with slurm as cluster manager)

Assuming that you are in the mapache directory, you can use the next lines to create an alias to copy or symlink those directories. You can adjust it by including only the directories that you want to copy/symlink.

```
# alias to cp directories
 echo 'alias copy_mapache="cp -r $(pwd -P)/{config,workflow,test_data,slurm} ." ' >> ~/.bash_profile
# alias to symlink directories
 echo 'alias symlink_mapache="ln -s $(pwd -P)/{config,workflow,test_data,slurm} ." ' >> ~/.bash_profile

source ~/.bash_profile

```

<p align="center">

<img src="https://github.com/sneuensc/mapache/wiki/images/symlink_mapache_loop.gif" width="600" height="400"/>

</p>


Make sure that the paths to 
- the reference genome(s) in the config files
- the file listing the FASTQ files (samples.tsv)

as well as the paths listed in samples.tsv  are properly specified


## Dry run

We highly recommend to make use of dry runs to get an idea of the jobs that will be executed.
```
# print jobs to be executed
snakemake -pn
# Visualization of the pipeline in a directed acyclic graph (DAG).
snakemake dag --dag | dot -Tpng > dag.png               
```

The command below will produce a figure with the rules that will be run. 

```
snakemake --rulegraph |dot -Tpng > rulegraph.png
```


## Run mapping pipeline 


mapache can be run locally by indicating the number of cores available or with a high-performance computing system, by configuring a profile.

Example of execution using only one CPU:

```
snakemake --cores 1
```

If you work on an HPC system managed by slurm, you can use the slurm profile in the repository (`mapache/slurm/`) by symlinking it to your working directory. We recommend to start a screen session prior to the job submission.

Example of a submission of 999 jobs (simoultaneously) with the slurm profile:

```
snakemake --jobs 999 --profile slurm
```

## Create report with mapping statistics

```
snakemake --report report.zip
```


USAGE:

```
## recommended workflow on a server
snakemake dag --dag | dot -Tsvg > dag.svg                Visualization of the pipeline in a directed acyclic graph (DAG).
snakemake -n                                             Dry run
snakemake --report report.html                           Create a html report (after execution)
sed -i "s/'runtime': 0.0/'runtime': 0.1/g" report.html   Correct report


## recommended workflow on a cluster (e.g. wally)
snakemake dag --dag | dot -Tsvg > dag.svg                Visualization of the pipeline in a directed acyclic graph (DAG).
snakemake -n 		                                 Dry run
snakemake mapping --profile wally                        Run all mappings (note that the profile has to be the last argument)
snakemake stats --profile wally                          Run all statistics (and all mapping steps which are required to get them...
snakemake --report report.html                           Create a html report (after execution)
sed -i "s/'runtime': 0.0/'runtime': 0.1/g" report.html   Correct report


## minimal set of commands to execute the workflow
snakemake                                                Run all (mappings & stats)
snakemake --report report.html                           Create a html report (after execution)
sed -i "s/'runtime': 0.0/'runtime': 0.1/g" report.html   Correct report

## commands to test if the workflow makes sense
snakemake dag --rulegraph | dot -Tsvg > rulegraph.svg Visualization the interplay of the rules.
snakemake dag --dag | dot -Tsvg > dag.svg  Visualization of the pipeline in a directed acyclic graph (DAG). 
snakemake -p -n                            Print out the commands

## helpfull commands to controll the execution
--profile axiom                           To send it on the cluster (must be the last argument).
-R RULE_NAME                              Force a start at at least the given rule.
--until RULE_NAME                         Run until the given rule (included).
--rerun-incomplete                        Re-run all jobs the output of which is recognized as incomplete.
--configfile FILE                         Define the config file (default config.yaml, this should one of the first arguments).
-t                                        Reset the timestamp that the output is not re-computed.
-p                                        Print out the shell commands that will be executed. 
```

# Citing mapache
We are preparing a manuscript describing mapache. In the meantime, if you use mapache for your study, please refer to mapache's repository on github (https://github.com/sneuensc/mapache) and cite the tools that you used within mapache. See the table below for a list of tools used at each step.



