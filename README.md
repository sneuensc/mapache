***!!! UNDER CONSTRUCTION !!!***

The pipeline is under development.

You are welcome to try it out. On the [Wiki](https://github.com/sneuensc/mapache/wiki), there is already a short manual describing how to get started. Happy to get any feedback/propositions.

***mapache*** is a lightweighted mapping pipeline for ancient DNA using the workflow manager *snakemake*

SUMMARY:

Ancient DNA is degraded and contaminated. Most standard bioinformatics tools to align sequenced reads have been designed for modern data and cannot be used “out of the box” to accommodate the features typical of ancient DNA. In this work, we propose a robust pipeline to align ancient DNA data and to have a first rough idea of the authenticity of the data. The implemented modules consist of the steps needed to go from a simple ‘fastq' file to a final ‘bam' file. The steps include quality control, mapping, filtering, duplicate removal, damage pattern inference, realignment, statistics to assess authenticity and inference of the sex of the organism. A final graphical report summarizes the statistics of the different modules allowing a quick overview of the data. The pipeline is implemented in the workflow manager snakemake providing the flexibility to run the pipeline on a workstation, a high memory server or a cluster. The installation of the pipeline and the underlying programs is made easy using conda. The pipeline may be used out of the box, or may be adapted with little knowledge of python and snakemake.

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=sneuensc/mapache).


# Quick start
To install mapache, you need to clone the repository, prepare a plain text file with the descriptions of the FASTQ files to be mapped, and adapt the configuration of the pipeline.


## Install

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

## Prepare samples file

A samples file for the test dataset is provided (`config/samples.tsv`). If you want to run it on your own datasets, you can prepare a tab-separated file containing the same columns as the example.

Example of a sample file for single-end libraries:

```
SM          LB        ID           Data                        MAPQ  PL
ind1        a_L2      a_L2_R1_001  reads/a_L2_R1_001.fastq.gz  30    ILLUMINA
ind1        a_L2      a_L2_R1_002  reads/a_L2_R1_002.fastq.gz  30    ILLUMINA
ind1        b_L2      b_L2_R1_001  reads/b_L2_R1_001.fastq.gz  30    ILLUMINA
ind1        b_L2      b_L2_R1_002  reads/b_L2_R1_002.fastq.gz  30    ILLUMINA
```

Example of a sample file for paired-end and single-end libraries:

```
SM          LB        ID           Data1                       Data2                       MAPQ  PL
ind1        a_L2      a_L2_R1_001  reads/a_L2_R1_001.fastq.gz  reads/a_L2_R2_001.fastq.gz  30    ILLUMINA
ind1        a_L2      a_L2_R1_002  reads/a_L2_R1_002.fastq.gz  reads/a_L2_R2_001.fastq.gz  30    ILLUMINA
ind1        b_L2      b_L2_R1_002  reads/b_L2_R1_002.fastq.gz  NULL                        30    ILLUMINA
```

In the first example, four fastq files will be mapped. They were generated from two different libraries (here, labelled as `a_L2` and `b_L2`) from a single sample (`ind1`). The reads will be mapped and retained if the mapping quality is above 30 (`MAPQ` column).

In the second example, there is still only one sample (`ind1`), and two libraries, sequenced in paired-end (`a_L2`) and single-end (`b_L2`) mode.

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



## Dry run
We highly recommend to make use of dry runs to get an idea of the jobs that will be executed.
```
# print jobs to be executed
snakemake -n
# Visualization of the pipeline in a directed acyclic graph (DAG).
snakemake dag --dag | dot -Tsvg > dag.svg               
```

This command below will produce a figure with the rules that will be run. 
```
snakemake --rulegraph |dot -Tpng > rulegraph.png
```


## Adapt config file
By default, mapache is configured to map ancient data to a human reference genome. You need to specify the path to the reference genome in FASTA format in the configuration file (`config/config.yaml`) provided by mapache.

If you need to further modify the pipeline (for instance, to ommit one step), see this link.

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

## Run mapping pipeline 

```
snakemake 
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
