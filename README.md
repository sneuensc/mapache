# Mapache

**Mapache** ([maˈpa.t͡ʃe]) is a lightweight mapping pipeline for ancient DNA using the workflow manager *Snakemake*.

Visit the [Wiki](https://github.com/sneuensc/mapache/wiki) for extensive documentation on how to use mapache and follow the [tutorial](https://github.com/sneuensc/mapache/wiki/4.-Tutorial:-Running-mapache) to map and impute the test dataset. 

If you already have some experience with DNA mapping and/or Snakemake, you can follow the quick guide below.


The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=sneuensc/mapache).

The following main steps are included in **mapache**:

```
for each fastq file:
    - subsampling of a small subset of reads (options, not run by default; useful for quick tests)
    - adapter removal with AdapterRemoval (optional; run by default)
    - mapping with bwa aln (default), bwa mem or bowtie2
for each library:
    - marking duplicates with MarkDuplicates (optional; duplicates are removed by default)
    - rescaling of base qualitiy scores with mapDamage (optional, not run by default)
for each sample:
    - re-alignment of reads with GATK (optional; run by default)
    - re-computation of the md flag with samtools (optional; run by default)
    - imputation of low-coverage genomes with GLIMPSE (optional, not run by default)
```


Moreover, mapache allows you to map large datasets to a single reference genome or to multiple genomes. 

Mapache is designed having in mind that one or multiple DNA libraries are generated per sample, and such libraries are sequenced at least once (e.g., for screening), but usually multiple times (normally prioritizing the highest-quality libraries). This results in several FASTQ files, which have to be mapped several times over the course of a project in order to generate a single BAM file.

The goal of mapache is to make the mapping process as easy and transparent as possible, whether you are mapping to a single or multiple genomes, mapping for the first time, needing to update a BAM file or even impute low-coverage genomes. 

<p align="center">

<img src="https://github.com/sneuensc/mapache/wiki/images/problem.001.jpeg"/>

</p>

# Quick guide
<p align="center">

<img src="https://github.com/sneuensc/mapache/wiki/images/how_to_run.001.jpeg" />

</p>


To install mapache, you need to clone the repository and create a conda environment following the instructions below. This will automatically install the software needed to map FASTQ files to a reference genome and create a report with the mapping statistics.

To execute mapache, you can move to the cloned directory (`mapache`) or symlink its content (the directories `config/`,  `results/`, `workflow/`, and `test_data/` if you want to run the test) to your working directory.

You will mostly interact with mapache through its **configuration file** (`config/config.yaml`), where you can tweak the parameteres of the pipeline, and the **samples file** (`config/samples.tsv`), in which you can list all the input FASTQ files.



## Installation

```bash
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

```yaml
SM      LB    ID       Data
ind1    L1    L1_1     reads/ind1.L1_R1_001.fastq.gz
ind1    L1    L1_2     reads/ind1.L1_R1_002.fastq.gz
ind1    L2    L2_1     reads/ind1.L2_R1_001.fastq.gz
ind1    L2    L2_2     reads/ind1.L2_R1_002.fastq.gz
```

Example of a sample file for paired-end and single-end libraries:

```yaml
SM      LB    ID       Data1                            Data2
ind2    L1    L1_1     reads/ind2.L1_R1_001.fastq.gz    reads/ind2.L1_R2_001.fastq.gz
ind2    L1    L1_2     reads/ind2.L1_R1_002.fastq.gz    reads/ind2.L1_R2_001.fastq.gz
ind2    L2    L2_1     reads/ind2.L2_R1_002.fastq.gz    NULL
```

In the first example, four fastq files will be mapped. They were generated from two different libraries (here, labelled as `L1` and `L2`) from a single sample (`ind1`).

In the second example, there is still only one sample (`ind2`), and two libraries, sequenced in paired-end (`L1`) and single-end (`L2`) mode.

You can add more samples/libraries/fastq files in the input dataset by adding a new row including the 4 or 5 fields indicated above.


The columns of the sample file are:
- **SM:** Sample name. Libraries are merged according to this name.
- **LB:** Library name. Fastq files are merged according to this name.
- **ID:** An ID for the fastq library (examples: id1, fq_1, ind1_lib1_fq2, etc.)
- **Data** (single-end format)**:** Path to the fastq file. The file may be gzipped or not. Path may be absolute or relative to the working directory.
- **Data1** (paired-end format)**:** Path to the forward fastq file (R1) for paired-end data or the fastq file for single-end data. The file may be gzipped or not. Path may be absolute or relative to the working directory.
- **Data2** (paired-end format)**:** Path to the reverse fastq file (R2) for paired-end data or `NULL` for single-end data. The file may be gzipped or not. Path may be absolute or relative to the working directory.

👆🏿 Note that you can select any label you want for the columns SM, LB, and ID, but they may not contain points '.'. That is, they do not need to match any substring or the name of your FASTQ files (FASTQ file names may in contrast contain points). However, we recommend that you stick to meaningful names (e.g.: Denisova, Mota, UstIshim, Anzick, lib1, lib2, lib3_USER, etc.) and ACII characters.

## Adapt config file
By default, mapache is configured to map ancient data to a human reference genome. You need to specify the path to the reference genome in FASTA format in the configuration file (`config/config.yaml`) provided by mapache.

### Specify reference genome

Mapache will map the input FASTQ files to the reference genome(s) indicated in the config file under the keyword `genome`.
The path to the genome in FASTA format should be indicated with the keyword `fasta`.

In the example below, all the samples will be mapped to two versions of the human genome. The final BAM files will be named with the format `{sample_ID}.{genome_name}.bam` (e.g., ind1.hg19.bam, ind1.GRCh38.bam).

```yaml
genome: 
    GRCh37: path_to_reference/GRCh37.fasta
    # GRCh38: path_to_reference/GCA_000001405.15_GRCh38.fa
```

### Enable/disable steps and specify memory and runtime
Mapache will perform different steps (see below) in order to produce a BAM file. Most of the steps are optional, but run by default.

The config file contains entries corresponding to each of the steps that can be run and customized. For example, the entry with the options to clean the fastq files looks like this:

```yaml
# fastq cleaning (optional)
cleaning:
    run: 'adapterremoval'   # options: adapterremoval (default), fastp, False
    params_adapterremoval: '--minlength 30 --trimns --trimqualities'
    params_fastp: ''
    threads: 4 
    mem: 4 ## in GB
    time: 2
```

If you need to modify the pipeline, for instance to ommit the step step, you can set in the config file its option to `run: False`, or `run: 'adapterremovel'` and `run: 'fastp'` to run `AdapterRemoval2` and `fastp`, respectively. Additionally, if you intend to run mapache on an HPC system, you can specify the number of `threads`, memory (`mem` in GB), and runtime (`time` in hours) to be allocated to each step. Finally, you can pass additional parameters to the tool to be executed (e.g. `--minlength 30` to AdapterRemoval2) with the keyword `params_adapterremoval`.



# Running mapache

To run mapache on your working directory, you will need to copy or create symbolic links to the following directories:
- workflow (mandatory)
- config (mandatory)
- test_data (optional, only if you want to run the tutorial)
- slurm (optional, useful to submit jobs in HPC systems with slurm as cluster manager)

Assuming that you are in the mapache directory, you can use the next lines to create an alias to copy or symlink those directories. You can adjust it by including only the directories that you want to copy/symlink.

```bash
# alias to cp directories
 echo "alias copy_mapache='cp -r $(pwd -P)/{config,workflow,test_data,slurm} . ' " >> ~/.bash_profile
# alias to symlink directories
 echo "alias symlink_mapache='ln -s $(pwd -P)/{config,workflow,test_data,slurm} . ' ">> ~/.bash_profile

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
```bash
# print jobs to be executed
snakemake -pn
# Visualization of the pipeline in a directed acyclic graph (DAG).
snakemake dag --dag | dot -Tpng > dag.png               
```

The command below will produce a figure with the rules that will be run. 

```bash
snakemake --rulegraph |dot -Tpng > rulegraph.png
```

<p align="center">

<img src="https://github.com/sneuensc/mapache/wiki/images/dag_rulegraph_loop.gif" width="600" height="400"/>

</p>

Examples for the DAG and rulegraph can be found [here](https://github.com/sneuensc/mapache/wiki/images/dag.png) and [here](https://github.com/sneuensc/mapache/wiki/images/rulegraph.png), respectively.


## Optional: configure a profile to automatically submit jobs to a queue

This is a very handy utility that will allow mapache to automatically schedule jobs for submission to a queuing system. If you work with _a system managed by slurm_, you can *skip this step*, as mapache comes with a slurm profile (`mapache/slurm`). 


Other profiles are available at: https://github.com/Snakemake-Profiles

Note however that the exact configuration might vary depending on your system (i.e., some users might have different accounts on a server for billing purposes). We thus highly encourage you to seek for help with your IT team if you are in doubt about the configuration that suits your system the best. Please note that we are not in charge of developing these profiles. 

## Run mapping pipeline 


mapache can be run locally by indicating the number of cores available or with a high-performance computing system, by configuring a profile.

Example of execution using only one CPU:

```bash
snakemake --cores 1
```


<p align="center">

<img src="https://github.com/sneuensc/mapache/wiki/images/run_mapache_loop.gif" width="600" height="400"/>

</p>


If you work on an HPC system managed by slurm, you can use the slurm profile in the repository (`mapache/slurm/`) by symlinking it to your working directory. We recommend to start a screen session prior to the job submission.

Example of a submission of 999 jobs (simoultaneously) with the slurm profile:

```bash
snakemake --jobs 999 --profile slurm
```

## Create report with mapping statistics

The template using to produce the html report will change depending on the version of Snakemake that you are using with `mapache`. Although we try to include one of the most recent distributions of Snakemake, we think that the template used in older versions is more useful to navigate through the results in the report. We thus recommend you to create a new environment with conda or mamba to create such report:


```bash
mamba create -n snakemake610 snakemake=6.10.0
mamba activate snakemake610
```

Once you have activated the new environment, creating the report is straightforward:

```bash
snakemake --report report.zip
```

or

```bash
snakemake --report report.html
```

We recommend creating the zip version of the report, as it contains the html report in it and it allows you to download any of the output tables or plots by clicking on the links of the report, making it easier to share with your colleagues.

Explore the [zipped report](https://github.com/sneuensc/mapache/wiki/report/report.zip) produced with Snakemake 6.10.0.


## Summary of useful commands

```bash
#------------------------------------------
# Get an idea of the jobs that will be executed. Not mandatory but very useful to spot possible mistakes in the configuration or input files

snakemake dag --dag | dot -Tsvg > dag.svg                Visualization of the pipeline in a directed acyclic graph (DAG).
snakemake dag --rulegraph | dot -Tsvg > rulegraph.svg    Visualization the interplay of the rules.
snakemake -n                                             Dry run
snakemake -p -n                                          Print out the commands in a dry run

#------------------------------------------
# recommended workflow on a cluster (e.g. slurm)
snakemake --jobs 999 --profile slurm                     Run all mappings and imputation (note that the profile has to be the last argument)

#------------------------------------------
# recommended workflow on a single machine without a queuing system
snakemake --cores 16                                     Run all mappings and imputation (assuming that you have 16 CPUs available)

#------------------------------------------
# Create report (after execution)

snakemake --report report.html                           Create a html report 
snakemake --report report.zip                            Create a zip report; useful if you want to download the files by clicking on the HTML report

#------------------------------------------
# Options to control the execution.
# Visit https://snakemake.readthedocs.io/en/stable/ for full documentation

--profile slurm                           To send it on the cluster (must be the last argument).
-R RULE_NAME                              Force a start at least at the given rule.
--until RULE_NAME                         Run until the given rule (included).
--rerun-incomplete                        Re-run all jobs the output of which is recognized as incomplete.
--configfile FILE                         Define the config file (default config/config.yaml).
-t                                        Reset the timestamp that the output is not re-computed.
-p                                        Print out the shell commands that will be executed.
--rerun-triggers mtime                    for a rerun, consider only the time stamp (previously the default 
                                          setting, now any change (code, config,..) lead to a rerun of 
                                          the rule)
```

## Citing mapache
If you use mapache in your study, please cite **Neuenschwander S, Cruz Dávalos DI, et al. (2023) Mapache: a flexible pipeline to map ancient DNA. _Bioinformatics_** (https://doi.org/10.1093/bioinformatics/btad028) and the tools that you used within mapache. See the table below for a list of tools used at each step.


## Tools included in mapache
<table class=" lightable-minimal" style="font-family: &quot;Trebuchet MS&quot;, verdana, sans-serif; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">  </th>
   <th style="text-align:left;"> Version </th>
   <th style="text-align:left;"> Reference </th>
   <th style="text-align:left;"> Link </th>
  </tr>
 </thead>
<tbody>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Workflow manager</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Snakemake </td>
   <td style="text-align:left;"> 7.18.2 </td>
   <td style="text-align:left;"> Mölder, et al. (2021) </td>
   <td style="text-align:left;"> https://github.com/snakemake/snakemake </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Subsample</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> seqtk </td>
   <td style="text-align:left;"> 1.3 </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> https://github.com/lh3/seqtk </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Clean</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> AdapterRemoval2 </td>
   <td style="text-align:left;"> 2.3.2 </td>
   <td style="text-align:left;"> Schubert, et al. (2016) </td>
   <td style="text-align:left;"> https://github.com/MikkelSchubert/adapterremoval </td>
  </tr>
<tr>
  <td style="text-align:left;padding-left: 2em;" indentlevel="1"> fastp </td>
  <td style="text-align:left;"> 0.23.2 </td>
  <td style="text-align:left;"> Chen, et al. (2018) </td>
  <td style="text-align:left;"> https://github.com/OpenGene/fastp </td>
 </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Map</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> BWA aln </td>
   <td style="text-align:left;"> 0.7.17 </td>
   <td style="text-align:left;"> Li and Durbin (2009) </td>
   <td style="text-align:left;"> https://github.com/lh3/bwa </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> BWA mem </td>
   <td style="text-align:left;"> 0.7.17 </td>
   <td style="text-align:left;"> Li and Durbin (2009) </td>
   <td style="text-align:left;"> https://github.com/lh3/bwa </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Bowtie2 </td>
   <td style="text-align:left;"> 2.4.4 </td>
   <td style="text-align:left;"> Langmead and Salzberg (2012) </td>
   <td style="text-align:left;"> https://github.com/BenLangmead/bowtie2 </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Sort</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> SAMtools </td>
   <td style="text-align:left;"> 1.14 </td>
   <td style="text-align:left;"> Danecek, et al. (2021) </td>
   <td style="text-align:left;"> https://github.com/samtools/samtools </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Filter</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> SAMtools </td>
   <td style="text-align:left;"> 1.14 </td>
   <td style="text-align:left;"> Danecek, et al. (2021) </td>
   <td style="text-align:left;"> https://github.com/samtools/samtools </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Merge lanes</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> SAMtools </td>
   <td style="text-align:left;"> 1.14 </td>
   <td style="text-align:left;"> Danecek, et al. (2021) </td>
   <td style="text-align:left;"> https://github.com/samtools/samtools </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Remove duplicates</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Picard MarkDuplicates </td>
   <td style="text-align:left;"> 2.25.5 </td>
   <td style="text-align:left;"> Broad Institute (2019) </td>
   <td style="text-align:left;"> http://broadinstitute.github.io/picard </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> dedup </td>
   <td style="text-align:left;"> 0.12.8 </td>
   <td style="text-align:left;"> Peltzer, et al.(2016) </td>
   <td style="text-align:left;"> https://github.com/apeltzer/DeDup </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Rescale damage</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> mapDamage2 </td>
   <td style="text-align:left;"> 2.2.1 </td>
   <td style="text-align:left;"> Jonsson, et al. (2013) </td>
   <td style="text-align:left;"> https://github.com/ginolhac/mapDamage </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Trim alignments</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> BamUtil </td>
   <td style="text-align:left;"> 1.0.15 </td>
   <td style="text-align:left;"> Jun, et al. (2015) </td>
   <td style="text-align:left;"> https://github.com/statgen/bamUtil </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Merge libraries</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> SAMtools </td>
   <td style="text-align:left;"> 1.14 </td>
   <td style="text-align:left;"> Danecek, et al. (2021) </td>
   <td style="text-align:left;"> https://github.com/samtools/samtools </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Realign indels</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> GATK IndelRealigner </td>
   <td style="text-align:left;"> 3.8 </td>
   <td style="text-align:left;"> DePristo, et al. (2011) </td>
   <td style="text-align:left;"> https://gatk.broadinstitute.org </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Recompute md flag</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> SAMtools </td>
   <td style="text-align:left;"> 1.14 </td>
   <td style="text-align:left;"> Danecek, et al. (2021) </td>
   <td style="text-align:left;"> https://github.com/samtools/samtools </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Imputation</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> GLIMPSE </td>
   <td style="text-align:left;"> 1.1.1 </td>
   <td style="text-align:left;"> Rubinacci, et al. (2021) </td>
   <td style="text-align:left;"> https://github.com/odelaneau/GLIMPSE </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> BCFtools </td>
   <td style="text-align:left;"> 1.15 </td>
   <td style="text-align:left;"> Danecek, et al. (2021) </td>
   <td style="text-align:left;"> https://github.com/samtools/bcftools </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Reports</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> FastQC </td>
   <td style="text-align:left;"> 0.11.9 </td>
   <td style="text-align:left;"> Andrews (2010) </td>
   <td style="text-align:left;"> https://www.bioinformatics.babraham.ac.uk/projects/fastqc </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Qualimap </td>
   <td style="text-align:left;"> 2.2.2d </td>
   <td style="text-align:left;"> Okonechnikov, et al. (2016) </td>
   <td style="text-align:left;"> http://qualimap.conesalab.org </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> MultiQC </td>
   <td style="text-align:left;"> 1.13 </td>
   <td style="text-align:left;"> Ewels, et al. (2016) </td>
   <td style="text-align:left;"> https://multiqc.info </td>
  </tr>
  <tr grouplength="1"><td colspan="4" style="border-bottom: 1px solid;"><strong>Statistics</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> BEDTools </td>
   <td style="text-align:left;"> 2.30.0 </td>
   <td style="text-align:left;"> Quinlan and Hall (2010) </td>
   <td style="text-align:left;"> https://github.com/arq5x/bedtools2 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> bamdamage </td>
   <td style="text-align:left;"> modified </td>
   <td style="text-align:left;"> Malaspinas, et al. (2014) </td>
   <td style="text-align:left;"> https://savannah.nongnu.org/projects/bammds </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> R </td>
   <td style="text-align:left;"> 4.0 </td>
   <td style="text-align:left;"> R Core Team (2022) </td>
   <td style="text-align:left;"> https://www.r-project.org </td>
  </tr>
</tbody>
</table>
