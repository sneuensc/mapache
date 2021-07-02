***!!! UNDER CONSTRUCTION !!!***

The pipeline is under strong development. Not all commits are working perfectly.

You are welcome to try it out. On the Wiki https://github.com/sneuensc/mapache/wiki, there is already a short manual describing how to get started. Happy to get any feedback/propositions.

***mapache*** is a lightweighted mapping pipeline for ancient DNA using the workflow manager *snakemake*

SUMMARY:

Ancient DNA is degraded and contaminated. Most standard bioinformatics tools to align sequenced reads have been designed for modern data and cannot be used “out of the box” to accommodate the features typical of ancient DNA. In this work, we propose a robust pipeline to align ancient DNA data and to have a first rough idea of the authenticity of the data. The implemented modules consist of the steps needed to go from a simple ‘fastq' file to a final ‘bam' file. The steps include quality control, mapping, filtering, duplicate removal, damage pattern inference, realignment, statistics to assess authenticity and inference of the sex of the organism. A final graphical report summarizes the statistics of the different modules allowing a quick overview of the data. The pipeline is implemented in the workflow manager snakemake providing the flexibility to run the pipeline on a workstation, a high memory server or a cluster. The installation of the pipeline and the underlying programs is made easy using conda. The pipeline may be used out of the box, or may be adapted with little knowledge of python and snakemake.


Following main steps are included in ***mapache***:
```
for each fastq file:
    1. adapter removal with AdapterRemoval (optional)
    2. mapping with bwa aln, bwa mem or bowtie2
for each library:
    3. marking duplicates with MarkDuplicates (optional)
    4. mapDamage is run, but by default not used for the next steps (optional)
for each sample:
    5. re-alignment of reads with GATK (optional)
    6. re-computation of the md flag with samtools (optional)
```

USAGE:
```
## recommended workflow on a server
snakemake dag --dag | dot -Tsvg > dag.svg                Visualization of the pipeline in a directed acyclic graph (DAG).
snakemake -n                                             Dry run
snakemake mapping                                        Run all mappings
snakemake stats                                          Run all statistics (and all mapping steps which are required to get them...)
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
