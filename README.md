Readme is under construction.

!!! Please consult the Wiki on how to install and use MAPache !!!

How to use **MAPache**

Do not continue here


==========================================================================


This pipeline maps (ancient) DNA libraries to the (human) genome as follows.
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
REQUIRED:
-	to be adapted for the current project
    -   Config file (config.yaml): 
        This file contains all settings and has to be adapted to the data.
    -   Sample file (samples.txt):
        This file lists all dependencies between libraries and samples and has to be generated for each run.
-	scripts needed (nothing to be adapted)
	-   Snakefile
	    This file contains the main rules to launch the pipeline.
	-   scripts/
	    This folder contains scripts used for the mapping pipeline
	-   report/
	    This folder contains files used for the reporting.
-	axiom/
	this folder is used to submit the pipeline to the Axiom cluster:
		- config.yaml
		- scheduler.py
		- slurm.yaml
		- status.py
		- submit.py


All these files can be cloned from the git lab at "https://github.com/sneuensc/snakemake_aDNA_mapping.git".
```
git clone https://github.com/sneuensc/snakemake_aDNA_mapping.git
```


USAGE:
```
## requirement on wally/axiom - snakemake has to be installed
module load Bioinformatics/Software/vital-it;
module add Development/snakemake/5.17.0;


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

## correct report for rules with 0s runtime (figure for runtimes are not correctly shown due to log scale)
sed -i "s/'runtime': 0.0/'runtime': 0.1/g" report.html
sed -i.bk "s/'runtime': 0.0/'runtime': 0.1/g" report.html    ## for MacOS


## temporal files
--notemp, --nt        Ignore temp() declarations. This is useful when
                        running only a part of the workflow, since temp()
                        would lead to deletion of probably needed files by
                        other parts of the workflow.
--delete-temp-output  Remove all temporary files generated by the workflow.
                        Use together with --dry-run to list files without
                        actually deleting anything. Note that this will not
                        recurse into subworkflows.
```
TYPICAL USAGE
```bash
snakemake --configfile config_wally.yaml -p --nt --rerun-incomplete --profile wally #runs the pipeline on wally, using the config_wally.yaml printing out the commands without removing temporal files reruning incomplete steps (this should actually happen by default) 		
```

GENERAL PARAMETERS:

*sample file*

The sample file lists the mapping information in the following format. The order of the columns is free, but the column names have to be exact. **ID, LB and SM names should not contain points ('.')**:

single-end libraries:
```
    ID           Data                        MAPQ  LB    PL        SM
    a_L2_R1_001  reads/a_L2_R1_001.fastq.gz  30    a_L2  ILLUMINA  ind1
    a_L2_R1_002  reads/a_L2_R1_002.fastq.gz  30    a_L2  ILLUMINA  ind1
    b_L2_R1_001  reads/b_L2_R1_001.fastq.gz  30    b_L2  ILLUMINA  ind1
    b_L2_R1_002  reads/b_L2_R1_002.fastq.gz  30    b_L2  ILLUMINA  ind1
```

paired-end libraries (or mix of PE and SE libraries):
```
    ID           Data1                       Data2                       MAPQ  LB    PL        SM
    a_L2_R1_001  reads/a_L2_R1_001.fastq.gz  reads/a_L2_R2_001.fastq.gz  30    a_L2  ILLUMINA  ind1
    a_L2_R1_002  reads/a_L2_R1_002.fastq.gz  reads/a_L2_R2_001.fastq.gz  30    a_L2  ILLUMINA  ind1
    b_L2_R1_001  reads/b_L2_R1_001.fastq.gz  reads/a_L2_R2_001.fastq.gz  30    b_L2  ILLUMINA  ind1
    b_L2_R1_002  reads/b_L2_R1_002.fastq.gz  NULL                        30    b_L2  ILLUMINA  ind1

```

The following code snippet allows to genrate  the sample file
```
readFolder=../reads
echo -e "ID\tData\tMAPQ\tLB\tPL\tSM" > samples.txt
paste <(ls $readFolder/*gz | rev | cut -d'/' -f1 | rev | cut -d'.' -f1) \
      <(ls $readFolder/*gz) \
      <(ls $readFolder/*gz | rev | cut -d'/' -f1 | rev | cut -d'_' -f-3) \
      <(ls $readFolder/*gz | rev | cut -d'/' -f1 | rev | cut -d'_' -f1 | sed 's/./&0/3') | \
awk '{print $1"_fq",$2,"30",$3"_lb","ILLUMINA",$4}' OFS='\t' >> samples.txt
```

*GENOME*

The path to the genome fasta file (for bwa) and to the bowtie index (for bowtie2).
```
SAMPLES: sample file [samples.txt]
GENOME:  reference genome (example: hg19: /archive/unibe/eg/amalaspi/group/genomes/reference_human/hs.build37.1/hs.build37.1.fa)
```


**adapter_removal (optional)**

Run AdapterRemoval which allows to
- remove adapters (default),
- remove low quality forgoing and trailing bases (default),
- remove a fixed number of forgoing and trailing bases (parameters --trim5p and --trim3p)
- collapse paired-end reads (parameter --collapse)

```
run_adapter_removal:        should AdapterRemoval be run [True]
adaptrem_params:            any additional parameter [--minlength 30]
threads_adaptrem:           number of threads [4]
```


**mapping**

The following mappers are available (default bwa_aln)
Available mappers:
- bwa_aln
- bwa_mem
- bowtie2

```
mapper:                     which mapper to use [bwa_aln] (available: bwa_aln, bwa_mem, bowtie2)
bwa_aln_params:             any additional parameter [-l 1024]
bwa_samse_params:           any additional parameter [-n 3]
bwa_mem_params:             any additional parameter []
bowtie2_params:             any additional parameter []
threads_mapping:            number of threads [4]
```


**samtools_sort**

Bam files are sorted.
```
threads_sort:               number of threads [4]
```


**samtools_filter**

Bam files are filtered for mapping quality (qualiuty is taken formthe sample file. 
Non-mapping reads and low quality mappings are stored in a separate bam file (nothing is lost).
```
threads_filter:             number of threads [4]
```


**merge_bam_fastq2library**

Merge the fastq based bam files to the library level. If only a single bam file is present it is symlinked.
```
threads_merge: number of threads [4]
```


**remove_duplicates (optional)**

Remove duplicates or store them in a separate bam file.
```
run_remove_duplicates:      should RemoveDuplicates be run [True] 
rmdup_params:               any additional parameter [REMOVE_DUPLICATES=true]
threads_rmdup:              number of threads [4]

## to store the duplicates in a separate bam file:
## in this case parameter 'rmdup_params' has NOT to be set to 'REMOVE_DUPLICATES=true'
run_extract_duplicates:     True
threads_extract_duplicates: 4
```


**mapDamage_stats and mapDamage_rescale (optional)**

Compute deamination pattern (default) and rescale the bam file if desired (not done by default).
```
run_mapDamage_rescale:      should the bam file be rescaled [False]
mapDamage_params:           any additional parameter []
```


**merge_bam_library2sample**

Merge the library based bam files to the sample level. If only a single one is present it is symlinked.
```
threads_merge:              number of threads [4] (same as for rule merge_bam_fastq2library)
```


**realign and samtools_index (optional)**

Index the bam file realing reads around indels
```
run_realign:                should GATK IndelRealigner be run [True]

```

**samtools_calmd (optional)**

Re-compute the md flag
```
run_compute_md:             should samtools calmd be run [True]
threads_calmd:              number of threads used  [4]
```

**stats**
These parameters allow to define which statistics should be computed and the time and memory allocation
```
run_damage: 'False'         infer the deamination pattern ['False']; Options: 'False' (or anything), 'mapDamage', 'bamdamage'
run_depth: True             run the script deptyh.py and report coverage statistics [False]
run_bammds: False           run bammds (not yet implemented)

## bamdamage 
bamdamage_fraction: 0       fraction/number of alignments to consider to infer the damage pattern
                            0: use all (default)
                            <1:  a fraction of the total number of alignments to consider 
                            >=1: absolute number of alignments to consider
                            
## memory (default 2GB)
fastqc_mem
samtools_flagstat_mem
multiqc_mem
depth_mem

## runtime (default 1h)
fastqc_time
samtools_flagstat_time
multiqc_time
depth_time
```

**Resource allocations for cluster**

If a job fails it can be re-launched automatically. The number of trials can be defined by the parameter 
*--restart-times* (In the profile axiom it is set to 3 trails). If the job fails due to resource allocations 
it is wise to increase the resource allocations iteratively for each trail:

*Memory allocation*

The following parameters may be used to specify the memory allocation. All allocations are in GB.
```
## memory allocation
foo_mem:                    memory allocation in GB for the inital run (default 4)

## specifyint increment for memory allocation
foo_mem_increment:          amount of memory to increase by each failure of the job (default same as foo_mem)
memory_increment_ratio:     this is a global parameter and specifies the ratio of the inital memory allocation 'foo_mem' 
                            to use as increment if the increment is not directly specifies with 'foo_mem_increment':
                            foo_mem_increment = memory_increment_ratio * foo_mem
                            set to 1: foo_mem_increment = foo_mem  (default)
                            set to 0: foo_mem_increment = 0        (no change)
```

*Runtime allocation*

The following parameters may be used to specify the runtime allocation. All allocations are in hours.
```
## memory allocation
foo_time:                   runtime allocation in hours for the inital run (default 12)

## specifyint increment for runtime allocation
foo_time_increment:         runtime in hours to increase by each failure of the job (default same as foo_time)
runtime_increment_ratio:    this is a global parameter and specifies the ratio of the inital runtime allocation 'foo_time' 
                            to use as increment if the increment is not directly specifies with 'foo_time_increment':
                            foo_time_increment = runtime_increment_ratio * foo_time
                            set to 1: foo_time_increment = foo_time  (default)
                            set to 0: foo_time_increment = 0         (no change)
```


