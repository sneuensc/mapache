# Preparation

In the config/config.yaml file you can find the parameters and options specified for each of the rules of mapache. These can be changed depending on specific needs. 
Running mapache requires that the path(s) to the reference genome(s) to be specified in the config/config.yaml file. For example, for hg19 human reference genome: 

```
genome: 
  hg19:             ## prefix for the reference genome (free to choose)
    fasta: path_to_reference/foo.fasta            ## reference genome (required)
```

It is possible to run the pipeline with a config/config.yaml file containing only the path to the reference genome and no other parameters or options, in which case it runs with the default parameters and options.  


# Sample file
Specify the location of the samples in the sample file, `config/samples.txt`, by default.

```
sample_file: config/samples.txt
```

The sample file is a table listing the relationship between `fastq files`, `libraries` and `samples`. The order of the columns is free, but the column names have to be exact. **ID, LB and SM names should not contain points ('.')**:

sample file for single-end libraries:
```
    ID           Data                        MAPQ  LB    PL        SM
    a_L2_R1_001  reads/a_L2_R1_001.fastq.gz  30    a_L2  ILLUMINA  ind1
    a_L2_R1_002  reads/a_L2_R1_002.fastq.gz  30    a_L2  ILLUMINA  ind1
    b_L2_R1_001  reads/b_L2_R1_001.fastq.gz  30    b_L2  ILLUMINA  ind1
    b_L2_R1_002  reads/b_L2_R1_002.fastq.gz  30    b_L2  ILLUMINA  ind1
```

sample file for paired-end libraries (or mix of PE and SE libraries):
```
    ID           Data1                       Data2                       MAPQ  LB    PL        SM
    a_L2_R1_001  reads/a_L2_R1_001.fastq.gz  reads/a_L2_R2_001.fastq.gz  30    a_L2  ILLUMINA  ind1
    a_L2_R1_002  reads/a_L2_R1_002.fastq.gz  reads/a_L2_R2_001.fastq.gz  30    a_L2  ILLUMINA  ind1
    b_L2_R1_001  reads/b_L2_R1_001.fastq.gz  reads/a_L2_R2_001.fastq.gz  30    b_L2  ILLUMINA  ind1
    b_L2_R1_002  reads/b_L2_R1_002.fastq.gz  NULL                        30    b_L2  ILLUMINA  ind1

```
