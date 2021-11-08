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
