# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs. 

To run ``mapache`` at least the following two parameters have to be set in the file ``config/config.yaml``:
* a reference genome (e.g., hg19 human reference genome):
```
genome: 
  hg19:                                  ## prefix for the reference genome (free to choose)
    fasta: path_to_reference/foo.fasta   ## file location
```
* the path to the sample file ``samples.tsv``:
```
sample_file: config/samples.tsv
```


The sample file ``samples.tsv`` lists all ``fastq files`` and their relationship to ``library`` and ``sample`` in a tab seperated format. For single-end libraries the format is:
|SM   |LB      |Data                           |ID             |
|:----|:-------|:------------------------------|:--------------|
|ind1 |lib1_lb |lib1_R1_001.fastq.gz |lib1_R1_001_fq |
|ind1 |lib1_lb |lib1_R1_002.fastq.gz |lib1_R1_002_fq |
|ind1 |lib2_lb |lib2_R1_001.fastq.gz |lib2_R1_001_fq |
|ind1 |lib2_lb |lib2_R1_002.fastq.gz |lib2_R1_002_fq |
|ind2 |lib3_lb |lib3_R1_001.fastq.gz |lib3_R1_001_fq |
|ind2 |lib3_lb |lib3_R1_002.fastq.gz |lib3_R1_002_fq |

where
* **SM**: sample name
* **LB**: library name
* **ID**: uniq identifier for the fastq file (row)
* **Data**: path to the ``fastq file`` 

For more details, please read the [Wiki](https://github.com/sneuensc/mapache/wiki). The code is available on [GitHub](https://github.com/sneuensc/mapache).

# Run mapache
After configuering the files ``config/config.yaml`` and the ``samples.tsv`` you can run ``mapache`` using the command
```
snakemake -j1
``` 