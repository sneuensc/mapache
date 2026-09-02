
<img src="images/mapache_logo.png" width="500">


Welcome to the wiki page of mapache! 🦝 🎉 

[mapache](https://github.com/sneuensc/mapache) is a pipeline designed to map sequencing reads to one or more reference genomes. 
Initially, the pipeline was aiming at mapping ancient DNA obtained from humans. Along the way, mapache has become a more general-purpose mapping pipeline.

In this page, you will find some of the nomenclature used to map sequencing reads with mapache. We hope that the users get familiar with terms such as "sample", "library", "FASTQ" used in mapache. More precise definitions can be found [here](https://github.com/sneuensc/mapache/wiki/2.-Sample-file).

If you are already familiar with these terms, you can try to [install and run the mapache tutorial.](https://github.com/sneuensc/mapache/wiki/4.-Tutorial:-Running-mapache) 

If you are using mapache, please cite
> Neuenschwander, Samuel, Diana I. Cruz Dávalos, Lucas Anchieri, Bárbara Sousa da Mota, Davide Bozzi, Simone Rubinacci, Olivier Delaneau, Simon Rasmussen, and Anna-Sapfo Malaspinas. 2023. “Mapache: A Flexible Pipeline to Map Ancient DNA.” Bioinformatics , 2023. https://doi.org/10.1093/bioinformatics/btad028.

as well as the software used to produce the mapping and/or imputed files (listed at the bottom of the README file in https://github.com/sneuensc/mapache).

# The problem(s)

You've got one or more samples (`SM`).

For each sample, one or more libraries (`LB`) were built.

Each of the libraries was sequenced once or more.

In each sequencing run, you received the `Data` in one or more FASTQ (`FASTQ`) files.

Before having your morning coffee ☕ , you organised this information nicely enough in a table like this:

|SM   |LB      |Data                           |ID             |
|:----|:-------|:------------------------------|:--------------|
|ind1 |lib1_lb |lib1_R1_001.fastq.gz |lib1_R1_001_fq |
|ind1 |lib1_lb |lib1_R1_002.fastq.gz |lib1_R1_002_fq |
|ind1 |lib2_lb |lib2_R1_001.fastq.gz |lib2_R1_001_fq |
|ind1 |lib2_lb |lib2_R1_002.fastq.gz |lib2_R1_002_fq |
|ind2 |lib3_lb |lib3_R1_001.fastq.gz |lib3_R1_001_fq |
|ind2 |lib3_lb |lib3_R1_002.fastq.gz |lib3_R1_002_fq |


Now, you have to map these data to one or more genomes (`genome`). Easy-PC, it's just a matter of trimming adapters, aligning to the reference(s), filtering for mapping quality, sorting, indexing, merging the right files with proper IDs, identify and remove duplicated reads, realign around indels, and fix the md flag. Do not forget to create a nice summary of what happened during the mapping. Now you can go ahead with your project.

👂🏽 you wanted to mark duplicates instead of removing them?

Easy, it's just a matter of finding the files right before you removed duplicates. Then identify and mark duplicated reads, realign around indels, fix the md flag. Do not forget to create a nice summary of what happened during the mapping. Now you can go ahead with your project.

👂🏽 the adapters you specified at the beginning were wrong?

Easy, it's just a matter of setting the right adapters, trimming adapters, aligning to the reference(s), filtering for mapping quality, sorting, indexing, merging the right files with proper IDs, identify and mark (remember: do not remove this time) duplicated reads, realign around indels, fix the md flag. Do not forget to create a nice summary of what happened during the mapping. Now you can go ahead with your project.

👂🏽 you've got more money to sequence more of that amazing sample?

Good to hear! You may repeat the steps above and merge the corresponding files, but also update your report.

However, if you sequenced more of the same library, it is better if you merge the BAM files before marking duplicates. Then it is just a matter of identifying and marking duplicated reads, realign around indels, fix the md flag. Do not forget to update that nice summary of what happened during the mapping. Now you can go ahead with your project.


# mapache's solution
We use the workflow management system [snakemake](https://snakemake.readthedocs.io/en/stable/) to automatise the mapping pipeline and generate statistics related to the mapping.

You can use an input table similar to the one shown above, and a config file to specify the parameters used at each of the mapping steps.

👂🏽 you wanted to mark duplicates instead of removing them?

Alright. Please adapt the config file and re-launch mapache.

```
# example of a section in the config file
markduplicates:
    run: True  
    params: 'REMOVE_DUPLICATES=false'
    save_duplicates: mark  ## remove, mark, extract
```


👂🏽 the adapters you specified at the beginning were wrong?

Just set the right adaptors in the config file and re-launch mapache. All steps that need to be re-run will be re-run. Take a break 🏓 

```
# example of a section in the config file
adapterremoval:
    run: True
    params: ' --adapter1 ACCCGTTTGTAAACTT --adapter2 TTTGCCCAATCCGATC'
```

👂🏽 you've got more money to sequence more of that amazing sample?

Congratulations! Just extend the input file with the info of the new data and re-launch mapache. All steps that need to be re-run will be re-run.


|SM   |LB      |Data                           |ID             |
|:----|:-------|:------------------------------|:--------------|
|ind1 |lib1_lb |lib1_R1_001.fastq.gz |lib1_R1_001_fq |
|ind1 |lib1_lb |lib1_R1_002.fastq.gz |lib1_R1_002_fq |
|ind1 |lib2_lb |lib2_R1_001.fastq.gz |lib2_R1_001_fq |
|ind1 |lib2_lb |lib2_R1_002.fastq.gz |lib2_R1_002_fq |
|ind2 |lib3_lb |lib3_R1_001.fastq.gz |lib3_R1_001_fq |
|ind2 |lib3_lb |lib3_R1_002.fastq.gz |lib3_R1_002_fq |
|ind2 |lib3_lb |lib3_R1_003.fastq.gz |I_got_money |

What about the summary of my mappings?

Every time you launch mapache, the tables with statistics are updated. One of such tables would look like this:

|genome |SM   | reads_raw| reads_trim| trim_prop| mapped_raw| duplicates| duplicates_prop| mapped_unique| length_reads_raw| length_reads_trimmed| length_mapped_raw| length_mapped_unique| endogenous_raw| endogenous_unique|Sex                                      | read_depth|
|:------|:----|---------:|----------:|---------:|----------:|----------:|---------------:|-------------:|----------------:|--------------------:|-----------------:|--------------------:|--------------:|-----------------:|:----------------------------------------|----------:|
|hg19   |ind1 |      7485|       7485|         1|         57|         10|       0.1754386|            47|         47.37214|             47.37214|          41.75439|             42.40426|      0.0076152|         0.0062792|XX |      6e-07|
|hg19   |ind2 |      3672|       3672|         1|         21|          4|       0.1904762|            17|         47.00000|             47.00000|          38.38095|             38.88235|      0.0057190|         0.0046296|XY |      2e-07|


# Workflow
The image below shows a graph for an example mapping pipeline, including the steps to get the mapping statistics.
Each node is a "rule" or a step that needs to be executed. Rules are connected if the output of one rule is used as input of another rule, and the direction is indicated by the edges.


[[images/rulegraph.png| ALT TEXT|width=700px]]

The image above gives you a simplified or summarised view of the steps that will run.
Note that, for example, the rule to remove adapter ("adapter_removal_se", green bubble on the top right) actually needs to be executed once for every input file.

The extended graph for an example dataset with five input FASTQ files would look like this: 

For the initial input table that we had, the directed acyclic graph with the rules to be run would look like [this](https://github.com/sneuensc/mapache/wiki/images/dag.png): [[images/dag.png| a dag]]

