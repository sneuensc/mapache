## functions to generate indexes.
## Rules to index the reference check first if the indexes are next to the original reference
## file (input.orig). If present the file are just symlinked and no index is generated.


## indexing fasta
rule genome_index_bwa:
    """
    Indexing the genome for bwa
    """
    input:
        fasta="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        orig=lambda wildcards: recursive_get(["genome", wildcards.GENOME, "fasta"], ""),
    output:
        multiext(
            "{folder}/00_reference/{GENOME}/{GENOME}.fasta",
            ".sa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
        ),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12),
    params:
        params=recursive_get(["indexing", "bwa_params"], ""),
        cmd="f'bwa index {snakemake.params[0]} {snakemake.input.fasta} 2> {snakemake.log}'",
    log:
        "{folder}/00_reference/{GENOME}/{GENOME}_bwa_index.log",
    threads: 1
    conda:
        "../envs/bwa.yaml"
    envmodules:
        module_bwa,
    message:
        "--- BWA INDEX  {input.fasta}"
    script:
        "../scripts/fasta_indexing.py"


rule genome_index_bowtie2:
    """
    Indexing the genome for bowtie2
    """
    input:
        fasta="{folder}/{GENOME}.fasta",
        orig=lambda wildcards: recursive_get(["genome", wildcards.GENOME, "fasta"], ""),
    output:
        multiext(
            "{folder}/{GENOME}.fasta",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12),
    params:
        recursive_get(["indexing", "bowtie2_params"], ""),
        cmd="f'bowtie2-build {snakemake.params[0]} --threads {snakemake.threads} {snakemake.input.fasta} > {snakemake.log}'",
    log:
        "{folder}/{GENOME}_bowtie2_build.log",
    threads: 1
    conda:
        "../envs/bowtie2.yaml"
    envmodules:
        module_bowtie2,
    message:
        "--- BOWTIE2-BUILD {input.fasta}"
    script:
        "../scripts/fasta_indexing.py"


rule samtools_index_fasta:
    """
    Indexing the genome with samtools faidx
    """
    input:
        fasta="{folder}/{GENOME}.fasta",
        orig=lambda wildcards: recursive_get(["genome", wildcards.GENOME, "fasta"], ""),
    output:
        "{folder}/{GENOME}.fasta.fai",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12),
    params:
        recursive_get(["indexing", "samtools_params"], ""),
        cmd="f'samtools faidx {snakemake.params[0]}  {snakemake.input.fasta} > {snakemake.log}'",
    log:
        "{folder}/{GENOME}_samtools_faidx.log",
    threads: 1
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS_FAIDX {input.fasta}"
    script:
        "../scripts/fasta_indexing.py"


rule genome_index_picard:
    """
    Indexing the genome with Picard CreateSequenceDictionary
    """
    input:
        fasta="{folder}/{GENOME}.fasta",
        orig=lambda wildcards: recursive_get(["genome", wildcards.GENOME, "fasta"], ""),
    output:
        "{folder}/{GENOME}.dict",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12),
    params:
        picard_cmd=get_picard_bin(),
    #     cmd=get_picard_bin() + f' CreateSequenceDictionary --REFERENCE {input.fasta} --OUTPUT {output}',
    log:
        "{folder}/{GENOME}_picard_index.log",
    threads: 1
    conda:
        "../envs/picard.yaml"
    envmodules:
        module_picard,
    message:
        "--- PICARD CreateSequenceDictionary {input.fasta}"
    script:
        "../scripts/fasta_indexing.py"


## indexing bam
rule samtools_index_bam:
    """
    Index bam file with samtools
    """
    input:
        "{folder}/{file}.bam",
    output:
        "{folder}/{file}.bai",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 24),
    threads: 1
    log:
        "{folder}/{file}_samtools_index.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS INDEX {input}"
    shell:
        "samtools index {input} {output} 2> {log};"