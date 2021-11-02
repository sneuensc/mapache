## functions to generate indexes.
## Rules to index the reference check first if the indexes are next to the original reference
## file (input.orig). If present the file are just symlinked and no index is generated.


## indexing fasta
rule genome_index_bwa:
    """
    Indexing the genome for bwa
    """
    input:
        fasta="{folder}/{GENOME}.fasta",
        orig=lambda wildcards: get_param3("genome", wildcards.GENOME, "fasta", ""),
    output:
        multiext("{folder}/{GENOME}.fasta", ".sa", ".amb", ".ann", ".bwt", ".pac"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12),
    params:
        get_param2("indexing", "bwa_params", ""),
    log:
        "{folder}/bwa_index_{GENOME}.log",
    threads: 1
    conda:
        "../envs/bwa.yaml"
    envmodules:
        module_bwa,
    message:
        "--- BWA INDEX  {input.fasta}"
    script:
        "../scripts/bwa_indexing.py"


rule genome_index_bowtie2:
    """
    Indexing the genome for bowtie2
    """
    input:
        fasta="{folder}/{GENOME}.fasta",
        orig=lambda wildcards: get_param3("genome", wildcards.GENOME, "fasta", ""),
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
        get_param2("indexing", "bowtie2_params", ""),
    log:
        "{folder}/bowtie2_build_{GENOME}.log",
    threads: 1
    conda:
        "../envs/bowtie2.yaml"
    envmodules:
        module_bowtie2,
    message:
        "--- BOWTIE2-BUILD  {input.fasta}"
    script:
        "../scripts/bowtie2_indexing.py"


rule samtools_index_fasta:
    """
    Indexing the genome with samtools faidx
    """
    input:
        fasta="{folder}/{GENOME}.fasta",
        orig=lambda wildcards: get_param3("genome", wildcards.GENOME, "fasta", ""),
    output:
        "{folder}/{GENOME}.fasta.fai",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12),
    params:
        get_param2("indexing", "samtools_params", ""),
    log:
        "{folder}/samtools_fasta_{GENOME}.log",
    threads: 1
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS_FAIDX {input.fasta}"
    script:
        "../scripts/samtools_indexing.py"


rule genome_index_picard:
    """
    Indexing the genome with Picard CreateSequenceDictionary
    """
    input:
        fasta="{folder}/{GENOME}.fasta",
        orig=lambda wildcards: get_param3("genome", wildcards.GENOME, "fasta", ""),
    output:
        "{folder}/{GENOME}.dict",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12),
    log:
        "{folder}/picard_{GENOME}.log",
    threads: 1
    conda:
        "../envs/picard.yaml"
    envmodules:
        module_picard,
    message:
        "--- PICARD CreateSequenceDictionary {input.fasta} {input.orig}"
    script:
        "../scripts/picard_indexing.py"


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
        "{folder}/samtools_bam_{file}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS INDEX {input}"
    shell:
        "samtools index {input} {output} 2> {log};"
