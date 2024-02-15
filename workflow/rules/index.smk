##########################################################################################
## all rules to generate indexes
##########################################################################################
## Rules to index the reference check first if the indexes are next to the original reference
## file (input.orig). If present the file are just symlinked and no index is generated.


## indexing fasta
rule genome_index_bwa:
    """
    Indexing the genome for bwa
    """
    input:
        fasta="{folder}/00_reference/{genome}/{genome}.fasta",
        orig=lambda wildcards: get_param(["genome", wildcards.genome], ""),
    output:
        temp(
            multiext(
                "{folder}/00_reference/{genome}/{genome}.fasta",
                ".sa",
                ".amb",
                ".ann",
                ".bwt",
                ".pac",
            )
        ),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12),
    params:
        params=get_param(["indexing", "bwa_params"], ""),
        cmd="f'bwa index {snakemake.params[0]} {snakemake.input.fasta} 2> {snakemake.log}'",
    log:
        "{folder}/00_reference/{genome}/{genome}_bwa_index.log",
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
        fasta="{folder}/{genome}.fasta",
        orig=lambda wildcards: get_param(["genome", wildcards.genome], ""),
    output:
        temp(
            multiext(
                "{folder}/{genome}.fasta",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            )
        ),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12),
    params:
        get_param(["indexing", "bowtie2_params"], ""),
        cmd="f'bowtie2-build {snakemake.params[0]} --threads {snakemake.threads} {snakemake.input.fasta} {snakemake.input.fasta}> {snakemake.log}'",
    log:
        "{folder}/{genome}_bowtie2_build.log",
    threads: 1
    conda:
        "../envs/bowtie2.yaml"
    envmodules:
        module_bowtie2,
    message:
        "--- BOWTIE2-BUILD {input.fasta} {input.fasta}"
    script:
        "../scripts/fasta_indexing.py"

rule genome_index_bowtie2_long:
    """
    Indexing the genome for bowtie2
    """
    input:
        fasta="{folder}/{genome}.fasta",
        orig=lambda wildcards: get_param(["genome", wildcards.genome], ""),
    output:
        temp(
            multiext(
                "{folder}/{genome}.fasta",
                ".1.bt2l",
                ".2.bt2l",
                ".3.bt2l",
                ".4.bt2l",
                ".rev.1.bt2l",
                ".rev.2.bt2l",
            )
        ),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12),
    params:
        get_param(["indexing", "bowtie2_params"], ""),
        cmd="f'bowtie2-build {snakemake.params[0]} --threads {snakemake.threads} {snakemake.input.fasta} {snakemake.input.fasta} > {snakemake.log}'",
    log:
        "{folder}/{genome}_bowtie2_build.log",
    threads: 1
    conda:
        "../envs/bowtie2.yaml"
    envmodules:
        module_bowtie2,
    message:
        "--- BOWTIE2-BUILD {input.fasta} {input.fasta}"
    script:
        "../scripts/fasta_indexing.py"


rule samtools_index_fasta:
    """
    Indexing the genome with samtools faidx
    """
    input:
        fasta="{folder}/{genome}.fasta",
        orig=lambda wildcards: get_param(["genome", wildcards.genome], ""),
    output:
        temp("{folder}/{genome}.fasta.fai"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12),
    params:
        get_param(["indexing", "samtools_params"], ""),
        cmd="f'samtools faidx {snakemake.params[0]}  {snakemake.input.fasta} > {snakemake.log}'",
    log:
        "{folder}/{genome}_samtools_faidx.log",
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
        fasta="{folder}/{genome}.fasta",
        orig=lambda wildcards: get_param(["genome", wildcards.genome], ""),
    output:
        temp("{folder}/{genome}.dict"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12),
    params:
        java_mem_overhead_factor = float(get_param(["software", "java_mem_overhead_factor"], 0.2)),
        cmd="f'picard CreateSequenceDictionary -Xmx{round(snakemake.resources.memory * (1.0 - snakemake.params.java_mem_overhead_factor))}M --REFERENCE {snakemake.input.fasta} --OUTPUT {snakemake.output}'",
    log:
        "{folder}/{genome}_picard_index.log",
    threads: 1
    conda:
        "../envs/picard.yaml"
    envmodules:
        module_picard,
    message:
        "--- PICARD CreateSequenceDictionary {input.fasta}"
    script:
        "../scripts/fasta_indexing.py"


ruleorder: samtools_index_bam > samtools_index_bam_temp


## indexing bam
rule samtools_index_bam_temp:
    """
    Index bam file with samtools
    """
    input:
        "{file}.bam",
    output:
        temp("{file}.bai"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 24),
    threads: 1
    log:
        "{file}_samtools_index.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS INDEX {input}"
    shell:
        "samtools index {input} {output} 2> {log};"


## same es the previous one, but output is not temp
rule samtools_index_bam:
    """
    Index last bam file with samtools (should remain)
    """
    input:
        "{folder}/03_sample/03_final_sample/{type}/{sm}.{genome}.bam",
    output:
        "{folder}/03_sample/03_final_sample/{type}/{sm}.{genome}.bai",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 24),
    threads: 1
    log:
        "{folder}/03_sample/03_final_sample/{type}/{sm}.{genome}_samtools_index.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS INDEX {input}"
    shell:
        "samtools index {input} {output} 2> {log};"
