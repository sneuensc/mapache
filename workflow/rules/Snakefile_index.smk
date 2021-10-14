## functions to generate indexes.
## Rules to index the reference check first if the indexes are next to the original reference
## file (input.orig). If present the file are just symlinked and no index is generated.


rule check_chromosomes:
    """
    Check if the passed chromosome names are correct
    """
    input:
        fai="results/00_reference/{GENOME}/{GENOME}.fasta.fai",
    output:
        touch("results/00_reference/{GENOME}/{GENOME}.ok"),
    params:
        femaleChr=lambda wildcards: get_param3(
            "genome", wildcards.GENOME, "femaleChr", "X"
        ),
        maleChr=lambda wildcards: get_param3("genome", wildcards.GENOME, "maleChr", "Y"),
        mtChr=lambda wildcards: get_param3("genome", wildcards.GENOME, "mtChr", "MT"),
        autosomeChr=lambda wildcards: eval_to_list(
            get_param3("genome", wildcards.GENOME, "autosomeChr", "")
        ),
    log:
        "results/00_reference/{GENOME}/{GENOME}.ok.log",
    message:
        "--- TESTING CHROMOSOME NAMES"
    run:
        import numpy as np

        ## get all chromsome names from the reference genome
        allChr = pd.read_csv(input.fai, header=None, sep="\t")[0].tolist()
        if params.femaleChr not in allChr:
            set_param3("genome", wildcards.genome, "femaleChr", "")
            print(
                f"WARNING: In parameter 'genome:{wildcards.GENOME}:femaleChr' the argument '{params.mtChr}' is unknown, assuming no female chromosome!"
            )

        if params.maleChr not in allChr:
            set_param3("genome", wildcards.genome, "maleChr", "")
            print(
                f"WARNING: In parameter 'genome:{wildcards.GENOME}:maleChr' the argument '{params.mtChr}' is unknown, assuming no male chromosome!"
            )

        if params.mtChr not in allChr:
            set_param3("genome", wildcards.genome, "mtChr", "")
            print(
                f"WARNING: In parameter 'genome:{wildcards.GENOME}:mtChr' the argument '{params.mtChr}' is unknown, assuming no MT chromosome!"
            )

        if params.autosomeChr == "":
            autosomeChr = np.intersect1d(
                allChr, [params.femaleChr, params.maleChr, params.mtChr]
            )
            set_param3(
                "genome", wildcards.GENOME, "autosomeChr", list_to_csv(autosomeChr)
            )
        else:
            unknown = list(set(params.autosomeChr) - set(allChr))
            if len(unknown) > 0:
                autosomeChr = np.intersect1d(allChr, params.autosomeChr)
                set_param3(
                    "genome",
                    wildcards.GENOME,
                    "autosomeChr",
                    list_to_csv(autosomeChr),
                )
                if len(unknown) == 1:
                    print(
                        f"WARNING: In parameter 'genome:{wildcards.GENOME}:autosomeChr' the name {unknown} is unknown, ignoring it!"
                    )
                else:
                    print(
                        f"WARNING: In parameter 'genome:{wildcards.GENOME}:autosomeChr' the names {unknown} are unknown, ignoring them!"
                    )


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
