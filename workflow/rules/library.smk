##########################################################################################
## all rules for libraries
##########################################################################################


rule merge_bam_fastq2library:
    """
    Merge the bam files of the fastq step
    """
    input:
        get_bam_4_merge_bam_fastq2library,
    output:
        temp("{folder}/02_library/00_merged_fastq/01_bam/{sm}/{lb}.{genome}.bam"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
    threads: get_threads("merging", 4)
    log:
        "{folder}/02_library/00_merged_fastq/01_bam/{sm}/{lb}.{genome}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS MERGE {output}"
    shell:
        """
        samtools merge -f --threads {threads} {output} {input} 2> {log};
        """


rule merge_bam_low_qual_fastq2library:
    """
    Merge the bam files of the fastq step
    """
    input:
        get_bam_4_merge_bam_low_qual_fastq2library,
    output:
        temp(
            "{folder}/02_library/00_merged_fastq/01_bam_low_qual/{sm}/{lb}.{genome}.bam"
        ),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
    threads: get_threads("merging", 4)
    log:
        "{folder}/02_library/00_merged_fastq/01_bam_low_qual/{sm}/{lb}.{genome}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS MERGE {output}"
    shell:
        """
        samtools merge -f --threads {threads} {output} {input} 2> {log};
        """


rule markduplicates:
    """
    Remove duplicated mappings with pricards markduplicates
    """
    input:
        get_bam_4_markduplicates,
    output:
        bam=temp(
            "{folder}/02_library/01_duplicated/01_markduplicates/{sm}/{lb}.{genome}.bam"
        ),
        stats="{folder}/02_library/01_duplicated/01_markduplicates/{sm}/{lb}.{genome}.stats",
    resources:
        ## Java: there is an overhead: so reduce slightly the amount of memory given to the tool comapred to the job
        memory=lambda wildcards, attempt: get_memory_alloc(
            "remove_duplicates", attempt, 4
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "remove_duplicates", attempt, 24
        ),
    params:
        params=lambda wildcards: get_paramGrp(
            ["remove_duplicates", "params_markduplicates"],
            "--REMOVE_DUPLICATES true",
            wildcards,
        ),
        PICARD=get_picard_bin(),
    threads: get_threads("remove_duplicates", 4)
    log:
        "{folder}/02_library/01_duplicated/01_markduplicates/{sm}/{lb}.{genome}.log",
    conda:
        "../envs/picard.yaml"
    envmodules:
        module_picard,
    message:
        "--- MARKDUPLICATES {output.bam}"
    shell:
        """
        {params.PICARD} MarkDuplicates --INPUT {input} --OUTPUT {output.bam} --METRICS_FILE {output.stats} \
            {params.params} --ASSUME_SORT_ORDER coordinate --VALIDATION_STRINGENCY LENIENT 2> {log};
        """


rule dedup:
    """
    Remove duplicated mappings with dedup (-m only for collapsed PE reads!!!)
    """
    input:
        get_bam_4_markduplicates,
    output:
        json="{folder}/02_library/01_duplicated/01_dedup/{sm}/{lb}.{genome}.dedup.json",
        hist="{folder}/02_library/01_duplicated/01_dedup/{sm}/{lb}.{genome}.hist",
        log="{folder}/02_library/01_duplicated/01_dedup/{sm}/{lb}.{genome}.log",
        bam=temp("{folder}/02_library/01_duplicated/01_dedup/{sm}/{lb}.{genome}.bam"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc(
            "remove_duplicates", attempt, 4
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "remove_duplicates", attempt, 24
        ),
    params:
        params=lambda wildcards: get_paramGrp(
            ["remove_duplicates", "params_dedup"], "-m", wildcards
        ),
    log:
        "{folder}/02_library/01_duplicated/01_dedup/{sm}/{lb}.{genome}.log",
    conda:
        "../envs/dedup.yaml"
    envmodules:
        module_dedup,
    message:
        "--- DEDUP {output.bam}"
    shell:
        """
        bam={output.bam}
        dedup -i {input} $params  -o $(dirname $bam);
        mv ${{bam%%.bam}}_rmdup.bam $bam
        """


rule mapDamage_stats:
    """
    Run mapDamage to quantify the deamination pattern
    """
    input:
        ref="{folder}/00_reference/{genome}/{genome}.fasta",
        bam=get_bam_4_damage_rescale,
    output:
        #directory("{folder}/02_library/02_rescaled/01_mapDamage/{sm}/{lb}.{genome}_results_mapDamage"),
        deamination=report(
            "{folder}/02_library/02_rescaled/01_mapDamage/{sm}/{lb}.{genome}_results_mapDamage/Fragmisincorporation_plot.pdf",
            category="Damage pattern",
        ),
        length=report(
            "{folder}/02_library/02_rescaled/01_mapDamage/{sm}/{lb}.{genome}_results_mapDamage/Length_plot.pdf",
            category="Read length distribution",
        ),
    log:
        "{folder}/02_library/02_rescaled/01_mapDamage/{sm}/{lb}.{genome}_stats.log",
    threads: 1
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapdamage", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapdamage", attempt, 24),
    params:
        params=lambda wildcards: get_paramGrp(["mapdamage", "params"], "", wildcards),
    conda:
        "../envs/mapdamage.yaml"
    envmodules:
        module_mapdamage,
    message:
        "--- MAPDAMAGE {input.bam}"
    shell:
        """
        mapDamage -i {input.bam} -r {input.ref} -d $(dirname {output.deamination}) \
        {params.params} --merge-reference-sequences 2> {log};
        """


rule mapDamage_rescale:
    """
    Run mapDamage to rescale bam file
    """
    input:
        ref="{folder}/00_reference/{genome}/{genome}.fasta",
        bam=get_bam_4_damage_rescale,
        deamination="{folder}/02_library/02_rescaled/01_mapDamage/{sm}/{lb}.{genome}_results_mapDamage/Fragmisincorporation_plot.pdf",
    output:
        bam=temp("{folder}/02_library/02_rescaled/01_mapDamage/{sm}/{lb}.{genome}.bam"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapdamage", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapdamage", attempt, 24),
    log:
        "{folder}/02_library/02_rescaled/01_mapDamage/{sm}/{lb}.{genome}_rescale.log",
    threads: 1
    params:
        params=lambda wildcards: get_paramGrp(["mapdamage", "params"], "", wildcards),
    conda:
        "../envs/mapdamage.yaml"
    envmodules:
        module_mapdamage,
    message:
        "--- MAPDAMAGE RESCALE {output.bam}"
    shell:
        """
        mapDamage -i {input.bam} -r {input.ref} -d $(dirname {input.deamination}) \
        {params.params} --merge-reference-sequences --rescale-only --rescale-out {output} 2>> {log};
        """


## params cannot be empty
rule bamutil:
    "Run BamUtil to trim the end of reads in the BAM file"
    input:
        bam=get_bam_4_bamutil,
    output:
        bam=temp("{folder}/02_library/03_trim/01_bamutil/{sm}/{lb}.{genome}.bam"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("bamutil", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("bamutil", attempt, 24),
    log:
        "{folder}/02_library/03_trim/01_bamutil/{sm}/{lb}.{genome}.log",
    threads: 1
    params:
        lambda wildcards: get_paramGrp(["bamutil", "params"], "", wildcards),
    conda:
        "../envs/bamutil.yaml"
    envmodules:
        module_samtools,
        module_bamutil,
    message:
        "--- TRIM BAM {output.bam}"
    shell:
        """
        bam trimBam {input} {output} {params} 2>> {log};
        """
