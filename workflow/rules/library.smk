##########################################################################################
## all rules for libraries
##########################################################################################


rule merge_bam_fastq2library:
    """
    Merge the bam files of the fastq step
    """
    input:
        mapped=lambda wildcards: expand(
            "{folder}/01_fastq/04_final_fastq/{type}/{SM}/{LB}/{ID}.{GENOME}.bam",
            ID=samples[wildcards.SM][wildcards.LB],
            allow_missing=True,
        ),
    output:
        "{folder}/02_library/00_merged_fastq/{type}/{SM}/{LB}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
    threads: get_threads("merging", 4)
    log:
        "{folder}/02_library/00_merged_fastq/{type}/{SM}/{LB}.{GENOME}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS MERGE {input}"
    script:
        "../scripts/merge_bam.py"


rule markduplicates:
    """
    Remove duplicated mappings with pricards markduplicates
    """
    input:
        "{folder}/02_library/00_merged_fastq/01_bam/{SM}/{LB}.{GENOME}.bam",
    output:
        bam="{folder}/02_library/01_duplicated/01_markduplicates/{SM}/{LB}.{GENOME}.bam",
        stats="{folder}/02_library/01_duplicated/01_markduplicates/{SM}/{LB}.{GENOME}.stats",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("remove_duplciates", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "remove_duplciates", attempt, 24
        ),
    params:
        params=recursive_get(["remove_duplciates", "params_markduplicates"], "--REMOVE_DUPLICATES true"),
        PICARD=get_picard_bin(),
    threads: get_threads("remove_duplciates", 4)
    log:
        "{folder}/02_library/01_duplicated/01_markduplicates/{SM}/{LB}.{GENOME}.log",
    conda:
        "../envs/picard.yaml"
    envmodules:
        module_picard,
    message:
        "--- MARKDUPLICATES {input}"
    shell:
        """
        {params.PICARD} MarkDuplicates --INPUT {input} --OUTPUT {output.bam} --METRICS_FILE {output.stats} \
            {params.params} --ASSUME_SORT_ORDER coordinate --VALIDATION_STRINGENCY LENIENT 2> {log};
        """


rule markduplicates_extract:
    """
    Extract duplicates of bam file
    """
    input:
        "{folder}/02_library/01_duplicated/01_markduplicates/{SM}/{LB}.{GENOME}.bam",
    output:
        mapped="{folder}/02_library/01_duplicated/01_markduplicates/{SM}/{LB}.{GENOME}_mapped.bam",
        dup="{folder}/02_library/01_duplicated/01_markduplicates/{SM}/{LB}.{GENOME}_duplicates.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("remove_duplciates", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "remove_duplciates", attempt, 24
        ),
    threads: get_threads("remove_duplciates", 4)
    log:
        "{folder}/02_library/01_duplicated/01_markduplicates/{SM}/{LB}.{GENOME}_extract_duplicates.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS EXTRACT DUPLICATES {input}"
    shell:
        """
        samtools view -b --threads {threads} -F 1024 -U {output.dup} {input} > {output.mapped} 2> {log}
        """


rule dedup:
    """
    Remove duplicated mappings with dedup (only for collapsed PE reads!!!)
    """
    input:
        "{folder}/02_library/00_merged_fastq/01_bam/{SM}/{LB}.{GENOME}.bam",
    output:
        json="{folder}/02_library/01_duplicated/01_dedup/{SM}/{LB}.{GENOME}.dedup.json",
        hist="{folder}/02_library/01_duplicated/01_dedup/{SM}/{LB}.{GENOME}.hist",
        log="{folder}/02_library/01_duplicated/01_dedup/{SM}/{LB}.{GENOME}.log",
        bam="{folder}/02_library/01_duplicated/01_dedup/{SM}/{LB}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("remove_duplciates", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "remove_duplciates", attempt, 24
        ),
    params:
        params=recursive_get(["markduplicates", "params"], "--REMOVE_DUPLICATES true"),
    #threads: get_threads("remove_duplciates", 1)
    log:
        "{folder}/02_library/01_duplicated/01_dedup/{SM}/{LB}.{GENOME}.log",
    conda:
        "../envs/dedup.yaml"
    envmodules:
        module_picard,
    message:
        "--- DEDUP {input}"
    shell:
        """
        bam={output$bam}
        dedup -i {input} -m  -o $(dirname $bam);
        mv ${bam%%.bam}_rmdup.bam $bam
        """


rule mapDamage_stats:
    """
    Run mapDamage to quantify the deamination pattern
    """
    input:
        ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        bam=get_mapDamage_bam,
    output:
        #directory("{folder}/02_library/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_results_mapDamage"),
        deamination=report(
            "{folder}/02_library/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_results_mapDamage/Fragmisincorporation_plot.pdf",
            category="Damage pattern",
        ),
        length=report(
            "{folder}/02_library/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_results_mapDamage/Length_plot.pdf",
            category="Read length distribution",
        ),
    log:
        "{folder}/02_library/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_stats.log",
    threads: 1
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapdamage", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapdamage", attempt, 24),
    params:
        params=recursive_get(["mapdamage", "params"], ""),
    conda:
        "../envs/mapdamage.yaml"
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
        ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        bam=get_mapDamage_bam,
        deamination="{folder}/02_library/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_results_mapDamage/Fragmisincorporation_plot.pdf",
    output:
        bam="{folder}/02_library/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapdamage", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapdamage", attempt, 24),
    log:
        "{folder}/02_library/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_rescale.log",
    threads: 1
    params:
        params=recursive_get(["mapdamage", "params"], ""),
    conda:
        "../envs/mapdamage.yaml"
    message:
        "--- MAPDAMAGE {input.bam}"
    shell:
        """
        mapDamage -i {input.bam} -r {input.ref} -d $(dirname {input.deamination}) \
        {params.params} --merge-reference-sequences --rescale-only --rescale-out {output} 2>> {log};
        """


##########################################################################################
rule get_final_library:
    """
    Get the final bam file of the library part
    """
    input:
        get_final_bam_library,
    output:
        "{folder}/02_library/03_final_library/{type}/{SM}/{LB}.{GENOME}.bam",
    message:
        "--- GET FINAL BAM {input} (LIBRARY LEVEL)"
    log:
        "{folder}/02_library/03_final_library/{type}/{SM}/{LB}.{GENOME}.bam.log",
    run:
        symlink_rev(input, output)
