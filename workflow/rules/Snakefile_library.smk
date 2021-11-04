##########################################################################################
## all rules for libraries
##########################################################################################


rule merge_bam_fastq2library:
    """
    Merge the bam files of the fastq step
    """
    input:
        lambda wildcards: [
            f"results/01_fastq/04_final_fastq/01_bam/{wildcards.SM}/{wildcards.LB}/{ID}.{wildcards.GENOME}.bam"
            for ID in samples[wildcards.SM][wildcards.LB]
        ],
    output:
        "{folder}/00_merged_fastq/01_bam/{SM}/{LB}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
    threads: get_threads("merging", 4)
    log:
        "{folder}/00_merged_fastq/01_bam/{SM}/{LB}.{GENOME}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS MERGE {input}"
    script:
        "../scripts/merge_files.py"


rule merge_bam_fastq2library_low_qual:
    """
    Merge the bam files of the fastq step
    """
    input:
        lambda wildcards: [
            f"results/01_fastq/04_final_fastq/01_bam_low_qual/{wildcards.SM}/{wildcards.LB}/{ID}.{wildcards.GENOME}.bam"
            for ID in samples[wildcards.SM][wildcards.LB]
        ],
    output:
        "{folder}/00_merged_fastq/01_bam_low_qual/{SM}/{LB}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
    threads: get_threads("merging", 4)
    log:
        "{folder}/00_merged_fastq/01_bam_low_qual/{SM}/{LB}.{GENOME}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS MERGE {input}"
    script:
        "../scripts/merge_files.py"


rule remove_duplicates:
    """
    Remove duplicated mappings
    """
    input:
        "{folder}/00_merged_fastq/01_bam/{SM}/{LB}.{GENOME}.bam",
    output:
        bam="{folder}/01_duplicated/01_rmdup/{SM}/{LB}.{GENOME}.bam",
        stats="{folder}/01_duplicated/01_rmdup/{SM}/{LB}.{GENOME}.stats",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("markduplicates", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "markduplicates", attempt, 24
        ),
    params:
        rmdup_params=get_param2("markduplicates", "params", "REMOVE_DUPLICATES=true"),
        PICARD=get_param2("software", "picard_jar", "picard.jar"),
    threads: get_threads("markduplicates", 4)
    log:
        "{folder}/01_duplicated/01_rmdup/{SM}/{LB}.{GENOME}.log",
    conda:
        "../envs/picard.yaml"
    envmodules:
        module_picard,
    message:
        "--- MARKDUPLICATES {input}"
    shell:
        """
        ## get binary
        jar={params.PICARD}
        if [ "${{jar: -4}}" == ".jar" ]; then
            bin="java -XX:ParallelGCThreads={threads} -XX:+UseParallelGC -XX:-UsePerfData \
                    -Xms{resources.memory}m -Xmx{resources.memory}m -jar {params.PICARD}"
           else
               bin={params.PICARD}
           fi

           ## run MarkDuplicates
        $bin MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.stats} \
            {params.rmdup_params} ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT 2> {log};
        """


rule samtools_extract_duplicates:
    """
    Extract duplicates of bam file
    """
    input:
        "{folder}/01_duplicated/01_rmdup/{SM}/{LB}.{GENOME}.bam",
    output:
        mapped="{folder}/01_duplicated/01_rmdup/{SM}/{LB}.{GENOME}_mapped.bam",
        dup="{folder}/01_duplicated/01_rmdup/{SM}/{LB}.{GENOME}_duplicates.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("markduplicates", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "markduplicates", attempt, 24
        ),
    threads: get_threads("markduplicates", 4)
    log:
        "{folder}/01_duplicated/01_rmdup/{SM}/{LB}.{GENOME}_extract_duplicates.log",
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


rule mapDamage_stats:
    """
    Run mapDamage to quantify the deamination pattern
    """
    input:
        ref="results/00_reference/{GENOME}/{GENOME}.fasta",
        bam=get_mapDamage_bam,
    output:
        #directory("{folder}/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_results_mapDamage"),
        deamination=report(
            "{folder}/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_results_mapDamage/Fragmisincorporation_plot.pdf",
            category="Damage pattern",
        ),
        length=report(
            "{folder}/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_results_mapDamage/Length_plot.pdf",
            category="Read length distribution",
        ),
    log:
        "{folder}/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_stats.log",
    threads: 1
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapdamage", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapdamage", attempt, 24),
    params:
        params=get_param2("mapdamage", "params", ""),
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
        ref="results/00_reference/{GENOME}/{GENOME}.fasta",
        bam=get_mapDamage_bam,
        deamination="{folder}/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_results_mapDamage/Fragmisincorporation_plot.pdf",
    output:
        bam="{folder}/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapdamage", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapdamage", attempt, 24),
    log:
        "{folder}/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_rescale.log",
    threads: 1
    params:
        params=get_param2("mapdamage", "params", ""),
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
        "{folder}/03_final_library/01_bam/{SM}/{LB}.{GENOME}.bam",
    message:
        "--- GET FINAL BAM {input} (LIBRARY LEVEL)"
    log:
        "{folder}/03_final_library/01_bam/{SM}/{LB}.{GENOME}.bam.log",
    run:
        symlink_rev(input, output)


rule get_final_library_low_qual:
    """
    Get the final bam file of the library part
    """
    input:
        "{folder}/00_merged_fastq/01_bam_low_qual/{SM}/{LB}.{GENOME}.bam",
    output:
        "{folder}/03_final_library/01_bam_low_qual/{SM}/{LB}.{GENOME}.bam",
    message:
        "--- GET FINAL LOW_QUAL BAM {input} (LIBRARY LEVEL)"
    log:
        "{folder}/03_final_library/01_bam_low_qual/{SM}/{LB}.{GENOME}.bam.log",
    run:
        symlink_rev(input, output)


rule get_final_library_lduplicate:
    """
    Get the final duplicate bam file of the library part
    """
    input:
        "{folder}/01_duplicated/01_rmdup/{SM}/{LB}.{GENOME}_duplicates.bam",
    output:
        "{folder}/03_final_library/01_bam_duplicate/{SM}/{LB}.{GENOME}.bam",
    message:
        "--- GET FINAL DUPLICATE BAM {input} (LIBRARY LEVEL)"
    log:
        "{folder}/03_final_library/01_bam_duplicate/{SM}/{LB}.{GENOME}.bam.log",
    run:
        symlink_rev(input, output)
