##########################################################################################
## all rules for the samples
##########################################################################################
localrules:
    get_final_bam,


rule merge_bam_library2sample:
    """
    Merge the bam files of the library step
    """
    input:
        get_bam_4_merge_bam_library2sample,
    output:
        temp("{folder}/03_sample/00_merged_library/01_bam/{SM}.{GENOME}.bam"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
    threads: get_threads("merging", 4)
    log:
        "{folder}/03_sample/00_merged_library/01_bam/{SM}.{GENOME}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS MERGE merge_bam_library2sample {output}"
    script:
        "../scripts/merge_bam.py"


rule merge_bam_low_qual_library2sample:
    """
    Merge the bam files of the library step
    """
    input:
        get_bam_4_merge_bam_low_qual_library2sample,
    output:
        temp("{folder}/03_sample/00_merged_library/01_bam_low_qual/{SM}.{GENOME}.bam"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
    threads: get_threads("merging", 4)
    log:
        "{folder}/03_sample/00_merged_library/01_bam_low_qual/{SM}.{GENOME}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS MERGE merge_bam_library2sample {output}"
    script:
        "../scripts/merge_bam.py"


rule realign:
    """
    Realign sequence around indels.
    """
    input:
        ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        fai="{folder}/00_reference/{GENOME}/{GENOME}.fasta.fai",
        dict="{folder}/00_reference/{GENOME}/{GENOME}.dict",
        bam=get_bam_4_realign,
        bai=lambda wildcards: bam2bai(get_bam_4_realign(wildcards)),
    output:
        bam=temp("{folder}/03_sample/01_realigned/01_realign/{SM}.{GENOME}.bam"),
        intervals=temp(
            "{folder}/03_sample/01_realigned/01_realign/{SM}.{GENOME}.intervals"
        ),
        bai=temp("{folder}/03_sample/01_realigned/01_realign/{SM}.{GENOME}.bai"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("realign", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("realign", attempt, 24),
    threads: get_threads("realign", 4)
    params:
        GATK=get_gatk_bin(),
    log:
        "{folder}/03_sample/01_realigned/01_realign/{SM}.{GENOME}.log",
    conda:
        "../envs/gatk3.yaml"
    envmodules:
        module_gatk3,
    message:
        "--- GATK INDELREALIGNER {output.bam}"
    shell:
        """
        {params.GATK} -I {input.bam} -R {input.ref} -T RealignerTargetCreator -o {output.intervals} 2> {log}; \
        {params.GATK} -I {input.bam} -T IndelRealigner -R {input.ref} -targetIntervals \
                {output.intervals} -o {output.bam} 2>> {log};         
        """


rule samtools_calmd:
    """
    Recompute the md flag.
    """
    input:
        ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        bam=get_bam_4_samtools_calmd,
    output:
        temp("{folder}/03_sample/02_md_flag/01_md_flag/{SM}.{GENOME}.bam"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("calmd", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("calmd", attempt, 24),
    threads: get_threads("calmd", 4)
    log:
        "{folder}/03_sample/02_md_flag/01_md_flag/{SM}.{GENOME}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS CALMD {output}"
    shell:
        """
        samtools calmd --threads {threads} {input.bam} {input.ref} 2> {log} | samtools view -bS - > {output}
        """


rule get_final_bam:
    """
    Get the final bam files 
    """
    input:
        get_bam_4_final_bam,
    output:
        "{folder}/03_sample/03_final_sample/01_bam/{SM}.{GENOME}.bam",
    threads: 1
    log:
        "{folder}/03_sample/03_final_sample/01_bam/{SM}.{GENOME}.bam.log",
    message:
        "--- SIMLINKK FINAL BAM {output}"
    run:
        symlink_rev(input, output)


rule get_final_bam_low_qual:
    """
    Get the final bam files 
    """
    input:
        get_bam_4_final_bam_low_qual,
    output:
        "{folder}/03_final_sample/01_bam_low_qual/{SM}.{GENOME}.bam",
    threads: 1
    log:
        "{folder}/03_final_sample/01_bam_low_qual/{SM}.{GENOME}.bam.log",
    message:
        "--- SIMLINKK FINAL LOW_QUAL BAM {output}"
    run:
        symlink_rev(input, output)


##########################################################################################
