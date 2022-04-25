##########################################################################################
## all rules for the samples
##########################################################################################
localrules:
    get_final_bam,

def get_bam_4_merge_bam_library2sample(wildcards):
    return [get_final_bam_library(wildcards.folder, wildcards.SM, LB, wildcards.GENOME) for LB in samples[wildcards.SM]]


if FINAL_BAM_FOLDER_SAMPLE == f"{RESULT_DIR}/03_sample/00_merged_library":
    rule merge_bam_library2sample:
        """
        Merge the bam files of the library step
        """
        input:
            get_bam_4_merge_bam_library2sample 
        output:
            "{folder}/03_sample/00_merged_library/{type}/{SM}.{GENOME}.bam"
        resources:
            memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
            runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
        threads: get_threads("merging", 4)
        log:
            "{folder}/03_sample/00_merged_library/{type}/{SM}.{GENOME}.log",
        conda:
            "../envs/samtools.yaml"
        envmodules:
            module_samtools,
        message:
            "--- SAMTOOLS MERGE merge_bam_library2sample {input}"
        script:
            "../scripts/merge_bam.py"
else:
    rule merge_bam_library2sample:
        """
        Merge the bam files of the library step
        """
        input:
            get_bam_4_merge_bam_library2sample
        output:
            "{folder}/03_sample/00_merged_library/{type}/{SM}.{GENOME}.bam"
        resources:
            memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
            runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
        threads: get_threads("merging", 4)
        log:
            "{folder}/03_sample/00_merged_library/{type}/{SM}.{GENOME}.log",
        conda:
            "../envs/samtools.yaml"
        envmodules:
            module_samtools,
        message:
            "--- SAMTOOLS MERGE merge_bam_library2sample {input}"
        script:
            "../scripts/merge_bam.py"


rule merge_bam_low_qual_library2sample:
    """
    Merge the bam files of the library step
    """
    input:
        mapped=lambda wildcards: expand(
            "{folder}/02_library/03_final_library/01_bam_low_qual/{SM}/{LB}.{GENOME}.bam",
            LB=samples[wildcards.SM],
            allow_missing=True,
        ),    
    output:
        "{folder}/03_sample/00_merged_library/01_bam_low_qual/{SM}.{GENOME}.bam",
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
        "--- SAMTOOLS MERGE merge_bam_library2sample {input}"
    script:
        "../scripts/merge_bam.py"

def get_bam_4_realign(wildcards):
    bam = get_bam_4_merge_bam_library2sample(wildcards)
    if len(bam) > 1: ## sample consits of more than one library return 00_merged_library
        return f"{wildcards.folder}/02_library/00_merged_library/01_bam/{wildcards.SM}.{wildcards.GENOME}.bam"
    else:              ## library consits of one fastq file: return return the location of the final library bam file
        return bam[0]

## get the corresponding bai file:
def get_bai_4_realign(wildcards):
    bam = get_bam_4_realign(wildcards)
    return f"{bam[:len(bam) - 4]}.bai"


if FINAL_BAM_FOLDER_SAMPLE == f"{RESULT_DIR}/03_sample/01_realigned/01_realign":
    rule realign:
        """
        Realign sequence around indels.
        """
        input:
            ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
            fai="{folder}/00_reference/{GENOME}/{GENOME}.fasta.fai",
            dict="{folder}/00_reference/{GENOME}/{GENOME}.dict",
            bam=get_bam_4_realign,
            bai=get_bai_4_realign,
        output:
            bam="{folder}/03_sample/01_realigned/01_realign/{SM}.{GENOME}.bam",
            intervals=temp("{folder}/03_sample/01_realigned/01_realign/{SM}.{GENOME}.intervals"),
            bai="{folder}/03_sample/01_realigned/01_realign/{SM}.{GENOME}.bai",
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
            "--- GATK INDELREALIGNER {input.bam}"
        shell:
            """
            {params.GATK} -I {input.bam} -R {input.ref} -T RealignerTargetCreator -o {output.intervals} 2> {log}; \
            {params.GATK} -I {input.bam} -T IndelRealigner -R {input.ref} -targetIntervals \
                    {output.intervals} -o {output.bam} 2>> {log};         
            """
else:
    rule realign:
        """
        Realign sequence around indels.
        """
        input:
            ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
            fai="{folder}/00_reference/{GENOME}/{GENOME}.fasta.fai",
            dict="{folder}/00_reference/{GENOME}/{GENOME}.dict",
            bam=get_bam4realign,
            bai=get_bai4realign,
        output:
            bam=temp("{folder}/03_sample/01_realigned/01_realign/{SM}.{GENOME}.bam"),
            intervals=temp("{folder}/03_sample/01_realigned/01_realign/{SM}.{GENOME}.intervals"),
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
            "--- GATK INDELREALIGNER {input.bam}"
        shell:
            """
            {params.GATK} -I {input.bam} -R {input.ref} -T RealignerTargetCreator -o {output.intervals} 2> {log}; \
            {params.GATK} -I {input.bam} -T IndelRealigner -R {input.ref} -targetIntervals \
                    {output.intervals} -o {output.bam} 2>> {log};         
            """


def get_bam_4_samtools_calmd(wildcards):
    if run_realign:
        return f"{wildcards.folder}/03_sample/01_realigned/01_realign/{wildcards.SM}.{wildcards.GENOME}.bam"
    else:
        return get_bam_4_realign(wildcards)


rule samtools_calmd:
    """
    Recompute the md flag.
    """
    input:
        ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        bam=get_md_flag_bam,
    output:
        "{folder}/03_sample/02_md_flag/01_md_flag/{SM}.{GENOME}.bam",
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
        "--- SAMTOOLS CALMD {input.bam}"
    shell:
        """
        samtools calmd --threads {threads} {input.bam} {input.ref} 2> {log} | samtools view -bS - > {output}
        """


def get_bam_4_final_bam(wildcards):
    if run_compute_md:
        return f"{wildcards.folder}/03_sample/02_md_flag/01_md_flag/{wildcards.SM}.{wildcards.GENOME}.bam"
    else:
        return get_bam_4_samtools_calmd(wildcards)


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
        "--- SIMLINKK FINAL BAM"
    run:
        symlink_rev(input, output)


rule get_final_bam_low_qual:
    """
    Get the final bam files 
    """
    input:
        "{folder}/00_merged_library/01_bam_low_qual/{SM}.{GENOME}.bam",
    output:
        "{folder}/03_final_sample/01_bam_low_qual/{SM}.{GENOME}.bam",
    threads: 1
    log:
        "{folder}/03_final_sample/01_bam_low_qual/{SM}.{GENOME}.bam.log",
    message:
        "--- SIMLINKK FINAL LOW_QUAL BAM"
    run:
        symlink_rev(input, output)


##########################################################################################
