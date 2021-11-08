##########################################################################################
## all rules for the samples
##########################################################################################
localrules:
    get_final_bam,  ## executed locally on a cluster


ruleorder: merge_bam_library2sample_low_qual > get_final_bam


rule merge_bam_library2sample:
    """
    Merge the bam files of the library step
    """
    input:
        mapped=lambda wildcards: [
            f"results/02_library/03_final_library/01_bam/{wildcards.SM}/{LB}.{wildcards.GENOME}.bam"
            for LB in samples[wildcards.SM]
        ],
    output:
        "{folder}/00_merged_library/01_bam/{SM}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
    threads: get_threads("merging", 4)
    log:
        "{folder}/00_merged_library/01_bam/{SM}.{GENOME}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS MERGE merge_bam_library2sample {input}"
    script:
        "../scripts/merge_files.py"


rule merge_bam_library2sample_low_qual:
    """
    Merge the low quality bam files of the library step
    """
    input:
        lambda wildcards: [
            f"results/02_library/03_final_library/01_bam_low_qual/{wildcards.SM}/{LB}.{wildcards.GENOME}.bam"
            for LB in samples[wildcards.SM]
        ],
    output:
        "{folder}/00_merged_library/01_bam_low_qual/{SM}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
    threads: get_threads("merging", 4)
    log:
        "{folder}/00_merged_library/01_bam_low_qual/{SM}.{GENOME}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS MERGE merge_bam_library2sample_low_qual {input}"
    script:
        "../scripts/merge_files.py"


rule merge_bam_library2sample_duplicates:
    """
    Merge the duplicate bam files of the library step
    """
    input:
        low_qual=lambda wildcards: get_bams_of_sample_low_cov(
            "library_rmdup",
            "_duplicates.bam",
            wildcards.SM,
            wildcards.GENOME,
        ),
    output:
        "{folder}/00_merged_library/01_bam_duplicate/{SM}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
    threads: get_threads("merging", 4)
    log:
        "{folder}/00_merged_library/01_bam_duplicate/{SM}.{GENOME}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS MERGE merge_bam_library2sample_duplicates {input}"
    script:
        "../scripts/merge_files.py"


rule realign:
    """
    Realign sequence around indels.
    """
    input:
        ref="results/00_reference/{GENOME}/{GENOME}.fasta",
        fai="results/00_reference/{GENOME}/{GENOME}.fasta.fai",
        dict="results/00_reference/{GENOME}/{GENOME}.dict",
        bam="{folder}/00_merged_library/01_bam/{SM}.{GENOME}.bam",
        bai="{folder}/00_merged_library/01_bam/{SM}.{GENOME}.bai",
    output:
        bam="{folder}/01_realigned/01_realign/{SM}.{GENOME}.bam",
        intervals="{folder}/01_realigned/01_realign/{SM}.{GENOME}.intervals",
        bai="{folder}/01_realigned/01_realign/{SM}.{GENOME}.bai",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("realign", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("realign", attempt, 24),
    threads: get_threads("realign", 4)
    params:
        GATK=recursive_get(["software", "gatk3_jar"], "GenomeAnalysisTK.jar"),
    log:
        "{folder}/01_realigned/01_realign/{SM}.{GENOME}.log",
    conda:
        "../envs/gatk3.yaml"
    envmodules:
        module_gatk3,
    message:
        "--- GATK INDELREALIGNER {input.bam}"
    shell:
        """
        ## get binary
         jar={params.GATK};
        if [ "${{jar: -4}}" == ".jar" ]; then
            bin="java -Djava.io.tmpdir=/tmp/ -XX:ParallelGCThreads={threads} -XX:+UseParallelGC \
                -XX:-UsePerfData -Xms15000m -Xmx15000m -jar {params.GATK}"
           else
               bin={params.GATK}
           fi

           ## run GATK
           $bin -I {input.bam} -R {input.ref} -T RealignerTargetCreator -o {output.intervals} 2> {log}; \
        $bin -I {input.bam} -T IndelRealigner -R {input.ref} -targetIntervals \
                {output.intervals} -o {output.bam} 2>> {log};         
        """


rule samtools_calmd:
    """
    Recompute the md flag.
    """
    input:
        ref="results/00_reference/{GENOME}/{GENOME}.fasta",
        bam=get_md_flag_bam,
    output:
        "results/03_sample/02_md_flag/01_md_flag/{SM}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("calmd", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("calmd", attempt, 24),
    threads: get_threads("calmd", 4)
    log:
        "results/03_sample/02_md_flag/01_md_flag/{SM}.{GENOME}.log",
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


rule get_final_bam:
    """
    Get the final bam files 
    """
    input:
        get_final_bam,
    output:
        "{folder}/03_final_sample/01_bam/{SM}.{GENOME}.bam",
    threads: 1
    log:
        "{folder}/03_final_sample/01_bam/{SM}.{GENOME}.bam.log",
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


rule move_final_bam_duplicate:
    """
    Get the final bam files
    """
    input:
        "{folder}/00_merged_library/01_bam_duplicate/{SM}.{GENOME}.bam",
    output:
        "{folder}/03_final_sample/01_bam_duplicate/{SM}.{GENOME}.bam",
    threads: 1
    log:
        "{folder}/03_final_sample/01_bam_duplicate/{SM}.{GENOME}.bam.log",
    message:
        "--- SIMLINKK FINAL DUPLICATE BAM"
    run:
        symlink_rev(input, output)


##########################################################################################
