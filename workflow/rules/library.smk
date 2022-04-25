##########################################################################################
## all rules for libraries
##########################################################################################
ruleorder: merge_bam_low_qual_fastq2library > merge_bam_fastq2library


def get_bam_4_merge_bam_fastq2library(wildcards):
    return [get_final_bam_fastq(wildcards.folder, wildcards.SM, wildcards.LB, ID, wildcards.GENOME) for ID in samples[wildcards.SM][wildcards.LB]]


rule merge_bam_fastq2library:
    """
    Merge the bam files of the fastq step
    """
    input:
        get_bam_4_merge_bam_fastq2library
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


rule merge_bam_low_qual_fastq2library:
    """
    Merge the bam files of the fastq step
    """
    input:
        mapped=lambda wildcards: expand(
            "{folder}/01_fastq/04_final_fastq/01_bam_low_qual/{SM}/{LB}/{ID}.{GENOME}.bam",
            ID=samples[wildcards.SM][wildcards.LB],
            allow_missing=True,
        ),
    output:
        "{folder}/02_library/00_merged_fastq/01_bam_low_qual/{SM}/{LB}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24),
    threads: get_threads("merging", 4)
    log:
        "{folder}/02_library/00_merged_fastq/01_bam_low_qual/{SM}/{LB}.{GENOME}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS MERGE {input}"
    script:
        "../scripts/merge_bam.py"


def get_bam_4_markduplicates(wildcards):
    bam = get_bam_4_merge_bam_fastq2library(wildcards)
    if len(bam) > 1: ## library consits of more than one fastq file: return 00_merged_fastq
        return f"{wildcards.folder}/02_library/00_merged_fastq/01_bam/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    else:              ## library consits of one fastq file: return return the location of the final library bam file
        return bam[0]


rule markduplicates:
    """
    Remove duplicated mappings with pricards markduplicates
    """
    input:
        get_bam_4_markduplicates,
    output:
        bam=temp("{folder}/02_library/01_duplicated/01_markduplicates/{SM}/{LB}.{GENOME}.bam"),
        stats="{folder}/02_library/01_duplicated/01_markduplicates/{SM}/{LB}.{GENOME}.stats",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc(
            "remove_duplicates", attempt, 4
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "remove_duplicates", attempt, 24
        ),
    params:
        params=recursive_get(
            ["remove_duplicates", "params_markduplicates"], "--REMOVE_DUPLICATES true"
        ),
        PICARD=get_picard_bin(),
    threads: get_threads("remove_duplicates", 4)
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


rule dedup:
    """
    Remove duplicated mappings with dedup (-m only for collapsed PE reads!!!)
    """
    input:
        get_bam_4_markduplicates,
    output:
        json="{folder}/02_library/01_duplicated/01_dedup/{SM}/{LB}.{GENOME}.dedup.json",
        hist="{folder}/02_library/01_duplicated/01_dedup/{SM}/{LB}.{GENOME}.hist",
        log="{folder}/02_library/01_duplicated/01_dedup/{SM}/{LB}.{GENOME}.log",
        bam=temp("{folder}/02_library/01_duplicated/01_dedup/{SM}/{LB}.{GENOME}.bam"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc(
            "remove_duplicates", attempt, 4
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "remove_duplicates", attempt, 24
        ),
    params:
        params=recursive_get(["remove_duplicates", "params_dedup"], "-m"),
        collapsed=collapse
        and "nan"
        not in [
            i["Data2"]
            for i in recursive_get(
                [wildcards.SM, wildcards.LB, wildcards.ID, "Data2"],
                [],
                my_dict=samples,
            ).values()
        ],
    log:
        "{folder}/02_library/01_duplicated/01_dedup/{SM}/{LB}.{GENOME}.log",
    conda:
        "../envs/dedup.yaml"
    envmodules:
        module_dedup,
    message:
        "--- DEDUP {input}"
    shell:
        """
        ## remove -m or --merged from $params (needed for SE or not collapsed PE reads)
        if [ {params.collapsed} ]; then
            params={params.params}
        else
            params=$(echo {params.params} | sed  's/ -m/ /g' | sed  's/ --merged/ /g')
        fi

        bam={output.bam}
        dedup -i {input} $params  -o $(dirname $bam);
        mv ${{bam%%.bam}}_rmdup.bam $bam
        """

def get_bam_4_damage(wildcards):
    if remove_duplicates == "markduplicates":
        bam = f"{wildcards.folder}/02_library/01_duplicated/01_markduplicates/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    elif remove_duplicates == "dedup":
        bam = f"{wildcards.folder}/02_library/01_duplicated/01_dedup/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    else:
        bam = get_bam_4_markduplicates(wildcards)
    return bam   

rule mapDamage_stats:
    """
    Run mapDamage to quantify the deamination pattern
    """
    input:
        ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        bam=get_bam_4_damage,
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
        bam=get_bam_4_damage,
        deamination="{folder}/02_library/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}_results_mapDamage/Fragmisincorporation_plot.pdf",
    output:
        bam=temp("{folder}/02_library/02_rescaled/01_mapDamage/{SM}/{LB}.{GENOME}.bam"),
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
