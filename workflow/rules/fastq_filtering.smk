##########################################################################################
## sorting


rule samtools_sort:
    """
    Sort bam file with samtools
    """
    input:
        get_bam_4_sorting,
    output:
        temp("{folder}/01_fastq/02_mapped/03_bam_sort/{sm}/{lb}/{id}.{genome}.bam"),
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("sorting", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("sorting", attempt, 24),
    log:
        "{folder}/01_fastq/02_mapped/03_bam_sort/{sm}/{lb}/{id}.{genome}.log",
    threads: get_threads("sorting", 4)
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS SORT {input}"
    shell:
        """
        samtools sort --threads {threads} {input} > {output} 2> {log}
        """


##########################################################################################
## filtering

if save_low_qual:

    rule samtools_filter:
        """
        Filter mappings following quality and keeping the low quality mappings
        """
        input:
            "{folder}/01_fastq/02_mapped/03_bam_sort/{sm}/{lb}/{id}.{genome}.bam",
        output:
            mapped=temp(
                "{folder}/01_fastq/03_filtered/01_bam_filter/{sm}/{lb}/{id}.{genome}.bam"
            ),
            low_qual=temp(
                "{folder}/01_fastq/03_filtered/01_bam_filter_low_qual/{sm}/{lb}/{id}.{genome}.bam"
            ),
        params:
            lambda wildcards: get_paramGrp(
                ["filtering", "params"], "-F 4 -q 30", wildcards
            ),
        resources:
            mem_mb=lambda wildcards, attempt: get_memory_alloc("filtering", attempt, 4),
            runtime=lambda wildcards, attempt: get_runtime_alloc(
                "filtering", attempt, 24
            ),
        log:
            "{folder}/01_fastq/03_filtered/01_bam_filter/{sm}/{lb}/{id}.{genome}.log",
        threads: get_threads("filtering", 4)
        conda:
            "../envs/samtools.yaml"
        envmodules:
            module_samtools,
        message:
            "--- SAMTOOLS FILTER {input}"
        shell:
            """
            samtools view -b --threads {threads} {params} \
            -U {output.low_qual} {input} > {output.mapped} 2> {log}
            """

    rule gam_filter:
        """
        Filter mappings following quality and keeping the low quality mappings
        """
        input:
            "{folder}/01_fastq/02_mapped/02_vg_giraffe/{sm}/{lb}/{id}.{genome}.gam"
        output:
            mapped=temp(
                "{folder}/01_fastq/03_filtered/01_gam_filter/{sm}/{lb}/{id}.{genome}.gam"
            ),
            low_qual=temp(
                "{folder}/01_fastq/03_filtered/01_gam_filter_low_qual/{sm}/{lb}/{id}.{genome}.gam"
            ),
        params:
            lambda wildcards: get_paramGrp(
                ["filtering", "params_gam"], "-q 30", wildcards
            ),
            bin = get_param(["software", "vg"], "vg")
        resources:
            mem_mb=lambda wildcards, attempt: get_memory_alloc("filtering", attempt, 4),
            runtime=lambda wildcards, attempt: get_runtime_alloc(
                "filtering", attempt, 24
            ),
        log:
            "{folder}/01_fastq/03_filtered/01_gam_filter/{sm}/{lb}/{id}.{genome}.log",
        threads: get_threads("filtering", 4)
        conda:
            "../envs/samtools.yaml"
        envmodules:
            module_samtools,
        message:
            "--- VG FILTER {input}"
        shell:
            """
            (
                {params.bin} filter {input} --threads {threads} {params} > {output.mapped} &&
                {params.bin} filter {input} --threads {threads} -U {params} > {output.low_qual} 
            ) 2> {log}
            """

else:

    rule samtools_filter:
        """
        Filter mappings following quality and discard low quality mappings
        """
        input:
            "{folder}/01_fastq/02_mapped/03_bam_sort/{sm}/{lb}/{id}.{genome}.bam",
        output:
            mapped=temp(
                "{folder}/01_fastq/03_filtered/01_bam_filter/{sm}/{lb}/{id}.{genome}.bam"
            ),
        params:
            lambda wildcards: get_paramGrp(
                ["filtering", "params"], "-F 4 -q 30", wildcards
            ),
        resources:
            mem_mb=lambda wildcards, attempt: get_memory_alloc("filtering", attempt, 4),
            runtime=lambda wildcards, attempt: get_runtime_alloc(
                "filtering", attempt, 24
            ),
        log:
            "{folder}/01_fastq/03_filtered/01_bam_filter/{sm}/{lb}/{id}.{genome}.log",
        threads: get_threads("filtering", 4)
        message:
            "--- SAMTOOLS FILTER {input}"
        shell:
            """
            samtools view -b --threads {threads} {params} {input} > {output.mapped} 2> {log}
            """

    rule gam_filter:
        """
        Filter mappings following quality and discard low quality mappings
        """
        input:
            "{folder}/01_fastq/02_mapped/02_vg_giraffe/{sm}/{lb}/{id}.{genome}.gam"
        output:
            mapped=temp(
                "{folder}/01_fastq/03_filtered/01_gam_filter/{sm}/{lb}/{id}.{genome}.gam"
            ),
        params:
            lambda wildcards: get_paramGrp(
                ["filtering", "params"], "-q 30", wildcards
            ),
            bin = get_param(["software", "vg"], "vg")
        resources:
            mem_mb=lambda wildcards, attempt: get_memory_alloc("filtering", attempt, 4),
            runtime=lambda wildcards, attempt: get_runtime_alloc(
                "filtering", attempt, 24
            ),
        log:
            "{folder}/01_fastq/03_filtered/01_gam_filter/{sm}/{lb}/{id}.{genome}.log",
        threads: get_threads("filtering", 4)
        message:
            "--- VG FILTER {input}"
        shell:
            """
            {params.bin} filter {input} --threads {threads} {params} > {output.mapped} 2> {log}
            """
