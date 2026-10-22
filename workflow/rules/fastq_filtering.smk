##########################################################################################
## sorting


rule samtools_sort:
    """
    Sort bam file with samtools
    """
    input:
        get_bam_4_sorting,
    output:
        maybe_temp("{folder}/01_fastq/02_mapped/03_bam_sort/{sm}/{lb}/{id}{chunk}.{genome}.bam", "sorting"),
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("sorting", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("sorting", attempt, 24),
    log:
        "{folder}/01_fastq/02_mapped/03_bam_sort/{sm}/{lb}/{id}{chunk}.{genome}.log",
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
            "{folder}/01_fastq/02_mapped/03_bam_sort/{sm}/{lb}/{id}{chunk}.{genome}.bam",
        output:
            mapped=maybe_temp(
                "{folder}/01_fastq/03_filtered/01_bam_filter/{sm}/{lb}/{id}{chunk}.{genome}.bam", "filtering"
            ),
            low_qual=maybe_temp(
                "{folder}/01_fastq/03_filtered/01_bam_filter_low_qual/{sm}/{lb}/{id}{chunk}.{genome}.bam", "filtering"
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
            "{folder}/01_fastq/03_filtered/01_bam_filter/{sm}/{lb}/{id}{chunk}.{genome}.log",
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
            "{folder}/01_fastq/02_mapped/02_vg_giraffe/{sm}/{lb}/{id}{chunk}.{genome}.gam"
        output:
            mapped=maybe_temp(
                "{folder}/01_fastq/03_filtered/01_gam_filter/{sm}/{lb}/{id}{chunk}.{genome}.gam", "filtering"
            ),
            low_qual=maybe_temp(
                "{folder}/01_fastq/03_filtered/01_gam_filter_low_qual/{sm}/{lb}/{id}{chunk}.{genome}.gam", "filtering"
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
            "{folder}/01_fastq/03_filtered/01_gam_filter/{sm}/{lb}/{id}{chunk}.{genome}.log",
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
            "{folder}/01_fastq/02_mapped/03_bam_sort/{sm}/{lb}/{id}{chunk}.{genome}.bam",
        output:
            mapped=maybe_temp(
                "{folder}/01_fastq/03_filtered/01_bam_filter/{sm}/{lb}/{id}{chunk}.{genome}.bam", "filtering"
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
            "{folder}/01_fastq/03_filtered/01_bam_filter/{sm}/{lb}/{id}{chunk}.{genome}.log",
        threads: get_threads("filtering", 4)
        conda:
            "../envs/samtools.yaml"
        envmodules:
            module_samtools,
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
            "{folder}/01_fastq/02_mapped/02_vg_giraffe/{sm}/{lb}/{id}{chunk}.{genome}.gam"
        output:
            mapped=maybe_temp(
                "{folder}/01_fastq/03_filtered/01_gam_filter/{sm}/{lb}/{id}{chunk}.{genome}.gam", "filtering"
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
            "{folder}/01_fastq/03_filtered/01_gam_filter/{sm}/{lb}/{id}{chunk}.{genome}.log",
        threads: get_threads("filtering", 4)
        message:
            "--- VG FILTER {input}"
        shell:
            """
            {params.bin} filter {input} --threads {threads} {params} > {output.mapped} 2> {log}
            """

##########################################################################################
## get final fastq files (if it was split)

rule final_fastq:
    """
    Get final bam file (if it was split)
    {kind}: bam, bam_low_qual, gam or gam_low_qual; 
    {ext}: bam (samtools merge) or gam (cat)
    """
    input:
        chunks = get_chunk_filename_of_fastq_split,
        bam = get_split_4_final,
    output:
        maybe_temp("{folder}/01_fastq/04_final/01_{kind}/{sm}/{lb}/{id}.{genome}.{ext}", "final_fastq"),
    wildcard_constraints:
        kind="bam|bam_low_qual|gam|gam_low_qual",
        ext="bam|gam",
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("final_fastq", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("final_fastq", attempt, 24),
    log:
        "{folder}/01_fastq/04_final/01_{kind}/{sm}/{lb}/{id}.{genome}.{ext}.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS MERGE OF SPLIT FASTQ FILES {output}"
    shell:
        """
        ## same header as without splitting:
        ##  - header of the first chunk, without the chunk names in the @PG command lines
        ##  - all chunks have the same @RG/@PG IDs: keep them once (-c/-p), instead of adding suffixes to the IDs
        ##  - no @PG line for samtools merge (--no-PG)
        if [ {wildcards.ext} = bam ]; then
            header={output}.header.sam
            samtools view -H --no-PG {input.bam[0]} | sed 's/_chunk[0-9][0-9]*//g' > $header
            samtools merge -f -c -p --no-PG -h $header --threads {threads} {output} {input.bam} 2> {log};
            rm $header
        else
            cat {input.bam} > {output} 2> {log};
        fi
        """
