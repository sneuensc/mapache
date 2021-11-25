##########################################################################################
## all rules for fastq files
##########################################################################################
## get/rename reference and fastq files


localrules:
    get_fastq,
    get_fasta,  ## executed locally on a cluster


## all rules for fastq files
rule get_fastq:
    """
    Subsample or symlink and rename all fastq files to a common folder (makes the DAG readable)
    """
    input:
        get_fastq_of_ID,
    output:
        "{folder}/00_reads/01_files_orig/{SM}/{LB}/{ID}.fastq.gz",
    threads: 1
    params:
        run=recursive_get(["subsampling","run"], False),
        number=recursive_get(["subsampling","number"], 100),
        params=recursive_get(["subsampling","params"], "-s1"),
    conda:
        "../envs/seqtk.yaml"
    envmodules:
        module_seqtk,
    message:
        "--- GET FASTQ FILES  {input}"
    log:
        "{folder}/00_reads/01_files_orig/{SM}/{LB}/{ID}.fastq.gz.log",
    script:
        "../scripts/subsample_fastq.py"


## all rules for fastq files
rule get_fasta:
    """
    Symlink and rename the reference (.fasta/.fa) to a new folder.
    """
    input:
        lambda wildcards: recursive_get(["genome",wildcards.GENOME,"fasta"], ""),
    output:
        "{folder}/00_reference/{GENOME}/{GENOME}.fasta",
    threads: 1
    message:
        "--- GET REFERENCE  {input}"
    log:
        "{folder}/00_reference/{GENOME}/{GENOME}.fasta.log",
    shell:
        """
        ln -srf {input} {output}
        """


##########################################################################################
## trimming
if paired_end:
    if collapse:
        ruleorder: adapter_removal_collapse > adapter_removal_pe > adapter_removal_se
    else:
        ruleorder:  adapter_removal_pe > adapter_removal_collapse  > adapter_removal_se
else:
    ruleorder: adapter_removal_se > adapter_removal_collapse > adapter_removal_pe


rule adapter_removal_collapse:
    """
    Remove adapter and low quality bases at the edges and collapse paired-end reads
    """
    input:
        R1="{folder}/01_fastq/00_reads/01_files_orig/{SM}/{LB}/{ID}_R1.fastq.gz",
        R2="{folder}/01_fastq/00_reads/01_files_orig/{SM}/{LB}/{ID}_R2.fastq.gz",
    output:
        R=       "{folder}/01_fastq/01_trimmed/01_adapter_removal_collapsed/{SM}/{LB}/{ID}.fastq.gz",
        trunc=   "{folder}/01_fastq/01_trimmed/01_adapter_removal_collapsed/{SM}/{LB}/{ID}_truncated.fastq.gz",
        R1=      "{folder}/01_fastq/01_trimmed/01_adapter_removal_collapsed/{SM}/{LB}/{ID}_R1.fastq.gz",
        R2=      "{folder}/01_fastq/01_trimmed/01_adapter_removal_collapsed/{SM}/{LB}/{ID}_R2.fastq.gz",
        strunc=  "{folder}/01_fastq/01_trimmed/01_adapter_removal_collapsed/{SM}/{LB}/{ID}.singleton.truncated.gz",
        disc=    "{folder}/01_fastq/01_trimmed/01_adapter_removal_collapsed/{SM}/{LB}/{ID}.discarded.gz",
        settings="{folder}/01_fastq/01_trimmed/01_adapter_removal_collapsed/{SM}/{LB}/{ID}.settings",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("adapterremoval", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "adapterremoval", attempt, 24
        ),
    params:
        recursive_get(
            ["adapterremoval","params_paired_end"],
            "--minlength 30 --trimns --trimqualities --collapse"
            
        ),
    log:
        "{folder}/01_fastq/01_trimmed/01_adapter_removal_collapsed/{SM}/{LB}/{ID}.log",
    threads: get_threads("adapterremoval", 4)
    conda:
        "../envs/adapterremoval.yaml"
    envmodules:
        module_adapterremoval,
    message:
        "--- ADAPTERREMOVAL PAIRED-END COLLAPSED {input.R1} {input.R2}"
    shell:
        """
        out={output.R};
        AdapterRemoval --threads {threads} {params} --file1 {input.R1} \
                --file2 {input.R2} --basename ${{out%%.fastq.gz}} --gzip \
                --output1 {output.R1} --output2 {output.R2} \
                --outputcollapsed {output.R} \
                --outputcollapsedtruncated {output.trunc} 2> {log};
        """


rule adapter_removal_pe:
    """
    Remove adapter and low quality bases at the edges
    """
    input:
        R1="{folder}/01_fastq/00_reads/01_files_orig/{SM}/{LB}/{ID}_R1.fastq.gz",
        R2="{folder}/01_fastq/00_reads/01_files_orig/{SM}/{LB}/{ID}_R2.fastq.gz",
    output:
        R1="{folder}/01_fastq/01_trimmed/01_adapter_removal/{SM}/{LB}/{ID}_R1.fastq.gz",
        R2="{folder}/01_fastq/01_trimmed/01_adapter_removal/{SM}/{LB}/{ID}_R2.fastq.gz",
        singleton="{folder}/01_fastq/01_trimmed/01_adapter_removal/{SM}/{LB}/{ID}.singleton.truncated.gz",
        discarded="{folder}/01_fastq/01_trimmed/01_adapter_removal/{SM}/{LB}/{ID}.discarded.gz",
        settings="{folder}/01_fastq/01_trimmed/01_adapter_removal/{SM}/{LB}/{ID}.settings",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("adapterremoval", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "adapterremoval", attempt, 24
        ),
    params:
        recursive_get(
            ["adapterremoval", "params_paired_end"],
            "--minlength 30 --trimns --trimqualities"
        ),
    log:
        "{folder}/01_fastq/01_trimmed/01_adapter_removal/{SM}/{LB}/{ID}.log",
    threads: get_threads("adapterremoval", 4)
    conda:
        "../envs/adapterremoval.yaml"
    envmodules:
        module_adapterremoval,
    message:
        "--- ADAPTERREMOVAL PAIRED-END {input.R1} {input.R2}"
    shell:
        """
        out={output.R1};
        AdapterRemoval --threads {threads} {params} --file1 {input.R1} \
                --file2 {input.R2} --basename ${{out%%_R1.fastq.gz}} --gzip \
                --output1 {output.R1} --output2 {output.R2} 2> {log};
        """



rule adapter_removal_se:
    """
    Remove adapter and low quality bases at the edges
    """
    input:
        "{folder}/01_fastq/00_reads/01_files_orig/{SM}/{LB}/{ID}.fastq.gz",
    output:
        fastq="{folder}/01_fastq/01_trimmed/01_adapter_removal/{SM}/{LB}/{ID}.fastq.gz",
        discard="{folder}/01_fastq/01_trimmed/01_adapter_removal/{SM}/{LB}/{ID}.discarded.gz",
        settings="{folder}/01_fastq/01_trimmed/01_adapter_removal/{SM}/{LB}/{ID}.settings",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("adapterremoval", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "adapterremoval", attempt, 24
        ),
    params:
        recursive_get(
            ["adapterremoval","params_single_end"],
            "--minlength 30 --trimns --trimqualities"
        ),
    log:
        "{folder}/01_fastq/01_trimmed/01_adapter_removal/{SM}/{LB}/{ID}.log",
    threads: get_threads("adapterremoval", 4)
    conda:
        "../envs/adapterremoval.yaml"
    envmodules:
        module_adapterremoval,
    message:
        "--- ADAPTERREMOVAL SINGLE-END {input}"
    shell:
        """
        out={output.fastq};
        AdapterRemoval --threads {threads} {params} --file1 {input} \
                --basename ${{out%%.fastq.gz}} --gzip \
                --output1 {output.fastq} 2> {log};
        """






##########################################################################################
## mapping
ruleorder: mapping_bwa_aln_pe > mapping_bwa_aln_se


rule mapping_bwa_aln_se:
    """
    Align reads to the reference
    """
    input:
        multiext(
            "{folder}/00_reference/{GENOME}/{GENOME}.fasta",
            ".sa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
        ),
        ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        fastq=lambda wildcards: get_fastq_for_mapping(wildcards, run_adapter_removal=run_adapter_removal),
    output:
        "{folder}/01_fastq/02_mapped/01_bwa_aln/{SM}/{LB}/{ID}.{GENOME}.sai",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        recursive_get(["mapping","bwa_aln_params"], "-l 1024"),
    log:
        "{folder}/01_fastq/02_mapped/01_bwa_aln/{SM}/{LB}/{ID}.{GENOME}.log",
    threads: get_threads("mapping", 4)
    conda:
        "../envs/bwa.yaml"
    envmodules:
        module_bwa,
    message:
        "--- BWA ALN  {input.fastq}"
    shell:
        """
        bwa aln {params} -t {threads} {input.ref} -f {output} {input.fastq} 2> {log}
        """


rule mapping_bwa_aln_pe:
    """
    Align reads to the reference
    """
    input:
        multiext(
            "{folder}/00_reference/{GENOME}/{GENOME}.fasta",
            ".sa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
        ),
        ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        fastq=lambda wildcards: get_fastq_for_mapping(wildcards, run_adapter_removal=run_adapter_removal),
    output:
        "{folder}/01_fastq/02_mapped/01_bwa_aln/{SM}/{LB}/{ID}.{GENOME}_R{id_read}.sai",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        recursive_get(["mapping","bwa_aln_params"], "-l 1024"),
    log:
        "{folder}/01_fastq/02_mapped/01_bwa_aln/{SM}/{LB}/{ID}.{GENOME}_R{id_read}.log",
    threads: get_threads("mapping", 4)
    conda:
        "../envs/bwa.yaml"
    envmodules:
        module_bwa,
    message:
        "--- BWA ALN  {input.fastq}"
    shell:
        """
        bwa aln {params} -t {threads} {input.ref} -f {output} {input.fastq} 2> {log}
        """


rule mapping_bwa_samse:
    """
    Creates bam file from sai file for SE reads
    """
    input:
        multiext(
            "{folder}/00_reference/{GENOME}/{GENOME}.fasta",
            ".sa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
        ),
        ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        fastq=lambda wildcards: get_fastq_for_mapping(wildcards, run_adapter_removal=run_adapter_removal),
        sai="{folder}/01_fastq/02_mapped/01_bwa_aln/{SM}/{LB}/{ID}.{GENOME}.sai",
    output:
        "{folder}/01_fastq/02_mapped/02_bwa_samse/{SM}/{LB}/{ID}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        PL=lambda wildcards: samples[wildcards.SM][wildcards.LB][wildcards.ID]["PL"],
        bwa_samse_params=recursive_get(["mapping","bwa_samse_params"], "-n 3"),
    log:
        "{folder}/01_fastq/02_mapped/02_bwa_samse/{SM}/{LB}/{ID}.{GENOME}.log",
    threads: 1
    conda:
        "../envs/bwa.yaml"
    envmodules:
        module_bwa,
        module_samtools,
    message:
        "--- BWA SAMSE  {input.fastq}"
    shell:
        """
        (bwa samse {params.bwa_samse_params} \
         -r \"@RG\\tID:{wildcards.ID}\\tLB:{wildcards.LB}\\tSM:{wildcards.SM}\\tPL:{params.PL}\" \
         {input.ref} {input.sai} {input.fastq} | samtools view -Sb > {output}) 2> {log}
        """


rule mapping_bwa_sampe:
    """
    Creates bam file from sai file for PE reads
    """
    input:
        multiext(
            "{folder}/00_reference/{GENOME}/{GENOME}.fasta",
            ".sa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
        ),
        ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        fastq=lambda wildcards: get_fastq_for_mapping(wildcards, run_adapter_removal=run_adapter_removal),  ## should get both pairs
        sai1="{folder}/01_fastq/02_mapped/01_bwa_aln/{SM}/{LB}/{ID}.{GENOME}_R1.sai",
        sai2="{folder}/01_fastq/02_mapped/01_bwa_aln/{SM}/{LB}/{ID}.{GENOME}_R2.sai",
    output:
        "{folder}/01_fastq/02_mapped/02_bwa_sampe/{SM}/{LB}/{ID}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        PL=lambda wildcards: samples[wildcards.SM][wildcards.LB][wildcards.ID]["PL"],
        bwa_samse_params=recursive_get(["mapping","bwa_samse_params"], "-n 3"),
    log:
        "{folder}/01_fastq/02_mapped/02_bwa_sampe/{SM}/{LB}/{ID}.{GENOME}.log",
    threads: 1
    conda:
        "../envs/bwa.yaml"
    envmodules:
        module_bwa,
        module_samtools,
    message:
        "--- BWA SAMPE {input.fastq}"
    shell:
        """
        (bwa sampe {params.bwa_samse_params} \
             -r \"@RG\\tID:{wildcards.ID}\\tLB:{wildcards.LB}\\tSM:{wildcards.SM}\\tPL:{params.PL}\" \
             {input.ref} {input.sai1} {input.sai2} {input.fastq} | \
             samtools view -Sb > {output}) 2> {log}
        """


rule mapping_bwa_mem:
    """
    Map reads to genome using bwa mem
    """
    input:
        multiext(
            "{folder}/00_reference/{GENOME}/{GENOME}.fasta",
            ".sa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
        ),
        ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        fastq=lambda wildcards: get_fastq_for_mapping(wildcards, run_adapter_removal=run_adapter_removal),
    output:
        "{folder}/01_fastq/02_mapped/02_bwa_mem/{SM}/{LB}/{ID}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        PL=lambda wildcards: samples[wildcards.SM][wildcards.LB][wildcards.ID]["PL"],
        bwa_mem_params=recursive_get(["mapping","bwa_mem_params"], ""),
    log:
        "{folder}/01_fastq/02_mapped/02_bwa_mem/{SM}/{LB}/{ID}.{GENOME}.log",
    threads: get_threads("mapping", 4)
    conda:
        "../envs/bwa.yaml"
    envmodules:
        module_bwa,
    message:
        "--- BWA MEM {input.fastq}"
    shell:
        """
        bwa mem {params.bwa_mem_params} -t {threads} \
            -R \"@RG\\tID:{wildcards.ID}\\tLB:{wildcards.LB}\\tSM:{wildcards.SM}\\tPL:{params.PL}\" \
            {input.ref} ${input.fastq} > {output} 2> {log};
        """


rule mapping_bowtie2:
    """
    Map reads to genome using bowtie2
    """
    input:
        multiext(
            "{folder}/00_reference/{GENOME}/{GENOME}.fasta",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        ref="{folder}/00_reference/{GENOME}/{GENOME}.fasta",
        fastq=lambda wildcards: get_fastq_for_mapping(wildcards, run_adapter_removal=run_adapter_removal),
    output:
        "{folder}/01_fastq/02_mapped/02_bwa_bowtie2/{SM}/{LB}/{ID}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        bowtie2_params=recursive_get(["mapping","bowtie2_params"], ""),
    log:
        "{folder}/01_fastq/02_mapped/02_bwa_bowtie2/{SM}/{LB}/{ID}.{GENOME}.log",
    threads: get_threads("mapping", 4)
    conda:
        "../envs/bowtie2.yaml"
    envmodules:
        module_bowtie2,
        module_samtools,
    message:
        "--- BOWTIE2 {input.fastq}"
    script:
        "../scripts/mapping_bowtie2.py"


##########################################################################################
## sorting


rule samtools_sort:
    """
    Sort bam file with samtools
    """
    input:
        get_bam_for_sorting,
    output:
        "{folder}/01_fastq/02_mapped/03_bam_sort/{SM}/{LB}/{ID}.{GENOME}.bam",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("sorting", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("sorting", attempt, 24),
    log:
        "{folder}/01_fastq/02_mapped/03_bam_sort/{SM}/{LB}/{ID}.{GENOME}.log",
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
            "{folder}/01_fastq/02_mapped/03_bam_sort/{SM}/{LB}/{ID}.{GENOME}.bam",
        output:
            mapped="{folder}/01_fastq/03_filtered/01_bam_filter/{SM}/{LB}/{ID}.{GENOME}.bam",
            low_qual="{folder}/01_fastq/03_filtered/01_bam_filter_low_qual/{SM}/{LB}/{ID}.{GENOME}.bam",
        params:
            q=lambda wildcards: samples[wildcards.SM][wildcards.LB][wildcards.ID][
                "MAPQ"
            ],
        resources:
            memory=lambda wildcards, attempt: get_memory_alloc("filtering", attempt, 4),
            runtime=lambda wildcards, attempt: get_runtime_alloc(
                "filtering", attempt, 24
            ),
        log:
            "{folder}/01_fastq/03_filtered/01_bam_filter/{SM}/{LB}/{ID}.{GENOME}.log",
        threads: get_threads("filtering", 4)
        conda:
            "../envs/samtools.yaml"
        envmodules:
            module_samtools,
        message:
            "--- SAMTOOLS FILTER {input}"
        shell:
            """
            samtools view -b --threads {threads} -F 4 -q {params.q} \
            -U {output.low_qual} {input} > {output.mapped} 2> {log}
            """


else:

    rule samtools_filter:
        """
        Filter mappings following quality and discard low quality mappings
        """
        input:
            "{folder}/01_fastq/02_mapped/03_bam_sort/{SM}/{LB}/{ID}.{GENOME}.bam",
        output:
            mapped="{folder}/01_fastq/03_filtered/01_bam_filter/{SM}/{LB}/{ID}.{GENOME}.bam",
        params:
            q=lambda wildcards: samples[wildcards.SM][wildcards.LB][wildcards.ID][
                "MAPQ"
            ],
        resources:
            memory=lambda wildcards, attempt: get_memory_alloc("filtering", attempt, 4),
            runtime=lambda wildcards, attempt: get_runtime_alloc(
                "filtering", attempt, 24
            ),
        log:
            "{folder}/01_fastq/03_filtered/01_bam_filter/{SM}/{LB}/{ID}.{GENOME}.log",
        threads: get_threads("filtering", 4)
        conda:
            "../envs/samtools.yaml"
        envmodules:
            module_samtools,
        message:
            "--- SAMTOOLS FILTER {input}"
        shell:
            """
            samtools view -b --threads {threads} -F 4 -q {params.q} {input} > {output.mapped} 2> {log}
            """


##########################################################################################
rule get_final_fastq:
    """
    Get the final bam file from the fastq part
    """
    input:
        get_final_bam_fastq,
    output:
        "{folder}/01_fastq/04_final_fastq/{type}/{SM}/{LB}/{ID}.{GENOME}.bam",
    message:
        "--- GET FINAL BAM {input} (FASTQ LEVEL)"
    log:
        "{folder}/01_fastq/04_final_fastq/{type}/{SM}/{LB}/{ID}.{GENOME}.bam.log",
    run:
        symlink_rev(input, output)


rule get_final_fastq_low_qual:
    """
    Get the final bam file from the fastq part
    """
    input:
        "{folder}/03_filtered/01_bam_filter_low_qual/{SM}/{LB}/{ID}.{GENOME}.bam",
    output:
        "{folder}/04_final_fastq/01_bam_low_qual/{SM}/{LB}/{ID}.{GENOME}.bam",
    message:
        "--- GET FINAL LOW_QUAL BAM {input} (FASTQ LEVEL)"
    log:
        "{folder}/04_final_fastq/01_bam_low_qual/{SM}/{LB}/{ID}.{GENOME}.bam.log",
    run:
        symlink_rev(input, output)
