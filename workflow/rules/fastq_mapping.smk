## mapping
ruleorder: mapping_bwa_aln_pe > mapping_bwa_aln_se


rule mapping_bwa_aln_se:
    """
    Align reads to the reference
    """
    input:
        multiext(
            "{folder}/00_reference/{genome}/{genome}.fasta",
            ".sa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
        ),
        ref="{folder}/00_reference/{genome}/{genome}.fasta",
        fastq=get_fastq_4_mapping,
    output:
        temp("{folder}/01_fastq/02_mapped/01_bwa_aln/{sm}/{lb}/{id}.{genome}.sai"),
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        lambda wildcards: get_paramGrp(
            ["mapping", "bwa_aln_params"], "-l 1024", wildcards
        ),
    log:
        "{folder}/01_fastq/02_mapped/01_bwa_aln/{sm}/{lb}/{id}.{genome}.log",
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
            "{folder}/00_reference/{genome}/{genome}.fasta",
            ".sa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
        ),
        ref="{folder}/00_reference/{genome}/{genome}.fasta",
        fastq=get_fastq_4_mapping,
    output:
        temp(
            "{folder}/01_fastq/02_mapped/01_bwa_aln/{sm}/{lb}/{id}.{genome}_R{id_read}.sai"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        lambda wildcards: get_paramGrp(
            ["mapping", "bwa_aln_params"], "-l 1024", wildcards
        ),
    log:
        "{folder}/01_fastq/02_mapped/01_bwa_aln/{sm}/{lb}/{id}.{genome}_R{id_read}.log",
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
            "{folder}/00_reference/{genome}/{genome}.fasta",
            ".sa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
        ),
        ref="{folder}/00_reference/{genome}/{genome}.fasta",
        fastq=get_fastq_4_mapping,
        sai="{folder}/01_fastq/02_mapped/01_bwa_aln/{sm}/{lb}/{id}.{genome}.sai",
    output:
        temp("{folder}/01_fastq/02_mapped/02_bwa_samse/{sm}/{lb}/{id}.{genome}.bam"),
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        PL=lambda wildcards: get_paramGrp(
            ["mapping", "platform"], "ILLUMINA", wildcards
        ),
        params=lambda wildcards: get_paramGrp(
            ["mapping", "bwa_samse_params"], "-n 3", wildcards
        ),
    log:
        "{folder}/01_fastq/02_mapped/02_bwa_samse/{sm}/{lb}/{id}.{genome}.log",
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
        (bwa samse {params.params} \
         -r \"@RG\\tID:{wildcards.id}\\tLB:{wildcards.lb}\\tSM:{wildcards.sm}\\tPL:{params.PL}\" \
         {input.ref} {input.sai} {input.fastq} | samtools view -Sb > {output}) 2> {log}
        """


rule mapping_bwa_sampe:
    """
    Creates bam file from sai file for PE reads
    """
    input:
        multiext(
            "{folder}/00_reference/{genome}/{genome}.fasta",
            ".sa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
        ),
        ref="{folder}/00_reference/{genome}/{genome}.fasta",
        fastq=get_fastq_4_mapping,
        ## should get both pairs
        sai1="{folder}/01_fastq/02_mapped/01_bwa_aln/{sm}/{lb}/{id}.{genome}_R1.sai",
        sai2="{folder}/01_fastq/02_mapped/01_bwa_aln/{sm}/{lb}/{id}.{genome}_R2.sai",
    output:
        temp("{folder}/01_fastq/02_mapped/02_bwa_sampe/{sm}/{lb}/{id}.{genome}.bam"),
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        PL=lambda wildcards: get_paramGrp(
            ["mapping", "platform"], "ILLUMINA", wildcards
        ),
        params=lambda wildcards: get_paramGrp(
            ["mapping", "bwa_sampe_params"], "-n 3", wildcards
        ),
    log:
        "{folder}/01_fastq/02_mapped/02_bwa_sampe/{sm}/{lb}/{id}.{genome}.log",
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
        (bwa sampe {params.params} \
             -r \"@RG\\tID:{wildcards.id}\\tLB:{wildcards.lb}\\tSM:{wildcards.sm}\\tPL:{params.PL}\" \
             {input.ref} {input.sai1} {input.sai2} {input.fastq} | \
             samtools view -Sb > {output}) 2> {log}
        """


rule mapping_bwa_mem:
    """
    Map reads to GENOMES using bwa mem
    """
    input:
        multiext(
            "{folder}/00_reference/{genome}/{genome}.fasta",
            ".sa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
        ),
        ref="{folder}/00_reference/{genome}/{genome}.fasta",
        fastq=get_fastq_4_mapping,
    output:
        temp("{folder}/01_fastq/02_mapped/02_bwa_mem/{sm}/{lb}/{id}.{genome}.bam"),
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        PL=lambda wildcards: get_paramGrp(
            ["mapping", "platform"], "ILLUMINA", wildcards
        ),
        params=lambda wildcards: get_paramGrp(
            ["mapping", "bwa_mem_params"], "", wildcards
        ),
    log:
        "{folder}/01_fastq/02_mapped/02_bwa_mem/{sm}/{lb}/{id}.{genome}.log",
    threads: get_threads("mapping", 4)
    conda:
        "../envs/bwa.yaml"
    envmodules:
        module_bwa,
        module_samtools,
    message:
        "--- BWA MEM {input.fastq}"
    shell:
        """
        bwa mem {params.params} -t {threads} \
            -R \"@RG\\tID:{wildcards.id}\\tLB:{wildcards.lb}\\tSM:{wildcards.sm}\\tPL:{params.PL}\" \
            {input.ref} {input.fastq} 2> {log} | samtools view -bS --threads {threads} - > {output};
        """


rule mapping_bowtie2:
    """
    Map reads to GENOMES using bowtie2
    """
    input:
        multiext(
            "{folder}/00_reference/{genome}/{genome}.fasta",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        ref="{folder}/00_reference/{genome}/{genome}.fasta",
        fastq=get_fastq_4_mapping,
    output:
        temp("{folder}/01_fastq/02_mapped/02_bwa_bowtie2/{sm}/{lb}/{id}.{genome}.bam"),
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        bowtie2_params=lambda wildcards: get_paramGrp(
            ["mapping", "bowtie2_params"], "", wildcards
        ),
    log:
        "{folder}/01_fastq/02_mapped/02_bwa_bowtie2/{sm}/{lb}/{id}.{genome}.log",
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


if config["mapping"]["mapper"] == "vg_giraffe":
    ruleorder: mapping_vg_giraffe > convert_gam2bam
else:
    ruleorder: convert_gam2bam > mapping_vg_giraffe


rule mapping_vg_giraffe:
    """
    Map reads to variation graph using vg giraffe
    """
    input:
        gbz="{folder}/00_reference/{genome}/{genome}.fasta.graph.gbz",
        dist="{folder}/00_reference/{genome}/{genome}.fasta.graph.dist",
        min="{folder}/00_reference/{genome}/{genome}.fasta.graph.min",
        zipcode="{folder}/00_reference/{genome}/{genome}.fasta.graph.zipcodes",
        #xg="{folder}/00_reference/{genome}/{genome}.fasta.graph.xg",
        dict="{folder}/00_reference/{genome}/{genome}.dict",
        fastq=get_fastq_4_mapping,
    output:
        temp("{folder}/01_fastq/02_mapped/02_vg_giraffe/{sm}/{lb}/{id}.{genome}.bam"),
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        PL=lambda wildcards: get_paramGrp(
            ["mapping", "platform"], "ILLUMINA", wildcards
        ),
        params=lambda wildcards: get_paramGrp(
            ["mapping", "giraffe_params"], "", wildcards
        ),
        bin = get_param(["software", "vg"], "vg")
    log:
        "{folder}/01_fastq/02_mapped/02_vg_giraffe/{sm}/{lb}/{id}.{genome}.log",
    threads: get_threads("mapping", 4)
    conda:
        "../envs/samtools.yaml"
    message:
        "--- VG GIRAFFE {input.fastq}"
    shell:
        """
        set -o pipefail
        (
            {params.bin} giraffe \
                -Z {input.gbz} \
                -d {input.dist} \
                -m {input.min} \
                -z {input.zipcode} \
                {params.params} \
                -t {threads} \
                -f {input.fastq} \
                -o BAM \
            | samtools addreplacerg \
                -O BAM \
                -r "@RG\\tID:{wildcards.id}\\tLB:{wildcards.lb}\\tSM:{wildcards.sm}\\tPL:{params.PL}" \
                -o - \
                - \
        ) 2> {log} > {output}.tmp;

        ## correct the order of the chromosomes (a bug in vg giraffe)
        samtools view -H {output}.tmp | grep '^@HD' > {output}.header;
        grep '^@SQ' {input.dict} | cut -f-3 >> {output}.header;
        samtools view -H {output}.tmp | grep  -Ev '^@(HD|SQ)' >> {output}.header;

        samtools reheader {output}.header {output}.tmp > {output};

        rm -f {output}.tmp {output}.header;
        """

rule mapping_vg_giraffe_gam:
    """
    Map reads to variation graph using vg giraffe (via gam file)
    """
    input:
        gbz="{folder}/00_reference/{genome}/{genome}.fasta.graph.gbz",
        dist="{folder}/00_reference/{genome}/{genome}.fasta.graph.dist",
        min="{folder}/00_reference/{genome}/{genome}.fasta.graph.min",
        zipcode="{folder}/00_reference/{genome}/{genome}.fasta.graph.zipcodes",
        xg="{folder}/00_reference/{genome}/{genome}.fasta.graph.xg",
        dict="{folder}/00_reference/{genome}/{genome}.dict",
        fastq=get_fastq_4_mapping,
    output:
        gam=temp("{folder}/01_fastq/02_mapped/02_vg_giraffe/{sm}/{lb}/{id}.{genome}.gam"),
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        PL=lambda wildcards: get_paramGrp(
            ["mapping", "platform"], "ILLUMINA", wildcards
        ),
        params=lambda wildcards: get_paramGrp(
            ["mapping", "giraffe_params"], "", wildcards
        ),
        bin = get_param(["software", "vg"], "vg")
    log:
        "{folder}/01_fastq/02_mapped/02_vg_giraffe/{sm}/{lb}/{id}.{genome}.log",
    threads: get_threads("mapping", 4)
    message:
        "--- VG GIRAFFE {input.fastq}"
    shell:
        """
        {params.bin} giraffe \
            -Z {input.gbz} \
            -d {input.dist} \
            -m {input.min} \
            -z {input.zipcode} \
            {params.params} \
            -t {threads} \
            -f {input.fastq} > {output.gam};
        """


rule convert_gam2bam:
    """
    Convert gam to bam
    """
    input:
        xg="{folder}/00_reference/{genome}/{genome}.fasta.graph.xg",
        dict="{folder}/00_reference/{genome}/{genome}.dict",
        gam="{folder}/01_fastq/02_mapped/02_vg_giraffe/{sm}/{lb}/{id}.{genome}.gam"
    output:
        bam=temp("{folder}/01_fastq/02_mapped/02_vg_giraffe/{sm}/{lb}/{id}.{genome}.bam")
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        PL=lambda wildcards: get_paramGrp(
            ["mapping", "platform"], "ILLUMINA", wildcards
        ),
        params=lambda wildcards: get_paramGrp(
            ["mapping", "giraffe_params"], "", wildcards
        ),
        bin = get_param(["software", "vg"], "vg")
    log:
        "{folder}/01_fastq/02_mapped/02_vg_giraffe/{sm}/{lb}/{id}.{genome}_convert.log",
    threads: get_threads("mapping", 4)
    conda:
        "../envs/samtools.yaml"
    message:
        "--- CONVERTING GAM TO BAM {input.gam}"
    shell:
        """
        set -euo pipefail

        exec > {log} 2>&1

        {params.bin} surject -b \
            -x {input.xg} \
            {input.gam} \
        | samtools addreplacerg \
            -O BAM \
            -r "@RG\\tID:{wildcards.id}\\tLB:{wildcards.lb}\\tSM:{wildcards.sm}\\tPL:{params.PL}" \
            -o {output.bam}.tmp \
            -;

        ## correct the order of the chromosomes (a bug in vg giraffe)
        samtools view -H {output.bam}.tmp | grep '^@HD' > {output.bam}.header;
        grep '^@SQ' {input.dict} | cut -f-3 >> {output.bam}.header;
        samtools view -H {output.bam}.tmp | grep  -Ev '^@(HD|SQ)' >> {output.bam}.header;

        samtools reheader {output.bam}.header {output.bam}.tmp > {output.bam};
        rm -f {output.bam}.tmp {output.bam}.header;
        """