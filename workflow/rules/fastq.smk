##########################################################################################
## all rules for fastq files
##########################################################################################
## get/rename reference and fastq files


localrules:
    get_fastq,
    get_fasta,  ## executed locally on a cluster


## all rules for fastq files
rule get_fastq_remote:
    """
    Download a remote fastq file from an anonymous ftp server and check the md5sum
    """
    output:
        "{folder}/00_reads/00_files_remote/{sm}/{lb}/{idd}.fastq.gz",
    threads: 1
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("download", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("download", attempt, 10),
    params:
        ftp=get_fastq_of_ID_0,
        md5=get_md5_of_ID,
    message:
        "--- GET FASTQ FILES REMOTELY {input}"
    log:
        "{folder}/00_reads/00_files_remote/{sm}/{lb}/{idd}.fastq.gz.log",
    shell:
        """
        ## download file
        wget -O {output} {params.ftp} > {log};

        ## test md5sum if available
        if [ "{params.md5}" == "''" ] || [ "{params.md5}" == "nan" ] ; then
            echo "WARNING: Downloaded fastq file '{params.ftp}' has no md5sum to verify the download!";
        else
            if [ $(md5sum {output} | cut -d' ' -f1) != "{params.md5}" ]; then
                echo "ERROR: Downloaded fastq file '{params.ftp}' has a wrong md5sum!";
                exit 1;
            fi
        fi
        """


## all rules for fastq files
rule get_fastq:
    """
    Subsample or symlink and rename all fastq files to a common folder (makes the DAG readable)
    """
    input:
        get_fastq_of_ID,
    output:
        temp("{folder}/00_reads/01_files_orig/{sm}/{lb}/{idd}.fastq.gz"),
    threads: 1
    params:
        run=lambda wildcards: get_paramGrp(["subsampling", "run"], False, wildcards),
        number=lambda wildcards: get_paramGrp(
            ["subsampling", "number"], 1000, wildcards
        ),
        params=lambda wildcards: get_paramGrp(
            ["subsampling", "params"], "-s1", wildcards
        ),
        ftp=lambda wildcards: get_fastq_of_ID_0(wildcards),
    conda:
        "../envs/seqtk.yaml"
    envmodules:
        module_seqtk,
    message:
        "--- GET FASTQ FILES  {input}"
    log:
        "{folder}/00_reads/01_files_orig/{sm}/{lb}/{idd}.fastq.gz.log",
    script:
        "../scripts/subsample_fastq.py"


## all rules for fastq files
rule get_fasta:
    """
    Symlink and rename the reference (.fasta/.fa) to a new folder.
    """
    input:
        lambda wildcards: get_param(["genome", wildcards.genome], ""),
    output:
        temp("{folder}/00_reference/{genome}/{genome}.fasta"),
    threads: 1
    message:
        "--- GET REFERENCE  {input}"
    log:
        "{folder}/00_reference/{genome}/{genome}.fasta.log",
    shell:
        """
        ln -srf {input} {output}
        """


##########################################################################################
## trimming


rule adapterremoval_collapse:
    """
    Remove adapter and low quality bases at the edges and collapse paired-end reads
    """
    input:
        get_fastq_4_cleaning,
    output:
        R=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}.fastq.gz"
        ),
        col_R=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}_collapsed.fastq.gz"
        ),
        col_trunc=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}_collapsed_truncated.fastq.gz"
        ),
        R1=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}_R1.fastq.gz"
        ),
        R2=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}_R2.fastq.gz"
        ),
        strunc=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}.singleton.truncated.gz"
        ),
        disc=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}.discarded.gz"
        ),
        settings="{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}.settings",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("cleaning", attempt, 24),
    params:
        params=lambda wildcards: get_paramGrp(
            ["cleaning", "params_adapterremoval"],
            "--minlength 30 --trimns --trimqualities",
            wildcards,
        ),
        collapsed=lambda wildcards: get_paramGrp(
            ["cleaning", "collapse_opt"],
            ["only_collapse", "collapse_trunc", "all"],
            wildcards,
        ),
    log:
        "{folder}/01_fastq/01_trimmed/01_adapter_removal_collapse/{sm}/{lb}/{id}.log",
    threads: get_threads("cleaning", 4)
    conda:
        "../envs/adapterremoval.yaml"
    envmodules:
        module_adapterremoval,
    message:
        "--- ADAPTERREMOVAL PAIRED-END COLLAPSED {input[0]} {input[1]}"
    shell:
        """
        out={output.R};
        AdapterRemoval --threads {threads} {params.params} --file1 {input[0]} \
                --file2 {input[1]} --basename ${{out%%.fastq.gz}} --gzip \
                --output1 {output.R1} --output2 {output.R2} \
                --outputcollapsed {output.col_R} \
                --outputcollapsedtruncated {output.col_trunc} 2> {log};

        ## what should be retained
        options={params.collapsed};
        if [[ "$options" == "only_collapse" ]]; then
            ln -srf {output.col_R} {output.R};
        elif [[ "$options" == "collapse_trunc" ]]; then
            cat {output.col_R} {output.col_trunc} > {output.R};
        elif [[ "$options" == "all" ]]; then
            cat {output.col_R} {output.col_trunc} {output.R1} {output.R2} {output.strunc} > {output.R};
        fi;
        """


rule adapter_removal_pe:
    """
    Remove adapter and low quality bases at the edges
    """
    input:
        get_fastq_4_cleaning,
    output:
        R1=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_pe/{sm}/{lb}/{id}_R1.fastq.gz"
        ),
        R2=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_pe/{sm}/{lb}/{id}_R2.fastq.gz"
        ),
        singleton=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_pe/{sm}/{lb}/{id}.singleton.truncated.gz"
        ),
        discarded=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_pe/{sm}/{lb}/{id}.discarded.gz"
        ),
        settings="{folder}/01_fastq/01_trimmed/01_adapterremoval_pe/{sm}/{lb}/{id}.settings",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("cleaning", attempt, 24),
    params:
        lambda wildcards: get_paramGrp(
            ["cleaning", "params_adapterremoval"],
            "--minlength 30 --trimns --trimqualities",
            wildcards,
        ),
    log:
        "{folder}/01_fastq/01_trimmed/01_adapterremoval_pe/{sm}/{lb}/{id}.log",
    threads: get_threads("cleaning", 4)
    conda:
        "../envs/adapterremoval.yaml"
    envmodules:
        module_adapterremoval,
    message:
        "--- ADAPTERREMOVAL PAIRED-END {input[0]} {input[1]}"
    shell:
        """
        out={output.R1};
        AdapterRemoval --threads {threads} {params} --file1 {input[0]} \
                --file2 {input[1]} --basename ${{out%%_R1.fastq.gz}} --gzip \
                --output1 {output.R1} --output2 {output.R2} 2> {log};
        """


rule adapterremoval_se:
    """
    Remove adapter and low quality bases at the edges
    """
    input:
        get_fastq_4_cleaning,
    output:
        R=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_se/{sm}/{lb}/{id}.fastq.gz"
        ),
        discard=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_se/{sm}/{lb}/{id}.discarded.gz"
        ),
        settings="{folder}/01_fastq/01_trimmed/01_adapterremoval_se/{sm}/{lb}/{id}.settings",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("cleaning", attempt, 24),
    params:
        lambda wildcards: get_paramGrp(
            ["cleaning", "params_dapterremoval"],
            "--minlength 30 --trimns --trimqualities",
            wildcards,
        ),
    log:
        "{folder}/01_fastq/01_trimmed/01_adapterremoval_se/{sm}/{lb}/{id}.log",
    threads: get_threads("cleaning", 4)
    conda:
        "../envs/adapterremoval.yaml"
    envmodules:
        module_adapterremoval,
    message:
        "--- ADAPTERREMOVAL SINGLE-END {input}"
    shell:
        """
        ## remove --collapse from $params (needed for SE libs in a paired-end setting)
        params=$(echo {params} | sed  's/ --collapse[^ ]*//g')

        out={output.R};
        AdapterRemoval --threads {threads} $params --file1 {input} \
                --basename ${{out%%.fastq.gz}} --gzip \
                --output1 {output.R} 2> {log};
        """


rule adapterremoval_infer_adapters:
    """
    Remove adapter and low quality bases at the edges and collapse paired-end reads
    """
    input:
        get_fastq_4_cleaning,
    output:
        adapters="{folder}/01_fastq/01_trimmed/01_adapterremoval_infer_adapters/{sm}/{lb}/{id}.txt",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("cleaning", attempt, 24),
    params:
        params=lambda wildcards: get_paramGrp(
            ["cleaning", "params_adapterremoval"],
            "--minlength 30 --trimns --trimqualities",
            wildcards,
        ),
        collapsed=lambda wildcards: get_paramGrp(
            ["cleaning", "collapse_opt"],
            ["only_collapse", "collapse_trunc", "all"],
            wildcards,
        ),
    log:
        "{folder}/01_fastq/01_trimmed/01_adapterremoval_infer_adapters/{sm}/{lb}/{id}.log",
    threads: get_threads("cleaning", 1)
    conda:
        "../envs/adapterremoval.yaml"
    envmodules:
        module_adapterremoval,
    message:
        "--- ADAPTERREMOVAL INFER ADAPTERTS {input}"
    shell:
        """
        set +e;
        AdapterRemoval --threads {threads} {params.params} --file1 {input[0]} \
                --file2 {input[1]} --identify-adapters > {output};
        """


rule fastp_collapse:
    """
    Clean fastq files with fastp (collapse paired-end)
    """
    input:
        get_fastq_4_cleaning,
    output:
        R=temp("{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}.fastq.gz"),
        R_merged=temp(
            "{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}_merged.fastq.gz"
        ),
        R1=temp("{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}_R1.fastq.gz"),
        R2=temp("{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}_R2.fastq.gz"),
        json="{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}.json",
        html="{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}.html",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("cleaning", attempt, 24),
    params:
        lambda wildcards: get_paramGrp(
            ["cleaning", "params_fastp"],
            "",
            wildcards,
        ),
    log:
        "{folder}/01_fastq/01_trimmed/02_fastp_pe/{sm}/{lb}/{id}.log",
    threads: get_threads("cleaning", 4)
    conda:
        "../envs/fastp.yaml"
    envmodules:
        module_fastp,
    message:
        "--- FASTP PAIRED-END {input}"
    shell:
        """
        fastp --in1 {input[0]} --in2 {input[1]} --out1 {output.R1} --out2 {output.R2} \
              --merged_out {output.R_merged} \
              --json {output.json} --html {output.html} --thread {threads} {params} \
              -R "Fastp report of {wildcards.sm}/{wildcards.lb}/{wildcards.id}" 2> {log};

        ## what should be retained
        options={params.collapsed};
        if [[ "$options" == "only_collapse" ]]; then
            ln -srf {output.R_merged} {output.R};
        elif [[ "$options" == "all" ]]; then
            cat {output.R_merged} {output.R1} {output.R2} > {output.R};
        fi;
        """


rule fastp_pe:
    """
    Clean fastq files with fastp (paired-end)
    """
    input:
        get_fastq_4_cleaning,
    output:
        R1=temp("{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}_R1.fastq.gz"),
        R2=temp("{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}_R2.fastq.gz"),
        json="{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}.json",
        html="{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}.html",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("cleaning", attempt, 24),
    params:
        lambda wildcards: get_paramGrp(
            ["cleaning", "params_fastp"],
            "",
            wildcards,
        ),
    log:
        "{folder}/01_fastq/01_trimmed/02_fastp_pe/{sm}/{lb}/{id}.log",
    threads: get_threads("cleaning", 4)
    conda:
        "../envs/fastp.yaml"
    envmodules:
        module_fastp,
    message:
        "--- FASTP PAIRED-END {input}"
    shell:
        """
        fastp --in1 {input[0]} --in2 {input[1]} --out1 {output.R1} --out2 {output.R2} \
              --json {output.json} --html {output.html} --thread {threads} {params} \
              -R "Fastp report of {wildcards.sm}/{wildcards.lb}/{wildcards.id}" 2> {log};
        """


rule fastp_se:
    """
    Clean fastq files with fastp (single-end)
    """
    input:
        get_fastq_4_cleaning,
    output:
        R=temp("{folder}/01_fastq/01_trimmed/01_fastp_se/{sm}/{lb}/{id}.fastq.gz"),
        json="{folder}/01_fastq/01_trimmed/01_fastp_se/{sm}/{lb}/{id}.json",
        html="{folder}/01_fastq/01_trimmed/01_fastp_se/{sm}/{lb}/{id}.html",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("cleaning", attempt, 24),
    params:
        lambda wildcards: get_paramGrp(
            ["cleaning", "params_fastp"],
            "",
            wildcards,
        ),
    log:
        "{folder}/01_fastq/01_trimmed/01_fastp_se/{sm}/{lb}/{id}.log",
    threads: get_threads("cleaning", 4)
    conda:
        "../envs/fastp.yaml"
    envmodules:
        module_fastp,
    message:
        "--- FASTP SINGLE-END {input}"
    shell:
        """
        fastp --in1 {input} --out1 {output.R} \
               --json {output.json} --html {output.html} --thread {threads} {params} \
              -R "Fastp report of {wildcards.sm}/{wildcards.lb}/{wildcards.id}" 2> {log};
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
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
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
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
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
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
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
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
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
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
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
        temp("{folder}/01_fastq/02_mapped/02_bowtie2/{sm}/{lb}/{id}.{genome}.bam"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24),
    params:
        bowtie2_params=lambda wildcards: get_paramGrp(
            ["mapping", "bowtie2_params"], "", wildcards
        ),
    log:
        "{folder}/01_fastq/02_mapped/02_bowtie2/{sm}/{lb}/{id}.{genome}.log",
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
        get_bam_4_sorting,
    output:
        temp("{folder}/01_fastq/02_mapped/03_bam_sort/{sm}/{lb}/{id}.{genome}.bam"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("sorting", attempt, 4),
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
            memory=lambda wildcards, attempt: get_memory_alloc("filtering", attempt, 4),
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
            memory=lambda wildcards, attempt: get_memory_alloc("filtering", attempt, 4),
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
            samtools view -b --threads {threads} {params} {input} > {output.mapped} 2> {log}
            """
