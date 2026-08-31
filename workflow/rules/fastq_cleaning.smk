##########################################################################################
## all rules for fastq files
##########################################################################################
## get/rename reference and fastq files




             

## all rules for fastq files
rule get_fastq_remote:
    """
    Download a remote fastq file from an anonymous ftp server and check the md5sum
    """
    output:
        temp("{folder}/01_fastq/00_reads/00_files_orig/{sm}/{lb}/{idd}.fastq.gz"),
    wildcard_constraints:
        idd="({})".format("|".join(remote_idd))
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("download", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("download", attempt, 10),
    params:
        ftp=get_fastq_name_original,
        md5=get_md5_of_ID,
    message:
        "--- GET REMOTE FASTQ FILE {input}"
    log:
        "{folder}/01_fastq/00_reads/00_files_orig/{sm}/{lb}/{idd}.fastq.gz.log",
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


rule get_fastq_local:
    """
    Symlink and rename local fastq file
    """
    input:
        get_fastq_name_original,
    output:
        temp("{folder}/01_fastq/00_reads/00_files_orig/{sm}/{lb}/{idd}.fastq.gz"),
    wildcard_constraints:
        idd="({})".format("|".join(local_idd))
    localrule: True
    threads: 1
    message:
        "--- GET LOCAL FASTQ FILE  {input}"
    log:
        "{folder}/01_fastq/00_reads/00_files_orig/{sm}/{lb}/{idd}.fastq.gz.log",
    shell:
        """
        ln -srf {input} {output}
        """


## all rules for fastq files
rule get_fasta:
    """
    Symlink and rename the reference (.fasta/.fa) to a new folder.
    """
    input:
        lambda wildcards: get_param(["genome", wildcards.genome], ""),
    output:
        temp("{folder}/00_reference/{genome}/{genome}.fasta"),
    localrule: True
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
## subsampling

rule fastq_subsample:
    """
    Subsample fastq
    """
    input:
        "{folder}/01_fastq/00_reads/00_files_orig/{sm}/{lb}/{idd}.fastq.gz"
    output:
        temp("{folder}/01_fastq/00_reads/01_subsample/{sm}/{lb}/{idd}.fastq.gz"),
    threads: 1
    params:
        run=lambda wildcards: get_paramGrp(["subsampling", "run"], False, wildcards),
        number=lambda wildcards: get_paramGrp(
            ["subsampling", "number"], 1000, wildcards
        ),
        params=lambda wildcards: get_paramGrp(
            ["subsampling", "params"], "-s1", wildcards
        ),
    conda:
        "../envs/seqtk.yaml"
    envmodules:
        module_seqtk,
    message:
        "--- SUBSAMPLE FASTQ FILE  {input}"
    log:
        "{folder}/01_fastq/00_reads/01_subsample/{sm}/{lb}/{idd}.fastq.gz.log",
    shell:
        """
        seqtk sample {params.params} {input} {params.number} | gzip > {output}
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
        mem_mb=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
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
        "{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}.log",
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
            mv {output.col_R} {output.R}; ## symlink is not working with temp(): 'invert' it !!!
            ln -srf {output.R} {output.col_R};
        elif [[ "$options" == "collapse_trunc" ]]; then
            cat {output.col_R} {output.col_trunc} > {output.R};
        elif [[ "$options" == "all" ]]; then
            cat {output.col_R} {output.col_trunc} {output.R1} {output.R2} {output.strunc} > {output.R};
        fi;
        """


rule adapterremoval_pe:
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
        mem_mb=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
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
        mem_mb=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("cleaning", attempt, 24),
    params:
        lambda wildcards: get_paramGrp(
            ["cleaning", "params_adapterremoval"],
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
        mem_mb=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
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
        mem_mb=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
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
        mem_mb=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
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
        mem_mb=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
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
