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
        maybe_temp("{folder}/01_fastq/00_reads/00_files_orig/{sm}/{lb}/{idd}.fastq.gz", "download"),
    wildcard_constraints:
        idd="({})".format("|".join(remote_idd))
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("download", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("download", attempt, 10),
    params:
        ftp=get_fastq_name_original,
        md5=get_md5_of_ID,
        tries=lambda wildcards: get_paramGrp(
            ["download", "tries"], "3", wildcards
        ),
    message:
        "--- GET REMOTE FASTQ FILE {input}"
    log:
        "{folder}/01_fastq/00_reads/00_files_orig/{sm}/{lb}/{idd}.fastq.gz.log",
    shell:
        """
        ## download file
        wget -c --tries={params.tries} -O {output} {params.ftp} > {log} 2>&1;

        ## test md5sum if available
        if [ "{params.md5}" == "''" ] ; then
            echo "WARNING: Downloaded fastq file '{params.ftp}' has no md5sum to verify the download!" | tee -a {log};
        else
            md5_actual=$(md5sum {output} | cut -d' ' -f1);
            if [ "$md5_actual" != "{params.md5}" ]; then
                echo "ERROR: Downloaded fastq file '{params.ftp}' has a wrong md5sum ($md5_actual != {params.md5})!" | tee -a {log};
                exit 1;
            fi
        fi
        """


rule get_fastq_local:
    """
    Symlink and rename local fastq file
    """
    input:
        lambda wildcards: ancient(get_fastq_name_original(wildcards)),   ## raw data: rerun only if missing, ignore its time stamp
    output:
        maybe_temp("{folder}/01_fastq/00_reads/00_files_orig/{sm}/{lb}/{idd}.fastq.gz"),
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
        maybe_temp("{folder}/00_reference/{genome}/{genome}.fasta"),
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
        ancient("{folder}/01_fastq/00_reads/00_files_orig/{sm}/{lb}/{idd}.fastq.gz"),   ## symlink to the raw data: ignore its time stamp
    output:
        maybe_temp("{folder}/01_fastq/00_reads/01_files_subsample/{sm}/{lb}/{idd}.fastq.gz", "subsampling"),
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
        "{folder}/01_fastq/00_reads/01_files_subsample/{sm}/{lb}/{idd}.fastq.gz.log",
    shell:
        """
        seqtk sample {params.params} {input} {params.number} | gzip > {output}
        """

##########################################################################################
# Split fastq into chunks
## checkpoint: output may not be temporal!
checkpoint fastq_split:
    """
    Split a single fastq file into chunks of 'nb_reads' reads to be separately computed in parallel.
    R1 and R2 of paired-end data are split in separate jobs; as both use the same 'nb_reads', chunks stay in sync.
    Output: {idd}/{id}_chunk0000_R1.fastq.gz, chunk0001..., (nb_reads <= 0: a single chunk symlinked to the input)
    """
    input:
        get_fastq_4_splitting,
    output:
        dir=temp(directory(
            "{folder}/01_fastq/00_reads/01_files_split/{sm}/{lb}/{idd}"
        )),
        chunks="{folder}/01_fastq/00_reads/01_files_split/{sm}/{lb}/{idd}_chunks.txt",
    params:
        nb_reads=lambda wildcards: int(
            get_paramGrp(["fastq_split", "nb_reads"], 0, wildcards)
        ),
    threads: get_threads("fastq_split", 1)
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("fastq_split", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("fastq_split", attempt, 4),
    message:
        "--- SPLIT FASTQ FILE(S) {input}"
    log:
        "{folder}/01_fastq/00_reads/01_files_split/{sm}/{lb}/{idd}.log",
    shell:
        """
        mkdir -p {output.dir}
        if [ {params.nb_reads} -le 0 ]; then
            echo {input} > {output.chunks}
        else
            ext=$(echo {wildcards.idd} | sed -n 's/.*\\(_R[12]\\)$/\\1/p')
            id=$(echo {wildcards.idd} | sed 's/_R[12]$//')
            zcat {input} | split -l $(( 4 * {params.nb_reads} )) -d -a 4 \
                --additional-suffix=$ext.fastq.gz \
                --filter='gzip > $FILE' - {output.dir}/${{id}}_chunk 2> {log}

            ## list of chunk files (one per line)
            (cd {output.dir} && ls -1 *_chunk*.fastq.gz) | sed 's/.*_chunk\\([0-9]*\\).*/\\1/' > {output.chunks}
        fi
        """


def get_fastq_split_files(wildcards):
    chunk_file = checkpoints.fastq_split.get(**wildcards).output.chunks
    with open(chunk_file) as f:
        files = [os.path.basename(file) for file in f.read().splitlines()]
    return files


##########################################################################################
## trimming

rule adapterremoval_collapse:
    """
    Remove adapter and low quality bases at the edges and collapse paired-end reads
    """
    input:
        get_fastq_4_cleaning,
    output:
        R=maybe_temp("{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}{chunk}.fastq.gz", "cleaning"),
        col_R=temp("{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}{chunk}_collapsed.fastq.gz"),
        col_trunc=temp("{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}{chunk}_collapsed_truncated.fastq.gz"),
        R1=temp("{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}{chunk}_R1.fastq.gz"),
        R2=temp("{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}{chunk}_R2.fastq.gz"),
        strunc=temp("{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}{chunk}.singleton.truncated.gz"),
        disc=temp("{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}{chunk}.discarded.gz"),
        settings="{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}{chunk}.settings",
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
        "{folder}/01_fastq/01_trimmed/01_adapterremoval_collapse/{sm}/{lb}/{id}{chunk}.log",
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
        R1=maybe_temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_pe/{sm}/{lb}/{id}{chunk}_R1.fastq.gz", "cleaning"
        ),
        R2=maybe_temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_pe/{sm}/{lb}/{id}{chunk}_R2.fastq.gz", "cleaning"
        ),
        singleton=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_pe/{sm}/{lb}/{id}{chunk}.singleton.truncated.gz"
        ),
        discarded=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_pe/{sm}/{lb}/{id}{chunk}.discarded.gz"
        ),
        settings="{folder}/01_fastq/01_trimmed/01_adapterremoval_pe/{sm}/{lb}/{id}{chunk}.settings",
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
        "{folder}/01_fastq/01_trimmed/01_adapterremoval_pe/{sm}/{lb}/{id}{chunk}.log",
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
        R=maybe_temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_se/{sm}/{lb}/{id}{chunk}.fastq.gz", "cleaning"
        ),
        discard=temp(
            "{folder}/01_fastq/01_trimmed/01_adapterremoval_se/{sm}/{lb}/{id}{chunk}.discarded.gz"
        ),
        settings="{folder}/01_fastq/01_trimmed/01_adapterremoval_se/{sm}/{lb}/{id}{chunk}.settings",
    resources:
        mem_mb=lambda wildcards, attempt: get_memory_alloc("cleaning", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("cleaning", attempt, 24),
    params:
        lambda wildcards: get_paramGrp(
            ["cleaning", "params_adapterremoval"],
            "--minlength 30 --trimns --trimqualities",
            wildcards,
        ),
    log: "{folder}/01_fastq/01_trimmed/01_adapterremoval_se/{sm}/{lb}/{id}{chunk}.log",
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
        adapters="{folder}/01_fastq/01_trimmed/01_adapterremoval_infer_adapters/{sm}/{lb}/{id}{chunk}.txt",
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
        "{folder}/01_fastq/01_trimmed/01_adapterremoval_infer_adapters/{sm}/{lb}/{id}{chunk}.log",
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
        R=maybe_temp("{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}{chunk}.fastq.gz", "cleaning"),
        R_merged=temp(
            "{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}{chunk}_merged.fastq.gz"
        ),
        R1=temp("{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}{chunk}_R1.fastq.gz"),
        R2=temp("{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}{chunk}_R2.fastq.gz"),
        json="{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}{chunk}.json",
        html="{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}{chunk}.html",
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
        "{folder}/01_fastq/01_trimmed/02_fastp_pe/{sm}/{lb}/{id}{chunk}.log",
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
        R1=maybe_temp("{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}{chunk}_R1.fastq.gz", "cleaning"),
        R2=maybe_temp("{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}{chunk}_R2.fastq.gz", "cleaning"),
        json="{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}{chunk}.json",
        html="{folder}/01_fastq/01_trimmed/01_fastp_pe/{sm}/{lb}/{id}{chunk}.html",
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
        "{folder}/01_fastq/01_trimmed/02_fastp_pe/{sm}/{lb}/{id}{chunk}.log",
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
        R=maybe_temp("{folder}/01_fastq/01_trimmed/01_fastp_se/{sm}/{lb}/{id}{chunk}.fastq.gz", "cleaning"),
        json="{folder}/01_fastq/01_trimmed/01_fastp_se/{sm}/{lb}/{id}{chunk}.json",
        html="{folder}/01_fastq/01_trimmed/01_fastp_se/{sm}/{lb}/{id}{chunk}.html",
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
        "{folder}/01_fastq/01_trimmed/01_fastp_se/{sm}/{lb}/{id}{chunk}.log",
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
