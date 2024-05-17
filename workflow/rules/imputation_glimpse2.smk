# This function will be useful later to know how many chunks will be merged per chromosome
def get_num_chunks2(wildcards, return_str=False):
    chunk_file = checkpoints.glimpse2_chunk.get(**wildcards).output[0]
    myCommand = f"wc -l {chunk_file}"
    proc = subprocess.Popen(myCommand, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    n_chunks = int(out.decode().split()[0]) + 1
    if return_str:
        n_chunks = str(n_chunks)
    return n_chunks


## 3. Split the genome into chunks
checkpoint glimpse2_chunk:
    """
    Split GENOMES into chunks to be separately imputed
    """
    input:
        #map="{folder}/03_sample/04_imputed/02_map/{genome}/chr{chr}.gz",
        map=lambda wildcards: get_param(
            ["imputation", wildcards.genome, "path_map"], ""
        ),
        vcf="{folder}/03_sample/04_imputed/01_panel/01_panel/{genome}/chr{chr}.vcf.gz",
        csi="{folder}/03_sample/04_imputed/01_panel/01_panel/{genome}/chr{chr}.vcf.gz.csi",
    output:
        chunks="{folder}/03_sample/04_imputed/04_glimpse2_chunked/{genome}/chunks_chr{chr}.txt",
    params:
        params=lambda wildcards: get_param(
            ["imputation", wildcards.genome, "glimpse2_chunk_params"], ""
        ),
    threads: lambda wildcards: get_threads2(["imputation", wildcards.genome, "glimpse2_chunk"], 1)
    message:
        "--- GLIMPSE2_CHUNK: split GENOMES (genome {wildcards.genome}; chr {wildcards.chr})"
    log:
        "{folder}/log/03_sample/04_imputed/04_glimpse2_chunked/{genome}/chunks_chr{chr}.log",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["imputation2", wildcards.genome, "glimpse2_chunk"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["imputation2", wildcards.genome, "glimpse2_chunk"], attempt, 4
        ),
    conda:
        "../envs/glimpse2.yaml"
    envmodules:
        module_glimpse2,
    shell:
        """
        GLIMPSE2_chunk {params.params} \
            --input {input.vcf} \
            --threads {threads} \
            --region {wildcards.chr} \
            --output {output.chunks} \
            --map {input.map} \
            --log {log}
        """


## 4. Create binary reference panel
rule glimpse2_split_reference:
    """
    Create a binary reference panel for quick reading time
    """
    input:
        vcf_ref="{folder}/03_sample/04_imputed/01_panel/01_panel/{genome}/chr{chr}.vcf.gz",
        csi_ref="{folder}/03_sample/04_imputed/01_panel/01_panel/{genome}/chr{chr}.vcf.gz.csi",
        chunks="{folder}/03_sample/04_imputed/04_glimpse2_chunked/{genome}/chunks_chr{chr}.txt",
        map=lambda wildcards: get_param(
            ["imputation", wildcards.genome, "path_map"], ""
        ),
    output:
        temp(
            "{folder}/03_sample/04_imputed/05_glimpse2_splitted_reference/{genome}/chr{chr}/chunk{n}.n"
        ),
    message:
        "--- IMPUTATION GLIMPSE2: create binary reference panel (genome {wildcards.genome}; chr {wildcards.chr}; chunk: {wildcards.n})"
    threads: lambda wildcards: get_threads2(
    ["imputation", wildcards.genome, "glimpse2_split_reference"], 1
)
    params:
        params=lambda wildcards: get_param(
            ["imputation", wildcards.genome, "GLIMPSE2_split_reference_params"], ""
        ),
    log:
        "{folder}/log/03_sample/04_imputed/05_glimpse2_splitted_reference/{genome}/chr{chr}/chunk{n}.log",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["imputation2", wildcards.genome, "glimpse2_split_reference"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["imputation2", wildcards.genome, "glimpse2_split_reference"], attempt, 4
        ),
    conda:
        "../envs/glimpse2.yaml"
    envmodules:
        module_glimpse2,
    shell:
        """
        ## get chunk
        LINE=$( head -n {wildcards.n} {input.chunks} |tail -n1 );

        ## run glimpse
        GLIMPSE2_split_reference {params.params} \
            --reference {input.vcf_ref} \
            --threads {threads} \
            --map {input.map}  \
            --input-region $(echo $LINE | cut -d" " -f3) \
            --output-region $(echo $LINE | cut -d" " -f4) \
            --output {output};

        ## rename output file
        mv {output}_$(echo $LINE | cut -d" " -f3 | sed 's/:/_/g' | sed 's/-/_/g').bin {output};
        """


## 5. Impute and phase a whole chromosome
rule glimpse2_phase:
    input:
        bam="{folder}/03_sample/03_final_sample/01_bam/{sm}.{genome}.bam",
        bai="{folder}/03_sample/03_final_sample/01_bam/{sm}.{genome}.bai",
        ref="{folder}/03_sample/04_imputed/05_glimpse2_splitted_reference/{genome}/chr{chr}/chunk{n}.n",
    output:
        bcf=temp(
            "{folder}/03_sample/04_imputed/05_glimpse2_phased/{sm}.{genome}/chr{chr}/chunk{n}.bcf"
        ),
        csi=temp(
            "{folder}/03_sample/04_imputed/05_glimpse2_phased/{sm}.{genome}/chr{chr}/chunk{n}.bcf.csi"
        ),
    message:
        "--- IMPUTATION GLIMPSE2: impute phase (sample {wildcards.sm}; genome {wildcards.genome}; chr {wildcards.chr}; chunk: {wildcards.n})"
    params:
        params=lambda wildcards: get_param(
            ["imputation", wildcards.genome, "glimpse2_phase_params"], ""
        ),
    log:
        "{folder}/03_sample/04_imputed/05_glimpse2_phased/{sm}.{genome}/chr{chr}/chunk{n}.log",
    threads: lambda wildcards: get_threads2(["imputation", wildcards.genome, "glimpse2_phase"], 1)
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["imputation", wildcards.genome, "glimpse2_phase"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["imputation", wildcards.genome, "glimpse2_phase"], attempt, 2
        ),
    conda:
        "../envs/glimpse2.yaml"
    envmodules:
        module_glimpse2,
        module_bcftools,
    shell:
        """
        ## run glimpse
        GLIMPSE2_phase {params.params} \
            --bam-file {input.bam} \
            --reference {input.ref} \
            --output {output.bcf} \
            --threads {threads};

        ## index bcf file
        bcftools index -f {output.bcf};
        """


# 6. Ligate chunks of the the same chromosome
rule glimpse2_ligate:
    input:
        bcf=lambda wildcards: [
            f"{wildcards.folder}/03_sample/04_imputed/05_glimpse2_phased/{wildcards.sm}.{wildcards.genome}/chr{wildcards.chr}/chunk{g}.bcf"
            for g in range(1, get_num_chunks2(wildcards))
        ],
        csi=lambda wildcards: [
            f"{wildcards.folder}/03_sample/04_imputed/05_glimpse2_phased/{wildcards.sm}.{wildcards.genome}/chr{wildcards.chr}/chunk{g}.bcf.csi"
            for g in range(1, get_num_chunks2(wildcards))
        ],
        chunks="{folder}/03_sample/04_imputed/04_glimpse2_chunked/{genome}/chunks_chr{chr}.txt",
    output:
        ligated_bcf=temp(
            "{folder}/03_sample/04_imputed/06_glimpse2_ligated/{sm}.{genome}/chr{chr}.bcf"
        ),
        ligated_csi=temp(
            "{folder}/03_sample/04_imputed/06_glimpse2_ligated/{sm}.{genome}/chr{chr}.bcf.csi"
        ),
    message:
        "--- IMPUTATION: ligate chunks (sample {wildcards.sm}; genome: {wildcards.genome}; chr {wildcards.chr})"
    log:
        "{folder}/log/03_sample/04_imputed/06_glimpse2_ligated/{sm}.{genome}/chr{chr}.log",
    params:
        params=lambda wildcards: get_param(
            ["imputation", wildcards.genome, "glimpse2_ligate_params"], ""
        ),
    threads: lambda wildcards: get_threads2(["imputation", wildcards.genome, "glimpse2_ligate"], 1)
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["imputation", wildcards.genome, "glimpse2_ligate"], attempt, 4
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["imputation", wildcards.genome, "glimpse2_ligate"], attempt, 6
        ),
    conda:
        "../envs/glimpse2.yaml"
    envmodules:
        module_bcftools,
        module_glimpse2,
    shell:
        """
        list_files={output.ligated_bcf}.list_files;
        for file in {input.bcf}; do echo $file; done > ${{list_files}};
        GLIMPSE2_ligate {params} \
            --threads {threads} \
            --input ${{list_files}} \
            --output {output.ligated_bcf};
        rm ${{list_files}};
        bcftools index -f {output.ligated_bcf};
        """
