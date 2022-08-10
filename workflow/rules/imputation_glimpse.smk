# import subprocess, os.path


# -----------------------------------------------------------------------------#
# if run_imputation:
# Getting some values specified in the config file
# ref_genome is the first defined GENOMES (recursive_get(["GENOMES", GENOMES[0], "fasta"], "")

# This string contains a wildcard where we will place the name of the chromosome
# something like "path/to/my/panel_chr{chr}.vcf.gz"
path_panel = recursive_get(["imputation", "path_panel"], "")
if path_panel == "":
    LOGGER.error(f"ERROR: Parameter config[imputation][path_panel] is not specified!")
    sys.exit(1)

path_map = recursive_get(["imputation", "path_map"], "")
if path_map == "":
    LOGGER.error(f"ERROR: Parameter config[imputation][path_map] is not specified!")
    sys.exit(1)

# Each GLIMPSE_phase step runs in about 1 minute (human GENOMES and defaults at least).
# Run at most n = num_imputations commands in one job
# Reasoning:
#   The human GENOMES is broken in approx. 2000 chunks.
#   GLIMPSE_phase is run for each block.
#   Unfortunately, snakemake might take a long time to infer the DAG with so many jobs.
#   2000 jobs is still fine for 1 individual, but if imputating more individuals it might be worth
#   to group a few GLIMPSE_phase commands in a single job, as they are usually fast (1-2 minutes each)
num_imputations = int(recursive_get(["imputation", "num_imputations"], 1))

# Imputation will be run by default on all chromosomes. The paramter below allow to select a subset of chromosomes.
chromosomes = to_list(recursive_get(["imputation", "chromosomes"], []))
if not chromosomes:
    chromosomes = get_chromosome_nams_of_genome(GENOMES[0])
else:
    if valid_chromosome_names(GENOMES[0], chromosomes):
        LOGGER.error(
            f"ERROR: In config[imputation][chromosomes], the following chromsome names are not recognized: {valid_chromosome_names(GENOMES[0], chromosomes)}!"
        )
        os._exit(1)

# This string contains a wildcard where we will place the name of the chromosome
# something like "path/to/my/panel_chr{chr}.vcf.gz"
path_panel = recursive_get(["imputation", "path_panel"], "")
for chr in chromosomes:
    file = path_panel.format(chr=chr)
    if not os.path.isfile(file):
        LOGGER.error(
            f"ERROR: Panel file config[imputation][path_panel] ({path_panel}) does not exist for 'chr={chr}'!"
        )
        sys.exit(1)


# This string contains a wildcard where we will place the name of the chromosome
# something like "path/to/my/panel_chr{chr}.vcf.gz"
path_map = recursive_get(["imputation", "path_map"], "")
for chr in chromosomes:
    file = path_map.format(chr=chr)
    if not os.path.isfile(file):
        LOGGER.error(
            f"ERROR: Map file config[imputation][path_map] ({path_panel}) does not exist for 'chr={chr}'!"
        )
        sys.exit(1)


# -----------------------------------------------------------------------------#
# Some useful functions

# This function will be useful later to know how many chunks will be merged per chromosome
def get_num_chunks(wildcards, return_str=False):
    chunk_file = checkpoints.glimpse_chunk.get(**wildcards).output[0]
    myCommand = f"wc -l {chunk_file}"
    proc = subprocess.Popen(myCommand, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    n_chunks = int(out.decode().split()[0]) + 1
    if return_str:
        n_chunks = str(n_chunks)
    return n_chunks


# -----------------------------------------------------------------------------#
# rules definitions
localrules:
    get_panel,
    get_map,


# if run_imputation:
#    gp=[0.8]
#    chromosomes=[20]
wildcard_constraints:
    #        gp="|".join([str(value) for value in gp]),
    chr="|".join([str(chr) for chr in chromosomes]),
    sm="|".join(
        [sm for sm in SAMPLES]
        + [sm for gen in EXTERNAL_SAMPLES for sm in EXTERNAL_SAMPLES[gen]]
    ),


rule get_panel:
    """
    Symlink panel and its index (or create index if not present)
    """
    input:
        recursive_get(["imputation", "path_panel"], ""),
    output:
        panel_vcf=temp(
            "{folder}/03_sample/04_imputed/01_panel/01_panel/chr{chr}.vcf.gz"
        ),
        panel_vcf_csi=temp(
            "{folder}/03_sample/04_imputed/01_panel/01_panel/chr{chr}.vcf.gz.csi"
        ),
    message:
        "--- IMPUTATION: symlink the panel (chr {wildcards.chr})"
    conda:
        "../envs/bcftools.yaml"
    envmodules:
        module_bcftools,
    shell:
        """
        ln -srf {input} {output.panel_vcf};

        ## if index does not exist create it
        if [ -f "{input}.csi" ]; then
            ln -srf {input}.csi {output.panel_vcf_csi};
        else
            bcftools index -f {output.panel_vcf} -o {output.panel_vcf_csi};
        fi
        """


rule get_map:
    """
    Symlink the map
    """
    input:
        recursive_get(["imputation", "path_map"], ""),
    output:
        temp("{folder}/03_sample/04_imputed/02_map/chr{chr}.gz"),
    message:
        "--- IMPUTATION: get the map (chr {wildcards.chr})"
    shell:
        """
        ln -srf {input} {output};
        """


#
# extract variable positions
rule extract_positions:
    """
    Extract variable positions
    """
    input:
        vcf="{folder}/03_sample/04_imputed/01_panel/01_panel/chr{chr}.vcf.gz",
        index="{folder}/03_sample/04_imputed/01_panel/01_panel/chr{chr}.vcf.gz.csi",
    output:
        sites=temp("{folder}/03_sample/04_imputed/01_panel/02_sites/chr{chr}.vcf.gz"),
        tsv=temp("{folder}/03_sample/04_imputed/01_panel/02_sites/chr{chr}.tsv.gz"),
        csi=temp("{folder}/03_sample/04_imputed/01_panel/02_sites/chr{chr}.vcf.gz.csi"),
        tbi=temp("{folder}/03_sample/04_imputed/01_panel/02_sites/chr{chr}.tsv.gz.tbi"),
    message:
        "--- IMPUTATION: extract positions (chr {wildcards.chr})"
    log:
        "{folder}/log/03_sample/04_imputed/01_panel/02_sites/{chr}.log",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["imputation", "extract_positions"], attempt, 4
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["imputation", "extract_positions"], attempt, 6
        ),
    conda:
        "../envs/bcftools.yaml"
    envmodules:
        module_bcftools,
    shell:
        """
        bcftools view -G -m 2 -M 2 -v snps \
            {input.vcf} \
            -Oz -o {output.sites};

        bcftools index -f {output.sites};
        bcftools query -f'%CHROM\\t%POS\\t%REF,%ALT\\n' {output.sites} | \
            bgzip -c > {output.tsv};
        tabix -s1 -b2 -e2 {output.tsv};
        """


rule bcftools_mpileup:
    input:
        ref=f"{{folder}}/00_reference/{GENOMES[0]}/{GENOMES[0]}.fasta",
        bam=f"{{folder}}/03_sample/03_final_sample/01_bam/{{sm}}.{GENOMES[0]}.bam",
        bai=f"{{folder}}/03_sample/03_final_sample/01_bam/{{sm}}.{GENOMES[0]}.bai",
        sites="{folder}/03_sample/04_imputed/01_panel/02_sites/chr{chr}.vcf.gz",
        tsv="{folder}/03_sample/04_imputed/01_panel/02_sites/chr{chr}.tsv.gz",
    output:
        final_vcf=temp("{folder}/03_sample/04_imputed/03_vcf/{sm}_chr{chr}.vcf.gz"),
        final_csi=temp("{folder}/03_sample/04_imputed/03_vcf/{sm}_chr{chr}.vcf.gz.csi"),
    threads: get_threads2(["imputation", "bcftools_mpileup"], 1)
    log:
        "{folder}/log/03_sample/04_imputed/03_vcf/{sm}_chr{chr}.log",
    message:
        "--- IMPUTATION: bcftools mpileup (sample {wildcards.sm}; chr {wildcards.chr})"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["imputation", "bcftools_mpileup"], attempt, 4
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["imputation", "bcftools_mpileup"], attempt, 6
        ),
    conda:
        "../envs/bcftools.yaml"
    envmodules:
        module_bcftools,
    shell:
        """
        header={output.final_vcf}.header;
        vcf_tmp={output.final_vcf}.tmp;
        bcftools mpileup -f {input.ref} \
            -I -E -a 'FORMAT/DP' \
            -T {input.sites} -r {wildcards.chr} \
            -Ob --threads {threads} \
            --ignore-RG {input.bam}  | \
            bcftools call -Aim -C alleles -T {input.tsv} -Oz -o ${{vcf_tmp}};
        echo {input.bam} {wildcards.sm} > ${{header}};
        bcftools reheader -s ${{header}} -o {output.final_vcf} ${{vcf_tmp}};
        bcftools index -f {output.final_vcf};
        rm -f ${{header}} ${{vcf_tmp}};
        """


# Split the GENOMES into chunks
## checkpoint: output may not be temporal!
checkpoint glimpse_chunk:
    """
    Split GENOMES into chunks to be separately imputed
    """
    input:
        "{folder}/03_sample/04_imputed/01_panel/01_panel/chr{chr}.vcf.gz",
    output:
        chunks="{folder}/03_sample/04_imputed/04_glimpse_chunked/chunks_chr{chr}.txt",
    params:
        params=recursive_get(["imputation", "glimse_chunk_params"], ""),
    message:
        "--- GLIMSE_CHUNK: split GENOMES (chr {wildcards.chr})"
    log:
        "{folder}/log/03_sample/04_imputed/04_glimpse_chunked/chunks_chr{chr}.log",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["imputation", "split_genome"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["imputation", "split_genome"], attempt, 4
        ),
    conda:
        "../envs/glimpse.yaml"
    envmodules:
        module_glimpse,
    shell:
        """
        GLIMPSE_chunk --input {input} \
            {params.params} \
            --output {output.chunks}
        """


# # 5. Impute and phase a whole chromosome
# as each job runs for ~1 minute, we can group ~10 jobs and send them as one
rule glimpse_phase:
    input:
        vcf_ref="{folder}/03_sample/04_imputed/01_panel/01_panel/chr{chr}.vcf.gz",
        csi_ref="{folder}/03_sample/04_imputed/01_panel/01_panel/chr{chr}.vcf.gz.csi",
        chunks="{folder}/03_sample/04_imputed/04_glimpse_chunked/chunks_chr{chr}.txt",
        vcf_sample="{folder}/03_sample/04_imputed/03_vcf/{sm}_chr{chr}.vcf.gz",
        csi_sample="{folder}/03_sample/04_imputed/03_vcf/{sm}_chr{chr}.vcf.gz.csi",
        map="{folder}/03_sample/04_imputed/02_map/chr{chr}.gz",
    output:
        bcf=temp(
            "{folder}/03_sample/04_imputed/05_glimpse_phased/{sm}/chr{chr}/chunk{n}.bcf"
        ),
        csi=temp(
            "{folder}/03_sample/04_imputed/05_glimpse_phased/{sm}/chr{chr}/chunk{n}.bcf.csi"
        ),
    message:
        "--- IMPUTATION: impute phase (sample {wildcards.sm}; chr {wildcards.chr}; chunk: {wildcards.n})"
    params:
        params=recursive_get(["imputation", "glimse_phase_params"], ""),
    log:
        "{folder}/log/03_sample/04_imputed/05_glimpse_phased/{sm}/chr{chr}/chunk{n}.log",
    threads: get_threads2(["imputation", "impute_phase"], 1)
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["imputation", "impute_phase"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["imputation", "impute_phase"], attempt, 2
        ),
    conda:
        "../envs/glimpse.yaml"
    envmodules:
        module_glimpse,
        module_bcftools,
    shell:
        """
        ## get chunk
        LINE=$( head -n {wildcards.n} {input.chunks} |tail -n1 );

        ## run glimpse
        GLIMPSE_phase {params.params} \
            --input {input.vcf_sample} \
            --reference {input.vcf_ref} \
            --map {input.map} \
            --input-region $(echo $LINE | cut -d" " -f3) \
            --output-region $(echo $LINE | cut -d" " -f4) \
            --output {output.bcf} \
            --thread {threads};

        ## index bcf file
        bcftools index -f {output.bcf};
        """


# #  # 6. Ligate multiple chunks together
rule glimpse_ligate:
    input:
        bcf=lambda wildcards: [
            f"{wildcards.folder}/03_sample/04_imputed/05_glimpse_phased/{wildcards.sm}/chr{wildcards.chr}/chunk{g}.bcf"
            for g in range(1, get_num_chunks(wildcards))
        ],
        csi=lambda wildcards: [
            f"{wildcards.folder}/03_sample/04_imputed/05_glimpse_phased/{wildcards.sm}/chr{wildcards.chr}/chunk{g}.bcf.csi"
            for g in range(1, get_num_chunks(wildcards))
        ],
        chunks="{folder}/03_sample/04_imputed/04_glimpse_chunked/chunks_chr{chr}.txt",
    output:
        ligated_bcf=temp(
            "{folder}/03_sample/04_imputed/06_glimpse_ligated/{sm}/chr{chr}.bcf"
        ),
        ligated_csi=temp(
            "{folder}/03_sample/04_imputed/06_glimpse_ligated/{sm}/chr{chr}.bcf.csi"
        ),
    message:
        "--- IMPUTATION: ligate chunks (sample {wildcards.sm}; chr {wildcards.chr})"
    log:
        "{folder}/log/03_sample/04_imputed/06_glimpse_ligated/{sm}/chr{chr}.log",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["imputation", "ligate_chunks"], attempt, 4
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["imputation", "ligate_chunks"], attempt, 6
        ),
    conda:
        "../envs/glimpse.yaml"
    envmodules:
        module_bcftools,
        module_glimpse,
    shell:
        """
        list_files={output.ligated_bcf}.list_files;
        for file in {input.bcf}; do echo $file; done > ${{list_files}};
        GLIMPSE_ligate --input ${{list_files}} --output {output.ligated_bcf};
        rm ${{list_files}};
        bcftools index -f {output.ligated_bcf};
        """


# #  # 6.1
# concatenate imputed genotypes per sample
rule concat:
    """
    Concatenate imputed genotypes per sample
    """
    input:
        bcf=lambda wildcards: expand(
            "{folder}/03_sample/04_imputed/06_glimpse_ligated/{sm}/chr{chr}.bcf",
            chr=chromosomes,
            sm=wildcards.sm,
            allow_missing=True,
        ),
        index=lambda wildcards: expand(
            "{folder}/03_sample/04_imputed/06_glimpse_ligated/{sm}/chr{chr}.bcf.csi",
            chr=chromosomes,
            sm=wildcards.sm,
            allow_missing=True,
        ),
    output:
        bcf="{folder}/03_sample/04_imputed/06_glimpse_ligated/{sm}.bcf",
        csi="{folder}/03_sample/04_imputed/06_glimpse_ligated/{sm}.bcf.csi",
    message:
        "--- IMPUTATION: concatenate chromosomes (sample {wildcards.sm})"
    log:
        "{folder}/log/03_sample/04_imputed/06_glimpse_ligated/{sm}.log",
    conda:
        "../envs/bcftools.yaml"
    envmodules:
        module_bcftools,
    shell:
        """
        bcftools concat {input.bcf} -Oz -o {output.bcf};
        bcftools index -f {output.bcf};
        """


# # # 6.2 Filter by GP
rule filter_gp_split_sm:
    input:
        bcf="{folder}/03_sample/04_imputed/06_glimpse_ligated/{sm}.bcf",
        index="{folder}/03_sample/04_imputed/06_glimpse_ligated/{sm}.bcf.csi",
    output:
        bcf=temp(
            "{folder}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp{gp}.bcf"
        ),
        csi=temp(
            "{folder}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp{gp}.bcf.csi"
        ),
    message:
        "--- IMPUTATION: GP filtering (sample {wildcards.sm}; GP {wildcards.gp})"
    log:
        "{folder}/log/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp{gp}.log",
    conda:
        "../envs/bcftools.yaml"
    envmodules:
        module_bcftools,
    shell:
        """
        bcftools view -i'FORMAT/GP[0:*]>{wildcards.gp}' -Ob -o {output.bcf} {input.bcf};
        bcftools index -f {output.bcf};
        """


rule get_gp:
    input:
        bcf="{folder}/03_sample/04_imputed/06_glimpse_ligated/{sm}.bcf",
        index="{folder}/03_sample/04_imputed/06_glimpse_ligated/{sm}.bcf.csi",
    output:
        table=report(
            "{folder}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp.txt",
            caption="../report/imputation_gp_table.rst",
            category="Imputation",
            subcategory="Tables",
        ),
    message:
        "--- IMPUTATION: get GP values (sample {wildcards.sm})"
    log:
        "{folder}/log/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp.log",
    conda:
        "../envs/bcftools.yaml"
    envmodules:
        module_bcftools,
    shell:
        """
        bcftools query -f '[%GT]\\t[%GP]\\n' {input.bcf} > {output.table};
        """


rule plot_gp:
    input:
        "{folder}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp.txt",
        to_trigger_rerun=[
            f"{{folder}}/03_sample/04_imputed/07_glimpse_sampled/unphased/{{sm}}_gp{gp}.bcf"
            for gp in str2list(recursive_get(["imputation", "gp_filter"], [0.8]))
        ],
    output:
        report(
            "{folder}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp.svg",
            caption="../report/imputation_gp_plot.rst",
            category="Imputation",
            subcategory="Plots",
        ),
    message:
        "--- IMPUTATION: plot GP values (sample {wildcards.sm})"
    params:
        script=workflow.source_path("../scripts/plot_imputation_gp.R"),
        width=recursive_get(["plots", "width"], 11),
        height=recursive_get(["plots", "height"], 7),
        gp=",".join(str2list(recursive_get(["imputation", "gp_filter"], [0.8]))),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["imputation", "plot_gp"], attempt, 4
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["imputation", "plot_gp"], attempt, 1
        ),
    log:
        "{folder}/log/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp_plot.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    shell:
        """
        Rscript {params.script} \
            --input={input} \
            --output={output} \
            --width={params.width} \
            --height={params.height} \
            --gp={params.gp} \
            --sample={wildcards.sm}
        """


# #  # 7. Sample haplotypes
rule glimpse_sample:
    input:
        ligated_bcf="{folder}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp{gp}.bcf",
        ligated_csi="{folder}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp{gp}.bcf.csi",
    output:
        phased_bcf="{folder}/03_sample/04_imputed/07_glimpse_sampled/{sm}_gp{gp}.bcf",
        phased_csi="{folder}/03_sample/04_imputed/07_glimpse_sampled/{sm}_gp{gp}.bcf.csi",
    message:
        "--- IMPUTATION: sample haplotypes (sample {wildcards.sm}; GP {wildcards.gp})"
    log:
        "{folder}/log/03_sample/04_imputed/07_glimpse_sampled/{sm}_gp{gp}.log",
    conda:
        "../envs/glimpse.yaml"
    envmodules:
        module_glimpse,
    shell:
        """
        GLIMPSE_sample --input {input.ligated_bcf} \
            --solve \
            --output {output.phased_bcf};
        bcftools index -f {output.phased_bcf};
        """
