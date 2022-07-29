# import subprocess, os.path


# -----------------------------------------------------------------------------#
# if run_imputation:
# Getting some values specified in the config file
# ref_genome is the first defined genome (recursive_get(["genome", genome[0], "fasta"], "")

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

# Each GLIMPSE_phase step runs in about 1 minute (human genome and defaults at least).
# Run at most n = num_imputations commands in one job
# Reasoning:
#   The human genome is broken in approx. 2000 chunks.
#   GLIMPSE_phase is run for each block.
#   Unfortunately, snakemake might take a long time to infer the DAG with so many jobs.
#   2000 jobs is still fine for 1 individual, but if imputating more individuals it might be worth
#   to group a few GLIMPSE_phase commands in a single job, as they are usually fast (1-2 minutes each)
num_imputations = int(recursive_get(["imputation", "num_imputations"], 1))

# Imputation will be run by default on all chromosomes. The paramter below allow to select a subset of chromosomes.
chromosomes = list(
    map(
        str,
        eval_to_list(
            recursive_get(
                ["imputation", "chromosomes"],
                [],
            )
        ),
    )
)
if not chromosomes:
    chromosomes = get_chromosome_nams_of_genome(genome[0])
else:
    if valid_chromsome_names(genome[0], chromosomes):
        LOGGER.error(
            f"ERROR: In config[imputation][chromosomes], the following chromsome names are not recognized: {valid_chromsome_names(GENOME, chromosomes)}!"
        )
        os._exit(1)


# -----------------------------------------------------------------------------#
# Some useful functions

# This function will be useful later to know how many chunks will be merged per chromosome
def get_num_chunks(wildcards, return_str=False):
    chunk_file = checkpoints.split_genome.get(**wildcards).output[0]
    myCommand = f"wc -l {chunk_file}"
    proc = subprocess.Popen(myCommand, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    n_chunks = int(out.decode().split()[0]) + 1
    if return_str:
        n_chunks = str(n_chunks)
    return n_chunks


# Return either
# a list with formatted paths of the imputed groups (with num_imputations per group) or
# the start and end indexes of jobs for imputation
def group_chunks(wildcards, num_imputations, start_end=False):
    n_chunks = get_num_chunks(wildcards)
    chunk_list = [i for i in range(1, n_chunks)]
    groups = [
        chunk_list[i : i + num_imputations]
        for i in range(0, len(chunk_list), num_imputations)
    ]
    if start_end:
        group_number = int(wildcards.g) - 1
        start = groups[group_number][0]
        end = groups[group_number][-1]
        return [start, end]
    else:
        n_groups = len(groups) + 1
        return n_groups


# -----------------------------------------------------------------------------#
# rules definitions
localrules:
    symlink_panel,


if run_imputation:

    wildcard_constraints:
        gp="|".join([str(value) for value in gp]),
        chr="|".join([str(chr) for chr in chromosomes]),
        sm="|".join([sm for sm in samples]),


# rule aLL:
#    input:
#        phased_bcf=[
#            f"GLIMPSE_phased/{sm}.GP{GP}.phased.{ext}"
#            for sm in samples
#            for GP in gp
#            for ext in ["bcf", "bcf.csi"]
#        ],


# Split the genome into chunks
checkpoint split_genome:
    """
    Split genome into chunks to be separately imputed
    """   
    input:
        recursive_get(["imputation", "path_panel"], ""),  # some_phased_haplotypes_chr{chr}.vcf.gz
    output:
        chunks=temp("{folder}/03_sample/04_imputed/02_chunks/chunks_chr{chr}.txt"),
    params:
        params=recursive_get(["imputation", "glimse_chunk_params"], ""),
    message:
        "--- IMPUTATION: split genome (chr {wildcards.chr})"
    log:
        "{folder}/03_sample/04_imputed/log/02_chunks/chunks_chr{chr}.log",
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


rule symlink_panel:
    """
    Symlink panel and its index (or create index if not present)
    """   
    input:
        recursive_get(["imputation", "path_panel"], "")
    output:
        panel_vcf=temp("{folder}/03_sample/04_imputed/01_panel/chr{chr}.vcf.gz"),
        panel_vcf_csi=temp("{folder}/03_sample/04_imputed/01_panel/chr{chr}.vcf.gz.csi"),
    message:
        "--- IMPUTATION: symlink panel (chr {wildcards.chr})"
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



# extract variable positions
rule extract_positions:
    """
    Extract variable positions
    """
    input:
        vcf="{folder}/03_sample/04_imputed/01_panel/chr{chr}.vcf.gz",
        index="{folder}/03_sample/04_imputed/01_panel/chr{chr}.vcf.gz.csi",
    output:
        sites=temp("{folder}/03_sample/04_imputed/03_sites/chr{chr}_sites.vcf.gz"),
        tsv=temp("{folder}/03_sample/04_imputed/03_sites/chr{chr}_sites.tsv.gz"),
        csi=temp("{folder}/03_sample/04_imputed/03_sites/chr{chr}_sites.vcf.gz.csi"),
        tbi=temp("{folder}/03_sample/04_imputed/03_sites/chr{chr}_sites.tsv.gz.tbi"),
    message:
        "--- IMPUTATION: extract positions (chr {wildcards.chr})"
    log:
        "{folder}/03_sample/04_imputed/log/03_sites/{chr}_sites.log",
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


# rule get_bam_stats:
#    input:
#        bam=lambda wildcards: samples_files[wildcards.sm],
#        bai=lambda wildcards: bam2bai(samples_files[wildcards.sm]),
#    output:


rule bcftools_mpileup:
    input:
        genome=f"{{folder}}/00_reference/{genome[0]}/{genome[0]}.fasta",
        bam=f"{{folder}}/03_sample/03_final_sample/01_bam/{{sm}}.{genome[0]}.bam",
        bai=f"{{folder}}/03_sample/03_final_sample/01_bam/{{sm}}.{genome[0]}.bai",
        sites="{folder}/03_sample/04_imputed/03_sites/chr{chr}_sites.vcf.gz",
        tsv="{folder}/03_sample/04_imputed/03_sites/chr{chr}_sites.tsv.gz",
    output:
        final_vcf=temp("{folder}/03_sample/04_imputed/04_vcf/{sm}_chr{chr}.vcf.gz"),
        final_csi=temp("{folder}/03_sample/04_imputed/04_vcf/{sm}_chr{chr}.vcf.gz.csi"),
    threads: 
        get_threads2(["imputation", "bcftools_mpileup"], 1)
    log:
        "{folder}/03_sample/04_imputed/log/04_vcf/{sm}_chr{chr}.log",
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
        bcftools mpileup -f {input.genome} \
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


# # 5. Impute and phase a whole chromosome
# as each job runs for ~1 minute, we can group ~10 jobs and send them as one
rule impute_phase:
    input:
        chunks="{folder}/03_sample/04_imputed/02_chunks/chunks_chr{chr}.txt",
        vcf_sample="{folder}/03_sample/04_imputed/04_vcf/{sm}_chr{chr}.vcf.gz",
        csi_sample="{folder}/03_sample/04_imputed/04_vcf/{sm}_chr{chr}.vcf.gz.csi",
        vcf_ref="{folder}/03_sample/04_imputed/01_panel/chr{chr}.vcf.gz",
        csi_ref="{folder}/03_sample/04_imputed/01_panel/chr{chr}.vcf.gz.csi",
        map=recursive_get(["imputation", "path_map"], ""),
    output:
        txt=temp(
            "{folder}/03_sample/04_imputed/05_GLIMPSE_imputed/done_{sm}_chr{chr}_group{g}.txt"
        ),
    message:
        "--- IMPUTATION: impute phase (sample {wildcards.sm}; chr {wildcards.chr})"
    params:
        imputed_bcf=lambda wildcards: f"{wildcards.folder}/03_sample/04_imputed/05_GLIMPSE_imputed/{wildcards.sm}_chr{wildcards.chr}.${{n}}.bcf",
        start_end=lambda wildcards: group_chunks(
            wildcards, num_imputations=num_imputations, start_end=True
        ),
        params=recursive_get(["imputation", "glimse_phase_params"], ""),
    log:
        "{folder}/03_sample/04_imputed/log/05_GLIMPSE_imputed/{sm}_chr{chr}_group{g}.log",
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
        for n in $(seq {params.start_end[0]} {params.start_end[1]})
        do
            LINE=$( head -n ${{n}} {input.chunks} |tail -n1 );
            IRG=$(echo $LINE | cut -d" " -f3);
            ORG=$(echo $LINE | cut -d" " -f4);

            GLIMPSE_phase {params.params} \
                --input {input.vcf_sample} \
                --reference {input.vcf_ref} \
                --map {input.map} \
                --input-region ${{IRG}} \
                --output-region ${{ORG}} \
                --output {params.imputed_bcf} \
                --thread {threads};

            bcftools index -f {params.imputed_bcf};
        done

        touch {output.txt};
        """


# #  # 6. Ligate multiple chunks together
## input.chunks is needed for the function get_num_chunks()
rule ligate_chunks:
    input:
        txt=lambda wildcards: [
            f"{wildcards.folder}/03_sample/04_imputed/05_GLIMPSE_imputed/done_{wildcards.sm}_chr{wildcards.chr}_group{g}.txt"
            for g in range(1, group_chunks(wildcards, num_imputations))
        ],
        chunks="{folder}/03_sample/04_imputed/02_chunks/chunks_chr{chr}.txt",
    output:
        ligated_bcf=temp(
            "{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}_chr{chr}.merged.bcf"
        ),
        ligated_csi=temp(
            "{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}_chr{chr}.merged.bcf.csi"
        ),
        list_files=temp(
            "{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/list_{sm}_chr{chr}.txt"
        ),
    message:
        "--- IMPUTATION: ligate chunks (sample {wildcards.sm}; chr {wildcards.chr})"
    log:
        "{folder}/03_sample/04_imputed/log/06_GLIMPSE_ligated/{sm}_chr{chr}.log",
    params:
        phased=lambda wildcards: [
            f"{wildcards.folder}/03_sample/04_imputed/05_GLIMPSE_imputed/{wildcards.sm}_chr{wildcards.chr}.{n}.bcf"
            for n in range(1, get_num_chunks(wildcards))
        ],
        phased_csi=lambda wildcards: [
            f"{wildcards.folder}/03_sample/04_imputed/05_GLIMPSE_imputed/{wildcards.sm}_chr{wildcards.chr}.{n}.bcf.csi"
            for n in range(1, get_num_chunks(wildcards))
        ],
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
        list_files={output.ligated_bcf}.list_files
        for file in {params.phased} ; do echo $file ; done > ${{list_files}};
        GLIMPSE_ligate --input ${{list_files}}; --output {output.ligated_bcf};
        rm {params.phased} {params.phased_csi} ${{list_files}};;
        bcftools index -f {output.ligated_bcf};
        """


# #  # 6.1
# concatenate imputed genotypes per sample
rule concat:
    input:
        bcf=lambda wildcards: expand(
            "{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}_chr{chr}.merged.bcf",
            chr=chromosomes,
            sm=wildcards.sm,
            allow_missing=True,
        ),
        index=lambda wildcards: expand(
            "{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}_chr{chr}.merged.bcf.csi",
            chr=chromosomes,
            sm=wildcards.sm,
            allow_missing=True,
        ),
    output:
        bcf=temp("{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}.bcf"),
        csi=temp("{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}.bcf.csi"),
    message:
        "--- IMPUTATION: concatenate (sample {wildcards.sm})"
    log:
        "{folder}/03_sample/04_imputed/log/06_GLIMPSE_ligated/{sm}.log",
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
        bcf="{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}.bcf",
        index="{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}.bcf.csi",
    output:
        bcf=temp("{folder}/03_sample/04_imputed/07_unphased/{sm}.GP{gp}.bcf"),
        csi=temp("{folder}/03_sample/04_imputed/07_unphased/{sm}.GP{gp}.bcf.csi"),
    message:
        "--- IMPUTATION: GP filtering (sample {wildcards.sm}; GP {wildcards.gp})"
    log:
        "{folder}/03_sample/04_imputed/log/07_unphased/{sm}_GP{gp}.log",
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
        bcf="{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}.bcf",
        index="{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}.bcf.csi",
    output:
        table=report(
            "{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}_gp.txt",
            caption="../report/imputation_gp_table.rst",
            category="Imputation",
            subcategory="Tables",
        ),
    message:
        "--- IMPUTATION: get GP values (sample {wildcards.sm})"
    log:
        "{folder}/03_sample/04_imputed/log/06_GLIMPSE_ligated/{sm}_GP_values.log",
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
        "{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}_gp.txt",
        to_trigger_rerun=[
            f"{{folder}}/03_sample/04_imputed/07_unphased/{{sm}}.GP{gp}.bcf"
            for gp in str2list(
                recursive_get(["imputation", "gp_filter"], [0.8])
            )
        ],
    output:
        report(
            "{folder}/03_sample/04_imputed/06_GLIMPSE_ligated/{sm}_gp.svg",
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
        gp=",".join(
            str2list(recursive_get(["imputation", "gp_filter"], [0.8]))
        ),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["imputation", "plot_gp"], attempt, 4
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["imputation", "plot_gp"], attempt, 1
        ),
    log:
        "{folder}/03_sample/04_imputed/log/06_GLIMPSE_ligated/{sm}_GP_values.log",
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
rule sample_haplotypes:
    input:
        ligated_bcf="{folder}/03_sample/04_imputed/07_unphased/{sm}.GP{gp}.bcf",
        ligated_csi="{folder}/03_sample/04_imputed/07_unphased/{sm}.GP{gp}.bcf.csi",
    output:
        phased_bcf="{folder}/03_sample/04_imputed/{sm}.GP{gp}.phased.bcf",
        phased_csi="{folder}/03_sample/04_imputed/{sm}.GP{gp}.phased.bcf.csi",
    message:
        "--- IMPUTATION: sample haplotypes (sample {wildcards.sm}; GP {wildcards.gp})"
    log:
        "{folder}/03_sample/04_imputed/log/{sm}.GP{gp}.log",
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
