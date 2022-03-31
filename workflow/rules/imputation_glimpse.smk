import subprocess, os.path

# -----------------------------------------------------------------------------#
# Getting some values specified in the config file
ref_genome = config["imputation_GLIMPSE"]["ref_genome"]

# This string contains a wildcard where we will place the name of the chromosome
# something like "path/to/my/panel_chr{chr}.vcf.gz"
path_panel = config.get("imputation_GLIMPSE", {}).get("path_panel")

# Each GLIMPSE_phase step runs in about 1 minute (human genome and defaults at least).
# Run at most n = num_imputations commands in one job
# Reasoning:
#   The human genome is broken in approx. 2000 chunks.
#   GLIMPSE_phase is run for each block.
#   Unfortunately, snakemake might take a long time to infer the DAG with so many jobs.
#   2000 jobs is still fine for 1 individual, but if imputating more individuals it might be worth
#   to group a few GLIMPSE_phase commands in a single job, as they are usually fast (1-2 minutes each)
num_imputations = config.get("imputation_GLIMPSE", {}).get("num_imputations", 1)

# gp = 1 means that all imputed genotypes are kept; otherwise genotypes with a genotype probability < gp will be filtered out
gp = config.get("imputation_GLIMPSE", {}).get("GP_filter", 0)

# Imputation will be run on the chromosomes present in the panel
# The version of the panel and chromosome names should match those of
# the reference genome to which individual samples were mapped.
chromosomes = [str(chr) for chr in range(1, 3)]

# Either a file listing BAM file or individual paths should be specified
if eval(config.get("imputation_GLIMPSE", {}).get("bam_list", "None")):
    bam_list = config["imputation_GLIMPSE"]["bam_list"]
    with open(bam_list, "r") as file:
        samples_files = {
            os.path.basename(line)
            .replace(".bam\n", ""): os.path(line)
            .replace(".bam\n", "")
            for line in file.readlines()
        }
else:
    samples_files = config.get("imputation_GLIMPSE", {}).get("paths", None)

samples = list(samples_files.keys())

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


wildcard_constraints:
    gp="|".join([str(value) for value in gp]),
    chr="|".join([str(chr) for chr in chromosomes]),
    sm="|".join([sm for sm in samples]),


rule all:
    input:
        phased_bcf=[
            f"GLIMPSE_phased/{sm}.GP{GP}.phased.{ext}"
            for sm in samples
            for GP in gp
            for ext in ["bcf", "bcf.csi"]
        ],


# Split the genome into chunks
#    resources:
#        runtime = 4 * 60,
#        mem = 2 * 1024
checkpoint split_genome:
    input:
        vcf_ref_panel=path_panel,  # some_phased_haplotypes_chr{chr}.vcf.gz
    output:
        chunks="chunks/chunks.chr{chr}.txt",
    params:
        window_size=config["imputation_GLIMPSE"]["window_size"]
        if "window_size" in config["imputation_GLIMPSE"].keys()
        else 1000000,
        buffer_size=config["imputation_GLIMPSE"]["buffer_size"]
        if "buffer_size" in config["imputation_GLIMPSE"].keys()
        else 200000,
    shell:
        """
        GLIMPSE_chunk --input {input.vcf_ref_panel} \
        --region {wildcards.chr} \
        --window-size {params.window_size} \
        --buffer-size {params.buffer_size} \
        --output {output.chunks}
        """


rule symlink_panel:
    input:
        path_panel,
    output:
        panel_vcf="panel/chr{chr}.vcf.gz",
        panel_vcf_csi="panel/chr{chr}.vcf.gz.csi",
    shell:
        """
        ln -s {input} {output.panel_vcf}
        ln -s {input}.csi {output.panel_vcf_csi}
        """


# extract variable positions
#    resources:
#        runtime = 6 * 60
rule extract_positions:
    input:
        vcf="panel/chr{chr}.vcf.gz",
        index="panel/chr{chr}.vcf.gz.csi",
    output:
        sites="sites/{chr}.sites.vcf.gz",
        tsv="sites/{chr}.sites.tsv.gz",
    log:
        "logs/extract_pos_{chr}.log",
    shell:
        """
        bcftools view -G -m 2 -M 2 -v snps \
        {input.vcf} \
        -Oz -o {output.sites}

        bcftools index -f {output.sites}
        bcftools query -f'%CHROM\\t%POS\\t%REF,%ALT\\n' {output.sites} | \
            bgzip -c > {output.tsv}
        tabix -s1 -b2 -e2 {output.tsv}
        """


#    resources:
# memory=lambda wildcards, attempt: get_memory_alloc("mpileup_mem", attempt, 4),
# runtime=lambda wildcards, attempt: get_runtime_alloc("mpileup_time", attempt, 6)
rule bcftools_mpileup:
    input:
        ref_genome=config["ref_genome"],
        sites="sites/{chr}.sites.vcf.gz",
        bam=lambda wildcards: samples_files[wildcards.sm],
        tsv="sites/{chr}.sites.tsv.gz",
    output:
        header=temp("vcf/individual/{sm}_chr{chr}.header"),
        temp_vcf=temp("vcf/individual/{sm}_chr{chr}.temp.vcf.gz"),
        final_vcf=temp("vcf/individual/{sm}_chr{chr}.vcf.gz"),
        final_csi=temp("vcf/individual/{sm}_chr{chr}.vcf.gz.csi"),
    threads: 1
    log:
        "logs/bcftools_mpileup_{sm}_{chr}.txt",
    shell:
        """
        bcftools mpileup -f {input.ref_genome} \
            -I -E -a 'FORMAT/DP' \
            -T {input.sites} -r {wildcards.chr} \
            -Ob --threads {threads} \
            --ignore-RG {input.bam}  | \
            bcftools call -Aim -C alleles -T {input.tsv} -Oz -o {output.temp_vcf}
        echo {input.bam} {wildcards.sm} > {output.header}
        bcftools reheader -s {output.header} -o {output.final_vcf} {output.temp_vcf}
        bcftools index -f {output.final_vcf}
        """


# # 5. Impute and phase a whole chromosome
# as each job runs for ~1 minute, we can group ~10 jobs and send them as one
#    resources:
# memory=lambda wildcards, attempt: get_memory_alloc("phase_mem", attempt, 2),
# runtime=lambda wildcards, attempt: get_runtime_alloc("phase_time", attempt, 2)
rule impute_phase:
    input:
        chunks="chunks/chunks.chr{chr}.txt",
        vcf_sample="vcf/individual/{sm}_chr{chr}.vcf.gz",
        csi_sample="vcf/individual/{sm}_chr{chr}.vcf.gz.csi",
        vcf_ref="panel/chr{chr}.vcf.gz",
        csi_ref="panel/chr{chr}.vcf.gz.csi",
        map=config["imputation_GLIMPSE"]["path_map"],
    output:
        txt=temp("GLIMPSE_imputed/done_{sm}_chr{chr}_group{g}.txt"),
    params:
        imputed_bcf=lambda wildcards: f"GLIMPSE_imputed/{wildcards.sm}_chr{wildcards.chr}.${{n}}.bcf",
        start_end=lambda wildcards: group_chunks(
            wildcards, num_imputations=num_imputations, start_end=True
        ),
    log:
        "logs/impute_phase_{sm}_chr{chr}_group{g}.txt",
    threads: 1
    shell:
        """
        for n in $(seq {params.start_end[0]} {params.start_end[1]})
        do
            LINE=$( head -n ${{n}} {input.chunks} |tail -n1 )
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)

            GLIMPSE_phase \
            --input {input.vcf_sample} \
                --reference {input.vcf_ref} \
                --map {input.map} \
                --input-region ${{IRG}} \
                --output-region ${{ORG}} \
                --output {params.imputed_bcf} \
                --thread {threads}

            bcftools index -f {params.imputed_bcf}
        done

        touch {output.txt}
        """


# #  # 6. Ligate multiple chunks together
rule ligate_chunks:
    input:
        txt=lambda wildcards: [
            f"GLIMPSE_imputed/done_{wildcards.sm}_chr{wildcards.chr}_group{g}.txt"
            for g in range(1, group_chunks(wildcards, num_imputations))
        ],
    output:
        ligated_bcf=temp("GLIMPSE_ligated/{sm}_chr{chr}.merged.bcf"),
        ligated_csi=temp("GLIMPSE_ligated/{sm}_chr{chr}.merged.bcf.csi"),
        list_files=temp("GLIMPSE_ligated/list_{sm}_chr{chr}.txt"),
    log:
        "logs/ligate_chunks_{sm}_chr{chr}.txt",
    params:
        phased=lambda wildcards: [
            f"GLIMPSE_imputed/{wildcards.sm}_chr{wildcards.chr}.{n}.bcf"
            for n in range(1, get_num_chunks(wildcards))
        ],
        phased_csi=lambda wildcards: [
            f"GLIMPSE_imputed/{wildcards.sm}_chr{wildcards.chr}.{n}.bcf.csi"
            for n in range(1, get_num_chunks(wildcards))
        ],
    shell:
        """
        for file in {params.phased} ; do echo $file ; done > {output.list_files}
        GLIMPSE_ligate --input {output.list_files} --output {output.ligated_bcf}
        rm {params.phased} {params.phased_csi}
        bcftools index -f {output.ligated_bcf}
        """


# #  # 6.1
# concatenate imputed genotypes per sample
rule concat:
    input:
        bcf=lambda wildcards: expand(
            "GLIMPSE_ligated/{sm}_chr{chr}.merged.bcf",
            chr=chromosomes,
            sm=wildcards.sm,
        ),
        index=lambda wildcards: expand(
            "GLIMPSE_ligated/{sm}_chr{chr}.merged.bcf.csi",
            chr=chromosomes,
            sm=wildcards.sm,
        ),
    output:
        bcf=temp("GLIMPSE_ligated/{sm}.bcf"),
        csi=temp("GLIMPSE_ligated/{sm}.bcf.csi"),
    log:
        "logs/concat_{sm}.txt",
    shell:
        """
        bcftools concat {input.bcf} -Oz -o {output.bcf}
        bcftools index -f {output.bcf}
        """


# # # 6.2 Filter by GP
rule filter_gp_split_sm:
    input:
        bcf="GLIMPSE_ligated/{sm}.bcf",
        index="GLIMPSE_ligated/{sm}.bcf.csi",
    output:
        bcf=temp("unphased/{sm}.GP{gp}.bcf"),
        csi=temp("unphased/{sm}.GP{gp}.bcf.csi"),
    log:
        "logs/filterGP_{sm}_GP{gp}.txt",
    shell:
        """
        bcftools view -i'FORMAT/GP[0:*]>{wildcards.gp}'  -Ob -o {output.bcf} {input.bcf}
        bcftools index -f {output.bcf}
        """


# #  # 7. Sample haplotypes
rule sample_haplotypes:
    input:
        ligated_bcf="unphased/{sm}.GP{gp}.bcf",
        ligated_csi="unphased/{sm}.GP{gp}.bcf.csi",
    output:
        phased_bcf="GLIMPSE_phased/{sm}.GP{gp}.phased.bcf",
        phased_csi="GLIMPSE_phased/{sm}.GP{gp}.phased.bcf.csi",
    log:
        "logs/sample_haplotype_{sm}.GP{gp}.txt",
    shell:
        """
        GLIMPSE_sample --input {input.ligated_bcf} \
        --solve \
        --output {output.phased_bcf}
        bcftools index -f {output.phased_bcf}
        """
