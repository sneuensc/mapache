# -----------------------------------------------------------------------------#
import pandas as pd


# -----------------------------------------------------------------------------#
## get sparse stats stats


rule fastqc:
    """
    Quality control of fastq file by fastqc (SE or R1)
    """
    input:
        "results/{file}.fastq.gz",
    output:
        html="results/04_stats/01_sparse_stats/{file}_fastqc.html",
        zip="results/04_stats/01_sparse_stats/{file}_fastqc.zip",
    log:
        "results/04_stats/01_sparse_stats/{file}_fastqc.log",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("fastqc_mem", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("fastqc_time", attempt, 1),
    conda:
        "../envs/fastqc.yaml"
    envmodules:
        module_fastqc,
    message:
        "--- FASTQC {input}"
    shell:
        "fastqc --quiet --outdir $(dirname {output.html}) {input} 2> {log}"


rule samtools_flagstat:
    """
    Compute samtools flagstat on bam file
    """
    input:
        bam="results/{file}.bam",
    output:
        "results/04_stats/01_sparse_stats/{file}_flagstat.txt",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc(
            "samtools_flagstat_mem", attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "samtools_flagstat_time", attempt, 1
        ),
    log:
        "results/04_stats/01_sparse_stats/{file}_flagstat.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS FLAGSTAT {input}"
    shell:
        "samtools flagstat --threads {threads} {input} > {output} 2> {log};"


rule multiqc_fastqc:
    """
    Merging fastqc output with multiqc
    """
    input:
        lambda wildcards: [
            f"results/04_stats/01_sparse_stats/{wildcards.folder}/{SM}/{LB}/{ID}_fastqc.zip"
            for SM in samples
            for LB in samples[SM]
            for ID in samples[SM][LB]
        ],
    output:
        html="results/04_stats/01_sparse_stats/{folder}/multiqc_fastqc.html",
        txt="results/04_stats/01_sparse_stats/{folder}/multiqc_fastqc_data/multiqc_fastqc.txt",
        html=report(
            "results/{dir_stats}/{folder}/multiqc_fastqc.html",
            category=" Quality control",
        ),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("multiqc_mem", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("multiqc_time", attempt, 1),
    log:
        "results/04_stats/01_sparse_stats/{folder}/multiqc_fastqc.log",
    conda:
        "../envs/multiqc.yaml"
    envmodules:
        module_multiqc,
    message:
        "--- MULTIQC fastqc: {input}"
    shell:
        """
        multiqc -n $(basename {output.html}) -f -d -o $(dirname {output.html}) {input}  2> {log}
        """


rule bedtools_genomecov:
    input:
        bam="results/{dir}/{file}.bam",
    output:
        genomecov="results/04_stats/01_sparse_stats/{dir}/{file}.genomecov",
    log:
        "results/04_stats/01_sparse_stats/{dir}/{file}.genomecov.log",
    conda:
        "../envs/bedtools.yaml"
    envmodules:
        module_bedtools,
    message:
        "--- BEDTOOLS GENOMECOV of {input.bam}"
    shell:
        """
        bedtools genomecov -ibam {input.bam} > {output.genomecov}
        """


rule read_length:
    input:
        bam="results/{file}.bam",
    output:
        length="results/04_stats/01_sparse_stats/{file}.length",
    log:
        "results/04_stats/01_sparse_stats/{file}.length.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- READ LENGTH of {input}"
    shell:
        """
        samtools view {input.bam} | workflow/scripts/read_length.pl -o {output.length}
        """


def is_quick(file_name, dict):
    if "quick" in dict.keys() and dict["quick"]:
        file_name.replace(".bam", ".downsampled.bam")
    return file_name


# sex_params = config["genome"]["sex"]["sex_params"] if "sex_params" in config["genome"]["sex"].keys() else {}


def get_sex_params(name):
    sex_params = get_param2("genome", name, {})
    x = " ".join(
        [f"--{key}='{eval_to_csv(sex_params[key])}'" for key in sex_params.keys()]
    )
    return x


rule assign_sex:
    input:
        genomecov="results/04_stats/01_sparse_stats/{file}.{GENOME}.genomecov",
        test="results/00_reference/{GENOME}/{GENOME}.ok",
    output:
        sex="results/04_stats/01_sparse_stats/{file}.{GENOME}.sex",
    params:
        sex_params=lambda wildcards: get_sex_params(wildcards.GENOME),
    log:
        "results/04_stats/01_sparse_stats/{file}.{GENOME}.sex.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- SEX ASSIGNEMENT {input}"
    shell:
        """
        Rscript workflow/scripts/sex_assignation.r \
            --genomecov={input.genomecov} \
            --out={output.sex} \
            {params.sex_params}
        """


# -----------------------------------------------------------------------------#
## merging individual stats
# path_multiqc_orig = "results/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/multiqc_fastqc_data/multiqc_fastqc.txt"  # raw sequenced reads
# path_multiqc_trim = "results/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_files_trim/multiqc_fastqc_data/multiqc_fastqc.txt" # raw trimmed reads
# path_flagstat_mapped_highQ = "results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/ind1/lib1_lb/lib1_R1_002_fq.hg19_flagstat.txt"       # mapped and high-qual reads
# path_length_mapped_highQ = "results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/ind1/lib1_lb/lib1_R1_002_fq.hg19.length"


rule merge_stats_per_fastq:
    input:
        multiqc_orig="results/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/multiqc_fastqc_data/multiqc_fastqc.txt",  # raw sequenced reads
        multiqc_trim="results/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_files_trim/multiqc_fastqc_data/multiqc_fastqc.txt",  # raw trimmed reads
        flagstat_mapped_highQ="results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/{SM}/{LB}/{ID}.{GENOME}_flagstat.txt",  # mapped and high-qual reads
        length_fastq_mapped_highQ="results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/{SM}/{LB}/{ID}.{GENOME}.length",
    output:
        "results/04_stats/02_separate_tables/{GENOME}/{SM}/{LB}/{ID}/fastq_stats.csv",
    log:
        "results/04_stats/02_separate_tables/{GENOME}/{SM}/{LB}/{ID}/fastq_stats.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- MERGE FASTQ LEVEL STATS"
    shell:
        """
        Rscript workflow/scripts/merge_stats_per_fastq.R \
            --ID={wildcards.ID} \
            --LB={wildcards.LB} \
            --SM={wildcards.SM} \
            --genome={wildcards.GENOME} \
            --output_file={output} \
            --path_multiqc_orig={input.multiqc_orig} \
            --path_multiqc_trim={input.multiqc_trim} \
            --path_flagstat_mapped_highQ={input.flagstat_mapped_highQ} \
            --path_length_mapped_highQ={input.length_fastq_mapped_highQ}
        """


## get all chromosome names given in the stats part and check if the name is given in the genome part. In this case overtake it
## Example config:
## genome:
##    GRCh38:
##        maleChr: chrY
##        femaleChr: chrX
##        mtChr: chrM
##        autosomeChr: '[f"chr{x}" for x in range(1,5)]'
## stats:
##    sample:
##        depth_chromosomes: "femaleChr, maleChr, mtChr, autosomeChr"
## ==> 'chrY,chrX,chr1,chr2,chr3,chr4,chrM'


def get_chrom(wildcards):
    chr = eval_list(
        "".join(get_param3("stats", "sample", "depth_chromosomes", "").split()).split(
            ","
        )
    )
    gen = get_param2("genome", wildcards.GENOME, {})
    chr_uniq = list(set(chr) - set(list(gen)))
    chr_def = list(set(chr).intersection(set(gen)))
    return ",".join(eval_list(chr_uniq) + [eval_if_possible(gen[c]) for c in chr_def])


rule merge_stats_per_lb:
    input:
        fastq_stats=lambda wildcards: [
            f"results/04_stats/02_separate_tables/{wildcards.GENOME}/{wildcards.SM}/{wildcards.LB}/{ID}/fastq_stats.csv"
            for ID in samples[wildcards.SM][wildcards.LB]
        ],
        flagstat_raw="results/04_stats/01_sparse_stats/02_library/00_merged_fastq/01_bam/{SM}/{LB}.{GENOME}_flagstat.txt",
        flagstat_unique="results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}_flagstat.txt",
        length_unique="results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}.length",
        genomecov_unique="results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}.genomecov",
        sex_unique="results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}.sex",
    output:
        "results/04_stats/02_separate_tables/{GENOME}/{SM}/{LB}/library_stats.csv",
    params:
        chrs_selected=get_chrom,
    log:
        "results/04_stats/02_separate_tables/{GENOME}/{SM}/{LB}/library_stats.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- MERGE LIBRARY LEVEL STATS"
    shell:
        """
        list_fastq_stats=$(echo {input.fastq_stats} |sed 's/ /,/g');
        Rscript workflow/scripts/merge_stats_per_LB.R \
            --LB={wildcards.LB} \
            --SM={wildcards.SM} \
            --genome={wildcards.GENOME} \
            --output_file={output} \
            --path_list_stats_fastq=${{list_fastq_stats}} \
            --path_flagstat_raw={input.flagstat_raw} \
            --path_flagstat_unique={input.flagstat_unique} \
            --path_length_unique={input.length_unique} \
            --path_genomecov_unique={input.genomecov_unique} \
            --path_sex_unique={input.sex_unique} \
            --chrs_selected={params.chrs_selected}
        """


rule merge_stats_per_sm:
    input:
        lb_stats=lambda wildcards: [
            f"results/04_stats/02_separate_tables/{wildcards.GENOME}/{wildcards.SM}/{LB}/library_stats.csv"
            for LB in samples[wildcards.SM]
        ],
        flagstat_unique="results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{GENOME}_flagstat.txt",
        length_unique="results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{GENOME}.length",
        genomecov_unique="results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{GENOME}.genomecov",
        sex_unique="results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{GENOME}.sex",
    output:
        "results/04_stats/02_separate_tables/{GENOME}/{SM}/sample_stats.csv",
    params:
        chrs_selected=get_chrom,
    log:
        "results/04_stats/02_separate_tables/{GENOME}/{SM}/sample_stats.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- MERGE SAMPLE LEVEL STATS of  of {wildcards.SM} / {wildcards.GENOME}"
    shell:
        """
        list_lb_stats=$(echo {input.lb_stats} |sed 's/ /,/g');
        Rscript workflow/scripts/merge_stats_per_SM.R \
            --SM={wildcards.SM} \
            --genome={wildcards.GENOME} \
            --output_file={output} \
            --path_list_stats_fastq=${{list_lb_stats}} \
            --path_flagstat_unique={input.flagstat_unique} \
            --path_length_unique={input.length_unique} \
            --path_genomecov_unique={input.genomecov_unique} \
            --path_sex_unique={input.sex_unique} \
            --chrs_selected={params.chrs_selected}
        """


# -----------------------------------------------------------------------------#
# merge all the stats in the same level


def path_stats_by_level(level):
    if level == "FASTQ":
        paths = [
            f"results/04_stats/02_separate_tables/{gen}/{SM}/{LB}/{ID}/fastq_stats.csv"
            for gen in genome
            for SM in samples
            for LB in samples[SM]
            for ID in samples[SM][LB]
        ]
    elif level == "LB":
        paths = [
            f"results/04_stats/02_separate_tables/{gen}/{SM}/{LB}/library_stats.csv"
            for gen in genome
            for SM in samples
            for LB in samples[SM]
        ]
    elif level == "SM":
        paths = [
            f"results/04_stats/02_separate_tables/{gen}/{SM}/sample_stats.csv"
            for gen in genome
            for SM in samples
        ]
    return paths


rule merge_stats_by_level:
    input:
        paths=lambda wildcards: path_stats_by_level(wildcards.level),
    output:
        "results/04_stats/03_final_tables/{level}.csv",
    log:
        "results/04_stats/03_final_tables/{level}.log",
    message:
        "--- MERGE STATS by {wildcards.level}"
    run:
        import pandas as pd

        df_list = [pd.read_csv(file) for file in input]
        df = pd.concat(df_list)
        df.to_csv(str(output), index=False)


# -----------------------------------------------------------------------------#
# plotting


rule plot_summary_statistics:
    """
    Plot summary statistics
    """
    input:
        fastq_stats="results/04_stats/03_final_tables/FASTQ.csv",
        library_stats="results/04_stats/03_final_tables/LB.csv",
        sample_stats="results/04_stats/03_final_tables/SM.csv",
    output:
        plot_1_nb_reads=report(
            "results/04_stats/03_final_tables/04_final_plots/1_nb_reads.png",
            caption="../report/1_nb_reads.rst",
            category="Mapping statistics plots",
        ),
        plot_2_mapped=report(
            "results/04_stats/03_final_tables/04_final_plots/2_mapped.png",
            caption="../report/2_mapped.rst",
            category="Mapping statistics plots",
        ),
        plot_3_endogenous=report(
            "results/04_stats/03_final_tables/04_final_plots/3_endogenous.png",
            caption="../report/3_endogenous.rst",
            category="Mapping statistics plots",
        ),
        plot_4_duplication=report(
            "results/04_stats/03_final_tables/04_final_plots/4_duplication.png",
            caption="../report/4_duplication.rst",
            category="Mapping statistics plots",
        ),
    log:
        "results/04_stats/03_final_tables/04_final_plots/plot.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- PLOT SUMMARY STATISTICS"
    shell:
        """
        Rscript workflow/scripts/sex_assignation.r \
            --genomecov={input.genomecov} \
            --out={output.sex} \
            --larger_chr={params.sex_params}
        """


rule plot_depth_statistics:
    input:
        sample_depth="{folder}/depth_stats_{GENOME}.csv",
    output:
        plot_5_AvgReadDepth=report(
            "{folder}/5_AvgReadDepth.{GENOME}.svg",
            caption="../report/5_AvgReadDepth.rst",
            category="Mapping statistics plots",
        ),
        plot_6_AvgReadDepth_MT=report(
            "{folder}/6_AvgReadDepth_MT.{GENOME}.svg",
            caption="../report/6_AvgReadDepth_MT.rst",
            category="Mapping statistics plots",
        ),
        plot_7_Sex=report(
            "{folder}/7_Sex.{GENOME}.svg",
            caption="../report/7_Sex.rst",
            category="Mapping statistics plots",
        ),
    threads: 1
    log:
        "{folder}/depth_stats_plot.{GENOME}.csv.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- PLOT DEPTH STATISTICS OF {input}"
    script:
        "../scripts/plot_depth.R"
