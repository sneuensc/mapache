# -----------------------------------------------------------------------------#
import pandas as pd


localrules:
    samtools_idxstats,
    plot_summary_statistics,
    merge_DoC_chr,
    assign_sex,
    merge_stats_per_fastq,
    merge_stats_per_lb,
    merge_stats_per_sm,
    merge_stats_by_level_and_genome,
    merge_stats_all_genomes,
    DoC_chr_SM,
    merge_DoC_chr,


# -----------------------------------------------------------------------------#
## get sparse stats stats


rule fastqc_for_mapping:
    """
    Quality control of fastq file by fastqc (SE or R1)
    """
    input:
        lambda wildcards: inputs_fastqc(wildcards, run_adapter_removal=run_adapter_removal),
    output:
        html="results/04_stats/01_sparse_stats/01_fastq/{folder}/{SM}/{LB}/{ID}_fastqc.html",
        zip="results/04_stats/01_sparse_stats/01_fastq/{folder}/{SM}/{LB}/{ID}_fastqc.zip",
    log:
        "results/04_stats/01_sparse_stats/01_fastq/{folder}/{SM}/{LB}/{ID}_fastqc.log",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("fastqc_mem", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("fastqc_time", attempt, 1),
    conda:
        "../envs/fastqc.yaml"
    envmodules:
        module_fastqc,
    message:
        "--- FASTQC {output.html}"
    shell:
        """
        ## symlink file to remove '_R1' of paired reads
        html={output.html}
        symlink=${{html%%_fastqc.html}}.fastq.gz
        ln -srf {input[0]} $symlink

        ## run fastqc
        fastqc --quiet --outdir $(dirname $html) $symlink 2> {log};

        ## remove symlink again
        rm $symlink
        """


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


rule samtools_index:
    input:
        "{file}.bam",
    output:
        "{file}.bam.bai",
    log:
        "{file}.bam.bai.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS INDEX of {input}"
    shell:
        """
        samtools index {input}
        """


rule samtools_idxstats:
    input:
        bam="results/{dir}/{file}.bam",
        bai="results/{dir}/{file}.bam.bai",
    output:
        idxstats="results/04_stats/01_sparse_stats/{dir}/{file}.idxstats",
    log:
        "results/04_stats/01_sparse_stats/{dir}/{file}.idxstats.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS IDXSTATS of {input.bam}"
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats}
        """


rule assign_sex:
    input:
        idxstats="results/04_stats/01_sparse_stats/{file}.{GENOME}.idxstats",
    output:
        sex="results/04_stats/01_sparse_stats/{file}.{GENOME}.sex",
    params:
        run_sex=str2bool(
            lambda wildcards: recursive_get(config, [ 
                ["genome", {}],
                [wildcards.GENOME, {}],
                ["sex_inference", {}],
                ["run", False]
            ]
            )
        ),
        #sex_params=get_sex_params,
        sex_params=lambda wildcards: " ".join(
            [
                f"--{key}='{value}'"
                for key, value in recursive_get(config, [ 
                    ["genome", {}], 
                    [wildcards.GENOME, {}],
                    ["sex_inference", {}],
                    ["params", {}]
                ]
                ).items()
            ]
        ),
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
        if [ {params.run_sex} == False ] 
        then 
            echo "Sex" > {output.sex}
            echo "Not requested in config" >> {output.sex}
        else
            Rscript workflow/scripts/assign_sex.R \
                --idxstats={input.idxstats} \
                --out={output.sex} \
                {params.sex_params}
        fi
        """


# -----------------------------------------------------------------------------#
## merging individual stats
# path_multiqc_orig = "results/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/multiqc_fastqc_data/multiqc_fastqc.txt"  # raw sequenced reads
# path_multiqc_trim = "results/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_files_trim/multiqc_fastqc_data/multiqc_fastqc.txt" # raw trimmed reads
# path_flagstat_mapped_highQ = "results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/ind1/lib1_lb/lib1_R1_002_fq.hg19_flagstat.txt"       # mapped and high-qual reads
# path_length_mapped_highQ = "results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/ind1/lib1_lb/lib1_R1_002_fq.hg19.length"


rule merge_stats_per_fastq:
    input:
    # adapterremoval settings:
    # "{folder}/01_trimmed/01_files_trim/{SM}/{LB}/{ID}.settings"
    # "{folder}/01_trimmed/01_files_trim/{SM}/{LB}/{ID}.settings",
    # "{folder}/01_trimmed/01_files_trim_collapsed/{SM}/{LB}/{ID}.settings"
        settings_stats="results/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_files_trim/{SM}/{LB}/{ID}.settings",
        fastqc_orig="results/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/{SM}/{LB}/{ID}_fastqc.zip",  # raw sequenced reads
        fastqc_trim="results/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_files_trim/{SM}/{LB}/{ID}_fastqc.zip" if run_adapter_removal else "Not trimmed",  # raw trimmed reads
        flagstat_mapped_highQ="results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/{SM}/{LB}/{ID}.{GENOME}_flagstat.txt",  # mapped and high-qual reads
        length_fastq_mapped_highQ="results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/{SM}/{LB}/{ID}.{GENOME}.length",
    output:
        temp(
            #"results/04_stats/02_separate_tables/{GENOME}/{SM}/{LB}/{ID}/fastq_stats.csv"
            "{folder}/{GENOME}/{SM}/{LB}/{ID}/fastq_stats.csv"
        ),
    log:
        #"results/04_stats/02_separate_tables/{GENOME}/{SM}/{LB}/{ID}/fastq_stats.log",
        "{folder}/{GENOME}/{SM}/{LB}/{ID}/fastq_stats.log",
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
            --path_fastqc_orig={input.fastqc_orig} \
            --path_fastqc_trim={input.fastqc_trim} \
            --path_flagstat_mapped_highQ={input.flagstat_mapped_highQ} \
            --path_length_mapped_highQ={input.length_fastq_mapped_highQ}
        """


rule merge_stats_per_lb:
    input:
        fastq_stats=lambda wildcards: [
            f"results/04_stats/02_separate_tables/{wildcards.GENOME}/{wildcards.SM}/{wildcards.LB}/{ID}/fastq_stats.csv"
            for ID in samples[wildcards.SM][wildcards.LB]
        ],
        flagstat_raw="results/04_stats/01_sparse_stats/02_library/00_merged_fastq/01_bam/{SM}/{LB}.{GENOME}_flagstat.txt",
        flagstat_unique="results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}_flagstat.txt",
        length_unique="results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}.length",
        #genomecov_unique="results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}.genomecov",
        idxstats_unique="results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}.idxstats",
        sex_unique="results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}.sex",
    output:
        temp("results/04_stats/02_separate_tables/{GENOME}/{SM}/{LB}/library_stats.csv"),
    params:
        chrs_selected=lambda wildcards: recursive_get(config, [ 
            ["genome", {}],
            [wildcards.GENOME, {}], 
            ["depth_chromosomes",  "not requested"]
        ]
        ),
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

        if [ {params.chrs_selected} == "not requested" ] 
        then
            chrsSelected=""
        else
            chrsSelected="--chrs_selected={params.chrs_selected}"
        fi

        Rscript workflow/scripts/merge_stats_per_LB.R \
            --LB={wildcards.LB} \
            --SM={wildcards.SM} \
            --genome={wildcards.GENOME} \
            --output_file={output} \
            --path_list_stats_fastq=${{list_fastq_stats}} \
            --path_flagstat_raw={input.flagstat_raw} \
            --path_flagstat_unique={input.flagstat_unique} \
            --path_length_unique={input.length_unique} \
            --path_idxstats_unique={input.idxstats_unique} \
            --path_sex_unique={input.sex_unique} \
            $chrsSelected
        """


rule merge_stats_per_sm:
    input:
        lb_stats=lambda wildcards: [
            f"results/04_stats/02_separate_tables/{wildcards.GENOME}/{wildcards.SM}/{LB}/library_stats.csv"
            for LB in samples[wildcards.SM]
        ],
        flagstat_unique="results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{GENOME}_flagstat.txt",
        length_unique="results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{GENOME}.length",
        idxstats_unique="results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{GENOME}.idxstats",
        sex_unique="results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{GENOME}.sex",
    output:
        temp("results/04_stats/02_separate_tables/{GENOME}/{SM}/sample_stats.csv"),
    params:
        chrs_selected=lambda wildcards: recursive_get(config, [ 
            ["genome", {}], 
            [wildcards.GENOME, {}], 
            ["depth_chromosomes", "not requested"]
        ]
        ),
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

        if [ {params.chrs_selected} == "not requested" ] 
        then
            chrsSelected=""
        else
            chrsSelected="--chrs_selected={params.chrs_selected}"
        fi

        Rscript workflow/scripts/merge_stats_per_SM.R \
            --SM={wildcards.SM} \
            --genome={wildcards.GENOME} \
            --output_file={output} \
            --path_list_stats_fastq=${{list_lb_stats}} \
            --path_flagstat_unique={input.flagstat_unique} \
            --path_length_unique={input.length_unique} \
            --path_idxstats_unique={input.idxstats_unique} \
            --path_sex_unique={input.sex_unique} \
            $chrsSelected
        """


# -----------------------------------------------------------------------------#
# merge all the stats in the same level


# Here level can be: SM, LB, FASTQ
rule merge_stats_by_level_and_genome:
    input:
        paths=path_stats_by_level,
    output:
        temp(
            #"results/04_stats/03_summary/{level}_stats.{GENOME}.csv"
            "{folder}/{level}_stats.{GENOME}.csv"
            ),
    log:
        #"results/04_stats/03_summary/{level}_stats.{GENOME}.log",
        "{folder}/{level}_stats.{GENOME}.log",
    message:
        "--- MERGE STATS by {wildcards.level}"
    run:
        import pandas as pd

        df_list = [pd.read_csv(file) for file in input]
        df = pd.concat(df_list)
        df.to_csv(str(output), index=False)


# Here level can be: SM, LB, FASTQ
rule merge_stats_all_genomes:
    input:
        lambda wildcards: expand(
            "results/04_stats/03_summary/{level}_stats.{GENOME}.csv",
            level=wildcards.level,
            GENOME=genome,
        ),
    output:
        report(
            "results/04_stats/03_summary/{level}_stats.csv",
            category="Mapping statistics",
            subcategory="Tables",
        ),
    log:
        "results/04_stats/03_summary/{level}_stats.log",
    message:
        "--- MERGE STATS by {wildcards.level}"
    run:
        import pandas as pd
        import re


        def get_csv(file):
            my_csv = pd.read_csv(file)
            genome = re.sub(".*_stats.", "", file).replace(".csv", "")

            my_csv["genome"] = genome
            return my_csv


        df_list = [get_csv(file) for file in input]
        df = pd.concat(df_list)
        df.to_csv(str(output), index=False)


##########################################################################################
# read depth by chromosome
rule DoC_chr_SM:
    input:
        "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{genome}.genomecov",
    output:
        "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{genome}_DoC_chrs.csv",
    params:
        SM="{SM}",
    log:
        "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{genome}_DoC_chrs.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    shell:
        """
        Rscript workflow/scripts/depth_by_chr.R \
            --path_genomecov={input} \
            --SM={params.SM} \
            --output_file={output}
        """


rule merge_DoC_chr:
    input:
        lambda wildcards: expand(
            "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{genome}_DoC_chrs.csv",
            genome=wildcards.genome,
            SM=samples,
        ),
    output:
        "results/04_stats/03_summary/DoC_by_chrs.{genome}.csv",
    log:
        "results/04_stats/03_summary/DoC_by_chrs.{genome}.log",
    run:
        import pandas as pd

        df_list = [pd.read_csv(file) for file in input]
        # this should work even if the columns are not in the same order
        # across tables
        df = pd.concat(df_list)
        df.to_csv(str(output), index=False)


##########################################################################################
#
# bamdamage


rule bamdamage:
    """
    Run bamdamage to quantify the deamination pattern
    """
    input:
        ref="results/00_reference/{id_genome}/{id_genome}.fasta",
        bam=lambda wildcards: get_mapDamage_bam(wildcards),
        bai=lambda wildcards: get_mapDamage_bam(wildcards, index=True),
    output:
        damage_pdf="results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}.dam.pdf",
        length_pdf="results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}.length.pdf",
        length_table=report(
            "results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}.length.csv",
            category="Read length",
            subcategory="Tables",
        ),
        dam_5prime_table=report(
            "results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}.dam_5prime.csv",
            category="Damage pattern",
            subcategory="Tables",
        ),
        dam_3prime_table=report(
            "results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}.dam_3prime.csv",
            category="Damage pattern",
            subcategory="Tables",
        ),
    log:
        "results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}_bamdamage.log",
    threads: 1
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapdamage_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc(
            "bamdamage_time", attempt, 24
        ),
    message:
        "--- RUN BAMDAMAGE {input.bam}"
    params:
        prefix="results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}/{id_library}.{id_genome}",
        bamdamage_params=recursive_get(config, [ ["bamdamage_params", ""] ]),
        fraction=recursive_get(config, [ ["bamdamage_fraction", 0] ]),
    log:
<<<<<<< HEAD
<<<<<<< HEAD
        "results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}_bamdamage.log",
=======
>>>>>>> Added a log to all rules.
=======
        "results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}_bamdamage.log",
>>>>>>> Corrected log statement.
    conda:
        "../envs/bamdamage.yaml"
    envmodules:
        module_r,
        module_samtools,
    shell:
        """
        nb=$(awk '{{sum += $3}} END {{print sum}}' {input.bai}); 
        nth_line=1; 

        # this part is to take a fraction of the total reads
        # to run bamdamage
        if [ {params.fraction} -eq 0 ]; then
            nth_line=1; 
        elif [ {params.fraction} -lt 1 ]; then
            nth_line=$(( {params.fraction} * $nb )); 
        elif [ {params.fraction} -lt "$nb" ]; then
           nth_line=$(( $nb / {params.fraction} )); 
        fi;

        workflow/scripts/bamdamage {params.bamdamage_params} \
            --nth_read $nth_line --output {output.damage_pdf} \
            --output_length {output.length_pdf} {input.bam} 2> {log};
        """


rule plot_bamdamage:
    """
    Run bamdamage to quantify the deamination pattern
    """
    input:
        length_table="results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}.length.csv",
        dam_5prime_table="results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}.dam_5prime.csv",
        dam_3prime_table="results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}.dam_3prime.csv",
    output:
        length=report(
            "results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}.length.svg",
            category="Read length",
            subcategory="Plots",
        ),
        damage=report(
            "results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}.dam.svg",
            category="Damage pattern",
            subcategory="Plots",
        ),
    message:
        "--- PLOT DAMAGE"
    log:
        "results/04_stats/01_sparse_stats/02_library/04_bamdamage/{id_sample}/{id_library}.{id_genome}_plot.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    shell:
        """
        Rscript workflow/scripts/plot_bamdamage.R \
            --length={input.length_table} \
            --five_prime={input.dam_5prime_table} \
            --three_prime={input.dam_3prime_table} \
            --sample={wildcards.id_sample} \
            --library={wildcards.id_library} \
            --genome={wildcards.id_genome} \
            --length_svg={output.length} \
            --damage_svg={output.damage}
        """


##########################################################################################
# plots
# -----------------------------------------------------------------------------#
# plotting


rule plot_depth_statistics:
    input:
        sample_depth="{folder}/depth_stats_{GENOME}.csv",
    output:
        plot_5_AvgReadDepth=report(
            "{folder}/5_AvgReadDepth.{GENOME}.svg",
            caption="../report/5_AvgReadDepth.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
        plot_6_AvgReadDepth_MT=report(
            "{folder}/6_AvgReadDepth_MT.{GENOME}.svg",
            caption="../report/6_AvgReadDepth_MT.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
        plot_7_Sex=report(
            "{folder}/7_Sex.{GENOME}.svg",
            caption="../report/7_Sex.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
    threads: 1
    log:
        "{folder}/depth_stats_plot.{GENOME}.csv.log",
    message:
        "--- PLOT DEPTH"
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- PLOT DEPTH STATISTICS OF {input}"
    script:
        "../scripts/plot_depth.R"


rule plot_summary_statistics:
    """
    Plot summary statistics
    """
    input:
        sample_stats="results/04_stats/03_summary/SM_stats.csv",
    output:
        plot_1_nb_reads=report(
            "results/04_stats/04_plots/1_nb_reads.svg",
            caption="../report/1_nb_reads.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
        plot_2_mapped=report(
            "results/04_stats/04_plots/2_mapped.svg",
            caption="../report/2_mapped.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
        plot_3_endogenous=report(
            "results/04_stats/04_plots/3_endogenous.svg",
            caption="../report/3_endogenous.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
        plot_4_duplication=report(
            "results/04_stats/04_plots/4_duplication.svg",
            caption="../report/4_duplication.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
        plot_5_AvgReadDepth=report(
            "results/04_stats/04_plots/5_AvgReadDepth.svg",
            caption="../report/5_AvgReadDepth.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
    log:
        "results/04_stats/04_plots/plot_summary_statistics.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- PLOT SUMMARY STATISTICS"
    params:
        samples=recursive_get(config, [ ["sample_file", ""] ]),
        x_axis=recursive_get(config, [ ["stats",{}], ["plots",{}], ["x_axis", "sample"] ]),
        split_plot=recursive_get(config, [ ["stats", {}], ["plots", {}], ["split_plot", "F"] ]),
        n_col=recursive_get(config, [ ["stats", {}], ["plots", {}], ["n_col", 1] ]),
        n_row=recursive_get(config, [ ["stats", {}], ["plots", {}], ["n_row", 1] ]),
    shell:
        """
        Rscript workflow/scripts/plot_stats.R \
            --samples={params.samples} \
            --SM={input.sample_stats}  \
            --out_1_reads={output.plot_1_nb_reads} \
            --out_2_mapped={output.plot_2_mapped} \
            --out_3_endogenous={output.plot_3_endogenous} \
            --out_4_duplication={output.plot_4_duplication} \
            --out_5_AvgReadDepth={output.plot_5_AvgReadDepth} \
            --x_axis={params.x_axis} \
            --split_plot={params.split_plot} \
            --n_col={params.n_col} \
            --n_row={params.n_row}
        """
