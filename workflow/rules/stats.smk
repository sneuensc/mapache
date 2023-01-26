##########################################################################################
## all statistic rules
##########################################################################################


localrules:
    samtools_idxstats,
    plot_summary_statistics,
    merge_DoC_chr,
    assign_sex,
    assign_no_sex,
    merge_stats_per_fastq,
    merge_stats_per_lb,
    merge_stats_per_sm,
    merge_stats_by_level_and_genome,
    merge_stats_all_genomes,
    DoC_chr_SM,
    merge_DoC_chr,


# -----------------------------------------------------------------------------#
## get sparse stats stats
rule fastqc:
    """
    Quality control of fastq file by fastqc (SE or R1)
    """
    input:
        fastq=inputs_fastqc,
    output:
        html="{folder}/04_stats/01_sparse_stats/01_fastq/{type}/{sm}/{lb}/{id}_fastqc.html",
        zip="{folder}/04_stats/01_sparse_stats/01_fastq/{type}/{sm}/{lb}/{id}_fastqc.zip",
    log:
        "{folder}/04_stats/01_sparse_stats/01_fastq/{type}/{sm}/{lb}/{id}_fastqc.log",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["stats", "fastqc"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["stats", "fastqc"], attempt, 1
        ),
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
        ln -srf {input.fastq[0]} $symlink

        ## run fastqc (-t 2 is a trick to increase the memory allocation / multitreathing is per file, so still only 1CPU used)
        fastqc --quiet -t 2 --outdir $(dirname $html) $symlink 2> {log};

        ## remove symlink again
        rm $symlink
        """


rule samtools_flagstat:
    """
    Compute samtools flagstat on bam file
    """
    input:
        bam=get_bam_file,
    output:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_flagstat.txt",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["stats", "samtools_flagstat"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["stats", "samtools_flagstat"], attempt, 1
        ),
    log:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_flagstat.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS FLAGSTAT {input}"
    shell:
        "samtools flagstat --threads {threads} {input} > {output} 2> {log};"


rule samtools_stats:
    """
    Compute samtools stats on bam file
    """
    input:
        get_bam_file,
    output:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_stats.txt",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["stats", "samtools_stats"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["stats", "samtools_stats"], attempt, 1
        ),
    log:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_stats.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS FLAGSTAT {input}"
    shell:
        "samtools stats --threads {threads} {input} > {output} 2> {log};"


rule bedtools_genomecov:
    input:
        "{folder}/{dir}/{file}.bam",
    output:
        "{folder}/04_stats/01_sparse_stats/{dir}/{file}_genomecov",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["stats", "bedtools_genomecov"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["stats", "bedtools_genomecov"], attempt, 1
        ),
    log:
        "{folder}/04_stats/01_sparse_stats/{dir}/{file}_genomecov.log",
    conda:
        "../envs/bedtools.yaml"
    envmodules:
        module_bedtools,
    message:
        "--- BEDTOOLS GENOMECOV of {input.bam}"
    shell:
        """
        bedtools genomecov -ibam {input} > {output}
        """


rule read_length:
    input:
        get_bam_file,
    output:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_length.txt",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["stats", "read_length"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["stats", "read_length"], attempt, 2
        ),
    log:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_length.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    params:
        script=workflow.source_path("../scripts/read_length.pl"),
    message:
        "--- READ LENGTH of {input}"
    shell:
        """
        samtools view {input} | perl {params.script} -o {output} >> {log}
        """


rule samtools_idxstats:
    input:
        bam=get_bam_file,
        bai=lambda wildcards: bam2bai(get_bam_file(wildcards)),
    output:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_idxstats.txt",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["stats", "samtools_idxstats"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["stats", "samtools_idxstats"], attempt, 1
        ),
    log:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_idxstats.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools,
    message:
        "--- SAMTOOLS IDXSTATS of {input.bam}"
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """


rule assign_sex:
    input:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_idxstats.txt",
    output:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_sex.txt",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["stats", "assign_sex"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["stats", "assign_sex"], attempt, 1
        ),
    params:
        sex_params=lambda wildcards: " ".join(
            [
                f"--{key}='{value}'"
                for key, value in get_param(
                    ["sex_inference", wildcards.genome, "params"], {}
                ).items()
            ]
        ),
        script=workflow.source_path("../scripts/assign_sex.R"),
    log:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_sex.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- SEX ASSIGNEMENT {input}"
    shell:
        """
        Rscript {params.script} \
                --idxstats={input} \
                --out={output} \
                {params.sex_params}
        """


rule assign_no_sex:
    output:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_nosex.txt",
    log:
        "{folder}/04_stats/01_sparse_stats/{file}.{genome}_sex.log",
    message:
        "--- NO SEX ASSIGNEMENT"
    shell:
        """
        echo "Sex,Rx,CI,signif_set,reads_autosomes,reads_X" > {output}
        echo "NaN,NaN,NaN,NaN,NaN,NaN" >> {output}
        """


# -----------------------------------------------------------------------------#
## merging individual stats
rule merge_stats_per_fastq:
    input:
        fastqc_orig="{folder}/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/{sm}/{lb}/{id}_fastqc.zip",  # raw sequenced reads
        fastqc_trim=lambda wildcards: "{folder}/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_adapter_removal/{sm}/{lb}/{id}_fastqc.zip"
        if get_paramGrp(["adapterremoval", "run"], ["True", "False"], wildcards)
        else "{folder}/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/{sm}/{lb}/{id}_fastqc.zip",
        stats_mapped_highQ="{folder}/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/{sm}/{lb}/{id}.{genome}_stats.txt",  # mapped and high-qual reads
        length_fastq_mapped_highQ="{folder}/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/{sm}/{lb}/{id}.{genome}_length.txt",
    output:
        "{folder}/04_stats/02_separate_tables/{genome}/{sm}/{lb}/{id}/fastq_stats.csv",
    log:
        "{folder}/04_stats/02_separate_tables/{genome}/{sm}/{lb}/{id}/fastq_stats.log",
    params:
        script=workflow.source_path("../scripts/merge_stats_per_fastq.R"),
        script_parse_fastqc=workflow.source_path("../scripts/parse_fastqc.R"),
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- MERGE FASTQ LEVEL STATS"
    shell:
        """

        Rscript {params.script} \
            --id={wildcards.id} \
            --lb={wildcards.lb} \
            --sm={wildcards.sm} \
            --genome={wildcards.genome} \
            --output_file={output} \
            --path_fastqc_orig={input.fastqc_orig} \
            --path_fastqc_trim={input.fastqc_trim} \
            --path_stats_mapped_highQ={input.stats_mapped_highQ} \
            --path_length_mapped_highQ={input.length_fastq_mapped_highQ} \
            --script_parse_fastqc={params.script_parse_fastqc}
        """


rule merge_stats_per_lb:
    input:
        fastq_stats=lambda wildcards: expand(
            "{folder}/04_stats/02_separate_tables/{genome}/{sm}/{lb}/{id}/fastq_stats.csv",
            id=SAMPLES[wildcards.sm][wildcards.lb],
            allow_missing=True,
        ),
        stats_unique="{folder}/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{sm}/{lb}.{genome}_stats.txt",
        length_unique="{folder}/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{sm}/{lb}.{genome}_length.txt",
        idxstats_unique="{folder}/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{sm}/{lb}.{genome}_idxstats.txt",
    output:
        "{folder}/04_stats/02_separate_tables/{genome}/{sm}/{lb}/library_stats.csv",
    params:
        chrs_selected=lambda wildcards: ",".join(
            to_list(
                get_param(["depth", wildcards.genome, "chromosomes"], "not_requested")
            )
        ),
        script=workflow.source_path("../scripts/merge_stats_per_LB.R"),
    log:
        "{folder}/04_stats/02_separate_tables/{genome}/{sm}/{lb}/library_stats.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- MERGE LIBRARY LEVEL STATS"
    shell:
        """
        list_fastq_stats=$(echo {input.fastq_stats} |sed 's/ /,/g');

        if [ {params.chrs_selected} == "not_requested" ] 
        then
            chrsSelected=""
        else
            chrsSelected="--chrs_selected={params.chrs_selected}"
        fi

        Rscript {params.script} \
            --lb={wildcards.lb} \
            --sm={wildcards.sm} \
            --genome={wildcards.genome} \
            --output_file={output} \
            --path_list_stats_fastq=${{list_fastq_stats}} \
            --path_stats_unique={input.stats_unique} \
            --path_length_unique={input.length_unique} \
            --path_idxstats_unique={input.idxstats_unique} \
            $chrsSelected
        """


def get_chroms(wildcards):
    chrs_selected = (
        to_list(get_param(["depth", wildcards.genome, "chromosomes"], "not_requested")),
    )
    # print(chrs_selected)
    return chrs_selected


rule merge_stats_per_sm:
    input:
        lb_stats=get_lb_stats,
        length_unique="{folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_length.txt",
        idxstats_unique="{folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_idxstats.txt",
        sex_unique=lambda wildcards: get_sex_file(wildcards),
    output:
        "{folder}/04_stats/02_separate_tables/{genome}/{sm}/sample_stats.csv",
    params:
        chrs_selected=lambda wildcards: ",".join(
            to_list(
                get_param(["depth", wildcards.genome, "chromosomes"], "not_requested")
            )
        ),
        #chrs_selected=get_chroms,
        script=workflow.source_path("../scripts/merge_stats_per_SM.R"),
    log:
        "{folder}/04_stats/02_separate_tables/{genome}/{sm}/sample_stats.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- MERGE SAMPLE LEVEL STATS of {wildcards.sm} / {wildcards.genome}"
    shell:
        """
        list_lb_stats=$(echo {input.lb_stats} |sed 's/ /,/g');

        if [ {params.chrs_selected} == "not_requested" ] 
        then
            chrsSelected=""
        else
            chrsSelected="--chrs_selected={params.chrs_selected}"
        fi

        Rscript {params.script} \
            --sm={wildcards.sm} \
            --genome={wildcards.genome} \
            --output_file={output} \
            --path_list_stats_lb=$list_lb_stats \
            --path_length_unique={input.length_unique} \
            --path_idxstats_unique={input.idxstats_unique} \
            --path_sex_unique={input.sex_unique} \
            $chrsSelected
        """


# -----------------------------------------------------------------------------#
# merge all the stats in the same level


# Here level can be: sm, lb, FASTQ
rule merge_stats_by_level_and_genome:
    input:
        paths=path_stats_by_level,
    output:
        "{folder}/04_stats/03_summary/{level}_stats.{genome}.csv",
    log:
        "{folder}/04_stats/03_summary/{level}_stats.{genome}.log",
    message:
        "--- MERGE STATS by {wildcards.level}"
    run:
        import pandas as pd

        df_list = [pd.read_csv(file) for file in input]
        df = pd.concat(df_list)
        df.to_csv(str(output), index=False)


# Here level can be: sm, lb, FASTQ
rule merge_stats_all_genomes:
    input:
        expand(
            "{folder}/04_stats/03_summary/{level}_stats.{genome}.csv",
            genome=GENOMES,
            allow_missing=True,
        ),
    output:
        report(
            "{folder}/04_stats/03_summary/{level}_stats.csv",
            category="Mapping statistics",
            subcategory="Tables",
        ),
    log:
        "{folder}/04_stats/03_summary/{level}_stats.log",
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
        "{folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_genomecov",
    output:
        "{folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_DoC_chrs.csv",
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["stats", "DoC_chr_SM"], attempt, 2
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["stats", "DoC_chr_SM"], attempt, 1
        ),
    params:
        script=workflow.source_path("../scripts/depth_by_chr.R"),
    log:
        "{folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_DoC_chrs.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    shell:
        """
        Rscript {params.script} \
            --path_genomecov={input} \
            --sm={wildcards.sm} \
            --output_file={output}
        """


rule merge_DoC_chr:
    input:
        lambda wildcards: expand(
            "{folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_DoC_chrs.csv",
            sm=SAMPLES,
            allow_missing=True,
        ),
    output:
        "{folder}/04_stats/03_summary/DoC_by_chrs.{genome}.csv",
    log:
        "{folder}/04_stats/03_summary/DoC_by_chrs.{genome}.log",
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
        ref="{folder}/00_reference/{genome}/{genome}.fasta",
        bam=get_final_bam_LB,
        idxstats="{folder}/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{sm}/{lb}.{genome}_idxstats.txt",
    output:
        damage_pdf="{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}.dam.pdf",
        length_pdf="{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}.length.pdf",
        length_table=report(
            "{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}.length.csv",
            category="Read length",
            subcategory="Tables",
        ),
        dam_5prime_table=report(
            "{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}.dam_5prime.csv",
            category="Damage pattern",
            subcategory="Tables",
        ),
        dam_3prime_table=report(
            "{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}.dam_3prime.csv",
            category="Damage pattern",
            subcategory="Tables",
        ),
    log:
        "{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}_bamdamage.log",
    threads: 1
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc(["bamdamage"], attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc(["bamdamage"], attempt, 24),
    message:
        "--- RUN BAMDAMAGE {input.bam}"
    params:
        params=get_param(["damage", "bamdamage_params"], ""),
        fraction=get_param(["damage", "bamdamage_fraction"], 0),
        script=workflow.source_path("../scripts/bamdamage"),
    log:
        "{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}_bamdamage.log",
    conda:
        "../envs/bamdamage.yaml"
    envmodules:
        module_r,
        module_samtools,
    shell:
        """
        ## get the subsampling interval
        nb=$(awk '{{sum += $3}} END {{print sum}}' {input.idxstats}); 
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

        perl {params.script} {params.params} \
            --nth_read $nth_line --output {output.damage_pdf} \
            --output_length {output.length_pdf} {input.bam} 2>> {log};
        """


rule plot_bamdamage:
    """
    Run bamdamage to quantify the deamination pattern
    """
    input:
        length_table="{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}.length.csv",
        dam_5prime_table="{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}.dam_5prime.csv",
        dam_3prime_table="{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}.dam_3prime.csv",
    output:
        length=report(
            "{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}.length.svg",
            category="Read length",
            subcategory="Plots",
        ),
        damage=report(
            "{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}.dam.svg",
            category="Damage pattern",
            subcategory="Plots",
        ),
    message:
        "--- PLOT DAMAGE"
    params:
        script=workflow.source_path("../scripts/plot_bamdamage.R"),
        params=get_param(["stats", "bamdamage_params"], ""),
    log:
        "{folder}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}_plot.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    shell:
        """
        ## extract the plot_length
        plot_length=$(echo {params.params} | sed 's/=/ /g' | awk -F '--plot_length' '{{print $2}}'  | awk '{{print $1}}')
        if [ "$plot_length" != "" ]; then plot_length=--plot_length=$plot_length; fi

        Rscript {params.script} \
            --length={input.length_table} \
            --five_prime={input.dam_5prime_table} \
            --three_prime={input.dam_3prime_table} \
            --sample={wildcards.sm} \
            --library={wildcards.lb} \
            --genome={wildcards.genome} \
            --length_svg={output.length} \
            --damage_svg={output.damage} \
            $plot_length

        ## delete the unwanted created Rplots.pdf...
        rm -f Rplots.pdf
        """


##########################################################################################
# plots
# -----------------------------------------------------------------------------#
# plotting


rule plot_summary_statistics:
    """
    Plot summary statistics
    """
    input:
        sample_stats="{folder}/04_stats/03_summary/SM_stats.csv",
    output:
        plot_1_nb_reads=report(
            "{folder}/04_stats/04_plots/1_nb_reads.svg",
            caption="../report/1_nb_reads.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
        plot_2_mapped=report(
            "{folder}/04_stats/04_plots/2_mapped.svg",
            caption="../report/2_mapped.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
        plot_3_endogenous=report(
            "{folder}/04_stats/04_plots/3_endogenous.svg",
            caption="../report/3_endogenous.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
        plot_4_duplication=report(
            "{folder}/04_stats/04_plots/4_duplication.svg",
            caption="../report/4_duplication.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
        plot_5_AvgReadDepth=report(
            "{folder}/04_stats/04_plots/5_AvgReadDepth.svg",
            caption="../report/5_AvgReadDepth.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
        plot_6_Sex=report(
            "{folder}/04_stats/04_plots/6_Sex.svg",
            caption="../report/6_Sex.rst",
            category="Mapping statistics",
            subcategory="Plots",
        ),
    log:
        "{folder}/04_stats/04_plots/plot_summary_statistics.log",
    conda:
        "../envs/r.yaml"
    envmodules:
        module_r,
    message:
        "--- PLOT SUMMARY STATISTICS"
    params:
        script=workflow.source_path("../scripts/plot_stats.R"),
        x_axis=get_param(["stats", "plots", "x_axis"], "auto"),
        n_col=get_param(["stats", "plots", "n_col"], 1),
        width=get_param(["stats", "plots", "width"], 11),
        height=get_param(["stats", "plots", "height"], 7),
        color=get_param(["stats", "plots", "color"], "blue"),
        sex_ribbons=get_param(
            ["stats", "plots", "sex_ribbons"], 'c("XX"="red","XY"="blue")'
        ).replace("=", "?"),
        sex_thresholds=get_sex_threshold_plotting(),
        show_numbers = get_param(["stats", "plots", "show_numbers"], "10"),
    shell:
        """
        Rscript {params.script} \
            --sm={input.sample_stats}  \
            --out_1_reads={output.plot_1_nb_reads} \
            --out_2_mapped={output.plot_2_mapped} \
            --out_3_endogenous={output.plot_3_endogenous} \
            --out_4_duplication={output.plot_4_duplication} \
            --out_5_AvgReadDepth={output.plot_5_AvgReadDepth} \
            --out_6_Sex={output.plot_6_Sex} \
            --x_axis={params.x_axis} \
            --n_col={params.n_col} \
            --color={params.color} \
            --thresholds='{params.sex_thresholds}' \
            --sex_ribbons='{params.sex_ribbons}' \
            --width={params.width} \
            --height={params.height} \
            --show_numbers={params.show_numbers}
        """


rule qualimap:
    """ 
    Run qualimap
    """
    input:
        bam="{folder}/{file}.bam",
    output:
        directory("{folder}/04_stats/01_sparse_stats/{file}_qualimap"),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["stats", "qualimap"], attempt, 4
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["stats", "qualimap"], attempt, 1
        ),
    threads: get_threads2(["stats", "qualimap"], 4)
    log:
        "{folder}/04_stats/01_sparse_stats/{file}_qualimap.log",
    message:
        "--- QUALIMAP on {input}"
    conda:
        "../envs/qualimap.yaml"
    envmodules:
        module_qualimap,
    shell:
        """
        qualimap bamqc -c -bam {input} -outdir {output} -nt {threads} --java-mem-size={resources.memory}M > {log}
        """


rule multiqc:
    """
    Running multiqc
    """
    input:
        files=get_files_4_multiqc,
    output:
        html=report(
            "{folder}/04_stats/02_separate_tables/{genome}/multiqc_mapache.html",
            category="MultiQC",
        ),
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc2(
            ["stats", "multiqc"], attempt, 4
        ),
        runtime=lambda wildcards, attempt: get_runtime_alloc2(
            ["stats", "multiqc"], attempt, 1
        ),
    params:
        config="workflow/report/multiqc_config.yaml",
        resultdir=RESULT_DIR,
    log:
        "{folder}/04_stats/02_separate_tables/{genome}/multiqc_mapache.log",
    conda:
        "../envs/multiqc.yaml"
    envmodules:
        module_multiqc,
    message:
        "--- MULTIQC of {genome}"
    shell:
        """
        ## run mutliqc
        multiqc -c {params.config} \
                -f \
                -n $(basename {output.html}) \
                -o $(dirname {output.html}) \
                --title 'Mapache report (genome {wildcards.genome})' \
                --cl-config "extra_fn_clean_trim: ['{params.resultdir}']" \
                {input.files} 2> {log};
        """
