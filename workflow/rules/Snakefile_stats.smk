# -----------------------------------------------------------------------------#
import pandas as pd


# -----------------------------------------------------------------------------#
## get sparse stats stats


rule fastqc:
    """
    Quality control of fastq file by fastqc (SE or R1)
    """
    input:
        get_fastq_for_mapping
    output:
        html = "results/04_stats/01_sparse_stats/01_fastq/{group1}/{group2}/{SM}/{LB}/{ID}_fastqc.html",
        zip =  "results/04_stats/01_sparse_stats/01_fastq/{group1}/{group2}/{SM}/{LB}/{ID}_fastqc.zip",
    log:
        "results/04_stats/01_sparse_stats/01_fastq/{group1}/{group2}/{SM}/{LB}/{ID}_fastqc.log",
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


# rule multiqc_fastqc:
#     """
#     Merging fastqc output with multiqc
#     """
#     input:
#         lambda wildcards: [f"results/04_stats/01_sparse_stats/{wildcards.folder}/{SM}/{LB}/{ID}_fastqc.zip"
#             for SM in samples
#             for LB in samples[SM]
#             for ID in samples[SM][LB]
#         ],
#     output:
#         html=report("results/04_stats/01_sparse_stats/{folder}/multiqc_fastqc.html",
#             "results/{dir_stats}/{folder}/multiqc_fastqc.html",
#             category=" Quality control",
#         ),
#         txt="results/04_stats/01_sparse_stats/{folder}/multiqc_fastqc_data/multiqc_fastqc.txt",
#     resources:
#         memory=lambda wildcards, attempt: get_memory_alloc("multiqc_mem", attempt, 2),
#         runtime=lambda wildcards, attempt: get_runtime_alloc("multiqc_time", attempt, 1),
#     log:
#         "results/04_stats/01_sparse_stats/{folder}/multiqc_fastqc.log",
#     conda:
#         "../envs/multiqc.yaml"
#     envmodules:
#         module_multiqc,
#     message:
#         "--- MULTIQC fastqc: {input}"
#     shell:
#         """
#         multiqc -n $(basename {output.html}) -f -d -o $(dirname {output.html}) {input}  2> {log}
#         """


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


rule assign_sex:
    input:
        genomecov="results/04_stats/01_sparse_stats/{file}.{GENOME}.genomecov",
        test="results/00_reference/{GENOME}/{GENOME}.ok",
    output:
        sex="results/04_stats/01_sparse_stats/{file}.{GENOME}.sex",
    params:
        sex_params = get_sex_params
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
        fastqc_orig="results/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/{SM}/{LB}/{ID}_fastqc.zip",  # raw sequenced reads
        fastqc_trim="results/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_files_trim/{SM}/{LB}/{ID}_fastqc.zip",  # raw trimmed reads
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
            --path_fastqc_orig={input.fastqc_orig} \
            --path_fastqc_trim={input.fastqc_trim} \
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


rule merge_stats_per_lb:
    input:
        fastq_stats=lambda wildcards: [
            f"results/04_stats/02_separate_tables/{wildcards.GENOME}/{wildcards.SM}/{wildcards.LB}/{ID}/fastq_stats.csv" for ID in samples[wildcards.SM][wildcards.LB]],
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
        lb_stats = lambda wildcards: [f"results/04_stats/02_separate_tables/{wildcards.GENOME}/{wildcards.SM}/{LB}/library_stats.csv"
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

# Here level can be: SM, LB, FASTQ 
rule merge_stats_by_level:
    input:
        paths = path_stats_by_level,
    output:
        "results/03_sample/04_stats/01_summary/{level}_stats.{GENOME}.csv",
    log:
        "results/03_sample/04_stats/01_summary/{level}_stats.{GENOME}.log",
    message:
        "--- MERGE STATS by {wildcards.level}"
    run:
        import pandas as pd

        df_list = [pd.read_csv(file) for file in input]
        df = pd.concat(df_list)
        df.to_csv(str(output), index=False)



##########################################################################################
# read depth by chromosome
rule DoC_chr_SM:
    input:
        "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{genome}.genomecov"
    output:
        "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{genome}_DoC_chrs.csv"
    params:
        SM = "{SM}"
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
            genome = wildcards.genome,
            SM = samples
            )
    output:
        "03_sample/04_stats/01_summary/DoC_by_chrs.{genome}.csv"
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
        bai=lambda wildcards: get_mapDamage_bam(wildcards, index = True)
    output:
        damage_pdf="results/02_library/04_stats/03_bamdamage/{id_sample}/{id_library}/{id_library}.{id_genome}.dam.pdf",
        length_pdf="results/02_library/04_stats/03_bamdamage/{id_sample}/{id_library}/{id_library}.{id_genome}.length.pdf",
        length_table=report("results/02_library/04_stats/03_bamdamage/{id_sample}/{id_library}/{id_library}.{id_genome}.length.csv", category="Read length table"),
        dam_5prime_table=report("results/02_library/04_stats/03_bamdamage/{id_sample}/{id_library}/{id_library}.{id_genome}.dam_5prime.csv", category="Damage pattern table"),
        dam_3prime_table=report("results/02_library/04_stats/03_bamdamage/{id_sample}/{id_library}/{id_library}.{id_genome}.dam_3prime.csv", category="Damage pattern table")      
    log:
        "logs/results/02_library/04_stats/03_bamdamage/{id_sample}/{id_library}/{id_library}.{id_genome}_stats.log"
    threads: 1
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapdamage_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("bamdamage_time", attempt, 24),        
    message: "--- BAMDAMAGE {input.bam}"
    params:
        prefix="results/02_library/04_stats/03_bamdamage/{id_sample}/{id_library}/{id_library}.{id_genome}",
        bamdamage_params = config["bamdamage_params"] if "bamdamage_params" in config.keys() else '',
        fraction = config["bamdamage_fraction"] if "bamdamage_fraction" in config.keys() else 0,
    envmodules:
        module_samtools,
        module_r
    shell:    
    	"""

    	nb=$(samtools idxstats {input.bam} | awk '{{sum += $3}} END {{print sum}}'); 
    	nth_line=1; 
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

##########################################################################################
# plots
# -----------------------------------------------------------------------------#
# plotting


# rule plot_summary_statistics:
#     """
#     Plot summary statistics
#     """
#     input:
#         fastq_stats="results/04_stats/03_final_tables/FASTQ.csv",
#         library_stats="results/04_stats/03_final_tables/LB.csv",
#         sample_stats="results/04_stats/03_final_tables/SM.csv",
#     output:
#         plot_1_nb_reads=report(
#             "results/04_stats/03_final_tables/04_final_plots/1_nb_reads.png",
#             caption="../report/1_nb_reads.rst",
#             category="Mapping statistics plots",
#         ),
#         plot_2_mapped=report(
#             "results/04_stats/03_final_tables/04_final_plots/2_mapped.png",
#             caption="../report/2_mapped.rst",
#             category="Mapping statistics plots",
#         ),
#         plot_3_endogenous=report(
#             "results/04_stats/03_final_tables/04_final_plots/3_endogenous.png",
#             caption="../report/3_endogenous.rst",
#             category="Mapping statistics plots",
#         ),
#         plot_4_duplication=report(
#             "results/04_stats/03_final_tables/04_final_plots/4_duplication.png",
#             caption="../report/4_duplication.rst",
#             category="Mapping statistics plots",
#         ),
#     log:
#         "results/04_stats/03_final_tables/04_final_plots/plot.log",
#     conda:
#         "../envs/r.yaml"
#     envmodules:
#         module_r,
#     message:
#         "--- PLOT SUMMARY STATISTICS"
#     shell:
#         """
#         Rscript workflow/scripts/sex_assignation.r \
#             --genomecov={input.genomecov} \
#             --out={output.sex} \
#             --larger_chr={params.sex_params}
#         """

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
    	module_r
    message: "--- PLOT DEPTH STATISTICS OF {input}"
    script:    	
    	"../scripts/plot_depth.R"

rule plot_summary_statistics:
    """
    Plot summary statistics
    """
    input:
        sample_stats = "results/03_sample/04_stats/01_summary/SM_stats.{id_genome}.csv"
        # fastq_stats = "results/01_fastq/05_stats/02_summary/fastq_stats.{id_genome}.csv",
        # library_stats = "results/02_library/04_stats/02_summary/library_stats.{id_genome}.csv",
    output: 
        plot_1_nb_reads = report("results/03_sample/04_stats/01_summary/1_nb_reads.{id_genome}.png", caption="../report/1_nb_reads.rst", category="Mapping statistics plots"),
        plot_2_mapped = report("results/03_sample/04_stats/01_summary/2_mapped.{id_genome}.png", caption="../report/2_mapped.rst", category="Mapping statistics plots"),        
        plot_3_endogenous = report("results/03_sample/04_stats/01_summary/3_endogenous.{id_genome}.png", caption="../report/3_endogenous.rst", category="Mapping statistics plots"),        
        # plot_4_duplication = report("results/03_sample/04_stats/01_summary/4_duplication.{id_genome}.png", caption="../report/4_duplication.rst", category="Mapping statistics plots")      
    log: 
        "results/03_sample/04_stats/01_summary/plot_summary_statistics_{id_genome}.log"
    conda:
    	"../envs/r.yaml"
    envmodules:
    	module_r
    message: "--- PLOT SUMMARY STATISTICS OF {wildcards.id_genome}"
    params:
        samples = config["sample_file"]
    shell:
    	"""
        Rscript workflow/scripts/plot_stats.R \
            --samples={params.samples} \
            --SM={input.sample_stats}  \
            --out_1_reads={output.plot_1_nb_reads} \
            --out_2_mapped={output.plot_2_mapped} \
            --out_3_endogenous={output.plot_3_endogenous}
        """
        
