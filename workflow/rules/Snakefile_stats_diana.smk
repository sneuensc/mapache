#-----------------------------------------------------------------------------#
import pandas as pd
configfile: "config/config.yaml"

genomes = list(config["GENOME"].keys())

input_data = pd.read_csv("samples.txt", sep = "\s+")
samples = list(input_data.SM.unique())

sm_lb_id = {
    sample: {
        lib: {
            fastq for fastq in input_data[input_data.SM == sample][input_data.LB == lib]["ID"]
        }
        for lib in input_data[input_data.SM == sample]["LB"]
    }
    for sample in samples
}


#-----------------------------------------------------------------------------#

rule all_stats_by_level:
    input:
        ["results/04_stats/03_final_tables/SM.csv",
        "results/04_stats/03_final_tables/LB.csv",
        "results/04_stats/03_final_tables/FASTQ.csv"]

#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
## get sparse stats stats

rule fastqc:
    """
    Quality control of fastq file by fastqc (SE or R1)
    """
    input:
        "results/{file}.fastq.gz"
    output:
        html="results/04_stats/01_sparse_stats/{file}_fastqc.html",
        zip="results/04_stats/01_sparse_stats/{file}_fastqc.zip"
    log:
        "results/logs/04_stats/01_sparse_stats/{file}_fastqc.log"
    # resources:
    #     memory=lambda wildcards, attempt: get_memory_alloc("fastqc_mem", attempt, 2),
    #     runtime=lambda wildcards, attempt: get_runtime_alloc("fastqc_time", attempt, 1)
    # conda:
    # 	"../envs/fastqc.yaml"
    # envmodules:
    # 	module_fastqc
    message: "--- FASTQC {input}"
    shell:
        "fastqc --quiet --outdir $(dirname {output.html}) {input} 2> {log}"

rule samtools_flagstat:
    """
    Compute samtools flagstat on bam file
    """
    input:
        bam="results/{file}.bam"
    output:
        "results/04_stats/01_sparse_stats/{file}_flagstat.txt"
    # resources:
    #     memory=lambda wildcards, attempt: get_memory_alloc("samtools_flagstat_mem", attempt, 2),
    #     runtime=lambda wildcards, attempt: get_runtime_alloc("samtools_flagstat_time", attempt, 1)
    log:
        "results/logs/04_stats/01_sparse_stats/{file}_flagstat.log"
    # conda:
    #     "../envs/samtools.yaml"
    # envmodules:
    #     module_samtools
    message: "--- SAMTOOLS FLAGSTAT {input}"        
    shell:
        "samtools flagstat --threads {threads} {input} > {output} 2> {log};"

rule multiqc_fastqc:
    """
    Merging fastqc output with multiqc
    """
    input:
        # {folder} is "01_fastq/00_reads/01_files_orig" or "01_fastq/01_trimmed/01_files_trim"
        lambda wildcards: [f"results/04_stats/01_sparse_stats/{wildcards.folder}/{SM}/{LB}/{ID}_fastqc.zip" for SM in samples for LB in sm_lb_id[SM] for ID in sm_lb_id[SM][LB]]
    output:
        html =  "results/04_stats/01_sparse_stats/{folder}/multiqc_fastqc.html",
        txt =  "results/04_stats/01_sparse_stats/{folder}/multiqc_fastqc_data/multiqc_fastqc.txt"
        #html = report("results/{dir_stats}/{folder}/multiqc_fastqc.html", category=" Quality control"),
    # resources:
    #     memory=lambda wildcards, attempt: get_memory_alloc("multiqc_mem", attempt, 2),
    #     runtime=lambda wildcards, attempt: get_runtime_alloc("multiqc_time", attempt, 1)
    log:
        "results/logs/04_stats/{folder}/multiqc_fastqc.log"
    # conda:
    #     "../envs/multiqc.yaml"
    #message: "--- MULTIQC fastqc: {input}"
    shell:
        """
        multiqc -n $(basename {output.html}) -f -d -o $(dirname {output.html}) {input}  2> {log}
        """

rule bedtools_genomecov:
    input:
        bam = "results/{dir}/{file}.bam"
    output:
        genomecov = "results/04_stats/01_sparse_stats/{dir}/{file}.genomecov"
    shell:
        """
        bedtools genomecov -ibam {input.bam} > {output.genomecov}
        """

rule read_length:
    input:
        bam = "results/{file}.bam"
    output:
        length = "results/04_stats/01_sparse_stats/{file}.length"
    shell:
        """
        samtools view {input.bam} | workflow/scripts/read_length.pl -o {output.length}
        """

def is_quick(file_name, dict):
    if "quick" in dict.keys() and dict["quick"]:
        file_name.replace(".bam", ".downsampled.bam")
    return file_name


sex_params = config["stats"]["sex"]["sex_params"] if "sex_params" in config["stats"]["sex"].keys() else {}

rule assign_sex:
    input:
        genomecov = "results/04_stats/01_sparse_stats/{file}.genomecov"
    output:
        sex = "results/04_stats/01_sparse_stats/{file}.sex"
    params:
        sex_params = " ".join([f"--{key}={sex_params[key]}" for key in sex_params.keys()])
    shell:
        """
        Rscript workflow/scripts/sex_assignation.r \
            --genomecov={input.genomecov} \
            --out={output.sex} \
            {params.sex_params}
        """

#-----------------------------------------------------------------------------#
## merging individual stats

rule merge_stats_per_fastq:
    input:
        multiqc_orig                = "results/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/multiqc_fastqc_data/multiqc_fastqc.txt",   # raw sequenced reads
        multiqc_trim                = "results/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_files_trim/multiqc_fastqc_data/multiqc_fastqc.txt", # raw trimmed reads
        flagstat_mapped_highQ       = "results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/{SM}/{LB}/{ID}.{genome}_flagstat.txt",       # mapped and high-qual reads
        length_fastq_mapped_highQ   = "results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/{SM}/{LB}/{ID}.{genome}.length",
    output:
        "results/04_stats/02_separate_tables/{genome}/{SM}/{LB}/{ID}/stats.csv"
    shell:
        """
        Rscript workflow/scripts/merge_stats_per_fastq.R \
            --ID={wildcards.ID} \
            --LB={wildcards.LB} \
            --SM={wildcards.SM} \
            --genome={wildcards.genome} \
            --output_file={output} \
            --path_multiqc_orig={input.multiqc_orig} \
            --path_multiqc_trim={input.multiqc_trim} \
            --path_flagstat_mapped_highQ={input.flagstat_mapped_highQ} \
            --path_length_mapped_highQ={input.length_fastq_mapped_highQ}
        """          

rule merge_stats_per_lb:
    input:
        fastq_stats         = lambda wildcards: [f"results/04_stats/02_separate_tables/{wildcards.genome}/{wildcards.SM}/{wildcards.LB}/{ID}/stats.csv" for ID in sm_lb_id[wildcards.SM][wildcards.LB]],
        flagstat_raw        = "results/04_stats/01_sparse_stats/02_library/00_merged_fastq/01_bam/{SM}/{LB}.{genome}_flagstat.txt",       
        flagstat_unique     = "results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{genome}_flagstat.txt",      
        length_unique       = "results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{genome}.length",      
        genomecov_unique    = "results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{genome}.genomecov",     
        sex_unique          = "results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{SM}/{LB}.{genome}.sex",     
    output:
        "results/04_stats/02_separate_tables/{genome}/{SM}/{LB}/stats.csv"
    params:
        chrs_selected = "--chrs_selected=" + config["stats"]["library"]["depth_chromosomes"] if "depth_chromosomes" in config["stats"]["library"] else ""
    shell:
        """
        list_fastq_stats=$(echo {input.fastq_stats} |sed 's/ /,/g')
        Rscript workflow/scripts/merge_stats_per_LB.R \
            --LB={wildcards.LB} \
            --SM={wildcards.SM} \
            --genome={wildcards.genome} \
            --output_file={output} \
            --path_list_stats_fastq=${{list_fastq_stats}} \
            --path_flagstat_raw={input.flagstat_raw} \
            --path_flagstat_unique={input.flagstat_unique} \
            --path_length_unique={input.length_unique} \
            --path_genomecov_unique={input.genomecov_unique} \
            --path_sex_unique={input.sex_unique} \
            {params.chrs_selected}
        """

rule merge_stats_per_sm:
    input:
        lb_stats           = lambda wildcards: [f"results/04_stats/02_separate_tables/{wildcards.genome}/{wildcards.SM}/{LB}/stats.csv" for LB in sm_lb_id[wildcards.SM]],
        flagstat_unique    = "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{genome}_flagstat.txt",      
        length_unique      = "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{genome}.length",
        genomecov_unique   = "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{genome}.genomecov",
        sex_unique         = "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{genome}.sex" 
    output:
        "results/04_stats/02_separate_tables/{genome}/{SM}/stats.csv"
    params:
        chrs_selected = "--chrs_selected=" + config["stats"]["sample"]["depth_chromosomes"] if "depth_chromosomes" in config["stats"]["sample"] else ""
    shell:
        """
        list_lb_stats=$(echo {input.lb_stats} |sed 's/ /,/g')
        Rscript workflow/scripts/merge_stats_per_SM.R \
            --SM={wildcards.SM} \
            --genome={wildcards.genome} \
            --output_file={output} \
            --path_list_stats_fastq=${{list_lb_stats}} \
            --path_flagstat_unique={input.flagstat_unique} \
            --path_length_unique={input.length_unique} \
            --path_genomecov_unique={input.genomecov_unique} \
            --path_sex_unique={input.sex_unique} \
            {params.chrs_selected}
        """


#-----------------------------------------------------------------------------#
# merge all the stats in the same level

def path_stats_by_level(level):
    if level == "FASTQ":
        paths = [f"results/04_stats/02_separate_tables/{genome}/{SM}/{LB}/{ID}/stats.csv" for genome in genomes for SM in samples for LB in sm_lb_id[SM] for ID in sm_lb_id[SM][LB]]
    elif level == "LB":
        paths = [f"results/04_stats/02_separate_tables/{genome}/{SM}/{LB}/stats.csv" for genome in genomes for SM in samples for LB in sm_lb_id[SM]]
    elif level == "SM":
        paths = [f"results/04_stats/02_separate_tables/{genome}/{SM}/stats.csv" for genome in genomes for SM in samples]
    return paths


rule merge_stats_by_level:
    input:
        paths = lambda wildcards: path_stats_by_level(wildcards.level)
    output:
        "results/04_stats/03_final_tables/{level}.csv"
    run:
        import pandas as pd
        df_list = [pd.read_csv(file) for file in input]
        df = pd.concat(df_list)
        df.to_csv(str(output), index = False)

