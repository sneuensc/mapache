import pandas as pd
import numpy as np
import itertools
import pathlib


##########################################################################################
# global variables

# autosomes = {} # will be set in check_chromosome_names()

##########################################################################################
## all functions for main snakemake file

## getter for config file
def get_param1(key, default):
    return config.get(key, default)


def get_param2(key1, key2, default):
    return config.get(key1, {}).get(key2, default)


def get_param3(key1, key2, key3, default):
    return config.get(key1, {}).get(key2, {}).get(key3, default)


## setter for config file
def set_param1(key, value):
    config[key] = value


def set_param2(key1, key2, value):
    if key1 not in config:
        config[key1] = {}
    config[key1][key2] = value


def set_param3(key1, key2, key3, value):
    if key1 not in config:
        config[key1] = {}
    if key2 not in config[key1]:
        config[key1][key2] = {}
    config[key1][key2][key3] = value


##########################################################################################
## functions to evaluate python code if necessary

## eval single element if needed
def eval_if_possible(x):
    try:
        return eval(x)
    except:
        return x


## eval single element if needed and return it as a list
def eval_to_list(x):
    a = eval_if_possible(x)
    if isinstance(a, list):
        return a
    else:
        return [a]


## eval single element if needed and return it as a comma separated string
def eval_to_csv(x):
    return ",".join(list(map(str, eval_to_list(x))))


## transform a list to a comma separated string
def list_to_csv(x):
    if isinstance(x, list):
        return ",".join(x)
    else:
        return x


## replace any element of the list
def eval_list(x):
    if isinstance(x, list):
        return list(map(eval_if_possible, x))
    else:
        return [eval_if_possible(x)]


## replace any element of the list
def eval_list_to_csv(x):
    return ",".join(eval_list(x))


##########################################################################################
## function to test the chromosome names
def check_chromosome_names(GENOME):

    global autosomes
    sex_chr = "X"
    # print(f"""
    # Evaluating GENOME:
    #     {GENOME}
    # """)

    ## get all chromsome names from the reference GENOME
    fasta = get_param3("genome", GENOME, "fasta", "")
    if pathlib.Path(f"{fasta}.fai").exists():
        allChr = list(
            map(str, pd.read_csv(f"{fasta}.fai", header=None, sep="\t")[0].tolist())
        )
    elif pathlib.Path(f"results/00_reference/{GENOME}/{GENOME}.fasta.fai").exists():
        fasta = f"results/00_reference/{GENOME}/{GENOME}.fasta"
        allChr = list(
            map(str, pd.read_csv(f"{fasta}.fai", header=None, sep="\t")[0].tolist())
        )
    elif pathlib.Path(fasta).exists():
        cmd = f"grep '^>' {fasta} | cut -c2- | awk '{{print $1}}'"
        allChr = subprocess.check_output(cmd, shell=True, text=True).split()
    else:
        # print(f"ERROR: Reference genome 'genome:{GENOME}:fasta' does not exist!")
        os._exit(1)

    # check if chromosomes for which DoC was requested exist
    depth_chromosomes = config["genome"][GENOME].get("depth_chromosomes", "")
    if len(depth_chromosomes):
        chromosomes = depth_chromosomes.split(",")

    else:
        chromosomes = []

    for chr in chromosomes:
        if chr not in allChr:
            print(
                f"""
            WARNING: at least one chromosome specified in config[genome][{GENOME}][depth_chromosomes] does not exist in FASTA reference file. 
            chr: {chr}
            """
            )
            os._exit(1)

    # print(f"""
    # chromosomes which will have their DoC reported in main table: {chromosomes}
    # """)

    # check if the chromosomes specified in sex determination exist
    # sex chromosome
    if config["genome"][GENOME].get("sex_inference", {}).get("run", False):
        # print(f"    Checking if chromosomes specified in config file for sex inference exist in genome {GENOME}.")
        sex_chr = str(config["genome"][GENOME]["sex_inference"].get("params", {}).get("sex_chr", "X"))
        if sex_chr not in allChr:

            print(f"""
            WARNING: sex chromosome specified in config[genome][{GENOME}][sex_inference][params][sex_chr] does not exist in FASTA reference file. 
            sex_chr: {sex_chr}
            Please set the right chromosome name for this reference genome.
            If you do not plan to infer sex from the mappings to this genome, set the option "run" to False for this genome.
            """
            )
            os._exit(1)

        # autosomes specified for sex inference
        chromosomes = list(
            map(
                str,
                eval_to_list(
                    config["genome"][GENOME]["sex_inference"]
                    .get("params", {})
                    .get("autosomes", [i for i in range(1, 23)])
                ),
            )
        )

        for chr in chromosomes:
            if chr not in allChr:
                print(
                    f"""
                WARNING: at least one chromosome specified in config[genome][{GENOME}][sex_inference][params][autosomes] does not exist in FASTA reference file. 
                chr: {chr}
                Please set the right chromosome names for this reference genome.
                If you do not plan to infer sex from the mappings to this genome, set the option "run" to False for this genome.
                """
                )
                os._exit(1)
        autosomes = 'c("' + '","'.join(chromosomes) + '")'
        config["genome"][GENOME]["sex_inference"]["params"]["autosomes"] = autosomes

        print(
            f"""
        sex_chr: {sex_chr}
        autosomes: {autosomes}
        """
        )
    print(f"WELL DONE. The chromosome names are well specified for genome {GENOME}.")


## convert string to boolean
def str2bool(v):
    return str(v).lower() in ("yes", "true", "t", "1")


## get incremental memory allocation when jobs fail
## 'startStr' is the first memory allocation in GB
## input is in GB; output in MB; default is 2GB, but can be changed by a rule

## get incremental memory allocation when jobs fail
## 'start' is the first memory allocation in GB (default 4GB)
## input is in GB; output is in MB;
## global variable memory_increment_ratio defines by how much (ratio) the memory is increased if not defined specifically
def get_memory_alloc(module, attempt, default=2):
    mem_start = int(config.get(module, {}).get("mem", default))
    mem_incre = int(
        config.get(module, {}).get("mem_increment", memory_increment_ratio * mem_start)
    )
    return int(1024 * ((attempt - 1) * mem_incre + mem_start))


def convert_time(seconds):
    day = seconds // (24 * 3600)
    seconds = seconds % (24 * 3600)
    hour = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    return "%d-%02d:%02d:%02d" % (day, hour, minutes, seconds)


## get incremental time allocation when jobs fail
## 'start' is the first time allocation in hours (default 12h)
## input is in hours; output is in minutes;
## global variable runtime_increment_ratio defines by how much (ratio) the time is increased if not defined specifically
def get_runtime_alloc(module, attempt, default=12):
    time_start = int(config.get(module, {}).get("time", default))
    time_incre = int(
        config.get(module, {}).get(
            "time_increment", runtime_increment_ratio * time_start
        )
    )
    return int(60 * ((attempt - 1) * time_incre + time_start))


# return convert_time(60*60 * ((attempt-1) * time_incre + time_start))


## get the number of threads of the given parameter
def get_threads(module, default=1):
    return int(config.get(module, {}).get("threads", default))


## define how to quantify the deamination pattern
def get_damage(run_damage):
    files = []
    if run_damage == "bamdamage":
        for GENOME in genome:
            files += [
                (
                    "results/04_stats/01_sparse_stats/02_library/04_bamdamage/{SM}/{LB}.{GENOME}.{type}.{ext}"
                ).format(SM=row["SM"], LB=row["LB"], GENOME=GENOME, type=type, ext=ext)
                for index, row in all_libraries.iterrows()
                for type in ["dam", "length"]
                for ext in ["pdf", "svg"]
            ]
    elif run_damage == "mapDamage":
        for GENOME in genome:
            files += [
                (
                    "results/04_stats/01_sparse_stats/02_library/04_mapDamage/{LB}.{GENOME}_results_mapDamage/Fragmisincorporation_plot.pdf"
                ).format(SM=row["SM"], LB=row["LB"], GENOME=GENOME)
                for index, row in all_libraries.iterrows()
            ]
    return files


##########################################################################################
##########################################################################################
## all functions for fastq files


def get_fastq_of_ID(wildcards):
    if "_R1" == wildcards.ID[-3:]:
        filename = samples[wildcards.SM][wildcards.LB][wildcards.ID[:-3]]["Data1"]
    elif "_R2" == wildcards.ID[-3:]:
        filename = samples[wildcards.SM][wildcards.LB][wildcards.ID[:-3]]["Data2"]
    elif paired_end:
    #elif paired_end != 0:  ## SE library in a paired-end sample file
        filename = samples[wildcards.SM][wildcards.LB][wildcards.ID]["Data1"]
    else:
        filename = samples[wildcards.SM][wildcards.LB][wildcards.ID]["Data"]
    return filename


def get_fastq_for_mapping(wildcards):
    if run_adapter_removal:
        if (
            # paired_end == 1
            collapse
            and str(samples[wildcards.SM][wildcards.LB][wildcards.ID]["Data2"]) != "nan"
        ):
            folder = f"results/01_fastq/01_trimmed/01_files_trim_collapsed/{wildcards.SM}/{wildcards.LB}"
        else:
            folder = f"results/01_fastq/01_trimmed/01_files_trim/{wildcards.SM}/{wildcards.LB}"
    else:
        folder = (
            f"results/01_fastq/00_reads/01_files_orig/{wildcards.SM}/{wildcards.LB}"
        )

    if not collapse:
    #if paired_end == 2:
        filename = [
            f"{folder}/{wildcards.ID}_R1.fastq.gz",
            f"{folder}/{wildcards.ID}_R2.fastq.gz",
        ]
    else:
        filename = [f"{folder}/{wildcards.ID}.fastq.gz"]
    return filename


def get_fastq_for_mapping_pe(wildcards):
    if run_adapter_removal:
        folder = f"results/01_fastq/01_trimmed/01_files_trim"
    else:
        folder = f"results/01_fastq/00_reads/01_files_orig"

    return f"{folder}/{wildcards.SM}/{wildcards.LB}/{wildcards.ID}_R{wildcards.id_read}.fastq.gz"


def get_bam_for_sorting(wildcards):
    if mapper == "bwa_aln":
        if (
            # paired_end == 2
            not collapse
            and str(samples[wildcards.SM][wildcards.LB][wildcards.ID]["Data2"]) != "nan"
        ):
            folder = "02_bwa_sampe"
        else:
            folder = "02_bwa_samse"
    elif mapper == "bwa_mem":
        folder = "02_bwa_mem"
    elif mapper == "bowtie2":
        folder = "02_bowtie2"
    else:
        print(
            f"ERROR: The parameter mapper is not correctly specified: {mapper} is unknown!"
        )
        os._exit(0)
    return f"results/01_fastq/02_mapped/{folder}/{wildcards.SM}/{wildcards.LB}/{wildcards.ID}.{wildcards.GENOME}.bam"


def get_final_bam_fastq(wildcards):
    if run_filtering:
        bam = f"{wildcards.folder}/03_filtered/01_bam_filter/{wildcards.SM}/{wildcards.LB}/{wildcards.ID}.{wildcards.GENOME}.bam"
    else:
        bam = f"{wildcards.folder}/02_mapped/03_bam_sort/{wildcards.SM}/{wildcards.LB}/{wildcards.ID}.{wildcards.GENOME}.bam"
    return bam


##########################################################################################
##########################################################################################
## all functions for libraries


def get_final_bam_library(wildcards):
    if run_mapDamage_rescale:
        bam = f"results/02_library/02_rescaled/01_mapDamage/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    elif run_mark_duplicates:
        if save_duplicates == "extract":
            bam = f"results/02_library/01_duplicated/01_rmdup/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}_mapped.bam"
        else:
            bam = f"results/02_library/01_duplicated/01_rmdup/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    else:
        bam = f"results/02_library/00_merged_fastq/01_bam/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    return bam


def get_mapDamage_bam(wildcards, index=False):
    if run_mark_duplicates:
        if save_duplicates == "extract":
            bam = f"results/02_library/01_duplicated/01_rmdup/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}_mapped.bam"
        else:
            bam = f"results/02_library/01_duplicated/01_rmdup/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}.bam"
    else:
        bam = f"results/02_library/00_merged_fastq/01_bam/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}.bam"
    if index:
        bam = bam.replace(".bam", ".bai")
    return bam


# if extract_duplicates:
#     bam = f"results/02_library/01_duplicated/01_rmdup/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}_mapped.bam"
# elif run_remove_duplicates:
#     bam = f"results/02_library/01_duplicated/01_rmdup/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}.bam"
# else:
#     bam = f"results/02_library/00_merged_fastq/01_bam/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}.bam"


##########################################################################################
##########################################################################################

##########################################################################################
## all functions for the samples


def get_final_bam(wildcards):
    if run_compute_md:
        bam = f"results/03_sample/02_md_flag/01_md_flag/{wildcards.SM}.{wildcards.GENOME}.bam"
    elif run_realign:
        bam = f"results/03_sample/01_realigned/01_realign/{wildcards.SM}.{wildcards.GENOME}.bam"
    else:
        bam = f"results/03_sample/00_merged_library/01_bam/{wildcards.SM}.{wildcards.GENOME}.bam"
    return bam


def get_md_flag_bam(wildcards):
    if run_realign:
        bam = f"results/03_sample/01_realigned/01_realign/{wildcards.SM}.{wildcards.GENOME}.bam"
    else:
        bam = f"results/03_sample/00_merged_library/01_bam/{wildcards.SM}.{wildcards.GENOME}.bam"
    return bam


## the function makes a reverse symlink, by moving the file to the new location and then symlink it back to the original location
def symlink_rev2(input, output):
    shell("mv {input} {output}")
    shell("ln -srf {output} {input}")
    shell("touch {output}")


def symlink_rev(input, output):
    shell("ln -srf {input} {output}")


##########################################################################################
##########################################################################################

##########################################################################################
## all functions for the stats


def is_quick(file_name, dict):
    if "quick" in dict.keys() and dict["quick"]:
        file_name.replace(".bam", ".downsampled.bam")
    return file_name


# sex_params = config["genome"]["sex"]["sex_params"] if "sex_params" in config["genome"]["sex"].keys() else {}


# def get_sex_params(wildcards):
#     #sex_params = get_param2("genome", wildcards.GENOME, {})
#     sex_dict = config["genome"][wildcards.GENOME].get("sex_inference", {}).get("params", {})

#     # autosomes are passed as a python expression, which was parsed previously in check_chromosome_names()
#     # the R script to assign sex needs the R expression that was stored in autosomes
#     if "autosomes" in sex_dict:
#         del sex_dict["autosomes"]
#         sex_dict["autosomes"] = autosomes[wildcards.GENOME]

#     sex_params = " ".join(
#         [f"--{key}='{sex_dict[key]}'" for key in sex_dict.keys()]
#     )
#     return sex_params


## get the individual depth files to combien them
def get_depth_files(wildcards):
    parts = pathlib.Path(wildcards.folder).parts
    if parts[1] == "02_library":
        filename = (
            [
                f"results/02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}_depth.txt"
                for GENOME in genome
                for SM in samples
                for LB in samples[SM]
            ],
        )
    else:
        filename = (
            [
                f"results/03_sample/03_final_sample/01_bam/{SM}.{GENOME}_depth.txt"
                for GENOME in genome
                for SM in samples
            ],
        )
    return filename


# def get_chrom(wildcards):
#     chr = eval_list(
#         "".join(get_param3("stats", "sample", "depth_chromosomes", "").split()).split(
#             ","
#         )
#     )
#     GENOME = get_param2("genome", wildcards.GENOME, {})
#     chr_uniq = list(set(chr) - set(list(GENOME)))
#     chr_def = list(set(chr).intersection(set(GENOME)))
#     return ",".join(
#         eval_list(chr_uniq) + [eval_if_possible(GENOME[c]) for c in chr_def]
#     )


def path_stats_by_level(wildcards):
    if wildcards.level == "FASTQ":
        paths = [
            f"results/04_stats/02_separate_tables/{wildcards.GENOME}/{SM}/{LB}/{ID}/fastq_stats.csv"
            # for GENOME in genome
            for SM in samples
            for LB in samples[SM]
            for ID in samples[SM][LB]
        ]
    elif wildcards.level == "LB":
        paths = [
            f"results/04_stats/02_separate_tables/{wildcards.GENOME}/{SM}/{LB}/library_stats.csv"
            # for GENOME in genome
            for SM in samples
            for LB in samples[SM]
        ]
    elif wildcards.level == "SM":
        paths = [
            f"results/04_stats/02_separate_tables/{wildcards.GENOME}/{SM}/sample_stats.csv"
            # for GENOME in genome
            for SM in samples
        ]
    return paths
