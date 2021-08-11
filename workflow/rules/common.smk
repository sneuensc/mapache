import pandas as pd
import numpy as np
import itertools
import pathlib


##########################################################################################
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
def check_chromsome_names(GENOME):
    ## get all chromsome names from the reference GENOME
    fasta = get_param3('genome', GENOME, 'fasta', '')
    if f"{fasta}.fai".exists():
        allChr = list(map(str, pd.read_csv(f"{fasta}.fai", header=None, sep="\t")[0].tolist()))
    elif f"results/00_reference/{GENOME}/{GENOME}.fasta.fai".exists():
        fasta = f"results/00_reference/{GENOME}/{GENOME}.fasta"
        allChr = list(map(str, pd.read_csv(f"{fasta}.fai", header=None, sep="\t")[0].tolist()))
    else: ## if the .fai file is not yet present
        cmd = f"grep '^>' {fasta} | cut -c2- | awk '{{print $1}}'"
        allChr = subprocess.check_output(cmd, shell=True, text=True).split()
    
    ## check female chromosome
    femaleChr = get_param3("genome", GENOME, "femaleChr", "X")
    if femaleChr not in allChr:
        set_param3("genome", GENOME, "femaleChr", "")
        print(f"WARNING: In parameter 'genome:{GENOME}:femaleChr' the chromosome name '{femaleChr}' is unknown, assuming no female chromosome!")
    
    ## check male chromosome
    maleChr = get_param3("genome", GENOME, "maleChr", "Y")
    if maleChr not in allChr:
        set_param3("genome", GENOME, "maleChr", "")
        print(f"WARNING: In parameter 'genome:{GENOME}:maleChr' the chromosome name '{maleChr}' is unknown, assuming no male chromosome!")
    
    ## check MT chromosome
    mtChr = get_param3("genome", GENOME, "mtChr", "MT")
    if mtChr not in allChr:
        set_param3("genome", GENOME, "mtChr", "")
        print(f"WARNING: In parameter 'genome:{GENOME}:mtChr' the chromosome name '{mtChr}' is unknown, assuming no MT chromosome!")
    
    ## check autosomes
    autosomeChr = list(map(str, eval_to_list(get_param3("genome", GENOME, "autosomeChr", ""))))
    if autosomeChr == "":
        ## if empty, autosomes are all not differently defined chromsomes
        autosomeChr = np.intersect1d(allChr, [femaleChr, maleChr, mtChr])
        set_param3("genome", GENOME, "autosomeChr", list_to_csv(autosomeChr))
    else:
        unknown = list(set(autosomeChr) - set(allChr))
        if len(unknown) > 0:
            autosomeChr = np.intersect1d(allChr, autosomeChr)
            set_param3("genome", GENOME,"autosomeChr",list_to_csv(autosomeChr),)
            if len(unknown) == 1:
                print(f"WARNING: In parameter 'genome:{GENOME}:autosomeChr' the chromosome name {unknown} is unknown, ignoring it!")
            else:
                print(f"WARNING: In parameter 'genome:{GENOME}:autosomeChr' the chromosome names {unknown} are unknown, ignoring them!")


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
# def get_damage(run_damage):
#     if run_damage == 'bamdamage':
#         files = [f"results/03_sample/{SM}/{LB}/library_bamdamage/{LB}.{GENOME}.dam.svg"
#            for GENOME in genome
#            for SM in samples
#            for LB in samples[SM]]
#     elif run_damage == 'mapDamage':
#         files = [f"results/03_sample/{SM}/{LB}/library_mapDamage/{LB}.{GENOME}_results_mapDamage/Fragmisincorporation_plot.pdf"
#            for GENOME in genome
#            for SM in samples
#            for LB in samples[SM]]
#
#     return (files)


##########################################################################################
##########################################################################################
## all functions for fastq files


def get_fastq_of_ID(wildcards):
    if "_R1" == wildcards.ID[-3:]:
        filename = samples[wildcards.SM][wildcards.LB][wildcards.ID[:-3]]["Data1"]
    elif "_R2" == wildcards.ID[-3:]:
        filename = samples[wildcards.SM][wildcards.LB][wildcards.ID[:-3]]["Data2"]
    elif paired_end != 0:  ## SE library in a paired-end sample file
        filename = samples[wildcards.SM][wildcards.LB][wildcards.ID]["Data1"]
    else:
        filename = samples[wildcards.SM][wildcards.LB][wildcards.ID]["Data"]
    return filename


def get_fastq_for_mapping(wildcards):
    if run_adapter_removal:
        if (
            paired_end == 1
            and str(samples[wildcards.SM][wildcards.LB][wildcards.ID]["Data2"]) != "nan"
        ):
            folder = f"results/01_fastq/01_trimmed/01_files_trim_collapsed/{wildcards.SM}/{wildcards.LB}"
        else:
            folder = f"results/01_fastq/01_trimmed/01_files_trim/{wildcards.SM}/{wildcards.LB}"
    else:
        folder = (
            f"results/01_fastq/00_reads/01_files_orig/{wildcards.SM}/{wildcards.LB}"
        )

    if paired_end == 2:
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
            paired_end == 2
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


def get_mapDamage_bam(wildcards):
    if run_mark_duplicates:
        if save_duplicates == "extract":
            bam = f"results/02_library/01_duplicated/01_rmdup/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}_mapped.bam"
        else:
            bam = f"results/02_library/01_duplicated/01_rmdup/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    else:
        bam = f"results/02_library/00_merged_fastq/01_bam/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    return bam


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


def get_sex_params(wildcards):
    sex_params = get_param2("genome", wildcards.GENOME, {})
    x = " ".join(
        [f"--{key}='{eval_to_csv(sex_params[key])}'" for key in sex_params.keys()]
    )
    return x


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


def get_chrom(wildcards):
    chr = eval_list(
        "".join(get_param3("stats", "sample", "depth_chromosomes", "").split()).split(
            ","
        )
    )
    GENOME = get_param2("genome", wildcards.GENOME, {})
    chr_uniq = list(set(chr) - set(list(GENOME)))
    chr_def = list(set(chr).intersection(set(GENOME)))
    return ",".join(
        eval_list(chr_uniq) + [eval_if_possible(GENOME[c]) for c in chr_def]
    )


def path_stats_by_level(wildcards):
    if wildcards.level == "FASTQ":
        paths = [
            f"results/04_stats/02_separate_tables/{GENOME}/{SM}/{LB}/{ID}/fastq_stats.csv"
            for GENOME in genome
            for SM in samples
            for LB in samples[SM]
            for ID in samples[SM][LB]
        ]
    elif wildcards.level == "LB":
        paths = [
            f"results/04_stats/02_separate_tables/{GENOME}/{SM}/{LB}/library_stats.csv"
            for GENOME in genome
            for SM in samples
            for LB in samples[SM]
        ]
    elif wildcards.level == "SM":
        paths = [
            f"results/04_stats/02_separate_tables/{GENOME}/{SM}/sample_stats.csv"
            for GENOME in genome
            for SM in samples
        ]
    return paths
