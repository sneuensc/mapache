import pandas as pd
import numpy as np
import itertools
import pathlib
import re


##########################################################################################
## all functions for main snakemake file

## get recursively the argument of the given keys in a nested dict
def recursive_get(keys, def_value, my_dict=config):
    key = keys[0]
    if len(keys) == 1:
        value = my_dict.get(key, def_value)
    else:
        value = recursive_get(keys[1:], def_value, my_dict=my_dict.get(key, {}))
    return value


## same as above, but the argument is tested if it is present in the 'available_args'
## first argument of 'available_args' is the default
def recursive_get_and_test(key, available_args, my_dict=config):
    arg = str(recursive_get(key, available_args[0], my_dict))
    if arg not in available_args:
        LOGGER.error(
            f"ERROR: The parameter '{':'.join(key)}' has no valid argument (currently '{arg}'; available {available_args})!"
        )
        sys.exit(1)
    return arg


## update the arguemnt in a nested dict (and create leaves if needed)
def update_value(keys, value, my_dict=config):
    key = keys[0]
    if len(keys) == 1:
        new_dict = {}
        new_dict[key] = value
        return new_dict
    else:
        if key in my_dict:
            new_dict = update_value(keys[1:], value, my_dict[key])
            my_dict[key].update(new_dict)
            return my_dict
        else:
            my_dict[key] = {}
            new_dict = update_value(keys[1:], value, my_dict[key])
            my_dict[key].update(new_dict)
            return my_dict


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
def check_chromosome_names(GENOME, logging=True):
    if logging:
        LOGGER.info(f"  - Genome '{GENOME}':")

    ## test if fasta is valid
    fasta = recursive_get(["genome", GENOME, "fasta"], "")
    if "fasta" == "":
        LOGGER.error(f"ERROR: Reference genome '{GENOME}' has no fasta file defined!")
    elif not os.path.isfile(fasta):
        LOGGER.error(
            f"ERROR: Reference genome '{GENOME}': Fasta file '{fasta}' does not exist!"
        )

    ## get all chromsome names from the reference GENOME
    if pathlib.Path(f"{fasta}.fai").exists():
        allChr = list(
            map(str, pd.read_csv(f"{fasta}.fai", header=None, sep="\t")[0].tolist())
        )
    elif pathlib.Path(
        f"{RESULT_DIR}/00_reference/{GENOME}/{GENOME}.fasta.fai"
    ).exists():
        fasta = f"{RESULT_DIR}/00_reference/{GENOME}/{GENOME}.fasta"
        allChr = list(
            map(str, pd.read_csv(f"{fasta}.fai", header=None, sep="\t")[0].tolist())
        )
    elif pathlib.Path(fasta).exists():
        cmd = f"grep '^>' {fasta} | cut -c2- | awk '{{print $1}}'"
        allChr = map(str, subprocess.check_output(cmd, shell=True, text=True).split())
    else:
        LOGGER.error(f"ERROR: Reference genome 'genome:{GENOME}:fasta' does not exist!")
        os._exit(1)

    ## try to recognise if it is the human hg19 or GRCh38 genome. If so apply default chromosome names
    hg19 = map(str, list(range(1, 23)) + ["X", "Y", "MT"])
    GRCh38 = [f"chr{x}" for x in list(range(1, 23)) + ["X", "Y", "M"]]
    if sorted(allChr) == sorted(hg19):
        detectedChromosomes = ["X", "Y", "MT"]
        detectedAutosomes = list(set(allChr) - set(detectedChromosomes))
        if logging:
            LOGGER.info(
                f"    - Detected genome as 'hg19': Appliyng default chromosome names for sex and mt chromsomes: {detectedChromosomes}."
            )
    elif sorted(allChr) == sorted(GRCh38):
        detectedChromosomes = ["chrX", "chrY", "chrM"]
        detectedAutosomes = list(set(allChr) - set(detectedChromosomes))
        if logging:
            LOGGER.info(
                f"    - Detected genome as 'GRCh38': Appliyng default chromosome names for sex and mt chromsomes {detectedChromosomes}."
            )
    else:
        detectedChromosomes = ["", "", ""]
        detectedAutosomes = []

    # check if the chromosomes specified in sex determination exist
    # sex chromosome
    if recursive_get(["genome", GENOME, "sex_inference", "run"], False):
        if logging:
            LOGGER.info(f"    - Inferring sex")
        ## X chromosome specified for the sex inference
        sex_chr = recursive_get(
            ["genome", GENOME, "sex_inference", "params", "sex_chr"],
            detectedChromosomes[0],
        )
        if sex_chr not in allChr:
            LOGGER.error(
                f"ERROR: Sex chromosome specified in 'config[genome][{GENOME}][sex_inference][params][sex_chr]' ({sex_chr}) does not exist in FASTA reference genome."
            )
            os._exit(1)

        # autosomes specified for sex inference
        autosomes = list(
            map(
                str,
                eval_to_list(
                    recursive_get(
                        ["genome", GENOME, "sex_inference", "params", "autosomes"],
                        detectedAutosomes,
                    )
                ),
            )
        )
        if list(set(autosomes) - set(allChr)):
            LOGGER.error(
                f"ERROR: In 'config[genome][{GENOME}][sex_inference][params][autosomes]', the following chromsome names are not recognized: {list(set(autosomes) - set(allChr))}!"
            )
            os._exit(1)

        config = update_value(
            ["genome", GENOME, "sex_inference", "params", "autosomes"],
            'c("' + '","'.join(autosomes) + '")',
        )

    # check if chromosomes for which DoC was requested exist
    depth_chromosomes = recursive_get(["genome", GENOME, "depth_chromosomes"], "")
    chromosomes = depth_chromosomes.split(",") if len(depth_chromosomes) else []
    if list(set(chromosomes) - set(allChr)):
        LOGGER.error(
            f"In 'config[genome][{GENOME}][depth_chromosomes]', the following chromsome names are not recognized: {list(set(chromosomes) - set(allChr))}!"
        )
        os._exit(1)
    if logging and depth_chromosomes:
        LOGGER.info(f"    - Computing depth of chromosomes {depth_chromosomes}.")


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
    moduleList = module
    if type(moduleList) is not list:
        moduleList = [module]
    mem_start = int(recursive_get(moduleList + ["mem"], default))
    mem_incre = int(
        recursive_get(
            moduleList + ["mem_increment"], memory_increment_ratio * mem_start
        )
    )
    return int(1024 * ((attempt - 1) * mem_incre + mem_start))


## in this second verion the 'mem' is added to the word of the last element
def get_memory_alloc2(module, attempt, default=2):
    moduleList = module
    if type(moduleList) is not list:
        moduleList = [moduleList]
    mem_start = int(recursive_get(moduleList[:-1] + [moduleList[-1] + "_mem"], default))
    mem_incre = int(
        recursive_get(
            moduleList[:-1] + [moduleList[-1] + "_mem_increment"],
            memory_increment_ratio * mem_start,
        )
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
    moduleList = module
    if type(moduleList) is not list:
        moduleList = [module]
    time_start = int(recursive_get(moduleList + ["time"], default))
    time_incre = int(
        recursive_get(
            moduleList + ["time_increment"], runtime_increment_ratio * time_start
        )
    )
    return int(60 * ((attempt - 1) * time_incre + time_start))


## in this second verion the 'time' is added to the word of the last element
def get_runtime_alloc2(module, attempt, default=12):
    moduleList = module
    if type(moduleList) is not list:
        moduleList = [module]
    time_start = int(
        recursive_get(moduleList[:-1] + [moduleList[-1] + "_time"], default)
    )
    time_incre = int(
        recursive_get(
            moduleList[:-1] + [moduleList[-1] + "_time_increment"],
            runtime_increment_ratio * time_start,
        )
    )
    return int(60 * ((attempt - 1) * time_incre + time_start))


# return convert_time(60*60 * ((attempt-1) * time_incre + time_start))


## get the number of threads of the given parameter
def get_threads(module, default=1):
    return int(config.get(module, {}).get("threads", default))


## define how to quantify the deamination pattern
def get_damage(run_damage):
    if run_damage == "bamdamage":
        files = [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/02_library/04_bamdamage/{SM}/{LB}.{GENOME}.{type}.{ext}"
            for SM in samples
            for LB in samples[SM]
            for type in ["dam", "length"]
            for ext in ["pdf", "svg"]
            for GENOME in genome
        ]
    elif run_damage == "mapDamage":
        files = [
            f"{{RESULT_DIR}}/04_stats/01_sparse_stats/02_library/04_mapDamage/{SM}/{LB}.{GENOME}_results_mapDamage/Fragmisincorporation_plot.pdf"
            for SM in samples
            for LB in samples[SM]
            for GENOME in genome
        ]
    else:
        LOGGER.error(f"ERROR: def get_damage({run_damage}): should never happen!")
        sys.exit(1)
    return files


##########################################################################################
## check if java is called by a .jar file or by a wrapper
def get_picard_bin():
    bin = recursive_get(["software", "picard_jar"], "picard.jar")
    if bin[-4:] == ".jar":
        bin = "java -XX:ParallelGCThreads={threads} -XX:+UseParallelGC -XX:-UsePerfData \
            -Xms{resources.memory}m -Xmx{resources.memory}m -jar bin"
    return bin


def get_gatk_bin():
    bin = recursive_get(["software", "gatk3_jar"], "GenomeAnalysisTK.jar")
    if bin[-4:] == ".jar":
        bin = "java -XX:ParallelGCThreads={threads} -XX:+UseParallelGC -XX:-UsePerfData \
            -Xms{resources.memory}m -Xmx{resources.memory}m -jar bin"
    return bin


##########################################################################################
## all functions for fastq files


def get_fastq_of_ID(wildcards):
    if "_R1" == wildcards.ID[-3:]:
        filename = samples[wildcards.SM][wildcards.LB][wildcards.ID[:-3]]["Data1"]
    elif "_R2" == wildcards.ID[-3:]:
        filename = samples[wildcards.SM][wildcards.LB][wildcards.ID[:-3]]["Data2"]
    elif paired_end:
        # elif paired_end != 0:  ## SE library in a paired-end sample file
        filename = samples[wildcards.SM][wildcards.LB][wildcards.ID]["Data1"]
    else:
        filename = samples[wildcards.SM][wildcards.LB][wildcards.ID]["Data"]
    return filename


def get_fastq_for_mapping(wildcards, run_adapter_removal):
    if run_adapter_removal:
        folder = f"{wildcards.folder}/01_fastq/01_trimmed/01_adapter_removal/{wildcards.SM}/{wildcards.LB}"
        if not paired_end:
            filename = [f"{folder}/{wildcards.ID}.fastq.gz"]
        elif collapse:
            filename = rules.adapter_removal_collapse.output.R
        else:
            # if str(samples[wildcards.SM][wildcards.LB][wildcards.ID]["Data2"]) == "nan":
            data2 = recursive_get(
                [wildcards.SM, wildcards.LB, wildcards.ID, "Data2"],
                "nan",
                my_dict=samples,
            )
            if data2 != data2:
                # single-end files
                filename = [f"{folder}/{wildcards.ID}.fastq.gz"]
            else:
                # paired-end files, not collapsing
                filename = [
                    f"{folder}/{wildcards.ID}_R1.fastq.gz",
                    f"{folder}/{wildcards.ID}_R2.fastq.gz",
                ]
    else:
        folder = f"{wildcards.folder}/01_fastq/00_reads/01_files_orig/{wildcards.SM}/{wildcards.LB}"
        if not paired_end:
            filename = [f"{folder}/{wildcards.ID}.fastq.gz"]
        else:
            data2 = recursive_get(
                [wildcards.SM, wildcards.LB, wildcards.ID, "Data2"],
                "nan",
                my_dict=samples,
            )
            # checking a single-end file
            if data2 != data2:
                filename = [f"{folder}/{wildcards.ID}.fastq.gz"]
            else:
                filename = [
                    f"{folder}/{wildcards.ID}_R1.fastq.gz",
                    f"{folder}/{wildcards.ID}_R2.fastq.gz",
                ]

    return filename


def inputs_fastqc(wildcards):
    if "trim" in wildcards.folder:
        return get_fastq_for_mapping(wildcards, True)
    else:
        return get_fastq_for_mapping(wildcards, False)


def get_bam_for_sorting(wildcards):
    if mapper == "bwa_aln":
        if (
            not collapse
            # and str(samples[wildcards.SM][wildcards.LB][wildcards.ID]["Data2"]) != "nan"
            and str(
                recursive_get(
                    [wildcards.SM, wildcards.LB, wildcards.ID, "Data2"],
                    "nan",
                    my_dict=samples,
                )
            )
            != "nan"
        ):
            folder = "02_bwa_sampe"
        else:
            folder = "02_bwa_samse"
    elif mapper == "bwa_mem":
        folder = "02_bwa_mem"
    elif mapper == "bowtie2":
        folder = "02_bowtie2"
    else:
        LOGGER.error(
            f"ERROR: The parameter mapper is not correctly specified: {mapper} is unknown!"
        )
        os._exit(1)
    return f"{wildcards.folder}/01_fastq/02_mapped/{folder}/{wildcards.SM}/{wildcards.LB}/{wildcards.ID}.{wildcards.GENOME}.bam"


def get_final_bam_fastq(wildcards):
    if wildcards.type == "01_bam":
        if run_filtering:
            bam = f"{wildcards.folder}/01_fastq/03_filtered/01_bam_filter/{wildcards.SM}/{wildcards.LB}/{wildcards.ID}.{wildcards.GENOME}.bam"
        else:
            bam = f"{wildcards.folder}/01_fastq/02_mapped/03_bam_sort/{wildcards.SM}/{wildcards.LB}/{wildcards.ID}.{wildcards.GENOME}.bam"
    elif wildcards.type == "01_bam_low_qual":
        bam = f"{wildcards.folder}/01_fastq/03_filtered/01_bam_filter_low_qual/{wildcards.SM}/{wildcards.LB}/{wildcards.ID}.{wildcards.GENOME}.bam"
    else:
        LOGGER.error(
            f"ERROR: def get_final_bam_library({wildcards.type}): should never happen!"
        )
        os._exit(1)
    return bam


##########################################################################################
##########################################################################################
## all functions for libraries


def get_final_bam_library(wildcards):
    if wildcards.type == "01_bam":
        if run_damage_rescale:
            bam = f"{wildcards.folder}/02_library/02_rescaled/01_mapDamage/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
        elif remove_duplicates == "markduplicates":
            bam = f"{wildcards.folder}/02_library/01_duplicated/01_markduplicates/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
        elif remove_duplicates == "dedup":
            bam = f"{wildcards.folder}/02_library/01_duplicated/01_dedup/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
        else:
            bam = f"{wildcards.folder}/02_library/00_merged_fastq/01_bam/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    elif wildcards.type == "01_bam_low_qual":
        bam = f"{wildcards.folder}/02_library/00_merged_fastq/01_bam_low_qual/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    else:
        LOGGER.error(
            f"ERROR: def get_final_bam_library({wildcards.type}): should never happen!"
        )
        os._exit(1)
    return bam


def get_mapDamage_bam(wildcards):
    if remove_duplicates == "markduplicates":
        bam = f"{wildcards.folder}/02_library/01_duplicated/01_markduplicates/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    elif remove_duplicates == "dedup":
        bam = f"{wildcards.folder}/02_library/01_duplicated/01_dedup/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    else:
        bam = f"{wildcards.folder}/02_library/00_merged_fastq/01_bam/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}.bam"
    return bam


##########################################################################################
##########################################################################################

##########################################################################################
## all functions for the samples


def get_final_bam(wildcards):
    if run_compute_md:
        bam = f"{wildcards.folder}/03_sample/02_md_flag/01_md_flag/{wildcards.SM}.{wildcards.GENOME}.bam"
    elif run_realign:
        bam = f"{wildcards.folder}/03_sample/01_realigned/01_realign/{wildcards.SM}.{wildcards.GENOME}.bam"
    else:
        bam = f"{wildcards.folder}/03_sample/00_merged_library/01_bam/{wildcards.SM}.{wildcards.GENOME}.bam"
    return bam


def get_md_flag_bam(wildcards):
    if run_realign:
        bam = f"{wildcards.folder}/03_sample/01_realigned/01_realign/{wildcards.SM}.{wildcards.GENOME}.bam"
    else:
        bam = f"{wildcards.folder}/03_sample/00_merged_library/01_bam/{wildcards.SM}.{wildcards.GENOME}.bam"
    return bam


def symlink_rev(input, output):
    shell("ln -srf {input} {output}")


def get_sex_file(wildcards, level):
    if level == "SM":
        folder=f"{wildcards.folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{wildcards.SM}.{wildcards.GENOME}"
    else:
        folder=f"{wildcards.folder}/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{wildcards.SM}/{wildcards.LB}.{wildcards.GENOME}"

    if str2bool(
        recursive_get_and_test(
            ["genome", wildcards.GENOME, "sex_inference", "run"], ["False", "True"]
        )
    ):
        return f"{folder}.sex"
    else:
        return f"{folder}.nosex"



##########################################################################################
##########################################################################################

##########################################################################################
## all functions for the stats


def path_stats_by_level(wildcards):
    if wildcards.level == "FASTQ":
        paths = [
            f"{wildcards.folder}/04_stats/02_separate_tables/{wildcards.GENOME}/{SM}/{LB}/{ID}/fastq_stats.csv"
            # for GENOME in genome
            for SM in samples
            for LB in samples[SM]
            for ID in samples[SM][LB]
        ]
    elif wildcards.level == "LB":
        paths = [
            f"{wildcards.folder}/04_stats/02_separate_tables/{wildcards.GENOME}/{SM}/{LB}/library_stats.csv"
            # for GENOME in genome
            for SM in samples
            for LB in samples[SM]
        ]
    elif wildcards.level == "SM":
        paths = [
            f"{wildcards.folder}/04_stats/02_separate_tables/{wildcards.GENOME}/{SM}/sample_stats.csv"
            # for GENOME in genome
            for SM in samples
        ]
    else:
        LOGGER.error(
            f"ERROR: def path_stats_by_level({wildcards.level}): should never happen!"
        )
        os._exit(1)
    return paths
