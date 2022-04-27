import pandas as pd
import numpy as np
import itertools
import pathlib
import re

from snakemake.io import Wildcards

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


def get_sex_threshold_plotting():

    thresholds = {
        GENOME: recursive_get(
            keys=["genome", GENOME, "sex_inference", "params", "thresholds"],
            def_value='list( "XX"=c(0.8, 1), "XY"=c(0, 0.6), "consistent with XX but not XY"=c(0.6, 1), "consistent with XY but not XX"=c(0, 0.8) )',
        )
        for GENOME in genome
    }

    sex_thresholds = "list({pair})".format(
        pair=",     ".join([f'"{GENOME}"={thresholds[GENOME]}' for GENOME in genome])
    )
    sex_thresholds = sex_thresholds.replace("=", "?")
    return sex_thresholds


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
def check_chromosome_names(GENOME, Logging=True):
    if Logging:
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
        allChr = list(
            map(str, subprocess.check_output(cmd, shell=True, text=True).split())
        )
    else:
        LOGGER.error(f"ERROR: Reference genome 'genome:{GENOME}:fasta' does not exist!")
        os._exit(1)

    ## try to recognise if it is the human hg19 or GRCh38 genome. If so apply default chromosome names
    hg19 = map(str, list(range(1, 23)) + ["X", "Y", "MT"])
    GRCh38 = [f"chr{x}" for x in list(range(1, 23)) + ["X", "Y", "M"]]
    if sorted(allChr) == sorted(hg19):
        detectedChromosomes = ["X", "Y", "MT"]
        detectedAutosomes = list(set(allChr) - set(detectedChromosomes))
        if Logging:
            LOGGER.info(
                f"    - Detected genome as 'hg19': Appliyng default chromosome names for sex and mt chromsomes: {detectedChromosomes}."
            )

        config = update_value(
            ["genome", GENOME, "sex_inference", "params", "sex_chr"],
            detectedChromosomes[0],
        )
        config = update_value(
            ["genome", GENOME, "sex_inference", "params", "autosomes"],
            detectedAutosomes,
        )
    elif sorted(allChr) == sorted(GRCh38):
        detectedChromosomes = ["chrX", "chrY", "chrM"]
        detectedAutosomes = list(set(allChr) - set(detectedChromosomes))
        if Logging:
            LOGGER.info(
                f"    - Detected genome as 'GRCh38': Appliyng default chromosome names for sex and mt chromsomes {detectedChromosomes}."
            )

        config = update_value(
            ["genome", GENOME, "sex_inference", "params", "sex_chr"],
            detectedChromosomes[0],
        )
        config = update_value(
            ["genome", GENOME, "sex_inference", "params", "autosomes"],
            detectedAutosomes,
        )

    # check if the chromosomes specified in sex determination exist
    # sex chromosome
    if recursive_get(["genome", GENOME, "sex_inference", "run"], False):
        if Logging:
            LOGGER.info(f"    - Inferring sex")
        ## X chromosome specified for the sex inferenceß
        sex_chr = recursive_get(
            ["genome", GENOME, "sex_inference", "params", "sex_chr"], []
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
                        [],
                    )
                ),
            )
        )
        if list(set(autosomes) - set(allChr)):
            LOGGER.error(
                f"ERROR: In 'config[genome][{GENOME}][sex_inference][params][autosomes]', the following chromsome names are not recognized: {list(set(autosomes) - set(allChr))}!"
            )
            os._exit(1)

        ## for the autosomes transfer python to R format
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
    if Logging and depth_chromosomes:
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
    moduleList = module
    if type(moduleList) is not list:
        moduleList = [module]
    return int(recursive_get(moduleList + ["threads"], default))


## in this second verion the 'threads' is added to the word of the last element
def get_threads2(module, default=1):
    moduleList = module
    if type(moduleList) is not list:
        moduleList = [module]
    return int(recursive_get(moduleList[:-1] + [moduleList[-1] + "_threads"], default))


def bam2bai(bam):
    # return bam.replace('.bam', '.bai')
    return f"{bam[:len(bam) - 4]}.bai"


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


def symlink_rev(input, output):
    shell("cp {input} {output}")
# shell("ln -srf {input} {output}")
