from tkinter import E
import pandas as pd
import numpy as np
import itertools
import pathlib
import re
import subprocess, os.path

from snakemake.io import Wildcards

##########################################################################################
## VERBOSE LOG
##########################################################################################
## print a summary of the contents
def write_log():
    ## reference genome
    if len(GENOMES) == 1:
        LOGGER.info(f"REFERENCE GENOME: 1 genome is specified:")
    elif len(GENOMES) < 4:
        LOGGER.info(f"REFERENCE GENOME: {len(GENOMES)} genomes are specified:")
    else:
        LOGGER.info(f"REFERENCE GENOME: {len(GENOMES)} genomes are specified.")

    ## get all chromosome names and store them in the dict config[chromosomes][genome][all] for later use
    for i, genome in enumerate(GENOMES):
        if i < 4:
            if len(recursive_get(["chromosome", genome], [])) > 1:
                LOGGER.info(f"  - {genome}:")
            else:
                LOGGER.info(f"  - {genome}")
            name = recursive_get(["chromosome", genome, "name"], "")
            if name:
                all_sex_chr = recursive_get(["chromosome", genome, "all_sex_chr"], "")
                LOGGER.info(
                    f"    - Detected genome as '{name}': Applying default chromosome names for sex and mt chromosomes"
                )
            sex_chr = recursive_get(["chromosome", genome, "sex_chr"], "")
            if sex_chr:
                LOGGER.info(f"    - Sex chromosome (X): {to_list(sex_chr)}")
            autosomes = recursive_get(["chromosome", genome, "autosomes"], "")
            if autosomes:
                if(len(autosomes)> 10):
                    LOGGER.info(f"    - Autosomes: {autosomes[:4] + ['...'] + autosomes[-4:]}")
                else:
                    LOGGER.info(f"    - Autosomes: {autosomes}")
        elif i == 4:
            LOGGER.info(f"    - ...")
        else:
            break

    # ----------------------------------------------------------------------------------------------------------
    ## sample file
    LOGGER.info(f"SAMPLES:"
        )
    if paired_end:
        LOGGER.info(
            f"  - Sample file '{sample_file}' in paired-end format:"
        )
        LOGGER.info(f"    - {len(SAMPLES)} SAMPLES")
        LOGGER.info(f"    - {len([l for s in SAMPLES.values() for l in s])} libraries")
        tmp = [
            i["Data2"] for s in SAMPLES.values() for l in s.values() for i in l.values()
        ]
        LOGGER.info(
            f"    - {len([x for x in tmp if str(x) != 'nan'])} paired-end fastq files"
        )
        nb = len([x for x in tmp if str(x) == "nan"])
        if nb:
            LOGGER.info(f"    - {nb} single-end fastq files")
    else:
        LOGGER.info(
            f"  - Sample file '{sample_file}' in single-end format:"
        )
        LOGGER.info(f"    - {len(SAMPLES)} SAMPLES")
        LOGGER.info(f"    - {len([l for s in SAMPLES.values() for l in s])} libraries")
        LOGGER.info(
            f"    - {len([i for s in SAMPLES.values() for l in s.values() for i in l])} single-end fastq files"
        )

    if len(EXTERNAL_SAMPLES) > 0:
        LOGGER.info(f"  - External samples (only stats are computed):")
        
        for i, genome in enumerate(EXTERNAL_SAMPLES):
            if i < 4:
                LOGGER.info(f"    - {len(EXTERNAL_SAMPLES[genome])} bam files {to_list(genome)}")
            elif i == 4:
                LOGGER.info(f"    - ...")
            else:
                break

    # ----------------------------------------------------------------------------------------------------------
    # print a summary of the workflow
    LOGGER.info("MAPPING WORKFLOW:")
    if run_subsampling:
        if subsampling_number < 1:
            LOGGER.info(f"  - Subsampling {100 * subsampling_number}% of the reads")
        else:
            LOGGER.info(f"  - Subsampling {subsampling_number} reads per fastq file")

    if run_adapter_removal:
        if collapse:
            LOGGER.info(
                f"  - Removing adapters with AdapterRemoval and collapsing paired-end reads"
            )
        else:
            LOGGER.info(f"  - Removing adapters with AdapterRemoval")

    LOGGER.info(f"  - Mapping with {mapper}")

    LOGGER.info(f"  - Sorting bam file")

    if run_filtering:
        if save_low_qual:
            LOGGER.info(
                f"  - Filtering and keeping separately low quality/unmapped reads"
            )
        else:
            LOGGER.info(f"  - Filtering and removing low quality/unmapped reads")

    if remove_duplicates != "False":
        if remove_duplicates == "markduplicates":
            LOGGER.info(f"  - Removing duplicates with MarkDuplicates")
        elif remove_duplicates == "dedup":
            LOGGER.info(f"  - Removing duplicates with DeDup")
        else:
            LOGGER.info(f"  - Removing duplicates with {remove_duplicates}")

    if run_damage_rescale:
        LOGGER.info(f"  - Rescaling damage with MapDamage2")

    if run_realign:
        LOGGER.info(f"  - Realigning indels with GATK v3.8")

    if run_compute_md:
        LOGGER.info(f"  - Recomputing md flag")

    # ----------------------------------------------------------------------------------------------------------
    LOGGER.info("FINAL BAM FILES:")
    LOGGER.info(
        f"  - BAM FILES ({len(final_bam)} files in folder '{os.path.dirname(final_bam[0])}'):"
    )
    for i, file in enumerate(final_bam):
        if i < 4:
            LOGGER.info(f"    - {file}")
        elif i == 4:
            LOGGER.info(f"    - ...")
        else:
            break

    if save_low_qual:
        LOGGER.info(
            f"  - LOW QUALITY AND UNMAPPED BAM FILES ({len(final_bam_low_qual)} files in folder '{os.path.dirname(final_bam_low_qual[0])}'):"
        )
        for i, file in enumerate(final_bam_low_qual):
            if i < 4:
                LOGGER.info(f"    - {file}")
            elif i == 4:
                LOGGER.info(f"    - ...")
            else:
                break

    # ----------------------------------------------------------------------------------------------------------
    LOGGER.info("ANALYSES:")
    if len(run_sex_inference) > 0:
        if len(genome) > 0:
            LOGGER.info(f"  - Inferring sex {run_sex_inference}")
        else:
            LOGGER.info(f"  - Inferring sex")

    if len(run_depth) > 0:
        if len(genome) > 0:
            LOGGER.info(
                f"  - Computing depth per given chromosomes {run_depth}"
            )
        else:
            LOGGER.info(f"  - Computing depth per given chromosomes ")

    if run_damage != "False":
        if run_damage == "mapDamage":
            LOGGER.info(
                f"  - Inferring damage and read length with MapDamage2 on all alignments"
            )

        fraction = recursive_get(["damage", "bamdamage_fraction"], 0)
        if fraction == 0:
            LOGGER.info(
                f"  - Inferring damage and read length with bamdamge on all alignments"
            )
        elif fraction < 1:
            LOGGER.info(
                f"  - Inferring damage and read length with bamdamge on {100 * fraction}% of the alignments"
            )
        else:
            LOGGER.info(
                f"  - Inferring damage and read length with bamdamge on {fraction} alignments"
            )

    if run_imputation:
        if len(genome) > 1:
            LOGGER.info(f"  - Imputing on genome {to_list(GENOMES[0])} (first genome in list)")
        else:
            LOGGER.info(f"  - Imputing")

    # ----------------------------------------------------------------------------------------------------------
    # print a summary of the stats
    LOGGER.info("REPORTS:")
    LOGGER.info(f"  - Mapping statistics tables:")
    for i, file in enumerate(stat_csv):
        LOGGER.info(f"    - {file}")

    if run_qualimap:
        LOGGER.info(
            f"  - Qualimap report: {len(qualimap_files)} report(s) on final bam files:"
        )
        for i, file in enumerate(qualimap_files):
            if i < 4:
                LOGGER.info(f"    - {file}")
            elif i == 4:
                LOGGER.info(f"    - ...")
            else:
                break

    if run_multiqc:
        LOGGER.info(
            f"  - Mulitqc report: {len(multiqc_files)} HTML report(s) combining statistics per genome:"
        )
        for i, file in enumerate(multiqc_files):
            if i < 4:
                LOGGER.info(f"    - {file}")
            elif i == 4:
                LOGGER.info(f"    - ...")
            else:
                break

    LOGGER.info(
        f"  => Snakemake report: run 'snakemake --report report.html' after finalization of the run to get the report"
    )


##########################################################################################################
def read_sample_file(file):
    ## if file does not exist or is not specified
    if not os.path.isfile(file):
        LOGGER.warning(f"WARNING: No sample file (config[sample_file]) is specified!")
        return {}, False, False

    ## read the file
    delim = recursive_get(["delim"], "\s+")
    db = pd.read_csv(file, sep=delim, comment="#", dtype=str)

    ## check number of columns and column names
    if len(db.columns) == 6:
        paired_end = False
        collapse = False
        colnames = ["ID", "Data", "MAPQ", "LB", "PL", "SM"]
        if not set(colnames).issubset(db.columns):
            LOGGER.error(
                f"ERROR: The column names in the sample file are wrong! Expected are for single-end reads {colnames}!"
            )
            sys.exit(1)

        ## test if all fastq files exist
        for fq in db["Data"].tolist():
            if not os.path.isfile(fq):
                LOGGER.error(f"ERROR: Fastq file '{fq}' does not exist!")
                sys.exit(1)

    elif len(db.columns) == 7:
        paired_end = True
        adaptrem_params = recursive_get(
            ["adapterremoval", "params"], "--minlength 30 --trimns --trimqualities"
        )
        if "--collapse" in adaptrem_params:
            collapse = True
        else:
            collapse = False
        colnames = ["ID", "Data1", "Data2", "MAPQ", "LB", "PL", "SM"]
        if not set(colnames).issubset(db.columns):
            LOGGER.error(
                f"ERROR: The column names in the sample file are wrong! Expected are for paired-end reads {colnames}!"
            )
            sys.exit(1)

        ## test if all fastq files exist
        for fq in db["Data1"].tolist():
            if fq != fq:  # test if NaN
                LOGGER.error(f"ERROR: Fastq file '{fq}' in column Data1 is missing!")
                sys.exit(1)
            if not fq and not os.path.isfile(fq):
                LOGGER.error(f"ERROR: Fastq file '{fq}' does not exist!")
                sys.exit(1)

        ## test if all fastq files exist or are NaN
        for fq in db["Data2"].tolist():
            if fq == fq and not os.path.isfile(fq):
                LOGGER.error(
                    f"ERROR: Fastq file '{fq}' in column Data2 does not exist!"
                )
                sys.exit(1)
    else:
        LOGGER.error(
            f"ERROR: The number of columns in the sample file is wrong ({len(db.columns)} columns)!"
        )
        sys.exit(1)
    ## test if SM, LB and ID names are valid
    # invalid_char = string.punctuation.replace("_", "").replace("-", "")
    # if bool(re.match("^[a-zA-Z0-9]*$", string)) == True:
    #    print("String does not contain any special characters.")
    # else:
    #    print("The string contains special characters.")

    ## --------------------------------------------------------------------------------------------------
    ## check if all IDs per LB and SM are unique
    all_fastq = db.groupby(["ID", "LB", "SM"])["ID"].agg(["count"]).reset_index()
    if max(all_fastq["count"]) > 1:
        LOGGER.error(
            f"ERROR: The ID's {all_fastq['ID']} are not uniq within library and sample!"
        )
        sys.exit(1)

    ## --------------------------------------------------------------------------------------------------
    ## dataframe to nested dict
    SAMPLES = {}
    cols = [e for e in list(db.columns) if e not in ("SM", "LB", "ID")]
    for index, row in db.iterrows():
        ## if key not present add new dict
        if row["SM"] not in SAMPLES:
            SAMPLES[row["SM"]] = {}
        if row["LB"] not in SAMPLES[row["SM"]]:
            SAMPLES[row["SM"]][row["LB"]] = {}
        if row["ID"] not in SAMPLES[row["SM"]][row["LB"]]:
            SAMPLES[row["SM"]][row["LB"]][row["ID"]] = {}

        ## add all remaining columns to this dict
        for col in cols:
            SAMPLES[row["SM"]][row["LB"]][row["ID"]][col] = row[col]

    return SAMPLES, paired_end, collapse


##########################################################################################
## all functions for main snakemake file

## get recursively the argument of the given keys in a nested dict (run eval())
def recursive_get(keys, def_value, my_dict=config):
    assert type(keys) is list
    first = keys[0]
    if len(keys) == 1:
        value = my_dict.get(first, def_value)
    else:
        value = recursive_get(keys[1:], def_value, my_dict=my_dict.get(first, {}))
    return eval_param(value)



## same as above, but the argument is tested if it is present in the 'available_args'
## first argument of 'available_args' is the default
def recursive_get_and_test(key, available_args, my_dict=config):
    arg = str(recursive_get(key, available_args[0], my_dict))
    if arg not in available_args:
        LOGGER.error(
            f"ERROR: The parameter config[{']['.join(key)}] has no valid argument (currently '{arg}'; available {available_args})!"
        )
        sys.exit(1)
    return arg


## update the argument in a nested dict (and create leaves if needed)
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
        genome: recursive_get(
            keys=["sex_inference", genome, "thresholds"],
            def_value='list( "XX"=c(0.8, 1), "XY"=c(0, 0.6), "consistent with XX but not XY"=c(0.6, 1), "consistent with XY but not XX"=c(0, 0.8) )',
        )
        for genome in GENOMES
    }

    sex_thresholds = "list({pair})".format(
        pair=",     ".join([f'"{genome}"={thresholds[genome]}' for genome in GENOMES])
    )
    sex_thresholds = sex_thresholds.replace("=", "?")
    return sex_thresholds


##########################################################################################
## functions to evaluate python code if necessary

## eval single element
def eval_elem(x):
    try:
        if type(x) is str:
            return eval(x)
        else:
            return x
    except:
        return x

def eval_param(x):
    if type(x) is list:
        return list(map(eval_elem, x))
    elif type(x) is str:
        return eval_elem(x)
    else:
        return x


def to_str(x):
    if type(x) is list:
        return list(map(str, x))
    elif type(x) is int:
        return str(x)        
    elif type(x) is float:
        return str(x)
    else:
        return x


def to_list(x):
    if type(x) is list:
        return x    
    return list(to_str(x).split(" "))

   
def is_external_sample(sample, genome):
    return ( genome in EXTERNAL_SAMPLES and sample in list(EXTERNAL_SAMPLES[genome]))

## eval a list with potential eval elements
#def eval_list(x):
##    return list(map(str, eval_to_list(x)))


## eval single element if needed and return it as a list
##def eval_to_list(x):
#    a = eval_if_possible(x)
#    if isinstance(a, list):
#        return a
#    else:
#        return [a]


## eval single element if needed and return it as a comma separated string
#def to_csv(x):
#    return ",".join(list(map(str, eval_to_list(x))))


## transform a list to a comma separated string
#def list_to_csv(x):
#    if isinstance(x, list):
#        return ",".join(x)
#    else:
#        return x


## replace any element of the list
#def eval_list(x):
#    if isinstance(x, list):
#        return list(map(eval_if_possible, x))
#    else:
#        return [eval_if_possible(x)]


## replace any element of the list
#def eval_list_to_csv(x):
#    return ",".join(eval_list(x))


## convert python list to R vector
def list_to_r_vector(x):
    return 'c("' + '","'.join(x) + '")'


##########################################################################################
## get all chromosome names of the given reference genome
def get_chromosome_names(genome):
    return to_str(recursive_get(["chromosome", genome, "all"], ""))


## set the chromosome names from fasta and store them in config
def set_chromosome_names(genome):
    ## test if fasta is valid
    fasta = recursive_get(["genome", genome], "")
    if not os.path.isfile(fasta):
        LOGGER.error(
            f"ERROR: Reference genome config[{genome}] does not exist ({fasta})!"
        )

    ## get all chromosome names from the reference genome
    if pathlib.Path(f"{fasta}.fai").exists():
        allChr = list(
            map(str, pd.read_csv(f"{fasta}.fai", header=None, sep="\t")[0].tolist())
        )
    elif pathlib.Path(
        f"{RESULT_DIR}/00_reference/{genome}/{genome}.fasta.fai"
    ).exists():
        fasta = f"{RESULT_DIR}/00_reference/{genome}/{genome}.fasta"
        allChr = list(
            map(str, pd.read_csv(f"{fasta}.fai", header=None, sep="\t")[0].tolist())
        )
    else:
        cmd = f"grep '^>' {fasta} | cut -c2- | awk '{{print $1}}'"
        allChr = list(
            map(str, subprocess.check_output(cmd, shell=True, text=True).split())
        )
    config = update_value(["chromosome", genome, "all"], allChr)


## return a list of the chromosome names which do not match
def valid_chromosome_names(genome, names):
    allChr = to_str(recursive_get(["chromosome", genome, "all"], []))
    if list(set(names) - set(allChr)):
        return list(set(names) - set(allChr))
    else:
        return []


## check the chromosome names if they are valid
def set_sex_inference(genome):
    ## get the chromosome names of the given genome
    allChr = get_chromosome_names(genome)

    ## try to recognize if it is the human hg19 or GRCh38 genome. If so apply default chromosome names
    hg19 = list(map(str, list(range(1, 23)) + ["X", "Y", "MT"]))
    GRCh38 = [f"str{x}" for x in list(range(1, 23)) + ["X", "Y", "M"]]
    if sorted(allChr) == sorted(hg19):
        name = "hg19"
        detectedSexChrom = ["X", "Y", "MT"]
        detectedAutosomes = list(set(allChr) - set(detectedSexChrom))
    elif sorted(allChr) == sorted(GRCh38):
        name = "GRCh38"
        detectedSexChrom = ["chrX", "chrY", "chrM"]
        detectedAutosomes = list(set(allChr) - set(detectedSexChrom))
    else:
        name = ""
        detectedSexChrom = []
        detectedAutosomes = []

    ## write the sex and autosome names
    config = update_value(["chromosome", genome, "name"], name)
    config = update_value(["chromosome", genome, "all_sex_chr"], detectedSexChrom)
    config = update_value(["chromosome", genome, "sex_chr"], detectedSexChrom[0])
    config = update_value(["chromosome", genome, "autosomes"], detectedAutosomes)

    # check if the chromosomes specified in sex determination exist
    ## X chromosome
    sex_chr = recursive_get(["sex_inference", genome, "sex_chr"], [])
    if len(sex_chr):
        if valid_chromosome_names(genome, sex_chr):
            LOGGER.error(
                f"ERROR: Sex chromosome specified in config[sex_inference][{genome}][sex_chr] ({sex_chr}) does not exist in the reference genome."
            )
            os._exit(1)
        config = update_value(["chromosome", genome, "sex_chr"], sex_chr)

    # autosomes
    autosomes = to_str(recursive_get(
                    ["sex_inference", genome, "autosomes"],
                    [],
                )
        )

    if len(autosomes):
        if valid_chromosome_names(genome, autosomes):
            LOGGER.error(
                f"ERROR: In config[sex_inference][{genome}][autosomes], the following chromosome names are not recognized: {valid_chromosome_names(genome, autosomes)}!"
            )
            os._exit(1)
        config = update_value(["chromosome", genome, "autosomes"], autosomes)


def read_depth(genome):
    # check if chromosomes for which DoC was requested exist
    allChr = recursive_get(["chromosome", genome, "all"], "")
    depth = recursive_get(["depth", genome, "chromosomes"], "")
    if valid_chromosome_names(genome, depth):
        LOGGER.error(
            f"ERROR: config[depth][{genome}][chromosomes] contains unrecognized chromosome names ({valid_chromosome_names(genome, depth)})!"
        )
        os._exit(1)
    config = update_value(["chromosome", genome, "depth"], depth)


## convert string to boolean
def str2bool(v):
    return str(v).lower() in ("yes", "true", "t", "1")


## convert any argument to a list of string(s)
def str2list(v):
    if type(v) is list:
        return [str(x) for x in v]
    else:
        return [str(v)]


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


## in this second version the 'mem' is added to the word of the last element
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


## in this second version the 'time' is added to the word of the last element
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


## get the number of threads of the given parameter
def get_threads(module, default=1):
    moduleList = module
    if type(moduleList) is not list:
        moduleList = [module]
    return int(recursive_get(moduleList + ["threads"], default))


## in this second version the 'threads' is added to the word of the last element
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
