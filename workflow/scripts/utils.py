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
            if len(get_param(["chromosome", genome], [])) > 1:
                LOGGER.info(f"  - {genome}:")
            else:
                LOGGER.info(f"  - {genome}")
            name = get_param(["chromosome", genome, "name"], "")
            if name:
                all_sex_chr = get_param(["chromosome", genome, "all_sex_chr"], "")
                LOGGER.info(
                    f"    - Detected genome as '{name}': Applying default chromosome names for sex and mt chromosomes"
                )
            sex_chr = get_param(["chromosome", genome, "sex_chr"], "")
            if sex_chr:
                LOGGER.info(f"    - Sex chromosome (X): {to_list(sex_chr)}")
            autosomes = get_param(["chromosome", genome, "autosomes"], "")
            if autosomes:
                if len(autosomes) > 10:
                    LOGGER.info(
                        f"    - Autosomes: {autosomes[:4] + ['...'] + autosomes[-4:]}"
                    )
                else:
                    LOGGER.info(f"    - Autosomes: {autosomes}")
        elif i == 4:
            LOGGER.info(f"    - ...")
        else:
            break

    # ----------------------------------------------------------------------------------------------------------
    ## sample file
    LOGGER.info(f"SAMPLES:")
    if len(SAMPLES):
        if PAIRED_END:
            if SAMPLE_FILE == "yaml":
                LOGGER.info(f"  - Samples (YAML input) in paired-end format:")
            else:
                LOGGER.info(f"  - Sample file ('{SAMPLE_FILE}') in paired-end format:")

            LOGGER.info(f"    - {len(SAMPLES_FINAL)} SAMPLES")
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
            if SAMPLE_FILE == "yaml":
                LOGGER.info(f"  - Samples (YAML input) in single-end format:")
            else:
                LOGGER.info(f"  - Sample file ('{SAMPLE_FILE}') in single-end format:")
            LOGGER.info(f"    - {len(SAMPLES)} SAMPLES")
            LOGGER.info(f"    - {len([l for s in SAMPLES.values() for l in s])} libraries")
            LOGGER.info(
                f"    - {len([i for s in SAMPLES.values() for l in s.values() for i in l])} single-end fastq files"
            )

    if len(EXTERNAL_SAMPLES):
        if EXTERNAL_SAMPLE_FILE == "yaml":
            LOGGER.info(f"  - External samples (YAML input; only stats are computed):")
        else:
            LOGGER.info(
                f"  - External sample file ('{EXTERNAL_SAMPLE_FILE}'; only stats are computed):"
            )

        for i, genome in enumerate(EXTERNAL_SAMPLES):
            if i < 4:
                LOGGER.info(
                    f"    - {len(EXTERNAL_SAMPLES[genome])} bam files {to_list(genome)}"
                )
            elif i == 4:
                LOGGER.info(f"    - ...")
            else:
                break

    # ----------------------------------------------------------------------------------------------------------
    if len(final_bam):
        # print a summary of the workflow
        LOGGER.info("MAPPING WORKFLOW:")
        if run_subsampling:
            if type(subsampling_number) is str:
                LOGGER.info(f"  - Subsampling with group specific settings")
            elif subsampling_number < 1:
                LOGGER.info(f"  - Subsampling {100 * subsampling_number}% of the reads")
            else:
                LOGGER.info(f"  - Subsampling {subsampling_number} reads per fastq file")

        run_adapter_removal = get_param(["cleaning", "run"], ["adapterremoval", "fastp", "False"]) ## get param 'simple'
        if run_adapter_removal != False:
            if type(run_adapter_removal) is dict:
                if COLLAPSE:
                    LOGGER.info(f"  - Removing adapters is variable and collapsing paired-end reads")
                else:
                    LOGGER.info(f"  - Removing adapters is variable")
            elif run_adapter_removal == 'adapterremoval':
                if COLLAPSE:
                    LOGGER.info(f"  - Removing adapters with AdapterRemoval and collapsing paired-end reads")
                else:
                    LOGGER.info(f"  - Removing adapters with AdapterRemoval")
            else:
                LOGGER.info(f"  - Removing adapters with fastp")

        LOGGER.info(f"  - Mapping with {mapper}")

        LOGGER.info(f"  - Sorting bam file")

        if run_filtering:
            if save_low_qual:
                LOGGER.info(
                    f"  - Filtering and keeping separately low quality/unmapped reads"
                )
            else:
                LOGGER.info(f"  - Filtering and discarding low quality/unmapped reads")

        rmduplicates = get_param(["remove_duplicates", "run"], "markduplicates") ## get param 'simple'
        if rmduplicates == "markduplicates":
            LOGGER.info(f"  - Removing duplicates with MarkDuplicates")
        elif rmduplicates == "dedup":
            LOGGER.info(f"  - Removing duplicates with DeDup")
            

        if run_damage_rescale:
            LOGGER.info(f"  - Rescaling damage with MapDamage2")

        if run_realign:
            LOGGER.info(f"  - Realigning indels with GATK v3.8")

        if run_compute_md:
            LOGGER.info(f"  - Recomputing md flag")

    # ----------------------------------------------------------------------------------------------------------
    final_bam_tmp = final_bam + final_external_bam
    if len(final_bam_tmp):
        LOGGER.info("FINAL BAM FILES:")
        LOGGER.info(
            f"  - BAM FILES ({len(final_bam_tmp)} files in folder '{os.path.dirname(final_bam_tmp[0])}'):"
        )
        for i, file in enumerate(final_bam_tmp):
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
    if len(run_sex_inference):
        if len(genome):
            LOGGER.info(f"  - Inferring sex {run_sex_inference}")
        else:
            LOGGER.info(f"  - Inferring sex")

    if len(run_depth):
        if len(genome):
            LOGGER.info(f"  - Computing depth per given chromosomes {run_depth}")
        else:
            LOGGER.info(f"  - Computing depth per given chromosomes ")
    if run_damage  != False:
        if run_damage == "mapDamage":
            LOGGER.info(
                f"  - Inferring damage and read length with MapDamage2 on all alignments"
            )

        fraction = get_param(["damage", "bamdamage_fraction"], 0)
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

    if len(run_imputation):
        if len(genome) > 1:
            LOGGER.info(f"  - Imputing {run_imputation}")
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
## helper function to deal with the SAMPLES nested dict


## get all values of the given key/column (multiple if it is at the lb or sm level)
def get_sample_column(col, wc=None):
    if wc is None:
        grp = [
            id[col]
            for sm in SAMPLES.values()
            for lb in sm.values()
            for id in lb.values()
            if col in id
        ]
    elif "id" in wc.keys():
        grp = [val for keys, val in SAMPLES[wc.sm][wc.lb][wc.id].items() if col in keys]
    elif "lb" in wc.keys():
        grp = [val[col] for val in SAMPLES[wc.sm][wc.lb].values() if col in val]
    elif "sm" in wc.keys():
        grp = [id[col] for lb in SAMPLES[wc.sm].values() for id in lb.values() if col in id]

    return grp


## get the number of fastq files of the given sample path
def get_sample_count(wc=None):
    if wc is None:
        grpNb = len(
            [id for sm in SAMPLES.values() for lb in sm.values() for id in lb.values()]
        )
    elif "id" in wc.keys():
        grpNb = 1
    elif "lb" in wc.keys():
        grpNb = len([val for val in SAMPLES[wc.sm][wc.lb].values()])
    elif "sm" in wc.keys():
        grpNb = len([id for lb in SAMPLES[wc.sm].values() for id in lb.values()])
    return grpNb


## get the path, e.g., 'SAMPLES[sm][lb][id]'
def get_sample_path(wc=None):
    if wc is None:
        return "SAMPLES"
    if "id" in wc.keys():
        return f"SAMPLES[{wc.sm}][{wc.lb}][{wc.id}]"
    if "lb" in wc.keys():
        return f"SAMPLES[{wc.sm}][{wc.lb}]"
    if "sm" in wc.keys():
        return f"SAMPLES[{wc.sm}]"


## return true if the data is paired-end and throw an error if it is a mixture
def is_paired_end(wc):
    ## get all Data2 elements
    data2 = get_sample_column("Data2", wc)

    ## if not present it is SE
    if len(data2) == 0:
        return False

    ## if it is not homogenous across the group, return an error
    nbItems = get_sample_count(wc)
    if len(data2) != nbItems:
        LOGGER.error(
            f"ERROR: {get_sample_path(wc)} does not contain key 'Data2' in all rows ({len(data2)} of {nbItems} present)!"
        )
        sys.exit(1)

    ## entry may be empty (SE): either all or none
    nbEmpty = data2.count("NULL")
    if nbEmpty == 0:
        return True

    if nbEmpty == nbItems:
        return False  ## all are empty

    LOGGER.error(
        f"ERROR: {get_sample_path(wc)} has a mixture of single-end ({nbEmpty}) and paired-end ({nbItems - nbEmpty}) fastq files!"
    )
    sys.exit(1)


## return true if the data is collapsed paired-end data and throw an error if it is a mixture
def is_collapse(wc):
    ## do we have paired reads at the start?
    paired = is_paired_end(wc)
    if not paired:
        return False

    ## does the cleaning collapses the paired reads?: check if collapse is mentioned in the param
    run_cleaning = get_paramGrp(["cleaning", "run"], ["adapterremoval", "fastp", "False"], wc)
    if run_cleaning == "adapterremoval":
        params = get_paramGrp(
            ["cleaning", "params_adapterremoval"],
            "--minlength 30 --trimns --trimqualities",
            wc,
        )
        return any(
            ext in params
            for ext in ["--collapse", "--collapse-deterministic", "--collapse-conservatively"]
        )
    elif run_cleaning == "fastp":
        params = get_paramGrp(
            ["cleaning", "params_fastp"], "--minlength 30 --trimns --trimqualities", wc
        )
        return any(ext in params for ext in ["--merge", "-m"])
    else:
        return False


#######################################################################################################################
## read sample file
def read_sample_file():
    file = get_param(["sample_file"], "")
    # print(file)
    if file == "":
        return {}, ""

    if type(file) is dict:  ## yaml input
        SAMPLES = file
        input = "yaml"

    else:  ## sample file
        ## read the file
        delim = get_param(["delim"], "\\s+")
        db = pd.read_csv(file, sep=delim, comment="#", dtype=str, keep_default_na=False)
        input = file

        ## check number of columns and column names
        colsSE = ["SM", "LB", "ID", "Data"]
        colsPE = ["SM", "LB", "ID", "Data1", "Data2"]
        if not set(colsSE).issubset(db.columns) and not set(colsPE).issubset(db.columns):
            LOGGER.error(
                f"ERROR: The column names in the sample file are wrong: single-end: {colsSE}; paired-end: {colsPE}"
            )
            sys.exit(1)

        ## if SM_FINAL does not exist, copy it from SM
        if "SM_FINAL" not in db.columns:
            db['SM_FINAL'] = db.loc[:, 'SM']

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
        cols = [e for e in list(db.columns) if e not in ("SM", "LB", "ID", "SM_FINAL")]
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

        ## create nested dict of the final sample name and the libraries (remove duplicates)
        SAMPLES_FINAL = {k: list(set(g["LB"])) for k,g in db.groupby("SM_FINAL")}

    return SAMPLES, SAMPLES_FINAL, input


def test_SAMPLES():
    ## do we have paired-end data or single-end data
    PAIRED_END = len(get_sample_column("Data1")) > 0

    ## test all fastq files
    if PAIRED_END:
        ## forward reads
        # print(get_sample_column("Data1"))
        for fq in get_sample_column("Data1"):
            if not os.path.isfile(fq) and fq[:3] != "ftp":
                LOGGER.error(f"ERROR: Fastq file '{fq}' does not exist!")
                sys.exit(1)

        ## reverse reads (may be NaN)
        for fq in get_sample_column("Data2"):
            if fq != "NULL" and not os.path.isfile(fq) and fq[:3] != "ftp":
                LOGGER.error(f"ERROR: Fastq file '{fq}' does not exist!")
                sys.exit(1)

    else:  ## single end
        for fq in get_sample_column("Data"):
            if not os.path.isfile(fq) and fq[:3] != "ftp":
                LOGGER.error(f"ERROR: Fastq file '{fq}' does not exist!")
                sys.exit(1)

    ## test if collapsing is correctly set
    COLLAPSE = "--collapse" in get_param(
        ["adapterremoval", "params"], "--minlength 30 --trimns --trimqualities"
    )

    ## do we have Group settings? test if 'default' is absent
    col = "Group"
    groups = get_sample_column(col)
    if len(groups) != 0:
        ## the word 'default' is reserved and may not be used as keyword
        grpList = [ii.split(",") for ii in groups]
        if [x.count("default") for x in grpList].count(0) != len(grpList):
            LOGGER.error(
                f"ERROR: Sample file column 'Group' contains the word 'default' which is not an allowed keyword!"
            )
            sys.exit(1)

    return PAIRED_END, COLLAPSE


##########################################################################################
## read config[external_sample]
def get_external_samples():
    external_sample = get_param(["external_sample"], "")
    if external_sample == "":
        return {}, ""

    if isinstance(external_sample, dict):
        samples_stats = external_sample
        EXTERNAL_SAMPLE_FILE = "yaml"

    elif os.path.isfile(external_sample):
        delim = get_param(["delim"], "\\s+")
        db_stats = pd.read_csv(external_sample, sep=delim, comment="#", dtype=str)
        EXTERNAL_SAMPLE_FILE = external_sample
        colnames = ["SM", "Bam", "Genome"]
        if not set(colnames).issubset(db_stats.columns):
            LOGGER.error(
                f"ERROR: The column names in the bam file (given by parameter config[external_sample]) are wrong! Expected are {colnames}!"
            )
            sys.exit(1)

        ## from dataframe to nested dict
        from collections import defaultdict

        d = defaultdict(dict)
        for row in db_stats.itertuples(index=False):
            d[row.Genome][row.SM] = row.Bam
        samples_stats = dict(d)

    else:
        LOGGER.error(
            f"ERROR: The argument of parameter config[external_sample] is not valide '{external_sample}'!"
        )
        sys.exit(1)

    ## test for each GENOMES separately
    for genome in list(samples_stats):
        ## does the GENOMES name exist?
        if genome not in GENOMES:
            LOGGER.error(
                f"ERROR: Genome name config[external_sample] does not exist ({genome})!"
            )
            sys.exit(1)

        ## test if there are duplicated sample names
        if len(list(samples_stats[genome])) != len(set(list(samples_stats[genome]))):
            LOGGER.error(
                f"ERROR: Parameter config[external_sample][{genome}] contains duplicated sample names!"
            )
            sys.exit(1)

        ## test each bam file
        for id, bam in list(samples_stats[genome].items()):
            if not os.path.isfile(bam):
                LOGGER.error(
                    f"ERROR: Bam file config[external_sample][genome][{id}][bam] does not exist ({bam})!"
                )
                sys.exit(1)

    return samples_stats, EXTERNAL_SAMPLE_FILE


##########################################################################################
## all functions for main snakemake file


## get the argument of the keys (recursively acorss teh keays)
## keys: list of keys
## def_value:
##     - the default value if the parameter is not specified
##     - if a list: possible arguments (throw error if different), first element is default value
def get_param(keys, def_value, my_dict=config):
    assert type(keys) is list

    ## check if only a default value is passed or a list of possible values (first one is default arg)
    if type(def_value) is list and len(def_value) > 1:
        def_valueI = def_value[0]
    else:
        def_valueI = def_value

    ## run recursively
    if len(keys) == 1:
        arg = my_dict.get(keys[0], def_valueI)
    else:
        arg = get_param(keys[1:], def_valueI, my_dict=my_dict.get(keys[0], {}))

    arg = eval_param(arg)

    ## check if the arg is within the list of possible values
    if type(def_value) is list and len(def_value) > 1:
        if str(arg) not in def_value:
            LOGGER.error(
                f"ERROR: The parameter config[{']['.join(keys)}] has no valid argument (currently '{arg}'; available {def_value})!"
            )
            sys.exit(1)

    return arg


## same as get_param(), but arguments may be a dict with group specific settings
## group names may be specified in the sample file and may be any word,
##   except 'default', which allows to define a default argument
def get_paramGrp(keys, def_value, wc, my_dict=config):
    if type(def_value) is list and len(def_value) > 1:
        arg = get_param(keys, def_value[0], my_dict)
    else:
        arg = get_param(keys, def_value, my_dict)

    ## if it is not a dict: return value and stop here
    if type(arg) is not dict:
        ## check if the arg is within the list of possible values
        if type(def_value) is list and len(def_value) > 1:
            if str(arg) not in def_value:
                LOGGER.error(
                    f"ERROR: The parameter config[{']['.join(keys)}] has no valid argument (currently '{arg}'; available {def_value})!"
                )
                sys.exit(1)
        param = arg

    else:
        ## get all 'Group' keywords (of the sample file of the given rank (sm, lb or id)
        col = "Group"
        if "id" in wc.keys():
            grp = [val for keys, val in SAMPLES[wc.sm][wc.lb][wc.id].items() if col in keys]
            grpName = f"SAMPLES[{wc.sm}][{wc.lb}][{wc.id}]"
        elif "lb" in wc.keys():
            grp = [val[col] for val in SAMPLES[wc.sm][wc.lb].values() if col in val]
            grpName = f"SAMPLES[{wc.sm}][{wc.lb}]"
        elif "sm" in wc.keys():
            grp = [
                id[col] for lb in SAMPLES[wc.sm].values() for id in lb.values() if col in id
            ]
            grpName = f"SAMPLES[{wc.sm}]"
        grpList = [ii.split(",") for ii in grp]
        if len(grpList) == 0:
            LOGGER.error(
                f"ERROR: The parameter config[{']['.join(keys)}] has group specific settings, but no keywords are available. Is the column '{col}' missing in the sample file!"
            )
            sys.exit(1)

        ## go through all keywords and get the correct argument
        if type(def_value) is list and len(def_value) > 1:
            param = def_value[0]  ## system default argument
        else:
            param = def_value  ## system default argument

        for keyword in arg.keys():
            ## if a default is defined take it, but continue searching
            if keyword == "default":
                param = arg[keyword]
                continue

            ## get the number of occurrences of the given group keyword
            argCount = [x.count(keyword) for x in grpList].count(0)

            ## if not present continue
            if argCount == len(grpList):
                continue

            ## if present in all rows, take it
            if argCount == 0:
                param = arg[keyword]
                break

            ## if not present in all rows throw an error
            LOGGER.error(
                f"ERROR: The parameter config[{']['.join(keys)}] has group specific settings, but the keyword '{keyword}' is not omnipresent in '{grpName}' (only {argCount} of {len(grpList)})!"
            )
            sys.exit(1)

        ## check if the arg is within the list of possible values
        if type(def_value) is list and len(def_value) > 1:
            if str(param) not in def_value:
                LOGGER.error(
                    f"ERROR: The parameter config[{']['.join(keys)}] has no valid argument (currently '{param}'; available {def_value})!"
                )
                sys.exit(1)

        # print(f"{param} <= {keys}  of  {grpList} | {grpName}")
    # print(f"{param} <= {keys}")
    return param


## same as above, but a boolean is returned
## if it has a dict (group specific setting) a True is returned
## used at teh beginning to get a global view what is used
def get_param_bool(key, def_value, my_dict=config):
    arg = get_param(
        key, def_value[0], my_dict
    )  ## search without predefined list (could be a dict...)
    if type(arg) is dict:
        return True
    if str(arg) not in def_value:
        LOGGER.error(
            f"ERROR: The parameter config[{']['.join(key)}] has no valid argument (currently '{arg}'; available {def_value})!"
        )
        sys.exit(1)
    return str2bool(arg)


def get_sex_threshold_plotting():
    thresholds = {
        genome: get_param(
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


##########################################################################################
## functions to evaluate python code if necessary


## eval single element
def eval_elem(x):
    try:
        if type(x) is str:
            return eval(x, globals=None, locals=None)
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
        # return [str(i) for i in x]
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
    return genome in EXTERNAL_SAMPLES and sample in list(EXTERNAL_SAMPLES[genome])


## convert python list to R vector
def list_to_r_vector(x):
    return 'c("' + '","'.join(x) + '")'


##########################################################################################
## get all chromosome names of the given reference genome
def get_chromosome_names(genome):
    return to_str(get_param(["chromosome", genome, "all"], ""))


## set the chromosome names from fasta and store them in config
def set_chromosome_names(genome):
    ## test if fasta is valid
    fasta = get_param(["genome", genome], "")
    if not os.path.isfile(fasta):
        LOGGER.error(f"ERROR: Reference genome config[{genome}] does not exist ({fasta})!")

    ## get all chromosome names from the reference genome
    if pathlib.Path(f"{fasta}.fai").exists():
        allChr = list(map(str, pd.read_csv(f"{fasta}.fai", header=None, sep="\t")[0].tolist()))
    elif pathlib.Path(f"{RESULT_DIR}/00_reference/{genome}/{genome}.fasta.fai").exists():
        fasta = f"{RESULT_DIR}/00_reference/{genome}/{genome}.fasta"
        allChr = list(map(str, pd.read_csv(f"{fasta}.fai", header=None, sep="\t")[0].tolist()))
    else:
        cmd = f"grep '^>' {fasta} | cut -c2- | awk '{{print $1}}'"
        allChr = list(map(str, subprocess.check_output(cmd, shell=True, text=True).split()))
    config = update_value(["chromosome", genome, "all"], allChr)


## return a list of the chromosome names which do not match
## empty list means that all matched
def valid_chromosome_names(genome, names):
    allChr = get_chromosome_names(genome)

    if type(names) is not list and names not in allChr:
        return [names]

    if list(set(names) - set(allChr)):
        return list(set(names) - set(allChr))

    return []


## check the chromosome names if they are valid (sex inference for this genome is set to true)
def set_sex_inference(genome):
    ## get the chromosome names of the given genome
    allChr = get_chromosome_names(genome)

    ## get the specified sex and autosome chromosome names
    sex_chr = to_str(get_param(["sex_inference", genome, "sex_chr"], ""))
    autosomes = to_str(get_param(["sex_inference", genome, "autosomes"], []))

    ## if the sex and autosome chromosome names are set, check if they make sense
    if sex_chr != "" and len(autosomes) > 0:
        # check if the chromosomes specified in sex determination exist
        ## X chromosome
        if len(sex_chr):
            # print(sex_chr)
            if valid_chromosome_names(genome, sex_chr):
                LOGGER.error(
                    f"ERROR: Sex chromosome specified in config[sex_inference][{genome}][sex_chr] ({sex_chr}) does not exist in the reference genome."
                )
                os._exit(1)
            config = update_value(["chromosome", genome, "sex_chr"], sex_chr)
        else:
            LOGGER.error(
                f"ERROR: No sex chromosome specified in config[sex_inference][{genome}][sex_chr]!"
            )
            os._exit(1)

        # autosomes
        if len(autosomes):
            if valid_chromosome_names(genome, autosomes):
                LOGGER.error(
                    f"ERROR: In config[sex_inference][{genome}][autosomes], the following chromosome names are not recognized: {valid_chromosome_names(genome , autosomes)}!"
                )
                os._exit(1)
            config = update_value(["chromosome", genome, "autosomes"], autosomes)
        else:
            LOGGER.error(
                f"ERROR: No autosomes specified in config[sex_inference][{genome}][autosomes]!"
            )
            os._exit(1)
    else:  ## if they are not set, try to infer the genome (hg19 or GRCh38)
        hg19 = list(map(str, list(range(1, 23)) + ["X", "Y", "MT"]))
        GRCh38 = [f"chr{x}" for x in list(range(1, 23)) + ["X", "Y", "M"]]
        if set(hg19).issubset(set(allChr)):
            name = "hg19"
            detectedSexChrom = ["X", "Y", "MT"]
            detectedAutosomes = list(set(hg19) - set(detectedSexChrom))
        elif set(GRCh38).issubset(set(allChr)):
            name = "GRCh38"
            detectedSexChrom = ["chrX", "chrY", "chrM"]
            detectedAutosomes = list(set(GRCh38) - set(detectedSexChrom))
        else:
            LOGGER.error(
                f"ERROR: For sex inference the parameters config[sex_inference][{genome}][sex_chr] and config[sex_inference][{genome}][autosomes] are required!"
            )
            os._exit(1)

        ## write the sex and autosome names
        config = update_value(["chromosome", genome, "name"], name)
        config = update_value(["chromosome", genome, "all_sex_chr"], detectedSexChrom)
        config = update_value(["chromosome", genome, "sex_chr"], detectedSexChrom[0])
        config = update_value(["chromosome", genome, "autosomes"], detectedAutosomes)


## check if the specified chromosomes to compute the depth on are available
def read_depth(genome):
    # check if chromosomes for which DoC was requested exist
    depth = to_str(get_param(["depth", genome, "chromosomes"], ""))
    if valid_chromosome_names(genome, depth):
        LOGGER.error(
            f"ERROR: config[depth][{genome}][chromosomes] contains unrecognized chromosome names ({valid_chromosome_names(genome , depth)})!"
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
    mem_start = int(get_param(moduleList + ["mem"], default))
    mem_incre = int(
        get_param(moduleList + ["mem_increment"], memory_increment_ratio * mem_start)
    )
    return int(1024 * ((attempt - 1) * mem_incre + mem_start))


## in this second version the 'mem' is added to the word of the last element
def get_memory_alloc2(module, attempt, default=2):
    moduleList = module
    if type(moduleList) is not list:
        moduleList = [moduleList]
    mem_start = int(get_param(moduleList[:-1] + [moduleList[-1] + "_mem"], default))
    mem_incre = int(
        get_param(
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
    time_start = int(get_param(moduleList + ["time"], default))
    time_incre = int(
        get_param(moduleList + ["time_increment"], runtime_increment_ratio * time_start)
    )
    return int(60 * ((attempt - 1) * time_incre + time_start))


## in this second version the 'time' is added to the word of the last element
def get_runtime_alloc2(module, attempt, default=12):
    moduleList = module
    if type(moduleList) is not list:
        moduleList = [module]
    time_start = int(get_param(moduleList[:-1] + [moduleList[-1] + "_time"], default))
    time_incre = int(
        get_param(
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
    return int(get_param(moduleList + ["threads"], default))


## in this second version the 'threads' is added to the word of the last element
def get_threads2(module, default=1):
    moduleList = module
    if type(moduleList) is not list:
        moduleList = [module]
    return int(get_param(moduleList[:-1] + [moduleList[-1] + "_threads"], default))


def bam2bai(bam):
    # return bam.replace('.bam', '.bai')
    return f"{bam[: len(bam) - 4]}.bai"


##########################################################################################
## check if java is called by a .jar file or by a wrapper


def get_gatk_bin():
    bin = get_param(["software", "gatk3_jar"], "GenomeAnalysisTK.jar")
    if bin[-4:] == ".jar":
        bin = f"java -XX:ParallelGCThreads={{snakemake.threads}} -XX:+UseParallelGC -XX:-UsePerfData \
            -Xms{{snakemake.resources.memory}}m -Xmx{{snakemake.resources.memory}}m -jar {bin}"
    return bin
