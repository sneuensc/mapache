##########################################################################################
## FASTQ LEVEL

## get the path to the fastq file given SM, LB and ID
def get_fastq_of_ID(wc):
    # print(f"get_fastq_of_ID: {wc}")
    if "_R1" == wc.ID[-3:]:
        filename = samples[wc.SM][wc.LB][wc.ID[:-3]]["Data1"]
    elif "_R2" == wc.ID[-3:]:
        filename = samples[wc.SM][wc.LB][wc.ID[:-3]]["Data2"]
    elif paired_end:
        # elif paired_end != 0:  ## SE library in a paired-end sample file
        filename = samples[wc.SM][wc.LB][wc.ID]["Data1"]
    else:
        filename = samples[wc.SM][wc.LB][wc.ID]["Data"]
    return filename


## get the fastq file(s) used for mapping
def get_fastq_4_mapping(wc, rem_adapt=""):
    # print(f"get_fastq_4_mapping: {wc}")
    if rem_adapt == "":  ## if not set take the global setting
        rem_adapt = run_adapter_removal

    if rem_adapt:
        folder = f"{wc.folder}/01_fastq/01_trimmed/01_adapter_removal/{wc.SM}/{wc.LB}"
        if not paired_end:
            filename = [f"{folder}/{wc.ID}.fastq.gz"]
        elif collapse:
            filename = rules.adapter_removal_collapse.output.R
        else:
            # if str(samples[wc.SM][wc.LB][wc.ID]["Data2"]) == "nan":
            data2 = recursive_get(
                [wc.SM, wc.LB, wc.ID, "Data2"],
                "nan",
                my_dict=samples,
            )
            if data2 != data2:
                # single-end files
                filename = [f"{folder}/{wc.ID}.fastq.gz"]
            else:
                # paired-end files, not collapsing
                filename = [
                    f"{folder}/{wc.ID}_R1.fastq.gz",
                    f"{folder}/{wc.ID}_R2.fastq.gz",
                ]
    else:
        folder = f"{wc.folder}/01_fastq/00_reads/01_files_orig/{wc.SM}/{wc.LB}"
        if not paired_end:
            filename = [f"{folder}/{wc.ID}.fastq.gz"]
        else:
            data2 = recursive_get(
                [wc.SM, wc.LB, wc.ID, "Data2"],
                "nan",
                my_dict=samples,
            )
            # checking a single-end file
            if data2 != data2:
                filename = [f"{folder}/{wc.ID}.fastq.gz"]
            else:
                filename = [
                    f"{folder}/{wc.ID}_R1.fastq.gz",
                    f"{folder}/{wc.ID}_R2.fastq.gz",
                ]

    return filename


## fastqc may be run on the original or trimmed fastq files
def inputs_fastqc(wc):
    if "trim" in wc.folder:
        return get_fastq_4_mapping(wc, True)
    else:
        return get_fastq_4_mapping(wc, False)


## get the bam file used for sorting
def get_bam_4_sorting(wc):
    # print(f"get_bam_for_sorting: {wc}")
    if mapper == "bwa_aln":
        if (
            not collapse
            # and str(samples[wc.SM][wc.LB][wc.ID]["Data2"]) != "nan"
            and str(
                recursive_get(
                    [wc.SM, wc.LB, wc.ID, "Data2"],
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
    return f"{wc.folder}/01_fastq/02_mapped/{folder}/{wc.SM}/{wc.LB}/{wc.ID}.{wc.GENOME}.bam"


## get the final bam file of at the FASTQ level
def get_final_bam_FASTQ(wc):
    if run_filtering:
        file = f"{wc.folder}/01_fastq/03_filtered/01_bam_filter/{wc.SM}/{wc.LB}/{wc.ID}.{wc.GENOME}.bam"
    else:
        file = f"{folder}/01_fastq/02_mapped/03_bam_sort/{wc.SM}/{wc.LB}/{wc.ID}.{wc.GENOME}.bam"
    # print(f"get_bam_4_final_fastq: {file}")
    return file


def get_final_bam_low_qual_FASTQ(wc):
    return f"{wc.folder}/01_fastq/03_filtered/01_bam_filter_low_qual/{wc.SM}/{wc.LB}/{wc.ID}.{wc.GENOME}.bam"


##########################################################################################
##########################################################################################
## LIBRARY LEVEL

## get the bam file(s) to be merged
def get_bam_4_merge_bam_fastq2library(wc):
    return [
        get_final_bam_FASTQ(Wildcards(wc, {"ID": ID})) for ID in samples[wc.SM][wc.LB]
    ]


def get_bam_4_merge_bam_low_qual_fastq2library(wc):
    return [
        get_final_bam_low_qual_FASTQ(Wildcards(wc, {"ID": ID}))
        for ID in samples[wc.SM][wc.LB]
    ]


## get the (merged) bam file
def get_merged_bam_LB(wc):
    bam = get_bam_4_merge_bam_fastq2library(wc)
    if (
        len(bam) > 1
    ):  ## library consits of more than one fastq file: return 00_merged_fastq
        return f"{wc.folder}/02_library/00_merged_fastq/01_bam/{wc.SM}/{wc.LB}.{wc.GENOME}.bam"
    else:  ## library consits of one fastq file: return return the location of the final library bam file
        return bam[0]


def get_merged_bam_low_qual_LB(wc):
    bam = get_bam_4_merge_bam_low_qual_fastq2library(wc)
    if (
        len(bam) > 1
    ):  ## library consits of more than one fastq file: return 00_merged_fastq
        return f"{wc.folder}/02_library/00_merged_fastq/01_bam_low_qual/{wc.SM}/{wc.LB}.{wc.GENOME}.bam"
    else:  ## library consits of one fastq file: return return the location of the final library bam file
        return bam[0]


## get the bam file used to remove duplicates
def get_bam_4_markduplicates(wc):
    return get_merged_bam_LB(wc)


## get the bam file used to compte the damage
def get_bam_4_damage(wc):
    if remove_duplicates == "markduplicates":
        bam = f"{wc.folder}/02_library/01_duplicated/01_markduplicates/{wc.SM}/{wc.LB}.{wc.GENOME}.bam"
    elif remove_duplicates == "dedup":
        bam = f"{wc.folder}/02_library/01_duplicated/01_dedup/{wc.SM}/{wc.LB}.{wc.GENOME}.bam"
    else:
        bam = get_bam_4_markduplicates(wc)
    return bam


## get the final bam files at the library level
def get_final_bam_LB(wc):
    if run_damage_rescale:
        file = f"{wc.folder}/02_library/02_rescaled/01_mapDamage/{wc.SM}/{wc.LB}.{wc.GENOME}.bam"
    else:
        file = get_bam_4_damage(wc)
    return file


def get_final_bam_low_qual_LB(wc):
    return get_merged_bam_low_qual_LB(wc)


##########################################################################################
##########################################################################################
##########################################################################################
## SAMPLE LEVEL

## get the bam file(s) to be merged
def get_bam_4_merge_bam_library2sample(wc):
    return [get_final_bam_LB(Wildcards(wc, {"LB": LB})) for LB in samples[wc.SM]]


def get_bam_4_merge_bam_low_qual_library2sample(wc):
    return [
        get_final_bam_low_qual_LB(Wildcards(wc, {"LB": LB})) for LB in samples[wc.SM]
    ]


## get the (merged) bam file
def get_merged_bam_SM(wc):
    bam = get_bam_4_merge_bam_library2sample(wc)
    if (
        len(bam) > 1
    ):  ## sample consits of more than one library return 00_merged_library
        return f"{wc.folder}/03_sample/00_merged_library/01_bam/{wc.SM}.{wc.GENOME}.bam"
    else:  ## library consits of one fastq file: return return the location of the final library bam file
        return bam[0]


def get_merged_bam_low_qual_SM(wc):
    # print(f"get_merged_bam_low_qual_SM: {wc}")
    bam = get_bam_4_merge_bam_low_qual_library2sample(wc)
    if (
        len(bam) > 1
    ):  ## sample consits of more than one library return 00_merged_library
        return f"{wc.folder}/03_sample/00_merged_library/01_bam_low_qual/{wc.SM}.{wc.GENOME}.bam"
    else:  ## library consits of one fastq file: return return the location of the final library bam file
        return bam[0]


## get the bam file used to realign indels
def get_bam_4_realign(wc):
    return get_merged_bam_SM(wc)


## get the bam file used to recompute the md flag
def get_bam_4_samtools_calmd(wc):
    if run_realign:
        return f"{wc.folder}/03_sample/01_realigned/01_realign/{wc.SM}.{wc.GENOME}.bam"
    else:
        return get_bam_4_realign(wc)


## get the final bam files at the sample level
def get_bam_4_final_bam(wc):
    if run_compute_md:
        return f"{wc.folder}/03_sample/02_md_flag/01_md_flag/{wc.SM}.{wc.GENOME}.bam"
    else:
        return get_bam_4_samtools_calmd(wc)


def get_bam_4_final_bam_low_qual(wc):
    return get_merged_bam_low_qual_SM(wc)


##########################################################################################
##########################################################################################
##########################################################################################
## STATS

## get all files generated by the deamintion inference
def get_damage_output():
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
    else:  ## if False
        files = []
    return files


## get all individual stat table files to concatenate
def path_stats_by_level(wc):
    if wc.level == "FASTQ":
        paths = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.GENOME}/{SM}/{LB}/{ID}/fastq_stats.csv"
            for SM in samples
            for LB in samples[SM]
            for ID in samples[SM][LB]
        ]
    elif wc.level == "LB":
        paths = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.GENOME}/{SM}/{LB}/library_stats.csv"
            for SM in samples
            for LB in samples[SM]
        ]
    elif wc.level == "SM":
        paths = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.GENOME}/{SM}/sample_stats.csv"
            for SM in samples
        ]
    else:
        LOGGER.error(
            f"ERROR: def path_stats_by_level({wc.level}): should never happen!"
        )
        os._exit(1)
    return paths


## get the corresponding bam file (generic, across all levels):
##  - if final bam file: get the corresponding final bam file at the different levels
##  - otherwise retake the path
def get_bam_file(wc):
    paths = list(pathlib.Path(wc.file).parts)
    if paths[1] == "04_final_fastq":  ## fastq
        file = get_final_bam_FASTQ(
            Wildcards(
                fromdict={
                    "folder": wc.folder,
                    "SM": paths[3],
                    "LB": paths[4],
                    "ID": paths[5],
                    "GENOME": wc.GENOME,
                }
            )
        )
    elif paths[1] == "03_final_library":  ## library
        file = get_final_bam_LB(
            Wildcards(
                fromdict={
                    "folder": wc.folder,
                    "SM": paths[3],
                    "LB": paths[4],
                    "GENOME": wc.GENOME,
                }
            )
        )
    else:  ## the path is the right one (symlink)
        file = f"{wc.folder}/{wc.file}.{wc.GENOME}.bam"
    # print(f"get_bam_file: {file}")
    return file


## sex may be infered at teh sample or/and library level
def get_sex_file(wc, level):
    if level == "SM":
        folder = f"{wc.folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{wc.SM}.{wc.GENOME}"
    else:
        folder = f"{wc.folder}/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{wc.SM}/{wc.LB}.{wc.GENOME}"

    if str2bool(
        recursive_get_and_test(
            ["genome", wc.GENOME, "sex_inference", "run"], ["False", "True"]
        )
    ):
        return f"{folder}_sex.txt"
    else:
        return f"{folder}_nosex.txt"


##########################################################################################
## stats

## return the sample names for the stats
def get_samples_stats():
    bam_list = recursive_get(["stats", "bam_list"], "None")
    if isinstance(bam_list, dict):
        samples_stats = bam_list
    # elif bam_list == "None":
    #    #samples_stats = samples
    elif os.path.isfile(bam_list):
        db_stats = pd.read_csv(bam_list, sep=delim, comment="#", dtype=str)
        colnames = ["SM", "Data"]
        if not set(colnames).issubset(db_stats.columns):
            LOGGER.error(
                f"ERROR: The column names in the bam file (given by paramter 'stats:bam_list') are wrong! Expected are {colnames}!"
            )
            sys.exit(1)
        samples_stats = dict(zip(db_stats.SM, db_stats.Data))
    else:
        LOGGER.error(
            "ERROR: The argument of parameter 'stats:bam_list' is not valide '{bam_list}'!"
        )
        sys.exit(1)

    ## test if all bam files exist
    for file in list(samples_stats.values()):
        if not os.path.isfile(file):
            LOGGER.error(f"ERROR: Bam file '{file}' does not exist!")
            sys.exit(1)

    return samples_stats


#    samplesGLIMPSE = list(samples_files.keys())


##########################################################################################
## get all final files to run snakemake
def get_final_bam_files():
    final_bam = [
        f"{RESULT_DIR}/03_sample/03_final_sample/01_bam/{SM}.{GENOME}.bam"
        for SM in samples
        for GENOME in genome
    ]
    return final_bam


def get_final_bam_low_qual_files():
    final_bam_low_qual = [
        f"{RESULT_DIR}/03_sample/03_final_sample/01_bam_low_qual/{SM}.{GENOME}.bam"
        for SM in samples
        for GENOME in genome
        if save_low_qual
    ]
    return final_bam_low_qual


def get_stat_csv_files():
    stat_csv = [
        f"{RESULT_DIR}/04_stats/03_summary/{level}_stats.csv"
        for level in ["SM", "LB", "FASTQ"]
    ]
    return stat_csv


def get_stat_plot_files():
    plots = [
        f"{RESULT_DIR}/04_stats/04_plots/{plot_type}.svg"
        for plot_type in [
            "1_nb_reads",
            "2_mapped",
            "3_endogenous",
            "4_duplication",
            "5_AvgReadDepth",
        ]
    ]
    return plots


def get_fastqc_files():
    fastqc = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/{SM}/{LB}/{ID}_fastqc.zip"
        for SM in samples
        for LB in samples[SM]
        for ID in samples[SM][LB]
    ] + [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_adapter_removal/{SM}/{LB}/{ID}_fastqc.zip"
        for SM in samples
        for LB in samples[SM]
        for ID in samples[SM][LB]
        if run_adapter_removal
    ]
    return fastqc


def get_samtools_stats_files():
    samtools_stats = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/{file}_stats.txt"
        for GENOME in genome
        for SM in samples
        for LB in samples[SM]
        for ID in samples[SM][LB]
        for file in [
            f"01_fastq/04_final_fastq/01_bam/{SM}/{LB}/{ID}.{GENOME}",
            f"02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}",
            f"03_sample/03_final_sample/01_bam/{SM}.{GENOME}",
        ]
    ]
    return list(set(samtools_stats))  ## remove duplicates


def get_length_files():
    lengths = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/{file}_length.txt"
        for GENOME in genome
        for SM in samples
        for LB in samples[SM]
        for ID in samples[SM][LB]
        for file in [
            f"01_fastq/04_final_fastq/01_bam/{SM}/{LB}/{ID}.{GENOME}",
            f"02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}",
            f"03_sample/03_final_sample/01_bam/{SM}.{GENOME}",
        ]
    ]
    return list(set(lengths))  ## remove duplicates


def get_length_files2(samples_):
    lengths = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{GENOME}_length.txt"
        for GENOME in genome
        for SM in samples_
    ]
    return list(set(lengths))  ## remove duplicates


def get_idxstats_files():
    idxstats = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/{file}_idxstats.txt"
        for GENOME in genome
        for SM in samples
        for LB in samples[SM]
        for file in [
            f"02_library/03_final_library/01_bam/{SM}/{LB}.{GENOME}",
            f"03_sample/03_final_sample/01_bam/{SM}.{GENOME}",
        ]
    ]
    return list(set(idxstats))  ## remove duplicates


def get_qualimap_files(samples_):
    qualimap_files = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{GENOME}_qualimap"
        for GENOME in genome
        for SM in samples_
        if run_qualimap
    ]
    return qualimap_files


def get_sex_files(samples_):
    sex_files = []
    for GENOME in genome:
        ext = (
            "sex"
            if str2bool(
                recursive_get(["genome", GENOME, "sex_inference", "run"], False)
            )
            else "nosex"
        )

        sex_files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{SM}.{GENOME}_{ext}.txt"
            for SM in samples_
        ]
    return list(set(sex_files))  ## remove duplicates


def get_multiqc_files():
    multiqc_files = [
        f"{RESULT_DIR}/04_stats/02_separate_tables/{GENOME}/multiqc_fastqc.html"
        for GENOME in genome
        if run_multiqc
    ]
    return multiqc_files


def get_imputation_files(samples_):
    files = [
        f"{RESULT_DIR}/03_sample/04_imputed/{SM}.GP{GP}.phased.{ext}"
        for SM in samples_
        for GP in str2list(recursive_get(["stats", "imputation", "gp_filter"], "[0.8]"))
        for ext in ["bcf", "bcf.csi"]
        if run_imputation
    ] + [
        f"{RESULT_DIR}/03_sample/04_imputed/06_GLIMPSE_ligated/{SM}_gp.txt"
        for SM in samples_
        if run_imputation
    ]
    return files


def get_imputation_plots(samples_):
    files = [
        f"{RESULT_DIR}/03_sample/04_imputed/06_GLIMPSE_ligated/{SM}_gp.svg"
        for SM in samples_
        if run_imputation
    ]
    return files
