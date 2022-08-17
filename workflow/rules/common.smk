##########################################################################################
## FASTQ LEVEL

## get the path to the fastq file given sm, lb and id
def get_fastq_of_ID(wc):
    # print(f"get_fastq_of_ID: {wc}")
    if "_R1" == wc.id[-3:]:
        filename = SAMPLES[wc.sm][wc.lb][wc.id[:-3]]["Data1"]
    elif "_R2" == wc.id[-3:]:
        filename = SAMPLES[wc.sm][wc.lb][wc.id[:-3]]["Data2"]
    elif PAIRED_END:
        # elif PAIRED_END != 0:  ## SE library in a paired-end sample file
        filename = SAMPLES[wc.sm][wc.lb][wc.id]["Data1"]
    else:
        filename = SAMPLES[wc.sm][wc.lb][wc.id]["Data"]
    return filename


## get the fastq file(s) used for mapping
def get_fastq_4_mapping(wc, rem_adapt=""):
    # print(f"get_fastq_4_mapping: {wc}")
    if rem_adapt == "":  ## if not set take the global setting
        rem_adapt = run_adapter_removal

    if rem_adapt:
        folder = f"{wc.folder}/01_fastq/01_trimmed/01_adapter_removal/{wc.sm}/{wc.lb}"
        if not PAIRED_END:
            filename = [f"{folder}/{wc.id}.fastq.gz"]
        elif COLLAPSE:
            filename = rules.adapter_removal_collapse.output.R
        else:
            # if str(SAMPLES[wc.sm][wc.lb][wc.id]["Data2"]) == "nan":
            data2 = recursive_get(
                [wc.sm, wc.lb, wc.id, "Data2"],
                "nan",
                my_dict=SAMPLES,
            )
            if data2 != data2:
                # single-end files
                filename = [f"{folder}/{wc.id}.fastq.gz"]
            else:
                # paired-end files, not collapsing
                filename = [
                    f"{folder}/{wc.id}_R1.fastq.gz",
                    f"{folder}/{wc.id}_R2.fastq.gz",
                ]
    else:
        folder = f"{wc.folder}/01_fastq/00_reads/01_files_orig/{wc.sm}/{wc.lb}"
        if not PAIRED_END:
            filename = [f"{folder}/{wc.id}.fastq.gz"]
        else:
            data2 = recursive_get(
                [wc.sm, wc.lb, wc.id, "Data2"],
                "nan",
                my_dict=SAMPLES,
            )
            # checking a single-end file
            if data2 != data2:
                filename = [f"{folder}/{wc.id}.fastq.gz"]
            else:
                filename = [
                    f"{folder}/{wc.id}_R1.fastq.gz",
                    f"{folder}/{wc.id}_R2.fastq.gz",
                ]

    return filename


## fastqc may be run on the original or trimmed fastq files
def inputs_fastqc(wc):
    if "trim" in wc.type:
        return get_fastq_4_mapping(wc, True)
    else:
        return get_fastq_4_mapping(wc, False)


## get the bam file used for sorting
def get_bam_4_sorting(wc):
    # print(f"get_bam_for_sorting: {wc}")
    if mapper == "bwa_aln":
        if (
            not COLLAPSE
            # and str(SAMPLES[wc.sm][wc.lb][wc.id]["Data2"]) != "nan"
            and str(
                recursive_get(
                    [wc.sm, wc.lb, wc.id, "Data2"],
                    "nan",
                    my_dict=SAMPLES,
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
            f"ERROR: The parameter config[mapping][mapper] is not correctly specified: {mapper} is unknown!"
        )
        os._exit(1)
    return f"{wc.folder}/01_fastq/02_mapped/{folder}/{wc.sm}/{wc.lb}/{wc.id}.{wc.genome}.bam"


## get the final bam file of at the FASTQ level
def get_final_bam_FASTQ(wc):
    if run_filtering:
        file = f"{wc.folder}/01_fastq/03_filtered/01_bam_filter/{wc.sm}/{wc.lb}/{wc.id}.{wc.genome}.bam"
    else:
        file = f"{folder}/01_fastq/02_mapped/03_bam_sort/{wc.sm}/{wc.lb}/{wc.id}.{wc.genome}.bam"
    # print(f"get_bam_4_final_fastq: {file}")
    return file


def get_final_bam_low_qual_FASTQ(wc):
    return f"{wc.folder}/01_fastq/03_filtered/01_bam_filter_low_qual/{wc.sm}/{wc.lb}/{wc.id}.{wc.genome}.bam"


##########################################################################################
##########################################################################################
## LIBRARY LEVEL

## get the bam file(s) to be merged
def get_bam_4_merge_bam_fastq2library(wc):
    return [
        get_final_bam_FASTQ(Wildcards(wc, {"id": id})) for id in SAMPLES[wc.sm][wc.lb]
    ]


def get_bam_4_merge_bam_low_qual_fastq2library(wc):
    return [
        get_final_bam_low_qual_FASTQ(Wildcards(wc, {"id": id}))
        for id in SAMPLES[wc.sm][wc.lb]
    ]


## get the (merged) bam file
def get_merged_bam_LB(wc):
    bam = get_bam_4_merge_bam_fastq2library(wc)
    if (
        len(bam) > 1
    ):  ## library consits of more than one fastq file: return 00_merged_fastq
        return f"{wc.folder}/02_library/00_merged_fastq/01_bam/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:  ## library consits of one fastq file: return return the location of the final library bam file
        return bam[0]


def get_merged_bam_low_qual_LB(wc):
    bam = get_bam_4_merge_bam_low_qual_fastq2library(wc)
    if (
        len(bam) > 1
    ):  ## library consits of more than one fastq file: return 00_merged_fastq
        return f"{wc.folder}/02_library/00_merged_fastq/01_bam_low_qual/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:  ## library consits of one fastq file: return return the location of the final library bam file
        return bam[0]


## get the bam file used to remove duplicates
def get_bam_4_markduplicates(wc):
    return get_merged_bam_LB(wc)


## get the bam file used to compte the damage
def get_bam_4_damage(wc):
    if remove_duplicates == "markduplicates":
        bam = f"{wc.folder}/02_library/01_duplicated/01_markduplicates/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    elif remove_duplicates == "dedup":
        bam = f"{wc.folder}/02_library/01_duplicated/01_dedup/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:
        bam = get_bam_4_markduplicates(wc)
    return bam


## get the final bam files at the library level
def get_final_bam_LB(wc):
    if run_damage_rescale:
        file = f"{wc.folder}/02_library/02_rescaled/01_mapDamage/{wc.sm}/{wc.lb}.{wc.genome}.bam"
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
    return [get_final_bam_LB(Wildcards(wc, {"lb": lb})) for lb in SAMPLES[wc.sm]]


def get_bam_4_merge_bam_low_qual_library2sample(wc):
    return [
        get_final_bam_low_qual_LB(Wildcards(wc, {"lb": lb})) for lb in SAMPLES[wc.sm]
    ]


## get the (merged) bam file
def get_merged_bam_SM(wc):
    bam = get_bam_4_merge_bam_library2sample(wc)
    if (
        len(bam) > 1
    ):  ## sample consits of more than one library return 00_merged_library
        return f"{wc.folder}/03_sample/00_merged_library/01_bam/{wc.sm}.{wc.genome}.bam"
    else:  ## library consits of one fastq file: return return the location of the final library bam file
        return bam[0]


def get_merged_bam_low_qual_SM(wc):
    # print(f"get_merged_bam_low_qual_SM: {wc}")
    bam = get_bam_4_merge_bam_low_qual_library2sample(wc)
    if (
        len(bam) > 1
    ):  ## sample consits of more than one library return 00_merged_library
        return f"{wc.folder}/03_sample/00_merged_library/01_bam_low_qual/{wc.sm}.{wc.genome}.bam"
    else:  ## library consits of one fastq file: return return the location of the final library bam file
        return bam[0]


## get the bam file used to realign indels
def get_bam_4_realign(wc):
    return get_merged_bam_SM(wc)


## get the bam file used to recompute the md flag
def get_bam_4_samtools_calmd(wc):
    if run_realign:
        return f"{wc.folder}/03_sample/01_realigned/01_realign/{wc.sm}.{wc.genome}.bam"
    else:
        return get_bam_4_realign(wc)


## get the final bam files at the sample level
def get_bam_4_final_bam(wc):
    if is_external_sample(wc.sm, wc.genome):
        return EXTERNAL_SAMPLES[wc.genome][wc.sm]
    if run_compute_md:
        return f"{wc.folder}/03_sample/02_md_flag/01_md_flag/{wc.sm}.{wc.genome}.bam"
    return get_bam_4_samtools_calmd(wc)


def get_bam_4_final_bam_low_qual(wc):
    return get_merged_bam_low_qual_SM(wc)


##########################################################################################
##########################################################################################

## get all files generated by the deamintion inference
def get_damage_output():
    if run_damage == "bamdamage":
        files = [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/02_library/04_bamdamage/{sm}/{lb}.{genome}.{type}.{ext}"
            for sm in SAMPLES
            for lb in SAMPLES[sm]
            for type in ["dam", "length"]
            for ext in ["pdf", "svg"]
            for genome in GENOMES
        ]
    elif run_damage == "mapDamage":
        files = [
            f"{{RESULT_DIR}}/04_stats/01_sparse_stats/02_library/04_mapDamage/{sm}/{lb}.{genome}_results_mapDamage/Fragmisincorporation_plot.pdf"
            for sm in SAMPLES
            for lb in SAMPLES[sm]
            for genome in GENOMES
        ]
    else:  ## if False
        files = []
    return files


##########################################################################################
## STATS
## get all individual stat table files to concatenate
def path_stats_by_level(wc):
    if wc.level == "FASTQ":
        paths = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.genome}/{sm}/{lb}/{id}/fastq_stats.csv"
            for sm in SAMPLES
            for lb in SAMPLES[sm]
            for id in SAMPLES[sm][lb]
        ]
    elif wc.level == "LB":
        paths = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.genome}/{sm}/{lb}/library_stats.csv"
            for sm in SAMPLES
            for lb in SAMPLES[sm]
        ]
    elif wc.level == "SM":
        paths = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.genome}/{sm}/sample_stats.csv"
            for sm in SAMPLES
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
                    "sm": paths[3],
                    "lb": paths[4],
                    "id": paths[5],
                    "genome": wc.genome,
                }
            )
        )
    elif paths[1] == "03_final_library":  ## library
        file = get_final_bam_LB(
            Wildcards(
                fromdict={
                    "folder": wc.folder,
                    "sm": paths[3],
                    "lb": paths[4],
                    "genome": wc.genome,
                }
            )
        )
    else:  ## the path is the right one (symlink)
        file = f"{wc.folder}/{wc.file}.{wc.genome}.bam"
    # print(f"get_bam_file: {file}")
    return file


## sex may be infered at teh sample or/and library level
def get_sex_file(wc):
    if str2bool(
        recursive_get_and_test(["sex_inference", wc.genome, "run"], ["False", "True"])
    ):
        return f"{wc.folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{wc.sm}.{wc.genome}_sex.txt"
    else:
        return f"{wc.folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{wc.sm}.{wc.genome}_nosex.txt"


##########################################################################################
## get all final files to run snakemake
def get_final_bam_files():
    final_bam = [
        f"{RESULT_DIR}/03_sample/03_final_sample/01_bam/{sm}.{genome}.bam"
        for sm in SAMPLES
        for genome in GENOMES
    ]
    return final_bam


def get_final_external_bam_files():
    final_bam = [
        f"{RESULT_DIR}/03_sample/03_final_sample/01_bam/{sm}.{genome}.bam"
        for genome in EXTERNAL_SAMPLES
        for sm in EXTERNAL_SAMPLES[genome]
    ]
    return final_bam


def get_final_bam_low_qual_files():
    final_bam_low_qual = [
        f"{RESULT_DIR}/03_sample/03_final_sample/01_bam_low_qual/{sm}.{genome}.bam"
        for sm in SAMPLES
        for genome in GENOMES
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
        f"{RESULT_DIR}/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/{sm}/{lb}/{id}_fastqc.zip"
        for sm in SAMPLES
        for lb in SAMPLES[sm]
        for id in SAMPLES[sm][lb]
    ] + [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_adapter_removal/{sm}/{lb}/{id}_fastqc.zip"
        for sm in SAMPLES
        for lb in SAMPLES[sm]
        for id in SAMPLES[sm][lb]
        if run_adapter_removal
    ]
    return fastqc


def get_samtools_stats_files():
    samtools_stats = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/{file}_stats.txt"
        for genome in GENOMES
        for sm in SAMPLES
        for lb in SAMPLES[sm]
        for id in SAMPLES[sm][lb]
        for file in [
            f"01_fastq/04_final_fastq/01_bam/{sm}/{lb}/{id}.{genome}",
            f"02_library/03_final_library/01_bam/{sm}/{lb}.{genome}",
            f"03_sample/03_final_sample/01_bam/{sm}.{genome}",
        ]
    ]
    return list(set(samtools_stats))  ## remove duplicates


def get_length_files():
    lengths = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/{file}_length.txt"
        for genome in GENOMES
        for sm in SAMPLES
        for lb in SAMPLES[sm]
        for id in SAMPLES[sm][lb]
        for file in [
            f"01_fastq/04_final_fastq/01_bam/{sm}/{lb}/{id}.{genome}",
            f"02_library/03_final_library/01_bam/{sm}/{lb}.{genome}",
            f"03_sample/03_final_sample/01_bam/{sm}.{genome}",
        ]
    ]
    return list(set(lengths))  ## remove duplicates


def get_length_files2():
    lengths = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_length.txt"
        for genome in GENOMES
        for sm in SAMPLES
    ]
    return list(set(lengths))  ## remove duplicates


def get_idxstats_files():
    idxstats = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/{file}_idxstats.txt"
        for genome in GENOMES
        for sm in SAMPLES
        for lb in SAMPLES[sm]
        for file in [
            f"02_library/03_final_library/01_bam/{sm}/{lb}.{genome}",
            f"03_sample/03_final_sample/01_bam/{sm}.{genome}",
        ]
    ]
    return list(set(idxstats))  ## remove duplicates


def get_qualimap_files():
    qualimap_files = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_qualimap"
        for genome in GENOMES
        for sm in SAMPLES
        if run_qualimap
    ]
    return qualimap_files


def get_sex_files():
    sex_files = []
    for genome in GENOMES:
        ext = (
            "sex"
            if str2bool(recursive_get(["sex_inference", genome, "run"], False))
            else "nosex"
        )

        sex_files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_{ext}.txt"
            for sm in SAMPLES
        ]
    return list(set(sex_files))  ## remove duplicates


def get_multiqc_files():
    multiqc_files = [
        f"{RESULT_DIR}/04_stats/02_separate_tables/{genome}/multiqc_fastqc.html"
        for genome in GENOMES
        if run_multiqc
    ]
    return multiqc_files


def get_imputation_files():
    files = [
        f"{RESULT_DIR}/03_sample/04_imputed/07_glimpse_sampled/{sm}_gp{GP}.{ext}"
        for sm in SAMPLES
        for GP in str2list(recursive_get(["imputation", "gp_filter"], "[0.8]"))
        for ext in ["bcf", "bcf.csi"]
        if run_imputation
    ] + [
        f"{RESULT_DIR}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp.txt"
        for sm in SAMPLES
        if run_imputation
    ]
    return files


def get_imputation_files_external():
    samples_ = [sm for gen in EXTERNAL_SAMPLES for sm in EXTERNAL_SAMPLES[gen]]
    files = [
        f"{RESULT_DIR}/03_sample/04_imputed/07_glimpse_sampled/{sm}_gp{GP}.{ext}"
        for sm in samples_
        for GP in str2list(recursive_get(["imputation", "gp_filter"], "[0.8]"))
        for ext in ["bcf", "bcf.csi"]
        if run_imputation
    ] + [
        f"{RESULT_DIR}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp.txt"
        for sm in samples_
        if run_imputation
    ]
    return files


def get_imputation_plots():
    files = [
        f"{RESULT_DIR}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp.svg"
        for sm in SAMPLES
        if run_imputation
    ]
    return files


def get_imputation_plots_external():
    samples_ = [sm for gen in EXTERNAL_SAMPLES for sm in EXTERNAL_SAMPLES[gen]]
    files = [
        f"{RESULT_DIR}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}_gp.svg"
        for sm in samples_
        if run_imputation
    ]
    return files
