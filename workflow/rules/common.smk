##########################################################################################
## FASTQ LEVEL

## get the path to the fastq file given sm, lb and id
def get_fastq_of_ID(wc):
    # print(f"get_fastq_of_ID: {wc}")
    if "_R1" == wc.idd[-3:]:
        filename = SAMPLES[wc.sm][wc.lb][wc.idd[:-3]]["Data1"]
    elif "_R2" == wc.idd[-3:]:
        filename = SAMPLES[wc.sm][wc.lb][wc.idd[:-3]]["Data2"]
    elif PAIRED_END:
        # elif PAIRED_END != 0:  ## SE library in a paired-end sample file
        filename = SAMPLES[wc.sm][wc.lb][wc.idd]["Data1"]
    else:
        filename = SAMPLES[wc.sm][wc.lb][wc.idd]["Data"]
    return filename


def get_fastq_4_cleaning(wc):
    folder = f"{wc.folder}/01_fastq/00_reads/01_files_orig/{wc.sm}/{wc.lb}"
    if is_paired_end(wc):
        filename = [
            f"{folder}/{wc.id}_R1.fastq.gz",
            f"{folder}/{wc.id}_R2.fastq.gz",
        ]
    else:
        filename = [f"{folder}/{wc.id}.fastq.gz"]
    return filename


def get_cleaning_folder_extension(wc):
    if is_collapse(wc):
        folder = "collapse"
    elif is_paired_end(wc):
        folder = "pe"
    else:
        folder = "se"
    return folder


## get the fastq file(s) used for mapping (output is a list)
def get_fastq_4_mapping(wc):
    cleaning = get_paramGrp(
        ["cleaning", "run"], ["adapterremoval", "fastp", "False"], wc
    )
    if cleaning == "adapterremoval":
        folder = f"{wc.folder}/01_fastq/01_trimmed/01_adapterremoval"
        if is_collapse(wc):
            filename = [f"{folder}_collapse/{wc.sm}/{wc.lb}/{wc.id}.fastq.gz"]
        elif is_paired_end(wc):
            filename = [
                f"{folder}_pe/{wc.sm}/{wc.lb}/{wc.id}_R1.fastq.gz",
                f"{folder}_pe/{wc.sm}/{wc.lb}/{wc.id}_R2.fastq.gz",
            ]
        else:
            filename = [f"{folder}_se/{wc.sm}/{wc.lb}/{wc.id}.fastq.gz"]

    elif cleaning == "fastp":
        folder = f"{wc.folder}/01_fastq/01_trimmed/02_fastp"
        if is_collapse(wc):
            filename = [f"{folder}_collapse/{wc.sm}/{wc.lb}/{wc.id}.fastq.gz"]
        elif is_paired_end(wc):
            filename = [
                f"{folder}_pe/{wc.sm}/{wc.lb}/{wc.id}_R1.fastq.gz",
                f"{folder}_pe/{wc.sm}/{wc.lb}/{wc.id}_R2.fastq.gz",
            ]
        else:
            filename = [f"{folder}_se/{wc.sm}/{wc.lb}/{wc.id}.fastq.gz"]

    else:  ## no cleaning
        filename = get_fastq_4_cleaning(wc)
    return filename


## fastqc may be run on the original or trimmed fastq files
def inputs_fastqc(wc):
    if "trim" in wc.type:
        return get_fastq_4_mapping(wc)
    else:
        return get_fastq_4_cleaning(wc)


## get the bam file used for sorting
def get_bam_4_sorting(wc):
    # print(f"get_bam_for_sorting: {wc}")
    if mapper == "bwa_aln":
        if is_paired_end(wc) and not is_collapse(wc):
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
    if str2bool(get_paramGrp(["filtering", "run"], ["True", "False"], wc)):
        file = f"{wc.folder}/01_fastq/03_filtered/01_bam_filter/{wc.sm}/{wc.lb}/{wc.id}.{wc.genome}.bam"
    else:
        file = f"{wc.folder}/01_fastq/02_mapped/03_bam_sort/{wc.sm}/{wc.lb}/{wc.id}.{wc.genome}.bam"
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
    ## library consits of more than one fastq file: return 00_merged_fastq
    if len(bam) > 1:
        return f"{wc.folder}/02_library/00_merged_fastq/01_bam/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:  ## library consits of one fastq file: return return the location of the final library bam file
        return bam[0]


def get_merged_bam_low_qual_LB(wc):
    bam = get_bam_4_merge_bam_low_qual_fastq2library(wc)
    ## library consits of more than one fastq file: return 00_merged_fastq
    if len(bam) > 1:
        return f"{wc.folder}/02_library/00_merged_fastq/01_bam_low_qual/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:  ## library consits of one fastq file: return return the location of the final library bam file
        return bam[0]


## get the bam file used to remove duplicates
def get_bam_4_markduplicates(wc):
    return get_merged_bam_LB(wc)


## get the bam file used to compte the damage
def get_bam_4_damage_rescale(wc):
    rm_duplicates = get_paramGrp(
        ["remove_duplicates", "run"], ["markduplicates", "dedup", "False"], wc
    )
    if rm_duplicates == "markduplicates":
        bam = f"{wc.folder}/02_library/01_duplicated/01_markduplicates/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    elif rm_duplicates == "dedup":
        ## use only on collapsed paired-end reads
        if not is_collapse(wc):
            LOGGER.warning(
                f"WARNING: DeDup should only be run on collapsed paired-end reads (SAMPLES[{wc.sm}][{wc.lb}])."
            )
        bam = f"{wc.folder}/02_library/01_duplicated/01_dedup/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:
        bam = get_bam_4_markduplicates(wc)
    return bam


## get the bam file used to trim the read ends with BamUtil
def get_bam_4_bamutil(wc):
    if str2bool(get_paramGrp(["damage_rescale", "run"], ["False", "True"], wc)):
        file = f"{wc.folder}/02_library/02_rescaled/01_mapDamage/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:
        file = get_bam_4_damage_rescale(wc)
    return file


## get the final bam files at the library level
def get_final_bam_LB(wc):
    if str2bool(get_paramGrp(["bamutil", "run"], ["False", "True"], wc)):
        file = (
            f"{wc.folder}/02_library/03_trim/01_bamutil/{wc.sm}/{wc.lb}.{wc.genome}.bam"
        )
    else:
        file = get_bam_4_bamutil(wc)
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
            for sm, smVals in SAMPLES.items()
            for lb in smVals
            for type in ["dam", "length"]
            for ext in ["pdf", "svg"]
            for genome in GENOMES
        ]
    elif run_damage == "mapDamage":
        files = [
            f"{{RESULT_DIR}}/04_stats/01_sparse_stats/02_library/04_mapDamage/{sm}/{lb}.{genome}_results_mapDamage/Fragmisincorporation_plot.pdf"
            for sm, smVals in SAMPLES.items()
            for lb in smVals
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
            for sm, smVals in SAMPLES.items()
            for lb, lbVals in smVals.items()
            for id in lbVals
        ]
    elif wc.level == "LB":
        paths = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.genome}/{sm}/{lb}/library_stats.csv"
            for sm, smVals in SAMPLES.items()
            for lb in smVals
        ]
    elif wc.level == "SM":
        paths = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.genome}/{sm}/sample_stats.csv"
            for sm in SAMPLES
        ] + [
            f"{wc.folder}/04_stats/02_separate_tables/{genome}/{sm}/sample_stats.csv"
            for g, gVal in EXTERNAL_SAMPLES.items()
            if g == wc.genome
            for sm in gVal
        ]
    else:
        LOGGER.error(
            f"ERROR: def path_stats_by_level({wc.level}): should never happen!"
        )
        os._exit(1)
    # print(wc)
    # print(paths)
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


## sex may be infered at the sample or/and library level
def get_sex_file(wc):
    if str2bool(get_param(["sex_inference", wc.genome, "run"], ["False", "True"])):
        return f"{wc.folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{wc.sm}.{wc.genome}_sex.txt"
    else:
        return f"{wc.folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{wc.sm}.{wc.genome}_nosex.txt"


def get_lb_stats(wc):
    if len(SAMPLES) and wc.sm in SAMPLES:
        return [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.genome}/{wc.sm}/{lb}/library_stats.csv"
            for lb in SAMPLES[wc.sm]
        ]
    return []


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
    if len(SAMPLES):
        ll = ["SM", "LB", "FASTQ"]
    else:
        ll = ["SM"]
    stat_csv = [f"{RESULT_DIR}/04_stats/03_summary/{level}_stats.csv" for level in ll]
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
        if len(SAMPLES)
    ]
    return plots


def get_fastqc_files():
    fastqc = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/{sm}/{lb}/{id}_fastqc.zip"
        for sm, smVals in SAMPLES.items()
        for lb, lbVals in smVals.items()
        for id in lbVals
    ] + [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_adapter_removal/{sm}/{lb}/{id}_fastqc.zip"
        for sm, smVals in SAMPLES.items()
        for lb, lbVals in smVals.items()
        for id in lbVals
        if get_paramGrp(
            ["adapterremoval", "run"],
            ["True", "False"],
            Wildcards(fromdict={"id": id, "lb": lb, "sm": sm}),
        )
    ]
    return fastqc


def get_samtools_stats_files():
    samtools_stats = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/{file}_stats.txt"
        for genome in GENOMES
        for sm, smVals in SAMPLES.items()
        for lb, lbVals in smVals.items()
        for id in lbVals
        for file in [
            f"01_fastq/04_final_fastq/01_bam/{sm}/{lb}/{id}.{genome}",
            f"02_library/03_final_library/01_bam/{sm}/{lb}.{genome}",
            f"03_sample/03_final_sample/01_bam/{sm}.{genome}"
        ]
    ] + [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/{file}_stats.txt"
            for genome, gVal in EXTERNAL_SAMPLES.items()
            for sm in gVal
            for file in [f"03_sample/03_final_sample/01_bam/{sm}.{genome}"]
    ]
    # print(samtools_stats)
    return list(set(samtools_stats))  ## remove duplicates


def get_length_files():
    lengths = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/{file}_length.txt"
        for genome in GENOMES
        for sm, smVals in SAMPLES.items()
        for lb, lbVals in smVals.items()
        for id in lbVals
        for file in [
            f"01_fastq/04_final_fastq/01_bam/{sm}/{lb}/{id}.{genome}",
            f"02_library/03_final_library/01_bam/{sm}/{lb}.{genome}",
            f"03_sample/03_final_sample/01_bam/{sm}.{genome}",
        ]
    ] + [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/{file}_length.txt"
            for genome, gVal in EXTERNAL_SAMPLES.items()
            for sm in gVal
            for file in [f"03_sample/03_final_sample/01_bam/{sm}.{genome}"]
    ]
    return list(set(lengths))  ## remove duplicates


def get_idxstats_files():
    idxstats = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/{file}_idxstats.txt"
        for genome in GENOMES
        for sm, smVals in SAMPLES.items()
        for lb in smVals
        for file in [
            f"02_library/03_final_library/01_bam/{sm}/{lb}.{genome}",
            f"03_sample/03_final_sample/01_bam/{sm}.{genome}",
        ]
    ] + [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/{file}_idxstats.txt"
        for genome, gVal in EXTERNAL_SAMPLES.items()
        if len(EXTERNAL_SAMPLES)
        for sm in gVal
        for file in [f"03_sample/03_final_sample/01_bam/{sm}.{genome}"]
    ]
    return list(set(idxstats))  ## remove duplicates


def get_qualimap_files():
    qualimap_files = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_qualimap"
        for genome in GENOMES
        for sm in SAMPLES
        if run_qualimap
    ] + [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_qualimap"
        for genome, gVal in EXTERNAL_SAMPLES.items()
        for sm in gVal
        if run_qualimap and len(EXTERNAL_SAMPLES)
    ]
    return qualimap_files


def get_multiqc_files():
    multiqc_files = [
        f"{RESULT_DIR}/04_stats/02_separate_tables/{genome}/multiqc_mapache.html"
        for genome in GENOMES
        if run_multiqc and len(SAMPLES)
    ] + [
        f"{RESULT_DIR}/04_stats/02_separate_tables/{genome}/multiqc_mapache.html"
        for genome, gVal in EXTERNAL_SAMPLES.items()
        for sm in gVal
        if run_multiqc and len(EXTERNAL_SAMPLES)
    ]
    return multiqc_files


##########################################################################################
## stats on final/external bam file
def get_sex_files():
    sex_files = []
    for genome in GENOMES:
        ext = (
            "sex"
            if str2bool(get_param(["sex_inference", genome, "run"], "False"))
            else "nosex"
        )

        ## SAMPLES
        sex_files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_{ext}.txt"
            for sm in SAMPLES
        ]

        ## EXTERNAL_SAMPLES
        sex_files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{genome}_{ext}.txt"
            for g, gVal in EXTERNAL_SAMPLES.items()
            if g == genome
            for sm in gVal
        ]
    return list(set(sex_files))  ## remove duplicates


def get_imputation_files():
    files = []
    for genome in GENOMES:
        if str2bool(get_param(["imputation", genome, "run"], ["False", "True"])):
            ## SAMPLES
            files += [
                f"{RESULT_DIR}/03_sample/04_imputed/07_glimpse_sampled/{sm}.{genome}_gp{GP}.{ext}"
                for sm in SAMPLES
                for GP in str2list(
                    get_param(["imputation", genome, "gp_filter"], "[0.8]")
                )
                for ext in ["bcf", "bcf.csi"]
            ]
            files += [
                f"{RESULT_DIR}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}.{genome}_gp.txt"
                for sm in SAMPLES
            ]

            ## EXTERNAL_SAMPLES
            files += [
                f"{RESULT_DIR}/03_sample/04_imputed/07_glimpse_sampled/{sm}.{genome}_gp{GP}.{ext}"
                for g, gVal in EXTERNAL_SAMPLES.items()
                if g == genome
                for sm in gVal
                for GP in str2list(
                    get_param(["imputation", genome, "gp_filter"], "[0.8]")
                )
                for ext in ["bcf", "bcf.csi"]
            ]
            files += [
                f"{RESULT_DIR}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}.{genome}_gp.txt"
                for g, gVal in EXTERNAL_SAMPLES.items()
                if g == genome
                for sm in gVal
            ]
    return files


def get_imputation_plots():
    files = []
    for genome in GENOMES:
        if str2bool(get_param(["imputation", genome, "run"], ["False", "True"])):
            ## SAMPLES
            files += [
                f"{RESULT_DIR}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}.{genome}_gp.svg"
                for sm in SAMPLES
            ]

            ## EXTERNAL_SAMPLES
            files += [
                f"{RESULT_DIR}/03_sample/04_imputed/07_glimpse_sampled/unphased/{sm}.{genome}_gp.svg"
                for g, gVal in EXTERNAL_SAMPLES.items()
                if g == genome
                for sm in gVal
            ]
    return files


#################################################################################################################
## multiqc input files
def get_files_4_multiqc(wc):
    files = []

    if len(SAMPLES):
        ## fastqc original
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/{sm}/{lb}/{id}_fastqc.zip"
            for sm, smVals in SAMPLES.items()
            for lb, lbVals in smVals.items()
            for id in lbVals
        ]

        ## adapterremoval
        files += [
            f"{RESULT_DIR}/01_fastq/01_trimmed/01_adapterremoval_{get_cleaning_folder_extension(Wildcards(fromdict={'id': id, 'lb': lb, 'sm': sm}))}/{sm}/{lb}/{id}.settings"
            for sm, smVals in SAMPLES.items()
            for lb, lbVals in smVals.items()
            for id in lbVals
            if get_paramGrp(
                ["cleaning", "run"],
                ["adapterremoval", "fastp", "False"],
                Wildcards(fromdict={"id": id, "lb": lb, "sm": sm}),
            )
            == "adapterremoval"
        ]

        ## fastp
        files += [
            f"{RESULT_DIR}/01_fastq/01_trimmed/01_fastp_{get_cleaning_folder_extension(Wildcards(fromdict={'id': id, 'lb': lb, 'sm': sm}))}/{sm}/{lb}/{id}.json"
            for sm, smVals in SAMPLES.items()
            for lb, lbVals in smVals.items()
            for id in lbVals
            if get_paramGrp(
                ["cleaning", "run"],
                ["adapterremoval", "fastp", "False"],
                Wildcards(fromdict={"id": id, "lb": lb, "sm": sm}),
            )
            == "fastp"
        ]

        ## fastqc trim
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_adapterremoval_{get_cleaning_folder_extension(Wildcards(fromdict={'id': id, 'lb': lb, 'sm': sm}))}/{sm}/{lb}/{id}_fastqc.zip"
            for sm, smVals in SAMPLES.items()
            for lb, lbVals in smVals.items()
            for id in lbVals
            if get_paramGrp(
                ["cleaning", "run"],
                ["adapterremoval", "fastp", "False"],
                Wildcards(fromdict={"id": id, "lb": lb, "sm": sm}),
            )
            == "adapterremoval"
        ] + [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_fastp_{get_cleaning_folder_extension(Wildcards(fromdict={'id': id, 'lb': lb, 'sm': sm}))}/{sm}/{lb}/{id}_fastqc.zip"
            for sm, smVals in SAMPLES.items()
            for lb, lbVals in smVals.items()
            for id in lbVals
            if get_paramGrp(
                ["cleaning", "run"],
                ["adapterremoval", "fastp", "False"],
                Wildcards(fromdict={"id": id, "lb": lb, "sm": sm}),
            )
            == "fastp"
        ]

        ## samtools_stats at final fastq bam
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/{sm}/{lb}/{id}.{wc.genome}_stats.txt"
            for sm, smVals in SAMPLES.items()
            for lb, lbVals in smVals.items()
            for id in lbVals
        ]

        ## picard markduplicates
        files += [
            f"{RESULT_DIR}/02_library/01_duplicated/01_markduplicates/{sm}/{lb}.{wc.genome}.stats"
            for sm, smVals in SAMPLES.items()
            for lb, lbVals in smVals.items()
            for id in lbVals
        ]

        ## samtools_stats at final library bam
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{sm}/{lb}.{wc.genome}_stats.txt"
            for sm, smVals in SAMPLES.items()
            for lb, lbVals in smVals.items()
            for id in lbVals
        ]

        ## samtools_stats at final bam file
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{wc.genome}_stats.txt"
            for sm, smVals in SAMPLES.items()
            for lb, lbVals in smVals.items()
            for id in lbVals
        ]

        ## qualimap at final bam file
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{wc.genome}_qualimap"
            for sm, smVals in SAMPLES.items()
            for lb, lbVals in smVals.items()
            for id in lbVals
            if str2bool(get_param(["stats", "qualimap"], False))
        ]

        files += [
            f"{RESULT_DIR}/04_stats/03_summary/SM_stats.{wc.genome}.csv",
            f"{RESULT_DIR}/04_stats/03_summary/LB_stats.{wc.genome}.csv",
            f"{RESULT_DIR}/04_stats/03_summary/FASTQ_stats.{wc.genome}.csv",
        ]

    if len(EXTERNAL_SAMPLES):
        ## samtools_stats at final bam file
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{wc.genome}_stats.txt"
            for genome, gVal in EXTERNAL_SAMPLES.items()
            for sm in gVal
        ]

        ## qualimap at final bam file
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{wc.genome}_qualimap"
            for genome, gVal in EXTERNAL_SAMPLES.items()
            for sm in gVal
            if str2bool(get_param(["stats", "qualimap"], False))
        ]

        files += [f"{RESULT_DIR}/04_stats/03_summary/SM_stats.{wc.genome}.csv"]

    return list(set(files))  ## remove duplicates


## return all files containing version information
def get_version_file_of_tools():
    ## extract tools form the conda environment file
    filename = "config/mapache-env.yaml"
    tools = []
    with open(filename) as file:
        for line in file:
            if "dependencies" in line:
                for line in file:
                    if (
                        "r-" not in line
                        and "charset-normalizer" not in line
                        and "mamba" not in line
                    ):
                        tools.append(
                            line.strip()
                            .replace("- ", "")
                            .replace("-bio", "")
                            .split("=")[0]
                        )

    files = [
        f"{RESULT_DIR}/04_stats/02_separate_tables/software/{tool}.txt"
        for tool in tools
    ]
    return files
