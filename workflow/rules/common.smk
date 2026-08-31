##########################################################################################
## FASTQ LEVEL

def get_fastq_name_original(wc):
    ## single or paired?
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

    
def idd_is_remote(sm, lb, idd):
    filename = get_fastq_name_original(Wildcards(fromdict={"sm": sm, "lb": lb, "idd": idd}))
    return not os.path.exists(filename) and filename.startswith(("ftp", "http://", "https://")) 



## check if the string as a valuable md5string 
def test_md5(md5, wc):
    MD5_RE = re.compile(r'^[a-fA-F0-9]{32}$')
    if bool(MD5_RE.match(md5)):
        return md5
    
    if  md5 != "" and md5 != "NA" and md5 != "NULL":
        LOGGER.warning(
            f"ERROR: The md5sum of {wc.sm}/{wc.lb}/{wc.idd} is not valid!"
        )
        os._exit(1)

    return ""


## get the md5 of the given ID (if available, otherwise return '')
def get_md5_of_ID(wc):
    if "_R1" == wc.idd[-3:]:
        if "MD5_1" in SAMPLES[wc.sm][wc.lb][wc.idd[:-3]]:
            md5 =   (SAMPLES[wc.sm][wc.lb][wc.idd[:-3]]["MD5_1"], wc)
        else:
            md5 = ""
    elif "_R2" == wc.idd[-3:]:
        if "MD5_2" in SAMPLES[wc.sm][wc.lb][wc.idd[:-3]]:
            md5 = test_md5(SAMPLES[wc.sm][wc.lb][wc.idd[:-3]]["MD5_2"], wc)
        else:
            md5 = ""
    elif PAIRED_END:
        # elif PAIRED_END != 0:  ## SE library in a paired-end sample file
        if "MD5_1" in SAMPLES[wc.sm][wc.lb][wc.idd]:
            md5 = test_md5(SAMPLES[wc.sm][wc.lb][wc.idd]["MD5_1"], wc)
        else:
            md5 = ""    
    else:
        if "MD5" in SAMPLES[wc.sm][wc.lb][wc.idd]:
            md5 = test_md5(SAMPLES[wc.sm][wc.lb][wc.idd]["MD5"], wc)
        else:
            md5 = ""
    
    #print(md5)
    return md5


def get_fastq_4_cleaning(wc):
    if run_subsampling:
        folder = f"{wc.folder}/01_fastq/00_reads/00_files_subsample/{wc.sm}/{wc.lb}"
    else:
        folder = f"{wc.folder}/01_fastq/00_reads/00_files_orig/{wc.sm}/{wc.lb}"

    if is_paired_end(wc):
        filename = [
            f"{folder}/{wc.id}_R1.fastq.gz",
            f"{folder}/{wc.id}_R2.fastq.gz"
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
    # print(wc)
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
    # print(filename)
    return filename


## fastqc may be run on the original or trimmed fastq files
def inputs_fastqc(wc):
    if "trim" in wc.type:
        return get_fastq_4_mapping(wc)
    else:
        return get_fastq_4_cleaning(wc)


## get the needed index files
def get_fasta_index(wc):
    ref = f"{wc.folder}/00_reference/{wc.genome}/{wc.genome}.fasta"
    if mapper == "bwa_aln" or mapper == "bwa_mem" or mapper == "bwa_mem":
        ext = ["", ".sa", ".amb", ".ann", ".bwt", ".pac"]
    elif mapper == "bowtie2":
        ext = ["", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
    elif mapper == "bowtie2l":
        ext = [
            "",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ]
    elif mapper == "giraffe":
        ext = [
            #".min",
            #".dist",
            ".gbz",
            #".min.gbz",
            #".dist.gbz",
        ]
    else:
        LOGGER.error(
            f"ERROR: The parameter config[mapping][mapper] is not correctly specified: {mapper} is unknown!"
        )
        os._exit(1)
    return [ref + x for x in ext]


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
    elif mapper == "bowtie2" or mapper == "bowtie2l":
        folder = "02_bowtie2"
    elif mapper == "vg_giraffe" or mapper == "vg_giraffe_gam":
        folder = "02_vg_giraffe"
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


def get_final_gam_FASTQ(wc):
    if str2bool(get_paramGrp(["filtering", "run"], ["True", "False"], wc)):
        file = f"{wc.folder}/01_fastq/03_filtered/01_gam_filter/{wc.sm}/{wc.lb}/{wc.id}.{wc.genome}.gam"
    else:
        file = f"{wc.folder}/01_fastq/02_mapped/02_vg_giraffe/{wc.sm}/{wc.lb}/{wc.id}.{wc.genome}.gam"
    # print(f"get_gam_4_final_fastq: {file}")
    return file


def get_final_gam_low_qual_FASTQ(wc):
    return f"{wc.folder}/01_fastq/03_filtered/01_gam_filter_low_qual/{wc.sm}/{wc.lb}/{wc.id}.{wc.genome}.gam"


##########################################################################################
##########################################################################################
## LIBRARY LEVEL


## get the bam file(s) to be merged
def get_bam_4_merge_bam_fastq2library(wc):
    return [
        get_final_bam_FASTQ(Wildcards(wc, {"id": id})) for id in SAMPLES[wc.sm][wc.lb]
    ]

def get_gam_4_merge_gam_fastq2library(wc):
    return [
        get_final_gam_FASTQ(Wildcards(wc, {"id": id})) for id in SAMPLES[wc.sm][wc.lb]
    ]


def get_bam_4_merge_bam_low_qual_fastq2library(wc):
    return [
        get_final_bam_low_qual_FASTQ(Wildcards(wc, {"id": id}))
        for id in SAMPLES[wc.sm][wc.lb]
    ]


def get_gam_4_merge_gam_low_qual_fastq2library(wc):
    return [
        get_final_gam_low_qual_FASTQ(Wildcards(wc, {"id": id}))
        for id in SAMPLES[wc.sm][wc.lb]
    ]


## get the (merged) bam file
def get_merged_bam_LB(wc):
    bam = get_bam_4_merge_bam_fastq2library(wc)
    ## library consists of more than one fastq file: return 00_merged_fastq
    if len(bam) > 1:
        return f"{wc.folder}/02_library/00_merged_fastq/01_bam/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:  ## library consists of one fastq file: return return the location of the final library bam file
        return bam[0]


def get_merged_bam_low_qual_LB(wc):
    bam = get_bam_4_merge_bam_low_qual_fastq2library(wc)
    ## library consists of more than one fastq file: return 00_merged_fastq
    if len(bam) > 1:
        return f"{wc.folder}/02_library/00_merged_fastq/01_bam_low_qual/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:  ## library consists of one fastq file: return return the location of the final library bam file
        return bam[0]


## get the (merged) gam file
def get_merged_gam_LB(wc):
    gam = get_gam_4_merge_gam_fastq2library(wc)
    ## library consists of more than one fastq file: return 00_merged_fastq
    if len(gam) > 1:
        return f"{wc.folder}/02_library/00_merged_fastq/01_gam/{wc.sm}/{wc.lb}.{wc.genome}.gam"
    else:  ## library consists of one fastq file: return return the location of the final library gam file
        return gam[0]
        

def get_merged_gam_low_qual_LB(wc):
    gam = get_gam_4_merge_gam_low_qual_fastq2library(wc)
    ## library consists of more than one fastq file: return 00_merged_fastq
    if len(gam) > 1:
        return f"{wc.folder}/02_library/00_merged_fastq/01_gam_low_qual/{wc.sm}/{wc.lb}.{wc.genome}.gam"
    else:  ## library consists of one fastq file: return the location of the final library gam file
        return gam[0]


## get the bam file used to remove duplicates
def get_bam_4_markduplicates(wc):
    return get_merged_bam_LB(wc)


##-------------------------------------------------------------------------------------------------------------------------------
## get the bam file after duplicate removal
def get_bam_4_after_rmDup(wc):
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
    # print(bam)
    return bam


## get all bam file(s) after the duplicates are removed
def get_all_bam_after_rmDup():
    files = [
        f"{RESULT_DIR}/02_library/02_bam_after_rmDup/01_bam/{sm}/{lb}.{genome}.bam"
        for sm in SAMPLES
        for lb in SAMPLES[sm]
        for genome in GENOMES
    ]
    # print(files)
    return files


def get_all_gam_after_rmDup():
    files = [
        f"{RESULT_DIR}/02_library/02_gam_after_rmDup/01_gam/{sm}/{lb}.{genome}.gam"
        for sm in SAMPLES
        for lb in SAMPLES[sm]
        for genome in GENOMES
    ]
    # print(files)
    return files
##-------------------------------------------------------------------------------------------------------------------------------


def get_bam_after_rmDup(wc):
    return get_bam_4_damage_rescale(wc)


## get the bam file for mapDamage2
def get_bam_4_damage_rescale(wc):
    if str2bool(get_paramGrp(["stats_after_rmDup", "run"], ["False", "True"], wc)):
        file = f"{wc.folder}/02_library/02_bam_after_rmDup/01_bam/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:
        file = get_bam_4_after_rmDup(wc)
    # print(file)
    return file


## get the bam file used to trim the read ends with BamUtil
def get_bam_4_bamutil(wc):
    if str2bool(get_paramGrp(["damage_rescale", "run"], ["False", "True"], wc)):
        file = f"{wc.folder}/02_library/02_rescaled/01_mapDamage/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:
        file = get_bam_4_damage_rescale(wc)
    # print(file)
    return file


## get the bam file used to mask the read ends with BamRefine
def get_bam_4_bamrefine(wc):
    if str2bool(get_paramGrp(["bamutil", "run"], ["False", "True"], wc)):
        file = (
            f"{wc.folder}/02_library/03_trim/01_bamutil/{wc.sm}/{wc.lb}.{wc.genome}.bam"
        )
    else:
        file = get_bam_4_bamutil(wc)
    # print(file)
    return file


## get the final bam files at the library level
def get_final_bam_LB(wc):
    if str2bool(get_paramGrp(["bamrefine", "run"], ["False", "True"], wc)):
        file = f"{wc.folder}/02_library/03_trim/02_bamrefine/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:
        file = get_bam_4_bamrefine(wc)
    return file


def get_final_bam_low_qual_LB(wc):
    if str2bool(get_paramGrp(["stats_after_rmDup", "run"], ["False", "True"], wc)):
        return f"{wc.folder}/02_library/02_bam_after_rmDup/01_bam_low_qual/{wc.sm}/{wc.lb}.{wc.genome}.bam"
    else:
        return get_merged_bam_low_qual_LB(wc)


def get_final_gam_LB(wc):
    return get_merged_gam_LB(wc)


def get_final_gam_low_qual_LB(wc):
    return get_merged_gam_low_qual_LB(wc)


##########################################################################################
##########################################################################################
##########################################################################################
## SAMPLE LEVEL


## get the corresponding 'old' sample name
def sm_final_2_sm(sm_final, lb):
    sm = SAMPLES_TABLE.loc[
        SAMPLES_TABLE[
            (SAMPLES_TABLE["SM_FINAL"] == sm_final) & (SAMPLES_TABLE["LB"] == lb)
        ].first_valid_index(),
        "SM",
    ]
    # print(f"{sm_final}/{lb} => {sm}/{lb}")
    return sm


## return all lines of the given sample as 'SM-LB' table
def sm_final_2_sm_table(sm_final):
    sm = (
        SAMPLES_TABLE[SAMPLES_TABLE["SM_FINAL"] == sm_final][["SM", "LB"]]
        .drop_duplicates()
        .reset_index()
    )
    return sm


## return true if the sample name has changed (and thus RG has to be adapted)
def sm_changed(sm_final):
    sm = SAMPLES_TABLE[SAMPLES_TABLE["SM_FINAL"] == sm_final]["SM"].tolist()
    return not all(element == sm_final for element in sm)


## get the bam file(s) to be merged
def get_bam_4_merge_bam_library2sample(wc):
    # print([Wildcards(wc, {"sm": sm_final_2_sm(wc.sm, lb), "lb": lb}) for lb in SAMPLES_FINAL[wc.sm]])
    return [
        get_final_bam_LB(Wildcards(wc, {"sm": sm_final_2_sm(wc.sm, lb), "lb": lb}))
        for lb in SAMPLES_FINAL[wc.sm]
    ]


def get_bam_4_merge_bam_low_qual_library2sample(wc):
    df = sm_final_2_sm_table(wc.sm)
    return [
        get_final_bam_low_qual_LB(Wildcards(wc, {"sm": sm, "lb": lb}))
        for sm, lb in zip(df["SM"], df["LB"])
    ]


def get_gam_4_merge_gam_library2sample(wc):
    # print([Wildcards(wc, {"sm": sm_final_2_sm(wc.sm, lb), "lb": lb}) for lb in SAMPLES_FINAL[wc.sm]])
    return [
        get_final_gam_LB(Wildcards(wc, {"sm": sm_final_2_sm(wc.sm, lb), "lb": lb}))
        for lb in SAMPLES_FINAL[wc.sm]
    ]


def get_gam_4_merge_gam_low_qual_library2sample(wc):
    df = sm_final_2_sm_table(wc.sm)
    return [
        get_final_gam_low_qual_LB(Wildcards(wc, {"sm": sm, "lb": lb}))
        for sm, lb in zip(df["SM"], df["LB"])
    ]


## get the (merged) bam file
def get_merged_bam_SM(wc):
    bam = get_bam_4_merge_bam_library2sample(wc)
    # print(bam)
    if len(bam) > 1 or sm_changed(
        wc.sm
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


def get_merged_gam_SM(wc):
    gam = get_gam_4_merge_gam_library2sample(wc)
    # print(gam)
    if len(gam) > 1 or sm_changed(
        wc.sm
    ):  ## sample consits of more than one library return 00_merged_library
        return f"{wc.folder}/03_sample/00_merged_library/01_gam/{wc.sm}.{wc.genome}.gam"
    else:  ## library consits of one fastq file: return return the location of the final library gam file
        return gam[0]


def get_merged_gam_low_qual_SM(wc):
    # print(f"get_merged_gam_low_qual_SM: {wc}")
    gam = get_gam_4_merge_gam_low_qual_library2sample(wc)
    if (
        len(gam) > 1
    ):  ## sample consits of more than one library return 00_merged_library
        return f"{wc.folder}/03_sample/00_merged_library/01_gam_low_qual/{wc.sm}.{wc.genome}.gam"
    else:  ## library consits of one fastq file: return return the location of the final library gam file
        return gam[0]


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


def get_gam_4_final_gam(wc):
    return get_merged_gam_SM(wc)


def get_gam_4_final_gam_low_qual(wc):
    return get_merged_gam_low_qual_SM(wc)


##########################################################################################
##########################################################################################


##########################################################################################
## STATS
## get all individual stat table files to concatenate
def path_stats_by_level(wc):
    # print(wc)
    if wc.level == "FASTQ":
        paths = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.genome}/{sm}/{lb}/{id}/FASTQ_stats.csv"
            for sm, smVals in SAMPLES.items()
            for lb, lbVals in smVals.items()
            for id in lbVals
        ]
    elif wc.level == "LB":
        paths = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.genome}/{sm}/{lb}/LB_stats.csv"
            for sm, smVals in SAMPLES.items()
            for lb in smVals
        ]
    elif wc.level == "LB_rmDup":
        paths = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.genome}/{sm}/{lb}/LB_rmDup_stats.csv"
            for sm, smVals in SAMPLES.items()
            for lb in smVals
        ]
    elif wc.level == "SM":
        paths = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.genome}/{sm}/SM_stats.csv"
            for sm in SAMPLES_FINAL
        ] + [
            f"{wc.folder}/04_stats/02_separate_tables/{genome}/{sm}/SM_stats.csv"
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


## sex may be inferred at the sample or/and library level
def get_sex_file_sample(wc):
    if str2bool(get_param(["sex_inference", wc.genome, "run"], ["False", "True"])):
        file = f"{wc.folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{wc.sm}.{wc.genome}_sex.txt"
    else:
        file = f"{wc.folder}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{wc.sm}.{wc.genome}_nosex.txt"
    return file


def get_sex_file_library_rmDup(wc):
    if str2bool(get_param(["sex_inference", wc.genome, "run"], ["False", "True"])):
        file = f"{wc.folder}/04_stats/01_sparse_stats/02_library/02_bam_after_rmDup/01_bam/{wc.sm}/{wc.lb}.{wc.genome}_sex.txt"
    else:
        file = f"{wc.folder}/04_stats/01_sparse_stats/02_library/02_bam_after_rmDup/01_bam/{wc.sm}/{wc.lb}.{wc.genome}_nosex.txt"
    # print(file)
    return file


def get_sex_file_library(wc):
    if str2bool(get_param(["sex_inference", wc.genome, "run"], ["False", "True"])):
        file = f"{wc.folder}/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{wc.sm}/{wc.lb}.{wc.genome}_sex.txt"
    else:
        file = f"{wc.folder}/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{wc.sm}/{wc.lb}.{wc.genome}_nosex.txt"
    # print(file)
    return file


def get_lb_stats(wc):
    df = sm_final_2_sm_table(wc.sm)
    if len(SAMPLES):
        files = [
            f"{wc.folder}/04_stats/02_separate_tables/{wc.genome}/{sm}/{lb}/LB_stats.csv"
            for sm, lb in zip(df["SM"], df["LB"])
        ]
    else:
        files = []
    # prinf(files)
    return files


##########################################################################################
## get all final files to run snakemake
def get_final_bam_files():
    final_bam = [
        f"{RESULT_DIR}/03_sample/03_final_sample/01_bam/{sm}.{genome}.bam"
        for sm in SAMPLES_FINAL
        for genome in GENOMES
    ]
    return final_bam


def get_final_gam_files():
    final_gam = [
        f"{RESULT_DIR}/03_sample/03_final_sample/01_gam/{sm}.{genome}.gam"
        for sm in SAMPLES_FINAL
        for genome in GENOMES
    ]
    return final_gam


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
        for sm in SAMPLES_FINAL
        for genome in GENOMES
        if save_low_qual
    ]
    return final_bam_low_qual


def get_final_gam_low_qual_files():
    final_gam_low_qual = [
        f"{RESULT_DIR}/03_sample/03_final_sample/01_gam_low_qual/{sm}.{genome}.gam"
        for sm in SAMPLES_FINAL
        for genome in GENOMES
        if save_low_qual
    ]
    return final_gam_low_qual
    

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


def get_fastqc_files2():
    ## orig
    fastqc = [
        f"{RESULT_DIR}/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/{sm}/{lb}/{id}_fastqc.zip"
        for sm, smVals in SAMPLES.items()
        for lb, lbVals in smVals.items()
        for id in lbVals
    ]

    ## adapterremoval
    prefix = (
        f"{RESULT_DIR}/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_adapterremoval"
    )
    for sm, smVals in SAMPLES.items():
        for lb, lbVals in smVals.items():
            for id in lbVals:
                wc = Wildcards(fromdict={"id": id, "lb": lb, "sm": sm})
                if get_paramGrp(["adapterremoval", "run"], ["True", "False"], wc):
                    if is_collapse(wc):
                        fastqc = fastqc + [
                            f"{prefix}_collapse/{sm}/{lb}/{id}_fastqc.zip"
                        ]
                    elif is_paired_end(wc):
                        fastqc = fastqc + [
                            f"{prefix}_pe/{wc.sm}/{wc.lb}/{wc.id}_R1.fastq.gz",
                            f"{prefix}_pe/{wc.sm}/{wc.lb}/{wc.id}_R2.fastq.gz",
                        ]
                    else:
                        fastqc = fastqc + [
                            f"{prefix}_se/{wc.sm}/{wc.lb}/{wc.id}.fastq.gz"
                        ]
    return fastqc


def get_multiqc_files(level="SM"):
    if not str2bool(get_param(["stats", "multiqc", "run"], ["False", "True"])):
        return []

    if level == "SM":
        files = [
            f"{RESULT_DIR}/04_stats/02_separate_tables/{genome}/multiqc_mapache.html"
            for genome in GENOMES + list(EXTERNAL_SAMPLES.keys())
        ]
    elif level == "LB_rmDup":
        files = [
            f"{RESULT_DIR}/04_stats/02_separate_tables/{genome}/multiqc_mapache_library_rmDup.html"
            for genome in GENOMES
        ]
    else:
        LOGGER.error(f"ERROR: def get_multiqc_files({level}): should never happen!")
        os._exit(1)

    # print(files)
    return files


##########################################################################################
## stats on final/external bam file
def get_imputation_files():
    files = []
    for genome in GENOMES:
        run_imputation = get_param(
            ["imputation", genome, "run"], ["False", "glimpse1", "glimpse2"]
        )
        if str(run_imputation) != "False":
            if run_imputation == "glimpse1":
                folder = "08_glimpse_sampled"
            else:
                folder = "07_gp_filtered"

            ## SAMPLES
            files += [
                f"{RESULT_DIR}/03_sample/04_imputed/{folder}/{sm}.{genome}_gp{GP}.{ext}"
                for sm in SAMPLES_FINAL
                for GP in str2list(
                    get_param(["imputation", genome, "gp_filter"], "[0.8]")
                )
                for ext in ["bcf", "bcf.csi"]
            ] + [
                f"{RESULT_DIR}/03_sample/04_imputed/07_gp_filtered/{sm}.{genome}_gp.txt"
                for sm in SAMPLES_FINAL
            ]

            files += [
                f"{RESULT_DIR}/03_sample/04_imputed/{folder}/{sm}.{genome}_gp{GP}.{ext}"
                for g, gVal in EXTERNAL_SAMPLES.items()
                if g == genome
                for sm in gVal
                for GP in str2list(
                    get_param(["imputation", genome, "gp_filter"], "[0.8]")
                )
                for ext in ["bcf", "bcf.csi"]
            ] + [
                f"{RESULT_DIR}/03_sample/04_imputed/07_gp_filtered/{sm}.{genome}_gp.txt"
                for g, gVal in EXTERNAL_SAMPLES.items()
                if g == genome
                for sm in gVal
            ]
    return files


def get_imputation_plots():
    files = []
    for genome in GENOMES:
        run_imputation = get_param(
            ["imputation", genome, "run"], ["False", "glimpse1", "glimpse2"]
        )
        if str(run_imputation) != "False":
            ## SAMPLES
            files += [
                f"{RESULT_DIR}/03_sample/04_imputed/07_gp_filtered/{sm}.{genome}_gp.svg"
                for sm in SAMPLES_FINAL
            ]

            ## EXTERNAL_SAMPLES
            files += [
                f"{RESULT_DIR}/03_sample/04_imputed/07_gp_filtered/{sm}.{genome}_gp.svg"
                for g, gVal in EXTERNAL_SAMPLES.items()
                if g == genome
                for sm in gVal
            ]
    return files


#################################################################################################################
## multiqc input files
def get_files_4_multiqc_library_rmDup(wc):
    files = []

    if len(SAMPLES):
        ## FASTQ ################################################
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

        ## LIBRARY ################################################
        ## preseq
        files += get_preseq_files_genome("tsv", wc.genome)

        ## picard markduplicates
        files += [
            f"{RESULT_DIR}/02_library/01_duplicated/01_markduplicates/{sm}/{lb}.{wc.genome}.stats"
            for sm, smVals in SAMPLES.items()
            for lb in smVals
            if get_paramGrp(
                ["remove_duplicates", "run"],
                ["False", "markduplicates", "dedup"],
                Wildcards(fromdict={"lb": lb, "sm": sm}),
            )
            == "markduplicates"
        ]

        ## samtools_stats after remDup library bam
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/02_library/02_bam_after_rmDup/01_bam/{sm}/{lb}.{wc.genome}_stats.txt"
            for sm, smVals in SAMPLES.items()
            for lb in smVals
            if str2bool(
                get_paramGrp(
                    ["stats_after_rmDup", "run"],
                    ["False", "True"],
                    Wildcards(fromdict={"lb": lb, "sm": sm}),
                )
            )
        ]

        ## qualimap at bam file (LB_rmDup)
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/02_library/02_bam_after_rmDup/01_bam/{sm}/{lb}.{wc.genome}_qualimap"
            for sm, smVals in SAMPLES.items()
            for lb in smVals
            if compute_qualimap("LB_rmDup", Wildcards(fromdict={"lb": lb, "sm": sm}))
        ]

        files += [
            f"{RESULT_DIR}/04_stats/02_separate_tables/{wc.genome}/LB_rmDup_stats.csv",
            f"{RESULT_DIR}/04_stats/02_separate_tables/{wc.genome}/FASTQ_stats.csv",
        ]

    # print(list(set(files)))
    return list(set(files))  ## remove duplicates


def get_files_4_multiqc(wc):
    files = get_files_4_multiqc_library_rmDup(wc)

    if len(SAMPLES):
        ## samtools_stats at final library bam
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{sm}/{lb}.{wc.genome}_stats.txt"
            for sm, smVals in SAMPLES.items()
            for lb in smVals
        ]

        ## qualimap at final bam file (LB)
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/{sm}/{lb}.{wc.genome}_qualimap"
            for sm, smVals in SAMPLES.items()
            for lb in smVals
            if compute_qualimap("LB", Wildcards(fromdict={"lb": lb, "sm": sm}))
        ]

        ## SAMPLE ################################################
        ## samtools_stats at final bam file
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{wc.genome}_stats.txt"
            for sm in list(SAMPLES_FINAL.keys())
            + get_external_samples_of_genome(wc.genome)
        ]

        ## qualimap at final bam file
        files += [
            f"{RESULT_DIR}/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/{sm}.{wc.genome}_qualimap"
            for sm in list(SAMPLES_FINAL.keys())
            + get_external_samples_of_genome(wc.genome)
            if compute_qualimap("SM", Wildcards(fromdict={"sm": sm}))
        ]

        files += [
            f"{RESULT_DIR}/04_stats/02_separate_tables/{wc.genome}/SM_stats.csv",
            f"{RESULT_DIR}/04_stats/02_separate_tables/{wc.genome}/LB_stats.csv",
        ]

    # print(list(set(files)))
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
    return tools
