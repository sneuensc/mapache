import pandas as pd
import numpy as np
import itertools
import pathlib


##########################################################################################
##########################################################################################
## all functions for main snakemake file

## functions to read the config file
def get_param1(key, default):
    return config.get(key, default)

def get_param2(key1, key2, default):
    return config.get(key1, {}).get(key2, default)

def get_param3(key1, key2, key3, default):
    return config.get(key1, {}).get(key2, {}).get(key3, default)    

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
    mem_incre = int(config.get(module, {}).get("mem_increment", memory_increment_ratio*mem_start))
    return int(1024 * ((attempt-1) * mem_incre + mem_start))
    
    
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
    time_incre = int(config.get(module, {}).get("time_increment", runtime_increment_ratio*time_start))
    return int(60 * ((attempt-1) * time_incre + time_start))
    #return convert_time(60*60 * ((attempt-1) * time_incre + time_start))
    
    
## get the number of threads of the given parameter
def get_threads(module, default=1):
    return int(config.get(module, {}).get("threads", default))


## define how to quantify the deamination pattern
# def get_damage(run_damage):
#     files=[]
#     if run_damage == 'bamdamage':
#         for genome in GENOME:
#             files+=[("results/03_sample/{SM}/{LB}/library_bamdamage/{LB}.{genome}.dam.svg").
#                 format(SM=row['SM'], LB=row['LB'], genome=genome) for index, row in all_libraries.iterrows()]
#     elif run_damage == 'mapDamage':
#         for genome in GENOME:
#             files+=[("{SM}/{LB}/library_mapDamage/{LB}.{genome}_results_mapDamage/Fragmisincorporation_plot.pdf").
#                 format(SM=row['SM'], LB=row['LB'], genome=genome) for index, row in all_libraries.iterrows()]
# 
#     return (files)
    

##########################################################################################
##########################################################################################
## all functions for fastq files

def get_from_sample_file(col, SM, LB, ID):
    filename = db[(db["ID"]==ID) & (db["LB"]==LB) & (db["SM"]==SM)][col].values
    return(filename)


def get_fastq_of_ID(wildcards):
	if "_R1" == wildcards.id_fastq[-3:]:
		filename = get_from_sample_file("Data1", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq[:-3])
	elif "_R2" == wildcards.id_fastq[-3:]:
		filename = get_from_sample_file("Data2", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq[:-3])
	elif paired_end != 0:   ## SE library in a paired-end sample file
		filename = get_from_sample_file("Data1", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)
	else:
		filename = get_from_sample_file("Data", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)
	return(filename)
	
	
def get_fastq_for_mapping(wildcards):
	if run_adapter_removal:
		if paired_end == 1 and str(get_from_sample_file("Data2", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0])!="nan":
			folder = f"results/01_fastq/01_trimmed/01_files_trim_collapsed/{wildcards.id_sample}/{wildcards.id_library}"
		else:
			folder = f"results/01_fastq/01_trimmed/01_files_trim/{wildcards.id_sample}/{wildcards.id_library}"
	else:
		folder = f"results/01_fastq/00_reads/01_files_orig/{wildcards.id_sample}/{wildcards.id_library}"

	if paired_end == 2:
		filename=[f"{folder}/{wildcards.id_fastq}_R1.fastq.gz", f"{folder}/{wildcards.id_fastq}_R2.fastq.gz"]
	else:
		filename=[f"{folder}/{wildcards.id_fastq}.fastq.gz"]
	return (filename)

def get_fastq_for_mapping_pe(wildcards):
	if run_adapter_removal:
		folder = f"results/01_fastq/01_trimmed/01_files_trim"
	else:
		folder = f"results/01_fastq/00_reads/01_files_orig"

	return (f"{folder}/{wildcards.id_sample}/{wildcards.id_library}/{wildcards.id_fastq}_R{wildcards.id_read}.fastq.gz")


def get_bam_for_sorting(wildcards):
	if mapper == "bwa_aln":
		if paired_end == 2 and str(get_from_sample_file("Data2", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0]) != "nan":
			folder = "02_bwa_sampe"
		else:
			folder = "02_bwa_samse"
	elif mapper == "bwa_mem":
		folder = "02_bwa_mem"
	elif mapper == "bowtie2":
		folder = "02_bowtie2"
	else:
		print(f"ERROR: The parameter mapper is not correctly specified: {mapper} is unknown!")
		os._exit(0)
	return (f"results/01_fastq/02_mapped/{folder}/{wildcards.id_sample}/{wildcards.id_library}/{wildcards.id_fastq}.{wildcards.id_genome}.bam")


	
##########################################################################################
##########################################################################################
## all functions for libraries

def get_final_bam_library(wildcards):
    if run_mapDamage_rescale:
        bam = f"results/02_library/02_rescaled/01_mapDamage/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}.bam"
    elif extract_duplicates:
        bam = f"results/02_library/01_duplicated/01_rmdup/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}_mapped.bam"
    elif run_remove_duplicates:
        bam = f"results/02_library/01_duplicated/01_rmdup/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}.bam"
    else:
        bam = f"results/02_library/00_merged_fastq/01_bam/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}.bam"
    return (bam)


def get_mapDamage_bam(wildcards):
	if extract_duplicates:
		bam = f"results/02_library/01_duplicated/01_rmdup/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}_mapped.bam"
	elif run_remove_duplicates:
		bam = f"results/02_library/01_duplicated/01_rmdup/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}.bam"
	else: 
		bam = f"results/02_library/00_merged_fastq/01_bam/{wildcards.id_sample}/{wildcards.id_library}.{wildcards.id_genome}.bam"
	return (bam)


##########################################################################################
##########################################################################################

##########################################################################################
## all functions for the samples

def get_final_bam(wildcards):
    if run_compute_md:
        bam = f"results/03_sample/02_md_flag/01_md_flag/{wildcards.id_sample}.{wildcards.id_genome}.bam"
    elif run_realign:
        bam = f"results/03_sample/01_realigned/01_realign/{wildcards.id_sample}.{wildcards.id_genome}.bam"
    else:
        bam = f"results/03_sample/00_merged_library/01_bam/{wildcards.id_sample}.{wildcards.id_genome}.bam"
    return (bam)


def get_md_flag_bam(wildcards):
    if run_realign:
        bam = f"results/03_sample/01_realigned/01_realign/{wildcards.id_sample}.{wildcards.id_genome}.bam"
    else:
        bam = f"results/03_sample/00_merged_library/01_bam/{wildcards.id_sample}.{wildcards.id_genome}.bam"
    return (bam)


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

## get the stat file names to run multiqc
def get_flagstat_for_multiqc(wildcards):
    parts = pathlib.Path(wildcards.folder).parts
    if parts[1] == "01_fastq":
        if wildcards.group == "final":
            filename = [("results/{level}/04_final_fastq/01_bam/{SM}/{LB}/{ID}.{genome}_flagstat.txt").format(level=parts[1], group=wildcards.group, ID=row['ID'], SM=row['SM'], LB=row['LB'], genome=wildcards.id_genome) for index, row in db.iterrows()]
        elif wildcards.group == "mapped":
            filename = [("results/{level}/02_mapped/03_bam_sort/{SM}/{LB}/{ID}.{genome}_flagstat.txt").format(level=parts[1], group=wildcards.group, ID=row['ID'], SM=row['SM'], LB=row['LB'], genome=wildcards.id_genome) for index, row in db.iterrows()]
        else:
            print(f"ERROR: This should never happen: error in def get_flagstat_for_multiqc ({wildcards.folder}, {wildcards.group})!")
            os._exit(0)
    elif parts[1] == "02_library":
        filename = [("results/{level}/03_final_library/01_bam/{SM}/{LB}.{genome}_flagstat.txt").format(level=parts[1], group=wildcards.group, SM=row['SM'], LB=row['LB'], genome=wildcards.id_genome) for index, row in all_libraries.iterrows()]
    elif parts[1] == "03_sample":	
        filename=expand("results/{level}/03_final_sample/01_bam/{id_sample}.{id_genome}_flagstat.txt", level=parts[1], id_sample=list(samples), group=wildcards.group, id_genome=wildcards.id_genome)
    else:
        print(f"ERROR: This should never happen: error in def get_flagstat_for_multiqc ({wildcards.folder}, {wildcards.group})!")
        os._exit(0)
    return(filename)


## get the individual depth files to combien them
def get_depth_files(wildcards):
    parts = pathlib.Path(wildcards.folder).parts
    if parts[1] == "02_library":
        filename = [("results/02_library/03_final_library/01_bam/{SM}/{LB}.{genome}_depth.txt").format(SM=row['SM'], LB=row['LB'], genome=wildcards.id_genome) for index, row in db.iterrows()]
    else:
        filename = [("results/03_sample/03_final_sample/01_bam/{SM}.{genome}_depth.txt").format(SM=SM, genome=wildcards.id_genome) for SM in list(samples)]
    return(filename)


def get_fastqc_for_multiqc(wildcards):
    folder = "results/01_fastq"
    if wildcards.group == "orig":
        folder = f"{folder}/00_reads/01_files_orig"
    else:
        folder = f"{folder}/01_trimmed/01_files_trim"
    fastqc = [("{folder}/{SM}/{LB}/{ID}_fastqc.zip").format(folder=folder, ID=row['ID'], SM=row['SM'], LB=row['LB']) for index, row in db.iterrows()]
    return (fastqc)


   
