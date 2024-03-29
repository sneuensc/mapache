#require(plyr)
# IMPORTANT!!!
#
# If you modify this script, you need to make sure that the following
# scripts adopt the right column names
# scripts/merge_stats_per_LB.R
# scripts/merge_stats_per_SM.R
# 

## The following assumes active Conda environment with `tzdata` installed
## This is needed to avoid the following error and thus an abort of snakemake
##    Error: Unknown TZ UTC
##    In addition: Warning message:
##    In OlsonNames() : no Olson database found
##    Execution halted
Sys.setenv("TZDIR"=paste0(Sys.getenv("CONDA_PREFIX"), "/share/zoneinfo"))

#library(fastqcr)

# reads_raw           SM,LB,ID    multiqc
# reads_trim          SM,LB,ID    multiqc
# trim_prop           SM,LB,ID    
# mapping_raw         SM,LB,ID    flagstat
# endo_prop_raw       SM,LB,ID    
# read_length         SM,LB,ID    python script



#-----------------------------------------------------------------------------#
args <- commandArgs(TRUE)
 
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      Get stats per fastq file
 
      Arguments:
      --am=sample                         - character, sample ID. Default: NA
      --lb=library                        - character, library ID. Default: NA
      --id=id                             - character, FASTQ ID. Default: NA
      --help                              - print this text
 
      Example:

      ")
 
  q(save="no")
}
 
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
#print(argsL)

get_args <- function(argsL, name){
    if(name %in% names(argsL)){
        value = argsL[[name]]
    }else{
        stop(paste0("Please specify ", name))
    }
    return(value)
}


# LB = "lib1_lb"
# SM = "ind1"
# ID = "lib1_R1_002_fq"
# genome = "hg19"
# output_file = "out.fq.stats"
# path_fastqc_orig = "results/04_stats/01_sparse_stats/01_fastq/00_reads/01_files_orig/ind1/lib1_lb/lib1_R1_002_fq_fastqc.zip"  # raw sequenced reads
# path_fastqc_trim = "results/04_stats/01_sparse_stats/01_fastq/01_trimmed/01_files_trim/ind1/lib1_lb/lib1_R1_002_fq_fastqc.zip" # raw trimmed reads
# path_stats_mapped_highQ = "results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/ind1/lib1_lb/lib1_R1_002_fq.hg19_stats.txt"       # mapped and high-qual reads
# path_length_mapped_highQ = "results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/ind1/lib1_lb/lib1_R1_002_fq.hg19.length"

ID = get_args(argsL, "id")
LB = get_args(argsL, "lb")
SM = get_args(argsL, "sm")
genome = get_args(argsL, "genome")
output_file = get_args(argsL, "output_file")
#paired = get_args(argsL, FALSE)
#collapsed = get_args(argsL, FALSE)
#path_adapterremoval = get_args(argsL, "path_adapterremoval")
path_fastqc_orig = get_args(argsL, "path_fastqc_orig")
path_fastqc_trim = get_args(argsL, "path_fastqc_trim")
path_stats_mapped_highQ = get_args(argsL, "path_stats_mapped_highQ")
path_length_mapped_highQ = get_args(argsL, "path_length_mapped_highQ")
script_parse_fastqc = get_args(argsL, "script_parse_fastqc")

#-----------------------------------------------------------------------------#

calc_avg_len <- function(l){ sum(l$Count * l$Length) / sum(l$Count) }

#-----------------------------------------------------------------------------#
length_mapped_highQ = read.table(path_length_mapped_highQ, sep = "\t", header = T, 
        col.names = c("Count", "Length"), 
        colClasses = c("numeric", "numeric"))

## get mapped reads
file <- readLines(path_stats_mapped_highQ)
mapped_raw = as.double(strsplit(grep("reads mapped:", file, value = T), "\t")[[1]][3])
is_paired = as.logical(as.double(strsplit(grep("reads paired:", file, value = T), "\t")[[1]][3]))
if(is_paired) mapped_raw = round(mapped_raw / 2)

#-----------------------------------------------------------------------------#
## original fastqc
#options(readr.show_col_types = FALSE)
source(script_parse_fastqc)


length_dist_raw <- parse_fastqc(file_name = path_fastqc_orig, df_name = "Sequence Length Distribution") 
length_dist_raw$Count <- as.double(length_dist_raw$Count)
length_dist_raw <- adjust_lengths(length_dist_raw)
reads_raw <- sum(length_dist_raw$Count)


length_reads_raw <- calc_avg_len(length_dist_raw)

## trimmed fastqc
if(path_fastqc_trim == "Not_trimmed"){
    reads_trim <- 0
    length_reads_trimmed <- 0
}else{
    length_dist_trimmed <- parse_fastqc(file_name = path_fastqc_trim, df_name = "Sequence Length Distribution")
    length_dist_trimmed <- adjust_lengths(length_dist_trimmed)
    length_reads_trimmed <- calc_avg_len(length_dist_trimmed)
    reads_trim <- sum(as.double(length_dist_trimmed$Count))
}


trim_prop = reads_trim / reads_raw
endogenous_raw = mapped_raw / reads_raw
# these are the mapped reads that passed the mapQ filter
# duplicates are still in the BAM at this point
length_mapped_raw = calc_avg_len(length_mapped_highQ)
#-----------------------------------------------------------------------------#


my_stats = data.frame(
    genome = genome, SM = SM, LB = LB, ID = ID, 
    reads_raw = reads_raw,
    reads_trim = reads_trim,
    trim_prop = trim_prop,
    mapped_raw = mapped_raw,
    length_reads_raw = length_reads_raw,
    length_reads_trimmed = length_reads_trimmed,
    length_mapped_raw = length_mapped_raw,
    endogenous_raw = endogenous_raw
)


#print(my_stats)

write.csv(my_stats, output_file, row.names = F)
