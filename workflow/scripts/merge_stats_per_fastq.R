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
      --SM=sample                         - character, sample ID. Default: NA
      --LB=library                        - character, library ID. Default: NA
      --ID=id                             - character, FASTQ ID. Default: NA
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
# path_flagstat_mapped_highQ = "results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/ind1/lib1_lb/lib1_R1_002_fq.hg19_flagstat.txt"       # mapped and high-qual reads
# path_length_mapped_highQ = "results/04_stats/01_sparse_stats/01_fastq/04_final_fastq/01_bam/ind1/lib1_lb/lib1_R1_002_fq.hg19.length"

ID = get_args(argsL, "ID")
LB = get_args(argsL, "LB")
SM = get_args(argsL, "SM")
genome = get_args(argsL, "genome")
output_file = get_args(argsL, "output_file")
#paired = get_args(argsL, FALSE)
#collapsed = get_args(argsL, FALSE)
#path_adapterremoval = get_args(argsL, "path_adapterremoval")
path_fastqc_orig = get_args(argsL, "path_fastqc_orig")
path_fastqc_trim = get_args(argsL, "path_fastqc_trim")
path_flagstat_mapped_highQ = get_args(argsL, "path_flagstat_mapped_highQ")
path_length_mapped_highQ = get_args(argsL, "path_length_mapped_highQ")
script_parse_fastqc = get_args(argsL, "script_parse_fastqc")

#-----------------------------------------------------------------------------#

calc_avg_len <- function(l){ sum(l$Count * l$Length) / sum(l$Count) }

#-----------------------------------------------------------------------------#
length_mapped_highQ = read.table(path_length_mapped_highQ, sep = "\t", header = T, col.names = c("Count", "Length"))
mapped_raw = as.numeric(strsplit(readLines(path_flagstat_mapped_highQ)[1], " ")[[1]][1])
        
#-----------------------------------------------------------------------------#
## original fastqc
#options(readr.show_col_types = FALSE)
source(script_parse_fastqc)

# lines_to_skip <-  as.integer(
#     system(
#         paste0(
#             "grep -n 'Length distribution' ",
#             path_adapterremoval,
#             " | cut -f 1 -d:"
#         ),
#         intern = T
#     )
# )

# length_dist_trimmed <- read.table(path_adapterremoval, skip = lines_to_skip, header = T)

length_dist_raw <- parse_fastqc(file_name = path_fastqc_orig, df_name = "Sequence Length Distribution") 
length_dist_raw <- adjust_lengths(length_dist_raw)
reads_raw <- sum(length_dist_raw$Count)
# raw reads

# grep_adapterremoval <- function(line){
#     as.integer(
#             system(
#                 paste0("grep ", line, path_adapterremoval, "|rev|cut -f1 -d ' ' |rev"), 
#                 intern = T
#             )
#     )
# }

# if(paired){
#     line <- 'Total number of read pairs:'
# }else{
#     line <- 'Total number of reads:'
# }

# reads_raw <- grep_adapterremoval(line)


length_reads_raw <- calc_avg_len(length_dist_raw)

## trimmed fastqc
if(path_fastqc_trim == "Not_trimmed"){
    reads_trim <- 0
    length_reads_trimmed <- 0
}else{
    length_dist_trimmed <- parse_fastqc(file_name = path_fastqc_trim, df_name = "Sequence Length Distribution")
    length_dist_trimmed <- adjust_lengths(length_dist_trimmed)
    length_reads_trimmed <- calc_avg_len(length_dist_trimmed)
    reads_trim <- sum(length_dist_trimmed$Count)
}



# if(paired){
#     if(collapsed){
#         line <- 'Number of full-length collapsed pairs:'
#         reads_trim <- grep_adapterremoval(line)

#         length_reads_trimmed <- sum(
#             length_dist_trimmed$Collapsed * length_dist_trimmed$Length
#         ) / sum(length_dist_trimmed$Collapsed)

#     }else{
#         well_aligned <- grep_adapterremoval('Number of well aligned read pairs:')
        
#         unaligned <- grep_adapterremoval('Number of unaligned read pairs:')

#         discarded <- grep_adapterremoval('Number of discarded mate 1 reads:') +
#             grep_adapterremoval('Number of discarded mate 2 reads:')

#         singleton <- grep_adapterremoval('Number of singleton mate 1 reads:') +
#             grep_adapterremoval('Number of singleton mate 2 reads:')
        
#         reads_trim <- well_aligned + unaligned - discarded - singleton

#         length_reads_trimmed <- sum(
#             length_dist_trimmed$Length * 
#             (length_dist_trimmed$Mate1 + length_dist_trimmed$Mate2)
#             ) / (2 * reads_trim)
#     }
# }else{
#     line <- 'Number of retained reads:'
#     reads_trim <- grep_adapterremoval(line)

#     length_reads_trimmed <- sum(
#             length_dist_trimmed$Length * 
#             length_dist_trimmed$Mate1 
#             ) / reads_trim
# }

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
