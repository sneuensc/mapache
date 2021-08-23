        #require(plyr)
        
## The following assumes active Conda environment with `tzdata` installed
## This is needed to avoid the following error and thus an abort of snakemake
##    Error: Unknown TZ UTC
##    In addition: Warning message:
##    In OlsonNames() : no Olson database found
##    Execution halted
Sys.setenv("TZDIR"=paste0(Sys.getenv("CONDA_PREFIX"), "/share/zoneinfo"))

library(fastqcr)

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
      Rscript sex_assignation.r --genomecov=file.genomecov [] \n\n")
 
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
path_fastqc_orig = get_args(argsL, "path_fastqc_orig")
path_fastqc_trim = get_args(argsL, "path_fastqc_trim")
path_flagstat_mapped_highQ = get_args(argsL, "path_flagstat_mapped_highQ")
path_length_mapped_highQ = get_args(argsL, "path_length_mapped_highQ")

#-----------------------------------------------------------------------------#



calc_avg_len <- function(m, nb){ 
	#print(m)
	d <- do.call(rbind, strsplit(m$Length, '-'))
	dd <- data.frame(apply(d, 1:2, function(x) as.numeric(as.character(x))))
	sum(apply(dd, 1, mean) * m$Count) / reads_raw
}

calc_avg_len2 <- function(l){ sum(l$n_reads * l$length) / sum(l$n_reads) }

#-----------------------------------------------------------------------------#
length_mapped_highQ = read.table(path_length_mapped_highQ, sep = "\t", header = T)
mapped_raw = as.numeric(strsplit(readLines(path_flagstat_mapped_highQ)[1], " ")[[1]][1])
        
#-----------------------------------------------------------------------------#
## original fastqc
data <- qc_read(path_fastqc_orig, "Sequence Length Distribution", F)
reads_raw <- sum(data$"sequence_length_distribution"$Count)
length_reads_raw <- calc_avg_len(data$"sequence_length_distribution", reads_raw)

## trimmed fastqc
data <- qc_read(path_fastqc_trim, "Sequence Length Distribution", F)
reads_trim <- sum(data$"sequence_length_distribution"$Count)
length_trim <- calc_avg_len(data$"sequence_length_distribution", reads_raw)

trim_prop = reads_trim / reads_raw
endogenous_raw = mapped_raw / reads_raw
length_mapped_raw = calc_avg_len2(length_mapped_highQ)
#-----------------------------------------------------------------------------#


my_stats = data.frame(
    genome = genome, SM = SM, LB = LB, ID = ID, 
    reads_raw = reads_raw,
    reads_trim = reads_trim,
    trim_prop = trim_prop,
    mapped_raw = mapped_raw,
    length_reads_raw = length_reads_raw,
    length_mapped_raw = length_mapped_raw,
    endogenous_raw = endogenous_raw
)


#print(my_stats)

write.csv(my_stats, output_file, row.names = F)