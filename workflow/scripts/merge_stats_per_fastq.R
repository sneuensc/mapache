        #require(plyr)

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
# path_multiqc_orig = "results/04_stats/01_fastq/00_reads/01_files_orig/multiqc_fastqc_data/multiqc_fastqc.txt"  # raw sequenced reads
# path_multiqc_trim = "results/04_stats/01_fastq/01_trimmed/01_files_trim/multiqc_fastqc_data/multiqc_fastqc.txt" # raw trimmed reads
# path_flagstat_mapped = "results/04_stats/01_fastq/02_mapped/03_bam_sort/ind1/lib1_lb/lib1_R1_002_fq.hg19_flagstat.txt"      # mapped and low-qual reads
# path_flagstat_mapped_highQ = "results/04_stats/01_fastq/04_final_fastq/01_bam/ind1/lib1_lb/lib1_R1_002_fq.hg19_flagstat.txt"       # mapped and high-qual reads
# path_length_mapped_highQ = "results/04_stats/01_fastq/04_final_fastq/01_bam/ind1/lib1_lb/lib1_R1_002_fq.hg19.length"
ID = get_args(argsL, "ID")
LB = get_args(argsL, "LB")
SM = get_args(argsL, "SM")
genome = get_args(argsL, "genome")
output_file = get_args(argsL, "output_file")
path_multiqc_orig = get_args(argsL, "path_multiqc_orig")
path_multiqc_trim = get_args(argsL, "path_multiqc_trim")
path_flagstat_mapped_highQ = get_args(argsL, "path_flagstat_mapped_highQ")
path_length_mapped_highQ = get_args(argsL, "path_length_mapped_highQ")

#-----------------------------------------------------------------------------#



parse_ids <- function(string, column, replace = F){
    words <- rev(strsplit(string, " \\| ")[[1]])
    if(column == "ID"){
        id <- words[1]
    }else if(column == "LB"){
        id <- words[2]
    }else if(column == "SM"){
        id <- words[3]
    }
    
    return(id)
}

calc_avg_len <- function(l){ sum(l$n_reads * l$length) / sum(l$n_reads) }
#-----------------------------------------------------------------------------#
multiqc_orig = read.table(path_multiqc_orig, sep = "\t", header = T)
multiqc_trim = read.table(path_multiqc_trim, sep = "\t", header = T)
length_mapped_highQ = read.table(path_length_mapped_highQ, sep = "\t", header = T)
mapped_raw = as.numeric(strsplit(readLines(path_flagstat_mapped_highQ)[1], " ")[[1]][1])
        
#-----------------------------------------------------------------------------#
for(column in c("SM", "LB", "ID")){
    multiqc_orig[, column] <- sapply(multiqc_orig$Sample, function(string) parse_ids(string, column))
    multiqc_trim[, column] <- sapply(multiqc_trim$Sample, function(string) parse_ids(string, column))
}

index = multiqc_orig$SM == SM & multiqc_orig$LB == LB & multiqc_orig$ID == ID
reads_raw =  multiqc_orig$Total.Sequences[index]
length_reads_raw = multiqc_orig$avg_sequence_length[index]

index = multiqc_trim$SM == SM & multiqc_trim$LB == LB & multiqc_trim$ID == ID
reads_trim =  multiqc_trim$Total.Sequences[index]
length_trim = multiqc_trim$avg_sequence_length[index]

trim_prop = reads_trim / reads_raw
endogenous_raw = mapped_raw / reads_raw
length_mapped_raw = calc_avg_len(length_mapped_highQ)
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