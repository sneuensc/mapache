# IMPORTANT!!!
#
# If you modify this script, you need to make sure that the following
# scripts adopt the right column names
# scripts/merge_stats_per_fastq.R
# scripts/merge_stats_per_SM.R

# reads_raw           SM,LB,ID    multiqc
# reads_trim          SM,LB,ID    multiqc
# trim_prop           SM,LB,ID    
# mapping_raw         SM,LB,ID    flagstat
# endo_prop_raw       SM,LB,ID    
# duplicates          SM,LB,ID    
# duplicates_prop     SM,LB,ID
# mapping_final       SM,LB,ID    flagstat
# endo_final_prop     SM,LB,ID    
# DoC                 SM,LB,ID    genomecov
# DoC_chrs_selected   SM,LB,ID    genomecov
# DoC_chrs_all        SM,LB,ID    genomecov
# Sex                 SM,LB,ID    R script
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
      --sm=sample                         - character, sample ID. Default: NA
      --lb=library                        - character, library ID. Default: NA
      --genome='hg19'                     - character, genome name
      
      --path_list_stats_fastq             - e.g., results/04_stats/02_separate_tables/GRCh38/ind1/lib2_lb/lib2_R1_001_fq/fastq_stats.csv
      --path_stats_unique                 - e.g., results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/ind1/lib2_lb.GRCh38_stats.txt
      --path_length_unique                - e.g., results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/ind1/lib2_lb.GRCh38_length.txt
      --path_idxstats_unique              - e.g., results/04_stats/01_sparse_stats/02_library/03_final_library/01_bam/ind1/lib2_lb.GRCh38_idxstats.txt

      --output_file                       - e.g., results/04_stats/02_separate_tables/GRCh38/ind1/lib2_lb/library_stats.csv
    
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

get_args <- function(argsL, name, default=""){
    if(name %in% names(argsL)){
        value = argsL[[name]]
    }else if("default" %in% objects()){
        value = default
    }else{
        stop(paste0("Please specify ", name))
    }
    return(value)
}



LB = get_args(argsL, "lb")
SM = get_args(argsL, "sm")
genome = get_args(argsL, "genome")
output_file = get_args(argsL, "output_file")

path_list_stats_fastq = get_args(argsL, "path_list_stats_fastq")
#path_flagstat_raw = get_args(argsL, "path_flagstat_raw")
path_stats_unique = get_args(argsL, "path_stats_unique")
path_length_unique = get_args(argsL, "path_length_unique")
#path_genomecov_unique = get_args(argsL, "path_genomecov_unique")
path_idxstats_unique = get_args(argsL, "path_idxstats_unique")
chrs_selected = get_args(argsL, "chrs_selected", NA)

# LB = "lib1_lb"
# SM = "ind1"
# genome = "hg19"
# output_file = "lib.stats"
# path_list_stats_fastq   = "results/04_stats/TABLES/hg19/ind1/lib1_lb/lib1_R1_001_fq/stats.csv,results/04_stats/TABLES/hg19/ind1/lib1_lb/lib1_R1_002_fq/stats.csv"
## path_flagstat_raw       = "results/04_stats/02_library/00_merged_fastq/01_bam/ind1/lib1_lb.hg19_flagstat.txt"
# path_stats_unique       = "results/04_stats/02_library/03_final_library/01_bam/ind1/lib1_lb.hg19_stats.txt"  
# path_length_unique      = "results/04_stats/02_library/03_final_library/01_bam/ind1/lib1_lb.hg19.length"
# path_genomecov_unique   = "results/04_stats/02_library/03_final_library/01_bam/ind1/lib1_lb.hg19_genomecov"    
# chrs_selected = c("X","Y", "MT")
# chrs_selected = unlist(strsplit("X,Y,MT", ","))
#-----------------------------------------------------------------------------#

## double is used instead of integer, as integers are limited in size:
##   "Note that current implementations of R use 32-bit integers for integer vectors, so the range
##   of representable integers is restricted to about +/-2*10^9: doubles can hold much larger integers exactly.""
stats_fastq = do.call(rbind, lapply(strsplit(path_list_stats_fastq, ",")[[1]], read.csv , colClasses = c(rep("character", 4), rep("numeric",8))))
#mapped_raw = as.double(strsplit(readLines(path_flagstat_raw)[1], " ")[[1]][1])
mapped_raw   = sum(stats_fastq$mapped_raw)

## get uniquely mapped reads
file <- readLines(path_stats_unique)
mapped_unique = as.double(strsplit(grep("reads mapped:", file, value = T), "\t")[[1]][3])
is_paired = as.logical(as.double(strsplit(grep("reads paired:", file, value = T), "\t")[[1]][3]))
if(is_paired) mapped_unique = round(mapped_unique / 2)

length_unique_table = read.table(path_length_unique, header = T, sep = "\t", colClasses = c("numeric", "numeric"))
#genomecov_unique = read.table(path_genomecov_unique, header = F, sep = "\t")
#colnames(genomecov_unique) =  c("chr", "depth", "counts", "length", "frac")
idxstats = read.table(
    path_idxstats_unique, 
    header=F,
    col.names=c("chr", "length", "mapped", "unmapped")
    )

#-----------------------------------------------------------------------------#
calc_avg_len <- function(l){ sum(l$n_reads * l$length) / sum(l$n_reads) }
calc_DoC <- function(genomecov, chr){
    genomecov = genomecov[genomecov$chr == chr,]
    chr_length = unique(genomecov$length)
    DoC = sum(genomecov$depth * genomecov$counts) / chr_length
    return(DoC)
}
calc_DoC_idxstats <- function(idxstats, read_length, chr){
    chr_length <- sum(idxstats$length[idxstats$chr %in% chr])
    reads_chr <- sum(idxstats$mapped[idxstats$chr %in% chr])
    DoC <- (reads_chr * read_length) / chr_length
    return(DoC)
}
#-----------------------------------------------------------------------------#

reads_raw = sum(stats_fastq$reads_raw) # check that this column is defined in scripts/merge_stats_per_fastq.R

reads_trim = sum(stats_fastq$reads_trim) # check scripts/merge_stats_per_fastq.R
trim_prop =       reads_trim / reads_raw 

duplicates = mapped_raw - mapped_unique
duplicates_prop = duplicates / mapped_raw

endogenous_raw = mapped_raw / reads_raw
endogenous_unique = mapped_unique / reads_raw

length_reads_raw = sum(stats_fastq$length_reads_raw * stats_fastq$reads_raw) / reads_raw # check scripts/merge_stats_per_fastq.R
length_reads_trimmed = sum(stats_fastq$length_reads_trimmed * stats_fastq$reads_trim) / reads_trim # check scripts/merge_stats_per_fastq.R

length_mapped_raw = sum(stats_fastq$length_mapped_raw * stats_fastq$mapped_raw) / mapped_raw # check scripts/merge_stats_per_fastq.R
length_mapped_unique = calc_avg_len(length_unique_table)

read_depth = calc_DoC_idxstats(idxstats = idxstats, read_length = length_mapped_unique, chr = idxstats$chr)

my_stats = data.frame(
    genome = genome, SM = SM, LB = LB, 
    reads_raw = reads_raw,
    reads_trim = reads_trim,
    trim_prop = trim_prop,
    mapped_raw = mapped_raw,
    duplicates = duplicates,
    duplicates_prop = duplicates_prop,
    mapped_unique = mapped_unique,
    length_reads_raw = length_reads_raw,
    length_reads_trimmed = length_reads_trimmed,
    length_mapped_raw = length_mapped_raw,
    length_mapped_unique = length_mapped_unique,
    endogenous_raw = endogenous_raw,
    endogenous_unique = endogenous_unique
)

if(!is.na(chrs_selected)){
    chrs_selected = unlist(strsplit(chrs_selected, ","))
    DoC_chrs_selected =  do.call(cbind, lapply(chrs_selected, function(chr) calc_DoC_idxstats(idxstats = idxstats, read_length = length_mapped_unique, chr = chr)   ))
    colnames(DoC_chrs_selected) = paste0("depth_", chrs_selected)
    my_stats = cbind(my_stats, DoC_chrs_selected)
}

write.csv(my_stats, output_file, row.names = F)