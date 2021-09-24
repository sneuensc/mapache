# IMPORTANT!!!
#
# If you modify this script, you need to make sure that the following
# scripts adopt the right column names
# scripts/merge_stats_per_LB.R
# scripts/merge_stats_per_fastq.R

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

get_args <- function(argsL, name, default){
    if(name %in% names(argsL)){
        value = argsL[[name]]
    }else if("default" %in% objects()){
        value = default
    }else{
        stop(paste0("Please specify ", name))
    }
    return(value)
}




SM = get_args(argsL, "SM")
genome = get_args(argsL, "genome")
output_file = get_args(argsL, "output_file")

path_list_stats_lb = get_args(argsL, "path_list_stats_fastq")
path_flagstat_unique = get_args(argsL, "path_flagstat_unique")
path_length_unique = get_args(argsL, "path_length_unique")
path_genomecov_unique = get_args(argsL, "path_genomecov_unique")
path_sex_unique = get_args(argsL, "path_sex_unique")
chrs_selected = get_args(argsL, "chrs_selected", NULL)


# SM                      = "ind2"
# genome                  = "GRCh38"
# output_file             = "results/04_stats/02_separate_tables/GRCh38/ind2/stats.csv"
# path_list_stats_lb      = "results/04_stats/02_separate_tables/GRCh38/ind2/lib3_lb/stats.csv"
# path_flagstat_unique    = "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/ind2.GRCh38_flagstat.txt"  
# path_length_unique      = "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/ind2.GRCh38.length"
# path_genomecov_unique   = "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/ind2.GRCh38.genomecov"    
# path_sex_unique         = "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/ind2.GRCh38.sex"
# chrs_selected           = "chrX,chrY,chrMT"
#-----------------------------------------------------------------------------#

stats_lb = do.call(rbind, lapply(strsplit(path_list_stats_lb, ",")[[1]], read.csv ))
# mapped_raw   = sum(stats_fastq$mapped_raw) # should give the same
mapped_unique = as.numeric(strsplit(readLines(path_flagstat_unique)[1], " ")[[1]][1])
length_unique_table = read.table(path_length_unique, header = T, sep = "\t")
genomecov_unique = read.table(path_genomecov_unique, header = F, sep = "\t")
colnames(genomecov_unique) =  c("chr", "depth", "counts", "length", "frac")
sex_unique = read.csv(path_sex_unique)
#-----------------------------------------------------------------------------#
calc_avg_len <- function(l){ sum(l$n_reads * l$length) / sum(l$n_reads) }
calc_DoC <- function(genomecov, chr){
    genomecov = genomecov[genomecov$chr == chr,]
    chr_length = unique(genomecov$length)
    DoC = sum(genomecov$depth * genomecov$counts) / chr_length
    return(DoC)
}
#-----------------------------------------------------------------------------#

reads_raw = sum(stats_lb$reads_raw)

reads_trim = sum(stats_lb$reads_trim)
trim_prop =       reads_trim / reads_raw 

mapped_raw = sum(stats_lb$mapped_raw)

duplicates = mapped_raw - mapped_unique 
duplicates_prop = duplicates / mapped_raw

endogenous_raw = mapped_raw / reads_raw
endogenous_unique = mapped_unique / reads_raw

length_reads_raw = sum(stats_lb$length_reads_raw * stats_lb$reads_raw) / reads_raw
length_reads_trimmed = sum(stats_lb$length_reads_trimmed * stats_lb$reads_trim) / reads_trim

length_mapped_raw = sum(stats_lb$length_mapped_raw * stats_lb$mapped_raw) / mapped_raw
length_mapped_unique = calc_avg_len(length_unique_table)

read_depth = calc_DoC(genomecov_unique, "genome")

Sex = sex_unique$Sex

my_stats = data.frame(
    genome = genome, SM = SM, 
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
    length_mapped_unique,
    endogenous_raw = endogenous_raw,
    endogenous_unique = endogenous_unique,
    Sex = Sex,
    read_depth = read_depth
)
if(!is.null(chrs_selected)){
    chrs_selected = unlist(strsplit(chrs_selected, ","))
    
	## warn if some chromosomes are unknown
	unknown = setdiff(chrs_selected, unique(genomecov_unique$chr))
	if(length(unknown)>0){
		print(paste0("Warning: The chromosome name '",unknown, "' is not valid in genome '", genome, "'. Ignoring it."))
		chrs_selected <- chrs_selected[!chrs_selected %in% unknown]
	}
	
    DoC_chrs_selected =  do.call(cbind, lapply(chrs_selected, function(chr) calc_DoC(genomecov_unique,chr)  ))
    colnames(DoC_chrs_selected) = paste0("depth_", chrs_selected)
    my_stats = cbind(my_stats, DoC_chrs_selected)
}

write.csv(my_stats, output_file, row.names = F)
