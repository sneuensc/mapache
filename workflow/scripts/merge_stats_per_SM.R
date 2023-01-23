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
      --sm=sample                 - character, sample ID. Default: NA

      --path_list_stats_lb        - e.g., results/04_stats/02_separate_tables/hg19/ind1/lib1_lb/library_stats.csv,results/04_stats/02_separate_tables/hg19/ind1/lib2_lb/library_stats.csv
      --path_length_unique        - e.g., results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/ind1.hg19_length.txt
      --path_idxstats_unique      - e.g., results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/ind1.hg19_idxstats.txt
      --path_sex_unique           - e.g., results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/ind1.hg19_sex.txt
      --chrs_selected             - e.g., X,Y,MT

      --output_file               - e.g., results/04_stats/02_separate_tables/hg19/ind1/sample_stats.csv

      --help                      - print this text"
      )
 
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




SM = get_args(argsL, "sm")
genome = get_args(argsL, "genome")
output_file = get_args(argsL, "output_file")

path_list_stats_lb = get_args(argsL, "path_list_stats_lb")
path_length_unique = get_args(argsL, "path_length_unique")
path_idxstats_unique = get_args(argsL, "path_idxstats_unique")
path_sex_unique = get_args(argsL, "path_sex_unique")
chrs_selected = get_args(argsL, "chrs_selected", NULL)


# SM                      = "ind2"
# genome                  = "GRCh38"
# output_file             = "results/04_stats/02_separate_tables/GRCh38/ind2/stats.csv"
# path_list_stats_lb      = "results/04_stats/02_separate_tables/GRCh38/ind2/lib3_lb/stats.csv"
## path_idxstats_unique    = "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/ind2.GRCh38_idxstats.txt"  
# path_length_unique      = "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/ind2.GRCh38.length"
# path_genomecov_unique   = "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/ind2.GRCh38_genomecov"    
# path_sex_unique         = "results/04_stats/01_sparse_stats/03_sample/03_final_sample/01_bam/ind2.GRCh38_sex"
# chrs_selected           = "chrX,chrY,chrMT"
#-----------------------------------------------------------------------------#
## get number of reported chromosome depths
nb_chrs_depth = 0
if(!is.null(chrs_selected)){
    nb_chrs_depth = length(unlist(strsplit(chrs_selected, ",")))
}

## double is used instead of integer, as integers are limited in size:
##   "Note that current implementations of R use 32-bit integers for integer vectors, so the range
##    of representable integers is restricted to about +/-2*10^9: doubles can hold much larger integers exactly.""
#-----------------------------------------------------------------------------#
calc_DoC_idxstats <- function(idxstats, read_length, chr){
    chr_length <- sum(idxstats$length[idxstats$chr %in% chr])
    reads_chr <- sum(idxstats$mapped[idxstats$chr %in% chr])
    DoC <- (reads_chr * read_length) / chr_length
    return(DoC)
} 
calc_DoC <- function(genomecov, chr){
    genomecov = genomecov[genomecov$chr == chr,]
    chr_length = unique(genomecov$length)
    DoC = sum(genomecov$depth * genomecov$counts) / chr_length
    return(DoC)
}
calc_avg_len <- function(l){ sum(l$n_reads * l$length) / sum(l$n_reads) }
#-----------------------------------------------------------------------------#
 
idxstats = read.table(
    path_idxstats_unique, 
    header=F,
    col.names=c("chr", "length", "mapped", "unmapped")
    )

length_unique_table = read.table(path_length_unique, header = T, sep = "\t", colClasses = c("numeric", "numeric"))
length_mapped_unique = calc_avg_len(length_unique_table)
read_depth = calc_DoC_idxstats(idxstats = idxstats, read_length = length_mapped_unique, chr = idxstats$chr)
sex_unique = read.csv(path_sex_unique)
mapped_unique = sum(idxstats$mapped)

if(path_list_stats_lb != "path_list_stats_lb"){
    stats_lb = do.call(rbind, lapply(strsplit(path_list_stats_lb, ",")[[1]], read.csv, 
            colClasses = c(rep("character", 3), rep("numeric", 13 + nb_chrs_depth))))
    #genomecov_unique = read.table(path_genomecov_unique, header = F, sep = "\t")
    #colnames(genomecov_unique) =  c("chr", "depth", "counts", "length", "frac")
    reads_raw = sum(stats_lb$reads_raw)
    reads_trim = sum(stats_lb$reads_trim)
    trim_prop = reads_trim / reads_raw 
    mapped_raw = sum(stats_lb$mapped_raw)
    duplicates = mapped_raw - mapped_unique 
    duplicates_prop = duplicates / mapped_raw
    endogenous_raw = mapped_raw / reads_raw
    endogenous_unique = mapped_unique / reads_raw
    length_reads_raw = sum(stats_lb$length_reads_raw * stats_lb$reads_raw) / reads_raw
    length_reads_trimmed = sum(stats_lb$length_reads_trimmed * stats_lb$reads_trim) / reads_trim
    length_mapped_raw = sum(stats_lb$length_mapped_raw * stats_lb$mapped_raw) / mapped_raw
} else {    ## if external sample and thus no library information
    reads_raw = NA
    reads_trim = NA
    trim_prop = NA
    mapped_raw = NA
    duplicates = NA
    duplicates_prop = NA
    endogenous_raw = NA
    endogenous_unique = NA
    length_reads_raw = NA
    length_reads_trimmed = NA
    length_mapped_raw = NA
}


my_stats = data.frame(
    genome = genome, 
    SM = SM, 
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
    length_mapped_uniquem = length_mapped_unique,
    endogenous_raw = endogenous_raw,
    endogenous_unique = endogenous_unique,
    read_depth = read_depth,
    Sex = sex_unique$Sex,
    Sex_Rx = sex_unique$Rx,
    Sex_CI = sex_unique$CI
)

if(!is.null(chrs_selected)){
    chrs_selected = unlist(strsplit(chrs_selected, ","))
    
	## warn if some chromosomes are unknown
	unknown = setdiff(chrs_selected, unique(idxstats$chr))
	if(length(unknown)>0){
		print(paste0("Warning: The chromosome name '",unknown, "' is not valid in genome '", genome, "'. Ignoring it."))
		chrs_selected <- chrs_selected[!chrs_selected %in% unknown]
	}
	
    DoC_chrs_selected =  do.call(cbind, lapply(chrs_selected, function(chr) calc_DoC_idxstats(idxstats = idxstats, read_length = length_mapped_unique, chr = chr)  ))
    colnames(DoC_chrs_selected) = paste0("depth_", chrs_selected)
    my_stats = cbind(my_stats, DoC_chrs_selected)
}

write.csv(my_stats, output_file, row.names = F)
