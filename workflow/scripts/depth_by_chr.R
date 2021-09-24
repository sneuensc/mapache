
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
      --path_genomecov=file.genomecov          - character, genomecov file. Default: NA
      --SM=sample                         - character, sample ID. Default: NA
      --output_file=out.txt                       - character, output file name.
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
path_genomecov = get_args(argsL, "path_genomecov")
output_file = get_args(argsL, "output_file")

genomecov = read.table(path_genomecov, header = F, sep = "\t")

colnames(genomecov) =  c("chr", "depth", "counts", "length", "frac")

calc_DoC <- function(genomecov, chr){
    genomecov = genomecov[genomecov$chr == chr,]
    chr_length = unique(genomecov$length)
    DoC = sum(genomecov$depth * genomecov$counts) / chr_length
    return(DoC)
}

chromosomes <- unique(genomecov$chr)

DoC_chrs <- cbind(
    data.frame(SM = SM),
    as.data.frame(t(sapply(chromosomes, calc_DoC, genomecov=genomecov)))
    )

write.csv(DoC_chrs, output_file,  row.names = F, quote = F)

