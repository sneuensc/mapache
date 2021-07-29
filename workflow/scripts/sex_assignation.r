# - XY
#   - X is larger than Y
#   - Males: X/genome = 0.5
#   - Females: X/genome = 1
# - ZW
#   - Z is larger than W
#   - Males: Z/genome = 1
#   - Females: Z/genome = 0.5
# - XO
#   - There is only one sex chromosome (X)
#   - Males: X/genome = 0.5
#   - Females: X/genome = 1
# - ZO
#   - There is only one sex chromosome (Z)
#   - Males: Z/genome = 1
#   - Females: Z/genome = 0.5

sex_system <- list(
  XY = list(
    larger_chr = "X",
    mappability_correction = 1,
    sex_with_larger_ratio = "Female",
    sex_to_autosomes_second_sex = 0.5,
    name_second_sex = "Male"  
  ),
  ZW = list(
    larger_chr = "Z",
    mappability_correction = 1,
    sex_with_larger_ratio = "Male",
    sex_to_autosomes_second_sex = 0.5,
    name_second_sex = "Female"  
  ),
  XO = list(
    larger_chr = "X",
    mappability_correction = 1,
    sex_with_larger_ratio = "Female",
    sex_to_autosomes_second_sex = 0.5,
    name_second_sex = "Male"  
  ),
  ZO = list(
    larger_chr = "Z",
    mappability_correction = 1,
    sex_with_larger_ratio = "Male",
    sex_to_autosomes_second_sex = 0.5,
    name_second_sex = "Female"  
  )
)

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      Sex determination
 
      Arguments:
      --genomecov=file.genomecov          - input file from the command 'bedtools genomecov -ibam file.bam > file.genomecov'
      --autosomeChr=1,2,..,22             - character, comma-separated list (chr1,chr2,chr3). Default: '1,2,..,22'
      --SM=sample                         - character, sample ID. Default: NA
      --LB=library                        - character, library ID. Default: NA
      --ID=id                             - character, FASTQ ID. Default: NA
      --out=sex.out                       - character, name of output file. Default: sex.out
      --system=XY                         - character, either of XY, ZW, XO or ZO. Default: XY
      --larger_chr=X                      - character, larger sex chromosome. Default: X
      --mappability_correction=1          - numeric, value to account for different mappabilities between sex and autosomal chromosomes. Default: 1 (no correction)
      --sex_with_larger_ratio=Female      - character, name assigned to the sex with the largest sex chromosome. Default: Female
      --sex_to_autosomes_second_sex=0.5   - numeric, expected DoC ratio between sex chromosome and autosomes in an organism of the second sex. Default: 0.5
      --name_second_sex=Male              - character, name assigned to the other sex
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
  }else{
    value = default
  }
  return(value)
}


# getting arguments
system                  = get_args(argsL, "system", "XY")
larger_chr              = get_args(argsL, "larger_chr", sex_system[[system]]$larger_chr)
mappability_correction  = get_args(argsL, "mappability_correction", sex_system[[system]]$mappability_correction)
sex_with_larger_ratio   = get_args(argsL, "sex_with_larger_ratio", sex_system[[system]]$sex_with_larger_ratio)
sex_to_autosomes        = get_args(argsL, "sex_to_autosomes", sex_system[[system]]$sex_to_autosomes)
name_second_sex         = get_args(argsL, "name_second_sex", sex_system[[system]]$name_second_sex)
autosomes               = get_args(argsL, "autosomeChr", paste0(1:22, collapse=','))

output_file             = get_args(argsL, "out", "sex.out")
SM                      = get_args(argsL, "SM", NA)
LB                      = get_args(argsL, "LB", NA)
ID                      = get_args(argsL, "ID", NA)

# get coverage on autosomes and sex chromosomes
genomecov <- read.table(argsL$genomecov)

# comma separated string to vector
autosomes <- unlist(strsplit(autosomes, ','))

## remove chromosomes which are not present


colnames(genomecov) <- c("chr", "depth", "counts", "length", "freq")
autosomes_length <- sum(
  sapply(
    autosomes, 
    function(chr) unique(genomecov$length[genomecov$chr == chr])
  )
)

genomecov_larger_chr <- genomecov[genomecov$chr == larger_chr,]
genomecov <- genomecov[genomecov$chr %in% autosomes,]

cov_autosomes <-  sum(genomecov$depth * genomecov$counts) / autosomes_length
cov_larger_chromosome <-  sum(genomecov_larger_chr$depth * genomecov_larger_chr$freq)


# likelihood for females (XY system)
probs <- dpois(
  genomecov_larger_chr$depth, 
  lambda = cov_autosomes * mappability_correction, 
  log = T
)
ll_first_sex <-  sum(probs * genomecov_larger_chr$counts)

# likelihood for males (XY system)
probs <- dpois(genomecov_larger_chr$depth, 
               lambda = cov_autosomes * sex_to_autosomes, 
               log = T)
ll_second_sex <-  sum(probs* genomecov_larger_chr$counts)

# assign sex
if(ll_first_sex > ll_second_sex){
  sex <- sex_system[[system]]$sex_with_larger_ratio
}else{
  sex <- sex_system[[system]]$name_second_sex
}

# make data frame with final info and save
sex_data <- data.frame(
  "SM" = SM,
  "LB" = LB,
  "ID" = ID,
  "Sex" = sex,
  "DoC_autosomes" = cov_autosomes,
  "DoC_larger_chr" = cov_larger_chromosome,
  "LogLikelihood_first_sex" = ll_first_sex,
  "LogLikelihood_second_sex" = ll_second_sex
)


for(column in c("SM", "LB", "ID")){
  if(is.na(sex_data[,column])){
    sex_data[,column] <- NULL
  }
}

colnames(sex_data) <- sub("larger_chr", sex_system[[system]]$larger_chr, colnames(sex_data))
colnames(sex_data) <- sub("first_sex", sex_system[[system]]$sex_with_larger_ratio, colnames(sex_data))
colnames(sex_data) <- sub("second_sex", sex_system[[system]]$name_second_sex, colnames(sex_data))
# print(sex_data)
write.csv(sex_data, output_file, row.names = F)