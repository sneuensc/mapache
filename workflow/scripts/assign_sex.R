# Reimplementation of Rx_identifier originally published by Mittnik et al., 2016
# Generalized to determine sex for organisms with different sex systems, 
# according to the thresholds specified by the user.
# Author: Diana I. Cruz Davalos

#-----------------------------------------------------------------------------#
# For each system, we specify:
# the sex chromosome which will be used for the inference
# the ratio of that sex chromosome to the autosomes
# the sex assigned according to this ratio 
# the thresholds to assign sex 

sex_system <- list(
  XY = list(
    sex_chr = "X",
    ratio = 1,
    sex = "Female",
    ratio_second_sex = 0.5,
    name_second_sex = "Male",
    thresholds = list(
      "XX" = c(0.8, 1),
      "XY" = c(0, 0.6),
      "consistent with XX but not XY" = c(0.6, 1),
      "consistent with XY but not XX" = c(0, 0.8)
    )
  ),
  ZW = list(
    sex_chr = "Z",
    ratio = 1,
    sex = "Male",
    ratio_second_sex = 0.5,
    name_second_sex = "Female",
    thresholds = list(
      "ZZ" = c(0.8, 1),
      "ZW" = c(0, 0.6),
      "consistent with ZZ but not ZW" = c(0.6, 1),
      "consistent with ZW but not ZZ" = c(0, 0.8)
    )
  ),
  XO = list(
    sex_chr = "X",
    ratio = 1,
    sex = "Female",
    ratio_second_sex = 0.5,
    name_second_sex = "Male",
    thresholds = list(
      "XO" = c(0.4, 0.5),
      "OO" = c(0, 0.2),
      "consistent with XO but not OO" = c(0.2, 0.5),
      "consistent with OO but not XO" = c(0, 0.4)
    )
  ),
  ZO = list(
    sex_chr = "Z",
    ratio = 1,
    sex = "Male",
    ratio_second_sex = 0.5,
    second_sex = "Female",
    thresholds = list(
      "ZO" = c(0.4, 0.5),
      "OO" = c(0, 0.2),
      "consistent with ZO but not OO" = c(0.2, 0.5),
      "consistent with ZZ but not ZO" = c(0, 0.4)
    )
  )
)

#-----------------------------------------------------------------------------#
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
      --idxstats=file.idxstats          - input file that can be obtained with 'samtools idxstats input.bam > file.idxstats'
      --autosomes=1:22             - character, quoted R expression. Examples: '1:22'; 'seq(1,22)'; 'paste0(\"chr\", 1:22)' 
      --signif=0.95
      --thresholds
      --SM=sample                         - character, sample ID. Default: NA
      --LB=library                        - character, library ID. Default: NA
      --ID=id                             - character, FASTQ ID. Default: NA
      --out=sex.out                       - character, name of output file. Default: sex.out
      --system=XY                         - character, either of XY, ZW, XO or ZO. Default: XY
      --sex_chr=X                         - character, larger sex chromosome. Default: X
      --help                              - print this text
 
      Example:
      Rscript assign_sex.R --idxstats=file.idxstats [] \n\n")
  
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(all_args){
  lapply(
    all_args, 
    function(str) {
      # removing -- from the name; this will go to the first column of argsDF
      new_str = sub("^--", "", str)
      # split string --x=y into x and y, on the first = sign
      # this works also for expressions like
      # --thresholds='list( "XX"=c(0.8, 1.3), "XY"=c(0, 0.6) )'
      # returning
      # thresholds
      # 'list( "XX"=c(0.8, 1.3), "XY"=c(0, 0.6) )'     which can be parsed later
      regmatches(new_str, regexpr("=", new_str), invert = TRUE)[[1]]
    }
    )
} 

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1


get_args <- function(argsL, name, default){
  if(name %in% names(argsL)){
    value = argsL[[name]]
  }else{
    value = default
  }
  return(value)
}


# getting arguments
idxstats_input          = get_args(argsL, "idxstats", NA)
system                  = get_args(argsL, "system", "XY")
sex_chr                 = get_args(argsL, "sex_chr", sex_system[[system]]$sex_chr)
autosomes_expression    = get_args(argsL, "autosomes", "1:22")
thresholds              = get_args(argsL, "thresholds", sex_system[[system]]$thresholds)
if("thresholds" %in% argsDF[,1]){

          thresholds =                       eval(parse(text=gsub("\\\\", "", thresholds)))
} 
signif                  = as.numeric(get_args(argsL, "signif", 0.95))

output_file             = get_args(argsL, "out", "sex.out")
SM                      = get_args(argsL, "SM", NA)
LB                      = get_args(argsL, "LB", NA)
ID                      = get_args(argsL, "ID", NA)
#-----------------------------------------------------------------------------#

# comma separated string to vector
autosomes <- eval(parse(text = autosomes_expression))

n_autosomes = length(autosomes)

idxstats <- read.table(
                    idxstats_input, 
                    header=F,
                    col.names=c("chr", "length", "reads", "unmapped")
                    )


idxstats = idxstats[idxstats$chr %in% c(autosomes, sex_chr), ]

if(sum(idxstats$reads == 0)){
# do not assign sex if there are chromosomes without reads mapping

  Sex <- "Not assigned - chromosomes without reads"
  Rx <- NA
  confinterval <- NA
  signif <- NA
  reads_autosomes <- NA
  reads_sex_chr <- NA 

}else{

  genome_length <- sum(idxstats$length)
  total_mapped <- sum(idxstats$reads)
  
  LM <- lm(reads ~ length, idxstats)
  # corresponding to f in Mittnik et al., 2016
  idxstats$fraction_mapped <- idxstats$reads / total_mapped

  idxstats$fraction_of_genome <- idxstats$length / genome_length

  # if the mappability is the same across chromosomes,
  # you expect this ratio to be 1 for every chromosome
  # if a chromosome is present in 1 copy instead of 2, this ratio would be 0.5
  # corresponding to rho in Mittnik et al., 2016
  idxstats$ratio_mappings <- idxstats$fraction_mapped / idxstats$fraction_of_genome

  sex_to_autosomes = idxstats$ratio_mappings[idxstats$chr == sex_chr] / idxstats$ratio_mappings[idxstats$chr != sex_chr]
  Rx = idxstats$ratio_mappings[idxstats$chr == sex_chr] / mean(idxstats$ratio_mappings[idxstats$chr != sex_chr])

  z <- qnorm(1 - ((1-signif)/2))
  confinterval <- z * (sd(sex_to_autosomes) / sqrt(n_autosomes))
  CI1 <- Rx-confinterval
  CI2 <- Rx+confinterval

  Sex <- "Not assigned"
  # loop through options for sex
  for(option in names(thresholds)){
    boundaries <- thresholds[[option]]
    if( CI1 > boundaries[1] & CI2 < boundaries[2] ){
      Sex <- option
      break
    }
  }

  reads_autosomes <- sum(idxstats$reads[idxstats$chr != sex_chr])
  reads_sex_chr <- idxstats$reads[idxstats$chr == sex_chr]
}

# make data frame with final info and save
sex_data <- data.frame(
  "SM" = SM,
  "LB" = LB,
  "ID" = ID,
  "Sex" = Sex,
  "Rx" = Rx,
  "CI" = confinterval,
  "signif_set" = signif,
  "reads_autosomes" = reads_autosomes,
  "reads_sex_chr" = reads_sex_chr
)


for(column in c("SM", "LB", "ID")){
  if(is.na(sex_data[,column])){
    sex_data[,column] <- NULL
  }
}

colnames(sex_data) <- sub("sex_chr", sex_system[[system]]$sex_chr, colnames(sex_data))

#print(sex_data)
write.csv(sex_data, output_file, row.names = F)
