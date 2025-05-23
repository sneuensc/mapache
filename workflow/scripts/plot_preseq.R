############################################################################
# to parse arguments
args <- commandArgs(TRUE)

# args <-c('--input=/Users/sneuensc/Documents/Vital-IT/Sapfo/lucas/EUB01_381-s1-d1-d2-e1-L1_expexted_yield.txt', '--libSize=100000000', '--xmax=1e9', '--depth=/Users/sneuensc/Documents/Vital-IT/Sapfo/lucas/depth.txt')

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      PLotting summary statistics
 
      Arguments:
        --input                                    - input csv filename
        --output                                   - output svg filename
        --name                                     - name to plot
        --libSize                                  - number of reads (library size)
        --xmax                                     - plotting: max x
        --help                                     - print this text
 
      Example:
      Rscript plot_preseq.R [] \n\n")
  
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

############################################################################
## read arguments
input = get_args(argsL, "input")
output = get_args(argsL, "output", gsub('.txt', '.svg', input))
name = get_args(argsL, "name", basename(gsub('.txt', '', input)))
libSize = as.numeric(get_args(argsL, "libSize"))
xmax = as.numeric(get_args(argsL, "xmax"), 0.8)
depth = get_args(argsL, "depth", 'NA')

############################################################################
library(tidyverse)

## read preseq results
data <- read_tsv(input, show_col_types=FALSE)
if(xmax>1){ ## absolute number
  data <- data[data$TOTAL_READS<=xmax,]
} else { ## relative number
  data <- data[data$TOTAL_READS<=xmax*data$TOTAL_READS[nrow((data))],]
}

## get the expected size of 'my' library
libSizeExpt <- round(approx(data$TOTAL_READS, data$EXPECTED_DISTINCT, xout = libSize)$y)


p1 <- ggplot(data, aes(x=TOTAL_READS, y=EXPECTED_DISTINCT)) +
  geom_ribbon(aes(x = TOTAL_READS, ymin =  LOWER_0.95CI, ymax = UPPER_0.95CI),
              fill = "lightblue", alpha = 0.6) +
  geom_line() +
  theme_bw() +
  ggtitle(paste("Expected distinct reads", name, sep='\n')) +
  geom_abline(intercept=0, slope=1, color = "grey", linetype = "longdash") +
  geom_segment(x = libSize, xend = libSize, y = -Inf, yend = libSizeExpt,
               color = "blue", linetype = "dotted") +
  geom_segment(x = -Inf, xend = libSize, y = libSizeExpt, yend = libSizeExpt,
               color = "blue", linetype = "dotted") +
  geom_point(aes(x=libSize, y=libSizeExpt))


if(depth != 'NA'){
  ## add coverage
  d <- read.table(depth, skip=3, header=T)
  cov <- d$SumReadLength[d$Chrom=='Total:']/d$ChromLength[d$Chrom=='Total:']
  
  p1 <- p1 +
    scale_y_continuous(sec.axis = sec_axis(transform=~./libSizeExpt*cov, name="COVERAGE")
    ) +
    geom_hline(yintercept=libSizeExpt, color = "blue", linetype = "dotted") +
    geom_segment(x = libSize, xend = Inf, y = libSizeExpt, yend = libSizeExpt,
                 color = "blue", linetype = "dotted") +
    geom_text(aes(x=libSize+0.01*data$TOTAL_READS[nrow(data)], y=libSizeExpt, 
                  label=paste0('x: ', formatC(libSize, format = "e", digits = 2), '\ny: ', formatC(libSizeExpt, format = "e", digits = 2), ' / ', round(cov,2), 'x')), 
              hjust=0, vjust=1.1, color="blue")
  
  
} else{
  p1 <- p1 +
    geom_text(aes(x=libSize+0.01*data$TOTAL_READS[nrow(data)], y=libSizeExpt, 
                  label=paste0('x: ', formatC(libSize, format = "e", digits = 2), '\ny: ', formatC(libSizeExpt, format = "e", digits = 2))), 
              hjust=-0.1, vjust=1.1, color="blue")
}


ggsave(output, p1, width = 11, height = 7)