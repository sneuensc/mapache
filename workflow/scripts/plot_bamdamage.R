############################################################################
# to parse arguments
args <- commandArgs(TRUE)
 
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      PLotting summary statistics
 
      Arguments:
        Input files:
        --length                                     - read length csv filename
        --five_prime                                 - substitution at 5' csv filename
        --three_prime                                - substitution at 3' csv filename
        --sample                                     - name of the sample
        --library                                    - name of the library
        --genome                                     - name of the genome
        
        Output files:  
        --length_svg                                 - read length plot filename
        --damage_svg                                 - damage plot filename
        --help                                       - print this text
 
      Example:
      Rscript plot_stats.R [] \n\n")
 
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

args = commandArgs(trailingOnly=TRUE)

length_csv = get_args(argsL, "length")
five_prime_csv = get_args(argsL, "five_prime")
three_prime_csv = get_args(argsL, "three_prime")
length_svg = get_args(argsL, "length_svg")
damage_svg = get_args(argsL, "damage_svg")

sample = get_args(argsL, "sample", "")
library = get_args(argsL, "library", "")
genome = get_args(argsL, "genome", "")

############################################################################
library(ggplot2)
library(RColorBrewer) 
library(reshape2)
library(gridExtra)

############################################################################
## plot length histogram

length <- read.csv(length_csv)

length$frequency<-length$counts/sum(length$counts)
my_plot <- ggplot(data=length, mapping=aes(x = length, y=frequency)) + 
  geom_bar(stat="identity") +
  ylab("frequency") +
  xlab("read length") +
  theme_classic() +
  ggtitle(paste0("Read length distribution of '", library, "' (sample: ", sample, ")"))

ggsave(length_svg, my_plot, width = 11, height = 7)

############################################################################
## plot damage pattern

##  read files
p5 <- read.csv(five_prime_csv)
p3 <- read.csv(three_prime_csv)
range.y <- c(0,max(p3[-1],p5[-1]))

## 5' plot
d <- melt(p5, 1)
p1 <- ggplot(d, aes(x=X, y=value, group=variable)) +
  geom_line(size=1, color="gray") +
  geom_line(data=subset(d, variable == 'G..A'), size=1.5, color="blue") +
  geom_line(data=subset(d, variable == 'C..T'), size=1.5, color="orange") +
  ylim(range.y) +
  xlab("position from 5' end") +
  ylab("frequency") +
  geom_text(x=mean(p5$X), y=0.95*max(range.y), label="C>T", color="orange", fontface="bold") +
  theme_classic() 

## 3' plot
d <- melt(p3, 1)
p2 <- ggplot(d, aes(x=X, y=value, group=variable)) +
  geom_line(size=1, color="gray") +
  geom_line(data=subset(d, variable == 'C..T'), size=1.5, color="orange") +
  geom_line(data=subset(d, variable == 'G..A'), size=1.5, color="blue") +
  xlab("position from 3' end") +
  ylab("frequency") +
  geom_text(x=-mean(p3$X), y=0.95*max(range.y), label="G>A", color="blue", fontface="bold") +
  theme_classic() +
  scale_y_continuous(lim=range.y, position = "right") +
  scale_x_reverse()

## combine plots
my_plot <- grid.arrange(grobs=list(p1, p2), ncol=2, top = paste("Damage pattern", sample, library, sep="\n"))   
ggsave(damage_svg, my_plot, width = 11, height = 7)

