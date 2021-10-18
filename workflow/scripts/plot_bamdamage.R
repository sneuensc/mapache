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
        --length                                     - csv to plot
        --five_prime                                 - csv to plot
        --three_prime                                - csv to plot
        --sample=ind1                                - name of the sample
        --sample=lib1                                - name of the library
        --genome=genome                              - name of the genome
        
        Output files:  
        --length_svg
        --damage_svg
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

args = commandArgs(trailingOnly=TRUE)

#path=args[1]
#genome=args[2]
#sample=args[3]
#lib=args[4]

length_csv = get_args(argsL, "length")
five_prime_csv = get_args(argsL, "five_prime")
three_prime_csv = get_args(argsL, "three_prime")
genome = get_args(argsL, "genome")
sample = get_args(argsL, "sample")
lib = get_args(argsL, "library")
length_svg = get_args(argsL, "length_svg")
damage_svg = get_args(argsL, "damage_svg")
# length<-read.csv(paste0(path,"/",lib,".",genome,".length.csv"))
length <- read.csv(length_csv)

# p5<-read.csv(paste0(path,"/",lib,".",genome,".dam_5prime.csv"), row.names=1)
# p3<-read.csv(paste0(path,"/",lib,".",genome,".dam_3prime.csv"), row.names=1)

p5 <- read.csv(five_prime_csv)
p3 <- read.csv(three_prime_csv)

############################################################################
library(ggplot2)
library(RColorBrewer) 

############################################################################
#svg(paste0(path,"/",lib,".",genome,".length.svg"), width = 11, height = 7)
svg(length_svg, width = 11, height = 7)

length$frequency<-length$counts/sum(length$counts)
ggplot(data=length, mapping=aes(x = length, y=frequency)) + 
  geom_bar(stat="identity") +
  ylab("frequency") +
  xlab("read length") +
  theme_classic() +
  ggtitle(paste0("Read length distribution of '", lib, "' (sample: ", sample, ")"))

dev.off()

############################################################################
#svg(paste0(path,"/",lib,".",genome,".dam.svg"), width = 11, height = 7)
svg(damage_svg, width = 11, height = 7)

par(mfrow=c(1,2)) 
par(oma=c(0,0,0,0)) 
par(mar=c(5,4,4,0))
yticks<-pretty(c(0,max(p3[-1],p5[-1])))
xticks<-pretty(1:nrow(p5))
plot(1:nrow(p5), p5[,1], ylim=c(0,max(p3,p5)), col.main="orange",
     xlab="position from 5' end", ylab="", axes=F,type="n")
axis(1, at = xticks, labels = xticks, tick = TRUE)
axis(2, at = yticks, labels = yticks, tick = TRUE)
mtext("frequency",side=2,col="black",line=2.5)
title("C>T", line = -1, col.main="orange")
for(s in 1:ncol(p5)){
  lines(1:nrow(p5), p5[,s], col="gray", lwd=2)
}
lines(1:nrow(p5), p5[,2], col="blue", lwd=4)
lines(1:nrow(p5), p5[,11], col="orange", lwd=4)

par(mar=c(5,0,4,4))
plot(rev(1:nrow(p3)), p3[,1], ylim=c(0,max(p3,p5)), col.main="blue",
     xlab="position from 3' end", ylab="", axes=F,  type="n")
axis(1, at = xticks, labels = rev(xticks), tick = TRUE)
axis(4, at = yticks, labels = yticks, tick = TRUE)
mtext("frequency",side=4,col="black",line=2.5)
title("G>A", line = -1, col.main="blue")
for(s in 1:ncol(p3)){
  lines(rev(1:nrow(p3)),p3[,s], col="gray", lwd=2)
}
lines(rev(1:nrow(p3)), p3[,11], col="orange", lwd=4)
lines(rev(1:nrow(p3)), p3[,2], col="blue", lwd=4)

mtext("Damage pattern", outer = TRUE, cex = 1.5, line = -1.2)
mtext(sample, outer = TRUE, cex = 1, line = -2.5)
mtext(lib, outer = TRUE, cex = 1, line = -3.8)


dev.off()

