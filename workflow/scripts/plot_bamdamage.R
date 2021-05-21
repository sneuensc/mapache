args = commandArgs(trailingOnly=TRUE)

path=args[1]
genome=args[2]
sample=args[3]
lib=args[4]

length<-read.csv(paste0(path,"/",lib,".",genome,".length.csv"))

p5<-read.csv(paste0(path,"/",lib,".",genome,".dam_5prime.csv"), row.names=1)
p3<-read.csv(paste0(path,"/",lib,".",genome,".dam_3prime.csv"), row.names=1)


library(ggplot2)
library(RColorBrewer) 

############################################################################
svg(paste0(path,"/",lib,".",genome,".length.svg"), width = 11, height = 7)

length$frequency<-length$counts/sum(length$counts)
ggplot(data=length, mapping=aes(x = length, y=frequency)) + 
  geom_bar(stat="identity") +
  ylab("frequency") +
  xlab("read length") +
  theme_classic() +
  ggtitle(paste0("Read length distribution of '", lib, "' (sample: ", sample, ")"))

dev.off()

############################################################################
svg(paste0(path,"/",lib,".",genome,".dam.svg"), width = 11, height = 7)

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

