library(ggplot2)
library(RColorBrewer)
library(reshape)

## input
fq <- read.csv(snakemake@input[['fastq_stats']])
lb <- read.csv(snakemake@input[['library_stats']])
sm <- read.csv(snakemake@input[['sample_stats']])

## output
plot_1_nb_reads <- snakemake@output[['plot_1_nb_reads']]
plot_2_mapped  <- snakemake@output[['plot_2_mapped']]
plot_3_endogenous  <- snakemake@output[['plot_3_endogenous']]
plot_4_duplication  <- snakemake@output[['plot_4_duplication']]

############################################################################

nb = sum(is.na(lb$reads_raw))
title <- "Total number of reads"
if(nb)title <- paste0(title, " (", nb, " values missing)")

mycolours<-colorRampPalette(brewer.pal(8, "Set2"))(nrow(lb))
p1 <- ggplot(data=lb, mapping=aes(x = SM, y=reads_raw, fill=LB)) + 
  geom_bar(stat="identity") +
  ylab("# reads") +
  xlab("") +
  theme_classic() +
  scale_fill_manual(values = mycolours) +
  ggtitle(title) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

if(nrow(lb)>10){
  p1 <- p1 + theme(legend.position = "none")
}

ggsave(plot_1_nb_reads, p1, width = 11, height = 7)

############################################################################

a<-sm[,c(1,7,9)]
colnames(a)[3] <- "endogenous"

## if one value is missing set the other one also to NaN 
a$endogenous[is.na(a$duplicates)] <-NaN
a$duplicates[is.na(a$endogenous)] <-NaN

nb = sum(is.na(sm$duplicates))
title <- "Mapped reads"
if(nb)title <- paste0(title, " (", nb, " values missing)")

sm2 <- melt(a, id.vars=1)
p2 <- ggplot(data=sm2, mapping=aes(x = SM, y=value, fill=variable)) + 
  geom_bar(stat="identity") +
  ylab("# reads") +
  xlab("") +
  theme_classic() +
  scale_fill_grey() +
  ggtitle(title) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

ggsave(plot_2_mapped, p2, width = 11, height = 7)

############################################################################

nb = sum(is.na(sm$endo_final_prop))
title <- "Endogenous content"
if(nb)title <- paste0(title, " (", nb, " values missing)")

p3 <- ggplot(data=sm, mapping=aes(x = SM, y=100*endo_final_prop)) + 
  geom_bar(stat="identity") +
  ylab("% reads") +
  xlab("") +
  theme_classic() +
  scale_fill_grey() +
  ggtitle(title) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

ggsave(plot_3_endogenous, p3, width = 11, height = 7)

############################################################################

nb = sum(is.na(sm$duplicates_prop))
title <- "Duplication level"
if(nb)title <- paste0(title, " (", nb, " values missing)")

p4 <- ggplot(data=sm, mapping=aes(x = SM, y=100*duplicates_prop)) + 
  geom_bar(stat="identity") +
  ylab("% reads") +
  xlab("") +
  theme_classic() +
  scale_fill_grey() +
  ggtitle(title) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

ggsave(plot_4_duplication, p4, width = 11, height = 7)


