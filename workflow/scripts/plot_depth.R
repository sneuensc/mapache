library(ggplot2)
library(RColorBrewer) 

## input
data<-read.csv(snakemake@input[['sample_depth']])

## output
plot_5_AvgReadDepth = snakemake@output[['plot_5_AvgReadDepth']]
plot_6_AvgReadDepth_MT = snakemake@output[['plot_6_AvgReadDepth_MT']]
plot_7_Sex = snakemake@output[['plot_7_Sex']]

############################################################################
svg(plot_5_AvgReadDepth, width = 11, height = 7)

nb = sum(is.na(data$AvgReadDepth_tot))
title <- "Average read depth"
if(nb)title <- paste0(title, " (", nb, " values missing)")

ggplot(data=data, mapping=aes(x = Sample, y=AvgReadDepth_tot)) +
  geom_bar(stat="identity") +
  ylab("Average read depth") +
  xlab("") +
  theme_classic() +
  ggtitle(title) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

dev.off()

############################################################################
svg(plot_6_AvgReadDepth_MT, width = 11, height = 7)

nb = sum(is.na(data$AvgReadDepth_MT))
title <- "Average read depth on MT"
if(nb)title <- paste0(title, " (", nb, " values missing)")

ggplot(data=data, mapping=aes(x = Sample, y=AvgReadDepth_MT)) +
  geom_bar(stat="identity") +
  ylab("Average read depth") +
  xlab("") +
  theme_classic() +
  ggtitle(title) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

dev.off()

############################################################################
svg(plot_7_Sex, width = 11, height = 7)

data$X025<-as.numeric(sapply(strsplit(as.character(data$X95.CI),"-"), `[`, 1))
data$X975<-as.numeric(sapply(strsplit(as.character(data$X95.CI),"-"), `[`, 2))

malelimit<-0.075                                                         
femalelimit<-0.016

nb = sum(is.na(data$R_y))
title <- "Infered sex"
if(nb)title <- paste0(title, " (", nb, " values missing)")

if(nrow(data) - nb > 0){ ## plot only if the values are not NaN
  ggplot(data=data, mapping=aes(x = Sample, y=R_y)) +
    ylim(0, 0.095) +
    geom_point() +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = malelimit, ymax = Inf),
              fill = "dodgerblue", alpha = 0.03) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = femalelimit),
              fill = "orangered", alpha = 0.03) +
    geom_text(aes(x = nrow(data)/2+0.5, y = femalelimit-0.01, label = "Female", fontface=2),
              size = 10, vjust = 0, hjust = 0.5, color = "white") +
    geom_text(aes(x = nrow(data)/2+0.5, y = malelimit+0.01, label = "Male", fontface=2),
              size = 10, vjust = 0, hjust = 0.5, color = "white") +
    geom_point() +
    geom_errorbar(aes(ymin = X025, ymax = X975, width = 0.15)) +
    ylab("NchrY/(NchrY+NchrX)") +
    xlab("") +
    theme_classic() +
    ggtitle(title) +
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
} else {
  ## empty plot
  ggplot(data=data, mapping=aes(x = Sample, y=R_y)) +
    ylim(0, 0.095) +
    theme_classic() +
    ggtitle("Infered sex") +
    geom_hline(yintercept = malelimit, linetype="dashed", color = "blue") +
    geom_hline(yintercept = femalelimit, linetype="dashed", color = "red") 
}

dev.off()

