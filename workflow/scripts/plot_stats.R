library(ggplot2)
library(RColorBrewer)
#library(reshape)


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
        --samples=samples.txt                        - input containing the samples names, libraries, etc.
        --SM=SM.csv                                  - summary statistics file name at the sample level. Default: SM.csv
        --FQ=FASTQ.csv                               - summary statistics file name at the fastq level. Default: FASTQ.csv
        --LB=LB.csv                                  - summary statistics file name at the library level. Default: LB.csv
        
        Output files:  
        --out_1_reads=plot_1_nb_reads.png            - plot file name of the 'Total number of reads'. Default: plot_1_nb_reads.png
        --out_2_mapped=plot_2_mapped.png             - plot file name of the 'Mapped reads'. Default: plot_2_mapped.png
        --out_3_endogenous=plot_3_endogenous.png     - plot file name of the 'Endogenous content'. Default: plot_3_endogenous.png
        --out_4_duplication=plot_4_duplication.png   - plot file name of the 'Duplication level'. Default: plot_4_duplication.png
          
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

# getting arguments
# fq                   = get_args(argsL, "FQ", "FASTQ.csv")
# lb                   = get_args(argsL, "LB", "LB.csv")
sm                   = get_args(argsL, "SM", "SM.csv")
samples_filename                   = get_args(argsL, "samples", "samples.txt")

out_1_reads      = get_args(argsL, "out_1_reads", "plot_1_nb_reads.png")
# out_2_mapped        = get_args(argsL, "out_2_mapped", "plot_2_mapped.png")
# out_3_endogenous    = get_args(argsL, "out_3_endogenous", "plot_3_endogenous.png")
# out_4_duplication   = get_args(argsL, "out_4_duplication", "plot_4_duplication.png")


############################################################################
#--------------------------------------------------------------------------#
# generic barplot with ggplot2

make_barplot <- function(
  data, x, y, color_by, fill_by, title = "", legend.position = "top"
  ){

  my_plot <- ggplot(
    data,
    aes_string(x = x, y = y, color = color_by, fill = fill_by) 
  ) +
  geom_bar(stat = "identity") +
  labs(title = title) +
  theme(legend.position = legend.position) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))


  return(my_plot)
}

#--------------------------------------------------------------------------#
# make barplots for some of the statistics

samples <- read.table(samples_filename, header = T)
n_samples <- length(unique(samples$SM))
colors_by_sample <- colorRampPalette(brewer.pal(8, "Set2"))(n_samples)

sample_stats <- read.csv(sm)

#--------------------------------------------------------------------------#
# "Total number of reads"
plot_nb_reads <- make_barplot(
  data = sample_stats, x = "SM", y = "reads_raw", color_by = "SM", 
  fill_by = "SM", title = "Total number of raw reads", legend.position = "top"
  )

plot_nb_reads <- plot_nb_reads +
  scale_color_manual(values = colors_by_sample) + 
  scale_fill_manual(values = colors_by_sample)

ggsave(out_1_reads, plot_nb_reads, width = 11, height = 7)

#--------------------------------------------------------------------------#

# title <- "Total number of reads"
# if(nb) title <- paste0(title, " (", nb, " values missing)")

# mycolours<-colorRampPalette(brewer.pal(8, "Set2"))(nrow(lb))

# p1 <- ggplot(data=lb, mapping=aes(x = SM, y=reads_raw, fill=LB)) + 
#   geom_bar(stat="identity") +
#   ylab("# reads") +
#   xlab("") +
#   theme_classic() +
#   scale_fill_manual(values = mycolours) +
#   ggtitle(title) +
#   theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

# if(nrow(lb)>10){
#   p1 <- p1 + theme(legend.position = "none")
# }
# if(length(unique(lb$genome))>1){
#   p1 <- p1 + facet_grid(cols = vars(genome))
# }

# ggsave(plot_1_nb_reads, p1, width = 11, height = 7)

############################################################################

# a<-sm[,c('genome','SM','duplicates', 'mapped_unique')]
# #a<-sm[,c(1,2,8,10)]
# colnames(a)[4] <- "endogenous"

# ## if one value is missing set the other one also to NaN 
# a$endogenous[is.na(a$duplicates)] <-NaN
# a$duplicates[is.na(a$endogenous)] <-NaN

# nb = sum(is.na(sm$duplicates))
# title <- "Mapped reads"
# if(nb)title <- paste0(title, " (", nb, " values missing)")

# sm2 <- melt(a, id.vars=c('genome','SM'))
# p2 <- ggplot(data=sm2, mapping=aes(x = SM, y=value, fill=variable)) + 
#   geom_bar(stat="identity") +
#   ylab("# reads") +
#   xlab("") +
#   theme_classic() +
#   scale_fill_grey() +
#   ggtitle(title) +
#   theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

# if(length(unique(sm2$genome))>1){
#   p2 <- p2 + facet_grid(cols = vars(genome))
# }

# ggsave(plot_2_mapped, p2, width = 11, height = 7)

# ############################################################################

# nb = sum(is.na(sm$endogenous_unique))
# title <- "Endogenous content"
# if(nb)title <- paste0(title, " (", nb, " values missing)")

# p3 <- ggplot(data=sm, mapping=aes(x = SM, y=100*endogenous_unique)) + 
#   geom_bar(stat="identity") +
#   ylab("% reads") +
#   xlab("") +
#   theme_classic() +
#   scale_fill_grey() +
#   ggtitle(title) +
#   theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

# if(length(unique(sm$genome))>1){
#   p3 <- p3 + facet_grid(cols = vars(genome))
# }

# ggsave(plot_3_endogenous, p3, width = 11, height = 7)

# ############################################################################

# nb = sum(is.na(sm$duplicates_prop))
# title <- "Duplication level"
# if(nb)title <- paste0(title, " (", nb, " values missing)")

# p4 <- ggplot(data=sm, mapping=aes(x = SM, y=100*duplicates_prop)) + 
#   geom_bar(stat="identity") +
#   ylab("% reads") +
#   xlab("") +
#   theme_classic() +
#   scale_fill_grey() +
#   ggtitle(title) +
#   theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

# if(length(unique(sm$genome))>1){
#   p4 <- p4 + facet_grid(cols = vars(genome))
# }

# ggsave(plot_4_duplication, p4, width = 11, height = 7)

