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
out_2_mapped        = get_args(argsL, "out_2_mapped", "plot_2_mapped.png")
out_3_endogenous    = get_args(argsL, "out_3_endogenous", "plot_3_endogenous.png")
out_4_duplication   = get_args(argsL, "out_4_duplication", "plot_4_duplication.png")


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
my_plot <- make_barplot(
  data = sample_stats, x = "SM", y = "reads_raw", color_by = "SM", 
  fill_by = "SM", title = "Total number of raw reads", legend.position = "none"
  )

my_plot <- my_plot +
  scale_color_manual(values = colors_by_sample) + 
  scale_fill_manual(values = colors_by_sample)

ggsave(out_1_reads, my_plot, width = 11, height = 7)

#--------------------------------------------------------------------------#
# "Mapped reads"
# plot with unique and duplicated reads
mapped_reads <- data.frame(
  SM = c(sample_stats$SM, sample_stats$SM),
  number_reads = c(sample_stats$mapped_unique, sample_stats$duplicates),
  read_type = rep(c("Unique", "Duplicates"), each = n_samples)
)

my_plot <- make_barplot(
  data = mapped_reads, x = "SM", y = "number_reads", color_by = "read_type", 
  fill_by = "read_type", title = "Mapped unique and duplicated reads",
  legend.position = "top"
  )


ggsave(out_2_mapped, my_plot, width = 11, height = 7)
#--------------------------------------------------------------------------#
# "Endogenous content"
require(scales)
my_plot <- make_barplot(
  data = sample_stats, x = "SM", y = "endogenous_unique", color_by = "SM", 
  fill_by = "SM", title = "Endogenous content (unique reads)",
  legend.position = "none"
  )


my_plot <- my_plot +
  scale_color_manual(values = colors_by_sample) + 
  scale_fill_manual(values = colors_by_sample) +
  scale_y_continuous(labels = percent)


ggsave(out_3_endogenous, my_plot, width = 11, height = 7)

#--------------------------------------------------------------------------#
# "Duplication level"
require(scales)
my_plot <- make_barplot(
  data = sample_stats, x = "SM", y = "duplicates_prop", color_by = "SM", 
  fill_by = "SM", title = "Duplicates per sample (percentage)",
  legend.position = "none"
  )


my_plot <- my_plot +
  scale_color_manual(values = colors_by_sample) + 
  scale_fill_manual(values = colors_by_sample) +
  scale_y_continuous(labels = percent)


ggsave(out_4_duplication, my_plot, width = 11, height = 7)

############################################################################


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

