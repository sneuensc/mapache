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
        --out_5_AvgReadDepth
        --split_plot=FALSE                           - TRUE or FALSE. If you are mapping many samples or to many genomes, break the long plot into subpanels.
        --n_col=1                                    - if split_plot=TRUE, the main plot will be broken into n_col number of columns
        --n_row=1                                    - if split_plot=TRUE, the main plot will be broken into n_row number of rows
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

out_1_reads         = get_args(argsL, "out_1_reads", "plot_1_nb_reads.png")
out_2_mapped        = get_args(argsL, "out_2_mapped", "plot_2_mapped.png")
out_3_endogenous    = get_args(argsL, "out_3_endogenous", "plot_3_endogenous.png")
out_4_duplication   = get_args(argsL, "out_4_duplication", "plot_4_duplication.png")
out_5_AvgReadDepth  = get_args(argsL, "out_5_AvgReadDepth", "plot_5_AvgReadDepth.png")
x_axis              = get_args(argsL, "x_axis", "sample")
split_plot          = eval( parse( text = get_args(argsL, "plit_plot", "FALSE") ) )

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
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = title) +
  theme(legend.position = legend.position) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))


  return(my_plot)
}

#--------------------------------------------------------------------------#
# make barplots for some of the statistics
samples_filename <- "../dev_mapache/Files/samples.txt"
sm <- "../dev_mapache/Files/SM_stats.csv"

samples <- read.table(samples_filename, header = T)

sample_stats <- read.csv(sm)

if(x_axis == "sample"){
  color_by = "genome"
  fill_by = "genome"
  x = "SM"
  n_x_bars <- length(unique(samples$SM))
}else if(x_axis == "genome"){
  color_by = "SM"
  fill_by = "SM"
  x = "genome"
  n_x_bars <- length(unique(samples$genome))
}

n_colors <- length(unique(sample_stats[,x]))

n_samples <- length(unique(samples$SM))

colors_by_sample <- colorRampPalette(brewer.pal(8, "Set2"))(n_colors)
#--------------------------------------------------------------------------#
# group stats in panels, if requested

break_into_panels <- function(my_plot, df, n_x_bars){
  n_col <- get_args(argsL, "n_col", round(sqrt(n_samples)) )
  n_row <- get_args(argsL, "n_row",ceiling(n_samples / n_col) )
  
  n_panels <- n_col * n_row
  
  df$group_to_plot <- cut_number(as.numeric(df[,x]), n_panels)
  new_plot <- my_plot + 
    facet_wrap(facets = df$group_to_plot, nrow = n_row, ncol = n_col, scales = "free_x") +
    theme(strip.text = element_blank())
  return(new_plot)
}

#--------------------------------------------------------------------------#
# "Total number of reads"
my_plot <- make_barplot(
  data = sample_stats, x = x, y = "reads_raw", color_by = color_by, 
  fill_by = fill_by, title = "Total number of raw reads", legend.position = "top"
  )

my_plot <- my_plot +
  scale_color_manual(values = colors_by_sample) + 
  scale_fill_manual(values = colors_by_sample)

if(split_plot){
  my_plot <- break_into_panels(my_plot = my_plot, df = sample_stats, n_x_bars = n_x_bars)
}

ggsave(out_1_reads, my_plot, width = 11, height = 7)

#--------------------------------------------------------------------------#
# "Mapped reads"
# plot with unique and duplicated reads
n_genomes <- length(unique(sample_stats$genome))

mapped_reads <- data.frame(
    SM = rep(sample_stats$SM, n_samples),
    number_reads = c(sample_stats$mapped_unique, sample_stats$mapped_unique+sample_stats$duplicates),
    read_type = rep(c("Unique", "Duplicates"), each = n_samples * n_genomes),
    genome = rep(sample_stats$genome, n_samples)
)


my_plot <- ggplot(
    mapped_reads,
    aes_string(x = x, y = "number_reads",  fill = fill_by, color = color_by,
    alpha = "read_type",
               group = fill_by) 
) +
    labs(title = "Mapped unique and duplicated reads") +
    theme(legend.position = "top") +
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
    scale_alpha_manual(values = c(0.5, 1)) +
    geom_bar(stat = "identity", position = position_dodge()) +
  scale_color_manual(values = colors_by_sample) + 
  scale_fill_manual(values = colors_by_sample) 

if(split_plot){
  my_plot <- break_into_panels(my_plot = my_plot, df = sampled_reads, n_x_bars = n_x_bars)
}

ggsave(out_2_mapped, my_plot, width = 11, height = 7)
#--------------------------------------------------------------------------#
# "Endogenous content"
require(scales)
my_plot <- make_barplot(
  data = sample_stats, x = "SM", y = "endogenous_unique", color_by = color_by, 
  fill_by = fill_by, title = "Endogenous content (unique reads)",
  legend.position = "top"
  )


my_plot <- my_plot +
  scale_color_manual(values = colors_by_sample) + 
  scale_fill_manual(values = colors_by_sample) +
  scale_y_continuous(labels = percent)

if(split_plot){
  my_plot <- break_into_panels(my_plot = my_plot, df = sample_stats, n_x_bars = n_x_bars)
  }

ggsave(out_3_endogenous, my_plot, width = 11, height = 7)

#--------------------------------------------------------------------------#
# "Duplication level"
require(scales)
my_plot <- make_barplot(
  data = sample_stats, x = "SM", y = "duplicates_prop", color_by = color_by, 
  fill_by = fill_by, title = "Duplicates per sample (percentage)",
  legend.position = "top"
  )


my_plot <- my_plot +
  scale_color_manual(values = colors_by_sample) + 
  scale_fill_manual(values = colors_by_sample) +
  scale_y_continuous(labels = percent)

if(split_plot){
  my_plot <- break_into_panels(my_plot = my_plot, df = sample_stats, n_x_bars = n_x_bars)
}

ggsave(out_4_duplication, my_plot, width = 11, height = 7)

#--------------------------------------------------------------------------#
# "Average read depth"
require(scales)
my_plot <- make_barplot(
  data = sample_stats, x = "SM", y = "read_depth", color_by = color_by, 
  fill_by = fill_by, title = "Average read depth",
  legend.position = "top"
  )


my_plot <- my_plot +
  scale_color_manual(values = colors_by_sample) + 
  scale_fill_manual(values = colors_by_sample) 

if(split_plot){
  my_plot <- break_into_panels(my_plot = my_plot, df = sample_stats, n_x_bars = n_x_bars)
}

ggsave(out_5_AvgReadDepth, my_plot, width = 11, height = 7)
############################################################################


# ############################################################################

