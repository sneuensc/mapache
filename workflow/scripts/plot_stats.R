library(ggplot2)
library(RColorBrewer)
library(scales)
require(plyr)

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
# print(argsL)

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
# fq                 = get_args(argsL, "FQ", "FASTQ.csv")
# lb                 = get_args(argsL, "LB", "LB.csv")
sm                  = get_args(argsL, "SM", "SM_stats.csv")
sex_path            = get_args(argsL, "sex_path", "sex.out")
sex_thresholds      = get_args(
                        argsL, 
                        "thresholds", 
                        'list("hg19"=list("XX"=c(0.8, 1), "XY"=c(0, 0.6), "consistent with XX but not XY"=c(0.6, 1),"consistent with XY but not XX"=c(0, 0.8) ), "GRCh38"=list("XX"=c(0.8, 1), "XY"=c(0, 0.6), "consistent with XX but not XY"=c(0.6, 1),"consistent with XY but not XX"=c(0, 0.8) ))'
                        )
#print(paste0("SEX THRESHOLDS: \n", sex_thresholds))

sex_thresholds      = eval(parse(text=gsub("\\?", "=", sex_thresholds)))


sex_ribbons         = eval(parse(text=get_args(argsL, "sex_ribbons", 'c("XX", "XY")')))

out_1_reads         = get_args(argsL, "out_1_reads", "plot_1_nb_reads.png")
out_2_mapped        = get_args(argsL, "out_2_mapped", "plot_2_mapped.png")
out_3_endogenous    = get_args(argsL, "out_3_endogenous", "plot_3_endogenous.png")
out_4_duplication   = get_args(argsL, "out_4_duplication", "plot_4_duplication.png")
out_5_AvgReadDepth  = get_args(argsL, "out_5_AvgReadDepth", "plot_5_AvgReadDepth.png")
out_6_Sex           = get_args(argsL, "out_6_Sex", "plot_6_Sex.png")
x_axis              = get_args(argsL, "x_axis", "sample")
split_plot          = eval( parse( text = get_args(argsL, "split_plot", "FALSE") ) )

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
sample_stats <- read.csv(sm)
sample_stats$SM <- factor(
  sample_stats$SM,
  levels = unique(sample_stats$SM),
  ordered = T
  )
sample_stats$genome <- factor(
  sample_stats$genome,
  levels = unique(sample_stats$genome),
  ordered = T
  )


if(x_axis == "sample"){
  color_by = "genome"
  fill_by = "genome"
  x = "SM"
  n_x_bars <- length(unique(sample_stats$SM))
}else if(x_axis == "genome"){
  color_by = "SM"
  fill_by = "SM"
  x = "genome"
  n_x_bars <- length(unique(sample_stats$genome))
}

n_colors <- length(unique(sample_stats[,color_by]))

n_samples <- length(unique(sample_stats$SM))

colors_by_sample <- colorRampPalette(brewer.pal(8, "Set2"))(n_colors)
#--------------------------------------------------------------------------#
# group stats in panels, if requested
break_into_panels <- function(my_plot, df, n_x_bars){
  if(n_x_bars == 1){
    return(my_plot)
  }else{
    n_col <- as.numeric( get_args(argsL, "n_col", 30 ) )
    n_row <- as.numeric( get_args(argsL, "n_row", ceiling(n_x_bars / n_col) ) )

    # we need as many panels as number of groups in the x-axis
    if( n_col * n_row < n_x_bars ){
      if(n_col > n_row){
        n_row <- min(ceiling(n_x_bars / n_col), 1) # at least one row
      }else{
        n_col <- min(ceiling(n_x_bars / n_row), 1) # at least one column
      }
    }

    n_panels <- min(c(n_col * n_row, n_x_bars))

    df$group_to_plot <- cut_number(as.numeric(df[,x]), n_panels)

    new_plot <- my_plot + 
      facet_wrap(facets = df$group_to_plot, nrow = n_row, ncol = n_col, scales = "free_x") +
      theme(strip.text = element_blank())

    return(new_plot)
  }
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
    SM = rep(sample_stats$SM, 2),
    number_reads = c(sample_stats$mapped_unique, sample_stats$mapped_unique+sample_stats$duplicates),
    read_type = rep(c("Unique", "Duplicates"), each = n_samples * n_genomes),
    genome = rep(sample_stats$genome, 2)
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
  my_plot <- break_into_panels(my_plot = my_plot, df = mapped_reads, n_x_bars = n_x_bars)
}

ggsave(out_2_mapped, my_plot, width = 11, height = 7)

#--------------------------------------------------------------------------#
# "Endogenous content"
my_plot <- make_barplot(
  data = sample_stats, x = x, y = "endogenous_unique", color_by = color_by, 
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
my_plot <- make_barplot(
  data = sample_stats, x = x, y = "duplicates_prop", color_by = color_by, 
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
my_plot <- make_barplot(
  data = sample_stats, x = x, y = "read_depth", color_by = color_by, 
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

#--------------------------------------------------------------------------#
# "Inferred molecular sex"
# Note that
# For each sample, sex inference was
# requested (file.sex)
#   - inferred
#   - not inferred if there were not reads mapping to sex chromosomes (Rx is NaN)
# not requested (file.nosex with NaN)
#sex_path <- "SM.GENOME"

# Build data frame from individual tables
sex <- data.frame()

for(i in 1:nrow(sample_stats)){
  
  sample_name <- sample_stats$SM[i]
  genome_name <- sample_stats$genome[i]
  path_ind <- sub( "GENOME", 
                   genome_name, 
                   sub( "SM", sample_name, sex_path)
                   )
  path_ind_sex <- paste0(path_ind, ".sex")
  
  if(file.exists(path_ind_sex)){
    if(nrow(sex)){
      sex_ind <- read.csv(path_ind_sex, stringsAsFactors = F)
      sex_ind$SM <- sample_name
      sex_ind$genome <- genome_name
      sex <- rbind(sex, sex_ind)
    }else{
      sex <- read.csv(path_ind_sex, stringsAsFactors = F)
      sex$SM <- sample_name
      sex$genome <- genome_name
    }
  }
}

# print(sex)
if(nrow(sex)){
  # This is necessary to plot all individuals, regardles of whether
  # sex inference was run or if the sex could not be assigned,
  # also considering all genomes present in the dataset
  sex <- join(sample_stats[,c("SM", "genome")], sex)
  sex$genome <- as.character(sex$genome)
  sex_with_data <- sex[!is.na(sex$Rx),]
  sex_no_data <- sex[is.na(sex$Rx),]

  # Get thresholds that will be plotted as ribbons/rectangles
  threshold_coord <- data.frame()
  for(genome in names(sex_thresholds)){
    
    for(sex_name in names(sex_thresholds[[genome]])){
      coord <- data.frame(
        Sex = sex_name,
        ymin = sex_thresholds[[genome]][[sex_name]][1],
        ymax = sex_thresholds[[genome]][[sex_name]][2],
        genome = genome
      )
      
      coord_not_assigned <- coord
      coord_not_assigned$genome <- paste0(coord_not_assigned$genome, "_not_assigned")
      if(nrow(threshold_coord)){
        threshold_coord <- rbind(threshold_coord, coord, coord_not_assigned)
      }else{
        threshold_coord <- rbind(coord, coord_not_assigned)
      }
    }
  }

  # add a point for the individuals/genomes
  # for which sex was not assigned/requested
  # (only when assignment was requested at least once; if it was not requested
  # at all, there will be a blank plot)
  not_assigned <- is.na(sex$Rx)
  #print(length(not_assigned))
  if(length(not_assigned)){
    sex$genome[not_assigned] <- paste0(sex$genome[not_assigned], "_not_assigned")
    sex$Rx[not_assigned] <- -.1
    sex$CI1[not_assigned] <- -.1
    sex$CI2[not_assigned] <- -.1
    sex$Sex[is.na(sex$Sex)] <- "Not requested"

  }

  # remove sex names that were not requested to be plotted in ribbons
  threshold_coord <- threshold_coord[threshold_coord$Sex %in% sex_ribbons,]

  all_sexes <- unique(
    c(
      sapply(names(sex_thresholds), function(genome) names(sex_thresholds[[genome]])),
      as.character(sex$Sex)
    )
  )


  
  ribbons_color <- colorRampPalette(brewer.pal(8, "Set3"))(length(all_sexes))
  sex_color <- colorRampPalette(brewer.pal(8, "Set2"))(length(all_sexes))
  my_plot <-     ggplot() +
    theme_bw() +
    geom_point(data=sex, 
               aes(
                 x = SM, y = Rx, 
                 color = Sex
               )) +
    geom_rect(data=threshold_coord,
              mapping=aes(xmin = 0, xmax=length(unique(sex$SM))+1, 
                          ymin=ymin,ymax=ymax, fill=Sex)
    ) +
    geom_point(data=sex, 
               aes(
                 x = SM, y = Rx, 
                 color = Sex
               )) +
    geom_errorbar(data=sex, aes( x = SM, y = Rx, ymin = CI1, ymax = CI2, color = Sex), width=0.25) + 
    scale_color_manual(breaks=all_sexes, values = sex_color, name="Sex inferred") +
    scale_fill_manual(breaks=all_sexes, 
                      values = alpha(ribbons_color, 0.5),
                      name="Thresholds") +
    facet_grid(genome~.) +
    theme(legend.position = "top")
}else{
  my_plot <- ggplot() +
    geom_text(data = data.frame(x=1,y=1), 
              aes(x = x, y = y, label = "Sex inference was not requested")) +
    theme_void()
}
  

ggsave(out_6_Sex, my_plot, width = 11, height = 7)
# ############################################################################

