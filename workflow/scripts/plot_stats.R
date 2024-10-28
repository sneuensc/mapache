library(ggplot2)
library(RColorBrewer)
library(scales)
library(plyr)

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
        --out_1_reads=plot_1_nb_reads.png            - plot file name of the 'Total number of reads'.
        --out_2_mapped=plot_2_mapped.png             - plot file name of the 'Mapped reads'.
        --out_3_endogenous=plot_3_endogenous.png     - plot file name of the 'Endogenous content'.
        --out_4_duplication=plot_4_duplication.png   - plot file name of the 'Duplication level'.
        --out_5_AvgReadDepth                         - plot file name of the 'Average read depth'.
        --out_6_Sex                                  - plot file name of the 'Sex determination' (or 'None' if not available).

        Parameters:
        --n_col=1                                    - number of figures per row if multiple samples and genomes are used
        --width=11                                   - width of the plot
        --height=7                                   - height of the plot
        --color=blue                                 - color of the plot
        --show_numbers                               - print numbers to columns if the number of samples is below the given number (10)
        --thresholds=list('XX'=c(0.8, 1), 'XY'=c(0, 0.6), \
                          'consistent with XX but not XY'=c(0.6, 1), \
                          'consistent with XY but not XX'=c(0, 0.8))  - sex determination thresholds
        --sex_ribbons=c('XX'='red', 'XY'='blue')     - what sex determination ribbons to color in the plot
        --help                                       - print this text
 
      Example:
      Rscsript plot_stats.R [] \n\n")
  
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
sm                  = get_args(argsL, "sm", "SM_stats.csv")
out_1_reads         = get_args(argsL, "out_1_reads", "plot_1_nb_reads.png")
out_2_mapped        = get_args(argsL, "out_2_mapped", "plot_2_mapped.png")
out_3_endogenous    = get_args(argsL, "out_3_endogenous", "plot_3_endogenous.png")
out_4_duplication   = get_args(argsL, "out_4_duplication", "plot_4_duplication.png")
out_5_AvgReadDepth  = get_args(argsL, "out_5_AvgReadDepth", "plot_5_AvgReadDepth.png")
out_6_Sex           = get_args(argsL, "out_6_Sex", "plot_6_Sex.png")
x_axis              = get_args(argsL, "x_axis", "auto")
n_col               = as.numeric( get_args(argsL, "n_col", 1 ) )
width               = as.numeric( get_args(argsL, "width", 11 ) )
height              = as.numeric( get_args(argsL, "height", 7 ) )
color               = get_args(argsL, "color", "blue")
show_numbers        = as.numeric( get_args(argsL, "show_numbers", 10 ) )


sex_thresholds      = get_args(
  argsL, 
  "thresholds", 
  'list("hg19"=list("XX"=c(0.8, 1), "XY"=c(0, 0.6), "consistent with XX but not XY"=c(0.6, 1),"consistent with XY but not XX"=c(0, 0.8) ), "GRCh38"=list("XX"=c(0.8, 1), "XY"=c(0, 0.6), "consistent with XX but not XY"=c(0.6, 1),"consistent with XY but not XX"=c(0, 0.8) ))'
)

sex_ribbons         = get_args(argsL, "sex_ribbons", 'c("XX"="red", "XY"="blue")')

############################################################################
#--------------------------------------------------------------------------#
# read the file
sample_stats <- read.csv(sm)

is_sample_file <- TRUE

if('LB' %in% colnames(sample_stats)){
  is_sample_file <- FALSE
  sample_stats$SM <- paste0(sample_stats$SM, '__', sample_stats$LB)
  sample_stats <- sample_stats[, names(sample_stats)!='LB']
}

if(is_sample_file){
  sample_stats$Sex <- factor(
    sample_stats$Sex,
    levels = unique(sample_stats$Sex),
    ordered = T
  )
}

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


n_genomes <- length(unique(sample_stats$genome))
n_samples <- length(unique(sample_stats$SM))

## define the x axis
if(x_axis == "auto"){
  if(n_genomes > n_samples){
    x_axis = "genome"
  } else {
    x_axis = "sample"
  }
}

if(x_axis == "sample"){
  x = "SM"
  group = "genome"
  n_figures = n_genomes
}else if(x_axis == "genome"){
  x = "genome"
  group = "SM"
  n_figures = n_samples
}

#--------------------------------------------------------------------------#
# generic barplot with ggplot2

make_barplot <- function(data, x, y, group, n_figures, cols, color="blue", title = "", percent=F){
  my_plot <- ggplot(data, aes_string(x = x, y = y)) +
    geom_bar(stat = "identity", position = position_dodge(), fill=color) +
    labs(title = title) +
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5),
          legend.position = "none") +
    xlab('')
  
  
  ## add numbers if there are not too many columns
  if(nrow(data)<=show_numbers){
    if(percent){
      meanVal = 100*mean(data[,y], na.rm=T)
      if(1e8 > meanVal & meanVal > 1e-3){
        vec = paste0(format(100*data[,y], scientific = FALSE, big.mark=",", digits = 3),'%')
      } else {
        vec = paste0(format(100*data[,y], scientific = TRUE, big.mark=",", digits = 3),'%')
      }
      my_plot <- my_plot + geom_text(aes(label=vec), vjust=-0.25, color="black") +
        scale_y_continuous(labels = scales::percent)
    }else{
      meanVal = mean(data[,y], na.rm=T)
      if(1e8 > meanVal & meanVal > 1e-3){
        vec = format(data[,y], scientific = F, big.mark=",", digits = 3)
      } else {
        vec = format(data[,y], scientific = T, big.mark=",", digits = 3)
      }
      my_plot <- my_plot + geom_text(aes(label=vec), vjust=-0.25, color="black")
    }
  }
  
  ## split figure if multidimensional
  if(n_figures > 1){
    my_plot <- my_plot + facet_wrap(as.formula(paste("~", group)), ncol = cols)
  }
  
  return(my_plot)
}


#--------------------------------------------------------------------------#
# "Total number of reads"
my_plot <- make_barplot(data = sample_stats, x = x, y = "reads_raw",
                        group = group, n_figures = n_figures, cols=n_col,
                        color=color, title = "Total number of raw reads")

ggsave(out_1_reads, my_plot, width = width, height = height)

#--------------------------------------------------------------------------#

# "Mapped reads"
# plot with unique and duplicated reads
mapped_reads <- data.frame(
  SM = rep(sample_stats$SM, 2),
  number_reads = c(sample_stats$mapped_unique, sample_stats$mapped_unique+sample_stats$duplicates),
  read_type = rep(c("Unique", "Duplicates"), each = n_samples * n_genomes),
  genome = rep(sample_stats$genome, 2)
)

my_plot <- ggplot(mapped_reads, aes_string(x = x, y = "number_reads",
                                           alpha = "read_type", group = group)) +
  labs(title = "Mapped unique and duplicated reads") +
  theme(legend.position = "top") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  geom_bar(stat = "identity", position = position_dodge(), fill=color) +
  guides(fill="none", color="none")

## add numbers if there are not too many columns
if(nrow(mapped_reads)/2<=show_numbers){
  my_plot <- my_plot + geom_text(aes(label=format(mapped_reads[,"number_reads"], big.mark=",")),
                                 vjust=-0.25, color="black", show.legend = FALSE)
}

## split figure if multidimensional
if(n_figures > 1){
  my_plot <- my_plot + facet_wrap(as.formula(paste("~", group)), ncol = n_col)
}

ggsave(out_2_mapped, my_plot, width = width, height = height)

#--------------------------------------------------------------------------#
# "Endogenous content"
my_plot <- make_barplot(data = sample_stats, x = x, y = "endogenous_unique",
                        group = group, n_figures = n_figures, cols=n_col,
                        color=color, title = "Endogenous content (unique reads)", percent=T) 

ggsave(out_3_endogenous, my_plot, width = width, height = height)

#--------------------------------------------------------------------------#
# "Duplication level"
my_plot <- make_barplot(data = sample_stats, x = x, y = "duplicates_prop",
                        group = group, n_figures = n_figures, cols=n_col,
                        color=color, title = "Duplicates per sample (percentage)", percent=T)

ggsave(out_4_duplication, my_plot, width = width, height = height)

#--------------------------------------------------------------------------#
# "Average read depth"
my_plot <- make_barplot(data = sample_stats, x = x, y = "read_depth",
                          group = group, n_figures = n_figures, cols=n_col,
                          color=color, title = "Average read depth")
  
ggsave(out_5_AvgReadDepth, my_plot, width = width, height = height)


#--------------------------------------------------------------------------#
# "Sex"
if('Sex' %in% colnames(sample_stats)){ ## not present (library file)
  if(sum(!is.na(sample_stats$Sex))){  ## not requested
    if(sum(!is.na(sample_stats$Sex_Rx))){  ## not possible to compute
      #sex_thresholds='list("hg19"?list( "XX"?c(0.8, 1), "XY"?c(0, 0.6), "consistent with XX but not XY"?c(0.6, 1), "consistent with XY but not XX"?c(0, 0.8) ),     "GRCh38"?list( "XX"?c(0.8, 1), "XY"?c(0, 0.6), "consistent with XX but not XY"?c(0.6, 1), "consistent with XY but not XX"?c(0, 0.8) ))'
      sex_thresholds      = eval(parse(text=gsub("\\?", "=", sex_thresholds)))
      
      ## get the ribbons as a matrix
      #sex_ribbons='c("XX"?"red", "XY"?"blue")'
      sex_ribbons         = eval(parse(text=gsub("\\?", "=", sex_ribbons)))
      m <-NULL
      for(i in names(sex_thresholds)){
        tmp <- do.call(rbind.data.frame, sex_thresholds[i])
        tmp <- as.data.frame(t(tmp[,colnames(tmp) %in% names(sex_ribbons)]))
        colnames(tmp) <- c("min","max")
        tmp$genome <- i
        tmp$Sex <- row.names(tmp)
        m <- rbind(m, tmp)
      }
      
      ## set all assignments not defined in 'sex_ribbons' to 'other'
      sample_stats$Sex2 <- as.character(sample_stats$Sex)
      sample_stats$Sex2[!sample_stats$Sex2 %in% names(sex_ribbons)] = "other"
      sex_ribbons2 <- c(sex_ribbons, "other"="darkgray")
      
      my_plot <- ggplot() +
        theme_bw() +
        geom_rect(data=m, aes(xmin = -Inf, xmax = Inf, ymin = min, ymax = max, fill = Sex), alpha = 0.2) +
        geom_errorbar(data=sample_stats, aes(x = SM, y = Sex_Rx, ymin = Sex_Rx-Sex_CI, ymax = Sex_Rx+Sex_CI, color = Sex2), width=0.25) +
        geom_point(data=sample_stats, aes(x=SM, y=Sex_Rx, color = Sex2)) +
        scale_color_manual(values = sex_ribbons2, name="Sex inferred") +
        scale_fill_manual(values = sex_ribbons, name="Sex range") +
        scale_y_continuous(breaks=seq(0,1,0.25)) +
        theme(legend.position = "top") +
        theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
      
      ## split figure if multidimensional
      if(n_figures > 1){
        my_plot <- my_plot + facet_wrap(as.formula(paste("~", group)), ncol = n_col)
      }
      
    }else{
      my_plot <- ggplot() +
        geom_text(data = data.frame(x=1,y=1),
                  aes(x = x, y = y, label = "Sex inference not possible")) +
        theme_void()
    }
  }else{
    my_plot <- ggplot() +
      geom_text(data = data.frame(x=1,y=1),
                aes(x = x, y = y, label = "Sex inference was not requested")) +
      theme_void()
  }
  
  ggsave(out_6_Sex, my_plot, width = width, height = height)
}
# ############################################################################