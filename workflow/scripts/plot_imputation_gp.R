library(ggplot2)
library(data.table)

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
        --input=GP.txt                - file with list of gp values (three columns)
        --output=plot.png             - output file
        --sample='Mota'               - sample name
        --gp=0.6,0.8                  - GP filter threshold, separated by comma
        --width=11                    - width of the plot
        --height=7                    - height of the plot
        --help                        - print this text

      Example:
      Rscsript plot_imputation.R [] \n\n")

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
input               = get_args(argsL, "input", "")
output              = get_args(argsL, "output", "")
sample              = get_args(argsL, "sample", "")
gp                  = as.numeric(unlist(strsplit(get_args(argsL, "gp", "0.8"),',')))
height              = as.numeric( get_args(argsL, "height", 12 ) )
width               = as.numeric( get_args(argsL, "width", 18 ) )

## read data
d <- fread(input, sep='\t', h=F)
d[, c('X1', 'X2', 'X3') :=tstrsplit(V2, ',')]

## get non-reference sites
d1 <- d[d$V1 != '0/0',][,-(1:2)]

## sort gp's
d1 <- data.frame(t(apply(d1,1,function(x)return(sort(x, decreasing=T)))))
d1 <- data.frame(lapply(d1,as.numeric))

## filter for gp
count <- unlist(lapply(gp, function(x)return(sum(d1$X1>x))))
rr <- data.frame(gp, 'str'=paste0(count, " (", round(100*count/nrow(d1)), "%) of the imputed non-reference sites are kept"))

colors <- c("X1" = "black", "X2" = gray(0.6), "X3" = gray(0.9))

my_plot <- ggplot(d1) +
  geom_histogram(aes(x=X3, fill="X3"), bins=100, color='black', alpha=0.5) +
  geom_histogram(aes(x=X2, fill="X2"), bins=100, color='black', alpha=0.5) +
  geom_histogram(aes(x=X1, fill="X1"), bins=100) +
  geom_text(data=rr, aes(x=gp+0.01, y=Inf, label=gp, color=factor(gp)), vjust=4, hjust=0) +
  labs(x = "genotype probabilities (GP)", y= "# sites", fill = "Genotype", color="GP filter") +
  geom_vline(xintercept=rr$gp, color = "blue") +
  ggtitle(paste0("Imputation genotype probabilities of sample '", sample, "' [", nrow(d1), " (", round(100*nrow(d1)/nrow(d), 1), "%) imputed non-reference sites]")) +
  scale_fill_manual(values = colors, labels=c('most probable','intermediate probable','least probable')) +
  scale_color_manual(values = rep("blue", nrow(rr)), labels=rr$str,
                     guide = guide_legend(override.aes = list(label = rr$gp))) +
  theme(legend.position = c(0.3, 0.7),
        legend.background = element_rect(fill='transparent'))

ggsave(output, my_plot, width = width, height = height)

