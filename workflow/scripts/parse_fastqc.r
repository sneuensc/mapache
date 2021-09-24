# function to parse fastqc output. It returns one of the 11 possible tables:
#  [1] ">>Adapter Content"                 ">>Basic Statistics"               
#  [3] ">>Overrepresented sequences"       ">>Per base N content"             
#  [5] ">>Per base sequence content"       ">>Per base sequence quality"      
#  [7] ">>Per sequence GC content"         ">>Per sequence quality scores"    
#  [9] ">>Per tile sequence quality"       ">>Sequence Duplication Levels"    
# [11] ">>Sequence Length Distribution"
parse_fastqc <- function(file_name, df_name){
    
    zip_dir <- sub(".zip" ,"", file_name)
    print(paste0("zip_dir: ", zip_dir))
    unzip(file_name, exdir = dirname(zip_dir))
    print(paste0("dirname(zip_dir): ", dirname(zip_dir)))
    print(paste0(zip_dir, "/fastqc_data.txt"))
    lines <- readLines(paste0(zip_dir, "/fastqc_data.txt"))
    # remove the very first line
    lines <- lines[-1]
    # remove the "fail" and "pass" tags from the lines that will become the names
    # of the data frames, and the apostrophe that brings some troubles
    lines <- gsub("\tfail", "", lines)
    lines <- gsub("\tpass", "", lines)
    lines <- gsub("\twarn", "", lines)
    lines <- gsub("\'", "prime", lines)
    # remove end of sections
    end_module <- grepl(">>END_MODULE", lines)
    lines <- lines[!end_module]

    # find headers and define groups
    is_header <- grepl(">>", lines)
    groups <- lines[is_header][cumsum(is_header)]

    # function to "read" the content from the lines
    Read <- function(x){
        if(length(x) > 1){
            t <- read.table(text = x, sep = "\t", fill = TRUE, comment = "x", header = T, skip = 1)
            return(t)
        }else{
            return(data.frame())
        }
    } 
        
    my_data_frames <- Map(Read, split(lines, groups))
    df <- my_data_frames[[paste0(">>", df_name)]]
    colnames(df) <- sub("^X.", "", colnames(df))
    return(df)
}



adjust_lengths <- function(length_dist){

    split_numbers <- function(sep, x){
        if(length(grep(sep, x))){
            as.numeric(unlist(strsplit(x, sep)))
        }else{
            return(as.numeric(c(x,x)))
        }
    }

    length_dist$length_1 <- sapply(length_dist$Length, function(i) split_numbers(sep="-", x=i))[1,]
    length_dist$length_2 <- sapply(length_dist$Length, function(i) split_numbers(sep="-", x=i))[2,]
    length_dist$Length <- apply(cbind(length_dist$length_1, length_dist$length_2), 1, function(x) mean(x))

    return(length_dist)

}
