args = commandArgs(trailingOnly = TRUE)
bwtool.file = args[1]

data <- read.delim(bwtool.file, as.is=T, sep = "\t", header = TRUE)
median <- median(sqrt(data$sum_of_squares/100000))
writeLines(paste(median,sep=" "))