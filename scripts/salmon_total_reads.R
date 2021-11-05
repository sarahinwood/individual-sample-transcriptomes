library(readr)
library(dplyr)
library(data.table)

##salmon log prints total # read pairs, counts are given in read pairs also (bbduk is in single reasd i.e. 2x number salmon and counts report)
##grep in snakemake for rows with "Observed " which should only be row with total reads
##read in txt with that

##find all files
salmon_total_reads <- list.files(path = "output/salmon_total_reads", recursive = TRUE, pattern = "(.+).csv", full.names = TRUE)
##read in
salmon_total_reads_list <- lapply(salmon_total_reads,
                             readr::read_csv, col_names = FALSE)
#change names of tables to sample name
names(salmon_total_reads_list) <- gsub(".*/(.+).csv", "\\1", salmon_total_reads)
##full results table
all_salmon_total_reads <- data.table(dplyr::bind_rows(salmon_total_reads_list, .id = "filename"))
##fix total reads column
all_salmon_total_reads$X1 <- tstrsplit(all_salmon_total_reads$X1, "Observed ", keep=c(2), fixed=TRUE)
all_salmon_total_reads$X1 <- tstrsplit(all_salmon_total_reads$X1, " total fragments", keep=c(1), fixed=TRUE)
setnames(all_salmon_total_reads, old=c("filename", "X1"), new=c("Sample_name", "Salmon_input_total_read_pairs"))
fwrite(all_salmon_total_reads, "output/salmon_total_reads/sample_total_read_pairs.csv")
