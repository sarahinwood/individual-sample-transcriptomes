library(data.table)

abundance_file <- "output/trinity_abundance/RSEM.isoforms.results"
isoform.list <- fread(abundance_file)
##subset into only those above 500bp
above_500_list <- subset(x=isoform.list, subset = isoform.list$length>500)
#sort by length to get longest isoform for each gene above 500bp
isoforms_by_length <- above_500_list[,.I[which.max(length)], by=gene_id][,V1]
fwrite(above_500_list[isoforms_by_length,list(transcript_id)], "output/length_filtered/transcript_ids_above_500.txt", col.names = FALSE)

#Effective length is the mean number of 5' start positions from which an RNA-Seq fragment could have been derived from
#this transcript, given the distribution of fragment lengths inferred by RSEM — the value is equal to
#(transcript_length - mean_fragment_length + 1)