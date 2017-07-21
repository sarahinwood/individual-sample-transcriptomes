library(data.table)
isoform.list <- fread('output/trinity_abundance/RSEM.isoforms.results')
exp.rows <- isoform.list[,.I[which.max(IsoPct)], by=gene_id][,V1]
fwrite(isoform.list[exp.rows,list(transcript_id)], 'output/trinity_abundance/isoforms_by_expression.txt', col.names = FALSE)


length.rows <- isoform.list[,.I[which.max(length)], by=gene_id][,V1]
isoform.list[length.rows,transcript_id]
fwrite(isoform.list[length.rows,list(transcript_id)], 'output/trinity_abundance/isoforms_by_length.txt', col.names = FALSE)
