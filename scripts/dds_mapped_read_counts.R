library(tximport)
library(data.table)
library(DESeq2)
library(textshape)
library(ggplot2)

##count reads that map onto asw and mh for each sample

##Import table describing samples
sample_data <- fread("data/sample_key.csv", header=TRUE)
setkey(sample_data, Sample_name)
##Find all salmon quant files
quant_files <- list.files(path="output/asw_mh_concat_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)

##asw dds
##asw tx2gene
asw_gene2tx <- fread("data/asw_edited_transcript_ids/Trinity.fasta.gene_trans_map", header = FALSE)
asw_tx2gene <- data.frame(asw_gene2tx[, .(V2, V1)])
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimaates e.g. salmon)
asw_txi <- tximport(quant_files, type = "salmon", tx2gene = asw_tx2gene, dropInfReps=TRUE)
##create asw_dds object and link to sample data  
asw_dds <- DESeqDataSetFromTximport(asw_txi, colData = sample_data[colnames(asw_txi$counts)], design = ~1)
saveRDS(asw_dds, file = "output/deseq2/asw_dds.rds")
#read back in
asw_dds <- readRDS("output/deseq2/asw_dds.rds")
##make asw counts matrix
asw_counts_matrix <- counts(asw_dds)
asw_counts_colsum <- data.table(data.frame(colSums(asw_counts_matrix)), keep.rownames = TRUE)

##check that counts returns un-normalised counts
##also need to get tabke of sample id vs total reads - needs to be total number provided to salmon not pre-cleanup


##mh dds
##mh tx2gene
mh_gene2tx <- fread("data/mh_edited_transcript_ids/Trinity.fasta.gene_trans_map", header = FALSE)
mh_tx2gene <- data.frame(mh_gene2tx[, .(V2, V1)])
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimaates e.g. salmon)
mh_txi <- tximport(quant_files, type = "salmon", tx2gene = mh_tx2gene, dropInfReps=TRUE)
##create mh_dds object and link to sample data  
mh_dds <- DESeqDataSetFromTximport(mh_txi, colData = sample_data[colnames(mh_txi$counts)], design = ~1)
saveRDS(mh_dds, file = "output/deseq2/mh_dds.rds")
##read back in
mh_dds <- readRDS("output/deseq2/mh_dds.rds")
##make mh counts matrix
mh_counts_matrix <- counts(mh_dds)
mh_counts_colsum <- data.table(data.frame(colSums(mh_counts_matrix)), keep.rownames = TRUE)

##merge asw and mh counts tables
asw_mh_mapped_reads <- merge(asw_counts_colsum, mh_counts_colsum, by="rn")
setnames(asw_mh_mapped_reads, old=c("rn", "colSums.asw_counts_matrix.", "colSums.mh_counts_matrix."), new=c("sample_id", "asw_counts_colsum", "mh_counts_colsum"))
##set sample_id column to rownames so can take rowSums
mapped_reads_df <- column_to_rownames(asw_mh_mapped_reads, "sample_id")
##sum of reads mapped to asw and mh transcriptomes - gives in read pairs
mapped_reads_df$total_mapped_reads <- rowSums(mapped_reads_df)
mapped_reads_dt <- setDT(mapped_reads_df, keep.rownames=TRUE)

##read in total read pairs from salmon log info
total_reads <- fread("output/salmon_total_reads/sample_total_read_pairs.csv")
sample_read_mapping <- merge(mapped_reads_dt, total_reads, by.x="rn", by.y="Sample_name")
##mapping rate
sample_read_mapping$`%_mapped` <- sample_read_mapping[,(((sample_read_mapping$total_mapped_reads)/(sample_read_mapping$Salmon_input_total_read_pairs))*100)]
##% of mapped reads that map to Mh
sample_read_mapping$`%_mapped_mh` <- sample_read_mapping[,(((sample_read_mapping$mh_counts_colsum)/(sample_read_mapping$total_mapped_reads))*100)]
##% of mapped reads that map to ASW
sample_read_mapping$`%_mapped_asw` <- sample_read_mapping[,(((sample_read_mapping$asw_counts_colsum)/(sample_read_mapping$total_mapped_reads))*100)]
fwrite(sample_read_mapping, "output/salmon_total_reads/sample_mapping.csv")

##merge with sample groups for plotting
sample_data <- fread("data/sample_key.csv")
sample_group <- sample_data[,c(2,5,6)]
sample_read_mapping_groups <- merge(sample_read_mapping, sample_group, by.x="rn", by.y="Sample_name")

##evasion subset
evasion_plot_data <- sample_read_mapping_groups[(sample_read_mapping_groups$Sample_group %in% c("E", "N")), ]
##plot
ggplot(evasion_plot_data, aes(x=Sample_group, y=`%_mapped_mh`, colour=Sample_group))+
  geom_boxplot(outlier.alpha = 0.2)+
  theme_light()+
  theme(legend.position = "none")+
  xlab("Sample group")+
  ylab("% of mapped reads")

##parasitism + heads subset
para_plot_data <- sample_read_mapping_groups[!(sample_read_mapping_groups$Sample_group %in% c("E", "N", "Abdo", "Thorax")), ]

##mh - order sample groups for graph - plotting para results
para_plot_data$Sample_group <- factor(para_plot_data$Sample_group, levels=c("Exposed Head", "NC Head", "Parasitised Abdomen", "NC Abdomen"))
##plot duplication per sample group
ggplot(para_plot_data, aes(x=Sample_group, y=`%_mapped_mh`, colour=Sample_group))+
  geom_boxplot(outlier.alpha = 0.2)+
  theme_light()+
  theme(legend.position = "none")+
  xlab("Sample group")+
  ylab("% of mapped reads")

##asw - plot duplication per sample group
ggplot(para_plot_data, aes(x=Sample_group, y=`%_mapped_asw`, colour=Sample_group))+
  geom_boxplot(outlier.alpha = 0.2)+
  theme_light()+
  theme(legend.position = "none")+
  xlab("Sample group")+
  ylab("% of mapped reads")
