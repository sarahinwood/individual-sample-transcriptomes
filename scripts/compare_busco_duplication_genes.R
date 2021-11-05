library(readr)
library(dplyr)
library(data.table)
library(ggplot2)
library(VennDiagram)

sample_data <- fread("data/sample_key.csv")
sample_group <- sample_data[,c(2,5,6)]

##find all busco length tables for evasive ASW samples - have been copied into one folder for ease
busco_fulltables <- list.files(path = "output/busco_full_tables", recursive = TRUE, pattern = "(.+)_full_table_length.tsv", full.names = TRUE)

##read in busco tables
busco_results_list <- lapply(busco_fulltables,
                             readr::read_tsv,
                             skip=4)

#change names of tables to sample name
names(busco_results_list) <- gsub(".*/(.+)_full_table_length.tsv", "\\1", busco_fulltables)

##full busco results table
busco_all <- dplyr::bind_rows(busco_results_list, .id = "filename")

##generate summary for each sample of each busco type
busco_summary <- busco_all %>%
  ##group by sample and busco status
  group_by(filename, Status) %>%
  ##summarise down to 1 value for each status
  summarise(n = n()) %>%
  ##make value a percentage
  mutate(percentage = n/sum(n)*100)

##filter for only duplicated
busco_duplicated <- dplyr::filter(busco_summary, Status == "Duplicated")
#write table of sample group with duplication level
busco_duplicated <- merge(busco_duplicated, sample_group, by.x="filename", by.y="Sample_name", all.x=TRUE)
fwrite(busco_duplicated, "output/busco_full_tables/sample_duplication_table.csv")

##order sample groups for graph - plotting para results
para_plot_data <- busco_duplicated[!(busco_duplicated$Sample_group %in% c("E", "N", "Abdo", "Thorax")), ]
##plot duplication per sample group - can change to parasitism status to show undetected abdo have low duplication
para_plot_data$Sample_group <- factor(para_plot_data$Sample_group, levels=c("Exposed Head", "NC Head", "Parasitised Abdomen", "NC Abdomen"))
ggplot(para_plot_data, aes(x=Sample_group, y=percentage, colour=Sample_group))+
  geom_boxplot(outlier.alpha = 0.2)+
  theme_light()+
  theme(legend.position = "none")+
  xlab("Sample group")+
  ylab("BUSCO Duplication %")

##for evasion samples
evasion_plot_data <- busco_duplicated[(busco_duplicated$Sample_group %in% c("E", "N")), ]
ggplot(evasion_plot_data, aes(x=Sample_group, y=percentage, colour=Sample_group))+
  geom_boxplot(outlier.alpha = 0.2)+
  theme_light()+
  theme(legend.position = "none")+
  xlab("Sample group")+
  ylab("BUSCO Duplication %")

##can use to save means or medians etc for each group
medians <- aggregate(percentage ~ Sample_group, busco_duplicated, median)
####################################################################################################



#make table of results for all samples
busco_all <- dplyr::bind_rows(busco_results_list, .id = "filename")
##3959 cases of duplicated BUSCO genes
busco_duplicated <- dplyr::filter(busco_all, Status == "Duplicated")
##drop columns I don't need
busco_duplicated <- busco_duplicated[-c(3,4,5,6)]
##only 663 BUSCO genes that are duplicated
length(unique(busco_duplicated$`# Busco id`))
##split table based on filename
busco_dup_split <- split(busco_duplicated, busco_duplicated$filename)

##compare duplicated BUSCO genes in para and ?para samples - cannot use more than 5 samples in venn
##- if I want to display this data it would be better to find a different method to do so
Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("L2_2h_A"=(busco_dup_split$`L2_2h_A`$`# Busco id`),
                            "L2_30m_A"=(busco_dup_split$`L2_30m_A`$`# Busco id`),
                            "L2_4h_A"=(busco_dup_split$`L2_4h_A`$`# Busco id`),
                            "R1_30m_A"=(busco_dup_split$`R1_30m_A`$`# Busco id`),
                            "R1_4h_A"=(busco_dup_split$`R1_4h_A`$`# Busco id`)),
                   filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd)
