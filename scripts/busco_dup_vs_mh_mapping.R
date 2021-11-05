library(ggplot2)
library(data.table)

##busco dup data
busco_duplication <- fread("output/busco_full_tables/sample_duplication_table.csv")
##mh mapping
sample_mapping <- fread("output/salmon_total_reads/sample_mapping.csv")

plot_data <- merge(busco_duplication, sample_mapping, by.x="filename", by.y="rn")
plot_data <- plot_data[,c(1,4,12)]
##sample key for grouping
sample_data <- fread("data/sample_key.csv")
sample_group <- sample_data[,c(2,5,6)]
##merge with plot data
plot_data_grouped <- merge(plot_data, sample_group, by.x="filename", by.y="Sample_name")
fwrite(plot_data_grouped, "output/busco_dup_vs_mapping.csv")

##para plot
para_plot_data <- plot_data_grouped[!(plot_data_grouped$Sample_group %in% c("E", "N", "Abdo", "Thorax")), ]
para_plot_data$Sample_group <- factor(para_plot_data$Sample_group, levels=c("Exposed Head", "NC Head", "Parasitized Abdomen", "NC Abdomen"))
##mean of groups
means <- aggregate(`%_mapped_mh` ~ Parasitism_status, para_plot_data, mean)
##plot
ggplot(para_plot_data, aes(x=percentage, y=`%_mapped_mh`, colour=Sample_group, size=2, alpha=1))+
  geom_point()+
  theme_light()+
  xlab("BUSCO duplication (%)")+
  ylab("% of mapped reads")+
  scale_alpha(guide="none")+
  scale_size(guide="none")
##by para status
ggplot(para_plot_data, aes(x=percentage, y=`%_mapped_mh`, colour=Parasitism_status, size=2, alpha=1))+
  geom_point()+
  theme_light()+
  xlab("BUSCO duplication (%)")+
  ylab("% of mapped reads")+
  scale_alpha(guide="none")+
  scale_size(guide="none")

##evasion plot - doesn't subset
ev_plot_data <- plot_data_grouped[(plot_data_grouped$Sample_group %in% c("E", "N")), ]
ggplot(ev_plot_data, aes(x=percentage, y=`%_mapped_mh`, colour=Sample_group, size=2, alpha=1))+
  geom_point()+
  theme_light()+
  xlab("BUSCO duplication (%)")+
  ylab("% of mapped reads")+
  scale_alpha(guide="none")+
  scale_size(guide="none")
