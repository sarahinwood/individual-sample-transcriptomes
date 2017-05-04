library(data.table)
library(ggplot2)
fread ("output/trinity_stats/xn50.out.txt")
data.raw<-fread ("output/trinity_stats/xn50.out.txt")
data.raw[,`#E`:=as.numeric(gsub("E", "", `#E`))]

ggplot(data.raw, aes(x=`#E`, y=`E-N50`)) + geom_point() +ggtitle("ASW Transcriptome") + scale_y_continuous(limits = c(0, 3000))
