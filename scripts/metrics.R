library(ggplot2)
library(reshape2)

stats<-read.csv("{INPUTFILE}",sep="\t")

stats_visualized <- stats
stats_visualized$TP<-NULL
stats_visualized$TN<-NULL
stats_visualized$FP<-NULL
stats_visualized$FN<-NULL

stats_visualized<-melt(stats_visualized, id="TF",variable.name = "Metric",value.name = "Percentage")

ggplot(stats_visualized, aes(fill=Metric, y=Percentage, x=TF)) +
    geom_bar(position="dodge", stat="identity")
ggsave("{OUTPUTFILE}", width=nrow(stats_visualized)/2, units="cm")