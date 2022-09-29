library(ggplot2)
library(reshape2)

safe_colorblind_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

stats<-read.csv("{INPUTFILE}",sep="\t")

stats_visualized <- stats
stats_visualized$TP<-NULL
stats_visualized$TN<-NULL
stats_visualized$FP<-NULL
stats_visualized$FN<-NULL

stats_visualized<-melt(stats_visualized, id="TF",variable.name = "Metric",value.name = "Percentage")

p <- ggplot(stats_visualized, aes(fill=Metric, y=Percentage, x=TF)) +
    geom_bar(position="dodge", stat="identity")

p + scale_fill_manual(values=safe_colorblind_palette) + theme_bw()

ggsave("{OUTPUTFILE}", width=nrow(stats_visualized)/2, units="cm")