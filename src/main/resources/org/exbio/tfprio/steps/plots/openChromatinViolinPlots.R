#!/usr/bin/env Rscript

library("ggplot2")
library("scales")
library("argparse")

parser <- ArgumentParser()
parser$add_argument("--input", required = TRUE)
parser$add_argument("--output", required = TRUE)

args <- parser$parse_args()

hm_to_lengths <- read.csv(args$input, sep = '\t')
p <- ggplot(hm_to_lengths, aes(HM, LENGTH))
p <- p + geom_violin() + geom_boxplot(width = 0.1)
p <- p + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))
p <- p + theme(axis.title.x = element_blank())
p <- p +
  ylab("log10 length of open chromatin") +
  theme(text = element_text(size = 20))
ggsave(filename = args$output, plot = p)
