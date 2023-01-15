#!/usr/bin/env Rscript

library(DESeq2)
library(pheatmap)
library(dplyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--counts", required = TRUE)
parser$add_argument("--tg1", required = TRUE)
parser$add_argument("--tg2", required = TRUE)
parser$add_argument("--heatmap", required = TRUE)
parser$add_argument("--outCounts", required = TRUE)
parser$add_argument("--groups", required = TRUE)
args <- parser$parse_args()

counts <- read.csv(args$counts, sep = "\t", header = TRUE, row.names = 1)
tg1 <- read.csv(args$tg1, sep = "\t", header = TRUE, row.names = 1)
tg2 <- read.csv(args$tg2, sep = "\t", header = TRUE, row.names = 1)
groups <- read.csv(args$groups, sep = "\t", header = TRUE, row.names = 1)

# Concatenate row names of tg1 and tg2
genes <- c(rownames(tg1), rownames(tg2))
genes <- intersect(genes, rownames(counts))
counts <- counts[genes,]

groups <- groups[, "group", drop = FALSE]

write.table(counts, file = args$outCounts, row.names = TRUE, col.names = TRUE, sep = "\t")
pheatmap(counts, scale = "row", filename = args$heatmap, legend = TRUE, annotation_legend = TRUE,
         annotation_col = groups, show_colnames = FALSE, show_col_dend = FALSE, show_row_dend = FALSE,
         treeheight_row = 0, treeheight_col = 0)