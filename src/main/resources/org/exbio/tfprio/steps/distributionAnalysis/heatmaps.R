#!/usr/bin/env Rscript

library(DESeq2)
library(pheatmap)
library(dplyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--affinity_file1", required = TRUE)
parser$add_argument("--affinity_file2", required = TRUE)
parser$add_argument("--expression_file", required = TRUE)
parser$add_argument("--heatmap_file", required = TRUE)
parser$add_argument("--metadata", required = TRUE)
args <- parser$parse_args()

tf <- tools::file_path_sans_ext(basename(args$target_genes_file))
pairing <- tools::file_path_sans_ext(group_pairing)

metaData <- read.csv(args$metadata, sep = "\t")

target_genes1 <- read.csv(args$affinity_file1, sep = "\t")
target_genes2 <- read.csv(args$affinity_file2, sep = "\t")
target_genes <- rbind(target_genes1, target_genes2)

expression <- read.csv(args$expression_file, sep = "\t")
target_gene_expression <- expression[expression$gene_id %in% target_genes$gene_id, ]

read_counts <- subset(chosen, select = -c(Geneid))
metaData_groupPairing <- metaData[colnames(read_counts), ]
metaData_groupPairing$sample_id <- NULL

dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = metaData_groupPairing,
                              design = ~{ DESIGN })
dds <- DESeq(dds, quiet = TRUE)

geneCounts_normalized <- counts(dds)

pheatmap(geneCounts_normalized, scale = "row", filename = plot_file, labels_row = chosen$geneSymbol, legend =
 TRUE, annotation_legend = TRUE, angle_col = 45, show_row_dend = FALSE, show_col_dend = FALSE, show_rownames
 = TRUE, show_colnames = FALSE, treeheight_row = 0, treeheight_col = 0, annotation_col = metaData_groupPairing)