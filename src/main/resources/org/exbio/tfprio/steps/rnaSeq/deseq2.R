#!/usr/bin/env Rscript

library("DESeq2")
library("argparse")
library("dplyr")

parser <- ArgumentParser()
parser$add_argument("--metadata", required = TRUE)
parser$add_argument("--input", required = TRUE)
parser$add_argument("--output", required = TRUE)

args <- parser$parse_args()

metadata_df <- read.csv(args$metadata, sep = "\t", header = T, row.names = NULL)
rownames(metadata_df) <- metadata_df$sample
metadata_df$sample <- NULL
metadata_df$group <- as.factor(metadata_df$group)
metadata_df$batch <- as.factor(metadata_df$batch)

include_batches <- (n_distinct(metadata_df$batch) > 1)

count_df <- read.csv(args$input, sep = "\t", header = T, row.names = NULL)
count_df <- count_df[!duplicated(count_df$gene_id), ]
row.names(count_df) <- count_df[, 1]
count_df$gene_id <- NULL

count_df <- round(count_df)

if (include_batches) {
  dds <- DESeqDataSetFromMatrix(
    countData = count_df,
    colData = metadata_df,
    design = ~ batch + group
  )
} else {
  dds <- DESeqDataSetFromMatrix(
    countData = count_df,
    colData = metadata_df,
    design = ~group
  )
}
threshold <- 50
keep <- rowSums(counts(dds)) >= threshold
dds <- dds[keep, ]
output_path <- args$output
dds <- DESeq(dds)
res <- results(dds)
summary(res)
write.table(res, file = output_path, sep = "\t")
