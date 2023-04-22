#!/usr/bin/env Rscript

library(pheatmap)
library(optparse)

option_list <- list(
    make_option(c("--counts"), action = "store", default = NA, type = "character", help = "Counts file"),
    make_option(c("--target_genes"), action = "store", default = NA, type = "character", help = "Target genes file"),
    make_option(c("--group1"), action = "store", default = NA, type = "character", help = "First group name"),
    make_option(c("--group2"), action = "store", default = NA, type = "character", help = "Second group name"),
    make_option(c("--hm"), action = "store", default = NA, type = "character", help = "Histone modification"),
    make_option(c("--samplesheet"), action = "store", default = NA, type = "character", help = "Sample sheet"),
    make_option(c("--map"), action = "store", default = NA, type = "character", help = "Map file")
)

args <- parse_args(OptionParser(option_list = option_list))

print("Reading counts")
counts <- read.csv(args$counts, sep = "\t", header = TRUE, row.names = 1)
print("Reading target genes")
target_genes <- read.csv(args$target_genes, sep = "\t", header = TRUE)
print("Reading samplesheet")
samplesheet <- read.csv(args$samplesheet, sep = ",", header = TRUE, row.names = 1)
print("Reading map")
map <- read.csv(args$map, sep = "\t", header = FALSE)
print("Done reading files")

# Keep only annotations where group is either group1 or group2
annotation <- samplesheet[(samplesheet$group == args$group1) | (samplesheet$group == args$group2), ]
annotation <- annotation[, "group", drop = FALSE]

# Keep only counts where rownames are in annotation
counts <- counts[, rownames(annotation)]
# Drop all rows containing only zeros
counts <- counts[rowSums(counts) != 0, ]

id_symbol <- setNames(map$V1, map$V2)

for (tf in colnames(target_genes)) {
    genes <- target_genes[, tf]
    genes <- intersect(rownames(counts), genes)
    tf_counts <- counts[genes, ]

    match_indices <- match(rownames(tf_counts), names(id_symbol))
    symbols <- ifelse(is.na(match_indices), rownames(tf_counts), id_symbol[match_indices])

    file <- paste0(args$group1, ":", args$group2, "_", args$hm, "_", tf, ".png")

    pheatmap(tf_counts, scale = "row", filename = file, labels_row = symbols, legend = TRUE, annotation_legend = TRUE,
         annotation_col = annotation, show_colnames = TRUE, show_col_dend = FALSE, show_row_dend = FALSE,
         treeheight_row = 0, treeheight_col = 0)
}
