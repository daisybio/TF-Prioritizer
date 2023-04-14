#!/usr/bin/env Rscript

library("DESeq2")
library("optparse")

option_list <- list(
    make_option(c("-m", "--metadata"), default = NULL, type = "character", metavar = "string", help = "metadata file path"),
    make_option(c("-i", "--input"), default = NULL, type = "character", metavar = "string", help = "input file path"),
    make_option(c("--group1"), default = NULL, type = "character", metavar = "string", help = "First group to consider"),
    make_option(c("--group2"), default = NULL, type = "character", metavar = "string", help = "Second group to consider")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$metadata) || is.null(opt$input)) {
    stop("You need to specify both options: metadata and input")
}

if (is.null(opt$group1) || is.null(opt$group2)) {
    stop("You need to specify both groups")
}

print("Reading metadata...")

metadata_df <- read.csv(opt$metadata, sep = ",", header = TRUE, row.names = 1)
metadata_df$group <- as.factor(metadata_df$group)

print("Reading input...")

count_df <- read.csv(opt$input, sep = "\t", header = TRUE, row.names = 1)
# count_df <- round(count_df)

# Make sure the order of the columns in the metadata and count dataframes match
metadata_df <- metadata_df[match(colnames(count_df), rownames(metadata_df)), ]

pairing <- c(opt$group1, opt$group2)
name <- paste(pairing, collapse = ":")

current_metadata <- metadata_df[metadata_df$group %in% pairing, ]
current_dataset <- count_df[, na.omit(match(colnames(count_df), rownames(current_metadata)))]

dds <- DESeqDataSetFromMatrix(
    countData = current_dataset,
    colData = current_metadata,
    design = ~group
)

dds <- DESeq(dds)
res <- results(dds)
write.table(res, file = paste0(name, ".tsv"), quote = FALSE)
