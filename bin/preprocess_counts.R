#!/usr/bin/env Rscript

library("DESeq2")
library("optparse")

option_list <- list(
  make_option(c("-m", "--metadata"), default=NULL, type="character", metavar="string", help="metadata file path"),
  make_option(c("-i", "--input"), default=NULL, type="character", metavar="string", help="input file path"),
  make_option(c("-o", "--output"), default=NULL, type="character", metavar="string", help="output file path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(is.null(opt$metadata) || is.null(opt$input) || is.null(opt$output)){
  stop("You need to specify all three options: metadata, input, and output")
}

print("Reading metadata...")

metadata_df <- read.csv(opt$metadata, sep=",", header=TRUE, row.names=1)
metadata_df$group <- as.factor(metadata_df$group)
metadata_df$batch <- as.factor(metadata_df$batch)

include_batches <- length(unique(metadata_df$batch)) > 1

print("Reading input...")

count_df <- read.csv(opt$input, sep="\t", header=TRUE, row.names=1)
count_df <- subset(count_df, select = -gene_name)
count_df <- round(count_df)

# Make sure the order of the columns in the metadata and count dataframes match
metadata_df <- metadata_df[match(colnames(count_df), rownames(metadata_df)), ]


print("Creating DESeq dataset...")
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

print("Running DESeq...")
dds <- DESeq(dds)
counts_normalized <- counts(dds)
write.table(counts_normalized, file="counts_normalized.tsv", sep="\t", row.names=TRUE, quote=FALSE)

print("Grouping samples...")

group_df <- data.frame(row.names=rownames(counts_normalized))

# Loop over each unique group
for (g in unique(metadata_df$group)) {
  # Subset the count_df to the samples in the current group
  # Attention: Requires counts_normalized and metadata_df to have the same order
  group_counts <- counts_normalized[, metadata_df$group == g]
  
  # Calculate the row means for the current group
  row_means <- rowMeans(group_counts)
  
  # Add a new column to group_df with the mean counts for the current group
  group_df[g] <- row_means
}

print("Writing output...")
write.table(group_df, file=opt$output, sep="\t", row.names=TRUE, quote=FALSE)