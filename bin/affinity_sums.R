#!/usr/bin/env Rscript

library(optparse)

# Define the options
option_list <- list(
    make_option(c("-a", "--input1"),
        type = "character", default = NULL,
        help = "First input file path"
    ),
    make_option(c("-b", "--input2"),
        type = "character", default = NULL,
        help = "Second input file path"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "Output file path"
    )
)

# Parse the options
opt <- parse_args(OptionParser(option_list = option_list))

# Check if input and output paths are provided
if (is.null(opt$input1) || is.null(opt$input2) || is.null(opt$output)) {
    stop("Input1, input2, and output paths are required.")
}

# Read input files
df1 <- read.delim(opt$input1, row.names = 1)
df2 <- read.delim(opt$input2, row.names = 1)

# Create intersection of genes
intersection <- intersect(row.names(df1), row.names(df2))
df1 <- df1[intersection, ]
df2 <- df2[intersection, ]

# Check if the data frames have the same dimensions
if (!all(dim(df1) == dim(df2))) {
    stop("The input files must have the same dimensions.")
}

if (!all(rownames(df1) == rownames(df2))) {
    stop("The input files must have the same row names.")
}

if (!all(colnames(df1) == colnames(df2))) {
    stop("The input files must have the same column names.")
}

# Calculate the sum between the two data frames
sum <- df1 + df2

# Write the sum to a file
write.table(sum, file = opt$output,
    sep = "\t", row.names = TRUE, quote = FALSE)
