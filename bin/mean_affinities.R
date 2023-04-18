#!/usr/bin/env Rscript

library(optparse)

# Define the options
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "List of input file paths"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "Output file path"
    )
)

# Parse the options
opt <- parse_args(OptionParser(option_list = option_list))

# Check if input and output paths are provided
if (is.null(opt$input) || is.null(opt$output)) {
    stop("Input and output paths are required.")
}

# Read all input files into a list
dfs <- lapply(opt$input, read.delim)

# Check if all data frames have the same dimensions
if (!all(sapply(dfs, function(df) all.equal(dim(df), dim(dfs[[1]]))))) {
    stop("The input files must have the same dimensions.")
}

# Check if all data frames have the same row names
if (!all(sapply(dfs, function(df) all.equal(rownames(df), rownames(dfs[[1]]))))) {
    stop("The input files must have the same row names.")
}

# Check if all data frames have the same column names
if (!all(sapply(dfs, function(df) all.equal(colnames(df), colnames(dfs[[1]]))))) {
    stop("The input files must have the same column names.")
}

# Calculate the mean of each cell
means <- Reduce("+", dfs) / length(dfs)

# Write the means to a file
write.table(means, file = opt$output, sep = "\t", row.names = FALSE)
