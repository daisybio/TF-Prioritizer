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

# Merge all data frames by common columns
merged <- Reduce(function(x, y) merge(x, y, all = TRUE), dfs)

# Calculate the mean of each cell
means <- apply(merged, c(1, 2), mean, na.rm = TRUE)

# Write the means to a file
write.table(means, opt$output, sep = "\t", row.names = FALSE)
