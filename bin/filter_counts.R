#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
    make_option(c("-i", "--input"), default = NULL, type = "character", metavar = "string", help = "input file path"),
    make_option(c("--output_tpm"), default = "tpm.tsv", type = "character", metavar = "string", help = "TPM output file"),
    make_option(c("--output_count"), default = "counts.tsv", type = "character", metavar = "string", help = "Count output file"),
    make_option(c("-c", "--min_count"), default = NULL, type = "integer", metavar = "integer", help = "Count cutoff"),
    make_option(c("-t", "--min_tpm"), default = NULL, type = "double", metavar = "double", help = "TPM cutoff")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    stop("Input file needs to be specified")
}

if (is.null(opt$output_tpm)) {
    stop("TPM output file needs to be specified")
}

if (is.null(opt$output_count)) {
    stop("Count output file needs to be specified")
}

if (is.null(opt$min_count) && is.null(opt$min_tpm)) {
    stop("At least one of min_count and min_tpm needs to be specified")
}

count_df <- read.csv(opt$input, sep = "\t", header = TRUE, row.names = 1)
count_df <- subset(count_df, select = -gene_name)
# count_df <- round(count_df)

if (!is.null(opt$min_count)) {
    qualifying_rows <- apply(count_df, 1, sum) > opt$min_count
    count_df <- count_df[qualifying_rows, ]
}

# Calculate TPM values
col_sums <- colSums(count_df)
rel_abundance <- t(count_df) / col_sums
scaling_factor <- 1e6 / apply(rel_abundance, 2, mean)
tpm_df <- t(rel_abundance * scaling_factor)

if (!is.null(opt$min_tpm)) {
    qualifying_rows <- apply(tpm_df, 1, sum) > opt$min_tpm
    count_df <- count_df[qualifying_rows, ]
    tpm_df <- tpm_df[qualifying_rows, ]
}

write.table(count_df, file = opt$output_count, sep = "\t", row.names = TRUE, quote = FALSE)
write.table(tpm_df, file = opt$output_tpm, sep = "\t", row.names = TRUE, quote = FALSE)
