library(rtracklayer)
library(GenomicRanges)
library(argparser)

args <- commandArgs(trailingOnly = T)
parser <- arg_parser("EHMM command line parser", name = "ehmm")
parser <- add_argument(parser, "-d", help = "Directory containing applyModel output")
parser <- add_argument(parser, "-o", help = "Output directory")
parser <- add_argument(parser, "-t", help = "Filtering ehmm score threshold", default = 0.9)
argv <- parse_args(parser, argv = args)
# find enhancer and promoter files
promoterFiles <- list.files(argv$d, pattern = "promoterRegions.bed", recursive = T, full.names = T)
enhancerFiles <- list.files(argv$d, pattern = "enhancerRegions.bed", recursive = T, full.names = T)
# import all ranges
promoterRegions <- sort(do.call(c, lapply(promoterFiles, import, format = "bed")))
enhancerRegions <- sort(do.call(c, lapply(enhancerFiles, import, format = "bed")))
# filter regions for score threshold
promoterRegions <- promoterRegions[promoterRegions$score >= argv$t]
enhancerRegions <- enhancerRegions[enhancerRegions$score >= argv$t]
message("predicted ", length(promoterRegions), " promoters and ", length(enhancerRegions), " enhancers.")
# write files to output
export(promoterRegions, file.path(argv$o, "predictedPromoters.bed"), format = "bed")
export(enhancerRegions, file.path(argv$o, "predictedEnhancers.bed"), format = "bed")
