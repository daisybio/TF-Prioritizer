library(rtracklayer)
library(GenomicRanges)
library(argparser)

args <- commandArgs(trailingOnly = T)
parser <- arg_parser("EHMM command line parser", name = "ehmm")
parser <- add_argument(parser, "-d", help = "Directory containing applyModel output")
parser <- add_argument(parser, "-o", help = "Output directory")
argv <- parse_args(parser, argv = args)
# find enhancer and promoter files
promoterFiles <- list.files(argv$d, pattern = "promoterRegions.bed", recursive = T, full.names = T)
enhancerFiles <- list.files(argv$d, pattern = "enhancerRegions.bed", recursive = T, full.names = T)
# import all ranges
promoterRegions <- sort(do.call(c, lapply(promoterFiles, import, format = "bed")))
enhancerRegions <- sort(do.call(c, lapply(enhancerFiles, import, format = "bed")))
# write files to output
export(promoterRegions, file.path(argv$o, "predictedPromoters.bed"), format = "bed")
export(enhancerRegions, file.path(argv$o, "predictedEnhancers.bed"), format = "bed")
