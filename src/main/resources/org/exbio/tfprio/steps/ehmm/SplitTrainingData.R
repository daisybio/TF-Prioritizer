library(ehmm)
library(argparser)
library(rtracklayer)
library(GenomicRanges)

filterRegions <- function(regions, w = 2000) {
  regions <- reduce(regions)
  message("filter enhancer and promoter regions for a length of < 2 kb")
  regions <- regions[width(regions) < w]
  d <- distance(regions[1:(length(regions)-1)], regions[2:length(regions)])
  message("remove neighboring regions with a pairwise distance of < 2 kb")
  regions[d >= w | is.na(d)]
}

args <- commandArgs(trailingOnly = T)
parser <- arg_parser("Filtering command line parser", name = "ehmm_training_filter")
parser <- add_argument(parser, "-e", help = "enhancer bed file")
parser <- add_argument(parser, "-p", help = "promoter bed file")
parser <- add_argument(parser, "-b", help = "background bed file")
parser <- add_argument(parser, "-o", help = "output directory", default = "./")
argv <- parse_args(parser, argv = args)

message("load background, enhancer, and promoter regions")
bgRegions <- import(argv$b, format = "bed")
names(bgRegions) <- 1:length(bgRegions)
eRegions <- import(argv$e, format = "bed")
pRegions <- import(argv$p, format = "bed")
message("extract enhancer and promoter regions within background regions")
eAndBg <- findOverlapPairs(bgRegions, eRegions)
pAndBg <- findOverlapPairs(bgRegions, pRegions)
eBgRegions <- filterRegions(eAndBg@first)
pBgRegions <- filterRegions(pAndBg@first)
message("removing enhancers and promoters from background regions")
bgRegions <- reduce(bgRegions[!names(bgRegions) %in% unique(c(names(eAndBg@first), names(pAndBg@first)))])
message("writing filtered bed files")
export.bed(bgRegions, file.path(argv$o, "background.bed"))
export.bed(eBgRegions, file.path(argv$o, "enhancers.bed"))
export.bed(pBgRegions, file.path(argv$o, "promoters.bed"))
