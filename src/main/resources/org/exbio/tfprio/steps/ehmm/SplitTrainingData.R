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

randomRanges <- function(regions, n, size) {
  if(!is.null(n) & n<=length(regions)) rRange <- regions[sample(1:length(regions), n)]
  else rRange <- regions
  rShift <- width(rRange)-size
  rShift <- sapply(rShift, function(shift) sample(0:shift, 1))
  rRange@ranges@start <- rRange@ranges@start+rShift
  width(rRange) <- size
  rRange
}

generateBackground <- function(gtf, e, p, size, n=100000) {
  gr <- import(gtf, format = "gtf")
  genic <- gr[gr$type=="gene"]
  intergenic <- gaps(genic)
  genic <- genic[width(genic)>=size]
  intergenic <- intergenic[width(intergenic)>=size]
  # get random ranges
  rIntergenic <- randomRanges(intergenic, n*0.9, size)
  rGenic <- randomRanges(genic, n*0.1, size)
  bg <- c(rIntergenic, rGenic)
  bg$score[is.na(bg$score)] <- 666
  # remove chr prefix from seqlevels
  seqlevels(bg) <- gsub("chr", "", seqlevels(bg))
  fg <- reduce(c(e, p))
  message("removing enhancer and promoter regions from background")
  subsetByOverlaps(bg, fg, invert = T)
}

splitRanges <- function(ranges, trainSplit, nSamples) {
  N <- length(ranges)
  if(nSamples == -1) {
    message("splitting data by a factor of ", trainSplit, " for training")
    nTrain <- N*trainSplit
  } else {
    message("Randomly selecting ", nSamples, " regions for training")
    nTrain <- nSamples
  }
  nTest <- N-nTrain
  R <- 1:N
  trainIdx <- sample(R, nTrain)
  testIdx <- R[-trainIdx]
  list(train=reduce(sort(ranges[trainIdx])), test=reduce(sort(ranges[testIdx])))
}

args <- commandArgs(trailingOnly = T)
parser <- arg_parser("Filtering command line parser", name = "ehmm_training_filter")
parser <- add_argument(parser, "-e", help = "enhancer bed file")
parser <- add_argument(parser, "-p", help = "promoter bed file")
parser <- add_argument(parser, "-b", help = "background bed file")
parser <- add_argument(parser, "-g", help = "path to gtf file")
parser <- add_argument(parser, "-o", help = "output directory", default = "./")
parser <- add_argument(parser, "-t", help = "Percent of data to use for training",
                       default = 0.8, type = "double")
parser <- add_argument(parser, "-f", help = "Fixed number of samples to consider for each training set",
                      default = -1, type = "integer")
parser <- add_argument(parser, "-s", help = "Size of random genomic regions",
                       default = 2000, type = "integer")
parser <- add_argument(parser, "-r", help = "seed",
                       default = 74726, type = "integer")
argv <- parse_args(parser, argv = args)
set.seed(argv$r)
message("load ChipAtlas, enhancer, and promoter regions")
bgRegions <- import(argv$b, format = "bed")
names(bgRegions) <- 1:length(bgRegions)
eRegions <- import(argv$e, format = "bed")
pRegions <- import(argv$p, format = "bed")
message("extract enhancer and promoter regions within ChipAtlas regions")
eAndBg <- findOverlapPairs(bgRegions, eRegions)
pAndBg <- findOverlapPairs(bgRegions, pRegions)
eBgRegions <- filterRegions(eAndBg@first)
pBgRegions <- filterRegions(pAndBg@first)
message("generating random background data from gft file")
bgNoCREs <- generateBackground(argv$g, eBgRegions, pBgRegions, argv$s)
export.bed(bgNoCREs, file.path(argv$o, "all_background.bed"))
export.bed(eBgRegions, file.path(argv$o, "all_enhancers.bed"))
export.bed(pBgRegions, file.path(argv$o, "all_promoters.bed"))
bgNoCREs <- splitRanges(bgNoCREs, argv$t, argv$f)
eBgRegions <- splitRanges(eBgRegions, argv$t, argv$f)
pBgRegions <- splitRanges(pBgRegions, argv$t, argv$f)
message("writing filtered train bed files")
export.bed(bgNoCREs$train, file.path(argv$o, "trainBackground.bed"))
export.bed(eBgRegions$train, file.path(argv$o, "trainEnhancers.bed"))
export.bed(pBgRegions$train, file.path(argv$o, "trainPromoters.bed"))
message("writing filtered test bed files")
export.bed(bgNoCREs$test, file.path(argv$o, "testBackground.bed"))
export.bed(eBgRegions$test, file.path(argv$o, "testEnhancers.bed"))
export.bed(pBgRegions$test, file.path(argv$o, "testPromoters.bed"))
