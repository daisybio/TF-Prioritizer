library(ehmm)
library(rtracklayer)
library(argparser)
library(RColorBrewer)
library(Matrix)

args <- commandArgs(trailingOnly = T)
parser <- arg_parser("EHMM command line parser", name = "ehmm")
parser <- add_argument(parser, "-b", help = "Background model RData")
parser <- add_argument(parser, "-e", help = "Enhancer model RData")
parser <- add_argument(parser, "-p", help = "Promoter model RData")
parser <- add_argument(parser, "-t", help = "Number of threads", default = 21, type = "integer")
parser <- add_argument(parser, "-o", help = "Output directory", default = "./")

argv <- parse_args(parser, argv = args)
cat("loading models and counts...\n")
load(argv$b)
bgModel <- model
bgCounts <- counts
bgRegions <- regions
load(argv$e)
eModel <- model
eCounts <- counts
eRegions <- regions
eStates <- selectStates(eModel)

load(argv$p)
pModel <- model
pCounts <- counts
pRegions <- regions
pStates <- selectStates(pModel)
cat("constructing combined model...\n")
nthreads = argv$t

model <- constructModelInspect(model.bg = bgModel, model.e = eModel, model.p = pModel,
                             counts.e = eCounts, counts.p = pCounts,
                             regions.e = eRegions, regions.p = pRegions,
                             nthreads = nthreads, outdir = argv$o)
