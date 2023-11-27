#!/usr/bin/env Rscript

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
print("Loading background model")
load(argv$b)
bgModel <- model
bgCounts <- counts
bgRegions <- regions

print("Loading promoter model")
load(argv$p)
pModel <- model
pCounts <- counts
pRegions <- regions

print("Selecting promoter states")
pStates <- selectStates(pModel)

print("Loading enhancer model")
load(argv$e)
eModel <- model
eCounts <- counts
eRegions <- regions

print("Selecting enhancer states")
eStates <- selectStates(eModel)

print("Constructing combined model...\n")
nthreads = argv$t

constructModel(model.bg = bgModel, model.e = eModel, model.p = pModel,
                             counts.e = eCounts, counts.p = pCounts,
                             regions.e = eRegions, regions.p = pRegions,
                             nthreads = nthreads, outdir = argv$o)
