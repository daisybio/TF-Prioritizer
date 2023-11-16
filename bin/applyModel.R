#!/usr/bin/env Rscript

library(ehmm)
library(rtracklayer)
library(argparser)
library(tools)
library(parallel)
library(dplyr)
library(bamsignals)


args <- commandArgs(trailingOnly = T)
parser <- arg_parser("EHMM command line parser", name = "ehmm")
parser <- add_argument(parser, "-r", help = "BED file of input regions")
parser <- add_argument(parser, "-g", help = "Genome sizes")
parser <- add_argument(parser, "-m", help = "Construced model")
parser <- add_argument(parser, "-b", help = "Path to directory containing the bam files of the input data")
parser <- add_argument(parser, "-o", help = "Output directory", default = "./")
parser <- add_argument(parser, "-s", help = "BinSize", default = 100, type = "integer")
parser <- add_argument(parser, "-t", help = "Number of threads to use", default = 21, type = "integer")
parser <- add_argument(parser, "-p", help = "Pseudocount", default = 1.0, type = "double")

argv <- parse_args(parser, argv = args)
message("loading input data")
binsize <- argv$s
nthreads <- argv$t
pseudocount <- argv$p
regions <- import(argv$r, format = "bed")
# remove scaffold chromosomes
levelsToDrop <- unique(unlist(lapply(c("Un", "M", "random",
                                       "hap", "alt", "GL", "NC", "hs"), function(x) which(grepl(x,
                                                                                          GenomeInfoDb:::seqlevels(regions))))))
if (length(levelsToDrop) > 0)
  GenomeInfoDb::seqlevels(regions, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(regions)[-levelsToDrop]
# load bam files
bamfiles <- list.files(path = argv$b, pattern = "*.bam$", full.names = T, recursive = T)
marks <- sapply(strsplit(file_path_sans_ext(basename(bamfiles)), "\\."), "[", 1)
marks[grep("ALL|RNA|ATC|DNASE|ATAC", marks, ignore.case = T)] <- "acc"
print(marks)
bamtab <- data.frame(mark = marks, path = bamfiles, check.names = F)
bamtab$shift <- sapply(tolower(bamtab$mark), function(m) ifelse(any(sapply(c("atac", "dhs", "dnase", "acc"),
                                                                           function(s) grepl(s, m))), 0, 75))
# generate counts from bam files
bamtab <- validateBamtab(bamtab)
npaths <- nrow(bamtab)
if (npaths == 0)
  stop("no marks provided")
if (!inherits(regions, "GRanges"))
  stop("regions must be a GRanges object")
if (length(binsize) != 1 || binsize <= 0)
  stop("invalid binsize")
if (any(width(regions)%%binsize!=0)) {
  # select invalid regions
  invalidRegions <- regions[width(regions)%%binsize!=0]
  message("INFO: ", length(invalidRegions), " regions are not multiples of the selected bin size ", binsize, ", adjusting...")
  deltaLengths <- width(invalidRegions)%%binsize
  # pad regions to binsize
  deltaLengths <- sapply(deltaLengths, function(d) ifelse(d>binsize/2, binsize-d, -d))
  # replace end points
  data <- as.data.frame(invalidRegions)
  data$end <- data$end+deltaLengths
  data$width <- data$width+deltaLengths
  data <- makeGRangesFromDataFrame(data, keep.extra.columns = T)
  regions[width(regions)%%binsize!=0] <- data
}
message("getting counts")
countsList <- safe_mclapply(mc.cores = nthreads, 1:npaths,
                            function(i) {
                              b <- bamtab[i, ]
                              b$pairedend <- ifelse(b$pairedend, "midpoint", "ignore")
                              bprof <- bamProfile(b$path, regions, binsize = binsize,
                                                  shift = b$shift, mapq = b$mapq, paired.end = b$pairedend)
                              unlist(as.list(bprof))
                            })
if (any(duplicated(bamtab$mark))) {
  message("putting replicate experiments together")
  countsList <- lapply(setNames(nm = unique(bamtab$mark)),
                       function(nm) {
                         subList <- countsList[bamtab$mark == nm]
                         if (length(subList) == 1)
                           return(subList[[1]])
                         vmat <- bindCols(subList)
                         defaultRepFun(vmat)
                       })
} else {
  names(countsList) <- bamtab$mark
}
message("pileup done, storing everything in a matrix")
counts <- t(bindCols(countsList)) + pseudocount
message("loading model")
model <- readModel(argv$m)
model$marks <- toupper(model$marks)
rownames(counts) <- toupper(rownames(counts))
commonMarks <- intersect(rownames(counts), model$marks)
if (length(commonMarks) < length(unique(c(rownames(counts), model$marks)))) {
  message("INFO: only using marks that are present in both training and testing set:\n",
          "-> ", paste(commonMarks, collapse = ", "))
  if(nrow(counts)>1) counts <- counts[commonMarks,]
  model$marks <- commonMarks
  if (length(model$emisP[[1]]$mus) > length(commonMarks)) {
    model$emisP <- lapply(model$emisP, lapply, "[", which(commonMarks==model$marks))
  }
}
message("loading genome sizes2")
genomeSizes <- import(argv$g, format = "bed")
genomeSizes <- setNames(genomeSizes@ranges@start, seqnames(genomeSizes))

applyModelInspect(regions, model = model, genomeSize = genomeSizes,
           bamtab = bamtab, counts = counts, nthreads = nthreads,
           outdir = argv$o)
message("ehmm finished successfully")
