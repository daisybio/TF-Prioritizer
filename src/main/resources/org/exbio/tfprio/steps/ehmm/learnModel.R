library(ehmm)
library(rtracklayer)
library(argparser)
library(tools)
library(parallel)
library(dplyr)
library(bamsignals)

statesToSegments_helper <- function (regions, states) {
  .Call("_ehmm_statesToSegments_helper", PACKAGE = "ehmm",
        regions, states)
}

statesToSegments <- function (states, regions) {
  h <- statesToSegments_helper(regions, states)
  GRanges(seqnames = h$chrs, IRanges(start = h$starts, end = h$ends,
                                     names = h$states))
}

validateVmat <- function (vmat)
{
  if (!is.matrix(vmat))
    stop("vmat must be a matrix")
  if (!is.numeric(vmat))
    stop("'vmat' must be numeric")
  if (length(vmat) <= 0)
    stop("'vmat' must be non-empty")
}

quantileNormalization <- function (vmat, ref = c("median", "min", "mean"), nthreads = 1)
{
  validateVmat(vmat)
  if (!is.numeric(ref)) {
    ref <- match.arg(ref)
    ref <- getRef(vmat, ref, nthreads)
  }
  if (length(ref) != nrow(vmat))
    stop("reference vector has the wrong length")
  quantileNorm(vmat, ref, nthreads = nthreads)
}

propagateErrors <- function (l)
{
  for (el in l) if (inherits(el, "try-error"))
    stop(el)
  l
}

safe_mclapply <- function (...)
  propagateErrors(mclapply(...))

colSummary <- function (mat, type, nthreads = 1L)
{
  .Call("_ehmm_colSummary", PACKAGE = "ehmm", mat, type, nthreads)
}

bindCols <- function (vlist, nthreads = 1L)
{
  .Call("_ehmm_bindCols", PACKAGE = "ehmm", vlist, nthreads)
}

defaultRepFun <- function (vmat, normFun = quantileNormalization, nthreads = 1,
          ...)
{
  normFunOpts <- list(...)
  normFunOpts$vmat <- vmat
  if ("nthreads" %in% names(formals(normFun))) {
    normFunOpts[["nthreads"]] <- nthreads
  }
  nvmat <- do.call(normFun, normFunOpts)
  colSummary(t(nvmat), "median", nthreads = nthreads)
}

bamtabDefaults <- list(75, 0, FALSE)
names(bamtabDefaults) <- c("shift", "mapq", "pairedend")

validateBamtab <- function (bamtab)
{
  if (!is.data.frame(bamtab))
    stop("'bamtab' must be a 'data.frame'")
  reqfields <- c("mark", "path")
  optfields <- "shift"
  if (!all(reqfields %in% names(bamtab)))
    stop("missing required fields")
  if (!all(names(bamtab) %in% c(reqfields, optfields)))
    stop("invalid fields")
  for (n in names(bamtabDefaults)) {
    if (!n %in% names(bamtab)) {
      bamtab[[n]] <- rep(bamtabDefaults[[n]], nrow(bamtab))
    }
  }
  if (!is.character(bamtab$path))
    stop("invalid path specification")
  if (any(!file.exists(bamtab$path)))
    stop("BAM file does not exist")
  if (!is.numeric(bamtab$mapq) || any(bamtab$mapq < 0 | bamtab$mapq >
                                      255)) {
    stop("invalid 'mapq'")
  }
  if (!is.numeric(bamtab$shift))
    stop("invalid 'shift'")
  if (!is.logical(bamtab$pairedend))
    stop("invalid 'pairedend'")
  bamtab
}

args <- commandArgs(trailingOnly = T)
parser <- arg_parser("EHMM command line parser", name = "ehmm")
parser <- add_argument(parser, "-r", help = "Regions of interest")
parser <- add_argument(parser, "-m", help = "Path to directory containing the bam files of the input data")
parser <- add_argument(parser, "-n", help = "Number of states", default = 10, type = "integer")
parser <- add_argument(parser, "-b", help = "Size of bins", default = 10, type = "integer")
parser <- add_argument(parser, "-t", help = "Number of threads", default = 1, type = "integer")
parser <- add_argument(parser, "-p", help = "Pseudocount", default = 1.0, type = "double")
parser <- add_argument(parser, "-o", help = "Output directory", default = "./")

argv <- parse_args(parser, argv = args)

regions <- import(argv$r, format = "bed")
# remove scaffold chromosomes
regions <- regions[!grepl("_", regions@seqnames)]
# change M to MT
regions@seqnames@values <- as.factor(gsub("M", "MT", regions@seqnames@values))
bamfiles <- list.files(path = argv$m, pattern = "*.bam$", full.names = T)
marks <- sapply(strsplit(file_path_sans_ext(basename(bamfiles)), "_"), "[", 1)
bamtab <- data.frame(mark = marks, path = bamfiles, check.names = F)
nstates = argv$n
binsize = argv$b
nthreads = argv$t
pseudocount = argv$p

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
  message("INFO: ", nrow(invalidRegions), " regions are multiples of the selected bin size ", binsize, ", adjusting regions...")
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
names(countsList) <- bamtab$mark
message("pileup done, storing everything in a matrix")
counts <- t(bind_cols(countsList)) + pseudocount

cat("learning model\n")
segmentation <- segment(counts = counts, regions = regions,
                        nstates = nstates, nthreads = nthreads,
                        verbose_kfoots = T,
                        nbtype = "lognormal")
cat("producing report\n")
viterbi_segments <- statesToSegments(segmentation$viterbi,
                                     segmentation$segments)
report(segments = viterbi_segments, model = segmentation$model,
       rdata = segmentation, outdir = argv$o)
model <- segmentation$model
cat("saving essential RData to file\n")
save(counts, model, regions, file = file.path(argv$o, "model.RData"))
