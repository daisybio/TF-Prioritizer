library(ehmm)
library(rtracklayer)
library(argparser)
library(tools)

args <- commandArgs(trailingOnly = T)
parser <- arg_parser("EHMM command line parser", name = "ehmm")
parser <- add_argument(parser, "-r", help = "Regions of interest")
parser <- add_argument(parser, "-m", help = "Path to directory containing the bam files of the input data")
parser <- add_argument(parser, "-n", help = "Number of states", type = "integer")
parser <- add_argument(parser, "-b", help = "Number of bins", default = 10, type = "integer")
parser <- add_argument(parser, "-t", help = "Number of threads", default = 1, type = "integer")
parser <- add_argument(parser, "-p", help = "Pseudocount", default = 1)
parser <- add_argument(parser, "-o", help = "Output directory", default = "./")

argv <- parse_args(parser, argv = args)

regions <- import(argv$r, format = "bed")
bamfiles <- list.files(path = argv$m, pattern = "*.bam$", full.names = T)
marks <- sapply(strsplit(file_path_sans_ext(basename(bamfiles)), "_"), "[", 1)
bamtab <- data.frame(mark = marks, path = bamfiles, check.names = F)
nstates = argv$n
binsize = argv$b
nthreads = argv$t
pseudocount = argv$p

statesToSegments_helper <- function (regions, states) {
  .Call("_ehmm_statesToSegments_helper", PACKAGE = "ehmm", 
        regions, states)
}

statesToSegments <- function (states, regions) {
  h <- statesToSegments_helper(regions, states)
  GRanges(seqnames = h$chrs, IRanges(start = h$starts, end = h$ends, 
                                     names = h$states))
}

bamtab$shift <- sapply(tolower(bamtab$mark), function(m) ifelse(any(sapply(c("atac", "dhs", "dnase", "acc"), 
                                                                           function(s) grepl(s, m))), 0, 75))
counts <- getcounts(regions, bamtab, binsize = binsize, 
                    nthreads = nthreads) + pseudocount

cat("learning model\n")
segmentation <- segment(counts = counts, regions = regions, 
                        nstates = nstates, nthreads = nthreads, verbose_kfoots = TRUE, 
                        nbtype = "lognormal")
cat("producing report\n")
viterbi_segments <- statesToSegments(segmentation$viterbi, 
                                     segmentation$segments)
report(segments = viterbi_segments, model = segmentation$model, 
       rdata = segmentation, outdir = argv$o)
model <- segmentation$model
cat("saving essential RData to file\n")
save(counts, model, regions, file = file.path(argv$o, "model.RData"))
