library(ehmm)
library(rtracklayer)
library(argparser)
library(tools)

args <- commandArgs(trailingOnly = T)
parser <- arg_parser("EHMM command line parser", name = "ehmm")
parser <- add_argument(parser, "-r", help = "BED file of input regions")
parser <- add_argument(parser, "-g", help = "Genome sizes")
parser <- add_argument(parser, "-m", help = "Construced model")
parser <- add_argument(parser, "-b", help = "Path to directory containing the bam files of the input data")
parser <- add_argument(parser, "-o", help = "Output directory", default = "./")

argv <- parse_args(parser, argv = args)
cat("loading input data")
regions <- import(argv$r, format = "bed")
model <- argv$m
bamfiles <- list.files(path = argv$m, pattern = "*.bam", full.names = T)
marks <- sapply(strsplit(file_path_sans_ext(basename(bamfiles)), "_"), "[", 1)
bamtab <- data.frame(mark = marks, path = bamfiles, check.names = F)

applyModel(regions, model = argv$m, genomeSize = argv$g, bamtab = bamtab, outdir = argv$o)
