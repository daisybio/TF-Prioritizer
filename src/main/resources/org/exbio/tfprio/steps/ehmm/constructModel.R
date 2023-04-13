library(ehmm)
library(rtracklayer)
library(argparser)

selectStates <- function(model){
  # assign marks
  mark.acc <- unlist(sapply(c('acc', 'atac', 'dhs', 'dnase'), function(pattern) which(grepl(pattern, tolower(model$marks)))))
  if (length(mark.acc) != 1){
    stop('Error: Exactly one of the given mark names must contain "atac", "dhs" or "dnase".\n
         If your chromatin accessibility assay is neither ATAC-seq nor DNase-seq, name it "ACC".\n
         Example: --mark ACC:/path/chromatin_accessibility_assay.bam')
  }
  mark.k27ac <- which(grepl('k27ac', tolower(model$marks)))
  mark.k4me1 <- which(grepl('k4me1', tolower(model$marks)))
  mark.k4me3 <- which(grepl('k4me3', tolower(model$marks)))
  if (any(sapply(c(mark.k27ac, mark.k4me1, mark.k4me3), function(mrk) length(mrk) != 1))){
    stop('Error: Automated state selection requires the three marks H3K27ac, H3K4me1, and H3K4me3 and their correct naming.\n
         If you work with different marks, please select the states manually.
         Example: --accStates.e 1,2 --nucStates.e 3,5 --accStates.p 2,3 --nucStates.p 1,4')
  }
  atac <- sapply(model$emisP, function(emis) emis$mu[mark.acc])
  k27ac <- sapply(model$emisP, function(emis) emis$mu[mark.k27ac])
  k4me1 <- sapply(model$emisP, function(emis) emis$mu[mark.k4me1])
  k4me3 <- sapply(model$emisP, function(emis) emis$mu[mark.k4me3])
  acc <- order(atac/k27ac, decreasing=T)[1:2]
  methyl.ratio <- order(k4me1/k4me3, decreasing=T)
  nuc.e <- methyl.ratio[!(methyl.ratio %in% acc)][1:2]
  nuc.p <- rev(methyl.ratio[!(methyl.ratio %in% acc)])[1:2]
  return(list(acc=acc, nuc.e=nuc.e, nuc.p=nuc.p))
}

args <- commandArgs(trailingOnly = T)
parser <- arg_parser("EHMM command line parser", name = "ehmm")
parser <- add_argument(parser, "-b", help = "Background model RData")
parser <- add_argument(parser, "-e", help = "Enhancer model RData")
parser <- add_argument(parser, "-p", help = "Promoter model RData")
parser <- add_argument(parser, "-t", help = "Number of threads", default = 1, type = "integer")
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

model <- constructModel(model.bg = bgModel, model.e = eModel, model.p = pModel,
               counts.e = eCounts, counts.p = pCounts,
               regions.e = eRegions, regions.p = pRegions,
               accStates.e = eStates$acc, accStates.p = pStates$acc,
               nucStates.e = eStates$nuc, nucStates.p = pStates$nuc,
               nthreads = nthreads, outdir = argv$o)
