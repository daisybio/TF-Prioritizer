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
  acc <- sapply(model$emisP, function(emis) emis$mu[mark.acc])
  hasK27ac <- grepl("k27ac", model$marks, ignore.case = T)
  if (any(hasK27ac)) {
    mark.k27ac <- which(hasK27ac)
    k27ac <- sapply(model$emisP, function(emis) emis$mu[mark.k27ac])
    acc <- acc/k27ac
  }
  acc <- order(acc, decreasing=T)[1:2]
  hasK3me1 <- grepl("k4me1", model$marks, ignore.case = T)
  hasK3me3 <- grepl("k4me3", model$marks, ignore.case = T)
  if (sum(hasK3me1) + sum(hasK3me3) == 0) {
    stop("Error: marks have to contain at least either H3K4me1 or H3K4me3")
  }
  if (any(hasK3me1) & any(hasK3me3)) {
    mark.k4me1 <- which(hasK3me1)
    mark.k4me3 <- which(hasK3me3)
    k4me1 <- sapply(model$emisP, function(emis) emis$mu[mark.k4me1])
    k4me3 <- sapply(model$emisP, function(emis) emis$mu[mark.k4me3])
    ratio <- k4me1/k4me3
  } else if (any(hasK3me1)) {
    mark.k4me1 <- which(hasK3me1)
    ratio <- sapply(model$emisP, function(emis) emis$mu[mark.k4me1])
  } else {
    mark.k4me3 <- which(hasK3me3)
    ratio <- sapply(model$emisP, function(emis) emis$mu[mark.k4me3])
  }
  methyl.ratio <- order(ratio, decreasing=T)
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
# eStates <- selectStates(eModel)
load(argv$p)
pModel <- model
pCounts <- counts
pRegions <- regions
# pStates <- selectStates(pModel)
cat("constructing combined model...\n")
nthreads = argv$t

model <- constructModel(model.bg = bgModel, model.e = eModel, model.p = pModel,
               counts.e = eCounts, counts.p = pCounts,
               regions.e = eRegions, regions.p = pRegions,
               nthreads = nthreads, outdir = argv$o)
