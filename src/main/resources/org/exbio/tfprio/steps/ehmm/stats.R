library(argparser)
library(rtracklayer)
library(GenomicRanges)


stats <- function(predictedRanges, referenceRanges, allRanges){
  N <- length(allRanges)
  predPos <- subsetByOverlaps(allRanges, predictedRanges)
  predNeg <- subsetByOverlaps(allRanges, predictedRanges, invert = T)

  refPos <- subsetByOverlaps(allRanges, referenceRanges)
  refNeg <- subsetByOverlaps(allRanges, referenceRanges, invert = T)

  tp <- length(findOverlaps(predPos, refPos))
  fp <- length(findOverlaps(predPos, refNeg))
  tn <- length(findOverlaps(predNeg, refNeg))
  fn <- length(findOverlaps(predNeg, refPos))
  acc <- (tp+tn)/N
  mcr <- 1-acc
  sens <- tp/length(refPos)
  spec <- tn/length(refNeg)
  prec <- tp/(tp+fp)
  prev <- length(refPos)/N
  return(c(TP=tp,FP=fp,TN=tn,FN=fn,N=N,Accuracy=acc,Missclassification=mcr,
           Sensitivity=sens,Specificity=spec,Precision=prec,Prevalence=prev))
}


args <- commandArgs(trailingOnly = T)
parser <- arg_parser("ehmm statistics command line parser", name = "ehmm_stats")
parser <- add_argument(parser, "--p_enhancers", help = "predicted enhancers bed file")
parser <- add_argument(parser, "--r_enhancers", help = "reference enhancers bed file")
parser <- add_argument(parser, "--p_promoters", help = "predicted promoters bed file")
parser <- add_argument(parser, "--r_promoters", help = "reference promoters bed file")
parser <- add_argument(parser, "-i", help = "input collapsed peak bed dir")
parser <- add_argument(parser, "-o", help = "output file", default = "./")
argv <- parse_args(parser, argv = args)

predicted_enhancers <- reduce(import(argv$p_enhancers, format = "bed"))
reference_enhancers <- reduce(import(argv$r_enhancers, format = "bed"))

predicted_promoters <- reduce(import(argv$p_promoters, format = "bed"))
reference_promoters <- reduce(import(argv$r_promoters, format = "bed"))

collapsedPeakFiles <- list.files(argv$i, pattern = ".bed", recursive = T, full.names = T)
peaks <- reduce(sort(do.call(c, lapply(collapsedPeakFiles, import, format = "bed"))))

enhancers_stats <- t(data.frame(round(stats(predicted_enhancers, reference_enhancers, peaks), 4)))
promoters_stats <- t(data.frame(round(stats(predicted_promoters, reference_promoters, peaks), 4)))
d <- data.frame(rbind(enhancers_stats, promoters_stats))
d <- data.frame(type=c("enhancers", "promoters"), d)

write.table(d, file.path(argv$o, "Statistics.tsv"), row.names = F, quote = F)
