library(argparser)
library(rtracklayer)
library(GenomicRanges)
library(caret)

stats <- function(predictedRanges, referenceRanges, allRanges){
  pred <- countOverlaps(allRanges, predictedRanges)
  ref <- countOverlaps(allRanges, referenceRanges)
  pred[pred > 0] <- 1
  ref[ref > 0] <- 1
  return(confusionMatrix(factor(pred), factor(ref)))
}


args <- commandArgs(trailingOnly = T)
parser <- arg_parser("ehmm statistics command line parser", name = "ehmm_stats")
parser <- add_argument(parser, "--p_enhancers", help = "predicted enhancers bed file")
parser <- add_argument(parser, "--r_enhancers", help = "reference enhancers bed file")
parser <- add_argument(parser, "--p_promoters", help = "predicted promoters bed file")
parser <- add_argument(parser, "--r_promoters", help = "reference promoters bed file")
parser <- add_argument(parser, "-i", help = "input collapsed peak bed dir")
parser <- add_argument(parser, "-o", help = "Output file")
argv <- parse_args(parser, argv = args)

predicted_enhancers <- reduce(import(argv$p_enhancers, format = "bed"))
reference_enhancers <- reduce(import(argv$r_enhancers, format = "bed"))

predicted_promoters <- reduce(import(argv$p_promoters, format = "bed"))
reference_promoters <- reduce(import(argv$r_promoters, format = "bed"))

collapsedPeakFiles <- list.files(argv$i, pattern = ".bed", recursive = T, full.names = T)
peaks <- reduce(sort(do.call(c, lapply(collapsedPeakFiles, import, format = "bed"))))

s_e <- stats(predicted_enhancers, reference_enhancers, peaks)
s_p <- stats(predicted_promoters, reference_promoters, peaks)
d <- data.frame(enhancers=s_e$byClass, promoters=s_p$byClass)

write.table(d, argv$o, row.names = F, quote = F)
