getDifferentialScoresOptions <- function(){
  opts <- list(
    list(arg="--scores1", type="character", required=TRUE, parser=readRegionsWithScores,
         help="Path to the first scores BED file."),
    list(arg="--scores2", type="character", required=TRUE, parser=readRegionsWithScores,
         help="Path to the second scores BED file."),
    list(arg="--outfile", type="character",
         help="Path to the output file.")
  )
  opts
}

differentialScoresCLI <- function(args, prog){
  differentialScoresOptions <- getDifferentialScoresOptions()
  #parse the options
  opt <- parseArgs(differentialScoresOptions, args, prog)
  #call 'differentialScores'
  segmentation <- do.call(differentialScores, opt)
}

#' Produce a segmentation based on a given model and extract enhancer / promoter elemenents.
#'
#' @param scores1 GRanges object containing scores.
#' @param scores2 GRanges object containing scores.
#' @param outfile path to the output file.
#' @return nothing.
#' 
#' @export
differentialScores <- function(scores1, scores2, outfile="differentialScores.bed"){
  # check arguments and define variables
  if (!all(scores1 == scores2)) stop('score files need to contain exactly the same regions.')

  cat("calculate differential scores\n")
  diffScores <- scores1
  diffScores$scores <- scores1$score - scores2$score
  export.bed(diffScores, outfile)
}
