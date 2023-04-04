getLearnModelOptions <- function(){
  opts <- list(
    list(arg="--regions", type="character", required=TRUE, parser=readRegions,
         help="Path to the BED file with the genomic regions of interest.
         These regions will be automatically partitioned into smaller, 
         consecutive bins. Only the first three fields in the file matter. 
         If the region lengths are not multiples of the given binsize
         a new bed file will be produced where each coordinate 
         is a multiple of binsize. Use this new file together with
         the count matrix for later analyses."),
    list(arg="--nstates", type="integer", required=TRUE,
         help="Number of states to use for training"),
    list(arg="--counts", type="character", parser=readCounts,
         help="Path to the count matrix. If not given, it will be calculated and written to the output directory."),
    list(arg="--mark", type="character", vectorial=TRUE, meta="label:path",
         help="Mark name and path to a bam file where to extract reads from. Only required if counts are not given.
         The bam files must be indexed and the chromosome names must match with
         those specified in the bed file. Entries with the same mark name will
         be treated as replicates and collapsed into one experiment.
         This option must be repeated for each mark, for example:
         `-m H3K4me3:/path1/foo1.bam -m H3K36me3:/path2/foo2.bam`"),
    list(arg="--outdir", type="character",
         help="Path to the output directory."),
    list(arg="--nthreads", type="integer", default=formals(learnModel)$nthreads,
         help="Number of threads to be used"),
    list(arg="--pseudoCount", type="integer",
         help="Pseudo-count added to read-counts for log-normal fit. Default = 1.")
    )
  opts
}

learnModelCLI <- function(args, prog){
  learnModelOptions <- getLearnModelOptions()
  #parse the options
  opt <- parseArgs(learnModelOptions, args, prog)
  ## make bamtab object
  # opt$bamtab <- makeBamtab(opt$mark)
  # opt <- opt[names(opt) != 'mark'] # remove 'mark' elements from opt
  #make bamtab object (only if opt$counts is missing)
  if ('mark' %in% names(opt)){
    opt$bamtab <- makeBamtab(opt$mark)
    opt <- opt[names(opt) != 'mark'] # remove 'mark' elements from opt
  } else {
    if (!('counts' %in% names(opt))) stop('either pass a count matrix (--counts) or specify bam-files to calculate it from (--mark)')
  }
  #call 'learnModel'
  segmentation <- do.call(learnModel, opt)
}

#' Learn a model and produce a segmentation
#'
#' @param regions GRanges object containing the genomic regions of interest.
#' @param nstates Number of states to learn.
#' @param bamtab Data frame describing how to get the counts from each file. 
#'     The following columns are required:
#'     'mark' (name of the histone mark),
#'     'path' (path to the bam file).
#' @param counts Count matrix matching with the \code{regions} parameter.
#' Each row of the matrix represents a mark and each column a bin resulting
#' from dividing the genomic regions into non-overlapping bins of equal size.
#' The rows of the matrix must be named with the name of the marks and these names must be unique.
#' If not given, it will be calculated and written to the output directory.

#' @param outdir path to the output directory.
#' @param nthreads number of threads used for learning.
#' @param pseudoCount pseudo-count to add to read counts.
#' This is necessary because log-counts are calculated in order to fit a log-normal distribution.
#' @return A list with the following arguments:
#' 
#' @export
learnModel <- function(regions, nstates, bamtab, counts=NULL, outdir=".", nthreads=1, pseudoCount=1){
  # if not given, calculate and save count matrix
  if (is.null(counts)){
    counts <- getCountMatrix(regions, bamtab, outdir, binsize=100, nthreads, pseudoCount)
  }

  # learn unsupervised model
  #TODO: implement fast learning on random 20mb subset..?
  cat("learning model\n")
  segmentation <- segment(counts=counts, regions=regions, nstates=nstates, nthreads=nthreads,
                   verbose_kfoots=TRUE, nbtype='lognormal')
  
  # produce 'report'
  cat("producing report\n")
  viterbi_segments <- statesToSegments(segmentation$viterbi, segmentation$segments) # create GRanges object with viterbi states
  report(segments=viterbi_segments, model=segmentation$model, rdata=segmentation, outdir=outdir)

  return(segmentation)
}

