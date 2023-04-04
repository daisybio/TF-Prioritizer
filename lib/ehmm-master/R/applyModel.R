getApplyModelOptions <- function(){
  opts <- list(
    list(arg="--regions", type="character", required=TRUE, parser=readRegions,
         help="Path to the BED file with the genomic regions of interest."),
    list(arg="--genomeSize", type="character", required=TRUE, parser=readGenomeSize,
         help="Path to a two-column file indicating chromosome names and sizes."),
    list(arg="--model", type="character", parser=readModel,
         help="Path to the file with the parameters of the HMM. Only required if --provideModel flag is not set."),
    list(arg="--provideModel", flag=TRUE,
         help="Whether or not to use the provided model that was learned on mouse embryonic stem cell data.
         If this flag is set, query data will be normalized to the data that was used during model training."),
    list(arg="--counts", type="character", parser=readCounts,
         help="Path to the count matrix. If not given, it will be calculated and written to the output directory."),
    list(arg="--mark", type="character", vectorial=TRUE, meta="label:path",
         help="Mark name and path to a bam file where to extract reads from. Only required if counts are not given.
         The bam files must be indexed and the chromosome names must match with
         those specified in the bed file. Entries with the same mark name will
         be treated as replicates and collapsed into one experiment.
         This option must be repeated for each mark, for example:
         `--mark H3K4me3:/path1/foo1.bam --mark H3K36me3:/path2/foo2.bam`"),
    list(arg="--outdir", type="character",
         help="Path to the output directory."),
    list(arg="--nthreads", type="integer", default=formals(applyModel)$nthreads,
         help="Number of threads to be used."),
    list(arg="--learnTrans", flag=TRUE,
         help="Whether or not to let the model learn the transition probabilities while fixing the emission parameters."),
    list(arg="--refCounts", type="character", parser=readCounts,
         help="Path to the count matrix of the reference model. If given, the counts of the query sample will be quantile-normalized to its distribution.
         refCounts should ideally represent a full genome and not just a subset.")
  )
  opts
}

applyModelCLI <- function(args, prog){
  applyModelOptions <- getApplyModelOptions()
  #parse the options
  opt <- parseArgs(applyModelOptions, args, prog)
  #make bamtab object (only if opt$counts is missing)
  if ('mark' %in% names(opt)){
    opt$bamtab <- makeBamtab(opt$mark)
    opt <- opt[names(opt) != 'mark'] # remove 'mark' elements from opt
  } else {
    if (!('counts' %in% names(opt))) stop('either pass a count matrix (--counts) or specify bam-files to calculate it from (--mark)')
  }
  #call 'applyModel'
  segmentation <- do.call(applyModel, opt)
}

#' Produce a segmentation based on a given model and extract enhancer / promoter elemenents.
#'
#' @param regions GRanges object containing the genomic regions of interest.
#' @param model A list with the parameters that describe the HMM.
#' @param bamtab Data frame describing how to get the counts from each file. 
#'     The following columns are required:
#'     'mark' (name of the histone mark),
#'     'path' (path to the bam file).
#' @param provideModel flag, whether or not to use the provided model that was learned on mouse embryonic stem cell data.
#' @param genomeSize vector with chromosome lengths.
#' @param counts Count matrix matching with the \code{regions} parameter.
#' Each row of the matrix represents a mark and each column a bin resulting
#' from dividing the genomic regions into non-overlapping bins of equal size.
#' The rows of the matrix must be named with the name of the marks and these names must be unique.
#' If not given, it will be calculated and written to the output directory.
#' @param outdir path to the output directory.
#' @param nthreads number of threads used for learning.
#' @param learnTrans flag, whether or not to let the model learn the transition probabilities while fixing the emission parameters.
#' @param refCounts Count matrix of the reference model.
#' Has to be given if \code{refCounts} and \code{counts} have different dimensions.
#' @return nothing.
#' 
#' @export
applyModel <- function(regions, model=NULL, provideModel=FALSE, genomeSize, counts=NULL, bamtab=NULL,
                       outdir=".", nthreads=1, learnTrans=FALSE, refCounts=NULL){
  # check arguments and define variables
  binsize <- 100
  if (!is.null(model)) provideModel <- FALSE
  if (is.null(model) && !provideModel){
    cat('No model specified. Provided mESC model will be used.')
    provideModel <- TRUE
  }
  # create outdir if it does not exist
  dir.create(file.path(outdir), showWarnings=FALSE)

  # remove unnamed, random and mitochondrial chromosomes from regions' seqlevels
  levelsToDrop <- unique(unlist(lapply(c('Un', 'M', 'random', 'hap', 'alt', 'GL', 'NC', 'hs'), function(x) which(grepl(x, GenomeInfoDb:::seqlevels(regions))))))
  if (length(levelsToDrop) > 0) GenomeInfoDb:::seqlevels(regions) <- GenomeInfoDb:::seqlevels(regions)[-levelsToDrop]
  
  # deal with the provideModel flag
  if (provideModel){
    # load provided rdata file with model and count table, which contains the counts as names and their numbers of occurrences as values.
    # reconstruct a refCounts matrix from the count table.
    # check the passed --mark arguments. They must fit to the mESC model, i.e. one accessibility containing 'atac', 'dhs', 'dnase' or 'acc',
    # as well as the three basic HM's H3K27ac and H3K4me1/3.
    mark.acc <- unlist(sapply(c('acc', 'atac', 'dhs', 'dnase'), function(pattern) which(grepl(pattern, tolower(bamtab$mark)))))
    if (length(mark.acc) != 1){
      stop('Error: Exactly one of the given mark names must contain "atac", "dhs" or "dnase".\n
         If your chromatin accessibility assay is neither ATAC-seq nor DNase-seq, name it "ACC".\n
         Example: --mark ACC:/path/chromatin_accessibility_assay.bam')
    }
    mark.k27ac <- which(grepl('k27ac', tolower(bamtab$mark)))
    mark.k4me1 <- which(grepl('k4me1', tolower(bamtab$mark)))
    mark.k4me3 <- which(grepl('k4me3', tolower(bamtab$mark)))
    
    load(system.file("extdata", "mESC.rdata", package="ehmm", mustWork=TRUE))
    model$marks <- bamtab$mark[c(mark.acc, mark.k27ac, mark.k4me1, mark.k4me3)] # rename provided model marks to the ones passed by the user.
    names(counts.tables) <- model$marks
    refCounts <- reconstructCountmatrix(counts.tables)
  }
  
  # if not given, calculate and save count matrix
  if (is.null(counts)){
    counts <- getCountMatrix(regions, bamtab, outdir, binsize=100, nthreads, pseudoCount=1)
  }
  
  # if reference count matrix is given, quantile normalize query count matrix
  if (!(is.null(refCounts))){
    # if regions seqnames don't match to the genomeSize names, the regions are probably in the annoying NCBI style ("1" instead of "chr1")
    # if so, change genomeSize format to NCBI.
    if (!all(as.vector(unique(seqnames(regions))) %in% names(genomeSize))){ 
      names(genomeSize) <- gsub('chr', '', names(genomeSize))
    }
    counts <- quantileNormalizeCounts(counts=counts, refCounts=refCounts, regions=regions, genomeSize=genomeSize,
                                      bamtab=bamtab, outdir=outdir, nthreads=nthreads)
  }

  # segment regions
  cat("apply model\n")
  segmentation <- segment(counts=counts, regions=regions, model=model, nstates=model$nstates, nthreads=nthreads,
                          verbose_kfoots=TRUE, nbtype='lognormal', notrain=!learnTrans, fix_emisP=learnTrans)
  
  # produce 'report'
  cat("producing report\n")
  viterbi_segments <- statesToSegments(segmentation$viterbi, segmentation$segments) # create GRanges object with viterbi states
  report(segments=viterbi_segments, model=model, rdata=segmentation, outdir=outdir, colors=model$colors, labels=model$labels)
  
  # extract enhancers / promoters
  cat("extract enhancer / promoter elements\n")
  extractRegions(segmentation, regions, genomeSize, outdir)
  
  cat("done\n\n")
}

extractRegions <- function(segmentation, regions, genomeSize, outdir){
  # this function tiles regions into 100 bp windows, assigns viterbi states and scores, extracts enhancer / promoter regions according to viterbi decoding,
  # and writes both scores and decoded elements to file
  labels <- segmentation$model$labels
  gr <- unlist(tile(regions, width=100))
  GenomeInfoDb:::seqlengths(gr) <- genomeSize[GenomeInfoDb:::seqlevels(gr)]
  gr$name <- labels[segmentation$viterbi]
  mapply(function(scores, elmName){
    label <- toupper(substr(elmName, 1, 1))
    gr$score <- scores
    export.bw(gr, sprintf('%s/%s.scores.bw', outdir, elmName))
    # export.bed(gr, sprintf('%s/%s.scores.bed', outdir, elmName))
    regionsBedfile <- sprintf('%s/%sRegions.bed', outdir, elmName)
    elms.tiled <- gr[startsWith(gr$name, paste0(label, '_A'))] # only save accessibility states
    file.create(regionsBedfile)
    if (length(elms.tiled) > 0){
      elms <- aggScore(reduce(elms.tiled), elms.tiled, 'max')
      export.bed(elms, regionsBedfile)
    }
  }, segmentation$score, c('enhancer', 'promoter'))
}

readGenomeSize <- function(genomeSize){
  # this function parses the genomeSize file that must contain one column of chromosome names and one column of chromosome sizes
  # unnamed, random and mitochondrial chromosomes are ignored
  df <- read.table(genomeSize)
  sizes <- df$V2
  names(sizes) <- df$V1
  levelsToDrop <- unique(unlist(lapply(c('Un', 'M', 'random', 'hap', 'alt'), function(x) which(grepl(x, names(sizes))))))
  if (length(levelsToDrop) > 1) sizes <- sizes[-levelsToDrop]
  return(sizes)
}

aggScore <- function(gr.reference, gr.tiled, func, aggName=F){
  # this function aggregates the scores of a tiled GRanges object (gr.tiled) to the overlapping regions of a GRanges object with broad regions (gr.reference)
  # by either taking the mean, max or the product given the passed function argument.
  # Often, gr.reference is equal to reduce(gr.tiled).
  ov <- findOverlaps(gr.reference, gr.tiled)
  aggSc <- aggregate(gr.tiled$score[subjectHits(ov)], list(queryHits(ov)), get(func))
  gr.reference$score <- 0
  gr.reference$score[unique(queryHits(ov))] <- aggSc$x
  if (aggName && !is.null(gr.tiled$name)) gr.reference$name <- aggregate(gr.tiled$name[subjectHits(ov)], list(queryHits(ov)), 'unique')$x
  return(gr.reference)
}

clipCounts <- function(cm, percentile){
  cm.clipped <- cm
  for (i in 1:nrow(cm)) {
    clipValue <- quantile(cm[i,], percentile)
    cm.clipped[i, (cm.clipped[i,] > clipValue)] <- clipValue
  }
  return(cm.clipped)
}

reconstructCountmatrix <- function(counts.tables){
  # This function reconstructs a sorted countmatrix based on passed contingency tables (counts with frequencies)
  # I use nested for loops over the preinitiated countmatrix instead of a one-liner sapply function for the sake of speed.
  refCounts <- matrix(nrow=length(counts.tables), ncol=sum(counts.tables[[1]]))
  row.names(refCounts) <- names(counts.tables)
  for (feature in names(counts.tables)) {
    tbl <- counts.tables[[feature]]
    end <- cumsum(tbl)
    start <- end-tbl+1
    for (count in names(tbl)) {
      freq <- tbl[count]
      refCounts[feature, start[count]:end[count]] <- as.integer(rep(count, freq))
    }
  }
  return(refCounts)
}
