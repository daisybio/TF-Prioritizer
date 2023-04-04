getConstructModelOptions <- function(){
  opts <- list(
    list(arg="--model.e", type="character", parser=readModel, required=TRUE,
         help="Path to the file with the parameters of the enhancer HMM."),
    list(arg="--model.p", type="character", parser=readModel, required=TRUE,
         help="Path to the file with the parameters of the promoter HMM."),
    list(arg="--counts.e", type="character", required=TRUE, parser=readCounts,
         help="Path to the countmatrix for the enhancer training regions."),
    list(arg="--counts.p", type="character", required=TRUE, parser=readCounts,
         help="Path to the countmatrix for the promoter training regions."),
    list(arg="--regions.e", type="character", required=TRUE, parser=readRegions,
         help="Path to the BED file with the enhancer training regions associated to the count matrix."),
    list(arg="--regions.p", type="character", required=TRUE, parser=readRegions,
         help="Path to the BED file with the promoter training regions associated to the count matrix."),
    list(arg="--model.bg", type="character", parser=readModel, required=FALSE,
         help="Path to the file with the parameters of the background HMM."),
    list(arg="--fg_to_bg", flag=TRUE,
         help="Flag for whether to create a background model from the foreground emission probs with unitized transition probs. Default: FALSE"),
    list(arg="--accStates.e", type="character", parser=readStates,
         help="String of comma-separated state-numbers for enhancer accessibility states. If not given, states selection will be automated."),
    list(arg="--nucStates.e", type="character", parser=readStates,
         help="String of comma-separated state-numbers for enhancer nucleosome states. If not given, states selection will be automated."),
    list(arg="--accStates.p", type="character", parser=readStates,
         help="String of comma-separated state-numbers for promoter accessibility states. If not given, states selection will be automated."),
    list(arg="--nucStates.p", type="character", parser=readStates,
         help="String of comma-separated state-numbers for promoter nucleosome states. If not given, states selection will be automated."),
    list(arg="--outdir", type="character",
         help="Path to the output directory."),
    list(arg="--nthreads", type="integer", default=formals(constructModel)$nthreads,
         help="Number of threads to be used")
  )
  opts
}

constructModelCLI <- function(args, prog){
  constructModelOptions <- getConstructModelOptions()
  #parse the options
  opt <- parseArgs(constructModelOptions, args, prog)
  #call 'constructModel'
  segmentation <- do.call(constructModel, opt)
}

#' Construct a total model by combining a background and two foreground models for enhancers and promoters.
#' Prior to combining, the foreground models are refined by relearning the transition parameters.
#'
#' @param model.e A list with the parameters that describe the enhancer HMM.
#' @param model.p A list with the parameters that describe the promoter HMM.
#' @param counts.e count matrix for enhancer training regions.
#' @param counts.p count matrix for promoter training regions.
#' @param regions.e GRanges object containing the enhancer training regions.
#' @param regions.p GRanges object containing the promoter training regions.
#' @param model.bg A list with the parameters that describe the background HMM.
#' @param fg_to_bg Flag for whether to create a background model from the foreground emission probs with unitized transition probs. Default: FALSE.
#' @param accStates.e String of comma-separated state-numbers for enhancer accessibility states. If not given, states selection will be automated.
#' @param nucStates.e String of comma-separated state-numbers for enhancer nucleosome states. If not given, states selection will be automated.
#' @param accStates.p String of comma-separated state-numbers for promoter accessibility states. If not given, states selection will be automated.
#' @param nucStates.p String of comma-separated state-numbers for promoter nucleosome states. If not given, states selection will be automated.
#' @param outdir path to the output directory
#' @param nthreads number of threads used for learning.
#' @return A list with the following arguments:
#' 
#' @export
constructModel <- function(model.e, model.p, counts.e, counts.p, regions.e, regions.p,
                           model.bg=NULL, fg_to_bg=FALSE, accStates.e=NULL, nucStates.e=NULL, accStates.p=NULL, nucStates.p=NULL, outdir=".", nthreads=1){
  # check arguments and define variables
  binsize <- 100
  if (any(is.null(accStates.e), is.null(nucStates.p), is.null(accStates.p), is.null(nucStates.p))){
    cat('Automated state selection\n')
    states.e <- selectStates(model.e)
    accStates.e <- states.e$acc
    nucStates.e <- states.e$nuc.e
    states.p <- selectStates(model.p)
    accStates.p <- states.p$acc
    nucStates.p <- states.p$nuc.p
  }
  if (all(is.null(model.bg), fg_to_bg==FALSE)) {
  	 cat('No background model specified. Creating FGtoBG background model using foreground emissions and unitized transitions.')
	 fg_to_bg <- TRUE
  }

  # construct model.bg from model.e and model.p if flag is set
  if (fg_to_bg == TRUE) {
  	 nstates.bg <- model.e$nstates+model.p$nstates
	 model.bg <- list(nstates=nstates.bg,
     		  	      marks=model.e$marks,
                 	  emisP=c(model.e$emisP, model.p$emisP),
                 	  transP=matrix(1/nstates.bg, nrow=nstates.bg, ncol=nstates.bg),
                 	  initP=rep(1/nstates.bg, nstates.bg),
                 	  labels=paste0('bg', 1:nstates.bg),
                 	  colors=tail(colorRampPalette(brewer.pal(9,'Greys'))(20), nstates.bg))
  }
  
  # construct initial enhancer / promoter models
  model.e.init <- initializeParams(model.e, accStates.e, nucStates.e)
  model.p.init <- initializeParams(model.p, accStates.p, nucStates.p)
  
  # refine enhancer / promoter models (relearn on training data with keeping emisP fixed)
  cat("Refine foreground models\n")
  segmentation.e.refinedTrans <- segment(counts=counts.e, regions=regions.e,  nstates=model.e.init$nstates, model=model.e.init, nthreads=nthreads,
                                         verbose_kfoots=TRUE, trainMode='viterbi', fix_emisP=TRUE, nbtype='lognormal', endstate=model.e.init$endstates)
  segmentation.p.refinedTrans <- segment(counts=counts.p, regions=regions.p, nstates=model.p.init$nstates, model=model.p.init, nthreads=nthreads,
                                         verbose_kfoots=TRUE, trainMode='viterbi', fix_emisP=TRUE, nbtype='lognormal', endstate=model.p.init$endstates)
  
  # add labels and colors to model objects
  segmentation.e.refinedTrans$model$labels <- c(paste0('E_N1.', 1:length(nucStates.e)), paste0('E_A.', 1:length(accStates.e)), paste0('E_N2.', 1:length(nucStates.e)))
  segmentation.p.refinedTrans$model$labels <- c(paste0('P_N1.', 1:length(nucStates.p)), paste0('P_A.', 1:length(accStates.p)), paste0('P_N2.', 1:length(nucStates.p)))
  model.bg$labels <- paste0('bg', 1:model.bg$nstates)
  segmentation.e.refinedTrans$model$colors <- c(rev(tail(brewer.pal(9,'Greens')[-c(8,9)], length(nucStates.e))),
                                                colorRampPalette(brewer.pal(9,'YlOrBr')[2:4])(length(accStates.e)),
                                                rev(tail(brewer.pal(9,'Greens')[-c(8,9)], length(nucStates.e))))
  segmentation.p.refinedTrans$model$colors <- c(rev(tail(brewer.pal(9,'Reds')[-9], length(nucStates.p))),
                                                colorRampPalette(brewer.pal(9,'YlOrBr')[2:4])(length(accStates.p)),
                                                rev(tail(brewer.pal(9,'Reds')[-9], length(nucStates.p))))
  model.bg$colors <- tail(colorRampPalette(brewer.pal(9,'Greys'))(20), model.bg$nstates)
  
  # combine fg / bg models, write to file
  model <- combineFgBgModels(model.bg, segmentation.e.refinedTrans$model, segmentation.p.refinedTrans$model)
  modelPath <- paste(outdir, 'model.txt', sep='/')
  writeModel(model, modelPath, type='lognormal')
}

readStates <- function(stateString){
  # This function parses a comma-separated string of state numbers to a integer-vector.
  return(as.integer(strsplit(stateString, ',')[[1]]))
}

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
