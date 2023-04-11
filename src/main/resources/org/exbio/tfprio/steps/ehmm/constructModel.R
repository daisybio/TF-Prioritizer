library(ehmm)
library(rtracklayer)
library(argparser)
library(RColorBrewer)
library(Matrix)

rowsToStr <- function (mat) {
  apply(mat, 1, paste, collapse = "\t")
}

selectStates <- function(model, n = 2){
  d <- unlist(sapply(model$emisP, "[", "mus"))
  acc <- order(d, decreasing = T)[1:n]
  nuc <- order(d)[1:n]
  return(list(acc=acc, nuc=nuc))
}

initializeParams <- function (model, states.a, states.n) {
  n.acc <- length(states.a)
  n.nuc <- length(states.n)
  model$nstates <- n.acc + 2 * n.nuc
  newStates.n1 <- 1:n.nuc
  newStates.a <- 1:n.acc + max(newStates.n1)
  newStates.n2 <- newStates.n1 + max(newStates.a)
  model$emisP <- model$emisP[c(states.n, states.a, states.n)]
  model$transP <- matrix(0, nrow = model$nstates, ncol = model$nstates)
  model$transP[newStates.n1, newStates.n1] <- 0.5/n.nuc
  model$transP[newStates.n1, newStates.a] <- 0.5/n.acc
  model$transP[newStates.a, newStates.a] <- 0.9/n.acc
  model$transP[newStates.a, newStates.n2] <- 0.1/n.nuc
  model$transP[newStates.n2, newStates.n2] <- 1/n.nuc
  model$initP <- matrix(c(rep(1/n.nuc, n.nuc), rep(0, n.acc +
                                                     n.nuc)))
  model$endstates <- newStates.n2
  return(model)
}

validateModel <- function (model, strict = FALSE, type = "dep") {
  if (is.null(model))
    stop("NULL is not a valid model")
  optArgs <- c("nstates", "marks", "labels", "colors", "emisP",
               "transP", "initP")
  if (is.null(names(model)))
    stop("model's element must be named")
  if (!strict && !all(optArgs %in% names(model)))
    stop(paste("a model requires the following fields:",
               paste(optArgs, collapse = ", ")))
  nstates <- model$nstates
  nmarks <- NULL

  if (!is.null(model$emisP)) {
    if (!is.list(model$emisP))
      stop("'model$emissP' must be a list'")
    if (!is.null(nstates)) {
      if (length(model$emisP) != nstates)
        stop("incorrect number of emission probabilities provided")
    }
    else nstates <- length(model$emisP)
    for (nm in model$emisP) {
      if (type == "lognormal") {
        needPars <- c("mus", "sigmasqs")
        if (!all(needPars %in% names(nm)))
          stop("missing fields from the emission probabilities")
        if (any(is.na(unlist(nm))))
          stop("NAs/NaNs in the emission probabilities not allowed")
        if (is.null(nmarks))
          nmarks <- length(nm$mus)
        if (length(nm$mus) != nmarks) {
          print(nm)
          stop("invalid parameters for the emission probabilities")
        }
      }
      else {
        needPars <- c("mu", "r", "ps")
        if (!all(needPars %in% names(nm)))
          stop("missing fields from the emission probabilities")
        if (any(is.na(unlist(nm))))
          stop("NAs/NaNs in the emission probabilities not allowed")
        if (is.null(nmarks))
          nmarks <- length(nm$ps)
        if (length(nm$ps) != nmarks)
          stop("invalid parameters for the emission probabilities")
        if (!isProbVector(nm$ps))
          stop("'model$emissP[[i]]$ps' must sum up to 1 for every i")
      }
    }
  }
  if (!is.null(model$transP)) {
    if (!is.matrix(model$transP) || ncol(model$transP) !=
        nrow(model$transP) || (!is.null(nstates) && nstates !=
                               ncol(model$transP)))
      stop("invalid trasition probabilities")
    if (!all(apply(model$transP, 1, isProbVector)))
      stop("rows of model$transP must sum up to 1")
    if (is.null(nstates))
      nstates <- nrow(model$transP)
  }
  if (!is.null(model$initP)) {
    if (!is.matrix(model$initP) || (!is.null(nstates) &&
                                    nrow(model$initP) != nstates))
      stop("invalid initial probabilities")
    if (!all(apply(model$initP, 2, isProbVector)))
      stop("columns of model$initP must sum up to 1")
    if (is.null(nstates))
      nstates <- nrow(model$initP)
  }
  list(nstates = nstates, nmarks = nmarks)
}

writeModel <- function (model, path, type) {
  if (all(names(model$emisP[[1]]) %in% c("mus", "sigmasqs")))
    type <- "lognormal"
  else type <- "negmultinom"
  # validateModel(model, strict = TRUE, type = type)
  txt <- list()
  if (!is.null(model$nstates))
    txt$nstates <- c("nstates", model$nstates)
  if (!is.null(model$labels))
    txt$labels <- c("labels", paste0(model$labels, collapse = ","))
  if (!is.null(model$colors))
    txt$colors <- c("colors", paste0(model$colors, collapse = ","))
  if (!is.null(model$marks))
    txt$marks <- c("marks", paste0(model$marks, collapse = "\t"))
  if (!is.null(model$emisP)) {
    if (type == "lognormal") {
      emisMat <- t(do.call(cbind, lapply(model$emisP,
                                         function(m) paste(m$mus, m$sigmasqs, sep = "|"))))
    }
    else {
      nbs <- getNBs(model$emisP)
      rs <- nbs[2, ]
      meansMat <- getMeanMatrix(model$emisP)
      emisMat <- t(rbind(rs, meansMat))
    }
    txt$emisP <- c("emisP", rowsToStr(emisMat))
  }
  if (!is.null(model$transP))
    txt$transP <- c("transP", rowsToStr(model$transP))
  if (!is.null(model$initP))
    txt$initP <- c("initP", rowsToStr(model$initP))
  writeLines(unlist(txt), path)
}

combineFgBgModels <- function (model.bg, model.e, model.p) {
  labels <- c(model.e$labels, model.p$labels, model.bg$labels)
  colors <- c(model.e$colors, model.p$colors, model.bg$colors)
  nstates <- sum(model.e$nstates, model.p$nstates, model.bg$nstates)
  marks <- unique(c(model.bg$marks, model.e$marks, model.p$marks))
  E <- c(model.e$emisP, model.p$emisP, model.bg$emisP)
  nEnhancers <- 399124
  nPromoters <- 70292
  nBins <- 3e+07
  enhancerFreq <- nEnhancers/nBins
  promoterFreq <- nPromoters/nBins
  states.n1.e <- which(startsWith(labels, "E_N1"))
  states.a.e <- which(startsWith(labels, "E_A"))
  states.n2.e <- which(startsWith(labels, "E_N2"))
  states.n1.p <- which(startsWith(labels, "P_N1"))
  states.a.p <- which(startsWith(labels, "P_A"))
  states.n2.p <- which(startsWith(labels, "P_N2"))
  states.bg <- which(startsWith(labels, "bg"))
  A <- as.matrix(bdiag(model.e$transP, model.p$transP, model.bg$transP))
  A[states.bg, states.n1.e] <- enhancerFreq/(length(states.bg) *
                                               length(states.n1.e))
  A[states.bg, states.n1.p] <- promoterFreq/(length(states.bg) *
                                               length(states.n1.p))
  A[states.bg, states.bg] <- A[states.bg, states.bg] * (1 -
                                                          rowSums(A[states.bg, c(states.n1.e, states.n1.p)]))
  A[states.n2.e, states.bg] <- rowSums(A[states.n1.e, states.a.e])/length(states.bg)
  A[states.n2.p, states.bg] <- rowSums(A[states.n1.p, states.a.p])/length(states.bg)
  A[states.n2.e, states.n2.e] <- A[states.n2.e, states.n2.e] *
    (1 - rowSums(A[states.n2.e, states.bg]))
  A[states.n2.p, states.n2.p] <- A[states.n2.p, states.n2.p] *
    (1 - rowSums(A[states.n2.p, states.bg]))
  I <- matrix(rep(0, nstates))
  I[c(states.n1.e, states.n1.p, states.bg)] <- 1/length(c(states.n1.e,
                                                          states.n1.p, states.bg))
  model <- list(nstates = nstates, marks = marks, emisP = E,
                transP = A, initP = I, labels = labels, colors = colors)
  return(model)
}

constructModel <- function (model.e, model.p, counts.e, counts.p, regions.e, regions.p,
          model.bg = NULL, fg_to_bg = FALSE, accStates.e = NULL, nucStates.e = NULL,
          accStates.p = NULL, nucStates.p = NULL, outdir = ".", nthreads = 1) {
  binsize <- 100
  if (any(is.null(accStates.e), is.null(nucStates.p), is.null(accStates.p),
          is.null(nucStates.p))) {
    cat("Automated state selection\n")
    states.e <- selectStates(model.e)
    accStates.e <- states.e$acc
    nucStates.e <- states.e$nuc.e
    states.p <- selectStates(model.p)
    accStates.p <- states.p$acc
    nucStates.p <- states.p$nuc.p
  }
  if (all(is.null(model.bg), fg_to_bg == FALSE)) {
    cat("No background model specified. Creating FGtoBG background model using foreground emissions and unitized transitions.")
    fg_to_bg <- TRUE
  }
  if (fg_to_bg == TRUE) {
    nstates.bg <- model.e$nstates + model.p$nstates
    model.bg <- list(nstates = nstates.bg, marks = model.e$marks,
                     emisP = c(model.e$emisP, model.p$emisP), transP = matrix(1/nstates.bg,
                                                                              nrow = nstates.bg, ncol = nstates.bg), initP = rep(1/nstates.bg,
                                                                                                                                 nstates.bg), labels = paste0("bg", 1:nstates.bg),
                     colors = tail(colorRampPalette(brewer.pal(9, "Greys"))(20),
                                   nstates.bg))
  }
  model.e.init <- initializeParams(model.e, accStates.e, nucStates.e)
  model.p.init <- initializeParams(model.p, accStates.p, nucStates.p)
  cat("Refine foreground models\n")
  segmentation.e.refinedTrans <- segment(counts = counts.e,
                                         regions = regions.e, nstates = model.e.init$nstates,
                                         model = model.e.init, nthreads = nthreads, verbose_kfoots = TRUE,
                                         trainMode = "viterbi", fix_emisP = TRUE, nbtype = "lognormal",
                                         endstate = model.e.init$endstates)
  segmentation.p.refinedTrans <- segment(counts = counts.p,
                                         regions = regions.p, nstates = model.p.init$nstates,
                                         model = model.p.init, nthreads = nthreads, verbose_kfoots = TRUE,
                                         trainMode = "viterbi", fix_emisP = TRUE, nbtype = "lognormal",
                                         endstate = model.p.init$endstates)
  segmentation.e.refinedTrans$model$labels <- c(paste0("E_N1.",
                                                       1:length(nucStates.e)), paste0("E_A.", 1:length(accStates.e)),
                                                paste0("E_N2.", 1:length(nucStates.e)))
  segmentation.p.refinedTrans$model$labels <- c(paste0("P_N1.",
                                                       1:length(nucStates.p)), paste0("P_A.", 1:length(accStates.p)),
                                                paste0("P_N2.", 1:length(nucStates.p)))
  model.bg$labels <- paste0("bg", 1:model.bg$nstates)
  segmentation.e.refinedTrans$model$colors <- c(rev(tail(brewer.pal(9,
                                                                    "Greens")[-c(8, 9)], length(nucStates.e))), colorRampPalette(brewer.pal(9,
                                                                                                                                            "YlOrBr")[2:4])(length(accStates.e)), rev(tail(brewer.pal(9,
                                                                                                                                                                                                      "Greens")[-c(8, 9)], length(nucStates.e))))
  segmentation.p.refinedTrans$model$colors <- c(rev(tail(brewer.pal(9,
                                                                    "Reds")[-9], length(nucStates.p))), colorRampPalette(brewer.pal(9,
                                                                                                                                    "YlOrBr")[2:4])(length(accStates.p)), rev(tail(brewer.pal(9,
                                                                                                                                                                                              "Reds")[-9], length(nucStates.p))))
  model.bg$colors <- tail(colorRampPalette(brewer.pal(9, "Greys"))(20),
                          model.bg$nstates)
  model <- combineFgBgModels(model.bg, segmentation.e.refinedTrans$model,
                             segmentation.p.refinedTrans$model)
  modelPath <- paste(outdir, "model.txt", sep = "/")
  writeModel(model, modelPath, type = "lognormal")
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
