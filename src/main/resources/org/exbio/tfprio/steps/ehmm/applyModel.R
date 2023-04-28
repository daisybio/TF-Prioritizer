library(ehmm)
library(rtracklayer)
library(argparser)
library(tools)
library(parallel)
library(dplyr)
library(bamsignals)

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
  if (!is.null(model$marks)) {
    if (anyDuplicated(model$marks))
      stop("model$marks cannot contain duplicate entries")
    nmarks <- length(model$marks)
  }
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
        if (length(nm$mus) != nmarks)
          stop("invalid parameters for the emission probabilities")
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

parseRows <- function (rows) {
  rowList <- strsplit(rows, "\t")
  rowLens <- sapply(rowList, length)
  if (any(rowLens != rowLens[1]))
    stop("matrix rows don't have the same length")
  rowList <- lapply(rowList, as.numeric)
  do.call(rbind, rowList)
}


parseNstates <- function (txt) {
  argname <- "nstates"
  matches <- grep(argname, txt, ignore.case = T)
  if (length(matches) != 1)
    stop(paste0("expecting one occurrence of the field '",
                argname, "', found ", length(matches)))
  if (matches + 1 > length(txt))
    stop(paste0("missing line below the argument '", argname,
                "'"))
  as.integer(txt[matches + 1])
}

searchOptArg <- function (txt, argname) {
  matches <- grep(argname, txt, ignore.case = T)
  if (length(matches) > 1)
    stop(paste0("expecting zero or one occurrences of the field '",
                argname, "', found ", length(matches)))
  matches + 1
}

ifHasField <- function (txt, field, nlines, f) {
  idx <- searchOptArg(txt, field)
  if (length(idx) == 1) {
    lastidx <- idx + nlines - 1
    if (lastidx > length(txt))
      stop(paste0("missing lines below the argument '",
                  field, "'"))
    f(txt[idx:lastidx])
  }
  else NULL
}

readModel <- function (path) {
  tryCatch({
    txt <- readLines(path)
    allowedIdxs <- 1:length(txt)
    nstates <- parseNstates(txt)
    marks <- ifHasField(txt, "marks", 1, function(lines) {
      strsplit(lines, "\t")[[1]]
    })
    labels <- NULL
    labels <- ifHasField(txt, "labels", 1, function(lines) {
      strsplit(lines, ",")[[1]]
    })
    colors <- NULL
    colors <- ifHasField(txt, "colors", 1, function(lines) {
      strsplit(lines, ",")[[1]]
    })
    firstline.E <- which(txt == "emisP") + 1
    if (grepl("\\|", strsplit(txt[firstline.E], "\t")[[1]][1]))
      modeltype <- "lognormal"
    else modeltype <- "negmultinom"
    if (modeltype == "negmultinom") {
      emisP <- ifHasField(txt, "emisP", nstates, function(lines) {
        emisMat <- parseRows(lines)
        lapply(1:nstates, function(i) {
          params <- emisMat[i, ]
          r = params[1]
          mus <- params[2:ncol(emisMat)]
          mu <- sum(mus)
          if (mu == 0) {
            ps <- rep(1/length(marks), length(marks))
          }
          else ps <- mus/mu
          list(mu = mu, r = r, ps = ps)
        })
      })
    }
    else if (modeltype == "lognormal") {
      emis.params <- as.numeric(unlist(strsplit(txt[firstline.E:(firstline.E +
                                                                   nstates - 1)], "[\t|]+")))
      mus <- matrix(emis.params[seq(1, length(emis.params),
                                    2)], nrow = nstates, byrow = T)
      sigmasqs <- matrix(emis.params[seq(2, length(emis.params),
                                         2)], nrow = nstates, byrow = T)
      emisP <- lapply(1:nstates, function(i) list(mus = mus[i,
      ], sigmasqs = sigmasqs[i, ]))
    }
    transP <- ifHasField(txt, "transP", nstates, parseRows)
    initP <- ifHasField(txt, "initP", nstates, parseRows)
    model <- list(nstates = nstates, marks = marks, emisP = emisP,
                  transP = transP, initP = initP)
    if (!is.null(labels))
      model$labels <- labels
    if (!is.null(colors))
      model$colors <- colors
    # validateModel(model, strict = TRUE, type = modeltype)
  }, error = function(e) {
    stop(paste0("Unable to parse the model file:\n", e$message))
  })
  model
}

getRef <- function (mat, type, nthreads = 1L) {
  .Call("_ehmm_getRef", PACKAGE = "ehmm", mat, type, nthreads)
}

quantileNorm <- function (mat, ref, nthreads = 1L, seed = 13L) {
  .Call("_ehmm_quantileNorm", PACKAGE = "ehmm", mat, ref,
        nthreads, seed)
}

statesToSegments_helper <- function (regions, states) {
  .Call("_ehmm_statesToSegments_helper", PACKAGE = "ehmm",
        regions, states)
}

statesToSegments <- function (states, regions) {
  h <- statesToSegments_helper(regions, states)
  GRanges(seqnames = h$chrs, IRanges(start = h$starts, end = h$ends,
                                     names = h$states))
}

validateVmat <- function (vmat) {
  if (!is.matrix(vmat))
    stop("vmat must be a matrix")
  if (!is.numeric(vmat))
    stop("'vmat' must be numeric")
  if (length(vmat) <= 0)
    stop("'vmat' must be non-empty")
}

quantileNormalization <- function (vmat, ref = c("median", "min", "mean"), nthreads = 1) {
  validateVmat(vmat)
  if (!is.numeric(ref)) {
    ref <- match.arg(ref)
    ref <- getRef(vmat, ref, nthreads)
  }
  if (length(ref) != nrow(vmat))
    stop("reference vector has the wrong length")
  quantileNorm(vmat, ref, nthreads = nthreads)
}

propagateErrors <- function (l) {
  for (el in l) if (inherits(el, "try-error"))
    stop(el)
  l
}

safe_mclapply <- function (...)
  propagateErrors(mclapply(...))

colSummary <- function (mat, type, nthreads = 1L) {
  .Call("_ehmm_colSummary", PACKAGE = "ehmm", mat, type, nthreads)
}

bindCols <- function (vlist, nthreads = 1L) {
  .Call("_ehmm_bindCols", PACKAGE = "ehmm", vlist, nthreads)
}

defaultRepFun <- function (vmat, normFun = quantileNormalization, nthreads = 1,
                           ...) {
  normFunOpts <- list(...)
  normFunOpts$vmat <- vmat
  if ("nthreads" %in% names(formals(normFun))) {
    normFunOpts[["nthreads"]] <- nthreads
  }
  nvmat <- do.call(normFun, normFunOpts)
  colSummary(t(nvmat), "median", nthreads = nthreads)
}

bamtabDefaults <- list(75, 0, FALSE)
names(bamtabDefaults) <- c("shift", "mapq", "pairedend")

validateBamtab <- function (bamtab) {
  if (!is.data.frame(bamtab))
    stop("'bamtab' must be a 'data.frame'")
  reqfields <- c("mark", "path")
  optfields <- "shift"
  if (!all(reqfields %in% names(bamtab)))
    stop("missing required fields")
  if (!all(names(bamtab) %in% c(reqfields, optfields)))
    stop("invalid fields")
  for (n in names(bamtabDefaults)) {
    if (!n %in% names(bamtab)) {
      bamtab[[n]] <- rep(bamtabDefaults[[n]], nrow(bamtab))
    }
  }
  if (!is.character(bamtab$path))
    stop("invalid path specification")
  if (any(!file.exists(bamtab$path)))
    stop("BAM file does not exist")
  if (!is.numeric(bamtab$mapq) || any(bamtab$mapq < 0 | bamtab$mapq >
                                      255)) {
    stop("invalid 'mapq'")
  }
  if (!is.numeric(bamtab$shift))
    stop("invalid 'shift'")
  if (!is.logical(bamtab$pairedend))
    stop("invalid 'pairedend'")
  bamtab
}

extractRegions <- function(segmentation, regions, genomeSize, outdir) {
  labels <- segmentation$model$labels
  gr <- unlist(tile(regions, width = 100))
  GenomeInfoDb:::seqlengths(gr) <- genomeSize[GenomeInfoDb:::seqlevels(gr)]
  gr$name <- labels[segmentation$viterbi]
  mapply(function(scores, elmName) {
    label <- toupper(substr(elmName, 1, 1))
    gr$score <- scores
    regionsBedfile <- sprintf("%s/%sRegions.bed", outdir,
                              elmName)
    elms.tiled <- gr[startsWith(gr$name, paste0(label, "_A"))]
    file.create(regionsBedfile)
    export.bed(elms.tiled, regionsBedfile)
  }, segmentation$score, c("enhancer", "promoter"))
}

applyModelLocal <- function (regions, model = NULL, provideModel = FALSE, genomeSize,
          counts = NULL, bamtab = NULL, outdir = ".", nthreads = 1,
          learnTrans = FALSE, refCounts = NULL) {
  binsize <- 100
  if (!is.null(model))
    provideModel <- FALSE
  if (is.null(model) && !provideModel) {
    cat("No model specified. Provided mESC model will be used.")
    provideModel <- TRUE
  }
  dir.create(file.path(outdir), showWarnings = FALSE)
  if (provideModel) {
    mark.acc <- unlist(sapply(c("acc", "atac", "dhs", "dnase"),
                              function(pattern) which(grepl(pattern, tolower(bamtab$mark)))))
    if (length(mark.acc) != 1) {
      stop("Error: Exactly one of the given mark names must contain \"atac\", \"dhs\" or \"dnase\".\n\n         If your chromatin accessibility assay is neither ATAC-seq nor DNase-seq, name it \"ACC\".\n\n         Example: --mark ACC:/path/chromatin_accessibility_assay.bam")
    }
    mark.k27ac <- which(grepl("k27ac", tolower(bamtab$mark)))
    mark.k4me1 <- which(grepl("k4me1", tolower(bamtab$mark)))
    mark.k4me3 <- which(grepl("k4me3", tolower(bamtab$mark)))
    load(system.file("extdata", "mESC.rdata", package = "ehmm",
                     mustWork = TRUE))
    model$marks <- bamtab$mark[c(mark.acc, mark.k27ac, mark.k4me1,
                                 mark.k4me3)]
    names(counts.tables) <- model$marks
    refCounts <- reconstructCountmatrix(counts.tables)
  }
  if (is.null(counts)) {
    counts <- getCountMatrix(regions, bamtab, outdir, binsize = 100,
                             nthreads, pseudoCount = 1)
  }
  if (!(is.null(refCounts))) {
    if (!all(as.vector(unique(seqnames(regions))) %in% names(genomeSize))) {
      names(genomeSize) <- gsub("chr", "", names(genomeSize))
    }
    counts <- quantileNormalizeCounts(counts = counts, refCounts = refCounts,
                                      regions = regions, genomeSize = genomeSize, bamtab = bamtab,
                                      outdir = outdir, nthreads = nthreads)
  }
  cat("applying model\n")
  segmentation <- segment(counts = counts, regions = regions,
                          model = model, nstates = model$nstates, nthreads = nthreads,
                          verbose_kfoots = TRUE, nbtype = "lognormal", notrain = !learnTrans,
                          fix_emisP = learnTrans)
  cat("producing report\n")
  viterbi_segments <- statesToSegments(segmentation$viterbi,
                                       segmentation$segments)
  report(segments = viterbi_segments, model = model, rdata = segmentation,
         outdir = outdir, colors = model$colors, labels = model$labels)
  cat("extract enhancer / promoter elements\n")
  extractRegions(segmentation, regions, genomeSize, outdir)
  cat("done\n\n")
}


args <- commandArgs(trailingOnly = T)
parser <- arg_parser("EHMM command line parser", name = "ehmm")
parser <- add_argument(parser, "-r", help = "BED file of input regions")
parser <- add_argument(parser, "-g", help = "Genome sizes")
parser <- add_argument(parser, "-m", help = "Construced model")
parser <- add_argument(parser, "-b", help = "Path to directory containing the bam files of the input data")
parser <- add_argument(parser, "-o", help = "Output directory", default = "./")
parser <- add_argument(parser, "-s", help = "BinSize", default = 100, type = "integer")
parser <- add_argument(parser, "-t", help = "Number of threads to use", default = 21, type = "integer")
parser <- add_argument(parser, "-p", help = "Pseudocount", default = 1.0, type = "double")

argv <- parse_args(parser, argv = args)
message("loading input data")
binsize <- argv$s
nthreads <- argv$t
pseudocount <- argv$p
regions <- import(argv$r, format = "bed")
# remove scaffold chromosomes
levelsToDrop <- unique(unlist(lapply(c("Un", "M", "random",
                                       "hap", "alt", "GL", "NC", "hs"), function(x) which(grepl(x,
                                                                                          GenomeInfoDb:::seqlevels(regions))))))
if (length(levelsToDrop) > 0)
  GenomeInfoDb::seqlevels(regions, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(regions)[-levelsToDrop]
# load bam files
bamfiles <- list.files(path = argv$b, pattern = "*.bam$", full.names = T, recursive = T)
marks <- basename(dirname(bamfiles))
marks[grep("ALL|RNA|ATC|POL", marks, ignore.case = T)] <- "acc"
bamtab <- data.frame(mark = marks, path = bamfiles, check.names = F)
bamtab$shift <- sapply(tolower(bamtab$mark), function(m) ifelse(any(sapply(c("atac", "dhs", "dnase", "acc"),
                                                                           function(s) grepl(s, m))), 0, 75))
# generate counts from bam files
bamtab <- validateBamtab(bamtab)
npaths <- nrow(bamtab)
if (npaths == 0)
  stop("no marks provided")
if (!inherits(regions, "GRanges"))
  stop("regions must be a GRanges object")
if (length(binsize) != 1 || binsize <= 0)
  stop("invalid binsize")
if (any(width(regions)%%binsize!=0)) {
  # select invalid regions
  invalidRegions <- regions[width(regions)%%binsize!=0]
  message("INFO: ", length(invalidRegions), " regions are not multiples of the selected bin size ", binsize, ", adjusting...")
  deltaLengths <- width(invalidRegions)%%binsize
  # pad regions to binsize
  deltaLengths <- sapply(deltaLengths, function(d) ifelse(d>binsize/2, binsize-d, -d))
  # replace end points
  data <- as.data.frame(invalidRegions)
  data$end <- data$end+deltaLengths
  data$width <- data$width+deltaLengths
  data <- makeGRangesFromDataFrame(data, keep.extra.columns = T)
  regions[width(regions)%%binsize!=0] <- data
}
message("getting counts")
countsList <- safe_mclapply(mc.cores = nthreads, 1:npaths,
                            function(i) {
                              b <- bamtab[i, ]
                              b$pairedend <- ifelse(b$pairedend, "midpoint", "ignore")
                              bprof <- bamProfile(b$path, regions, binsize = binsize,
                                                  shift = b$shift, mapq = b$mapq, paired.end = b$pairedend)
                              unlist(as.list(bprof))
                            })
if (any(duplicated(bamtab$mark))) {
  message("putting replicate experiments together")
  countsList <- lapply(setNames(nm = unique(bamtab$mark)),
                       function(nm) {
                         subList <- countsList[bamtab$mark == nm]
                         if (length(subList) == 1)
                           return(subList[[1]])
                         vmat <- bindCols(subList)
                         defaultRepFun(vmat)
                       })
} else {
  names(countsList) <- bamtab$mark
}
message("pileup done, storing everything in a matrix")
counts <- t(bindCols(countsList)) + pseudocount
message("loading model")
model <- readModel(argv$m)
model$marks <- toupper(model$marks)
rownames(counts) <- toupper(rownames(counts))
commonMarks <- intersect(rownames(counts), model$marks)
if (length(commonMarks) < length(unique(c(rownames(counts), model$marks)))) {
  message("INFO: only using marks that are present in both training and testing set:\n",
          "-> ", paste(commonMarks, collapse = ", "))
  if(nrow(counts)>1) counts <- counts[commonMarks,]
  model$marks <- commonMarks
  if (length(model$emisP[[1]]$mus) > length(commonMarks)) {
    model$emisP <- lapply(model$emisP, lapply, "[", which(commonMarks==model$marks))
  }
}
message("loading genome sizes")
genomeSizes <- import(argv$g, format = "bed")
genomeSizes <- setNames(genomeSizes@ranges@start, seqnames(genomeSizes))

applyModelLocal(regions, model = model, genomeSize = genomeSizes,
           bamtab = bamtab, counts = counts, nthreads = nthreads,
           outdir = argv$o)
message("ehmm finished successfully")
