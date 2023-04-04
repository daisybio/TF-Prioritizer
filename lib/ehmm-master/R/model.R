#a model object is an R list that encodes all the parameters of the HMM
#emission probabilities: a list of parameters for the negative multinomials
#transition probabilities: a square matrix, rows must sum up to one
#initial probabilities: a matrix where the rows corresponds to the states and the columns to the sequences
#columns must sum up to one
#marks: mark names, used as ids and needed to check consistency with the count matrix

#checks that the model object has the right format
#every field is optional if input=FALSE, otherwise every field is mandatory
#fields must be consistent
#return a list(nstates, nmarks)
validateModel <- function(model, strict=FALSE, type="dep"){
  if (is.null(model)) stop("NULL is not a valid model")
    optArgs <- c("nstates", "marks", "labels", "colors", "emisP", "transP", "initP")
    if (is.null(names(model))) stop("model's element must be named")
    # if (!all(names(model) %in% optArgs)) warning("some arguments will be ignored (incorrect naming?)")
    if (!strict && !all(optArgs %in% names(model))) stop(paste("a model requires the following fields:", paste(optArgs, collapse=", ")))
    
    nstates <- model$nstates #can be NULL
    nmarks <- NULL
    if (!is.null(model$marks)) {
        if (anyDuplicated(model$marks)) stop("model$marks cannot contain duplicate entries")
        nmarks <- length(model$marks)
    }
    #check emission probabilities
    if (!is.null(model$emisP)){
        if (!is.list(model$emisP)) stop("'model$emissP' must be a list'")
        if (!is.null(nstates)){
            if (length(model$emisP)!=nstates) stop("incorrect number of emission probabilities provided")
        } else nstates <- length(model$emisP)
	for (nm in model$emisP){
	    if (type == "lognormal") {
	       needPars <- c("mus", "sigmasqs")
	       if (!all(needPars %in% names(nm))) stop("missing fields from the emission probabilities")
               if (any(is.na(unlist(nm)))) stop("NAs/NaNs in the emission probabilities not allowed")
               if (is.null(nmarks)) nmarks <- length(nm$mus)
               if (length(nm$mus) != nmarks) stop("invalid parameters for the emission probabilities")
	    }
	    else {
	       needPars <- c("mu", "r", "ps")
	       if (!all(needPars %in% names(nm))) stop("missing fields from the emission probabilities")
               if (any(is.na(unlist(nm)))) stop("NAs/NaNs in the emission probabilities not allowed")
               if (is.null(nmarks)) nmarks <- length(nm$ps)
               if (length(nm$ps) != nmarks) stop("invalid parameters for the emission probabilities")
               if (!isProbVector(nm$ps)) stop("'model$emissP[[i]]$ps' must sum up to 1 for every i")
	    }
	}
    }
    #check transition probabilities
    if (!is.null(model$transP)){
        #validate transition probabilities (additional validation performed in kfoots)
        if (!is.matrix(model$transP) || ncol(model$transP)!=nrow(model$transP) ||
                (!is.null(nstates) && nstates != ncol(model$transP))) stop("invalid trasition probabilities")
        if (!all(apply(model$transP, 1, isProbVector))) stop("rows of model$transP must sum up to 1")
        if (is.null(nstates)) nstates <- nrow(model$transP)
    }
    
    #check initial probabilities
    if (!is.null(model$initP)){
        if (!is.matrix(model$initP) || (!is.null(nstates) && nrow(model$initP) != nstates)) stop("invalid initial probabilities")
        if (!all(apply(model$initP, 2, isProbVector))) stop("columns of model$initP must sum up to 1")
        if (is.null(nstates)) nstates <- nrow(model$initP)
    }
    list(nstates=nstates, nmarks=nmarks)
}


#reorder the marks in the model
reoderMarks <- function(model, mperm){
    if (length(model$marks) != length(mperm)) stop("invalid model permutation")
    model$marks <- model$marks[mperm]
    if (!is.null(model$emisP)){
        for (i in seq_along(model$emisP)){
            model$emisP[[i]]$ps <- model$emisP[[i]]$ps[mperm]
        }
    }
    model
}

isProbVector <- function(v, tol=1e-9){
    all(v >= 0) && abs(sum(v)-1) <= tol
}


#rows: marks, columns: states
getPS <- function(emisP){
    do.call(cbind, lapply(emisP, function(m) m$ps))
}

getmus <- function(emisP){
  do.call(cbind, lapply(emisP, function(m) m$mus))
}


#first row: mu
#second row: r
#columns: states
getNBs <- function(emisP){
    do.call(cbind, lapply(emisP, function(m) c(m$mu, m$r)))
}

#rows: marks, columns: states
getMeanMatrix <- function(emisP){
  ps <- getPS(emisP)
  mus <- getNBs(emisP)[1,]
  ps * mus[col(ps)]
}


plotNBs <- function(nbs, eps=0.2, xlab="count",...){
    nbs <- nbs[,ncol(nbs):1, drop=FALSE]
    ys <- seq(0,1, length.out=ncol(nbs))
    xs <- nbs[1,]
    sds <- sqrt(xs + (xs^2)/nbs[2,])
    xlim <- c(min(xs-sds), max(xs+sds))
    if (ncol(nbs) <= 1) { spacer <- 0.5
    } else                spacer <- 1/(2*(ncol(nbs)-1))
    ylim <- c(0-spacer,1+spacer)
    par(yaxs="i", las=1, mar=c(5, 7, 4, 1))
    plot(NULL, xlim=xlim, ylim=ylim, yaxt="n", xlab=xlab, ...)
    axis(side=2, tick=F, at=ys, labels=colnames(nbs))
    X <- as.numeric(rbind(xs-sds, xs+sds, NA))
    Y <- as.numeric(rbind(ys, ys, NA))
    lines(X, Y)
    X <- as.numeric(rbind(xs-sds, xs-sds, NA))
    Y <- as.numeric(rbind(ys - spacer*eps, ys + spacer*eps, NA))
    lines(X, Y)
    X <- as.numeric(rbind(xs+sds, xs+sds, NA))
    lines(X, Y)
    points(xs, ys)
}

parseNstates <- function(txt){
    argname <- "nstates"
    matches <- grep(argname, txt, ignore.case=T)
    if (length(matches)!=1) 
        stop(paste0("expecting one occurrence of the field '", argname, "', found ", length(matches)))
    if (matches+1 > length(txt)) stop(paste0("missing line below the argument '", argname, "'"))
    as.integer(txt[matches+1])
}

searchOptArg <- function(txt, argname){
    matches <- grep(argname, txt, ignore.case=T)
    if (length(matches)>1) 
        stop(paste0("expecting zero or one occurrences of the field '", argname, "', found ", length(matches)))
    matches+1
}

ifHasField <- function(txt, field, nlines, f){
    idx <- searchOptArg(txt, field)
    if (length(idx)==1) {
        #the argument is there, the data starts at line idx
        lastidx <- idx + nlines -1
        if (lastidx > length(txt)) stop(paste0("missing lines below the argument '", field, "'"))
        
        f(txt[idx:lastidx])
    } else NULL
}

parseRows <- function(rows){
    rowList <- strsplit(rows, "\t")
    rowLens <- sapply(rowList, length)
    if (any(rowLens != rowLens[1])) stop("matrix rows don't have the same length")
    rowList <- lapply(rowList, as.numeric)
    do.call(rbind, rowList)
}

rowsToStr <- function(mat){
    apply(mat, 1, paste, collapse="\t")
}


readModel <- function(path){
    tryCatch({
        txt <- readLines(path)
        allowedIdxs <- 1:length(txt)
        #parsing nstates
        nstates <- parseNstates(txt)
        #parsing marks
        marks <- ifHasField(txt, "marks", 1, function(lines){
            strsplit(lines, "\t")[[1]]
        })
        #parsing labels
        labels <- NULL
        labels <- ifHasField(txt, "labels", 1, function(lines){
            strsplit(lines, ",")[[1]]
        })
        #parsing colors
        colors <- NULL
        colors <- ifHasField(txt, "colors", 1, function(lines){
          strsplit(lines, ",")[[1]]
        })
        #parsing emisP
        firstline.E <- which(txt == 'emisP') + 1
        if (grepl('\\|', strsplit(txt[firstline.E], '\t')[[1]][1])) modeltype <- 'lognormal' else modeltype <- 'negmultinom'
        if (modeltype == 'negmultinom'){
          emisP <- ifHasField(txt, "emisP", nstates, function(lines){
            emisMat <- parseRows(lines)
            lapply(1:nstates, function(i){
                params <- emisMat[i,]
                r=params[1]
                mus <- params[2:ncol(emisMat)]
                mu <- sum(mus)
                if (mu==0) { ps <- rep(1/length(marks), length(marks))
                } else ps <- mus/mu
                list(mu=mu, r=r, ps=ps)
            })
          })
        } else if (modeltype == 'lognormal'){
          emis.params <- as.numeric(unlist(strsplit(txt[firstline.E:(firstline.E+nstates-1)], "[\t|]+")))
          mus <- matrix(emis.params[seq(1,length(emis.params),2)], nrow=nstates, byrow=T)
          sigmasqs <- matrix(emis.params[seq(2,length(emis.params),2)], nrow=nstates, byrow=T)
          emisP <- lapply(1:nstates, function(i) list(mus=mus[i,], sigmasqs=sigmasqs[i,]))
        }
        #parsing transP
        transP <- ifHasField(txt, "transP", nstates, parseRows)
        #parsing initP
        initP <- ifHasField(txt, "initP", nstates, parseRows)
        
        model <- list(nstates=nstates, marks=marks, emisP=emisP, transP=transP, initP=initP)
        if (!is.null(labels)) model$labels <- labels
        if (!is.null(colors)) model$colors <- colors
        validateModel(model, strict=TRUE, type=modeltype)
    }, error=function(e){
        stop(paste0("Unable to parse the model file:\n", e$message))
    })
    model
}

writeModel <- function(model, path, type){
    #determine fit type
    if (all(names(model$emisP[[1]]) %in% c("mus", "sigmasqs"))) type <- "lognormal" else type <- "negmultinom"
    validateModel(model, strict=TRUE, type=type)
    txt <- list()
    #nstates to strings
    if (!is.null(model$nstates)) txt$nstates <- c("nstates", model$nstates)
    #labels to strings
    if(!is.null(model$labels)) txt$labels <- c("labels", paste0(model$labels, collapse=","))
    #colors to strings
    if(!is.null(model$colors)) txt$colors <- c("colors", paste0(model$colors, collapse=","))
    #marks to strings
    if (!is.null(model$marks)) txt$marks <- c("marks", paste0(model$marks, collapse="\t"))
    #emisP to strings
    if (!is.null(model$emisP)){
	if (type == "lognormal"){
	   emisMat <- t(do.call(cbind, lapply(model$emisP, function(m) paste(m$mus, m$sigmasqs, sep="|"))))
    	}
    	else{
	   nbs <- getNBs(model$emisP)
           rs <- nbs[2,]
           meansMat <- getMeanMatrix(model$emisP)
           emisMat <- t(rbind(rs, meansMat))
    	}
    	txt$emisP <- c("emisP", rowsToStr(emisMat))
    }
    #transP to strings
    if (!is.null(model$transP)) txt$transP <- c("transP", rowsToStr(model$transP))
    #initP to strings
    if (!is.null(model$initP)) txt$initP <- c("initP", rowsToStr(model$initP))
    
    writeLines(unlist(txt), path)
}

initializeParams <- function(model, states.a, states.n){
  # this function extracts given accessibility / nucleosome states from a given model and constructs a model with a N1-A-N2 architecture
  n.acc <- length(states.a)
  n.nuc <- length(states.n)
  model$nstates <- n.acc + 2 * n.nuc
  newStates.n1 <- 1:n.nuc
  newStates.a <- 1:n.acc + max(newStates.n1)
  newStates.n2 <- newStates.n1 + max(newStates.a)
  
  # emission probabilities: adopted from the given model
  model$emisP <- model$emisP[c(states.n, states.a, states.n)]
  
  # transition probabilities: initial guess independent from the given model. will be refined (learned) in the next step.
  model$transP <- matrix(0, nrow=model$nstates, ncol=model$nstates)
  model$transP[newStates.n1, newStates.n1] <- 0.5 / n.nuc
  model$transP[newStates.n1, newStates.a] <- 0.5 / n.acc
  model$transP[newStates.a, newStates.a] <- 0.9 / n.acc
  model$transP[newStates.a, newStates.n2] <- 0.1 / n.nuc
  model$transP[newStates.n2, newStates.n2] <- 1 / n.nuc
  
  # initial probabilities: initial guess of equal probabilities for N1 states, 0 for the rest.
  model$initP <- matrix(c(rep(1/n.nuc, n.nuc), rep(0, n.acc+n.nuc)))

  # define endstates (N2*)
  model$endstates <- newStates.n2
  
  return(model)
}

combineFgBgModels <- function(model.bg, model.e, model.p){
  # this function combines a background, enhancer and promoter model, setting 'forbidden' transitions to zoer, e.g. bg -> accessible
  
  # combine labels, colors, set marks and nstates
  labels <- c(model.e$labels, model.p$labels, model.bg$labels)
  colors <- c(model.e$colors, model.p$colors, model.bg$colors)
  nstates <- sum(model.e$nstates, model.p$nstates, model.bg$nstates)
  marks <- model.e$marks
  
  # combine emisP lists
  E <- c(model.e$emisP, model.p$emisP, model.bg$emisP)
  
  # combine transP matrices
  nEnhancers <- 399124 # Bernstein et al. 2012, ENCODE, https://www.nature.com/articles/nature11247: Number of regions with enhancer-like features.
  nPromoters <- 70292 # Bernstein et al. 2012, ENCODE, https://www.nature.com/articles/nature11247: Number of regions with promoter-like features.
  nBins <- 30000000 # hg19: 30'956'738, mm10: 27'255'182
  enhancerFreq <- nEnhancers / nBins
  promoterFreq <- nPromoters / nBins
  states.n1.e <- which(startsWith(labels, 'E_N1')); states.a.e <- which(startsWith(labels, 'E_A')); states.n2.e <- which(startsWith(labels, 'E_N2'))
  states.n1.p <- which(startsWith(labels, 'P_N1')); states.a.p <- which(startsWith(labels, 'P_A')); states.n2.p <- which(startsWith(labels, 'P_N2'))
  states.bg <- which(startsWith(labels, 'bg'))
  A <- as.matrix(bdiag(model.e$transP, model.p$transP, model.bg$transP))
  A[states.bg, states.n1.e] <- enhancerFreq / (length(states.bg) * length(states.n1.e))
  A[states.bg, states.n1.p] <- promoterFreq / (length(states.bg) * length(states.n1.p))
  A[states.bg, states.bg] <- A[states.bg, states.bg] * (1 - rowSums(A[states.bg, c(states.n1.e, states.n1.p)]))
  A[states.n2.e, states.bg] <- rowSums(A[states.n1.e, states.a.e]) / length(states.bg)
  A[states.n2.p, states.bg] <- rowSums(A[states.n1.p, states.a.p]) / length(states.bg)
  A[states.n2.e, states.n2.e] <- A[states.n2.e, states.n2.e] * (1 - rowSums(A[states.n2.e, states.bg]))
  A[states.n2.p, states.n2.p] <- A[states.n2.p, states.n2.p] * (1 - rowSums(A[states.n2.p, states.bg]))
  
  # set uniform initial probs for 'allowed' states (bg, N1)
  I <- matrix(rep(0,nstates))
  I[c(states.n1.e, states.n1.p, states.bg)] <- 1 / length(c(states.n1.e, states.n1.p, states.bg))
  
  # create model object
  model <- list(nstates=nstates, marks=marks, emisP=E, transP=A, initP=I, labels=labels, colors=colors)
  return(model)
}
