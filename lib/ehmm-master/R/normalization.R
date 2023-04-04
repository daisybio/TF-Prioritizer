validateVmat <- function(vmat){
    if (!is.matrix(vmat)) stop("vmat must be a matrix")
    if (!is.numeric(vmat)) stop("'vmat' must be numeric")
    if (length(vmat) <= 0) stop("'vmat' must be non-empty")
}

quantileNormalization <- function(vmat, ref=c("median", "min", "mean"), nthreads=1){
    validateVmat(vmat)
    if (!is.numeric(ref)){
        ref <- match.arg(ref)
        ref <- getRef(vmat, ref, nthreads)
    }
    if (length(ref) != nrow(vmat)) stop("reference vector has the wrong length")
    
    quantileNorm(vmat, ref, nthreads=nthreads)
}

quantileNormalizeToReference <- function(cm.reference, cm.query){
  # This function quantile-normalizes a query count matrix to a reference count matrix,
  # i.e. it sorts both distributions and deploys the values of the reference distribution to the entries of the query distribution with respect to their rank.
  cm.query.normalized <- matrix(nrow=nrow(cm.query), ncol=ncol(cm.query))
  row.names(cm.query.normalized) <- row.names(cm.query)
  dict.list <- vector(mode="list", length=nrow(cm.query))
  for (i in 1:nrow(cm.query)){
    target <- preprocessCore::normalize.quantiles.determine.target(as.matrix(cm.reference[i,]))
    query <- as.matrix(cm.query[i,])
    query.normalized <- preprocessCore::normalize.quantiles.use.target(x=query, target=target)
    cm.query.normalized[i,] <- query.normalized
    dict <- as.vector(query.normalized)
    names(dict) <- as.vector(query)
    dict.list[[i]] <- dict[unique(names(dict))]
    dict.list[[i]] <- dict.list[[i]][order(as.integer(names(dict.list[[i]])))] # sort by names (= query counts)
  }
  return(list(cm.query.normalized=cm.query.normalized, dict.list=dict.list))
}

quantileNormalizeCounts <- function(counts, refCounts, regions, genomeSize, bamtab=NULL, outdir, nthreads){
  # Preparatory function for the actual quantile normalization function that is called within this function.
  # Checks wether the all the regions are full chromosomes, and if not, calculates a countmatrix on the full target genome
  # to enable a comparison to the reference. Counts are clipped to the 99.9 percentile.
  regions.full <- GRanges(seqnames=names(genomeSize), IRanges(start=101, end=as.integer(genomeSize/100)*100))
  if (any(is.na(match(regions, regions.full)))) { # this checks for every entry in regions whether there is a equal entry in regions.full, dealing with the possibility of different seqlevels.
    cat('Calculate counts for full genome\n') # always do normalization on countmatrix for full genome
    if(is.null(bamtab)) stop('Error: no directory with bam-files was passed. It is necessary to calculate count matrix for the full genome 
                             because the passed regions file does not span the full genome.')
    counts.full <- getCountMatrix(regions=regions.full, bamtab=bamtab, binsize=100, nthreads=nthreads, pseudoCount=1) # not written to file without passed outdir argument
    counts.full.clipped <- clipCounts(counts.full, .999)
    counts.clipped <- counts # here (in the next line) we just set any count greater than the rowwise max of counts.full.clipped to that exact max.
    for (i in 1:nrow(counts.clipped)) counts.clipped[i, (counts[i,] > max(counts.full.clipped[i,]))] <- max(counts.full.clipped[i,])
  } else {
    counts.full <- counts # if the test above fails, the given counts are treated as full genome (technically they could be only for one full chromosome, but that's also fine)
    counts.full.clipped <- clipCounts(counts.full, .999)
    counts.clipped <- counts.full.clipped # here, we actually just reassign counts.full.clipped because counts and counts.full are the same.
  }
  refCounts.clipped <- clipCounts(refCounts, .999)
  cat('normalizing count matrix to reference\n')
  res <- quantileNormalizeToReference(cm.reference=refCounts.clipped, cm.query=counts.full.clipped)
  rnames <- row.names(counts)
  # normalize the actual countmatrix. deal with counts not present in the dict (e.g. due to shifted regions compared to the full genome countmatrix):
  # interpolate with values closest two dict-keys and add to the dict.
  for (i in 1:nrow(counts.clipped)){
    countsNotInDict <- setdiff(unique(counts.full.clipped[i,]), names(res$dict.list[[i]]))
    for (count in countsNotInDict){
      nearestValues <- res$dict.list[[i]][order(abs(as.numeric(names(res$dict.list[[i]])) - countsNotInDict))[1:2]]
      interpolatedValue <- approxfun(c(names(nearestValues)), c(nearestValues))(count)
      res$dict.list[[i]][as.character(count)] <- interpolatedValue
    }
  }
  counts.normalized <- t(sapply(1:nrow(counts.clipped), function(i) as.vector(res$dict.list[[i]][as.character(counts.clipped[i,])])))
  row.names(counts.normalized) <- rnames
  # write normalized count matrix to file
  filename <- paste(outdir, 'countmatrix_normalized.txt', sep='/')
  cat(sep="", "writing normalized count matrix to the file '", filename, "'\n")
  writeCountsDouble(counts.normalized, filename)
  return(counts.normalized)
}

defaultSFFun <- function(vmat, method="RLE", ...){
    edgeR::calcNormFactors(vmat, method=method, ...)
}

linearNormalization <- function(vmat, sfFun=defaultSFFun, ...){
    validateVmat(vmat)
    #compute scaling factors
    sf <- sfFun(vmat, ...)
    
    #scale vectors
    res <- round(vmat*sf[col(vmat)])
    storage.mode(res) <- "integer"
    res
}

defaultRepFun <- function(vmat, normFun=quantileNormalization, nthreads=1, ...){
    normFunOpts <- list(...)
    normFunOpts$vmat <- vmat
    if ("nthreads" %in% names(formals(normFun))) {
        normFunOpts[["nthreads"]] <- nthreads
    }
    nvmat <- do.call(normFun, normFunOpts)
    colSummary(t(nvmat), "median", nthreads=nthreads)
}

# rpmNormalizeCounts <- function(counts, bamtab, binsize, pseudoCount){
#   # rpm normalization requires total number of reads obtained from the bam files.
#   # note: add (pseudoCount * nbins) to the total number of reads.
#   # nbins is the theoretical number of bins in the whole genome given by the total seqlength divided by the binsize.
#   nreads.total <- rep(0,length(bamtab$mark)); names(nreads.total) <- bamtab$mark
#   bamstats <- sapply(1:nrow(bamtab), function(i) Rsamtools:::idxstatsBam(bamtab$path[i]))
#   seqlengths.total <- sapply(1:nrow(bamtab), function(i) sum(as.numeric(bamstats['seqlength',][[i]])))
#   nbins <- seqlengths.total / binsize
#   nreads.total <- sapply(1:nrow(bamtab), function(i) sum(bamstats['mapped',][[i]])) + nbins
#   rpmCounts <- counts / nreads.total * 10^6
#   return(rpmCounts)
# }
