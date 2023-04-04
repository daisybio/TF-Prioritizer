validateCounts <- function(counts){
    if (is.null(counts) || !(is.matrix(counts))) stop("'counts' must be a matrix")
    if (is.null(rownames(counts))) stop("rows of the 'counts' matrix must be named")
    if (anyDuplicated(rownames(counts))) stop("row names of the 'counts' matrix must be unique")
    
}

validateCList <- function(clist){
    sapply(clist, validateCounts)
    if (length(clist) < 2) return(clist)
    marks <- rownames(clist[[1]])
    ncols <- ncol(clist[[1]])
    #make count matrices compatible (i-th row is always the i-th mark)
    for (i in 2:length(clist)){
        newmarks <- rownames(clist[[i]])
        #check compatibility
        if (!setequal(marks, newmarks) || ncol(clist[[i]]) != ncols) {
            stop("non-compatible count matrices") }
        if (any(marks != newmarks)){
            #rearrange rows
            clist[[i]] <- clist[[i]][marks,]
        }
    }
    clist
}

readCounts <- function(path){
    if (isRdata(path)){
        counts <- readRdata(path)
    } else {
        counts <- tryCatch({ # in case of rpm-normalized counts, the countmatrix is in decimals instead of integers
          read.table(path, header=TRUE, sep="\t", colClasses='integer', row.names=NULL, check.names=FALSE)
          counts <- t(bindCols(countsList))
          }, error = function(err) {
            countsList <- read.table(path, header=TRUE, sep="\t", colClasses='double', row.names=NULL, check.names=FALSE)
            return(t(bindCols_numeric(countsList)))
            }
        )
    }
    validateCounts(counts)
    counts
}

writeCounts <- function(counts, path){
    validateCounts(counts)
    if (isRdata(path)){
        save(counts, file=path)
    } else {
        writeCountsTXT(counts, rownames(counts), path)
    }
}

writeCountsDouble <- function(counts, path){
  validateCounts(counts)
  if (isRdata(path)){
    save(counts, file=path)
  } else {
    writeCountsTXT_double(counts, rownames(counts), path)
  }
}

getCountMatrix <- function(regions, bamtab, outdir=NULL, binsize=100, nthreads=1, pseudoCount=1){
  # This script calculates a count matrix for given regions based on all bam-files in a given bam-directory and writes it to file (only if outdir is given).
  
  # make bamtab object
  # default shift is 75. set to 0 in case of ATAC / DHS
  bamtab$shift <- sapply(tolower(bamtab$mark), function(m) ifelse(any(sapply(c('atac', 'dhs', 'dnase', 'acc'), function(s) grepl(s, m))), 0, 75))
  
  # create count matrix
  counts <- getcounts(regions, bamtab, binsize=binsize, nthreads=nthreads) + pseudoCount
  
  # write count matrix
  if(!is.null(outdir)){
    target <- paste(outdir, 'countmatrix.txt', sep='/')
    cat(sep="", "writing count matrix to the file '", target, "'\n")
    writeCountsDouble(counts, target)
  }
  
  return(counts)
}
