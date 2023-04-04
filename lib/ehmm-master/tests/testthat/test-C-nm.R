context("C: Negative multinomial distribution")

#R EQUIVALENTS TO THE COMPILED C++ FUNCTIONS
their_nmllikfun <- function(counts, model){
    dnbinom(log=T, colSums(counts), mu=model$mu, size=model$r) + apply(log=T, counts, 2, dmultinom, prob=model$ps)
}

#returns a braketing triple for
#function f with domain (0, +Inf)
#returns a triple such that:
#1: xlow < xmid < xup
#2: flow < fmid > fup
#if the maximum is at +Inf then xmax=Inf
#if the maximum is at 0 then xmin=0
their_braketMaximum <- function(initX, f, rate=1.6){
    if (rate <= 1) stop("invalid rate provided")
    if (!is.finite(initX) || initX/rate <= 0) stop("invalid initX provided")
    xmid <- initX; fmid <- f(xmid)
    xlow <- min(xmid/rate, xmid-1e-10); flow <- f(xlow)
    xup <- xmid*rate; fup <- f(xup)
    #extend braketing to the left
    while (flow >= fmid){
        xup <- xmid; fup <- fmid
        xmid <- xlow; fmid <- flow
        xlow <- min(xlow/rate, xlow-1e-10)
        if (xlow<=0) {xlow <- 0; break}
        flow <- f(xlow)
    }
    #extend braketing to the right
    while (fup >= fmid){
        xlow <- xmid; flow <- fmid;
        xmid <- xup; fmid <- fup;
        xup <- rate*xup
        if (xup==Inf) break
        fup <- f(xup)
    }
    
    list(xlow=xlow, xmid=xmid, xup=xup, flow=flow, fmid=fmid, fup=fup)
}

their_fitNB <- function(cs, wts=rep(1, length(cs)), tol=1e-9){
    mn <- sum(cs*wts)/sum(wts)
    vr <- sum(wts*((cs - mn)^2))/sum(wts)
    #poisson case
    if (vr < mn) return(list(mu=mn, r=Inf))
    
    #optimize the weighted log likelihood
    optimfun <- function(r) sum(dnbinom(cs, mu=mn, size=r, log=T)*wts)
    
    rinit <- (mn^2)/(vr-mn)#methods-of-moment estimate
    
    btrip <- their_braketMaximum(rinit, optimfun)
    if (btrip$xlow<=0||!is.finite(btrip$xup)) stop("unexpected result")
    
    #call brent algo
    r = optimize(optimfun, lower=btrip$xlow, upper=btrip$xup, maximum=T, tol=tol)$maximum
    list(mu=mn, r=r)
}

their_fitMultinom <- function(counts, wts=rep(1, ncol(counts))){
    mns <- rowSums(counts*wts[col(counts)])
    ps <- mns/sum(mns)
}

their_fitNM <- function(counts, wts=rep(1, ncol(counts)), tol=1e-9){
    cs <- colSums(counts)
    nb <- their_fitNB(cs, wts=wts, tol=tol)
    nb$ps <- their_fitMultinom(counts, wts)
    nb
}

their_fitNBs_1r <- function(cs, wtsmat, tol=1e-9){
    nmods <- nrow(wtsmat)
    mns <- sapply(1:nmods, function(i){
        sum(cs*wtsmat[i,])/sum(wtsmat[i,])
    })
    vrs <- sapply(1:nmods, function(i){
        sum(wtsmat[i,]*(cs-mns[i])^2)/sum(wtsmat[i,])
    })
    
    #optimize the weighted log likelihood for all models, using the same r
    optimfun <- function(r) {
        lliks <- t(sapply(1:nmods, function(i){
            if (r==Inf){
                dpois(cs, lambda=mns[i], log=T)
            } else {
                dnbinom(cs, mu=mns[i], size=r, log=T)
            }
        }))
        sum(lliks*wtsmat)
    }
    
    #initial guess for r
    valid <- vrs < mns
    if (sum(valid)==0) {
        rinit <- 1
    } else {
        initRs <- (mns[valid]^2)/(vrs[valid]-mns[valid])
        rinit <- median(initRs)
    }
    
    #find a bracketing triple
    btrip <- their_braketMaximum(rinit, optimfun)
    if (btrip$xlow<=0) stop("unexpected result")
    if (btrip$xup==Inf) {
        r = Inf
    } else{
        r <- optimize(optimfun, lower=btrip$xlow, upper=btrip$xup, maximum=T, tol=tol)$maximum
    }
    
    lapply(1:nmods, function(i) list(mu=mns[i], r=r))
}

their_fitNMs_1r <- function(counts, wtsmat, tol=1e-9){
    nmods <- nrow(wtsmat)
    nms <- their_fitNBs_1r(colSums(counts), wtsmat, tol=tol)
    mns <- lapply(1:nmods, function(i) {
        their_fitMultinom(counts, wtsmat[i,])
    })
    lapply(1:nmods, function(i){
        nms[[i]]$ps <- mns[[i]]
        nms[[i]]
    })
}

#START TESTING

test_that("log likelihood of the NM works",{
    #create some random data: a 11*1000 matrix of counts
    set.seed(17)
    counts <- exampleData(1000, indep=T)
    
    #test the negative binomial alone first
    cs <- colSums(counts)
    model <- list(mu=4, r=2, ps=1)
    their_csllik <- dnbinom(log=T, cs, mu=model$mu, size=model$r)
    #mapToUnique should already have been tested
    my_csllik <- numeric(length(cs))
    ucs <- mapToUnique(cs)
    mConst <- rep(0, length(cs))
    lLikMat(matrix(cs, nrow=1), list(model), ucs=ucs, mConst=mConst, lliks=my_csllik, nthreads=1)
    expect_equal(my_csllik, their_csllik)
    #try with more threads
    lLikMat(matrix(cs, nrow=1), list(model), ucs=ucs, mConst=mConst, lliks=my_csllik, nthreads=20)
    expect_equal(my_csllik, their_csllik)
    #check that if r==Inf you get a poisson
    pmodel <- list(mu=4, r=Inf, ps=1)
    their_csllik <- dpois(log=T, cs, lambda=pmodel$mu)
    lLikMat(matrix(cs, nrow=1), list(pmodel), ucs, mConst=mConst, lliks=my_csllik)
    expect_equal(my_csllik, their_csllik)
    #also the extreme-case/bug where mu*r==Inf (it should be fixed in R3.1.1)
    pmodel <- list(mu=5333.588, r=1.52916e+305, ps=1)
    their_csllik <- dpois(log=T, cs, lambda=pmodel$mu)
    lLikMat(matrix(cs, nrow=1), list(pmodel), ucs, mConst=mConst, lliks=my_csllik)
    expect_equal(my_csllik, their_csllik)
    
    #test multinomConst
    their_mConst <- lfactorial(colSums(counts)) - colSums(lfactorial(counts))
    my_mConst <- getMultinomConst(counts, nthreads=1)
    expect_equal(my_mConst, their_mConst)
    #with more threads
    my_mConst <- getMultinomConst(counts, nthreads=20)
    expect_equal(my_mConst, their_mConst)
    
    #test the negative multinomial with one model
    model$ps <- c(1,2,3,4,0.1,19,2,4,9,1,7)
    model$ps <- model$ps/sum(model$ps)
    my_nmllik <- numeric(length(cs))
    lLikMat(counts, list(model), ucs=ucs, mConst=my_mConst, lliks=my_nmllik, nthreads=1)
    their_nmllik <- their_nmllikfun(counts, model)
    expect_equal(my_nmllik, their_nmllik)
    #with more threads
    lLikMat(counts, list(model), ucs=ucs, mConst=my_mConst, lliks=my_nmllik, nthreads=20)
    expect_equal(my_nmllik, their_nmllik)
    #add some more models
    model2 = list(mu=40, r=0.4, ps=c(1,8,5,8,5,6,5,4,3,2,1))
    model3 = list(mu=20, r=2, ps=c(1,1,1,1,1,3,4,5,6,5,4))
    model2$ps = model2$ps/sum(model2$ps)
    model3$ps = model3$ps/sum(model3$ps)
    models <- list(model, model2, model3)
    their_nmllik <- t(sapply(models, their_nmllikfun, counts=counts))
    my_nmllik <- matrix(0, nrow=length(models), ncol=ncol(counts))
    lLikMat(counts, models, ucs=ucs, mConst=my_mConst, lliks=my_nmllik, nthreads=1)
    expect_equal(my_nmllik, their_nmllik)
    #with more threads
    lLikMat(counts, models, ucs=ucs, mConst=my_mConst, lliks=my_nmllik, nthreads=20)
    expect_equal(my_nmllik, their_nmllik)
})


test_that("fitting the NM works",{
    #create some random data: a 11*1000 matrix of counts
    set.seed(71)
    counts <- exampleData(1000, indep=T)
    
    #test the negative binomial alone first
    
    #we need to lower the error tolerance for the r parameter,
    #no way of getting errors less or equal to 1e-8...
    their_nbs <- apply(counts, 1, their_fitNB, tol=1e-9)
    my_nbs <- apply(counts, 1, fitNB, tol=1e-9)
    expect_equal(my_nbs, their_nbs, tol=1e-7)
    
    #repeat it with more threads
    my_nbs2 <- apply(counts, 1, fitNB, tol=1e-9, nthreads=20)
    expect_equal(my_nbs, my_nbs2, tol=1e-7)
    
    #now with some random weights
    wts <- abs(rnorm(ncol(counts)))
    
    their_nbs <- apply(counts, 1, their_fitNB, wts=wts, tol=1e-9)
    my_nbs <- apply(counts, 1, fitNB, posteriors=wts,  tol=1e-9)
    expect_equal(my_nbs, their_nbs, wts=wts,  tol=1e-7)
    
    #repeat it with more threads
    my_nbs2 <- apply(counts, 1, fitNB, posteriors=wts, tol=1e-9, nthreads=20)
    expect_equal(my_nbs, my_nbs2, tol=1e-7)
    
    #underdispersed data: we should get a poisson (r=Inf)
    cs <- rbinom(1000, size=100, prob=0.5)
    my_pois <- fitNB(cs)
    expect_true(my_pois$r==Inf)
    
    if (F){
        #this test should be written in C: you can't guarantee that scorefun 
        #is IDENTICAL to the one at the C level...
        scorefun <- function(cs, wts, model) sum(dnbinom(cs, mu=model$mu, size=model$r, log=T)*wts)
        #even with a good starting model, we need the guarantee that
        #the fitted one will not have a lower score than the starting one
        #use their_nms as starting models
        for (i in seq_along(their_nms)){
            cs <- counts[i,]
            wts <- rep(1, length(cs))
            ucs <- mapToUnique(cs)
            cs2 <- ucs$values
            wts2 <- sumAt(wts, ucs$map, length(cs2), zeroIdx=T)
            startmod <- their_nms[[i]]
            fitmod <- fitNB(cs, wts, old_r=startmod$r, tol=1e-8)
            expect_true(scorefun(cs2, wts2, startmod) <= scorefun(cs2, wts2, fitmod))
        }
    }
    
    #test the negative multinomial
    
    #with only one model
    ucs <- mapToUnique(colSums(counts))
    emptyModel <- list(mu=-1, r=-1, ps=-1)
    their_nm <- their_fitNM(counts, wts, tol=1e-9)
    my_nm <- fitModels(counts, matrix(wts, nrow=1), list(emptyModel), ucs=ucs, tol=1e-9)[[1]]
    expect_equal(my_nm, their_nm, tol=1e-7)
    #with many models
    wtsmat <- matrix(abs(rnorm(ncol(counts)*nrow(counts))), nrow=nrow(counts))
    their_nms <- lapply(1:nrow(counts), function(i) their_fitNM(counts, wtsmat[i,], tol=1e-9))
    my_nms <- fitModels(counts, wtsmat, lapply(1:nrow(counts), function(i) emptyModel), ucs=ucs, tol=1e-9)
    expect_equal(my_nms, their_nms, tol=1e-7)
    #with many threads
    my_nm2 <- fitModels(counts, matrix(wts, nrow=1), list(emptyModel), ucs=ucs, tol=1e-9, nthreads=20)[[1]]
    expect_equal(my_nm, my_nm2, tol=1e-7)
    my_nms2 <- fitModels(counts, wtsmat, lapply(1:nrow(counts), function(i) emptyModel), ucs=ucs, tol=1e-9, nthreads=20)
    expect_equal(my_nms, my_nms2, tol=1e-7)
    
    #try out the nbtype="dep" option
    their_nms <- their_fitNMs_1r(counts, wtsmat, tol=1e-9)
    my_nms <- fitModels(counts, wtsmat, lapply(1:nrow(counts), function(i) emptyModel), ucs=ucs, tol=1e-9, type="dep")
    expect_equal(my_nms, their_nms, tol=1e-7)
    #with many threads
    my_nms2 <- fitModels(counts, wtsmat, lapply(1:nrow(counts), function(i) emptyModel), ucs=ucs, tol=1e-9, type="dep", nthreads=20)
    expect_equal(my_nms, my_nms2, tol=1e-7)
    
    #try out the nbtype="pois" option
    my_poiss <- fitModels(counts, wtsmat, lapply(1:nrow(counts), function(i) emptyModel), ucs=ucs, tol=1e-9, type="pois")
    their_nms <- lapply(1:nrow(counts), function(i) their_fitNM(counts, wtsmat[i,], tol=1e-9))
    their_nms <- lapply(their_nms, function(nm) {nm$r <- Inf; nm})
    expect_equal(my_poiss, their_nms)
    
    #try out the nbtype="nofit" option
    r <- 53.1231
    fixedModel <- list(mu=-1, r=r, ps=-1)
    my_nms <- fitModels(counts, wtsmat, lapply(1:nrow(counts), function(i) fixedModel), ucs=ucs, tol=1e-9, type="nofit")
    their_nms <- lapply(1:nrow(counts), function(i) their_fitNM(counts, wtsmat[i,], tol=1e-9))
    their_nms <- lapply(their_nms, function(nm) {nm$r <- r; nm})
    expect_equal(my_nms, their_nms)
    
    #make sure that fitting models to vectors of 0s doesn't not produce NaNs
    zeros <- matrix(0, nrow=10, ncol=10)
    ones <- matrix(1, nrow=1, ncol=10)
    for (type in c("dep", "indep", "pois")){
        mod <- fitModels(zeros, ones, list(emptyModel), mapToUnique(colSums(zeros)), type=type)
        expect_false(any(is.na(unlist(mod))))
    }
    #see if a matrix with one column gives what you expect
    mat <- matrix(1:10, ncol=1)
    for (type in c("dep", "indep", "pois")){
        mod <- fitModels(mat, matrix(1, ncol=1), list(emptyModel), mapToUnique(colSums(mat)), type=type)[[1]]
        expect_equal(mod$ps, mat[,1]/sum(mat[,1]))
        expect_equal(mod$mu, sum(mat))
    }
})


