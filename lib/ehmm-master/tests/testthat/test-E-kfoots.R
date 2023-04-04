context("E: The main functions run as expected")
source("utils.R")

#the r parameters in the model object can vary considerably when it
#diverges to infinity, but this is no meaningful difference.
#Better to compare mean and variance (things that can be measured)
convertModels <- function(fit){
    fit$models <- lapply(fit$models, function(model){
        res <- model
        res$var <- res$mu + (res$mu^2)/res$r
        res$r <- NULL
    })
    fit
}

test_that("kfoots works", {
    counts <- exampleData(5000, indep=T)
    
    ms <- c(1, nrow(counts))
    ks <- list(1,10)
    inits <- list("rnd", "count", "pca")
    nbtypes <- list("dep", "indep", "pois")
    maxiters <- list(10)
    nthreadss <- list(1, 10)
    
    res <- list()
    mcres <- list()
    for (m in ms){
        newcounts <- counts[1:m,,drop=F]
        for (k in ks){
            for (init in inits){
                for (nbtype in nbtypes){
                    for (maxiter in maxiters){
                        for (nthreads in nthreadss){
                            label <- paste0("m_", m, ",k_", k, ",init_", init, ",nbtype_", nbtype)
                            foot <- kfoots(newcounts, k, maxiter=maxiter, nbtype=nbtype, init=init, nthreads=nthreads, framework="MM", verbose=F)
                            models <- foot$models
                            expect_true(modelsAreOk(models, k, nrow(newcounts), nbtype))
                            mix_coeff <- foot$mix_coeff
                            expect_equal(length(mix_coeff), k)
                            expect_equal(sum(mix_coeff), 1)
                            expect_true(foot$loglik <= 0)
                            expect_true(is.logical(foot$converged))
                            llhistory <- foot$llhistory
                            expect_true(length(llhistory) <= maxiter)
                            expect_true(all(llhistory <= 0))
                            posteriors <- foot$posteriors
                            expect_equal(dim(posteriors), c(k, ncol(newcounts)))
                            expect_true(all(posteriors >= 0 & posteriors <= 1))
                            expect_equal(colSums(posteriors), rep(1, ncol(newcounts)))
                            clusters <- foot$clusters
                            expect_equal(length(clusters), ncol(newcounts))
                            expect_true(all(clusters>=1&clusters<=k))
                            
                            if (init != "rnd" && m > 1){
                                if (nthreads == nthreadss[[1]]) res[[label]] <- convertModels(foot)
                                else if (nthreads == nthreadss[[2]]) mcres[[label]] <- convertModels(foot)
                            }
                        }
                    }
                }
            }
        }
    }
    expect_equal(res, mcres, 1e-6)
})



test_that("hmmfoots works", {
    counts <- exampleHMMData(c(1000, 2500, 1500))
    
    ms <- c(1, nrow(counts))
    ks <- list(1,10)
    seqlenss <- list(5000, c(1, 4999), c(1000, 2500, 1500))
    inits <- list("rnd", "count", "pca")
    nbtypes <- list("dep", "indep", "pois")
    maxiters <- list(10)
    nthreadss <- list(1, 10)
    
    res <- list()
    mcres <- list()
    for (m in ms){
        newcounts <- counts[1:m,,drop=F]
        for (k in ks){
            for (seqlens in seqlenss){
                for (init in inits){
                    for (nbtype in nbtypes){
                        for (maxiter in maxiters){
                            for (nthreads in nthreadss){
                                for (split4speed in c(TRUE, FALSE)){
                                    label <- paste0("m_", m, ",k_", k, ",init_", init, ",nbtype_", nbtype, ",slens_", paste(seqlens, collapse=";"))
                                    hmm <- kfoots(newcounts, k, maxiter=maxiter, nbtype=nbtype, init=init, seqlens=seqlens, nthreads=nthreads, framework="HMM", verbose=F, split4speed=split4speed)
                                    models <- hmm$models
                                    expect_true(modelsAreOk(models, k, nrow(newcounts), nbtype))
                                    trans <- hmm$trans
                                    expect_equal(dim(trans), c(k,k))
                                    expect_equal(rowSums(trans), rep(1, k))
                                    initP <- hmm$initP
                                    expect_equal(dim(initP), c(k, length(seqlens)))
                                    expect_equal(colSums(initP), rep(1, length(seqlens)))
                                    expect_true(hmm$loglik <= 0)
                                    expect_true(is.logical(hmm$converged))
                                    llhistory <- hmm$llhistory
                                    expect_true(length(llhistory) <= maxiter)
                                    expect_true(all(llhistory <= 0))
                                    posteriors <- hmm$posteriors
                                    expect_equal(dim(posteriors), c(k, ncol(newcounts)))
                                    expect_true(all(posteriors >= 0 & posteriors <= 1))
                                    expect_equal(colSums(posteriors), rep(1, ncol(newcounts)))
                                    clusters <- hmm$clusters
                                    expect_equal(length(clusters), ncol(newcounts))
                                    expect_true(all(clusters>=1&clusters<=k))
                                    viterbi <- hmm$viterbi
                                    expect_equal(length(viterbi$vpath), ncol(newcounts))
                                    expect_true(all(viterbi$vpath >=1 & viterbi$vpath <=k))
                                    eq <- all.equal(viterbi$vllik, hmm$loglik)
                                    if (eq != T) expect_true(viterbi$vllik < hmm$loglik)
                                    
                                    if (init != "rnd" && !split4speed  && m > 1){
                                        #we don't test that the viterbi paths coincide, because
                                        #there are many equivalent viterbi paths and the choice
                                        #depends a lot on the numerical fuzz
                                        hmm$viterbi <- NULL
                                        if (nthreads == nthreadss[[1]]) res[[label]] <- convertModels(hmm)
                                        else if (nthreads == nthreadss[[2]]) mcres[[label]] <- convertModels(hmm)
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    expect_equal(res, mcres, 1e-5)
    #bug <- list(res=res, mcres=mcres, counts=counts)
    #save(bug, file="/project/ale/home/data/kfoots_pkg/bug.Rdata")
})


test_that("maxiter=0 means no training", {
    counts <- exampleHMMData(c(1000, 2500, 1500))
    k <- 4
    init <- initAlgo(counts, k)
    trans <- t(sapply(1:k, function(i) init$mix_coeff))
    initP <- matrix(nrow=k, rep(init$mix_coeff, 1))
    models <- init$models
    mix_coeff <- init$mix_coeff

    fit <- kfoots(counts, models, "HMM", trans=trans, initP=initP, maxiter=0, verbose=F)
    expect_equal(models, fit$models)
    expect_equal(trans, fit$trans)
    expect_equal(initP, fit$initP)

    fit <- kfoots(counts, models, "MM", mix_coeff=mix_coeff, maxiter=0, verbose=F)
    expect_equal(models, fit$models)
    expect_equal(mix_coeff, fit$mix_coeff)
})

test_that("split4speed option runs", {
    seqlens <- c(1000, 2500, 1500)
    counts <- exampleHMMData(seqlens)
    for (slens in list(seqlens, sum(seqlens))){
        k <- 4
        init <- initAlgo(counts, k)
        trans <- t(sapply(1:k, function(i) init$mix_coeff))
        initP <- matrix(nrow=k, rep(init$mix_coeff, length(slens)))
        models <- init$models
        mix_coeff <- init$mix_coeff

        kfoots(counts, models, "HMM", trans=trans, initP=initP, maxiter=5,
               nthreads=5, split4speed=FALSE, seqlens=slens, verbose=F)

        expect_true(TRUE)

        kfoots(counts, models, "HMM", trans=trans, initP=initP, maxiter=5,
               nthreads=5, split4speed=TRUE, seqlens=slens, verbose=F)

        expect_true(TRUE)

        fit <- kfoots(counts, models, "HMM", trans=trans, initP=initP, maxiter=0,
                        nthreads=5, split4speed=TRUE, seqlens=slens, verbose=F)
        expect_equal(initP, fit$initP)
    }
})


    

    

