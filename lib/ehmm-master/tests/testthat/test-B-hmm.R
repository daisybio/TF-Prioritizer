context("B: HMM and mixture models")

#generate a random vector of length l
#such that the values (all non-negative) sum up to 1
rndPS <- function(l){
    ps <- runif(l)
    ps/sum(ps)
}

if (require("HMM")){
    #this is copy and pasted from the viterbi function 
    #in the package HMM. Only the last line is modified to return
    #also the loglikelihood of the viterbi path
    HMM__viterbi <- function(hmm, observation){
        hmm$transProbs[is.na(hmm$transProbs)] = 0
        hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
        nObservations = length(observation)
        nStates = length(hmm$States)
        v = array(NA, c(nStates, nObservations))
        dimnames(v) = list(states = hmm$States, index = 1:nObservations)
        for (state in hmm$States) {
                v[state, 1] = log(hmm$startProbs[state] * hmm$emissionProbs[state, 
                        observation[1]])
        }
        for (k in 2:nObservations) {
                for (state in hmm$States) {
                        maxi = NULL
                        for (previousState in hmm$States) {
                                temp = v[previousState, k - 1] + log(hmm$transProbs[previousState, 
                                    state])
                                maxi = max(maxi, temp)
                        }
                        v[state, k] = log(hmm$emissionProbs[state, observation[k]]) + 
                                maxi
                }
        }
        viterbiPath = rep(NA, nObservations)
        for (state in hmm$States) {
                if (max(v[, nObservations]) == v[state, nObservations]) {
                        viterbiPath[nObservations] = state
                        break
                }
        }
        for (k in (nObservations - 1):1) {
                for (state in hmm$States) {
                        if (max(v[, k] + log(hmm$transProbs[, viterbiPath[k + 
                                1]])) == v[state, k] + log(hmm$transProbs[state, 
                                viterbiPath[k + 1]])) {
                                viterbiPath[k] = state
                                break
                        }
                }
        }
        
        list(vpath=viterbiPath, vllik=max(v[,nObservations]))
    }
}

test_that("hmm functions work",{
    if (require("HMM")){
        nstat <- 10
        nsymb <- 10
        seqlen <- c(500, 100, 10)
        #make the HMM object
        States <- paste0("state_", 1:nstat)
        Symbols <- paste0("symbol_", 1:nsymb)
        transProbs <- t(sapply(rep(nstat, nstat), rndPS))
        emissionProbs <- t(sapply(rep(nsymb, nstat), rndPS))
        startProbs <- rndPS(nstat)
        hmm <- initHMM(States, Symbols, startProbs, transProbs, emissionProbs)
        #make the observations
        obs <- sample(nsymb, sum(seqlen), replace=T)
        #get their posteriors
        seqstart <- cumsum(c(1,seqlen))[1:length(seqlen)]
        postList <- lapply(1:length(seqlen), function(idx){
            s <- seqstart[idx]
            e <- s + seqlen[idx] - 1
            HMM:::posterior(hmm, obs[s:e])
        })
        their_post <- do.call(cbind, postList)
        #get their new init probabilities
        their_initP <- their_post[,seqstart]
        
        #get their new transition probabilities
        transList <- lapply(1:length(seqlen), function(idx){
            s <- seqstart[idx]
            e <- s + seqlen[idx] - 1
            HMM:::baumWelchRecursion(hmm, obs[s:e])$TransitionMatrix
        })
        their_exptrans <- prop.table(Reduce("+", transList),1)
        #get their total log-likelihood
        their_llik <- sum(sapply(1:length(seqlen), function(idx){
            s <- seqstart[idx]
            e <- s + seqlen[idx] - 1
            f <- HMM:::forward(hmm, obs[s:e])
            f_lastcol <- f[,ncol(f)]
            llik1 <- max(f_lastcol)
            llik2 <- log(sum(exp(f_lastcol-llik1)))
            llik1+llik2
        }))
        
        #get my posteriors, new initP and new trans
        llik <- log(t(sapply(1:nstat, function(s){ emissionProbs[s, obs] })))
        my_post <- matrix(0, nrow=nstat, ncol=sum(seqlen))
        initP <- matrix(rep(startProbs, length(seqlen)), ncol=length(seqlen))
        fb <- forward_backward(initP, transProbs, llik, seqlen, my_post)
        my_exptrans <- fb$new_trans
        my_llik <- fb$tot_llik
        my_initP <- fb$new_initP
        
        #compare
        dimnames(my_post) <- dimnames(their_post)
        expect_equal(my_post, their_post)
        dimnames(my_exptrans) <- dimnames(their_exptrans)
        expect_equal(my_exptrans, their_exptrans)
        dimnames(my_initP) <- dimnames(their_initP)
        expect_equal(my_initP, their_initP)
        expect_equal(my_llik, their_llik)
        
        #with more threads
        llik <- log(t(sapply(1:nstat, function(s){ emissionProbs[s, obs] })))
        fb <- forward_backward(initP, transProbs, llik, seqlen, my_post, nthreads=20)
        my_exptrans <- fb$new_trans
        my_llik <- fb$tot_llik
        my_initP <- fb$new_initP
        
        #compare
        dimnames(my_post) <- dimnames(their_post)
        expect_equal(my_post, their_post)
        dimnames(my_exptrans) <- dimnames(their_exptrans)
        expect_equal(my_exptrans, their_exptrans)
        dimnames(my_initP) <- dimnames(their_initP)
        expect_equal(my_initP, their_initP)
        expect_equal(my_llik, their_llik)
        
        #test the viterbi algorithm
        vitList <- lapply(1:length(seqlen), function(idx){
            s <- seqstart[idx]
            e <- s + seqlen[idx] - 1
            HMM__viterbi(hmm, obs[s:e])
        })
        their_viterbi <- list(vpath=c(), vllik=0)
        for (vit in vitList) {
            their_viterbi$vpath <- c(their_viterbi$vpath, vit$vpath)
            their_viterbi$vllik <- their_viterbi$vllik + vit$vllik
        }
        #convert characters to state numbers
        transf <- 1:nstat; names(transf) <- States
        their_viterbi$vpath <- transf[their_viterbi$vpath]; names(their_viterbi$vpath) <- NULL
        llik <- log(t(sapply(1:nstat, function(s){ emissionProbs[s, obs] })))
        my_viterbi <- viterbi(initP, transProbs, llik, seqlen)
        
        expect_equal(my_viterbi, their_viterbi)
    }
    
    #test that the underflow generates an error
    for (nthreads in c(1,4)){
        load(system.file(package="ehmm", "extdata/uflowdata.Rdata"))
        post <- uflowdata$post; initP <- uflowdata$initP; 
        trans <- uflowdata$trans; lliks <- uflowdata$lliks; 
        slens <- uflowdata$seqlens
        expect_error(
        forward_backward(posteriors=post, initP, trans, lliks, slens, 
            nthreads=nthreads))
    }
    
})

test_that("mixture models work", {
    nmod <- 10
    lliks <- matrix(rnorm(1e3*nmod), nrow=nmod)
    mix_coeff <- rndPS(nmod)
    their_post <- exp(log(mix_coeff) + lliks)
    cs <- colSums(their_post)
    their_llik <- sum(log(cs))
    their_post <- their_post/cs[col(their_post)]
    their_new_mix_coeff <- rowMeans(their_post)
    
    my_post <- matrix(0, nrow=nrow(lliks), ncol=ncol(lliks))
    l2p <- llik2posteriors(lliks, mix_coeff, my_post)
    my_llik <- l2p$tot
    my_new_mix_coeff <- l2p$new_mix_coeff
    expect_equal(my_post, their_post)
    expect_equal(my_llik, their_llik)
    expect_equal(my_new_mix_coeff, their_new_mix_coeff)
    #with more threads
    l2p <- llik2posteriors(lliks, mix_coeff, my_post, nthreads=20)
    my_llik <- l2p$tot
    my_new_mix_coeff <- l2p$new_mix_coeff
    expect_equal(my_post, their_post)
    expect_equal(my_llik, their_llik)
    expect_equal(my_new_mix_coeff, their_new_mix_coeff)
    
})
