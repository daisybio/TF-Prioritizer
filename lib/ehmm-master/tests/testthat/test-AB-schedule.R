context("AB: job scheduling works")

validateSchedule <- function(breaks, njobs, nthreads){
    if (length(breaks) != nthreads+1) stop("invalid breaks' length")
    if (breaks[1] != 0 || breaks[nthreads+1] != njobs) stop("invalid breaks endpoints")
    if (any(diff(breaks) < 0)) stop("breaks must be non-decreasing!")
}

#tiny wrapper around testSchedule which checks if the result makes sense
sched <- function(jobs, nthreads, algo){
    res <- testSchedule(jobs, nthreads, algo)
    validateSchedule(res$breaks, length(jobs), nthreads)
    res
}

test_that("job scheduling",{
    #test special cases
    #no jobs
    jobs <- integer(0)
    for (algo in 0:2){
        for (nthreads in c(1:10)){
            expect_true(sched(jobs, nthreads, algo)$makespan==0)
        }
    }
    #other special cases
    for (jobs in list(10, rep(0, 10), c(10, rep(0, 10)))){
        for (algo in 0:2){
            for (nthreads in c(1:10)){
                expect_true(sched(jobs, nthreads, algo)$makespan==jobs[1])
            }
        }
    }
    #another special case
    jobs <- rep(1/4, 8)
    for (algo in 0:2){
        expect_equal(sched(jobs, 4, algo)$makespan, .5)
    }
    
    #test the optimal schedule
    jobs <- rep(0, 100)
    jobs[c(10, 20, 21, 22, 40, 99)] <- 10
    expect_true(sched(jobs, 3, 2)$makespan == 20)
    expect_true(sched(jobs, 6, 2)$makespan == 10)
    #test on random inputs
    for (i in 1:50){
        jobs <- rnbinom(100, mu=50, size=1)
        if (i %% 2 == 0) jobs <- sort(jobs)
        for (nthreads in c(1, 2, 10)){
            spans <- sapply(0:2, function(algo) sched(jobs, 10, algo)$makespan)
            expect_true(spans[3] <= min(spans[1:2]))
        }
    }
    
    #test refineSplits
    rs <- refineSplits(100, 10)
    expect_equal(rs$newlens, rep(10,10))
    expect_equal(rs$origstarts, 1)
    for (i in 1:3){
    for (nseq in c(1, 10, 100)){
        for (nthreads in c(1, 10)){
            slens <- rnbinom(nseq, mu=1000, size=1)
            rs <- refineSplits(slens, nthreads)
            expect_equal(sum(rs$newlens), sum(slens))
            starts <- cumsum(c(0, slens[-length(slens)]))
            newstarts <- cumsum(c(0, rs$newlens[-length(rs$newlens)]))
            expect_true(setequal(rs$origstarts, match(starts, newstarts)))
        }
    }}
})
