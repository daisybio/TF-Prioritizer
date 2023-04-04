#include <Rcpp.h>
#include "array.cpp"

//the purpose of these functions is explained at the bottom of the document

inline int validateBreaks(std::vector<int>& breaks){
    int nthreads = breaks.size()-1;
    if (nthreads <= 0) Rcpp::stop("invalid breaks vector");
    return nthreads;
}

inline void validateJobs(Vec<double> jobSize){
    for (int i = 0; i < jobSize.len; ++i) if (jobSize[i] < 0) Rcpp::stop("negative job size");
}

//compute the makespan
inline double getMakespan(Vec<double> jobSize, std::vector<int>& breaks){
    int nthreads = validateBreaks(breaks);
    if (breaks[0] != 0 || breaks[nthreads] != jobSize.len) {
        Rcpp::stop("invalid breaks"); }
    double makespan = 0;
    for (int i = 1, j = 0; i <= nthreads; ++i){
        double groupSize = 0;
        while (j < breaks[i]){
            groupSize += jobSize[j++];
        }
        if (groupSize > makespan) makespan = groupSize;
    }
    return makespan;
}

//splits the jobs equally, without checking how big they are
inline void scheduleNaive(Vec<double> jobSize, std::vector<int>& breaks){
    validateJobs(jobSize);
    int nthreads = validateBreaks(breaks);
    //compute breaks
    double avgJobPerThread = jobSize.len/((double) nthreads);
    breaks[0] = 0;
    for (int i = 1; i <= nthreads; ++i){
        breaks[i] = round(i*avgJobPerThread);
    }
}

//tries to make the size of each partition as similar as possible to 
//the size of equal partitions
inline void scheduleGreedy(Vec<double> jobSize, std::vector<int>& breaks){
    validateJobs(jobSize);
    int njobs = jobSize.len;
    int nthreads = validateBreaks(breaks);
    //compute the size of a perfect partitioning
    double avgGroup = 0; 
    for (int i = 0; i < jobSize.len; ++i) avgGroup += jobSize[i];
    avgGroup /= nthreads;
    //compute breaks
    breaks[0] = 0;
    double cumSize = 0;
    for (int i = 1; i <= nthreads; ++i){
        double targetCumSize = avgGroup*i;
        int j = breaks[i-1];
        while (j < njobs && cumSize + jobSize[j] <= targetCumSize){
            cumSize += jobSize[j++];
        }
        //choose between before j and after j so as to minimize the 
        //absolute difference with the target
        if (j < njobs &&  targetCumSize - cumSize >= cumSize + jobSize[j] - targetCumSize){
            cumSize += jobSize[j++];
        }
        breaks[i] = j;
    }
}

//exact solution to the optimization problem using dynamic programming
//it could take time m*n*n, where n = jobSize.len and m = breaks.size(),
//even though in practice it will take much less
inline void scheduleOptimal(Vec<double> jobSize, std::vector<int>& breaks){
    validateJobs(jobSize);
    int njobs = jobSize.len;
    int nthreads = validateBreaks(breaks);
    //make sure the breaks vector starts out with all 0s
    for (int i = 0; i <= nthreads; ++i) breaks[i] = 0;
    //get rid of the case where jobSize is empty
    if (njobs == 0) return;
    //get the cumulative job sizes
    std::vector<double> cumJobSize(njobs+1); 
    cumJobSize[0] = 0;
    double acc = 0;
    for (int j = 1; j <= njobs; ++j){
        acc += jobSize[j-1];
        cumJobSize[j] = acc;
    }
    //element of the DP table, the first element is the score, the second
    //element is used for backtracking
    typedef std::pair<double, int> cell;
    //dp table, columns are threads, rows are jobs
    std::vector<cell> DP_storage(njobs*nthreads);
    Mat<cell> DP = asMat(DP_storage, nthreads);
    //use i or t for threads, j or k for jobs
    //DP(j,i): the optimal solution to the problem of scheduling
    //the first (j+1) jobs using (i+1) threads
    //DP(j,i).first: score of the solution (makespan)
    //DP(j,i).second: all jobs from DP(j,i).second (excluded) to j (included)
    //are done by the (i+1)-th thread. -1 is a valid value for DP(j,i).second, 
    //it means that all jobs are done by the last thread.
    //first iteration, first thread does the first job
    DP(0,0).first = jobSize[0]; DP(0,0).second = -1; //no more backtracking
    for (int i = 1; i < nthreads; ++i) {
        cell& c = DP(0, i);
        c.first = jobSize[0]; c.second = 0; //backtrack
    }
    //all other iterations
    for (int j = 1; j < njobs; ++j){
        //one thread does everything
        DP(j, 0).first = cumJobSize[j+1]; DP(j, 0).second = -1;
        for (int t = 1; t < nthreads; ++t){
            cell& c = DP(j, t);
            c.first = DP(j,t-1).first; 
            c.second = j;
            double thisThreadMakespan = 0;
            //here you could do binary search instead of trying all ks, 
            //because thisThreadMakespan always increases and 
            //otherThreadsMakespan always decreases, and you need to find the
            //crossing point. This would make the algorithm
            //n*m*log(n), but I don't want to debug that.
            for (int k = j; k >= 1; --k){
                thisThreadMakespan += jobSize[k];
                double otherThreadsMakespan = DP(k-1, t-1).first;
                double makeSpan = std::max(thisThreadMakespan, otherThreadsMakespan);
                if (makeSpan <= c.first){
                    c.first = makeSpan;
                    c.second = k-1;
                } else break;
            }
        }
    }
    //do backtracking
    int lastBreak = njobs; int t = nthreads;
    breaks[t] = lastBreak;
    while (t >= 0 && lastBreak > 0){
        cell& c = DP(lastBreak-1, t-1);
        lastBreak = c.second + 1;
        t -= 1;
        breaks[t] = lastBreak;
    }
}

//decides which jobs each thread gets. Each thread gets a set of consecutive jobs
//so that the makespan is minimized. The result is a vector of length 'nthreads+1',
//with the "breaks", such that the first element is 0 and the last is 'jobSize.len'.
inline std::vector<int> scheduleJobs(Vec<double> jobSize, int nthreads){
    double timeEstimate = nthreads*((double)jobSize.len)*jobSize.len;
    std::vector<int> breaks(nthreads+1);
    if (timeEstimate < 1e9){
        //use the exact algorithm
        scheduleOptimal(jobSize, breaks);
        return breaks;
    }
    //take the best of two fast heuristics
    scheduleNaive(jobSize, breaks);
    double t1 = getMakespan(jobSize, breaks);
    scheduleGreedy(jobSize, breaks);
    double t2 = getMakespan(jobSize, breaks);
    if (t1 < t2) scheduleNaive(jobSize, breaks);

    return breaks;
}
