#include "array.cpp"
#include <algorithm> 
#include <unordered_map>
#include <Rcpp.h>
#include <cmath>
#include "optim.cpp"
#ifdef SUPPORT_OPENMP
    #include <omp.h>
#endif


//it would take too long to organize these headers properly, but conceptually, you should use only the following functions:
//map2unique_core
//getMultinomConst_core
//fitNBs_core
//fitMultinoms_core
//lLikMat_core

/* MAP 2 UNIQUE ROUTINES */
struct Avatar {
    int count;
    int pos;
    Avatar(int _count, int _pos): count(_count), pos(_pos){}
};

struct Avatardouble {
    double count;
    double pos;
    Avatardouble(double _count, double _pos): count(_count), pos(_pos){}
};

static inline bool avatarSorter(const Avatar& a, const Avatar& b){
    return b.count > a.count;
}

static inline bool avatardoubleSorter(const Avatardouble& a, const Avatardouble& b){
    return b.count > a.count;
}

//values is the input variable, map and uvalues are the output variables
//map and values can also wrap the same pointer
static inline void map2unique_core(Vec<double> values, Vec<int> map, std::vector<int>& uvalues){
    //this function is not parallelized:
    //it should be done just once and the running time should be very low
    
    if (values.len <= 0){
        return;//otherwise the code breaks with vectors of length 0
    }
    
    
    //sort the counts keeping track of the original position
    Avatar empty(0,0);
    std::vector<Avatar> avatars(values.len, empty);
    for (int i = 0, e = values.len; i < e; ++i){
        avatars[i].count = values[i];
        avatars[i].pos = i;
        
    }
    std::sort(avatars.begin(), avatars.end(), avatarSorter);
    
    //fill in uniquevalues and map
    int lastVal = avatars[0].count;
    uvalues.push_back(lastVal);
    map[avatars[0].pos] = 0;
    
    for (int i = 1, e = avatars.size(); i < e; ++i){
        Avatar& a = avatars[i];
        if (a.count != lastVal) {
            lastVal = a.count;
            uvalues.push_back(lastVal);
        }
        map[a.pos] = uvalues.size()-1;
    }
}

static inline void map2unique_Double_core(Vec<double> values, Vec<double> map, std::vector<double>& uvalues){
    //this function is not parallelized:
    //it should be done just once and the running time should be very low
    if (values.len <= 0){
        return;//otherwise the code breaks with vectors of length 0
    }
    
    
    //sort the counts keeping track of the original position
    Avatardouble empty(0,0);
    std::vector<Avatardouble> avatars(values.len, empty);
    for (int i = 0, e = values.len; i < e; ++i){
        avatars[i].count = values[i];
        avatars[i].pos = i;
    }
    std::sort(avatars.begin(), avatars.end(), avatardoubleSorter);
    
    //fill in uniquevalues and map
    double lastVal = avatars[0].count;
    uvalues.push_back(lastVal);
    map[avatars[0].pos] = 0;
    
    for (int i = 1, e = avatars.size(); i < e; ++i){
        Avatardouble& a = avatars[i];
        if (a.count != lastVal) {
            lastVal = a.count;
            uvalues.push_back(lastVal);
        }
        map[a.pos] = uvalues.size()-1;
    }
}


/* MULTINOM CONST ROUTINES */

struct CachedLFact{
    std::unordered_map<int, double> cache;
    CachedLFact(double load_factor){
        cache.max_load_factor(load_factor);
    }
    
    inline double operator()(int n){
        if (n <= 1){ return 0; }
        double cachedValue = cache[n];
        if (cachedValue == 0){//value not present in cache
            cache[n] = cachedValue = Rf_lgammafn(n + 1);
        }
        return cachedValue;
    }
};

//it can be implemented specifically for SWMat. The calls to 
//the lgamma functions would be the same, but less calls to lfact
//and faster colSums
template<template <typename> class TMat>
static void getMultinomConst_core(TMat<double> counts, Vec<double> multinomConst, int nthreads){
    //compute the log multinomial for each column
    CachedLFact lfact(0.75);
    int ncol = counts.ncol;
    int nrow = counts.nrow;
    
    #pragma omp parallel for firstprivate(lfact) num_threads(std::max(1, nthreads))
    for (int col = 0; col < ncol; ++col){
        double* countsCol = counts.colptr(col);
        double tmp = 0;
        double colsum = 0;
        for (int row = 0; row < nrow; ++row){
            double c = countsCol[row];
            colsum += c;
            tmp -= lfact(c);
        }
        multinomConst[col] = tmp + lfact(colsum);
    }
}


/* FITTING ROUTINES */

//this will slow down a bit, but it's safe... hopefully they will fix it in R 3.1...
//the caller has the responsibility that mu and r are >= 0
//this thing can be broken...
#define lognbinom(c, mu, r) (std::isfinite(mu*r)?Rf_dnbinom_mu(c, r, mu, 1):Rf_dpois(c, mu, 1))

struct NMPreproc {
    Vec<double> uniqueCS;
    Vec<int> map;
    Vec<double> multinomConst;
    NMPreproc(){}
    NMPreproc(Vec<double> _uniqueCS, Vec<int> _map, Vec<double> _multinomConst){
        uniqueCS = _uniqueCS;
        map = _map;
        multinomConst = _multinomConst;
    }
};

//data for optimization function.
//It is always assumed that mu is also the mean count.
struct optimData {
    Vec<double> counts;
    Vec<double> posteriors;
    double mu;
    double spost;
    int nthreads;
    bool verbose;
};


//optim fun for Brent
static double fn1d(double logr, void* data){
    optimData* info = (optimData*) data;
    double mu = info->mu;
    double r = exp(logr);
    double* post = info->posteriors.ptr;
    double* counts = info->counts.ptr;
    int e = info->counts.len;
    int nthreads = info->nthreads;
    bool verbose = info->verbose;
    
    long double llik = 0;
    #pragma omp parallel for schedule(static) reduction(+:llik) num_threads(nthreads)
    for (int i = 0; i < e; ++i){
        if (post[i] > 0){
            llik += lognbinom(counts[i], mu, r)*post[i];
        }
    }
    if (verbose) {
        std::cerr << "[numeric optimization] nbinom param:"
                  << "  r=" << r
                  << ", mu=" << mu
                  << std::endl;
    }
    //std::cout << "fn1d(" << r << "):\t" << -llik << std::endl;
    return -llik;
}

//data for optimization function 2.
struct optimData2 {
    Vec<double> counts;
    Mat<double> posteriors;
    Vec<double> mus;
    bool verbose;
    int nthreads;
};

//optim fun for Brent , 1r function
static double fn1d_2(double logr, void* data){
    optimData2* info = (optimData2*) data;
    Vec<double> mus = info->mus;
    Mat<double> post = info->posteriors;
    double* counts = info->counts.ptr;
    int e = info->counts.len;
    int nmod = post.nrow;
    int nthreads = info->nthreads;
    bool verbose = info->verbose;
    
    double r = exp(logr);

    long double llik = 0;
    #pragma omp parallel for schedule(static) reduction(+:llik) num_threads(nthreads)
    for (int i = 0; i < e; ++i){
        double c = counts[i];
        double* postcol = post.colptr(i);
        for (int mod = 0; mod < nmod; ++mod){
            if (postcol[mod] > 0){
                llik += lognbinom(c, mus[mod], r)*postcol[mod];
            }
        }
    }
    if (verbose) {
        std::cerr << "[numeric optimization] nbinom param:"
                  << " r=" << r;
        std::cerr << ", mu=(";
        for (int i = 0; i < nmod; ++i) {
            if (i != 0) {
                std::cerr << ",";
            }
            std::cerr << mus[i];
        }
        std::cerr << ")" << std::endl;
    }
    return -llik;
}

static inline void meanAndVar(Vec<double> weights, Vec<double> counts, double* mean, double* var, double* wsum, int nthreads){
    //get weighted averages
    long double sum_c2p = 0;
    long double sum_cp = 0;
    long double sum_p = 0;
    double* c = counts.ptr;
    double* p = weights.ptr;
    int len = counts.len;
    #pragma omp parallel for reduction(+:sum_p,sum_cp,sum_c2p) num_threads(nthreads)
    for (int i = 0; i < len; ++i){
        double P = p[i], C = c[i];
        sum_p += P;
        sum_cp += P*C;
        sum_c2p += P*C*C;
    }
    *wsum = sum_p;
    *mean = sum_cp/sum_p;
    *var = sum_c2p/sum_p - (*mean)*(*mean);
}


//make sure both values are finite and different from each other
//takes the logarithm
static void validateAndLogBraketing(double* initR, double* guessR, double tol){
    //low sample variance -> poisson case
    if (*guessR < 0 || !std::isfinite(*guessR)){ 
        *guessR = DBL_MAX;
    }
    if (*initR < 0){
        *initR = *guessR;
    } else if (!std::isfinite(*initR)){
        //this is just to avoid to work with infinities, 
        //but it should give the same score as infinity
        *initR = DBL_MAX; 
    }
    //optimization is done in the log space to avoid boundary constraints
    *guessR = log(*guessR);
    *initR = log(*initR);
    //initial points too close to each other
    if (fabs(*guessR - *initR) < tol*fabs(*guessR + *initR)/2){
        *guessR = *initR*0.9;
    }
}



//you need the guarantee that the r that you get will be better or equal to initR
//if initR is acceptable (not negative and not infinity) it must be one of the two
//points used for brent
static inline void fitNB_core(Vec<double> counts, Vec<double> posteriors, double* mu, double* r, double initR, double tol=1e-8, bool verbose = false, int nthreads=1){
    //get the variance of the data
    double var, wsum;
    meanAndVar(posteriors, counts, mu, &var, &wsum, 1);
    if (var <= *mu){//underdispersion, ML maximum at r=Inf
        *r = std::numeric_limits<double>::infinity();
        return;
    }
    //the r that matches the sample variance
    double guessR = (*mu)*(*mu) / (var  -  (*mu));
    validateAndLogBraketing(&initR, &guessR, tol);
    
    optimData info;
    info.counts = counts;
    info.posteriors = posteriors;
    info.mu = *mu;
    info.spost = wsum;
    info.nthreads = nthreads;
    info.verbose = verbose;
    
    initR = brent_wrapper(initR, guessR, fn1d, &info, tol);
        
    *r = exp(initR);
}


template<template <typename> class TVec>
static inline void fitMeans_core(TVec<double> counts, Mat<double> posteriors, Vec<double> mus, int nthreads){
    int ncol = posteriors.ncol;
    int nmod = posteriors.nrow;
    if (mus.len != nmod || ncol != counts.len){
        throw std::invalid_argument("invalid parameters passed to fitNBs_core");
    }
    //fill the mus with zeros, just to be safe
    for (int i = 0; i < nmod; ++i) mus[i] = 0;
    std::vector<long double> tot_scp(nmod);
    std::vector<long double> tot_sp(nmod);
    
    #pragma omp parallel num_threads(nthreads)
    {
        std::vector<long double> scp(nmod); long double* SCP = scp.data();
        std::vector<long double> sp(nmod); long double* SP = sp.data();
        #pragma omp for schedule(static) nowait
        for (int i = 0; i < ncol; ++i){
            double* post = posteriors.colptr(i);
            double count = counts[i];
            for (int j = 0; j < nmod; ++j){
                SP[j] += post[j];
                SCP[j] += count*post[j];
            }
        }
        
        #pragma omp critical
        for (int j = 0; j < nmod; ++j){
            tot_scp[j] += scp[j];
            tot_sp[j] += sp[j];
        }
    }
    
    for (int j = 0; j < nmod; ++j){
        mus[j] = tot_scp[j]/tot_sp[j];
    }
}

static inline void fitNBs_core(Vec<double> values, Mat<double> tposteriors, Vec<double> mus, Vec<double> rs, double tol=1e-8, bool verbose = false, int nthreads=1){
    int nval = tposteriors.nrow;
    int nmod = tposteriors.ncol;
    if (mus.len != nmod || rs.len != nmod || nval != values.len){
        throw std::invalid_argument("invalid parameters passed to fitNBs_helper");
    }
    //for the fitNB function call, I didn't find anything better than
    //nested parallel regions... this happens when nmod < nthreads
    int nthreads_outer = nmod > nthreads? nthreads : nmod;
    int nthreads_inner = ceil(((double)nthreads)/nmod);
    
    #ifdef SUPPORT_OPENMP
        omp_set_nested(true);
        //make sure we're not using more than nthreads threads
        omp_set_dynamic(true);
        //fitting the negative binomials
    #endif
    
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads_outer)
    for (int mod = 0; mod < nmod; ++mod){
        //call the optimization procedure. This will use up to nthreads_inner threads
        fitNB_core(values, tposteriors.getCol(mod), &mus[mod], &rs[mod], rs[mod], tol=tol, verbose=verbose, nthreads=nthreads_inner);
    }
}

static inline void fitNBs_1r_core(Vec<double> values, Mat<double> posteriors, Vec<double> mus, double* r, double tol=1e-8, bool verbose = false, int nthreads=1){
    if (values.len != posteriors.ncol || posteriors.nrow != mus.len) {
        std::invalid_argument("invalid parameters passed to fitNBs_1r_helper");
    }
    //fit the means
    fitMeans_core(values, posteriors, mus, nthreads);
    //fit R
    //get some global posteriors (should be the same as the occurrences of each count)
    std::vector<double> totPost(values.len);
    colSums(posteriors, asVec(totPost), nthreads);
    //get the variance of the data
    double mean, var, wsum;
    meanAndVar(asVec(totPost), values, &mean, &var, &wsum, nthreads);
    //get two initial points
    //the initial r
    double initR = *r;
    //the r that matches the sample variance
    double guessR = mean*mean / (var  -  mean);
    validateAndLogBraketing(&initR, &guessR, tol);
    
    optimData2 info;
    info.counts = values;
    info.posteriors = posteriors;
    info.mus = mus;
    info.nthreads = nthreads;
    info.verbose = verbose;
    
    initR = brent_wrapper(initR, guessR, fn1d_2, &info, tol);
        
    *r = exp(initR);
}

static inline void colpost(Mat<double> cpost, Mat<double> post, Vec<int> map){
    int ncol = post.ncol;
    int nmod = post.nrow;
    if (cpost.nrow != post.nrow || post.ncol != map.len){
        throw std::invalid_argument("invalid parameters passed to colpost");
    }
    double* S = post.ptr;
    for (int c = 0; c < ncol; ++c, S += nmod){
        double* D = cpost.colptr(map[c]);
        for (int m = 0; m < nmod; ++m) D[m] += S[m];
    }
}

static inline void collapsePosteriors_core(Mat<double> cpost, Mat<double> post, NMPreproc& preproc, int nthreads=1){
    int ncol = post.ncol;
    int nmod = post.nrow;
    if (cpost.nrow != post.nrow || post.ncol != preproc.map.len || cpost.ncol != preproc.uniqueCS.len){
        throw std::invalid_argument("invalid parameters passed to collapsePosteriors_core");
    }
    Vec<int> map = preproc.map;
    Vec<double> values = preproc.uniqueCS;
    
    //make sure that cpost is zeros at the beginning
    memset(cpost.ptr, 0, values.len*nmod*sizeof(double));
    
    //heuristic to pick the best number of threads
    nthreads = std::max(1, nthreads);
    int threadlim = std::max(1, (int)(round(ncol/(double)(values.len))));
    int maxthreads = std::min(nthreads, threadlim);
    //std::cout << "using " << maxthreads << " threads " << std::endl;
    
    //no parallelization
    if (maxthreads <= 1) {colpost(cpost, post, map); return; }
    
    //allocate memory for each thread as a matrix
    std::vector<double> threads_cpost_mem(cpost.nrow*cpost.ncol*maxthreads);
    Mat<double> threads_cpost = asMat(threads_cpost_mem, maxthreads);
    
    //std::cout << "computing breaks" << std::endl;
    //decide which columns each thread works on
    std::vector<int> breaks(maxthreads + 1); double colperthread = ncol/((double)maxthreads);
    for (int i = 1; i <= maxthreads; ++i) breaks[i] = round(i*colperthread);
    breaks[0] = 0; breaks[maxthreads] = ncol; //shouldn't be necessary...
    //std::cout << "breaks:\n" << std::endl;
    //for (int i = 0; i < breaks.size(); ++i) std::cout << breaks[i] << ", ";
    //std::cout << std::endl;
    
    //gogo
    //std::cout << "main loop" << std::endl;
    #pragma omp parallel for num_threads(maxthreads)
    for (int t = 0; t < maxthreads; ++t){
        int start = breaks[t]; int end = breaks[t+1];
        Mat<double> thread_cpost(threads_cpost.colptr(t), cpost.nrow, cpost.ncol);
        colpost(thread_cpost, post.subsetCol(start, end), map.subset(start, end));
    }
    //std::cout << "sum up" << std::endl;
    //sum up
    for (int t = 0; t < maxthreads; ++t){
        double* thread_cpost = threads_cpost.colptr(t);
        for (int i = 0, e = cpost.nrow*cpost.ncol; i < e; ++i){
            cpost[i] += thread_cpost[i];
        }
    }
}

//ps matrix needs to be clean at the beginning
template<template <typename> class TMat>
static void fitMultinoms_core(TMat<double> counts, Mat<double> posteriors, Mat<double> ps, int nthreads){
    //std::cout << "\ncounts\n" << counts << "\nposteriors\n" << posteriors << "\nps\n" << ps << "\nnthreads\n" << nthreads;
    int nmod = posteriors.nrow;
    int nrow = counts.nrow;
    int ncol = counts.ncol;
    if (ncol <= 0 || nrow <= 0 || nmod <= 0 || posteriors.ncol != ncol || ps.nrow != nrow || ps.ncol != nmod){
        throw std::invalid_argument("invalid parameters passed to fitMultinoms_core");
    }
    
    #pragma omp parallel num_threads(nthreads)
    {
        std::vector<double> loc_ps_std(nrow*nmod, 0);
        Mat<double> loc_ps = asMat(loc_ps_std, nmod);
        if (nmod == 1){//separate case with one model for efficiency
            #pragma omp for schedule(static) nowait
            for (int col = 0; col < ncol; ++col){
                double post = posteriors[col];
                if (post > 0){
                    double* countsCol = counts.colptr(col);
                    for (int row = 0; row < nrow; ++row){
                        loc_ps[row] += countsCol[row]*post;
                    }
                }
            }
        } else {
            #pragma omp for schedule(static) nowait
            for (int col = 0; col < ncol; ++col){
                //This you could probably do better with BLAS...
                double* postCol = posteriors.colptr(col);
                double* countsCol = counts.colptr(col);
                for (int mod = 0; mod < nmod; ++mod, ++postCol){
                    double post = *postCol;
                    if (post > 0){
                        double* loc_psCol = loc_ps.colptr(mod);
                        for (int row = 0; row < nrow; ++row){
                            loc_psCol[row] += countsCol[row]*post;
                        }
                    }
                }
            }
        }
        
        #pragma omp critical
        {
            for (int i = 0, e = nrow*nmod; i < e; ++i){
                ps[i] += loc_ps[i];
            }
        }
        #pragma omp barrier
        //ps contains now the non-normalized optimal ps parameter, normalize!
        //parallelization is probably overkill, but can it be worse than without?
        
        #pragma omp for schedule(static)
        for (int mod = 0; mod < nmod; ++mod){
            double* psCol = ps.colptr(mod);
            double sum = 0;
            for (int row = 0; row < nrow; ++row){sum += psCol[row];}
            if (sum > 0){
                for (int row = 0; row < nrow; ++row){psCol[row] /= sum;}
            } else {
                //no training data... use uniform distr.
                double unif = 1.0/nrow;
                for (int row = 0; row < nrow; ++row){psCol[row] = unif;}
            }
        }
    }
}

/* EXPERIMENTAL LOG LIKELIHOOD ROUTINES 

static inline double dotprod(int* v1, double* v2, int len){
    int startup = len % 5;
    double sum = 0; int i = 0;
    for (; i < startup; ++i){ sum += v1[i]*v2[i]; }
    for (; i < len; i += 5){
        sum += v1[i]*v2[i] + v1[i+1]*v2[i+1] + v1[i+2]*v2[i+2] + v1[i+3]*v2[i+3] + v1[i+4]*v2[i+4]; 
    }
    return sum;
}
static void lLik_nbinom(double mu, double r, Vec<int> uniqueCS, Vec<double> tmpNB, int nthreads){
    if (tmpNB.len != uniqueCS.len) throw std::invalid_argument("the preprocessed data were not computed on the same count matrix");
    int nUCS = uniqueCS.len;
    #pragma omp parallel for schedule(static) num_threads(nthreads)
    for (int c = 0; c < nUCS; ++c){
        tmpNB[c] = lognbinom(uniqueCS[c], mu, r);
    }
}

template<template <typename> class TMat>
static void lLik_multinom(TMat<double> counts, Vec<double> ps, Vec<double> llik, Vec<int> mapp, Vec<double> tmpNB, Vec<double> mconst, int nthreads){
    if (ps.len != counts.nrow) throw std::invalid_argument("incoherent models provided");
    if (counts.ncol != mapp.len) throw std::invalid_argument("the preprocessed data were not computed on the same count matrix");
    
    int ncol = counts.ncol;
    int nrow = counts.nrow;
    int* map = mapp.ptr;
    double* llikptr = llik.ptr;
    double* mtnmConst = mconst.ptr;
    std::vector<double> logPsSTD(nrow);
    Vec<double> logPs = asVec(logPsSTD);
    
    //pre-compute all the log(p)
    for (int i = 0; i < nrow; ++i){ logPs[i] = log(ps[i]);}
    
    #pragma omp parallel for schedule(static) num_threads(nthreads)
    for (int col = 0; col < ncol; ++col){
        llikptr[col] = tmpNB[map[col]] + mtnmConst[col] + dotprod(counts.colptr(col), logPs.ptr, nrow);
    }
}
 */



/* LOG LIKELIHOOD ROUTINE  */

template<template <typename> class TMat>
static void lLikMat_core_lognormal(TMat<double> counts, Mat<double> mus, Mat<double> sigmasqs, Mat<double> llik, NMPreproc& preproc, Mat<double> tmpNB, int nthreads){
  if (mus.ncol != sigmasqs.ncol || mus.nrow != counts.nrow) throw std::invalid_argument("incoherent models provided");

  nthreads = std::max(1, nthreads);
  int nmodels = mus.ncol;
  int ncol = counts.ncol;
  int nrow = counts.nrow;
  int nUCS = preproc.uniqueCS.len;
  double* uniqueCS = preproc.uniqueCS.ptr;
  int* map = preproc.map.ptr;

  /* MAIN LOOP */
  if (nmodels==1){
    double mu = mus[0], sigmasq = sigmasqs[0]; double* llikptr = llik.ptr;
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp for schedule(static)
      for (int col = 0; col < ncol; ++col){
	double* currcount = counts.colptr(col);
	for (int row = 0; row < nrow; ++row){
	  double count = currcount[row];
	  double logcount = log(count);
	  // calculate probability of logcount given LN(mu,sigmasq)
	  llikptr[row,col] = -log(std::sqrt(2 * M_PI * sigmasq) * count) - (std::pow((logcount - mu),2) / (2 * sigmasq));
	}
      }
    }
  } else{
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp for schedule(static)
      for (int col = 0; col < ncol; ++col){
	double* currcount = counts.colptr(col);
	double* llikCol = llik.colptr(col);
	
	for (int mod = 0; mod < nmodels; ++mod){
	  llikCol[mod] = 0;
	  for (int row = 0; row < nrow; ++row){
	    double mu = mus.colptr(mod)[row], sigmasq = sigmasqs.colptr(mod)[row];
	    double count = currcount[row];
	    double logcount = log(count);
	    // calculate probability of logcount given N(mu,sigmasq)
	    // if sigmasq == 0: if logcount == mu; llik = inf; else llik = 0

		if (sigmasq == 0){
	      if (logcount == mu) {
			//llikCol[mod] += std::numeric_limits<double>::infinity();
			llikCol[mod] += 10000;
	      }
	      else{
			llikCol[mod] += 0;
	      }
	    }
	    else{
		  llikCol[mod] += -log(std::sqrt(2 * M_PI * sigmasq) * count) - (std::pow((logcount - mu),2) / (2 * sigmasq));
	    }
	  }
	  //std::cout << "mod,col: " << mod << "," << col << ", llik: " << llikCol[mod] << std::endl;
	}
      }
    }
  }
}

template<template <typename> class TMat>
static void lLikMat_core(TMat<double> counts, Vec<double> mus, Vec<double> rs, Mat<double> ps, Mat<double> llik, NMPreproc& preproc, Mat<double> tmpNB, int nthreads){
  if (rs.len != mus.len || mus.len != ps.ncol || ps.nrow != counts.nrow) throw std::invalid_argument("incoherent models provided");
  if (counts.ncol != preproc.map.len || counts.ncol != preproc.multinomConst.len) throw std::invalid_argument("the preprocessed data were not computed on the same count matrix");
  
  nthreads = std::max(1, nthreads);
  int nmodels = mus.len;
  int ncol = counts.ncol;
  int nrow = counts.nrow;
  int nUCS = preproc.uniqueCS.len;
  int logPsSize = nrow*nmodels;
  double* uniqueCS = preproc.uniqueCS.ptr;
  int* map = preproc.map.ptr;
  double* mtnmConst = preproc.multinomConst.ptr;
  std::vector<double> logPsSTD(nrow*nmodels);
  Mat<double> logPs = asMat(logPsSTD, nmodels);
  
  //pre-compute all the log(p)
  for (int i = 0; i < logPsSize; ++i){
    logPs[i] = log(ps[i]);
  }
  
  /* MAIN LOOP */
  if (nmodels==1){ //separate case where there is only one model, for efficiency
    double mu = mus[0], r = rs[0]; double* llikptr = llik.ptr; 
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp for schedule(static)
      for (int c = 0; c < nUCS; ++c){
	tmpNB[c] = lognbinom(uniqueCS[c], mu, r);
      }
      
      //contribution of the multinomial
      #pragma omp for schedule(static)
      for (int col = 0; col < ncol; ++col){
	double tmp = tmpNB[map[col]] + mtnmConst[col];
	double* currcount = counts.colptr(col);
	const double* currlogp = logPs.ptr;
	for (int row = 0; row < nrow; ++row, ++currcount, ++currlogp){
	  double c = *currcount;
	  if (c != 0){
	    tmp += c*(*currlogp);
	  }
	}
	llikptr[col] = tmp;
      }
    }
    
  } else { 
    #pragma omp parallel num_threads(nthreads)
    {
      //compute all the log neg binom on the unique column sums
      //Here the outer loop is on the models because I think it divides the
      //computation more evenly (even if we need to collapse): uniqueCS
      //is sorted...
      #pragma omp for schedule(static) collapse(2)
      for (int p = 0; p < nmodels; ++p){
	for (int c = 0; c < nUCS; ++c){
	  tmpNB(p,c) = lognbinom(uniqueCS[c], mus[p], rs[p]);
	}
      }
      
      //compute and add contribution of the multinomial 
      #pragma omp for schedule(static) 
      for (int col = 0; col < ncol; ++col){
	double* countsCol = counts.colptr(col);
	double* llikCol = llik.colptr(col);
	double mtnm = mtnmConst[col];
	double* tmpNBcol = tmpNB.colptr(map[col]);
        
	for (int mod = 0; mod < nmodels; ++mod){
	  double tmp = tmpNBcol[mod] + mtnm;
	  double* logp = logPs.colptr(mod);
	  for (int row = 0; row < nrow; ++row){
	    double c = countsCol[row];
	    if (c != 0){
	      tmp += c*logp[row];
	    }
	  }
	  llikCol[mod] = tmp;
	}
      }
    }
  }
}
