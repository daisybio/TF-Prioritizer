#include <Rcpp.h>
#include "core.cpp"
#include "schedule.cpp"
#include "logarithmetic.hpp"

template <typename T>
static inline bool allequal(const T& a, const T& b) {return a==b;}
template <typename T, typename... S>
static inline bool allequal(const T& a, const T& b, const T& c, const S&... d){ 
    return a==b && allequal(a, c, d...); 
}

using namespace Rcpp;
typedef IntegerVector::iterator iiter;

//temporary storage needed by fb_iter
struct FBtmp {
    int nstates;
    std::vector<double> mem;
    Vec<double> backward;
    Vec<double> new_backward;
    Mat<double> tmp;
    FBtmp(int ns) : nstates(ns), mem(ns*(ns+2)) {
        double* start = mem.data();
        backward = Vec<double>(start, nstates);
        new_backward = Vec<double>(start + nstates, nstates);
        tmp = Mat<double>(start + 2*nstates, nstates, nstates);
    }
};

//forward backward iteration for one sequence of observations
//posteriors: at the beginning it doesn't matter, at the end, the posterior probabilities
//eprobs: at the beginning and at the end, the emission probabilities
//eprobs and prosteriors have the same dimensions. 
//The rows are the states, the columns are the observations
//initP: at the beginning and at the end the initial probabilities
//new_initP: at the beginning it doesn't matter, at the end the fitted initial probabilitites
//trans: at the beginning and end, the transition probabilitites
//new_trans: at the beginning, it doesn't matter, at the end the new transition
//probabilities ARE ADDED to the initial values, unnormalized, to allow accumulation
//llik: at the beginning, it doesn't matter, at the end the log likelihood of
//this sequence is ADDED
//storage: temporary storage needed by the function. You can provided this arg
//for efficiency, but you don't need to.
//return codes:
//-2: everything fine
//-1: invalid dimensions of input arrays
//n>=0: at column n there was an underflow error
const int FB_OK = -2;
const int FB_WRONG_DIM = -1;
inline int fb_iter(Mat<double> eprobs, Mat<double> posteriors, Mat<double> alpha,
				   Vec<double> initP, Vec<double> new_initP, 
				   Mat<double> trans, Mat<double> new_trans, long double& llik,
				   FBtmp& storage){
  
  int nobs    = eprobs.ncol;
  int nstates = eprobs.nrow;
  if (nobs != posteriors.ncol || 
	  !allequal(nstates, posteriors.nrow, initP.len, new_initP.len, 
				trans.ncol, trans.nrow, new_trans.nrow, new_trans.ncol, 
                    storage.nstates)){
	return FB_WRONG_DIM; //invalid dimensions of input arrays
  }
  if (nobs == 0) return FB_OK;
  /* FORWARD LOOP */
  /* first iteration is from fictitious start state */
  {
	double cf = -std::numeric_limits<double>::infinity();//scaling factor
	double* emissprob = eprobs.colptr(0);
	double* forward   = posteriors.colptr(0);
	double* alpha_colptr = alpha.colptr(0);
	for (int r = 0; r < nstates; ++r){
	  double p = emissprob[r] + log(initP[r]);
	  forward[r] = p;
	  alpha_colptr[r] = p;
	  cf = logadd(cf, p);
	}
	if (cf==-std::numeric_limits<double>::infinity()) return 0;  //underflow error
	for (int r = 0; r < nstates; ++r){
	  forward[r] = forward[r] - cf;
	}
	llik += cf;
  }
  /* all other iterations */
  for (int i = 0 + 1; i < nobs; ++i){
	double cf = -std::numeric_limits<double>::infinity();//scaling factor
	double*    emissprob = eprobs.colptr(i);
	double*      forward = posteriors.colptr(i);
	double* last_forward = posteriors.colptr(i-1);
	double* alpha_colptr = alpha.colptr(i);
	double* last_alpha_colptr = alpha.colptr(i-1);
	for (int t = 0; t < nstates; ++t){
	  double* transcol = trans.colptr(t);
	  double acc = -std::numeric_limits<double>::infinity(); double acc_unnormalized = -std::numeric_limits<double>::infinity();
	  for (int s = 0; s < nstates; ++s){
		acc = logadd(acc, last_forward[s] + log(transcol[s]));
		acc_unnormalized = logadd(acc_unnormalized, last_alpha_colptr[s] + log(transcol[s]));
	  }
	  acc += emissprob[t];
	  acc_unnormalized += emissprob[t];
	  forward[t] = acc;
	  alpha_colptr[t] = acc_unnormalized;
	  cf = logadd(cf, acc);
	}
	if (cf==-std::numeric_limits<double>::infinity()) return i; //underflow error
	for (int t = 0; t < nstates; ++t){
	  forward[t] = forward[t] - cf;
	}
	llik += cf;
  }  
  
  /* BACKWARD LOOP */
  /* we don't keep the backward matrix, only a 'backward' column */
  /* this gets replaced by 'new_backward' at each iteration */
  /* first iteration set backward to 1/k, 
   * last column of posteriors is already ok */
  
  Vec<double>     backward = storage.backward;
  Vec<double> new_backward = storage.new_backward;
  Mat<double>          tmp = storage.tmp;

  /* Set backward and posteriors at last position */
  double* last_posterior = posteriors.colptr(nobs-1);
  double norm = -std::numeric_limits<double>::infinity();
  for (int r = 0; r < nstates; ++r){
	//backward[r] = log(1.0/nstates);
	backward[r] = log(1.0);
	norm = logadd(norm, backward[r] + last_posterior[r]);
  }
  for (int r = 0; r < nstates; ++r){
	last_posterior[r] = last_posterior[r] + backward[r] - norm;
  }

  /* all other iterations */
  if (nobs > 1) {
	for (int i = nobs-2; i >= 0; --i){
	  //std::cout << "\nposition: " << i << ":" << std::endl;
	  double* emissprob = eprobs.colptr(i+1);
	  double* posterior = posteriors.colptr(i);
	  double cf   = -std::numeric_limits<double>::infinity();
	  double norm = -std::numeric_limits<double>::infinity();
	  /* joint probabilities and backward vector */
	  for (int s = 0; s < nstates; ++s){
		//the forward variable is going to be overwritten with the posteriors
		double pc  = posterior[s];
		double acc = -std::numeric_limits<double>::infinity();
		
		for (int t = 0; t < nstates; ++t){
		  double p = log(trans(s, t)) + emissprob[t] + backward[t];
		  tmp(s, t) = pc + p;
		  acc = logadd(acc, p);
		}
		
		new_backward[s] = acc;
		cf = logadd(cf, acc);
	  }
	  //if (cf==-std::numeric_limits<double>::infinity()) return i; //underflow error
	  /* update backward vector */
	  for (int s = 0; s < nstates; ++s){
		//backward[s] = new_backward[s];
		backward[s] = new_backward[s] - cf;
		norm = logadd(norm, backward[s] + posterior[s]);
	  }
	  /* update transition probabilities */
	  for (int t = 0, e = nstates*nstates; t < e; ++t){
		//std::cout << t << ", " << new_trans[t] << ", " << tmp[t] << ", " << (norm+cf) << std::endl;
		//new_trans[t] = new_trans[t] + exp(tmp[t]-norm);
		new_trans[t] = new_trans[t] + exp(tmp[t] - (norm+cf));
	  }
	  /* get posteriors */
	  for (int s = 0; s < nstates; ++s){
		posterior[s] = posterior[s] + backward[s] - norm;
	  }
	}
  }
  /* exponentiate posterior */
  for (int i = 0; i < nobs; ++i){
	double* posterior = posteriors.colptr(i);
	for (int s = 0; s < nstates; ++s){
	  posterior[s] = exp(posterior[s]);
	}
  }
  /* set new_initP */
  double* posterior = posteriors.colptr(0);
  for (int r = 0; r < nstates; ++r){
	new_initP[r] = posterior[r];
  }
  
  return FB_OK;;
}

inline int fb_iter(Mat<double> eprobs, Mat<double> posteriors,
				   Mat<double> alpha,
                   Vec<double> initP, Vec<double> new_initP, 
                   Mat<double> trans, Mat<double> new_trans,
				   long double& llik){
  FBtmp storage(eprobs.nrow);
  return fb_iter(eprobs, posteriors, alpha, initP, new_initP, trans, new_trans, llik, storage);

}

inline void compute_escore_r(int r, int maxn, Vec<double> prior, Vec<double> escore, Mat<double> eprobs, Mat<double> alpha, Mat<double> trans_tilde, Vec<double> alpha_tilde, Vec<double> last_alpha_tilde){
  // this function calculates scores for enhancers of length n = {1,...,maxn} and logadds each score to the respective t = r+i (the middle of the enhancer)
  int nobs    = eprobs.ncol;
  int nstates = eprobs.nrow;

  for (int i = 0; i < nstates; i++) {
      alpha_tilde[i] = -std::numeric_limits<double>::infinity();
  }

  // bugfix for:
  // double* eprobs_colptr = eprobs.colptr(r+n-1);
  for (int n = 1; n <= maxn && r+n < eprobs.ncol; ++n){
	int i = ceil(n/2.) - 1;
	int t = r + i;

	// transition from background to enhancer, i.e. from alpha to alpha_tilde
	if (n == 1){
	  for (int y = 0; y < nstates; ++y){
		double* eprobs_colptr = eprobs.colptr(r);
		double* last_alpha_colptr = alpha.colptr(r-1);
		double* trans_colptr = trans_tilde.colptr(y);
		double p = -std::numeric_limits<double>::infinity();
		for (int x = 0; x < nstates; ++x){
		  p = logadd(p, last_alpha_colptr[x] + log(trans_colptr[x]));
		}
		alpha_tilde[y] = p + eprobs_colptr[y];
	  }
	}
	else {
	  // transition from enhancer to enhancer position, i.e. from alpha_tilde to alpha_tilde
	  memcpy(last_alpha_tilde.ptr, alpha_tilde.ptr, sizeof(double)*nstates); // copy alpha_tilde from last iteration to last_alpha_tilde

	  for (int y = 0; y < nstates; ++y){
		double* eprobs_colptr = eprobs.colptr(r+n-1);
		double* trans_colptr = trans_tilde.colptr(y);
		double p = -std::numeric_limits<double>::infinity();
		for (int x = 0; x < nstates; ++x){
		  p = logadd(p, last_alpha_tilde[x] + log(trans_colptr[x]));
		}
		alpha_tilde[y] = p + eprobs_colptr[y];
	  }
    }
	
	// for current enhancer length n, calculate score at current position r+n-1 and save it to t = r+i (middle of enhancer).
	// this is done by normalizing sum_y(alpha_tilde(y,r+n-1)) by sum_y(alpha(y,r+n-1)) (use logadd). values are in log-space, therefore subtract for normalization.
	// a respective prior is added to each log-probability, corresponding to the expected length of the enhancer elements.
	double* alpha_colptr = alpha.colptr(r+n-1);
	double alpha_sum = -std::numeric_limits<double>::infinity();
	double alpha_tilde_sum = -std::numeric_limits<double>::infinity();
	for (int y = 0; y < nstates; ++y){
	  alpha_tilde_sum = logadd(alpha_tilde_sum, alpha_tilde[y]);
	}
	for (int y = 0; y < nstates; ++y){
	  alpha_sum = logadd(alpha_sum, alpha_colptr[y]);
	}
	escore[t] = logadd(escore[t], prior[n-1] + alpha_tilde_sum - alpha_sum);
  }
}

inline void compute_escore_0(int maxn, Vec<double> prior, Vec<double> escore, Mat<double> eprobs, Mat<double> alpha, Mat<double> trans_tilde, Vec<double> initP_tilde, Vec<double> alpha_tilde, Vec<double> last_alpha_tilde){
  // this function calculates scores for enhancers of length n = {1,...,maxn} and logadds each score to the respective t = r+i (the middle of the enhancer)
  // special case r = 0: enhancers starting at the first position.
  int r = 0;
  int nobs    = eprobs.ncol;
  int nstates = eprobs.nrow;

  for (int i = 0; i < nstates; i++) {
      alpha_tilde[i] = -std::numeric_limits<double>::infinity();
  }

  for (int n = 1; n <= maxn && r+n < eprobs.ncol; ++n){
	int i = ceil(n/2.) - 1;
	int t = r + i;
	double* eprobs_colptr = eprobs.colptr(r+n-1);

	// transition from background to enhancer, i.e. from alpha to alpha_tilde
	if (n == 1){
	  for (int y = 0; y < nstates; ++y){
		alpha_tilde[y] = log(initP_tilde[y]) + eprobs_colptr[y];
	  }
	}
	else {
	  // transition from enhancer to enhancer position, i.e. from alpha_tilde to alpha_tilde
	  memcpy(last_alpha_tilde.ptr, alpha_tilde.ptr, sizeof(double)*nstates); // copy alpha_tilde from last iteration to last_alpha_tilde

	  for (int y = 0; y < nstates; ++y){
		double* trans_colptr = trans_tilde.colptr(y);
		double p = -std::numeric_limits<double>::infinity();
		for (int x = 0; x < nstates; ++x){
		  p = logadd(p, last_alpha_tilde[x] + log(trans_colptr[x]));
		}
		alpha_tilde[y] = p + eprobs_colptr[y];
	  }
	}

	// for current enhancer length n, calculate score at current position r+n-1 and save it to t = r+i (middle of enhancer).
	// this is done by normalizing sum_y(alpha_tilde(y,r+n-1)) by sum_y(alpha(y,r+n-1)) (use logadd). values are in log-space, therefore subtract for normalization.
	double* alpha_colptr = alpha.colptr(r+n-1);
	double alpha_sum = -std::numeric_limits<double>::infinity();
	double alpha_tilde_sum = -std::numeric_limits<double>::infinity();
	for (int y = 0; y < nstates; ++y){
	  alpha_tilde_sum = logadd(alpha_tilde_sum, alpha_tilde[y]);
	}
	for (int y = 0; y < nstates; ++y){
	  alpha_sum = logadd(alpha_sum, alpha_colptr[y]);
	}
	escore[t] = logadd(escore[t], prior[n-1] + alpha_tilde_sum - alpha_sum);
  }
}

inline void escore_core(Vec<double> escore, Mat<double> lliks, Mat<double> alpha, Mat<double> trans, Mat<double> initPs, IntegerVector enhancer_states, Vec<double> prior, Vec<int> seqlens, int maxn, int nthreads=1){
  // this function calculates scores for multiple enhancer lengths n, summing them up to a total enhancer score.
  // each enhancer score for length n is multiplied by a prior for enhancer length n.
  int nrow = lliks.nrow;
  int ncol = lliks.ncol;
  int nchunk = seqlens.len;

  //temporary objects
  std::vector<int> chunk_startsSTD(seqlens.len, 0);
  Vec<int> chunk_starts = asVec<int>(chunk_startsSTD);
  //get the start of each chunk
  for (int i = 0, acc = 0; i < nchunk; ++i){chunk_starts[i] = acc; acc += seqlens[i];}

  long double tot_llik = 0;

  //figure out how to assign the chromosomes to the threads
  //covert seqlens to double
  std::vector<double> jobSize(nchunk);
  for (int i = 0; i < nchunk; ++i) jobSize[i] = seqlens[i];
  //get the assignments
  std::vector<int> breaks = scheduleJobs(asVec(jobSize), nthreads);
  int gretcode = FB_OK; //return code from fb_iter, aggregated across threads
  int guflowchunk = -1; //chunk that caused the underflow
  
  // define initPs_tilde, a matrix with initial probabilities of each region (e.g. chromosome) as a column and its values for non-enhancer states set to zero.
  int nrow_initPs = initPs.nrow;
  int ncol_initPs = initPs.ncol;
  NumericMatrix tmpMat(nrow_initPs, ncol_initPs);
  Mat<double> initPs_tilde = asMat(tmpMat);
  for (int x = 0; x < nrow_initPs; ++x){
	if (std::find(enhancer_states.begin(), enhancer_states.end(), x) != enhancer_states.end()){
	  for (int y = 0; y < ncol_initPs; ++y){
		initPs_tilde(x,y) = initPs(x,y);
	  }
	}
	else{
	  for (int y = 0; y < ncol_initPs; ++y){
		initPs_tilde(x,y) = 0;
	  }
	}
  }

  // define trans_tilde, a transition matrix with transitions to non-enhancer states set to zero.
  NumericVector tmpVec(nrow, 0.);
  Vec<double> zeroVec = asVec(tmpVec);
  tmpMat = NumericMatrix(nrow, nrow);
  Mat<double> trans_tilde = asMat(tmpMat);
  for (int y = 0; y < nrow; ++y){
	if (std::find(enhancer_states.begin(), enhancer_states.end(), y) != enhancer_states.end()){
	  memcpy(trans_tilde.colptr(y), trans.colptr(y), sizeof(double)*nrow);
	}
	else{
	  memcpy(trans_tilde.colptr(y), zeroVec.ptr, sizeof(double)*nrow);
	}
  }

  #pragma omp parallel num_threads(nthreads)
  {
	NumericVector tmpVec1(nrow);
	Vec<double> alpha_tilde = asVec(tmpVec1);
	NumericVector tmpVec2(nrow);
	Vec<double> last_alpha_tilde = asVec(tmpVec2);
	// Do forward and backward loop for each chunk (defined by seqlens)                                                                                                                                   
	// Chunks might have very different lengths (that's why they have been scheduled).
    #pragma omp for schedule(static) nowait
	for (int thread = 0; thread < nthreads; ++thread){
	  for (int o = breaks[thread]; o < breaks[thread+1]; ++o){
		//o identifies the sequence/chunk
		int chunk_start = chunk_starts[o];
		int chunk_end =  chunk_start + seqlens[o];

		// first position
		int r = 0;
		compute_escore_0(maxn, prior,
						 escore.subset(chunk_start, chunk_end),
						 lliks.subsetCol(chunk_start, chunk_end),
						 alpha.subsetCol(chunk_start, chunk_end),
						 trans_tilde, initPs_tilde.getCol(o),
						 alpha_tilde,
						 last_alpha_tilde);

		// all other positions
        // bugfix:
		for (int r = 1; r < seqlens[o]; ++r){
		  maxn = std::min(maxn, ncol-r); // limit maximum enhancer towards the end
		  compute_escore_r(r, maxn, prior,
						   escore.subset(chunk_start, chunk_end),
						   lliks.subsetCol(chunk_start, chunk_end),
						   alpha.subsetCol(chunk_start, chunk_end),
						   trans_tilde,
                           alpha_tilde,
                           last_alpha_tilde);
		}
	  }
	}
  }
}

List enhancer_score(NumericMatrix initP, NumericMatrix trans, NumericMatrix lliks, IntegerVector seqlens,
					NumericMatrix alpha, IntegerVector enhancer_states, NumericVector escore, NumericVector prior, int nthreads=1){
    int nmod = initP.nrow();
    double totlen = Rcpp::sum(seqlens);
	int maxn = prior.size();

    if (nmod != trans.nrow() || nmod != trans.ncol() || nmod != lliks.nrow() || nmod != alpha.nrow()) Rcpp::stop("Unable to figure out the number of models");
    if (((double) lliks.ncol()) != totlen || ((double)alpha.ncol()) != totlen) Rcpp::stop("Seqence lengths don't match with the provided matrices");
    if (initP.ncol() != seqlens.length()) Rcpp::stop("'initP' must have as many columns as the number of sequences");

	escore_core(asVec(escore), asMat(lliks), asMat(alpha), asMat(trans), asMat(initP), enhancer_states, asVec(prior), asVec(seqlens), maxn, nthreads);
	List res = List::create(_("escore")=escore);
	return res;
}
 
static inline double fb_core(Mat<double> initPs, Mat<double> trans, Mat<double> lliks, Vec<int> seqlens, 
                             Mat<double> posteriors, Mat<double> alpha, Mat<double> new_trans, Mat<double> new_initPs, int nthreads){
    nthreads = std::max(1, nthreads);
    
    int nrow = lliks.nrow;
    int ncol = lliks.ncol;
    int nchunk = seqlens.len;
    
    //temporary objects 
    std::vector<int> chunk_startsSTD(seqlens.len, 0);
    Vec<int> chunk_starts = asVec<int>(chunk_startsSTD);
    //get the start of each chunk
    for (int i = 0, acc = 0; i < nchunk; ++i){chunk_starts[i] = acc; acc += seqlens[i];}
    
    long double tot_llik = 0;
    
    //figure out how to assign the chromosomes to the threads
    //covert seqlens to double
    std::vector<double> jobSize(nchunk);
    for (int i = 0; i < nchunk; ++i) jobSize[i] = seqlens[i];
    //get the assignments
    std::vector<int> breaks = scheduleJobs(asVec(jobSize), nthreads);
    int gretcode = FB_OK; //return code from fb_iter, aggregated across threads
    int guflowchunk = -1; //chunk that caused the underflow
    
    #pragma omp parallel num_threads(nthreads)
    {
        //each thread gets one copy of these temporaries
        std::vector<double> thread_new_transSTD(nrow*nrow, 0); 
        Mat<double> thread_new_trans = asMat(thread_new_transSTD, nrow);
        long double thread_llik = 0;
        FBtmp thread_tmp(nrow);
        int retcode = FB_OK; //return code from fb_iter
        int uflowchunk = -1; //chunk that caused the underflow

        /* Do forward and backward loop for each chunk (defined by seqlens)
         * Chunks might have very different lengths (that's why they have been scheduled). */
        #pragma omp for schedule(static) nowait
        for (int thread = 0; thread < nthreads; ++thread){
            for (int o = breaks[thread]; o < breaks[thread+1] && retcode == FB_OK; ++o){
                //o identifies the sequence/chunk
                int chunk_start = chunk_starts[o];
                int chunk_end =  chunk_start + seqlens[o];
				//std::cout << std::endl;
                retcode = fb_iter(lliks.subsetCol(chunk_start, chunk_end), 
                                posteriors.subsetCol(chunk_start, chunk_end),
								alpha.subsetCol(chunk_start, chunk_end),
								initPs.getCol(o), new_initPs.getCol(o), trans,
								thread_new_trans, thread_llik, thread_tmp);
                if (retcode >= 0) {//there was an underflow error
                    uflowchunk = o;
                }
            }
        }
        //protected access to the shared variables
        #pragma omp critical
        {
            //set gretcode to the smallest invalid retcode
            if (retcode != FB_OK && gretcode == FB_OK){
                gretcode = retcode;
                guflowchunk = uflowchunk;
            } else {
                tot_llik += thread_llik;
                for (int p = 0, q = nrow*nrow; p < q; ++p){
                    new_trans[p] += thread_new_trans[p];
                }
            }
        }
    }
	  
    //outside #pragma omp we can use exceptions :D
    if (gretcode != FB_OK) {
        if (gretcode == FB_WRONG_DIM) Rcpp::stop("Invalid array dimensions passed to 'fb_iter'");
        std::ostringstream msg_strm;
        msg_strm << "Underflow error at sequence:relpos:abspos @" <<
                    guflowchunk + 1 << ":" << gretcode + 1 << ":" << 
                    chunk_starts[guflowchunk] + gretcode + 1 << "@";
        std::string msg = msg_strm.str();
        // Rcpp::stop is broken in Rcpp 0.12.8
        // std::range_error is a workaround. I am changing it
        // only here because this is the only error message
        // that is part of an API with the epicseg package
        // (I know, bad practice to use error messages in an API).
        throw std::range_error(msg);
    }
    
    /* normalizing new_trans matrix */
    // The parallelization overhead might take longer than
    // this loop....
    for (int row = 0; row < nrow; ++row){
        double sum = 0;
        for (int col = 0; col < nrow; ++col){sum += new_trans(row, col);}
		if (sum == 0) {
		  for (int col = 0; col < nrow; ++col){
			if (col == row) {new_trans(row, col) = 1;}
			else {new_trans(row, col) = 0;}
		  }
		} else {
		  for (int col = 0; col < nrow; ++col){
			new_trans(row, col) /= sum;
		  }
		}
    }

    return (double) tot_llik;
}



using namespace Rcpp;

typedef NumericVector::iterator diter;
typedef IntegerVector::iterator iiter;

//' Forward-backward algorithm
//'
//' Forward-backward algorithm using the scaling technique.
//' That's more stable (and maybe even faster) than the method with the logarithm.
//' Warning: this function overwrites the lliks matrix. 
//' @param initP matrix of initial probabilities: each column corresponds to a sequence
//' @param trans transition matrix (rows are previous state, columns are next state)
//' @param lliks matrix with emission probabilities for each datapoint and each state.
//' Columns are datapoints and rows are states.
//' @param seqlens length of each subsequence of datapoints (set this to ncol(lliks)
//' if there is only one sequence).
//' @param posteriors the posteriors matrix where the posteriors will be written.
//' its value when the function is called does not matter, but it needs to have
//' the right dimensions (rows are states and columns are observations).
//' @param nthreads number of threads used. Sequences of observations are
//' processed independently by different threads (if \code{length(seqlens) > 1}).
//' @return a list with the following arguments:
//'    \item{posteriors}{posterior probability of being in a certain state for a certain datapoint.
//'     Same matrix used as input argument.}
//'    \item{tot_llik}{total log-likelihood of the data given the hmm model.}
//'    \item{new_trans}{update for the transition probabilities (it is already normalized).}
//' @export
// [[Rcpp::export]]
List forward_backward(NumericMatrix initP, NumericMatrix trans, NumericMatrix lliks, IntegerVector seqlens, NumericMatrix posteriors, int nthreads=1){
    int nmod = initP.nrow();
    double totlen = Rcpp::sum(seqlens);

    if (nmod != trans.nrow() || nmod != trans.ncol() || nmod != lliks.nrow() || nmod != posteriors.nrow()) Rcpp::stop("Unable to figure out the number of models");
    if (((double) lliks.ncol()) != totlen || ((double)posteriors.ncol()) != totlen) Rcpp::stop("Seqence lengths don't match with the provided matrices");
    if (initP.ncol() != seqlens.length()) Rcpp::stop("'initP' must have as many columns as the number of sequences");

    NumericMatrix newTrans(trans.nrow(), trans.ncol());
    NumericMatrix newInitP(initP.nrow(), initP.ncol());
	NumericMatrix alpha(posteriors.nrow(), posteriors.ncol());
    double tot_llik = fb_core(asMat(initP), asMat(trans), asMat(lliks), asVec(seqlens), asMat(posteriors), asMat(alpha), asMat(newTrans), asMat(newInitP), nthreads);
	List res = List::create(_("posteriors")=posteriors, _("alpha")=alpha, _("tot_llik")=tot_llik, _("new_trans")=newTrans, _("new_initP")=newInitP);
	/*
	while (alpha.length() > 0) {
	  // R has problem with garbage collection in the next line (or the destructor of alpha if this workaround is not used).
	  // alpha.~PreserveStorage() calls Rcpp_ReleaseObject, which in turn does a huge recursive search on all saved
	  // objects in RecursiveRelease. This eventually overflows the stack...
	  alpha.erase(alpha.end()-1);
	}
	*/
	return res;
}

//' Viterbi algorithm
//'
//' Standard viterbi algorithm in the log space
//' @param initP matrix of initial probabilities: each column corresponds to a sequence
//' @param trans transition matrix (rows are previous state, columns are next state)
//' @param lliks matrix with emission probabilities for each datapoint and each state.
//' Columns are datapoints and rows are states.
//' @param seqlens length of each subsequence of datapoints (set this to ncol(lliks)
//' if there is only one sequence).
//' @return a list with the following arguments:
//'    \item{vpath}{viterbi path}
//'    \item{vllik}{log-likelihood of the viterbi path}
//'    \item{vscores}{viterbi scores}
//' @export
// [[Rcpp::export]]
List viterbi(NumericMatrix vscores, NumericMatrix initP, NumericMatrix trans, NumericMatrix lliks, NumericVector seqlens, IntegerVector endstate){
  int nmod = initP.nrow();
  double totlen = Rcpp::sum(seqlens);
  if (nmod != trans.nrow() || nmod != trans.ncol() || nmod != lliks.nrow()) Rcpp::stop("Unable to figure out the number of models");
  if (((double) lliks.ncol()) != totlen) Rcpp::stop("Sequence lengths don't match with the provided matrix");
  
  int ncol = lliks.ncol();
  IntegerVector vpath(ncol);
  Mat<double> vs = asMat(vscores);
  IntegerMatrix backtrack(nmod, max(seqlens));
  std::vector<long double> scores(nmod);
  std::vector<long double> new_scores(nmod);
  
  /* avoid the Rcpp matrix object for now */
  Mat<double> myinitP = asMat(initP);
  Mat<double> mylliks = asMat(lliks);
  Mat<int> mybacktrack = asMat(backtrack);
  
  /* log-transform the transition probabilities */
  NumericMatrix ltrans(nmod,nmod);
  for (diter curr = ltrans.begin(), currt = trans.begin(); curr < ltrans.end(); ++curr, ++currt){
	*curr = log(*currt);
  }
  
  /* create transition matrix for the last step (transition to endstate) */
  NumericMatrix end_trans = trans;
  if (endstate.size() > 0) {
	// check endstate argument
	for(IntegerVector::iterator it = endstate.begin(); it != endstate.end(); ++it) {
	  if (*it < 1 || *it > trans.ncol()) {
		Rcpp::stop("invalid endstate argument");
	  }
	}
	end_trans = NumericMatrix(Dimension(trans.nrow(), trans.ncol()));
	for (int i = 0; i < trans.nrow(); i++) {
	  double sum = 0.0;
	  for(IntegerVector::iterator it = endstate.begin(); it != endstate.end(); ++it) {
		int j = *it-1;
		end_trans(i,j) = trans(i,j);
		sum += trans(i,j);
	  }
	  for(IntegerVector::iterator it = endstate.begin(); it != endstate.end(); ++it) {
		int j = *it-1;
		if (sum != 0){
		  end_trans(i,j) /= sum;
		}
	  }
	}
  }
  /* log-transform the end_transition probabilities */
  NumericMatrix end_ltrans(nmod,nmod);
  for (diter curr = end_ltrans.begin(), currt = end_trans.begin(); curr < end_ltrans.end(); ++curr, ++currt){
	*curr = log(*currt);
  }
  
  NumericMatrix ltrans_orig = ltrans; // save original log transition matrix (later it is adapted in case of endstate transitions)
  /* Viterbi independently on each chunk */
  double tot_maxscore = 0;
  for (int o = 0, chunk_start = 0; o < seqlens.length(); chunk_start += seqlens[o], ++o){
	int chunk_end = chunk_start + seqlens[o];
	/* dynamic programming */
	{
	  double* llikcol = mylliks.colptr(chunk_start);
	  double* curr_initP = myinitP.colptr(o);
	  double* vscorescol = vs.colptr(chunk_start);
	  for (int t = 0; t < nmod; ++t){
		scores[t] = llikcol[t] + log(curr_initP[t]);
		vscorescol[t] = scores[t];
	  }
	}
	
	for (int i = chunk_start + 1; i < (chunk_end); ++i){
	  /* use separate endstate transition matrix for last transition */
	  if (i == (chunk_end-1)){
		ltrans = end_ltrans;
	  } else{
		ltrans = ltrans_orig;
	  }
	  
	  double* llikcol = mylliks.colptr(i);
	  int* backtrackcol = mybacktrack.colptr(i-chunk_start);
	  double* vscorescol = vs.colptr(i);
	  for (int t = 0; t < nmod; ++t){
		int maxs = 0;
		long double maxscore = scores[0] + ltrans(0, t);
		//vscorescol[0] = maxscore;
		for (int s = 1; s < nmod; ++s){
		  long double currscore = scores[s] + ltrans(s,t);
		  //vscorescol[s] = currscore;
		  
		  if (currscore > maxscore){
			maxscore = currscore;
			maxs = s;
		  }
		}
		backtrackcol[t] = maxs;
		new_scores[t] = llikcol[t] + maxscore;
		vscorescol[t] = new_scores[t];
	  }
      
	  memcpy(scores.data(), new_scores.data(), sizeof(long double)*nmod);
	}
	
	/* backtracking */
	int maxp = 0;
	double maxscore = scores[0];
	for (int p = 1; p < nmod; ++p){
	  if (scores[p] > maxscore){
		maxscore = scores[p];
		maxp = p;
	  }
	}
	tot_maxscore += maxscore;
	vpath[chunk_end - 1] = maxp + 1;
	for (int i = chunk_end - 2; i >= chunk_start; --i){
	  maxp = backtrack(maxp, i - chunk_start + 1);
	  vpath[i] = maxp + 1; //in R indices are 1-based
	}
  }
  return List::create(_("vpath")=vpath, _("vllik")=tot_maxscore, _("vscores")=vscores);
}

/*
// [[Rcpp::export]]
Rcpp::IntegerVector orderColumns(Rcpp::IntegerMatrix mat){
    Rcpp::IntegerVector order(mat.ncol());
    orderColumns_core(asMat(mat), asVec(order));
    return order;
}
*/


// [[Rcpp::export]]
Rcpp::List testSchedule(Rcpp::NumericVector jobs, int nthreads, int type){
    std::vector<int> breaks(nthreads+1);
    Vec<double> jobSize = asVec(jobs);
    if (type == 0){
        scheduleNaive(jobSize, breaks);
    } else if (type == 1){
        scheduleGreedy(jobSize, breaks);
    } else if (type == 2){
        scheduleOptimal(jobSize, breaks);
    } else Rcpp::stop("invalid type");
    double makespan = getMakespan(jobSize, breaks);
    return Rcpp::List::create(  
        Rcpp::Named("makespan")=makespan, 
        Rcpp::Named("breaks")=Rcpp::wrap(breaks));
}

//static inline void collapsePosteriors_core(Mat<double> cpost, Mat<double> post, NMPreproc& preproc, int nthreads=1)

// [[Rcpp::export]]
Rcpp::NumericMatrix testColPost(Rcpp::NumericMatrix post, Rcpp::List m2u, int nthreads){
    Rcpp::NumericVector values = Rcpp::as<Rcpp::NumericVector>(m2u["values"]);
    Rcpp::IntegerVector map = Rcpp::as<Rcpp::IntegerVector>(m2u["map"]);
    if (post.ncol() != map.length()) Rcpp::stop("posteriors doesn't match with m2u");
    
    Rcpp::NumericMatrix smallerPost(post.nrow(), values.length());
    Vec<double> foo; NMPreproc preproc(asVec(values), asVec(map), foo);
    collapsePosteriors_core(asMat(smallerPost), asMat(post), preproc);
    return smallerPost;
}
