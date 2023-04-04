#include <Rcpp.h>
#include "array.cpp"
#include "tabbing.h"
#include <algorithm> 
#include <unordered_map>
#include <R_ext/BLAS.h>

//parallel sorting. SUPPORT_OPENMP is not a guarantee that
//the header parallel/algorithm is present, but I don't have better ideas
//other than using autoconf...
#ifdef SUPPORT_OPENMP
    #include <parallel/algorithm>
    template <class RandomAccessIterator>
    inline void parallelSort(RandomAccessIterator first, RandomAccessIterator last, int nthreads){
        __gnu_parallel::sort(first, last, __gnu_parallel::parallel_tag(nthreads));
    }
#else
    template <class RandomAccessIterator>
    inline void parallelSort(RandomAccessIterator first, RandomAccessIterator last, int /*nthreads*/){
        std::sort(first, last);
    }
#endif



//splits the range [0, len) into nthreads ranges of approximately same size
//the ith range is [breaks[i], breaks[i+1])
std::vector<int> makeBreaks(int len, int nthreads){
	std::vector<int> breaks(nthreads+1);
	double denom = (double)nthreads;
	for (int i = 0; i <= nthreads; ++i){ 
		breaks[i] = round((i/denom)*len); 
	}
	return breaks;
}

//computes the kth pair (i,j) such that j<=i (or j < i if noDiag is true)
//the sequence of pairs looks like this: (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), ...
inline void kthPair(int k, int& i, int& j, bool noDiag=false){
	i = (int)((-1 + sqrt(1 + 8*k))*0.5);
	j = k - (i*(i+1))/2;
	if (noDiag) ++i;
}

//computes the total number of unique unordered pairs (i,j)
//if noDiag==true, then it excludes the pairs s.t. i==j
inline int nPairs(int n, bool noDiag=false){
	if (noDiag) --n;
	return (n*(n+1))/2;
}

// [[Rcpp::export]]
Rcpp::IntegerVector tabFast(Rcpp::IntegerVector counts){
	std::vector<int> tab(100);
	tabFast_impl(counts.begin(), counts.end(), tab);
	return Rcpp::wrap(tab);
}

//equivalent to apply(counts, 1, tabFast)
std::vector<std::vector<int> > tabRows(Rcpp::IntegerMatrix counts, int nthreads=1){
	int nrow = counts.nrow();
	int ncol = counts.ncol();
	nthreads = std::min(nrow, nthreads);
	std::vector<std::vector<int> > result(nrow);
	for (int i = 0; i < nrow; ++i) result[i].resize(100);
	int chunkSize = 1e5/nrow;
	int chunkStart = 0;
	while (chunkStart < ncol){
		int chunkEnd = std::min(chunkStart + chunkSize, ncol);
		#pragma omp parallel for num_threads(nthreads)
		for (int i = 0; i < nrow; ++i){
			Rcpp::MatrixRow<INTSXP> row = counts.row(i);
			tabFast_impl(row.begin() + chunkStart, row.begin() + chunkEnd, result[i], false);
		}
		chunkStart = chunkEnd;
	}
	
	for (int i = 0; i < nrow; ++i) shrink(result[i]);
	
	return result;
}

//method for doing quantile normalization between an empirical and a theoretical distribution
//cfreq is the cumulative frequency of some empirical counts (position zero corresponds to count zero and so on)
//cdf is a callable function or class such that cdf(i) gives the cumulative prob
//ass is an assignment of the same length of cfreq (specified by flen) that assigns an observed count
//to a theoretical count
//there are many ways of dealing with discretization, here 
//we compute the assignment such that the norm ||CDF_empirical_after_assignment - CDF_theoretical||_1 is minimized.
template<typename Tfun>
void quantNormDiscrete_inner(double* cfreq, int* ass, int flen, Tfun cdf){
	int j = 0;
	bool lastSep = false;
	if (flen <= 0) return;
	if (abs(cfreq[flen-1]-1) >= 1e-10) Rcpp::stop("cumulative frequencies not tending to one");
	for (int i = 0; !lastSep; ++i){
		double P = cdf(i);
		if (P < 0 && P > 1) Rcpp::stop("Invalid cumulative distribution function");
		while (cfreq[j] < P && j < flen) ++j;
		//std::cout << "j: " << j << " i: " << i << std::endl;
		double diffMinus = j==0?P:(P-cfreq[j-1]);
		double diffPlus = cfreq[j] - P;
		//std::cout << "P: " << P << " cfreq[j]: " << cfreq[j] << std::endl;
		//std::cout << "diffMinus: " << diffMinus << " diffPlus: " << diffPlus << std::endl;
		if (diffMinus < diffPlus){
			++ass[j];
			//std::cout << "diffMinus" << std::endl;
		} else {
			//std::cout << "diffPlus" << std::endl;
			if (j >= flen-1) lastSep = true;
			else ++ass[j+1];
		}
	}
	
	int currPoint = 0;
	for (int i = 0; i < flen; ++i){
		currPoint += ass[i];
		ass[i] = currPoint;
	}
}



template<typename Titer>
std::vector<double> getCFreq(Titer begin, Titer end){
	int flen = end-begin;
	std::vector<double> cfreq(flen);
	long double acc = 0;
	for (int i = 0; begin < end; ++i, ++begin){
		acc += *begin;
		cfreq[i] = acc;
	}
	if (acc != 1){//make sure that frequencies sum up to 1
		for (int i = 0; i < flen; ++i) cfreq[i] /= acc;
	}
	
	return cfreq;
}


struct theoCDF {
	std::vector<double>& cdf;
	double operator() (int i){
		if (i < 0) return 0;
		if (i >= ((int)cdf.size())) return 1;
		else return cdf[i];
	}
	theoCDF(std::vector<double>& acdf) : cdf(acdf){}
};

struct poisCDF {
	double lambda;
	double operator() (int i){
		return Rf_ppois(i, lambda, true, false);
	}
};

//quantile normalization of the empirical distribution (just a vector of count frequencies)
//according to a theoretical distribution (also a vector of count frequenceis)

// [[Rcpp::export]]
Rcpp::IntegerVector labelCounts(Rcpp::NumericVector empirical, Rcpp::NumericVector theoretical){
	std::vector<double> cfreqE = getCFreq(empirical.begin(), empirical.end());
	std::vector<double> cfreqT = getCFreq(theoretical.begin(), theoretical.end());
	
	theoCDF cdf(cfreqT);
	int flen = empirical.length();
	Rcpp::IntegerVector assignment(flen);
	quantNormDiscrete_inner(cfreqE.data(), assignment.begin(), flen, cdf);
	return assignment;
}


//coords and counts have the same format
//clusters is an object to form clusters (seeds) based on the counts of a single histone mark
//typically clusters is produced by doing something like: 
//clusters <- apply(counts, 1, function(row) labelCounts(tabFast(row), theoretical=rep(0.1, 10)))
//each column of the counts matrix will be assigned to nrow(counts) different seeds,
//there are nrow(counts) independent sets of seeds
//column i, when looking at row j, will be assigned to the seed in the set j specified by clusters[[j]][coords[i,j] + 1]
// [[Rcpp::export]]
Rcpp::List clusterAverages2(Rcpp::NumericMatrix counts, Rcpp::NumericMatrix coords, Rcpp::List clusters, int nthreads=1){
  if (coords.ncol()!=counts.ncol()||coords.nrow()!=counts.nrow()) Rcpp::stop("counts and coords must have the same format");
  int ncomp = counts.nrow();
  int ncol = counts.ncol();
  if (ncomp != clusters.length()) Rcpp::stop("one set of clusters for each row of the count matrix is required");
  //switch to the array of arrays format
  std::vector<Rcpp::NumericVector> Cs;
  for (int i = 0; i < ncomp; ++i) {
    Rcpp::NumericVector v(clusters[i]);
    Cs.push_back(v);
  }
  //determine the maximum number of different clusters per row of count matrix
  int maxclust = 0;
  for (int i = 0; i < ncomp; ++i){
    int m = *(std::max_element(Cs[i].begin(), Cs[i].end()));
    maxclust = std::max(maxclust, m);
  }
  maxclust += 1;
  
  //allocate memory
  Rcpp::NumericMatrix mus(ncomp, maxclust*ncomp);
  std::vector<long> musvec(ncomp*maxclust*ncomp);
  Rcpp::NumericMatrix sigmasqs(ncomp, maxclust*ncomp);
  Rcpp::IntegerVector sizes(maxclust*ncomp);
  
  //avoid Rcpp::Matrix for now (buggy with long vectors)
  Mat<double> mycounts = asMat(counts);
  Mat<double> mycoords = asMat(coords);
  
  //parallelize on the columns of the matrix
  #pragma omp parallel num_threads(nthreads)
  {
    //accumulators local to this thread
    std::vector<double> vsum(ncomp*maxclust*ncomp);
    std::vector<int> vnum(ncomp*maxclust);
    #pragma omp for nowait
    for (int i = 0; i < ncol; ++i){
      double* tcol = mycounts.colptr(i);
      double* dcol = mycoords.colptr(i);
      for (int j = 0; j < ncomp; ++j){
	double c = dcol[j];
	if (c < 0 || c >= Cs[j].length()) Rcpp::stop("invalid clustering or invalid counts");
	int clust = Cs[j][c];
	if (clust < 0 || clust >= maxclust) Rcpp::stop("something wrong in detecting maxclust");
	clust = (j*maxclust + clust);
	int offset = clust*ncomp;
	for (int k = 0; k < ncomp; ++k){
	  vsum[offset + k] += tcol[k];
	}
	++vnum[clust];
      }
    }
    //threads accumulator update the global accumulators in a protected section
    #pragma omp critical
    {
      for (int i = 0, e = vsum.size(); i < e; ++i){
	mus[i] += vsum[i];
      }
      for (int i = 0, e = vnum.size(); i < e; ++i){
	sizes[i] += vnum[i];
      }
    }
  }
  for (int i = 0, e = mus.ncol(); i < e; ++i){
    int n = sizes[i];
    if (n > 0){
      for (int j = 0; j < ncomp; ++j){
	mus(j, i) /= n;
	musvec[j + ncomp*i] = mus(j, i); // for easy usage during calculation of sigmasq
      }
    }
  }

  // calculate sigmasqs
  #pragma omp parallel num_threads(nthreads)
  {
    //accumulators local to this thread
    std::vector<double> sig(ncomp*maxclust*ncomp);
    std::vector<int> vnum(ncomp*maxclust);
    Mat<double> musMat = asMat(mus);

    #pragma omp for nowait
    for (int i = 0; i < ncol; ++i){
      double* tcol = mycounts.colptr(i);
      double* dcol = mycoords.colptr(i);
      for (int j = 0; j < ncomp; ++j){
	double c = dcol[j];
	if (c < 0 || c >= Cs[j].length()) Rcpp::stop("invalid clustering or invalid counts");
	int clust = Cs[j][c];
	if (clust < 0 || clust >= maxclust) Rcpp::stop("something wrong in detecting maxclust");
	clust = (j*maxclust + clust);
	int offset = clust*ncomp;
	for (int k = 0; k < ncomp; ++k){
	  double mu = musvec[offset + k];
	  sig[offset + k] += std::pow((tcol[k] - mu), 2);
	}
	++vnum[clust];
      }
    }
    //threads accumulator update the global accumulators in a protected section
    #pragma omp critical
    {
      for (int i = 0, e = sig.size(); i < e; ++i){
	sigmasqs[i] += sig[i];
      }
      for (int i = 0, e = vnum.size(); i < e; ++i){
	sizes[i] += vnum[i];
      }
    }
  }
  for (int i = 0, e = sigmasqs.ncol(); i < e; ++i){
    int n = sizes[i];
    if (n > 0){
      for (int j = 0; j < ncomp; ++j){
	sigmasqs(j, i) /= n;
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("mus")=mus, Rcpp::Named("sigmasqs")=sigmasqs, Rcpp::Named("sizes")=sizes);
}


// [[Rcpp::export]]
Rcpp::List clusterAverages(Rcpp::NumericMatrix counts, Rcpp::List clusters, int nthreads=1){
	return clusterAverages2(counts, counts, clusters, nthreads);
}

//fill it in, written in c:
//for (i in 1:ncol(counts)){
//	for (j in 1:nrow(counts)){
//		clust <- multidclust[[j]][counts[j,i]+1] + 1
//		posteriors[clust,i] <- posteriors[clust,i] + 1
//	}
//}

// [[Rcpp::export]]
Rcpp::NumericMatrix fillPosteriors(Rcpp::IntegerMatrix coords, Rcpp::List clusters, int nclust, int nthreads=1){
	int ncomp = coords.nrow();
	int ncol = coords.ncol();
	if (ncomp != clusters.length()) Rcpp::stop("one set of clusters for each row of the count matrix is required");
	Rcpp::NumericMatrix posteriors(nclust, ncol);
	
	//switch to the std::vector of Rcpp::vector format
	std::vector<Rcpp::IntegerVector> Cs;
	for (int i = 0; i < ncomp; ++i) {
		Rcpp::IntegerVector v(clusters[i]);
		Cs.push_back(v);
	}
	
	//avoid Rcpp matrix for now
	Mat<double> myposteriors = asMat(posteriors);
	Mat<int> mycoords = asMat(coords);
	
	//parallelize on the columns of the matrix
	#pragma omp parallel for num_threads(nthreads)
	for (int i = 0; i < ncol; ++i){
		int* ccol = mycoords.colptr(i);
		double* pcol = myposteriors.colptr(i);
		for (int j = 0; j < ncomp; ++j){
			int c = ccol[j];
			if (c < 0 || c >= Cs[j].length()) Rcpp::stop("invalid clustering or invalid coords");
			int clust = Cs[j][c];
			if (clust < 0 || clust >= nclust) Rcpp::stop("count mapped to invalid cluster");
			++pcol[clust];
		}
	}
	
	return posteriors;
}

template<typename TVec1, typename TVec2, typename TCoeff>
inline void addCol(int n, TVec1 v1, TVec2 v2, TCoeff c){
	for (int i = 0; i < n; ++i) v1[i] += c*v2[i];
}

template<typename TVec1, typename TVec2>
inline void addAllCols(int n, TVec1 v1, TVec2 v2){
	for (int i = n; i > 0; --i) addCol(i, v1 + (i-1)*n, v2, v2[i-1]);
}


//computes the pairwise dot-product between any pair of rows and 
//divides it by the total number of rows (-1 if besselCorr==true)
//in R: 
//(counts %*% t(counts)) / (nrow(counts)-1)

// [[Rcpp::export]]
Rcpp::NumericMatrix rowdotprod(Rcpp::NumericMatrix counts, bool besselCorr=true, int nthreads=1){
	int ncol = counts.ncol();
	int nrow = counts.nrow();
	std::vector<long double> std_cov(nrow*nrow);
	Mat<long double> cov = asMat(std_cov, nrow);
	Mat<double> CNTS = asMat(counts);
	
	#pragma omp parallel num_threads(nthreads)
	{
		std::vector<long double> t_std_cov(nrow*nrow);
		long double* t_ptr_cov = t_std_cov.data();
		Mat<long double> t_cov = asMat(t_std_cov, nrow);
		#pragma omp for nowait
		for (int col = 0; col < ncol; ++col){
			addAllCols(nrow, t_ptr_cov, CNTS.colptr(col));
		}
		#pragma omp critical
		{
			for (int j = 0; j < nrow; ++j){
				for (int i = 0; i <= j; ++i){
					cov(i,j) += t_cov(i,j);
				}
			}
		}
	}
	Rcpp::NumericMatrix result(nrow, nrow);
	long double denom = besselCorr?(ncol-1):ncol;
	for (int j = 0; j < nrow; ++j){
		for (int i = 0; i <= j; ++i){
			result(i,j) = result(j,i) = cov(i,j)/denom;
		}
	}
	return result;
}

/*
//matrix multiplication. Parallelization is done on the columns of mat2 and mat3
//that works only if you don't have to allocate mat3, otherwise R is faster
// [[Rcpp::export]]
void matprod(Rcpp::NumericMatrix mat1, Rcpp::NumericMatrix mat2, Rcpp::NumericMatrix mat3, int nthreads=1){
	int ncol1 = mat1.ncol();
	int nrow1 = mat1.nrow();
	int ncol2 = mat2.ncol();
	int nrow2 = mat2.nrow();
	int ncol3 = mat3.ncol();
	int nrow3 = mat3.nrow();
	if (ncol1 != nrow2 || nrow1 != nrow3 || ncol2 != ncol3) Rcpp::stop("non conformable matrices");
	
	double* M1 = mat1.begin();
	double* M2 = mat2.begin();
	double* M3 = mat3.begin();
	const char *transa = "N", *transb = "N";
	double one = 1.0, zero = 0.0;
	std::vector<int> breaks(nthreads+1);
	for (int i = 0; i <= nthreads; ++i){
		breaks[i] = round(i*(((double)ncol3)/nthreads));
	}
	#pragma omp parallel for num_threads(nthreads)
	for (int i = 0; i < nthreads; ++i){
		double* TM2 = M2 + breaks[i]*nrow2;
		double* TM3 = M3 + breaks[i]*nrow3;
		int ncol23 = breaks[i+1]-breaks[i];
		F77_CALL(dgemm)(transa, transb, &nrow1, &ncol23, &ncol1, &one,
				M1, &nrow1, TM2, &nrow2, &zero, TM3, &nrow1);
	}
}
*/

template<typename TIter>
inline void rangeImpl(double* min, double* max, TIter C, TIter E){
	for (; C < E; ++C){
		double c = *C;
		if (c < *min) *min = c;
		if (c > *max) *max = c;
	}
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix discretizeRows(Rcpp::NumericMatrix scores, int nlevels, int nthreads=1){
	int nrow = scores.nrow();
	int ncol = scores.ncol();
	
	std::vector<double> mins(nrow, std::numeric_limits<double>::infinity());
	std::vector<double> maxs(ncol, -std::numeric_limits<double>::infinity());
	int chunkSize = 1e5/nrow; if (chunkSize < 1) chunkSize = 1;
	int nchunks = ceil(scores.ncol()/((double)chunkSize));
	#pragma omp parallel num_threads(nthreads)
	{	
		std::vector<double> tmins(mins);
		std::vector<double> tmaxs(maxs);
		#pragma omp for nowait
		for (int i = 0; i < nchunks; ++i){
			int chunkStart = chunkSize*i;
			int chunkEnd = std::min(chunkStart + chunkSize, ncol);
			for (int j = 0; j < nrow; ++j){
				Rcpp::MatrixRow<REALSXP>::iterator rowIter = scores.row(j).begin();
				rangeImpl(&mins[j], &maxs[j], rowIter + chunkStart, rowIter + chunkEnd);
			}
		}
		#pragma omp critical
		{
			for (int j = 0; j < nrow; ++j){
				mins[j] = std::min(mins[j], tmins[j]);
				maxs[j] = std::max(maxs[j], tmaxs[j]);
			}
		}
	}
	
	//overwrite the maxs with the deltas
	std::vector<double>& deltas = maxs;
	for (int j = 0; j < nrow; ++j){
		deltas[j] = (maxs[j]-mins[j]) / nlevels;
	}
	
	
	Rcpp::IntegerMatrix mat(nrow, ncol);
	Mat<double> orig = asMat(scores);
	Mat<int> neww = asMat(mat);
	#pragma omp parallel for num_threads(nthreads)
	for (int i = 0; i < ncol; ++i){
		double* O = orig.colptr(i);
		int* N = neww.colptr(i);
		for (int j = 0; j < nrow; ++j){
			N[j] = (O[j] - mins[j])/deltas[j];
		}
	}
	
	mat.attr("dimnames") = scores.attr("dimnames");
	return mat;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix splitAxes(Rcpp::NumericMatrix scores, int nsplit, int nthreads=1){
	int nrow = scores.nrow();
	int ncol = scores.ncol();
	
	//allocate the matrix to be returned
	Rcpp::IntegerMatrix mat(nrow, ncol);
	
	//allocate the necessary memory
	std::vector<std::pair<double,int> > avatars(ncol);
	for (int j = 0; j < nrow; ++j){
		Rcpp::MatrixRow<REALSXP> scrow = scores.row(j);
		Rcpp::MatrixRow<INTSXP> matrow = mat.row(j);
		//fill in the avatars vector
		#pragma omp parallel for num_threads(nthreads)
		for (int i = 0; i < ncol; ++i){
			avatars[i].first = scrow[i];
			avatars[i].second = i;
		}
		//sort it
		parallelSort(avatars.begin(), avatars.end(), nthreads);
		//fill in the matrix
		double mult = nsplit/((double)ncol);
		#pragma omp parallel for num_threads(nthreads)
		for (int i = 0; i < ncol; ++i){
			int a = avatars[i].second;
			matrow[a] = floor(i*mult);
		}
	}
	
	mat.attr("dimnames") = scores.attr("dimnames");
	return mat;
}

//same as before, but with count sort
// [[Rcpp::export]]
Rcpp::IntegerMatrix splitAxesInt(Rcpp::IntegerMatrix scores, int nsplit, int nthreads=1){
	int nrow = scores.nrow();
	int ncol = scores.ncol();
	
	//allocate the matrix to be returned
	Rcpp::IntegerMatrix mat(nrow, ncol);
	
	//tabulate the values for each row
	std::vector<std::vector<int> > tabs = tabRows(scores, nthreads=nthreads);
	
	#pragma omp parallel for num_threads(nthreads)
	for (int j = 0; j < nrow; ++j){
		Rcpp::MatrixRow<INTSXP> scrow = scores.row(j);
		Rcpp::MatrixRow<INTSXP> matrow = mat.row(j);
		//get abundance of each symbol
		std::vector<int>& tab = tabs[j];
		//we need the cumulative sums of it
		
		
		int acc = 0;
		for (int i = 0, e = tab.size(); i < e; ++i){
			int tmp = tab[i];
			tab[i] = acc;
			acc += tmp;
		}
		//fill in mat
		double mult = nsplit/((double)ncol);
		for (int i = 0; i < ncol; ++i){
			int pos = tab[scrow[i]]++;
			matrow[i] = floor(mult*pos);
		}
	}
	
	mat.attr("dimnames") = scores.attr("dimnames");
	return mat;
}

//symmetrized kullback-leibler divergence between independent multivariate poissons
//lmus1 and lmus2 are the logs of mus1 and mus2. In case mus1[i] is null, lmus1[i] must be a finite value (like -DBL_MAX).
inline double Symm_KLdiv_pois(double* mus1, double* mus2, double* lmus1, double* lmus2, int len){
	double logRatioSum = 0; //sum of (mus1[i]-mus2[i])*log(mus1[i]/mus2[i])
	bool retInf = false;
	for (int i = 0; i < len; ++i){
		retInf = retInf || ((mus1[i]>0)^(mus2[i]>0)); // if one of mus1 and mus2 is positive and the other negative, retInf becomes TRUE and 'inf' is returned
		logRatioSum += (mus1[i]-mus2[i])*(lmus1[i]-lmus2[i]);
	}
	if (retInf) return std::numeric_limits<double>::infinity();
	return logRatioSum*.5;
}

//symmetrized kullback-leibler divergence between two lognormal distributions
inline double Symm_KLdiv_LN(double* mus1, double* mus2, double* sigmasqs1, double* sigmasqs2, int len){
  double logRatioSum = 0; //sum of (mus1[i]-mus2[i])*log(mus1[i]/mus2[i])
  //bool retInf = false;
  for (int i = 0; i < len; ++i){
	if (mus1[i] == mus2[i] & sigmasqs1[i] == sigmasqs2[i]) logRatioSum += 0 ; //this deals with the case that sigmasq is equal to zero (but so far only if they both are)
    else logRatioSum += .5 * (sigmasqs1[i]/sigmasqs2[i] + sigmasqs2[i]/sigmasqs1[i] + (1/sigmasqs1[i] + 1/sigmasqs2[i]) * std::pow((mus1[i] - mus2[i]), 2)) - 1;
  }
  return logRatioSum;
}

bool allPos(Rcpp::NumericVector v){
  for (int i = 0, e = v.length(); i < e; ++i){
    if (v[i] < 0) return false;
  }
  return true;
}

//compute symmetric kullback leibler divergence between all pairs of lognormals
// [[Rcpp::export]]
Rcpp::NumericMatrix KL_dist_mat_LN(Rcpp::NumericMatrix mus, Rcpp::NumericMatrix sigmasqs, int nthreads=1){
  //check that all elements of sigmasqs are positive
  if (!allPos(sigmasqs)) Rcpp::stop("the sigmasqs matrix cannot contain negative numbers");
  
  int nrow = mus.nrow();
  int ncol = mus.ncol();

  //along the diagonal the distance is zero, because Symm_KL_div(musi, musi, sigmasqi, sigmasqi, nrow)==0
  Rcpp::NumericMatrix dmat(ncol, ncol);
  int npairs = nPairs(ncol, true);
  std::vector<int> breaks = makeBreaks(npairs, nthreads);
  
  #pragma omp parallel for num_threads(nthreads)
  for (int t = 0; t < nthreads; ++t){
    int k = breaks[t], lastK = breaks[t+1]; 
    int i, j; kthPair(k, i, j, true);
    while (k < lastK){
      
      //double* mus1 = mus.colptr(i); double* sigmasqs1 = sigmasqs.colptr(i);
      //double* mus2 = mus.colptr(j); double* sigmasqs2 = sigmasqs.colptr(j);
      double* mus1 = mus.column(i).begin(); double* sigmasqs1 = sigmasqs.column(i).begin();
      double* mus2 = mus.column(j).begin(); double* sigmasqs2 = sigmasqs.column(j).begin();

      dmat(i,j) = Symm_KLdiv_LN(mus1, mus2, sigmasqs1, sigmasqs2, nrow);
      dmat(j,i) = dmat(i,j);
      
      ++k; ++j; if (j == i) { j = 0; ++i; }
    }
  }
  return dmat;
}

//compute symmetric kullback leibler divergence between all pairs of negative multinomials
// [[Rcpp::export]]
Rcpp::NumericMatrix KL_dist_mat(Rcpp::NumericMatrix nbs, double r, int nthreads=1){
	//check that all elements are positive
	if (!allPos(nbs)) Rcpp::stop("the matrix cannot contain negative numbers");
	//check that r is positive
	if (r <= 0) Rcpp::stop("r must be strictly positive");
	
	int nrow = nbs.nrow();
	int ncol = nbs.ncol();

	//precompute all the logarithms
	int len = nrow*ncol;
	std::vector<double> lognbs(len); 
	#pragma omp parallel for num_threads(nthreads)
	for (int i = 0; i < len; ++i){
		lognbs[i] = nbs[i]>0?log(nbs[i]):(-DBL_MAX);
	}
	Mat<double> lnbs = asMat(lognbs, ncol);
	//precompute all column sums and their logarithms + r(not needed for the poisson case)
	std::vector<double> csums(ncol);
	std::vector<double> lcsumsr(ncol);
	if (R_FINITE(r)) {
		//column sums
		colSums(asMat(nbs), asVec(csums), nthreads);
		//logarithms of the column sums + r
		#pragma omp parallel for num_threads(nthreads)
		for (int i = 0; i < ncol; ++i){
			lcsumsr[i] = log(csums[i] + r);
		}
	}
	
	//along the diagonal the distance is zero, because Symm_KL_div(musi, musi, r, nrow)==0
	Rcpp::NumericMatrix dmat(ncol, ncol);
	int npairs = nPairs(ncol, true);
	std::vector<int> breaks = makeBreaks(npairs, nthreads);
	
	#pragma omp parallel for num_threads(nthreads)
	for (int t = 0; t < nthreads; ++t){
		int k = breaks[t], lastK = breaks[t+1]; 
		int i, j; kthPair(k, i, j, true);
		while (k < lastK){
			
			double* mus1 = nbs.column(i).begin(); double* lmus1 = lnbs.colptr(i);
			double* mus2 = nbs.column(j).begin(); double* lmus2 = lnbs.colptr(j);
			double mu1 = csums[i]; double lmur1 = lcsumsr[i];
			double mu2 = csums[j]; double lmur2 = lcsumsr[j];
			dmat(i,j) = Symm_KLdiv_pois(mus1, mus2, lmus1, lmus2, nrow) + //poisson contribution
									(mu1 - mu2)*(lmur2 - lmur1)*0.5; //neg binom contribution
			dmat(j,i) = dmat(i,j);
			
			 ++k; ++j; if (j == i) { j = 0; ++i; }
		}
	}
	
	return dmat;
}

unsigned hashVec(Vec<int> v){
	std::hash<int> hasher;
	unsigned hashSum = 0;
	int len = v.len;
	for (int i = 0; i < len; ++i){
		hashSum += hasher(v[i]);
	}
	return hashSum;
}

bool sameVec(Vec<int> v1, Vec<int> v2){
	if (v1.len != v2.len) return false;
	for (int i = 0, e = v1.len; i < e; ++i){
		if (v1[i] != v2[i]) return false;
	}
	return true;
}

//permutation should be a permutation of the numbers from 1 to ncol(counts), 
//it can be created in R by doing: sample(ncol(counts), ncol(counts))
// [[Rcpp::export]]
Rcpp::IntegerVector findUniqueSeeds(Rcpp::IntegerMatrix counts, Rcpp::IntegerVector permutation, int k){
	if (counts.ncol() != permutation.length()) Rcpp::stop("matrix and permutation don't match");
	int ncol = counts.ncol();
	Mat<int> mat = asMat(counts);
	//finding the seeds
	typedef std::unordered_map<Vec<int>, int, std::function<unsigned(Vec<int>)>,std::function<bool(Vec<int>,Vec<int>)> > Maptype; 
	Maptype map(2*k, hashVec, sameVec);
	for (int i = 0; i < ncol && map.size() < (unsigned)k; ++i){
		map[mat.getCol(permutation[i]-1)] = permutation[i];
	}
	if (map.size() < (unsigned)k) Rcpp::stop("unable to find enough distinct columns");
	//writing them in a new vector
	Rcpp::IntegerVector seeds(k);
	typedef Maptype::iterator Mapiter;
	int pos = 0;
	for (Mapiter s = map.begin(), e = map.end(); s != e; ++s, ++pos){
		seeds[pos] = s->second;
	}
	return seeds;
}
