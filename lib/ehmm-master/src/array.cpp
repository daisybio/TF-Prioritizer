#pragma once

#include <Rcpp.h>
#include <stdexcept> 
#include <iostream>
#ifdef SUPPORT_OPENMP
    #include <omp.h>
#endif
#include <csignal>

/* Matrix and vector classes.
 * All these classes do not own the memory (no alloc/dealloc), they
 * just wrap pointers to memory allocated/deallocated somewhere else */

template<typename T>
struct Vec {
	T* ptr;
	int len;
	
	Vec(){}
	
	Vec(T* _ptr, int _len):
		ptr(_ptr), len(_len){}
	
	inline T operator[] (int i) const {return ptr[i];}
	inline T& operator[] (int i){return ptr[i];}
	
	Vec subset(int start, int end){
	  if (start < 0 || start >= len) {
	    raise(SIGSEGV);
	  }
	  if (end < start || end > len) {
	    raise(SIGSEGV);
	  }
		return Vec(ptr + start, end - start);
	}
	
	Vec(SEXP x){
		const int RType = Rcpp::traits::r_sexptype_traits<T>::rtype;
		//no memory allocation
		if (TYPEOF(x) != RType) Rcpp::stop("incompatible types");

		Rcpp::Vector< RType > rcpp_vec(x);
		ptr = rcpp_vec.begin();
		len = rcpp_vec.length();
	}
};

//vector where elements correspond to non necessarily adjacent areas of memory 
template<typename T>
struct GapVec {
	T* ptr;
	int* set;
	int len;
	
	GapVec(){}
	
	GapVec(T* _ptr, int* _set, int _len):
		ptr(_ptr), set(_set), len(_len){}
	
	inline T operator[] (int i) const {return ptr[set[i]];}
	inline T& operator[] (int i){return ptr[set[i]];}
	
	GapVec subset(int start, int end){
		return GapVec(ptr, set + start, end - start);
	}
};

//normal matrix
template<typename T>
struct Mat {
	T* ptr;
	int nrow;
	int ncol;
	
	Mat(){}
	
	Mat(T* _ptr, int _nrow, int _ncol):
		ptr(_ptr), nrow(_nrow), ncol(_ncol){}
	
	inline T operator[] (int i) const {
	  if (i < 0 || i >= nrow*ncol) {
	    raise(SIGSEGV);
	  }
    return ptr[i];
	}
	inline T& operator[] (int i){
	  if (i < 0 || i >= nrow*ncol) {
	    raise(SIGSEGV);
	  }
	  return ptr[i];
	}
	
	inline T* colptr(int col){
	  if (col < 0 || col >= ncol) {
      raise(SIGSEGV);
	  }
		return ptr + ((long)col)*nrow;
	}

	inline T operator() (int row, int col) const {
	  if (row < 0 || row >= nrow) {
	    raise(SIGSEGV);
	  }
	  if (col < 0 || col >= ncol) {
	    raise(SIGSEGV);
	  }
	  return *(colptr(col) + row);
	}
	inline T& operator() (int row, int col) {
	  if (row < 0 || row >= nrow) {
	    raise(SIGSEGV);
	  }
	  if (col < 0 || col >= ncol) {
	    raise(SIGSEGV);
	  }
	  return *(colptr(col) + row);
	}
	
	Mat subsetCol(int colStart, int colEnd){
		return Mat(colptr(colStart), nrow, colEnd-colStart);
	}
	
	
	inline T get(int row, int col){
	  if (row < 0 || row >= nrow) {
	    raise(SIGSEGV);
	  }
	  if (col < 0 || col >= ncol) {
	    raise(SIGSEGV);
	  }
		return *(colptr(col) + row);
	}
	
	inline Vec<T> getCol(int col){
	  if (col >= ncol) {
	    raise(SIGSEGV);
	  }
		return Vec<T>(colptr(col), nrow);
	}
	
	Mat(SEXP x){
		const int RType = Rcpp::traits::r_sexptype_traits<T>::rtype;
		
		//no memory allocation
		if (TYPEOF(x) != RType) Rcpp::stop("incompatible types");
		Rcpp::Matrix< RType > rcpp_mat(x);
		ptr = rcpp_mat.begin();
		nrow = rcpp_mat.nrow();
		ncol = rcpp_mat.ncol();
	}
};

//matrix corresponding to a sliding window of length nrow on a sequence *ptr, the window slides step positions every time
template<typename T>
struct SWMat {
	T* ptr;
	int nrow;
	int ncol;
	int step;
	
	inline T* colptr(int col){
		return ptr + col*step;
	}
	
	//I never used this function so far, it's just a reminder 
	//about how long the pointer should be at least
	inline int ptrlen(){
		return nrow + (ncol-1)*step;
	}
	
	SWMat(){}
	
	SWMat(T* _ptr, int _nrow, int _ncol, int _step):
		ptr(_ptr), nrow(_nrow), ncol(_ncol), step(_step){}
	
	SWMat(Vec<T> vec, int _nrow, int _step){
		ptr = vec.ptr;
		nrow = _nrow;
		step = _step;
		
		if ((vec.len - nrow) % step != 0) throw std::invalid_argument("the window can be slid a fractional number of times...");
		
		ncol = (vec.len - nrow)/step + 1;
	}
	
	inline T operator[] (int i) const {return ptr[i];}
	inline T& operator[] (int i){return ptr[i];}
	
	inline T operator() (int row, int col) const {return *(row + colptr(col));}
	inline T& operator() (int row, int col) {return *(row + colptr(col));}
	
	SWMat subsetCol(int colStart, int colEnd){
		return SWMat(colptr(colStart), nrow, colEnd-colStart, step);
	}
	
	
	inline T get(int row, int col){
		return ptr[row + colptr(col)];
	}
	
	inline Vec<T> getCol(int col){
		return Vec<T>(colptr(col), nrow);
	}
};


//matrix where columns correspond to non necessarily adjacent areas of memory 
//the only difference between this and a pointer of pointers is that here 
//an int per column is stored (the offset from the start of the matrix), instead of a 64 bits pointer.
template<typename T>
struct GapMat {
  T* ptr;
  int* colset;
  int nrow;
  int ncol;
  
  
  inline T* colptr(int col){
    return ptr + colset[col];
  }
  
  GapMat(){}
  
  GapMat(T* _ptr, int* _colset, int _nrow, int _ncol):
    ptr(_ptr), colset(_colset), nrow(_nrow), ncol(_ncol){}
  
  inline T operator() (int row, int col) const {return *(row + colptr(col));}
  inline T& operator() (int row, int col) {return *(row + colptr(col));}
  
  GapMat subsetCol(int colStart, int colEnd){
    return GapMat(ptr, colset + colStart, nrow, colEnd-colStart);
  }
  
  
  inline T get(int row, int col){
    //return ptr[row + colptr(col)];
    return colptr(col)[row];
  }
  
  inline Vec<T> getCol(int col){
    return Vec<T>(colptr(col), nrow);
  }
};


/* from std and rcpp datastructures extract pointers and wrap them */

//so far I am using it only to translate at compile-time INTSXP to int and REALSXP to double
#define CType(RType) typename Rcpp::traits::storage_type<RType>::type
//the opposite is:  Rcpp::traits::r_sexptype_traits<T>::rtype

//this is not going to work with vector<bool>,
//because it doesn't contain a pointer to bools
template<typename T>
inline Vec<T> asVec(std::vector<T>& v){
	return Vec<T>(v.data(), v.size());
}

template<int RType>
inline Vec< CType(RType) > asVec(Rcpp::Vector<RType>& v){
	return Vec< CType(RType) >(v.begin(), v.length());
}

template<typename T>
inline Mat<T> asMat(std::vector<T>& v, int ncol){
	if (v.size() % ncol != 0) throw std::invalid_argument("number of columns must be a divisor of vector length");
	return Mat<T>(v.data(), v.size()/ncol, ncol);
}

template<int RType>
inline Mat< CType(RType) > asMat(Rcpp::Matrix<RType>& m){
	return Mat< CType(RType) >(m.begin(), m.nrow(), m.ncol());
}


//these converters start from a SEXP, you must make sure that the SEXP survives
template<typename T>
inline GapMat<T> asGapMat(SEXP x){
	if (!Rf_inherits(x, "gapmat")) Rcpp::stop("the given object does not inherit from gapmat");
	Rcpp::List list = Rcpp::as<Rcpp::List>(x);
	//you must make sure that no memory is allocated here. The pointers won't be valid otherwise
	Vec<T> vec((SEXP)list["vec"]);
	Vec<int> colset((SEXP)list["colset"]);
	int nrow = list["nrow"];
	
	return GapMat<T>(vec.ptr, colset.ptr, nrow, colset.len);
}

template<typename T>
inline SWMat<T> asSWMat(SEXP x){
	if (!Rf_inherits(x, "swmat")) Rcpp::stop("the given object does not inherit from swmat");
	Rcpp::List list = Rcpp::as<Rcpp::List>(x);
	//you must make sure that no memory is allocated here. The pointers won't be valid otherwise
	Vec<T> vec((SEXP)list["vec"]);
	int step = list["step"];
	int nrow = list["nrow"];
	
	return SWMat<T>(vec, nrow, step);
}

/* colsums and rowsums */

template<typename TNumMat, typename TNumVec, template <typename> class TMat>
static void colSums(TMat<TNumMat> mat, Vec<TNumVec> vec, int nthreads){
	if (mat.ncol != vec.len) throw std::invalid_argument("provided vector has invalid length");

	TNumVec*  cs = vec.ptr;
	int nrow = mat.nrow;
	int ncol = mat.ncol;
	
	#pragma omp parallel for schedule(static) num_threads(std::max(1, nthreads))
	for (int col = 0; col < ncol; ++col){
		TNumMat* ptr = mat.colptr(col);
		TNumMat tmp = 0;
		for (int row = 0; row < nrow; ++row){
			tmp += *ptr++;
		}
		cs[col] = tmp;
	}
}


//here vec must be clean at the beginning
template<typename TNumMat, typename TNumVec, template <typename> class TMat>
static void rowSums(TMat<TNumMat> mat, Vec<TNumVec> vec, int nthreads){
	if (mat.nrow != vec.len) throw std::invalid_argument("provided vector has invalid length");

	int nrow = mat.nrow;
	int ncol = mat.ncol;
	
	#pragma omp parallel num_threads(std::max(1, nthreads))
	{
		std::vector<TNumVec> acc(nrow, 0);
		TNumVec* accBegin = acc.data();
		#pragma omp for schedule(static) nowait
		for (int col = 0; col < ncol; ++col){
			TNumMat* matCol = mat.colptr(col);
			TNumVec* accIter = accBegin;
			for (int row = 0; row < nrow; ++row){//this loop should be unrolled...
				*accIter++ += *matCol++;
			}
		}
		#pragma omp critical
		{
			for (int row = 0; row < nrow; ++row){
				vec[row] += acc[row];
			}
		}
	}
}


//colsums specialized for sliding-window matrices
template<typename TNumMat, typename TNumVec>
static void colSums(SWMat<TNumMat> mat, Vec<TNumVec> vec, int nthreads){
	if (mat.ncol != vec.len) throw std::invalid_argument("provided vector has invalid length");

	nthreads = std::max(nthreads, 1);
	
	TNumVec*  cs = vec.ptr;
	const int nrow = mat.nrow;
	const int step = mat.step;
	//in this case the normal colSums is probably faster:
	if (step*4 > nrow){
		colSums<TNumMat, TNumVec, SWMat>(mat, vec, nthreads); return;
	}
	
	const int ncol = mat.ncol;
	const int iold = -step;
	const int inew = nrow - step;
	
	
	//splitting the computation manually for multithreading.
	//Even though it is really hard to think of a case where multithreading here would
	//make any difference...
	std::vector<int> breaks(nthreads + 1);
	double facc = 0, f = ncol/((double)nthreads);
	for (int i = 1; i < nthreads; ++i){ breaks[i] = round(facc); facc += f; }
	breaks[nthreads] = ncol;
	
	
	#pragma omp parallel for schedule(static) num_threads(nthreads)
	for (int chunk = 0; chunk < nthreads; ++chunk){
		int firstcol = breaks[chunk];
		int lastcol = breaks[chunk+1];
		if (lastcol > firstcol){
			//first iteration without recursion
			TNumMat* colValues = mat.colptr(firstcol);
			TNumMat tmp = 0;
			for (int row = 0; row < nrow; ++row){
				tmp += colValues[row];
			}
			cs[firstcol] = tmp;
			//other iterations with recursion
			
			//it would be much more readable to place this if inside the loop,
			//but there is no guarantee that the compiler puts it here
			// (or, even more efficient, before)
			// (or, even more efficient, before)
			if (step == 1){
				for (int col = firstcol + 1; col < lastcol; ++col){ 
					colValues += 1;
					tmp += colValues[inew] - colValues[-1];
					
					cs[col] = tmp;
				}
			} else if (step == 2){
				for (int col = firstcol + 1; col < lastcol; ++col){ 
					colValues += 2; 
					tmp += colValues[inew] + colValues[inew + 1] - colValues[-2] - colValues[-1];
					
					cs[col] = tmp;
				}
			} else {
				for (int col = firstcol + 1; col < lastcol; ++col){
					colValues += step;
					//subtracting old values
					for (int i = iold; i < 0; ++i){ tmp -= colValues[i]; }
					//adding new ones
					for (int i = inew; i < nrow; ++i){ tmp += colValues[i]; }
					
					cs[col] = tmp;
				}
			}
		}
	}
}

/* sorting columns of a matrix */

/*
struct colPtr {
	int* ptr;
	int colref;
	colPtr(int* _ptr, int _colref): ptr(_ptr), colref(_colref){}
};

struct colSorter{
	int nrow;
	
	inline bool operator() (const colPtr& cp1, const colPtr& cp2) {
		int* ptr1 = cp1.ptr;
		int* ptr2 = cp2.ptr;
		
		for (int i = 0; i < nrow; ++i){
			int diff = ptr2[i] - ptr1[i];
			if (diff!=0){
				return diff > 0;
			}
		}
		
		return false;
	}
	
	colSorter(int _nrow) : nrow(_nrow){}
};

static inline void orderColumns_core(Mat<int> mat, Vec<int> vec){
	int nc = mat.ncol;
	colSorter srtr(mat.nrow);
	colPtr foo(0,0);
	std::vector<colPtr> cols(nc, foo);
	for (int i = 0; i < nc; ++i){cols[i] = colPtr(mat.colptr(i), i);}
	std::sort(cols.begin(), cols.end(), srtr);
	for (int i = 0; i < nc; ++i){vec[i] = cols[i].colref+1;}
}
*/

/* getting the row index with the maximum element for every column */

static inline void pwhichmax_core(Mat<double> posteriors, Vec<int> clusters, int nthreads){
	int ncol = posteriors.ncol;
	int nrow = posteriors.nrow;
	#pragma omp parallel for num_threads(nthreads)
	for (int col = 0; col < ncol; ++col){
		double* postCol = posteriors.colptr(col);
		double max = *postCol;
		int whichmax = 1;
		for (int row = 1; row < nrow; ++row){
			++postCol;
			if (*postCol > max){
				max = *postCol;
				whichmax = row + 1;
			}
		}
		clusters[col] = whichmax;
	}
}


/*	starts and lengths represent the IRanges object passed in from R,
	so starts are 1-indexed. 
	n represents the number of threads that will operate on the ranges,
	the ranges need to be re-split and grouped so that we get n groups
	of ranges of approximately equal total span.
	newstarts and newlengths at the beginning are empty vectors, at the
	end they will represent an equivalent IRanges object, but intervals can
	be split in smaller adjacent ones, and starts become 0-indexed.
	breaks and scorebreaks at the beginning are empty vectors, at the end
	they represent the starts and end of each group.


static inline void splitForThreads(
										Vec<int> starts,
										Vec<int> lengths,
										int nthreads,
										std::vector<int>& newstarts,
										std::vector<int>& newlengths,
										std::vector<int>& breaks,
										std::vector<int>& spanbreaks){
	
	//compute tot span
	int totspan = 0;
	int nrng = starts.len;
	for (int r = 0; r < nrng; ++r){
		totspan += lengths[r];
	}
	
	//initialize loop variables
	int currspan = 0;//span processed so far
	int r = 0;//range index in the starts and lengths vectors
	int initsize = newstarts.size(); //should always be zero
	breaks.push_back(initsize);
	spanbreaks.push_back(currspan);
	// (currlen, currstart) == (starts[r], lengths[r]) or second half of the last splitted interval
	int currstart = starts[r];
	int currlen = lengths[r];
	double N = nthreads;
	
	//loop on the threads
	for (int n = 0; n < nthreads; ++n){
		int targetspan = round(totspan*((n+1)/N)); //targetspan for the new iteration
		
		while (currspan < targetspan){
			if (currspan + currlen <= targetspan){
				newstarts.push_back(currstart-1);//from 1-indexed to 0-indexed
				newlengths.push_back(currlen);
				currspan += currlen;
				if (r < nrng){
					++r;
					currstart = starts[r];
					currlen = lengths[r];
				} else break; //this should happen only when currspan == targetspan
			
			} else {
				int lenfirsthalf = targetspan - currspan;
				newstarts.push_back(currstart-1);//from 1-indexed to 0-indexed
				newlengths.push_back(lenfirsthalf);
				currstart += lenfirsthalf;
				currlen -= lenfirsthalf;//lensecondhalf
				currspan += lenfirsthalf;
			}
		}
		
		breaks.push_back(newstarts.size()-initsize);
		spanbreaks.push_back(currspan);
	}
}

*/
