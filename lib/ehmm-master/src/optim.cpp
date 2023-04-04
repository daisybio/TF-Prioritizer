#include <cfloat> 
#include <cmath> 
#include <string>
#include <stdexcept> 
#include <iostream>

#include <Rcpp.h>
#include <R_ext/Applic.h>

/*
	* Just a wrapper around the optim function call in R
	* This is a lot less general:
	* 1. it uses only the "L-BFGS-B" optim method
	* 2. it uses only functions of one variable
	* 3. it does not use derivative information
	* In case of need, points 2 and 3 can be easily extended.
	* 
	* 
	* I abandoned the L-BFGS-B approach because very often it didn't converge.
	* Both from R and from C sometimes it starts cycling through the same values.
	* It could be that it is not well suited to one-dimensional problems,
	* (as the optim help message says) or that it is just buggy.
*/

/* Type of the optimization function:
	
	The first argument is the number of parameters in the second argument (always one here...).
	The second argument is the point where the function is evaluated.
	The third argument is a pointer passed down from the calling routine, normally used to carry auxiliary information. 
	The return value is the value of the function.
*/
typedef double optimfn(int, double *, void *);

/* Type of the gradient of the optimization function:
	
	The first argument is the number of parameters in the second argument.
	The second argument is the point where the function is evaluated.
	The third argument at the end contains the computed gradient vector.
	The fourth argument is a pointer passed down from the calling routine, normally used to carry auxiliary information. 
*/
typedef void optimgr(int, double *, double *, void *);


/*
*  x is the starting parameter on entry and x the final parameter on exit
*  fn is the optimization function
*  gr is the gradient function
*  lb is a pointer to the lower bound, or 0 if there are no lower bounds
*  ub is a pointer to the upper bound, or 0 if there are no upper bounds
*  fx at the end contains the attained minimum
*  ex carries external information for the function fn
*/


static inline void lbfgsb_wrapper(optimfn fn, optimgr gr, void* ex, double* x, double* fx, double* lb, double* ub){
	
	//dimension of the parameter space
	int n = 1; 
	/*
	  'lmm' is an integer giving the number of BFGS updates retained in
          the '"L-BFGS-B"' method, It defaults to '5'.
	 */
	int lmm = 5;
	/*
	 nbd is an integer array of dimension n.
	 On entry nbd represents the type of bounds imposed on the
	   variables, and must be specified as follows:
	   nbd(i)=0 if x(i) is unbounded,
		  1 if x(i) has only a lower bound,
		  2 if x(i) has both lower and upper bounds, and
		  3 if x(i) has only an upper bound.
	 On exit nbd is unchanged.

	*/
	int nbd = 0;
	if (lb != 0){
		if (ub != 0){
			nbd = 2;
		} else {
			nbd = 1;
		}
	} else if (ub != 0){
		nbd = 3;
	}
	
	
	/* did it fail? */
	int fail = 0;
	/*
	 * 'factr' controls the convergence of the '"L-BFGS-B"' method.
      Convergence occurs when the reduction in the objective is
      within this factor of the machine tolerance. Default is
      '1e7', that is a tolerance of about '1e-8'.
   */
	double factr = 1e+07;
	/*
	 * 'pgtol' helps control the convergence of the '"L-BFGS-B"' method.
      It is a tolerance on the projected gradient in the current
      search direction. This defaults to zero, when the check is
      suppressed.
   */
	double pgtol =  0;
	/*
	 * 'trace' Non-negative integer. If positive, tracing information on
      the progress of the optimization is produced. Higher values
      may produce more tracing information: for method '"L-BFGS-B"'
      there are six levels of tracing.  (To understand exactly what
      these do see the source code: higher levels give more
      detail.)
   */
   int trace = 0;
   //count how many times the optimization function was called
   int fncount = 0;
   //count how many times the gradient of the optimization function was called
	int grcount = 0;
	//maximum allowed number of iterations
	int maxit = 100;
	/*
	 did't quite get what the 'msg' argument is, it could be this:
	 'task' is a working string of characters of length 60 indicating
	 the current job when entering and quitting this subroutine.
	 */
	char msg[60];
	/*
	 'REPORT' The frequency of reports for the '"BFGS"', '"L-BFGS-B"'
	 and '"SANN"' methods if 'control$trace' is positive. Defaults
	 to every 10 iterations for '"BFGS"' and '"L-BFGS-B"', or
	 every 100 temperatures for '"SANN"'.
	 */
	int nREPORT = 10;
	
	//finally call the optimization function
	lbfgsb(n, lmm, x, lb, ub, &nbd, fx, fn, gr, &fail, ex, factr,
			pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
	
	
	/*sometime it fails... nothing to do about it....
	if (fail){
		throw std::invalid_argument("The L-BFGS-B algorithm failed with this error message:\n" + std::string(msg));
	}
	*/
}

/*
static inline void bfgs_wrapper(optimfn fn, optimgr gr, void* ex, double* x, double* fx){
	int n = 1; 
	int trace = 0;
	//count how many times the optimization function was called
	int fncount = 0;
	//count how many times the gradient of the optimization function was called
	int grcount = 0;
	//maximum allowed number of iterations
	int maxit = 100;
	int nREPORT = 10;
	//R's defaults 
	double reltol = sqrt(DBL_EPSILON);
	double abstol = -INFINITY;
	//That's what R does in optim.c
	int mask = 1;
	// did it fail? 
	int fail = 0;
	
	vmmin(n, x, fx, fn, gr, maxit, trace, &mask, abstol, reltol, nREPORT, ex, &fncount, &grcount, &fail);
	
	if (fail){
		throw std::invalid_argument("The BFGS algorithm did not converge within the given number of iterations");
	}
}
*/


/*
 * Given two points ax and bx it finds a bracketing triple ax, bx, cx
 * such that bx is within [ax, cx] (or [cx, ax]) and fa > fb < fc
 * copied from "numerical recipies in C: routine for initially bracketing a Minimum"
 * only modifications: 
 * 1. from float to double
 * 2. function func now takes the parameter info
 * 3. replaced the pointers with temp variables
 */
#define GOLD 1.618034
#define GLIMIT 100
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define FMAX(a,b) ((a) > (b) ? (a) : (b))
//it must hold that f(x), where x is the middle of the output braketing triple,
// is less or equal to both f(ax) and f(bx)
static void mnbrak(double* _ax, double* _bx, double* _cx, double* _fa, double* _fb, double* _fc, double (*func)(double, void*), void* info){
	double ax = *_ax;
	double bx = *_bx;
	double ulim, u, r, q, fu, dum, cx, fa, fb, fc;
	fa = (*func)(ax, info);
	fb = (*func)(bx, info);
	//makes sure that fb <= fa by swapping b with a
	if (fb > fa){
		SHFT(dum, ax, bx, dum);
		SHFT(dum, fb, fa, dum);
	}
	//first guess for c
	cx = bx + GOLD*(bx-ax);
	fc = (*func)(cx, info);
	while (fb > fc){//keep returning here until we bracket
		r = (bx-ax)*(fb-fc);
		q = (bx-cx)*(fb-fa);
		u = bx - ((bx-cx)*q - (bx-ax)*r)/(2.0*SIGN(FMAX(fabs(q-r), TINY), q-r));
		ulim = bx + GLIMIT*(cx - bx);//we won't go further than this
		//test various possibilities
		if ((bx-u)*(u-cx) > 0.0){ //parabolic u is between b and c, try it!
			fu = (*func)(u, info);
			if (fu < fc){//got a minimum between b and c
				*_ax = bx; *_fa = fb;
				*_bx = u; *_fb = fu;
				*_cx = cx; *_fc = fc;
				return;
			} else if (fu > fb){//got a minimum between a and u
				*_ax = ax; *_fa = fa;
				*_bx = bx; *_fb = fb;
				*_cx = u; *_fc = fu;
				return;
			}
			//parabolic fit was no use, use default magnification
			u = cx + GOLD*(cx-bx);
			fu = (*func)(u, info);
		} else if ((cx-u)*(u-ulim) > 0.0){//parabolic fit is between c and its allowed limit
			fu = (*func)(u, info);
			if (fu < fc){
				SHFT(bx, cx, u, cx + GOLD*(cx-bx))
				SHFT(fb, fc, fu, (*func)(u, info))
			}
		} else if ((u-ulim)*(ulim-cx) >= 0.0){//limit parabolic u to maximum allowed value
			u = ulim;
			fu = (*func)(u, info);
		} else {//reject parabolic u, use default magnification
			u = cx + GOLD*(cx-bx);
			fu = (*func)(u, info);
		}
		//Eliminate oldest point and continue
		SHFT(ax, bx, cx, u)
		SHFT(fa, fb, fc, fu)
	}
	//fb < fc, we found the bracketing triple
	*_ax = ax; *_fa = fa;
	*_bx = bx; *_fb = fb;
	*_cx = cx; *_fc = fc;
	return;
}

//adapted from the R's library
//it must hold: ax  < vx < bx
//I modified only the initializations so that it can use
//all the information gained in mnbrak: the three initial
//points and the three initial values.
//it must hold that f(x), where x is the output of the algorithm,
// is less or equal to all of f(ax), f(bx) and f(xx)
static double Brent_fmin(double ax, double xx, double bx, double fai, double fxi, double fbi, double (*f)(double, void *), void *info, double tol)
{
    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

/*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    a = ax;
    b = bx;
    //done by me starts here -------
    if (ax > xx || xx > bx){
		 throw std::invalid_argument("the three initial points must be in ascending order");
	 }
    x = xx; fx = fxi;
    if (fai > fbi){
		 w = b; fw = fbi;
		 v = a; fv = fai;
	 } else {
		 v = b; fv = fbi;
		 w = a; fw = fai;
	 }

    d = DBL_MAX;
    e = DBL_MAX;
    //done by me ends here -------
    tol3 = tol / 3.;

/*  main loop starts here ----------------------------------- */

    for(;;) {
	xm = (a + b) * .5;
	tol1 = eps * fabs(x) + tol3;
	t2 = tol1 * 2.;

	/* check stopping criterion */

	if (fabs(x - xm) <= t2 - (b - a) * .5) break;
	p = 0.;
	q = 0.;
	r = 0.;
	if (fabs(e) > tol1) { /* fit parabola */

	    r = (x - w) * (fx - fv);
	    q = (x - v) * (fx - fw);
	    p = (x - v) * q - (x - w) * r;
	    q = (q - r) * 2.;
	    if (q > 0.) p = -p; else q = -q;
	    r = e;
	    e = d;
	}

	if (fabs(p) >= fabs(q * .5 * r) ||
	    p <= q * (a - x) || p >= q * (b - x)) {/* a golden-section step */

	    if (x < xm) e = b - x; else e = a - x;
	    d = c * e;
	}
	else { /* a parabolic-interpolation step */
	    d = p / q;
	    u = x + d;

	    /* f must not be evaluated too close to ax or bx */

	    if (u - a < t2 || b - u < t2) {
		d = tol1;
		if (x >= xm) d = -d;
	    }
	}

	/* f must not be evaluated too close to x */

	if (fabs(d) >= tol1)
	    u = x + d;
	else if (d > 0.)
	    u = x + tol1;
	else
	    u = x - tol1;

	fu = (*f)(u, info);

	/*  update  a, b, v, w, and x */

	if (fu <= fx) {
	    if (u < x) b = x; else a = x;
	    v = w;    w = x;   x = u;
	    fv = fw; fw = fx; fx = fu;
	} else {
	    if (u < x) a = u; else b = u;
	    if (fu <= fw || w == x) {
		v = w; fv = fw;
		w = u; fw = fu;
	    } else if (fu <= fv || v == x || v == w) {
		v = u; fv = fu;
	    }
	}
    }
    /* end of main loop */

    return x;
}

//x1 here comes from initR
static double brent_wrapper(double x1, double x2, double (*f)(double, void *), void *info, double tol){
	//find a bracketing triple using x1 and x2 as starting points
	double x3, f1, f2, f3, tmp;
	
	/*
					double initScore = (*f)(x1, info);
	*/
	//std::cout << "initial pair:\t" << exp(x1) << "\t" << exp(x2) << std::endl;
	//std::cout << "tolerance:\t" << tol << std::endl;
	//std::cout << "braketing the minimum" << std::endl;
	mnbrak(&x1, &x2, &x3, &f1, &f2, &f3, f, info);
	//reorder points from smallest to biggest
	if (x3 < x1){
		SHFT(tmp, x1, x3, tmp)
		SHFT(tmp, f1, f3, tmp)
	}
	//std::cout << "braketing triple (x):\t" << exp(x1) << "\t" << exp(x2) << "\t" << exp(x3) << std::endl;
	//std::cout << "braketing triple f(x):\t" << f1 << "\t" << f2 << "\t" << f3 << std::endl;
	/*
				double F1 = (*f)(x1, info);
				double F2 = (*f)(x2, info);
				double F3 = (*f)(x3, info);
				
				if (F1 != f1 || F2 != f2 || F3 != f3){
					std::cout << "Something wrong in the bookkeeping" << std::endl;
				}
				if (F2 > F1 || F2 > F3){
					std::cout << "That's not a braketing triple..." << std::endl;
				}
				
				if (initScore < F2){
					std::cout << "That's a bad braketing triple... " << std::endl;
				}
	*/
	//std::cout << "calling brent" << std::endl;
	double ret = Brent_fmin(x1, x2, x3, f1, f2, f3, f, info, tol);
	//std::cout << "brent returned:\t" << exp(ret) << std::endl;
	/*
				double finalScore = (*f)(ret, info);
				if (finalScore > F2){
					std::cout << "Something wrong in brent" << std::endl;
				}
	*/
	
	return ret;
}
