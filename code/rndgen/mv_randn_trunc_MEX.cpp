/*
 * mv_randn_trunc_MEX.cpp
 DESCRIPTION  -----------------------------------------------------
   MEX function for sampling from a multivariate truncated normal distribution
       y ~ TruncNormal_d( mu, eye(D) ) is a D dimensional vector restricted s.t.
            y[d] is the largest entry of y
            i.e. y[d] >= y[d'] for every d' != d
 USAGE
    Matlab Syntax:
      Y = mv_randn_trunc_MEX( mu, dim, Nsamps, seed )
 INPUT   ------------------------------------------------------
    mu     :  Dx1 column vector
    dim    :  integer within [1,2,...D] indicates dimension for truncation
    Nsamps :  integer number of samples to draw from trunc. normal distr.
    SEED   :  integer seed for random number generator
 OUTPUTS    -------------------------------------------------------
    Y     :  (DxNsamps) matrix
                each column Y(:,n) is a sample from the trunc. normal. distr.
 METADATA  --------------------------------------------------------
   Author: Mike Hughes ( mhughes@cs.brown.edu )
   Date:  11 January 2013
 RUNTIME ADVICE
   Weird matlab errors even after successful compilation?
   Restart matlab with a pointer to the system's standard c++ library:
      >> LD_PRELOAD=/usr/lib/libstdc++.so.6 matlab
   This overrides the default matlab library, and allow
 COMPILATION ------------------------------------------------------
   mex -I<path/to>/boost/    \     (not needed in some cases)
       -I<path/to>/Eigen/    \
       mv_randn_trunc_MEX.cpp
   relies on mersenneTwister2002.c in same directory when mexified.
 DEPENDENCIES       -----------------------------------------------
   Eigen : allows matrix manipulations for this multivariate distribution
           see http://eigen.tuxfamily.org/
   Boost : trustworthy implementation of math special functions
               erf/erfc/inverf for use in approximating NormalCDF
           see http://www.boost.org/
 ACKNOWLEDGEMENTS   -----------------------------------------------
   Makoto Matsumoto and Takuji Nishimura for excellent basic rand generator
   see  ( mersenneTwister2002.c ) for details
 */

#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <boost/math/special_functions/erf.hpp>
#include "mersenneTwister2002.c"
#include "Eigen/Dense"

using namespace std;

#define SQRT2  1.41421356237309504
#define EPS (.000000000001)  // 1e-12

/* Input Arguments */
#define    mu_IN     prhs[0]
#define    dim_IN    prhs[1]
#define    Nsamps_IN prhs[2]
#define    SEED_IN   prhs[3]

/* Output Arguments */
#define    y_OUT    plhs[0]

extern void _main();

const int numInputArgs  = 4;
const int numOutputArgs = 1;

typedef Eigen::ArrayXXd MatrixType;
typedef Eigen::ArrayXd VectorType;

/* =====================  Function Declarations =====================  */
void gibbs_sample_mv_randn_trunc( VectorType&, const Eigen::Map<VectorType>&, int );

double get_max_value_ignore_dim(  VectorType&, int);

double randn_trunc_above( double, double, double );
double randn_trunc_below( double, double, double );

double randn_trunc_tail_rejection( double );

double myinvnormcdf( double );
double mynormcdf( double );

void randperm( int, int*);

int Nburn = 3; // # samples to discard as we wander from initial configuration

/* =================================================================
                     GATEWAY FUNCTION
   Matlab Syntax:
      Y = mv_randn_trunc_MEX( mu, 3, Nsamps, seed );
   Requires 4 input args, 1 output arg
   This inputs are given macro names like "y_OUT" and "mu_IN"
     in the definitions at the top of this code file
   =================================================================
 */
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Check to see if we have the correct number of input args
  if (nrhs != numInputArgs) {
    mexPrintf("Incorrect number of input arguments. Required %d, got %d\n", numInputArgs, nrhs );
    mexErrMsgTxt("Exiting...\n");
  }

  int muR = mxGetM( mu_IN );
  int muC = mxGetN( mu_IN );

  if ( muC != 1 ) {
    mexErrMsgTxt("Mu needs to be a column vector.\n");
  }

  int dim    = mxGetScalar( dim_IN );
  dim = dim -1; // go to zero-based in index here in C++ land
  if (dim < 0 || dim >= muR) {
    mexErrMsgTxt("Target dimension outside valid range.  Needs to be within {1,2,...length(mu)}.\n");
  }

  int Nsamps = mxGetScalar( Nsamps_IN );
  int SEED   = mxGetScalar( SEED_IN );
  init_genrand( SEED ); // call to setup our PRNG: mersenneTwister2002.c

	Eigen::Map<VectorType> mu( mxGetPr(mu_IN), muR);

	MatrixType yKeep = MatrixType::Zero(muR, Nsamps);
	VectorType yCur = VectorType::Zero(muR);

  int rkeep=0;
  for (int rr=0; rr<Nburn+Nsamps; rr++) {
    gibbs_sample_mv_randn_trunc( yCur,  mu, dim );
    
    if (rr >= Nburn) {
      yKeep.col(rkeep) = yCur;
      rkeep++;
    }
    
  }

  // Translate the Eigen matrix collecting our output samples 
  //  into the right matlab format
	y_OUT = mxCreateDoubleMatrix( muR, Nsamps, mxREAL);
	memcpy( mxGetPr( y_OUT ), yKeep.data(), muR*Nsamps*sizeof(double));

}

/* =================================================================
   COMPUTE THE MAXIMUM VALUE OF ENTRIES IN VECTOR y
     IGNORING PARTICULAR DIMENSION targetDim
   Routine function used by the Gibbs sampler to 
 */
double get_max_value_ignore_dim( VectorType& y, int targetDim) {
  int D = (int) y.size();
  double maxval = -1 * numeric_limits<double>::max(); //very small value
  for (int dd=0; dd < D; dd++) {
    if (dd == targetDim) {
      continue;
    }
    if ( y(dd) > maxval ) {
      maxval = y(dd);
    }
  }
  return maxval;
}

/* =================================================================
   GIBBS SAMPLER ALGORITHM TO DRAW A NEW VECTOR "y" ~ TruncNormal of dimension D
       GIVEN PREVIOUS VECTOR "y" AS INPUT
   Visits each entry in y in random order, 
     and resamples that entry using efficient 1D truncated normal sampling routines.
   When we visit the *target* dimension,
     we enforce that it must be larger than the maximum of all other entries
   Otherwise, we visit other dimensions and enforce that they must be smaller than y[target]
 */
void gibbs_sample_mv_randn_trunc( VectorType& y, const Eigen::Map<VectorType>& mu, int targetDim) {
  int D = (int) mu.size();
  double maxOthers;
  int *perm = new int[D];
  randperm( D, perm );

  for (int dd=0; dd < D; dd++) {
    int d = perm[dd];
    if (d == targetDim) {
      maxOthers = get_max_value_ignore_dim( y, targetDim );
      y(d) = randn_trunc_below( mu(d), 1.0, maxOthers );
    } else {
      y(d) = randn_trunc_above( mu(d), 1.0, y(targetDim) );
    }
  }
  delete [] perm;
  perm = NULL;
}

/* =================================================================
 * GENERIC SAMPLING FROM THE 1D NORMAL DISTRIBUTION
 *   TRUNCATED FROM ABOVE BY GIVEN BOUND 'b'
   Draws y ~ Norm( mu, sigma^2 )  s.t.   -Inf <= y <= b
   Uses numerical methods based on inverting normal CDF
     unless truncated region is so far (>5*sigma) away from the mean
     that these methods become unstable 
   In these extreme cases, uses a rejection sampler with assymptotically 
     perfect acceptance rates as distance from mean increases. 
 * See http://web.michaelchughes.com/research/sampling-from-truncated-normal
    and references therein
 */
double randn_trunc_above( double mu, double sigma, double b ) {
  if ( mu - b > 5*sigma  ) {
    return mu - sigma*randn_trunc_tail_rejection(  (mu-b)/sigma  );
  } else {
    double bbar = mynormcdf(  (b-mu)/sigma  );
    double u    = bbar*genrand_double();  // u ~ Unif(0,CDF(bbar)) 
    return mu + sigma*myinvnormcdf( u );
  }
}

/* =================================================================
 * GENERIC SAMPLING FROM THE 1D NORMAL DISTRIBUTION
 *   TRUNCATED FROM BELOW BY GIVEN BOUND 'a'
   Draws y ~ Norm( mu, sigma^2 )  s.t.   a <= y <= +Inf
   Uses numerical methods based on inverting normal CDF
     unless truncated region is so far (>5*sigma) away from the mean
     that these methods become unstable 
   In these extreme cases, uses a rejection sampler with assymptotically 
     perfect acceptance rates as distance from mean increases. 
 * See http://web.michaelchughes.com/research/sampling-from-truncated-normal
    and references therein
 */
double randn_trunc_below( double mu, double sigma, double a ) {
  if ( a - mu > 5*sigma  ) {
    return mu + sigma*randn_trunc_tail_rejection(  (a-mu)/sigma  );
  } else {
    double c = mynormcdf(  (a-mu)/sigma  );
    double u = c + (1-c)*genrand_double();  // u ~ Unif(CDF(abar),1) 
    return mu + sigma*myinvnormcdf( u );
  }
}

/* =================================================================
 * REJECTION SAMPLING FROM TAIL OF THE 1D NORMAL DISTRIBUTION
 * see http://web.michaelchughes.com/research/sampling-from-truncated-normal
    and references therein
 */
double randn_trunc_tail_rejection( double a ) {
 double u, w, x;
  do {
    u = genrand_double();
    w = genrand_double();
    x = sqrt( a*a - 2*log(u) );
  } while ( w*x > a );
  return x;
}

/* =================================================================
 *                   NORMAL CUMULATIVE DISTRIBUTION FUNCTION
 */
double mynormcdf( double x ) {
  return 0.5*erfc( - x/SQRT2 );
}

/* =================================================================
 *         INVERSE  NORMAL CUMULATIVE DISTRIBUTION FUNCTION
 */
double myinvnormcdf( double p ) {
  return SQRT2 * boost::math::erf_inv( 2*p -1 );
}


void randperm( int N, int* perm) {
  int i,j,t;
  for (i=0; i<N; i++)
    perm[i] = i;
  for (i=0; i<N; i++) {
    j = genrand_int31() % (N-i) + i;
    t = perm[j];
    perm[j] = perm[i];
    perm[i] = t;
  }
}
