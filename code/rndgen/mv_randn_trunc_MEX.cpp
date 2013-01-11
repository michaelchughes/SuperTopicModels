/*
 * sampleTopicsForDoc_LDA_DirMult.cpp
 DESCRIPTION  -----------------------------------------------------
   MEX function for sampling from a multivariate truncated normal distribution
       y ~ Normal( mu, eye ) restricted s.t.
            y[dim] is the largest entry of y
            e.g. y[dim] > y[d] for every d != dim
 METADATA  --------------------------------------------------------
   Author: Mike Hughes ( mhughes@cs.brown.edu )
   Date:  11 January 2013
 INPUT   ------------------------------------------------------
 
    SEED   :  Integer seed for random number generator
 OUTPUTS    -------------------------------------------------------
    ys     :  (Nx1) Matlab vector of 
 COMPILATION ------------------------------------------------------
   mex -I<path/to>/SuperTopics/rndgen/ sampleTopicsForDoc_LDA_DirMult.cpp
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
#define pi (3.141592653589793)
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

int Nburn = 3;

/* =================================================================
 *                    GATEWAY FUNCTION
 * Syntax:

 * =================================================================
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
    mexErrMsgTxt("Target dimension outside valid range.  Needs to be in 1...length(mu).\n");
  }


  int Nsamps = mxGetScalar( Nsamps_IN );
  int SEED   = mxGetScalar( SEED_IN );

  // Initialize random seed --------------------------------------------
  init_genrand( SEED );

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

	y_OUT = mxCreateDoubleMatrix( muR, Nsamps, mxREAL);
	memcpy( mxGetPr( y_OUT ), yKeep.data(), muR*Nsamps*sizeof(double));

}

double get_max_value_ignore_dim( VectorType& y, int targetDim) {
  int D = (int) y.size();
  double maxval = -1000.0;
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

double randn_trunc_above( double mu, double sigma, double b ) {
  if ( mu - b > 5*sigma  ) {
    return mu - sigma*randn_trunc_tail_rejection(  (mu-b)/sigma  );
  } else {
    double bbar = mynormcdf(  (b-mu)/sigma  );
    double u    = bbar*genrand_double();  // u ~ Unif(0,CDF(bbar)) 
    return mu + sigma*myinvnormcdf( u );
  }
}

double randn_trunc_below( double mu, double sigma, double a ) {
  if ( a - mu > 5*sigma  ) {
    return mu + sigma*randn_trunc_tail_rejection(  (a-mu)/sigma  );
  } else {
    double c = mynormcdf(  (a-mu)/sigma  );
    double u = c + (1-c)*genrand_double();  // u ~ Unif(CDF(abar),1) 
    return mu + sigma*myinvnormcdf( u );
  }
}

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
