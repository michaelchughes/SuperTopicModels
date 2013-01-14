/*
 * randn_trunc_MEX.cpp
 DESCRIPTION  -----------------------------------------------------
   MEX function for sampling from a 1D truncated normal distribution
      y ~ Normal( mean=mu, var=sigma^2 ) restricted s.t.
             y <=0  if doA=0
             y  >0  if doA=1      where y is a scalar real
 INPUT   ------------------------------------------------------
     mu    :  scalar real mean for the distribution
                 (before truncation)
    sigma  :  scalar positive *standard deviation* for distribution
                 (before truncation)
    doA    :  binary flag
                1 := enforce y <=0 (truncate from above)
                0 := enforce y > 0 (truncate from below)
    Nsamps :  positive integer, number of desired samples for y
    SEED   :  integer seed for random number generator
 OUTPUTS    -------------------------------------------------------
    ys     :  (Nx1) Matlab vector of indep. draws from 1D truncated normal
 METADATA  --------------------------------------------------------
   Author: Mike Hughes ( mhughes@cs.brown.edu )
   Date:  11 January 2013
 COMPILATION ------------------------------------------------------
   mex -I<path/to>/boost/    \
        randn_trunc_MEX.cpp
   relies on mersenneTwister2002.c in same directory when mexified.
 REFERENCES          ----------------------------------------------
   http://web.michaelchughes.com/research/sampling-from-truncated-normal
 DEPENDENCIES       -----------------------------------------------
   Boost : trustworthy implementation of math special functions
               erf/erfc/inverf for use in approximating NormalCDF
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

using namespace std;

#define SQRT2  1.41421356237309504
#define pi (3.141592653589793)
#define EPS (.000000000001)  // 1e-12

/* Input Arguments */
#define    mu_IN     prhs[0]
#define    sigma_IN  prhs[1]
#define    doA_IN    prhs[2]
#define    Nsamps_IN prhs[3]
#define    SEED_IN   prhs[4]

/* Output Arguments */
#define    y_OUT    plhs[0]

extern void _main();

const int numInputArgs  = 5;
const int numOutputArgs = 1;


/* =====================  Function Declarations =====================  */
double randn_trunc_above( double, double, double );
double randn_trunc_below( double, double, double );

double randn_trunc_tail_rejection( double );

double myinvnormcdf( double );
double mynormcdf( double );


/* =================================================================
 *                    GATEWAY FUNCTION
 * Matlab Syntax:
      y ~ randn_trunc_MEX( mu, sigma, doA, Nsamps, seed )
   Requires 5 input args, 1 output arg
   This inputs are given macro names like "y_OUT" and "mu_IN"
     in the definitions at the top of this code file
     , which are much more readable than prhs[0] or plhs[0]
 * =================================================================
 */
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Check to see if we have the correct number of input args
  if (nrhs != numInputArgs) {
    mexPrintf("Incorrect number of input arguments. Required %d, got %d\n", numInputArgs, nrhs );
    mexErrMsgTxt("Exiting...\n");
  }

  double mu  = mxGetScalar( mu_IN );
  double sigma = mxGetScalar( sigma_IN );
  int doA    = mxGetScalar( doA_IN );
  int Nsamps = mxGetScalar( Nsamps_IN );
  int SEED   = mxGetScalar( SEED_IN );

  init_genrand( SEED ); // call to setup our PRNG: mersenneTwister2002.c

  y_OUT = mxCreateDoubleMatrix(Nsamps,1,mxREAL);
  double *ys = mxGetPr( y_OUT );

  if ( doA == 0 ) {
    for (int nn=0; nn<Nsamps; nn++) {      
      ys[nn] = randn_trunc_above( mu, sigma, 0 );
    }
  } else {
    for (int nn=0; nn<Nsamps; nn++) {      
      ys[nn] = randn_trunc_below( mu, sigma, 0 ); 
    }
  }
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
