/*
 * sampleTopicsForDoc_sLDA_DirMult.cpp
 * DESCRIPTION  -----------------------------------------------------
 *   MEX function for completing one Gibbs Sampler sweep 
 *     through a single document for Latent Dirichlet Allocation
 *   Assigns each word token in the document to a unique topic
 *     based on conditional probability given current sample of
 *           of all other topic assignments in the corpus
 *   Intended for use as the fast core of a Matlab-based Gibbs Sampler
 * METADATA  --------------------------------------------------------
 *   Author: Mike Hughes ( mhughes@cs.brown.edu )
 *   Date:  24 June 2011
 * INPUT   ------------------------------------------------------
 *    ws     :  (1xNd) Matlab row vector of vocab word IDs observed in this document
 *    y      :  (1xC)  Matlab row vector of observed/hidden real response variable
 *    zs     :  (1xNd) Matlab row vector of current topic assignments for this document
 *    eta    :  (KxC)  Regression weight vector
 *    lambda :  (1x1) scalar precision parameter [only used for regression]
 *    Ndk    :  (1xK) Matlab matrix storing document-topic counts
 *                 In Matlab,  Ndk(k)  yields # tokens in doc d asgnd to topic k
 *    Nkt    :  (KxV) Matlab matrix storing topic-term counts
 *                 In Matlab,  Nkt(k,t) yields # tokens asgnd to topic k in corpus
  *                In C/C++ mex code, access the same entry via Nkt[k+K*t]
 *    Nk     :  (1xK) Matlab row vector storing topic counts across corpus
 *                             Nk(k) yields #times topic k assigned to any token in corpus
 *    ALPHA  :  Scalar symmetric hyperparameter on Dirichlet document-topic distribution 
 *    BETA   :  Scalar symmetric hyperparameter on Dirichlet topic-word distribution
 *    SEED   :  Integer seed for random number generator
 * OUTPUTS    -------------------------------------------------------
 *    zs, Ndk, Nkt, Nk  returned with updated assignments
 * COMPILATION ------------------------------------------------------
 * mex -I<path/to>/SuperTopics/rndgen/ sampleTopicsForDoc_LDA_DirMult.cpp
 * ACKNOWLEDGEMENTS   -----------------------------------------------
 *   Makoto Matsumoto and Takuji Nishimura for excellent basic rand generator
 *   see  ( mersenneTwister2002.c ) for details
 */

#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

using namespace std;

#define pi (3.141592653589793)
#define EPS (.000000000001)  // 1e-12

/* Input Arguments */
#define    Terms_IN  prhs[0]
#define    y_IN      prhs[1]
#define    Topics_IN prhs[2]
#define    eta_IN    prhs[3] 
#define    lambda_IN prhs[4]
#define    Ndk_IN    prhs[5]
#define    Nkt_IN    prhs[6]
#define    Nk_IN     prhs[7]
#define    ALPHA_IN  prhs[8]
#define    BETA_IN   prhs[9]
#define    PERM_IN   prhs[10]
#define    RAND_IN   prhs[11]
#define    debug_IN  prhs[12]

/* Output Arguments */
#define    Topics_OUT    plhs[0]
#define    Ndk_OUT    plhs[1]
#define    Nkt_OUT    plhs[2]
#define    Nk_OUT    plhs[3]

extern void _main();

const int numInputArgs  = 13;
const int numOutputArgs = 4;


/* =====================  Function Declarations =====================  */
void sampleTopicsForDoc( mxArray*, mxArray*, mxArray*, \
                         mxArray*, double, \
                         mxArray*, mxArray*, mxArray*, \
                         double, double, int, int, int, \
                         mxArray*, mxArray*, int );

void randperm( int, int* );

/* =================================================================
 *                    GATEWAY FUNCTION
 * Syntax:

 * =================================================================
 */
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  //Declarations
  int Ntoken;
  int V, K;
  double ALPHA, BETA;
  double yd, lambda;

  // Check to see if we have the correct number of input args
  if (nrhs != numInputArgs) {
    mexPrintf("Incorrect number of input arguments. Required %d, got %d\n", numInputArgs, nrhs );
    mexErrMsgTxt("Exiting...\n");
  }

  // --------------------   parse constants  ------------------------------
  ALPHA = mxGetScalar(ALPHA_IN);
  BETA = mxGetScalar(BETA_IN); 
  
  lambda = mxGetScalar(lambda_IN);

  // --------------------   create Ndk,Nkt,Nk ------------------------------
  //  remember define statements up top set things up so
  //    Topics_OUT => plhs[0], etc.
  Topics_OUT = mxDuplicateArray( Topics_IN );
  Ndk_OUT = mxDuplicateArray( Ndk_IN );
  Nkt_OUT = mxDuplicateArray( Nkt_IN );
  Nk_OUT = mxDuplicateArray( Nk_IN );

  V = (int) mxGetN( Nkt_IN );
  K = (int) mxGetM( Nkt_IN);  
  Ntoken = (int) mxGetN( Terms_IN);

  int debugFlag = (int) mxGetScalar( debug_IN );

  sampleTopicsForDoc( (mxArray *)Terms_IN, (mxArray *) y_IN, Topics_OUT, \
                       (mxArray *)eta_IN, lambda, \
                       Ndk_OUT, Nkt_OUT, Nk_OUT, \
                       ALPHA, BETA, V, K, Ntoken, \
                       (mxArray *)PERM_IN, (mxArray *)RAND_IN, debugFlag );
}




/* =================================================================
                     GIBBS SAMPLER MAIN LOOP
  Purpose:
    Reassigns tokens in given document to potentially new topics
      via one-at-a-time Gibbs sampling via conditional distribution for
         token i's topic z_i when fixing other existing assignments z_{not i}
      p( z_i = k | w_i=t, z_{not i}, ... ) = 
                 (Ndk + ALPH) * ( Nkt + BETA )/(Nk + BETA*V)
            where all counts exclude token i                                      
            see eq. 5 in Griffiths and Steyvers (2004) for math details
    Modifies Ndk,Nkt,Nk counts in place
  Error checking:
    Several optional checks in place to ensure indices fall in proper ranges
       to avoid Segmentation Faults. Rarely needed once familiar with code.
  =================================================================
*/
void sampleTopicsForDoc( mxArray *Terms, mxArray *yArr, mxArray *Topics, \
                 mxArray *etaArr, double lambda, \
                 mxArray *NdkArr, mxArray *NktArr, mxArray *NkArr, \
                 double ALPHA, double BETA, \
                 int V, int K, int Nd, \
                 mxArray *permArr, mxArray *randArr, int doDebug) {

  double *ts =  mxGetPr( Terms );
  double *zs =  mxGetPr( Topics );
  double *Ndk = mxGetPr( NdkArr );
  double *Nkt = mxGetPr( NktArr );
  double *Nk = mxGetPr( NkArr );
  double *y  = mxGetPr( yArr );
  
  double *perm = mxGetPr( permArr );
  double *randChoices = mxGetPr( randArr );

  int C = mxGetN(  etaArr );  
  double *eta = mxGetPr( etaArr );

  double BETASUM = (double) (V*BETA);
  double *ps = new double[K];
  
  double mu = 0;
  double *muvec = new double[C];
  if (C == 1) {
    for (int kk=0; kk<K; kk++) {
      mu += eta[kk]*Ndk[kk];
    }
  }
  else {
    for (int cc=0; cc<C; cc++) {
      muvec[cc] = 0;
      for (int kk=0; kk<K; kk++) {
        muvec[cc] += eta[kk+K*cc]*Ndk[kk];
      }
    }
  }
  
  for (int nn=0; nn<Nd; nn++) {
    int n = (int) perm[nn]-1; // adjust for zero-based idx here in C++ land    
    int t = (int) ts[n]-1;  
    int k = (int) zs[n]-1;

    if ( k < -1 || k >= K ) {
      mexPrintf(  "  ERROR: k out of bounds.  zs[%d] =%d\n", n, k);
      mexErrMsgTxt( "   somethings not right here.");
      return;
    }
    if ( t < 0 || t >= V ) {
      mexPrintf(  "  ERROR: t out of bounds.  ts[%d]=%d\n", n, t);
      mexErrMsgTxt( "   somethings not right here.");
      return;
    }
    
    // --------------  decrement counts at token
    if ( k >= 0 ) {
      Ndk[k] -= 1;
      Nkt[k+K*t ] -= 1;
      Nk[k] -= 1;
      if (C==1) {
        mu = mu - eta[k];
      } else {
        for (int cc=0; cc<C; cc++) {
          muvec[cc] = muvec[cc] - eta[k+K*cc];
        }
      }
    }
      
    // --------------  generate probs for this token
    double total=0;
    double diff;
    for (int kk=0; kk<K; kk++) {        
      ps[kk] = (Ndk[kk] + ALPHA)*( Nkt[kk+K*t] + BETA ) / (Nk[kk] + BETASUM);
      if (C==1) {
        diff = y[0]  - ( (mu+eta[kk])/Nd );
        ps[kk] *= exp( -0.5*lambda* diff*diff );
      } else {
        diff = 0;
        for (int cc=0; cc<C; cc++) {
          diff += y[cc] - ( (muvec[cc] + eta[kk+K*cc])/Nd );   
        }
        ps[kk] *= exp( -0.5*diff*diff );
      }
      total += ps[kk];
    }
    
    // --------------  choose new topic for this token
    double r = total * randChoices[n];
    double cursum = ps[0];
    int newk = 0;
    while ( r >= cursum && newk < K-1) {
      newk++;
      cursum += ps[newk];
    }
    if ( newk < 0 || newk >= K ) {
      mexPrintf(  "  ERROR: newk out of bounds | zs[%d]=%d\n", n, newk);
      mexErrMsgTxt( "   somethings not right here.");
      return;
    }
    // --------------  increment counts for token
    zs[n] = newk+1; // return to 1-based idx in Matlab land
    Ndk[newk] += 1;
    Nkt[newk+K*t ] += 1;
    Nk[newk] += 1;
    if (C==1) {
      mu = mu + eta[newk];
    } else {
      for (int cc=0; cc<C; cc++) {
        muvec[cc] = muvec[cc] + eta[newk+K*cc];
      }
    }
    
    
    if (doDebug && (nn < 3 || nn > Nd-4) ) {
      for (int kk=0; kk<K; kk++) {
        ps[kk] = ps[kk]/total;
      }
      mexPrintf( "%2d %.3f | %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f--> %d \n", n, randChoices[n],\
                       ps[0], ps[1], ps[2], ps[3], ps[4], ps[5], ps[6], ps[7], ps[8], ps[9], newk+1 );
    }
    
  }
  
 
  delete [] ps;
  ps = NULL;
  delete [] muvec;
  muvec = NULL;
}
