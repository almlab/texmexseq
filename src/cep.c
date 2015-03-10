/*
 * author: Scott Olesen <swo@mit.edu>
 * date: May 2014
 * 
 */

#include <math.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>
#include <Rinternals.h>

void cep(int *n1, int *n2, double *mu1, double *mu2, double *sig1, double *sig2, double *rho, int *trunc, int *n_obs, double *val)
{
  /* compute ceps
   *
   * n1 : a single integer
   * n2 : an array of n2 values, with length n_obs
   * val : output spot
   */
  
  double cumprob=0; // cumulative probability
  int x=0; // counter for the sum
  int n2j=0; // position of the next n2 value where cumprob will be dumped
  
  // for every n2, accumulate its probability
  for(x=0; x<=*n_obs; x++) {
    
  }
}