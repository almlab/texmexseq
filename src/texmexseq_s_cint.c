/*
 * author: Scott Olesen <swo@mit.edu>
 * date: May 2014
 *
 * acknowledgement: This code borrowed and slightly refactored from the R package 'poilog'
 * written by Vidar Grotan and Steinar Engen.
 * 
 */

#include <math.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>
#include <Rinternals.h>

// integral utility functions
double argmax_log_integrand (int n, double mu, double sig);
double upper (int n, double mu, double sig, double argmax);
double lower (int n, double mu, double sig, double argmax);

// functions to populate individual slices in the integral
double my_f (double z, int n, double mu, double sig, double fac);

// functions to populate all the integral slices
void my_f_vec (double *z, int n_slices, void *p);

// pdfs for single observations
double poilog_singleton (int n, double mu, double sig);

// pdf wrappers for sequences of observations
void poilog (int *n, double *mu, double *sig, int *n_obs, double *val);


// -- utility functions --
double argmax_log_integrand (int n, double mu, double sig) {
  /* Find the maximum of the log of the integrand.
   *
   * n, mu, sig : poilog parameters
   *
   * returns : position of log(lambda) that maximizes the log of integrand
   */
  
  double d=100; // size of the next step
  double z=0; // log lambda
  
  // do bisection root-finding until step size is below abs tolerance
  while (d>0.00001) {
    if (n-1 - exp(z) - (z-mu)/pow(sig, 2)>0) {
      z=z+d;
    } else {
      z=z-d;
    }
    d /= 2;
  }
  
  return(z);
}

double upper (int n, double mu, double sig, double argmax) {
  /* Find an upper bound, where the integrand falls to some fold of its original value.
   *
   * n, mu, sig : poilog parameters
   * argmax : position of the maximum of the log of integrand
   *
   * returns : upper bound
   */
  
  double d=10; // size of the next step
  double z=argmax+20; // log lambda, starting above the maximum
  double mf; // value of the log integrand at the argmax
  
  mf = (n-1)*argmax - exp(argmax) - 0.5*pow((argmax-mu)/sig, 2);
  z = argmax+20;
  
  d = 10;
  while (d>0.000001) {
    if ((n-1)*z - exp(z) - 0.5*pow((z-mu)/sig, 2) - mf + log(1000000) > 0 ) {
      z=z+d;
    } else {
      z=z-d;
    }
    d=d/2;
  }
  
  return(z);
}

double lower (int n, double mu, double sig, double argmax) {
  double d;
  double z;
  double mf;
  
  mf = (n-1)*argmax - exp(argmax) - pow((argmax-mu)/sig, 2);
  z = argmax-20;
  
  d = 10;
  while (d>0.000001) {
    if ((n-1)*z - exp(z) - 0.5*pow((z-mu)/sig, 2) - mf + log(1000000) > 0 ) {
      z=z-d;
    } else {
      z=z+d;
    }
    d=d/2;
  }
  return(z);
}


// -- structs for parameters --
struct my_f_params {
  int n;
  double mu;  
  double sig;
  double fac;
};


double my_f (double z, int n, double mu, double sig, double fac) {
  // Create the slices in the poilog integral
  return exp(z*n - exp(z) - 0.5*pow((z-mu)/sig, 2) - fac);
}

void my_f_vec (double *z, int n_slices, void *p)
{
  // Fill the poilog integral slices
  int i;
  struct my_f_params *params = (struct my_f_params *)p;
  int n = (params->n);
  double sig = (params->sig);
  double mu = (params->mu);
  double fac = (params->fac);
  
  // populate the slices
  for (i=0; i<n_slices; i++) {
    z[i] = my_f(z[i], n, mu, sig, fac);
  }
  
  return;
}

// -- pdf functions --
double poilog_singleton (int n, double mu, double sig) {
  double lb, ub; // integral bounds
  double fac; // log(n!)
  double m; // position of the max of the log of the integrand
  double val; // output variable
  
  double result, abserr;
  int last, neval, ier;
  int lenw;
  int *iwork;
  double *work;
  int limit=100;
  double reltol=0.00001;
  double abstol=0.00001;
  lenw = 4*limit;
  iwork = (int *) Calloc(limit, int);
  work = (double *) Calloc(lenw, double);
  char message[80]; // error message

  m = argmax_log_integrand(n, mu, sig);
  lb = lower(n, mu, sig, m);
  ub = upper(n, mu, sig, m);
  fac = lgamma(n+1);
  
  struct my_f_params p = { n, mu, sig, fac };

  // numerically integrate using QAGS algorithm
  Rdqags(my_f_vec, (void *) &p, &lb, &ub,
         &abstol, &reltol, &result, &abserr, &neval, &ier,
         &limit, &lenw, &last, iwork, work);

  if (ier!=0) {
    sprintf(message, "error in integration (code %d)\n", ier);
    error(message);
  }

  // add the prefactor to the integration results
  val = result/(sqrt(2*M_PI)*sig);
  
  Free(iwork);
  Free(work);
  return(val);
}


// -- convenience wrapper --
void poilog (int *n, double *mu, double *sig, int *n_obs, double *val) {
  /* compute multiple bipoilog pdfs
   *
   * n, mu, sig : poilog parameters
   * n_obs : length of n
   * val : output variable
   */
  
  int i;
  for (i=0; i < *n_obs; i++){
    val[i] = poilog_singleton(n[i], *mu, *sig);
    R_CheckUserInterrupt();
  }
}