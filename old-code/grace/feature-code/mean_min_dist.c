/*
 * function: mean_min_dist
 *
 * summary : computes the mean of 
 *   the minimum distances from each point
 *   to the others
 *
 * function call :
 * 
 *   d = mean_min_dist(X,n,p);
 *     double d -> mean of minimum distances
 *     double *X -> data
 *     int n -> number of observations in X
 *     int p -> dimension of observations
 *
 * author: Ulisses Braga-Neto
 *   
 * last revision: 05/01/2003 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "err_est.h"

double mean_min_dist(double *X,
		     int n,
		     int p)
{
  register int i,j,k;
  double *min;
  double d,dm,v,mean;

  /* compute mimimum distances */

  dmem(min,n);

  for (i=0; i<n; i++) {

     dm = HUGE_VAL;

     for (j=0; j<i; j++) {
       for (d=k=0; k<p; k++) {
	 v = X[i*p+k]-X[j*p+k];
	 d += v*v;
       }
       if (d<dm) dm = d;
     }

     for (j=i+1; j<n; j++) {
       for (d=k=0; k<p; k++) {
	 v = X[i*p+k]-X[j*p+k];
	 d += v*v;
       }
       if (d<dm) dm = d;
     }

     min[i] = sqrt(dm);
  }

  /* find mean */
  dm = 0;
  for (i=0; i<n; i++)
    dm += min[i];
  mean = dm/n;
  
  free(min);

  /* return mean */

  return mean;
}
