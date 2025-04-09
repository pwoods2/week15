/*
 * file : bresub_lda.c
 *
 * summary : computes the bolstered resubstitution
 *   error estimate for LDA classifier
 *
 * function call :
 * 
 *   e = bresub_lda(X1,X2,n1,n2,p,a,m);
 *     double e -> bolstered resubstitution error
 *     double *X1,*X2 -> data
 *     int n1,n2 -> number of observations in X1,X2
 *     int p -> dimension of observations
 *     double *a,m -> LDA coefficients, y = a^T x + m
 *
 * author: Ulisses Braga-Neto
 *   
 * last revision: 07/10/2003 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "err_est.h"

double bresub_lda(double *X1,
		  double *X2,
		  int n1,
		  int n2,
		  int p,
		  double *a,
		  double m)
{
  register int i,j;
    double sig1,sig2,x,norm=0,v,e=0;
/*    double pi=4*atan(1);*/ // if needed, can be deactivated

  if (p>5) {
    fprintf(stderr,"Sorry, %d dimensions not supported.\n",p);
    exit(1);
  }

  /* correction factors versus p, 0.5/0.5 method */
  /* median of chi distribution with p degrees of freedom */

  double cp[] = {1,0.67448425,1.17741394,1.53816223,1.83213806,2.08601379};

  /* compute norm of a */

  for (i=0; i<p; i++)
    norm += a[i]*a[i];
  norm = sqrt(norm);

  /* get sig1 and sig2 */

  sig1 = mean_min_dist(X1,n1,p)/cp[p];
  sig2 = mean_min_dist(X2,n2,p)/cp[p];

  /* compute estimate */

  for (i=0; i<n1; i++) {
    v = m;
    for (j=0; j<p; j++)
      v += X1[i*p+j]*a[j];

    x = fabs(v)/(norm*sig1);
    
    if (v>0)
      e = e + ncdf(x); /* wrong side */
    else
      e = e + 1 - ncdf(x); /* right side */
  }

  for (i=0; i<n2; i++) {
    v = m;
    for (j=0; j<p; j++)
      v += X2[i*p+j]*a[j];

    x = fabs(v)/(norm*sig2);

    if (v<0)
      e = e + ncdf(x); /* wrong side */
    else
      e = e + 1 - ncdf(x); /* right side */
  }

  return (double) e/(n1+n2);
}
