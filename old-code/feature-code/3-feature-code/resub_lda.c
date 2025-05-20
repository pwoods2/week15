/*
 * file : resub_lda.c
 *
 * summary : computes the resubstitution
 *   error estimate for LDA classifier
 *
 * function call :
 * 
 *   e = resub_lda(X1,X2,n1,n2,p,a,m);
 *     double e -> resubstitution error
 *     double *X1,*X2 -> data
 *     int n1,n2 -> number of observations in X1,X2
 *     int p -> dimension of observations
 *     double *a,m -> LDA coefficients, y = a^T x + m
 *
 *
 * author: Ulisses Braga-Neto
 *   
 * last revision: 02/26/2003 
 */

#include "err_est.h"

double resub_lda(double *X1,
		 double *X2,
		 int n1,
		 int n2,
		 int p,
		 double *a,
		 double m)
{
  int i,e=0;

  /* test on X1 */

  for (i=0; i<n1; i++)
    if (test_lda(X1+i*p,p,a,m)==2) e++;

  /* test on X2 */

  for (i=0; i<n2; i++)
    if (test_lda(X2+i*p,p,a,m)==1) e++;

  return (double) e/(n1+n2);
}
