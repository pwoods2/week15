/*
 * file : test_lda.c
 *
 * summary : test an LDA classifier
 *   on a given data point
 *
 * function call :
 * 
 *   d = test_lda(x,p,a,m);
 *     int d -> class 1 or 2 
 *       (this is consistent with the LDA design in lda.c)
 *     double *x -> data point
 *     int p -> dimensionality
 *     double *a,*m -> LDA coefficients, y = a^T x + m
 *
 * author: Ulisses Braga-Neto
 *   
 * last revision: 02/26/2003 
 */

#include <stdlib.h>
#include <stdio.h>
#include "err_est.h"

int test_lda(double *x,
	     int p,
	     double *a,
	     double m)
{
  int d;
  double *l;

  l =  mmul(x,a,1,p,1);

  d = (*l+m <= 0) ? 1 : 2;
  
  free(l);
  
  return d;
}
