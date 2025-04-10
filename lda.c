/*
 * file : lda.c
 *
 * summary : computes coefficients for
 *   linear discriminant analysis
 *
 * function call :
 * 
 *   lda(X1,X2,n1,n2,p,a,m);
 *     double *X1,*X2 -> data
 *     int n1,n2 -> number of observations in X1,X2
 *     int p -> dimension of observations
 *     double *a,*m -> LDA coefficients, y = a^T x + m (OUTPUT)
 *
 * author: Ulisses Braga-Neto
 *   
 * last revision: 03/31/2003 
 */

#include <stdlib.h>
#include <stdio.h>
#include "err_est.h"
#include <math.h>

void lda(double *X1, // int add by chen
    double *X2,
    int n1,
    int n2,
    int p,
    double *a,
    double *m)
{
  double *m1,*m2,*md,*ms,
         *S1,*S2,*S,*Si;
  double *at,*mt;
  double ratio = (double)n2/(double)n1;
  register int i;

  /* find means and covariances */

  m1 = mean(X1,n1,p);
  m2 = mean(X2,n2,p);
  S1 = cov(X1,n1,p);
  S2 = cov(X2,n2,p);

  /* pooled covariance */

  dmem(S,p*p);
//  for (i=0; i<p*p; i++)
//    S[i] = 0.5*(S1[i]+S2[i]);
  
  for (i=0; i<p*p; i++)
    S[i] =((n1-1)*S1[i]+(n2-1)*S2[i])/(n1+n2-2);  // changed by chen, 07/08/08

  /* find inverse of covariance*/
dmem(Si,p*p);
  Si = minv(S,p);

  /* compute LDA coefficients */

  dmem(md,p);
  dmem(ms,p);
  for (i=0; i<p; i++) {
    md[i] = m2[i] - m1[i];
    ms[i] = m1[i] + m2[i];
  }

  at = mmul(Si,md,p,p,1);
//  double at0 = *(at+0); //
//  double at1 = *(at+1); //
  mt = mmul(at,ms,1,p,1);

  memcpy(a,at,p*sizeof(double));
//  a = mmul(Si,md,p,p,1);
//  double a0 = *(a+0);//
//  double a1 = *(a+1);//
  *m = -0.5*(*mt)+log(ratio); // changed by chen to consider un-equal prior prob of two classes  Nov.23.2008
//  *m = *m + log(ratio);
  free(m1);free(m2);
  free(S1);free(S2);free(S);free(Si);
  free(md);free(ms);
  free(at);free(mt);
  
//  printf("X1[1]=%f\n",X1[1]);
//  printf("X2[1]=%f\n",X2[1]);
  

//  return 0;
}
