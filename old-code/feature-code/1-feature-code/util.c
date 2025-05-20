/*
 * file : util.c
 *
 * summary : utility functions
 *   - mean vector
 * 
 *   - covariance matrix
 * 
 *   - normal CDF (translation of the MatLab
 *     implementation of an algorithm by 
 *     W. J. Cody; see his paper "Rational Chebyshev
 *     approximations for the error function",
 *     Math. Comp., 1969, pp. 631-638)
 *
 *   - matrix multiplication
 * 
 *   - matrix inversion (via adjoint formula;
 *     see Horn and Johnson, "Matrix Analysis,"
 *     Cambridge University Press, 1985)
 * 
 *   - determinant (via Laplace expansion;
 *     see Horn and Johnson's book) 
 *
 *   - random permutation (via an in-place
 *     algorithm, described in Ross,
 *     "A First Course in Probability,"
 *     Macmillan, 1994)
 *
 *   - generation of (uncorrelated) normal
 *     random points (via the polar method;
 *     see Ross' book)
 *
 * function calls :
 * 
 *   m = mean(X,n,p);
 *     double *m -> mean vector (px1)
 *     double *X -> array of input vectors (nxp)
 *     int n -> number of vectors
 *     int p -> vector dimension
 *
 *   S = cov(X,n,p);
 *     double *S -> covaraince matrix (pxp)
 *     double *X -> array of input vectors (nxp)
 *     int n -> number of vectors
 *     int p -> vector dimension
 * 
 *   e = ncdf(double x);
 *     double e -> normal CDF (0<=e<=1)
 *     double x -> real value
 * 
 *   M = mmul(A,B,n,p,q);
 *     double *M -> output matrix (nxq)
 *     double *A -> input matrix (nxp)
 *     double *B -> input matrix (pxq)
 *     int n,p,q -> matrix dimensions
 *
 *   I = minv(A,n);
 *     double *I -> output matrix (nxn)
 *     double *A -> input matrix (nxn)
 *     int n -> matrix dimension
 *     
 *   d = det(A,n);
 *     double d -> determinant
 *     double *A -> input matrix (nxn)
 *     int n -> matrix dimension
 *
 *   pmvec(I,n);
 *     int *I -> input AND output vector
 *     int n -> vector length
 *
 *   X = norm_points(m,s,n,p);
 *     double *X -> samples
 *     double *m -> mean vector
 *     double *s -> std vector 
 *     int n -> number of samples 
 *              (must be even)
 *     int p -> sample dimension
 *  
 * 
 * S = double *autoCor(double *X,
	    int n,
	    int p)
 * 
 * author: Ulisses Braga-Neto
 *   
 * last revision: 01/15/2004 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "err_est.h"

double *mean(double *X,
	     int n,
	     int p)
{
  double *m;
  register int i,j;

  dmem(m,p);

  for (i=0; i<p; i++) {
    m[i] = 0;
    for (j=0; j<n; j++)
      m[i] += X[j*p+i];
    m[i] /= n;
  }
  
  return m;
}

double *cov(double *X,
	    int n,
	    int p)
{
  double *S,*St,*m,*mt;
  double k;
  register int i,j;

  dmem(S,p*p);
  dmem(mt,p);
 
  k = 1.0/(n-1); 
  m = mean(X,n,p);

  for (i=0; i<p*p; i++) S[i] = 0;
 
  for (i=0; i<n; i++) {
    for (j=0; j<p; j++)
      mt[j] = X[i*p+j] - m[j];
    St = mmul(mt,mt,p,1,p);
    for (j=0; j<p*p; j++)
      S[j] += k*St[j];
    free(St);
  }
  
  free(m);free(mt);

  return S;
}

double ncdf(double x)
{
  const double a[] = {3.16112374387056560e00,
		      1.13864154151050156e02,
		      3.77485237685302021e02,
		      3.20937758913846947e03,
		      1.85777706184603153e-1};
  const double b[] = {2.36012909523441209e01,
		      2.44024637934444173e02,
		      1.28261652607737228e03,
		      2.84423683343917062e03};
  const double c[] = {5.64188496988670089e-1,
		      8.88314979438837594e00,
		      6.61191906371416295e01,
		      2.98635138197400131e02,
		      8.81952221241769090e02,
		      1.71204761263407058e03,
		      2.05107837782607147e03,
		      1.23033935479799725e03,
		      2.15311535474403846e-8};
  const double d[] = {1.57449261107098347e01,
		      1.17693950891312499e02,
		      5.37181101862009858e02,
		      1.62138957456669019e03,
		      3.29079923573345963e03,
		      4.36261909014324716e03,
		      3.43936767414372164e03,
		      1.23033935480374942e03};
  const double p[] = {3.05326634961232344e-1,
		      3.60344899949804439e-1,
		      1.25781726111229246e-1,
		      1.60837851487422766e-2,
		      6.58749161529837803e-4,
		      1.63153871373020978e-2};
  const double q[] = {2.56852019228982242e00,
		      1.87295284992346047e00,
		      5.27905102951428412e-1,
		      6.05183413124413191e-2,
		      2.33520497626869185e-3};
  const double xbreak = 0.46875;
  const double pi = acos(-1);
  double y,z,xnum,xden,r,del,e;
  int i;

  x = x/sqrt(2);
  y = fabs(x);
  
  /* evaluate erf for |x|<=0.46875 */

  if (fabs(x) <= xbreak) {
    z = y*y;
    xnum = a[4]*z;
    xden = z;
    for (i=0;i<3;i++) {
      xnum = (xnum+a[i])*z;
      xden = (xden+b[i])*z;
    }
    r = x*(xnum+a[3])/(xden+b[3]);
  }

  /* evaluate erf for 0.46875<=|x|<=4.0 */
  else if ((fabs(x) > xbreak) && (fabs(x) <= 4.0)) {
    xnum = c[8]*y;
    xden = y;
    for (i=0; i<7; i++) {
      xnum = (xnum+c[i])*y;
      xden = (xden+d[i])*y;
    }
    r = (xnum+c[7])/(xden+d[7]);
    z = floor(y*16.0)/16.0;
    del = (y-z)*(y+z);
    r = exp(-z*z)*exp(-del)*r;
    if (x>0)
      r = 1-r;
    else
      r = r-1;
  }

  /* evaluate erf for |x|>4.0 */
  else {
    z = 1.0/(y*y);
    xnum = p[5]*z;
    xden = z;
    for (i=0; i<4; i++) {
      xnum = (xnum+p[i])*z;
      xden = (xden+q[i])*z;
    }
    r = z*(xnum+p[4])/(xden+q[4]);
    r = ((1/sqrt(pi))-r)/y;
    z = floor(y*16.0)/16.0;
    del = (y-z)*(y+z);
    r = exp(-z*z)*exp(-del)*r;
    if (x>0)
      r = 1-r;
    else
      r = r-1;
  }
  
  e = 0.5*(1+r);
  if (e>1) e=1;
  
  return e;
}

double *mmul(double *A,
	     double *B,
	     int n,
	     int p,
	     int q)
{
  register int i,j,k;
  double *M;

  dmem(M,n*q);

  for (i=0; i<n*q; i++) M[i] = 0;
 
  for (i=0; i<n; i++)
    for (j=0; j<q; j++)
      for (k=0; k<p; k++)
	M[i*q+j] += A[i*p+k]*B[k*q+j];

  return M;
}

double *minv(double *A, int n) {
	// added by chen 08/07/08 to deal with 1 demision case
	double *I, *It;
	double d, di;
	register int i, j, k, c;
	int s;
	
	
	dmem(I,n*n);
			
	if (n==1) {
		
		I[0] = 1/A[0];
		return I;
	}

	else {
		

		if ((d = det(A, n)) == 0) {
			fprintf(stderr,"minv: matrix is singular");
			exit(2);
		}
		di = 1/d;

		dmem(It,(n-1)*(n-1));

		for (i=0; i<n; i++)
			for (j=0; j<n; j++) {
				c = 0;
				for (k=0; k<n*n; k++)
					if ((k/n!=j) && (k%n!=i))
						It[c++] = A[k];
				if ((i+j)%2==0)
					s = 1;
				else
					s = -1;
				I[i*n+j] = s*di*det(It, n-1);
			}

		free(It);

		return I;
	}
}

double det(double *A,
	   int n)
{
  double *It;
  double d;
  register int j,k,c;
  int s;

  if (n==1)
    return A[0];
  else if (n==2)
    return A[0]*A[3]-A[1]*A[2];

  dmem(It,(n-1)*(n-1));

  d = 0;
  for (j=0; j<n; j++) {
    c = 0;
    for (k=n; k<n*n; k++)
      if (k%n!=j)
	It[c++] = A[k];
    if (j%2==0)
      s = 1;
    else
      s = -1;
    d += s*A[j]*det(It,n-1);
  }

  free(It);

  return d;
}

int pmvec(int *I,
	  int n)
{
  register int i;
  int k,tmp;

  for (i=n-1; i>=1; i--) {
    k = rnd_int(i+1);
    tmp = I[k];
    I[k] = I[i];
    I[i] = tmp;
  }

  return 0;
}

double *norm_points(double *m,
		    double *s,
		    int n,
		    int p)
{
  double u1,u2,r,*X;
  register int i,j;
  dmem(X,n*p);

  for (i=0; i<n/2; i++)
    for (j=0; j<p; j++) {
      
      do {

	/* uniform r.v.'s on [-1,1] */ 

	u1 = -1 + 2*(double)rand()/RAND_MAX;
 	u2 = -1 + 2*(double)rand()/RAND_MAX;
	
      } while ((r = u1*u1 + u2*u2) > 1);

      r = sqrt(-2*log(r)/r);
      X[i*p+j] = m[j] + s[j]*r*u1;
      X[(i+n/2)*p+j] = m[j] - s[j]*r*u2; 
      
    }

  return X;
}


double *autoCor(double *X,
	    int n,
	    int p)
{
  double *S,*St,*mt;
  double k;
  register int i,j;

  dmem(S,p*p);
  dmem(mt,p);
 
  k = 1.0/n; 
  

  for (i=0; i<p*p; i++) S[i] = 0;
 
  for (i=0; i<n; i++) {
    for (j=0; j<p; j++)
      mt[j] = X[i*p+j];
    St = mmul(mt,mt,p,1,p);
    for (j=0; j<p*p; j++)
      S[j] += k*St[j];
    free(St);
  }
  
  free(mt);

  return S;
}
