/*
 * file : err_est.h
 *
 * summary : data type definition,
 *           function prototypes
 *           and macro definitions.
 *
 * author: Ulisses Braga-Neto
 *   
 * last revision: 07/10/2003 
 */

/* type definition */
#include <string.h> // added by chen
typedef struct node {

  int leaf;           /* if nonzero, leaf node */
  int class;          /* class 1 or 2, if node is a leaf */
  int var;            /* variable to split on */ 
  double split;       /* split value */
  struct node *left;  /* left child */
  struct node *right; /* right child */

} tree;

/* classifier design and test functions */

void lda(double *,double *,int,int,int,double *,double *);
void SF(double *,double *,int,int,int,double *,double *);
tree *cart(double *,double *,int,int,int,int);
int test_lda(double *,int,double *,double);
int test_knn(double *,int,int,double *,double *,int,int);
int test_cart(double *, tree *);
int free_cart(tree *);
 
/* error functions */

double err_lda(int,double,double,double,double,double *,double);
double err_knn(int,double,double,double,double,double *,double *,int,int,int);
double err_cart(int,double,double,double,double,tree *);

double resub_lda(double *,double *,int,int,int,double *,double);
double resub_knn(double *,double *,int,int,int,int);
double resub_cart(double *,double *,int,int,int,tree *);

double loo_lda(double *,double *,int,int,int);
double loo_knn(double *,double *,int,int,int,int);
double loo_cart(double *,double *,int,int,int,int);

double cv_lda(double *,double *,int,int,int,int,int);
double cv_knn(double *,double *,int,int,int,int,int,int);
double cv_cart(double *,double *,int,int,int,int,int,int);

double boot_corr_lda(double *,double *,int,int,int,int,double);
double boot_corr_knn(double *,double *,int,int,int,int,int,double);
double boot_corr_cart(double *,double *,int,int,int,int,int,double);

double boot_632_lda(double *,double *,int,int,int,int,double);
double boot_632_knn(double *,double *,int,int,int,int,int,double);
double boot_632_cart(double *,double *,int,int,int,int,int,double);

double bresub_lda(double *,double *,int,int,int,double *,double);
double bresub_knn(double *,double *,int,int,int,int);
double bresub_cart(double *,double *,int,int,int,tree *);
double bresub_SF(double *,double *,int,int,int,double *,double);


double sresub_lda(double *,double *,int,int,int,double *,double);
double sresub_knn(double *,double *,int,int,int,int);
double sresub_cart(double *,double *,int,int,int,tree *);

double bloo_lda(double *,double *,int,int,int);
double bloo_knn(double *,double *,int,int,int,int);
double bloo_cart(double *,double *,int,int,int,int);

/* utility functions */

double *mean(double *,int,int);
double *cov(double *,int,int);
double ncdf(double);
double *mmul(double *,double *,int,int,int);
double *minv(double *,int);
double det(double *,int);
int    pmvec(int *, int);
double *norm_points(double *,double *,int,int);
double min_dist(double *,double *,int,int);
double mean_min_dist(double *,int,int);
double *autoCor(double *,int,int); // added by chen 08/07/2008


/* macros */

/* imem: allocates N ints at the address pointed to by P */

#define imem(P,N) \
    if ((P = (int *) malloc((N)*sizeof(int))) == NULL) { \
    fprintf(stderr,"imem: Memory allocation failed"); \
    exit(1); \
  }  

/* dmem: allocates N doubles at the address pointed to by P */

#define dmem(P,N) \
    if ((P = (double *) malloc((N)*sizeof(double))) == NULL) { \
    fprintf(stderr,"dmem: Memory allocation failed"); \
    exit(1); \
  }  

/* tmem: allocates one node of a CART pointed to by P */

#define tmem(P) \
    if ((P = (tree *) malloc(sizeof(tree))) == NULL) { \
    fprintf(stderr,"tmem: Memory allocation failed"); \
    exit(1); \
  }  

/* rnd_int: generates a random integer
   between 0 and N-1 (inclusive) */

#define rnd_int(N) \
  rand()%(N)
