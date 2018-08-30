#include <catch.hpp>
#include <iostream>
#include <fenv.h>
#include "gsl/gsl_matrix.h"
#include <cblas.h>

#include <algorithm>
#include <limits>
#include <numeric>

#include "debug.h"
#include "mathfunc.h"
#include "fastblas.h"
#include "fastopenblas.h"

using namespace std;

TEST_CASE("calcPab", "[math]") {
// Hi_eval
// vector size: 49
  double Hi_eval[] = {1.0,0.00340397,0.00316875,0.00285475,0.000338947,0.000157865,0.000132286,9.72953e-05,8.30921e-05,7.85302e-05,7.50248e-05,6.77641e-05,6.74403e-05,5.97916e-05,3.66931e-05,3.56332e-05,3.26649e-05,2.94286e-05,2.78745e-05,2.57487e-05,2.34702e-05,2.23381e-05,2.17387e-05,2.11297e-05,1.99572e-05,1.92797e-05,1.79017e-05,1.77354e-05,1.76441e-05,1.69226e-05,1.6689e-05,1.62749e-05,1.56035e-05,1.53859e-05,1.50713e-05,1.46394e-05,1.43762e-05,1.35848e-05,1.29509e-05,9.6414e-06,9.07189e-06,8.36172e-06,7.36353e-06,7.28662e-06,6.93901e-06,4.85174e-06,3.93325e-06,3.44085e-06,1.68874e-06};
// Uab
// matrix size: 49 cols, 6 rows
  double Uab[] = {49,8.88476e-27,1.60576e-26,2.79763e-26,1.04467e-27,8.13814e-28,3.43994e-29,1.83964e-31,1.00874e-29,1.51195e-29,5.16275e-29,2.30463e-30,3.49053e-30,6.02949e-30,1.02633e-30,4.2725e-29,1.29895e-30,2.51337e-31,2.0631e-30,1.35086e-32,5.46109e-31,6.47883e-31,1.80278e-30,2.23342e-30,1.89715e-30,4.93038e-30,8.54899e-31,1.67604e-31,2.68166e-30,1.11242e-30,4.93038e-32,2.6366e-31,3.15544e-30,2.41589e-30,7.64402e-31,4.77751e-32,2.49601e-31,1.92593e-32,2.94045e-30,1.70175e-30,1.17174e-30,5.117e-30,1.30193e-31,5.69767e-30,7.70372e-32,8.49335e-32,3.64434e-31,9.53456e-32,1.05464e-30,// row 0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,// row 1
0.0829848,-8.68273e-17,7.90197e-16,5.28651e-17,-3.03831e-16,3.1215e-17,-4.35439e-17,3.44317e-18,-2.91753e-17,-1.41231e-17,5.76483e-17,8.20483e-18,-1.90286e-18,1.4879e-17,-2.41261e-18,-1.20099e-16,-1.50572e-17,1.22551e-17,-1.98056e-17,1.29155e-19,-1.93584e-17,-1.16457e-17,-1.01071e-17,4.18941e-18,-2.53464e-17,2.53682e-17,3.17125e-18,-1.39346e-17,3.12915e-17,2.82287e-18,-4.43288e-18,-2.77284e-18,6.938e-18,2.71779e-17,2.64918e-17,-4.61261e-18,-5.52614e-18,1.17298e-18,-6.61135e-18,6.76754e-17,-1.20832e-17,-4.93118e-17,-4.31897e-18,-2.83118e-17,-3.32353e-18,2.9943e-19,-6.21337e-18,-1.92731e-17,4.96634e-17,// row 2
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,// row 3
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,// row 4
0.00014054,8.48529e-07,3.88858e-05,9.98962e-08,8.83656e-05,1.19729e-06,5.51193e-05,6.44442e-05,8.43825e-05,1.31924e-05,6.43714e-05,2.92104e-05,1.03734e-06,3.67172e-05,5.67139e-06,0.000337593,0.00017454,0.000597554,0.000190132,1.23484e-06,0.000686213,0.000209331,5.66639e-05,7.85841e-06,0.000338635,0.000130526,1.17637e-05,0.00115853,0.000365131,7.1633e-06,0.000398558,2.91612e-05,1.52549e-05,0.000305741,0.000918125,0.00044534,0.000122349,7.14401e-05,1.4865e-05,0.00269132,0.000124604,0.000475211,0.000143276,0.000140681,0.000143383,1.05563e-06,0.000105934,0.00389584,0.00233867,// row 5
  };
// ab
// vector size: 6
  double ab[] = {6.91831e-310,6.91831e-310,0,0,0,0};
// Pab
// matrix size: 3 cols, 6 rows
  double Pab[] = {49,4.08372e-33,1.2538e-76,// row 0
0,1.72345e-47,4.65925e-33,// row 1
0.0829848,3.06246e-57,3.2526e-86,// row 2
0,0,1.36901e-71,// row 3
0,0,4.17633e-62,// row 4
0.000140915,-0.00014054,-0.00014054,// row 5
  };
}

TEST_CASE( "Math functions", "[math]" ) {
  debug_set_debug_mode(true);
  debug_set_no_check_mode(false);
  debug_set_strict_mode(true);
  double data[] = { 2,-1, 0,
                   -1, 2,-1,
                    0,-1, 2};
  gsl_matrix *m = gsl_matrix_alloc(3,3);
  copy(data, data+9, m->data);
  REQUIRE( isMatrixPositiveDefinite(m) );
  REQUIRE( isMatrixSymmetric(m) );
  // REQUIRE( checkMatrixEigen(m,0.001) );

  double data1[] = {1.0,0,0,
                    0,3.0,0,
                    0,0,2.0};
  copy(data1, data1+9, m->data);
  REQUIRE( isMatrixPositiveDefinite(m) );
  // REQUIRE( checkMatrixEigen(m) );

  double data2[] = {1,1,1,
                    1,1,1,
                    1,1,0.5};
  copy(data2, data2+9, m->data);
  REQUIRE( !isMatrixPositiveDefinite(m));
  // REQUIRE( !checkMatrixEigen(m) );

  double data3[] = {1.0,  0,  0,
                    3.0,3.0,  0,
                      0,  0,2.0};
  copy(data3, data3+9, m->data);
  REQUIRE( !isMatrixPositiveDefinite(m) );
  REQUIRE( !isMatrixSymmetric(m) );
  // REQUIRE( checkMatrixEigen(m) );

  // ---- NaN checks
  vector<double> v = {1.0, 2.0};
  REQUIRE (!std::isnan(std::accumulate(v.begin(), v.end(), 0)));
  vector<double> v2 = {1.0, 2.0, std::numeric_limits<double>::quiet_NaN()};
  REQUIRE (std::isnan(v2[2]));
  REQUIRE(has_nan(v2));
  // test minus nan
  vector<double> v3 = {1.0, 2.0, -std::numeric_limits<double>::quiet_NaN()};
  REQUIRE (std::isnan(v3[2]));
  REQUIRE(has_nan(v3));
}

TEST_CASE("cblas_dgemm", "[math]") {
   double *A, *B, *C;
   int m, n, k, i, j;
   double alpha, beta;

   printf ("\n This example computes real matrix C=alpha*A*B+beta*C using \n"
           " Intel(R) MKL function dgemm, where A, B, and  C are matrices and \n"
           " alpha and beta are double precision scalars\n\n");

   m = 2000, k = 200, n = 1000;
   printf (" Initializing data for matrix multiplication C=A*B for matrix \n"
           " A(%ix%i) and matrix B(%ix%i)\n\n", m, k, k, n);
   alpha = 1.0; beta = 0.0;

   printf (" Allocating memory for matrices aligned on 64-byte boundary for better \n"
           " performance \n\n");
   A = (double *)malloc( m*k*sizeof( double ));
   B = (double *)malloc( k*n*sizeof( double ));
   C = (double *)malloc( m*n*sizeof( double ));

   printf (" Intializing matrix data \n\n");
   for (i = 0; i < (m*k); i++) {
     A[i] = (double)(i+1);
   }

   for (i = 0; i < (k*n); i++) {
     B[i] = (double)(-i-1);
   }

   for (i = 0; i < (m*n); i++) {
     C[i] = 0.0;
   }

   printf (" Computing matrix product using Intel(R) MKL dgemm function via CBLAS interface \n\n");
   assert(m==2000);
   assert(k==200);
   assert(n==1000);
   //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
   //            m, n, k, alpha, A, k, B, n, beta, C, n);
   fast_cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    m, n, k, alpha, A, k, B, n, beta, C, n);

   REQUIRE(trunc(C[0]) == -2666620100.0 );
   REQUIRE(trunc(C[1]) == -2666640200.0 );
   REQUIRE(trunc(C[2003]) == -10627000400.0 );

}

TEST_CASE("fast_dgemm", "[math]") {
   double *A, *B, *C;
   int m, n, k, i, j;
   double alpha, beta;

   printf ("\n This example computes real matrix C=alpha*A*B+beta*C using \n"
           " Intel(R) MKL function dgemm, where A, B, and  C are matrices and \n"
           " alpha and beta are double precision scalars\n\n");

   m = 2000, k = 200, n = 1000;
   printf (" Initializing data for matrix multiplication C=A*B for matrix \n"
           " A(%ix%i) and matrix B(%ix%i)\n\n", m, k, k, n);
   alpha = 1.0; beta = 0.0;

   printf (" Allocating memory for matrices aligned on 64-byte boundary for better \n"
           " performance \n\n");
   A = (double *)malloc( m*k*sizeof( double ));
   B = (double *)malloc( k*n*sizeof( double ));
   C = (double *)malloc( m*n*sizeof( double ));

   printf (" Intializing matrix data \n\n");
   for (i = 0; i < (m*k); i++) {
     A[i] = (double)(i+1);
   }

   for (i = 0; i < (k*n); i++) {
     B[i] = (double)(-i-1);
   }

   for (i = 0; i < (m*n); i++) {
     C[i] = 0.0;
   }

   printf (" Computing matrix product using Intel(R) MKL dgemm function via CBLAS interface \n\n");
   // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
   //            m, n, k, alpha, A, k, B, n, beta, C, n);
   // eigenlib_dgemm(const char *TransA, const char *TransB, const double alpha,
   //                const gsl_matrix *A, const gsl_matrix *B, const double beta,
   //                gsl_matrix *C) {
   gsl_matrix *AM = gsl_matrix_safe_alloc(m,k); // rows x cols
   gsl_matrix *BM = gsl_matrix_safe_alloc(k,n);
   gsl_matrix *CM = gsl_matrix_calloc(m,n);

   fast_copy(AM,A);
   fast_copy(BM,B);
   fast_copy(CM,C);
   fast_dgemm("N","N",alpha,AM,BM,beta,CM);
   printf ("\n Computations completed.\n\n");
   A = AM->data;
   B = BM->data;
   C = CM->data;

   REQUIRE(trunc(C[0]) == -2666620100.0 );
   REQUIRE(trunc(C[1]) == -2666640200.0 );
   REQUIRE(trunc(C[2003]) == -10627000400.0 );

}

// The following code is normally disabled as it stops the unit tests!

#ifdef TEST_SIGFPE

#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>

void sighandler(int signum)
{
    printf("Process %d got signal %d\n", getpid(), signum);
    signal(signum, SIG_DFL);
    kill(getpid(), signum);
}

TEST_CASE("NaN handlers", "[math]") {

  signal(SIGFPE, sighandler);
  feenableexcept(FE_INVALID   |
                 FE_DIVBYZERO |
                 FE_OVERFLOW  |
                 FE_UNDERFLOW);
  double dirty = 0.0;
  double nanval= 0.0/dirty;
  printf("Succeeded! dirty=%lf, nanval=%lf\n",dirty,nanval);
}

#endif // TEST_SIGFPE
