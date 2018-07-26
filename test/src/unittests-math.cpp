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
