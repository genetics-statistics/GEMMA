/*
    Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright © 2011-2017, Xiang Zhou
    Copyright © 2017, Peter Carbonetto
    Copyright © 2017-2020, Pjotr Prins

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>    // std::min
#include <cmath>
#include <iomanip>
#include <vector>
#include "debug.h"
#include "mathfunc.h"
#include <string.h>
#include "fastblas.h"

const char *FastblasTrans = "T";
const char *FastblasNoTrans = "N";

using namespace std;

/*
   Reasonably fast function to copy data from standard C array into
   gsl_matrix. Avoid it for performance critical sections.
*/
gsl_matrix *fast_copy(gsl_matrix *m, const double *mem) {
  auto rows = m->size1;
  auto cols = m->size2;
  if (is_strict_mode()) { // slower correct version
    for (size_t r=0; r<rows; r++) {
      for (size_t c=0; c<cols; c++) {
        gsl_matrix_set(m,r,c,mem[r*cols+c]);
      }
    }
  } else { // faster goes by row
    auto v = gsl_vector_calloc(cols);
    enforce(v); // just to be sure
    for (size_t r=0; r<rows; r++) {
      assert(v->size == cols);
      assert(v->block->size == cols);
      assert(v->stride == 1);
      memcpy(v->block->data,&mem[r*cols],cols*sizeof(double));
      gsl_matrix_set_row(m,r,v);
    }
    gsl_vector_free(v);
  }
  return m;
}

/*
    Helper function fast_cblas_dgemm runs the local dgemm
*/
void fast_cblas_dgemm(const enum CBLAS_ORDER Order,
                      const enum CBLAS_TRANSPOSE TransA,
                      const enum CBLAS_TRANSPOSE TransB,
                      const size_t M,
                      const size_t N,
                      const size_t K,
                      const double alpha,
                      const double *A,
                      const size_t lda,
                      const double *B,
                      const size_t ldb,
                      const double beta,
                      double *C,
                      const size_t ldc) {
#ifndef NDEBUG
  if (is_debug_mode()) {
    #ifdef DISABLED
    size_t i,j;
    printf (" Top left corner of matrix A: \n");
    for (i=0; i<min(M,6); i++) {
      for (j=0; j<min(K,6); j++) {
        printf ("%12.0f", A[j+i*K]);
      }
      printf ("\n");
    }

    printf ("\n Top left corner of matrix B: \n");
    for (i=0; i<min(K,6); i++) {
      for (j=0; j<min(N,6); j++) {
        printf ("%12.0f", B[j+i*N]);
      }
      printf ("\n");
    }

    printf ("\n Top left corner of matrix C: \n");
    for (i=0; i<min(M,6); i++) {
      for (j=0; j<min(N,6); j++) {
        printf ("%12.5G", C[j+i*N]);
      }
      printf ("\n");
    }
    #endif

    cout << scientific << setprecision(3) << "* RowMajor " << Order << "\t" ;
    cout << "transA " << TransA << "\t" ;
    cout << "transB " << TransB << "\t" ;
    cout << "m " << M << "\t" ;
    cout << "n " << N << "\t" ;
    cout << "k " << K << "\n" ;
    cout << "* lda " << lda << "\t" ;
    cout << "ldb " << ldb << "\t" ;
    cout << "ldc " << ldc << "\t" ;
    cout << "alpha " << alpha << "\t" ;
    cout << "beta " << beta << "\n" ;
    cout << "* A03 " << A[3] << "\t" ;
    cout << "B03 " << B[3] << "\t" ;
    cout << "C03 " << C[3] << "\t" ;
    cout << "Asum " << sum(A,M,K) << "\t" ;
    cout << "Bsum " << sum(B,K,N) << "\n" ;
    cout << "Csum " << sum(C,M,N) << "\n" ;
  }
#endif // NDEBUG

  // Check for (integer) overflows
  enforce(M>0);
  enforce(N>0);
  enforce(K>0);

  // check_int_mult_overflow(560000,8000); // fails on default int (32-bits)
  check_int_mult_overflow(M,K);
  check_int_mult_overflow(N,K);
  check_int_mult_overflow(M,N);

  cblas_dgemm(Order,TransA,TransB,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);

#ifndef NDEBUG
  #ifdef DISABLED
  if (is_debug_mode()) {
    printf (" Top left corner of matrix A (cols=k %i, rows=m %i): \n",K,M);
    for (i=0; i<min(M,6); i++) {
      for (j=0; j<min(K,6); j++) {
        printf ("%12.0f", A[j+i*K]);
      }
      printf ("\n");
    }

    printf ("\n Top left corner of matrix B: \n");
    for (i=0; i<min(K,6); i++) {
      for (j=0; j<min(N,6); j++) {
        printf ("%12.0f", B[j+i*N]);
      }
      printf ("\n");
    }

    printf ("\n Top left corner of matrix C: \n");
    for (i=0; i<min(M,6); i++) {
      for (j=0; j<min(N,6); j++) {
      printf ("%12.5G", C[j+i*N]);
      }
      printf ("\n");
    }
  }
  #endif
#endif // NDEBUG
}

/*
    Helper function fast_cblas_dgemm converts a GEMMA layout to cblas_dgemm.
*/
static void fast_cblas_dgemm(const char *TransA, const char *TransB, const double alpha,
                             const gsl_matrix *A, const gsl_matrix *B, const double beta,
                             gsl_matrix *C) {
  // C++ is row-major
  auto transA = (*TransA == 'N' || *TransA == 'n' ? CblasNoTrans : CblasTrans);
  auto transB = (*TransB == 'N' || *TransB == 'n' ? CblasNoTrans : CblasTrans);
  const size_t M   = C->size1;
  const size_t N   = C->size2;
  const size_t MA  = (transA == CblasNoTrans) ? A->size1 : A->size2;
  const size_t NA  = (transA == CblasNoTrans) ? A->size2 : A->size1;
  const size_t MBx = (transB == CblasNoTrans) ? B->size1 : B->size2;
  const size_t NB  = (transB == CblasNoTrans) ? B->size2 : B->size1;

  if (M == MA && N == NB && NA == MBx) {  /* [MxN] = [MAxNA][MBxNB] */

    auto K = NA;

    // Check for (integer) overflows
    enforce(M>0);
    enforce(N>0);
    enforce(K>0);

    // check_int_mult_overflow(560000,8000);
    check_int_mult_overflow(M,K);
    check_int_mult_overflow(N,K);
    check_int_mult_overflow(M,N);

    cblas_dgemm (CblasRowMajor, transA, transB, M, N, NA,
                 alpha, A->data, A->tda, B->data, B->tda, beta,
                 C->data, C->tda);

  } else {
    fail_msg("Range error in dgemm");
  }
}


/*
   Use the fast/supported way to call BLAS dgemm
*/

void fast_dgemm(const char *TransA, const char *TransB, const double alpha,
                const gsl_matrix *A, const gsl_matrix *B, const double beta,
                gsl_matrix *C) {
  fast_cblas_dgemm(TransA,TransB,alpha,A,B,beta,C);

#ifdef DISABLE
  if (is_check_mode()) {
    // ---- validate with original implementation
    gsl_matrix *C1 = gsl_matrix_alloc(C->size1,C->size2);
    eigenlib_dgemm(TransA,TransB,alpha,A,B,beta,C1);
    enforce_msg(gsl_matrix_equal(C,C1),"dgemm outcomes are not equal for fast & eigenlib");
    gsl_matrix_free(C1);
  }
#endif
}

void fast_eigen_dgemm(const char *TransA, const char *TransB, const double alpha,
                      const gsl_matrix *A, const gsl_matrix *B, const double beta,
                      gsl_matrix *C) {
    fast_cblas_dgemm(TransA,TransB,alpha,A,B,beta,C);
}

/*
 *  Inverse in place
 */

#include <gsl/gsl_permutation.h>
// #include <gsl/gsl_linalg.h>

extern "C" {
  int gsl_linalg_LU_invert(const gsl_matrix * LU, const gsl_permutation * p, gsl_matrix * inverse);
  int gsl_linalg_LU_decomp(gsl_matrix * A, gsl_permutation * p, int * signum);
}

void gsl_matrix_inv(gsl_matrix *m)
{
    size_t n=m->size1;

    gsl_matrix *temp1=gsl_matrix_calloc(n,n);
    gsl_matrix_memcpy(temp1,m);

    gsl_permutation *p=gsl_permutation_calloc(n);
    int sign=0;
    gsl_linalg_LU_decomp(temp1,p,&sign);
    gsl_matrix *inverse=gsl_matrix_calloc(n,n);

    gsl_linalg_LU_invert(temp1,p,inverse);
    gsl_matrix_memcpy(m,inverse);

    gsl_permutation_free(p);
    gsl_matrix_free(temp1);
    gsl_matrix_free(inverse);

}

void fast_inverse(gsl_matrix *m) {
    gsl_matrix_inv(m);
}
