/*
    Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright (C) 2011-2017, Xiang Zhou

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

#include "gsl/gsl_matrix.h"
#include <cmath>
#include <iomanip>
#include <vector>
#include <cblas.h>
#include "debug.h"
#include "fastblas.h"
#include "mathfunc.h"
#ifndef NDEBUG
#include "eigenlib.h"
#endif

using namespace std;

/*
   Reasonably fast function to copy data from standard C array into
   gsl_matrix. Avoid it for performance critical sections.
*/
gsl_matrix *fast_copy(gsl_matrix *m, const double *mem) {
  auto rows = m->size1;
  auto cols = m->size2;
  if (is_strict_mode()) { // slower correct version
    for (auto r=0; r<rows; r++) {
      for (auto c=0; c<cols; c++) {
        gsl_matrix_set(m,r,c,mem[r*cols+c]);
      }
    }
  } else { // faster goes by row
    for (auto r=0; r<rows; r++) {
      auto v = gsl_vector_calloc(cols);
      assert(v->size == cols);
      assert(v->block->size == cols);
      assert(v->stride == 1);
      memcpy(v->block->data,&mem[r*cols],cols*sizeof(double));
      gsl_matrix_set_row(m,r,v);
    }
  }
  return m;
}

/*
    Helper function fast_cblas_dgemm runs the local dgemm
*/
void fast_cblas_dgemm(OPENBLAS_CONST enum CBLAS_ORDER Order,
                      OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                      OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB,
                      OPENBLAS_CONST blasint M,
                      OPENBLAS_CONST blasint N,
                      OPENBLAS_CONST blasint K,
                      OPENBLAS_CONST double alpha,
                      OPENBLAS_CONST double *A,
                      OPENBLAS_CONST blasint lda,
                      OPENBLAS_CONST double *B,
                      OPENBLAS_CONST blasint ldb,
                      OPENBLAS_CONST double beta,
                      double *C,
                      OPENBLAS_CONST blasint ldc) {
#ifndef NDEBUG
  size_t i,j;
  if (is_debug_mode()) {
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

  cblas_dgemm(Order,TransA,TransB,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);

#ifndef NDEBUG
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
  // A(m x k) * B(k x n) = C(m x n))
  auto rowsA = A->size1;
  auto colsA = A->size2;
  blasint M = A->size1;
  blasint K = B->size1;
  assert(K == colsA);
  blasint N = B->size2;
  // cout << M << "," << N "," << K << endl;
  // Layout = CblasRowMajor: Trans: K , NoTrans M
  blasint lda = (transA==CblasNoTrans ? K : M );
  blasint ldb = (transB==CblasNoTrans ? N : K );
  blasint ldc = N;

  fast_cblas_dgemm(CblasRowMajor, transA, transB, M, N, K, alpha,
              /* A */ A->data,
              /* lda */ lda,
              /* B */ B->data,
              /* ldb */ ldb,
              /* beta */ beta,
              /* C */ C->data, ldc);
}

/*
   Use the fasted/supported way to call BLAS dgemm
*/

void fast_dgemm(const char *TransA, const char *TransB, const double alpha,
                const gsl_matrix *A, const gsl_matrix *B, const double beta,
                gsl_matrix *C) {
  fast_cblas_dgemm(TransA,TransB,alpha,A,B,beta,C);

  #ifndef NDEBUG
  if (is_strict_mode() && !is_no_check_mode()) {
    // ---- validate with original implementation
    gsl_matrix *C1 = gsl_matrix_alloc(C->size1,C->size2);
    eigenlib_dgemm(TransA,TransB,alpha,A,B,beta,C1);
    enforce_msg(gsl_matrix_equal(C,C1),"dgemm outcomes are not equal for fast & eigenlib");
  }
  #endif
}
