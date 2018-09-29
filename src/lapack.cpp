/*
    Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright (C) 2011-2017 Xiang Zhou

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

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_sys.h" // for gsl_isnan, gsl_isinf, gsl_isfinite
#include <cmath>
#include <iostream>
#include <vector>

#include "debug.h"
#include "lapack.h"
#include "mathfunc.h"

using namespace std;

extern "C" void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
                       double *ALPHA, double *A, int *LDA, double *B, int *LDB,
                       double *BETA, double *C, int *LDC);
extern "C" void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFO);
extern "C" void dpotrs_(char *UPLO, int *N, int *NRHS, double *A, int *LDA,
                        double *B, int *LDB, int *INFO);
extern "C" void dsyev_(char *JOBZ, char *UPLO, int *N, double *A, int *LDA,
                       double *W, double *WORK, int *LWORK, int *INFO);
extern "C" void dsyevr_(char *JOBZ, char *RANGE, char *UPLO, int *N, double *A,
                        int *LDA, double *VL, double *VU, int *IL, int *IU,
                        double *ABSTOL, int *M, double *W, double *Z, int *LDZ,
                        int *ISUPPZ, double *WORK, int *LWORK, int *IWORK,
                        int *LIWORK, int *INFO);
extern "C" double ddot_(int *N, double *DX, int *INCX, double *DY, int *INCY);

// Cholesky decomposition, A is destroyed.
void lapack_cholesky_decomp(gsl_matrix *A) {
  int N = A->size1, LDA = A->size1, INFO;
  char UPLO = 'L';

  if (N != (int)A->size2) {
    cout << "Matrix needs to be symmetric and same dimension in "
         << "lapack_cholesky_decomp." << endl;
    return;
  }

  dpotrf_(&UPLO, &N, A->data, &LDA, &INFO);
  if (INFO != 0) {
    cout << "Cholesky decomposition unsuccessful in "
         << "lapack_cholesky_decomp." << endl;
    return;
  }

  return;
}

// Cholesky solve, A is decomposed.
void lapack_cholesky_solve(gsl_matrix *A, const gsl_vector *b, gsl_vector *x) {
  int N = A->size1, NRHS = 1, LDA = A->size1, LDB = b->size, INFO;
  char UPLO = 'L';

  if (N != (int)A->size2 || N != LDB) {
    cout << "Matrix needs to be symmetric and same dimension in "
         << "lapack_cholesky_solve." << endl;
    return;
  }

  gsl_vector_memcpy(x, b);
  dpotrs_(&UPLO, &N, &NRHS, A->data, &LDA, x->data, &LDB, &INFO);
  if (INFO != 0) {
    cout << "Cholesky solve unsuccessful in lapack_cholesky_solve." << endl;
    return;
  }

  return;
}

void lapack_dgemm(char *TransA, char *TransB, double alpha, const gsl_matrix *A,
                  const gsl_matrix *B, double beta, gsl_matrix *C) {
  int M, N, K1, K2, LDA = A->size1, LDB = B->size1, LDC = C->size2;

  if (*TransA == 'N' || *TransA == 'n') {
    M = A->size1;
    K1 = A->size2;
  } else if (*TransA == 'T' || *TransA == 't') {
    M = A->size2;
    K1 = A->size1;
  } else {
    cout << "need 'N' or 'T' in lapack_dgemm" << endl;
    return;
  }

  if (*TransB == 'N' || *TransB == 'n') {
    N = B->size2;
    K2 = B->size1;
  } else if (*TransB == 'T' || *TransB == 't') {
    N = B->size1;
    K2 = B->size2;
  } else {
    cout << "need 'N' or 'T' in lapack_dgemm" << endl;
    return;
  }

  if (K1 != K2) {
    cout << "A and B not compatible in lapack_dgemm" << endl;
    return;
  }
  if (C->size1 != (size_t)M || C->size2 != (size_t)N) {
    cout << "C not compatible in lapack_dgemm" << endl;
    return;
  }

  gsl_matrix *A_t = gsl_matrix_alloc(A->size2, A->size1);
  gsl_matrix_transpose_memcpy(A_t, A);
  gsl_matrix *B_t = gsl_matrix_alloc(B->size2, B->size1);
  gsl_matrix_transpose_memcpy(B_t, B);
  gsl_matrix *C_t = gsl_matrix_alloc(C->size2, C->size1);
  gsl_matrix_transpose_memcpy(C_t, C);

  check_int_mult_overflow(M,K1);
  check_int_mult_overflow(N,K1);
  check_int_mult_overflow(M,N);

  dgemm_(TransA, TransB, &M, &N, &K1, &alpha, A_t->data, &LDA, B_t->data, &LDB,
         &beta, C_t->data, &LDC);

  gsl_matrix_transpose_memcpy(C, C_t);

  gsl_matrix_free(A_t);
  gsl_matrix_free(B_t);
  gsl_matrix_free(C_t);
  return;
}

// Eigenvalue decomposition, matrix A is destroyed. Returns eigenvalues in
// 'eval'. Also returns matrix 'evec' (U).
void lapack_eigen_symmv(gsl_matrix *A, gsl_vector *eval, gsl_matrix *evec,
                        const size_t flag_largematrix) {
  if (flag_largematrix == 1) { // not sure this flag is used!
    int N = A->size1, LDA = A->size1, INFO, LWORK = -1;
    char JOBZ = 'V', UPLO = 'L';

    if (N != (int)A->size2 || N != (int)eval->size) {
      cout << "Matrix needs to be symmetric and same "
           << "dimension in lapack_eigen_symmv." << endl;
      return;
    }

    LWORK = 3 * N;
    double *WORK = new double[LWORK];
    dsyev_(&JOBZ, &UPLO, &N, A->data, &LDA, eval->data, WORK, &LWORK, &INFO);
    if (INFO != 0) {
      cout << "Eigen decomposition unsuccessful in "
           << "lapack_eigen_symmv." << endl;
      return;
    }

    gsl_matrix_view A_sub = gsl_matrix_submatrix(A, 0, 0, N, N);
    gsl_matrix_memcpy(evec, &A_sub.matrix);
    gsl_matrix_transpose(evec);

    delete[] WORK;
  } else {
    // entering here
    int N = A->size1, LDA = A->size1, LDZ = A->size1, INFO;
    int LWORK = -1, LIWORK = -1;
    char JOBZ = 'V', UPLO = 'L', RANGE = 'A';
    double ABSTOL = 1.0E-7;

    // VL, VU, IL, IU are not referenced; M equals N if RANGE='A'.
    double VL = 0.0, VU = 0.0;
    int IL = 0, IU = 0, M;

    if (N != (int)A->size2 || N != (int)eval->size) {
      cout << "Matrix needs to be symmetric and same "
           << "dimension in lapack_eigen_symmv." << endl;
      return;
    }

    int *ISUPPZ = new int[2 * N];

    double WORK_temp[1];
    int IWORK_temp[1];

    // disable fast NaN checking for now - dsyevr throws NaN errors,
    // but fixes them (apparently)
    if (is_check_mode()) disable_segfpe();

    // DSYEVR - computes selected eigenvalues and, optionally,
    // eigenvectors of a real symmetric matrix
    // Here compute both (JOBZ is V), all eigenvalues (RANGE is A)
    // Lower triangle is stored (UPLO is L)
    dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A->data, &LDA, &VL, &VU, &IL, &IU,
            &ABSTOL, &M, eval->data, evec->data, &LDZ, ISUPPZ, WORK_temp,
            &LWORK, IWORK_temp, &LIWORK, &INFO);
    // If info = 0, the execution is successful.
    // If info = -i, the i-th parameter had an illegal value.
    // If info = i, an internal error has occurred.

    if (INFO != 0) cerr << "ERROR: value of INFO is " << INFO;
    enforce_msg(INFO == 0, "lapack_eigen_symmv failed");
    LWORK = (int)WORK_temp[0];    // The dimension of the array work.
    LIWORK = (int)IWORK_temp[0];  // The dimension of the array iwork, lwork≥ max(1, 10n).

    double *WORK = new double[LWORK];
    int *IWORK = new int[LIWORK];

    dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A->data, &LDA, &VL, &VU, &IL, &IU,
            &ABSTOL, &M, eval->data, evec->data, &LDZ, ISUPPZ, WORK, &LWORK,
            IWORK, &LIWORK, &INFO);
    if (INFO != 0) cerr << "ERROR: value of INFO is " << INFO;
    enforce_msg(INFO == 0, "lapack_eigen_symmv failed");

    if (is_check_mode()) enable_segfpe(); // reinstate fast NaN checking

    gsl_matrix_transpose(evec);

    delete[] ISUPPZ;
    delete[] WORK;
    delete[] IWORK;
  }

  return;
}

// Does NOT set eigenvalues to be positive. G gets destroyed. Returns
// eigen trace and values in U and eval (eigenvalues).
double EigenDecomp(gsl_matrix *G, gsl_matrix *U, gsl_vector *eval,
                   const size_t flag_largematrix) {
  lapack_eigen_symmv(G, eval, U, flag_largematrix);
  assert(!has_nan(eval));
  // write(eval,"eval");

  // Calculate track_G=mean(diag(G)).
  double d = 0.0;
  for (size_t i = 0; i < eval->size; ++i)
    d += gsl_vector_get(eval, i);

  d /= (double)eval->size;

  return d;
}

// Does NOT set eigenvalues to be positive. G gets destroyed. Returns
// eigen trace and values in U and eval (eigenvalues).  Same as
// EigenDecomp but zeroes eigenvalues close to zero. When negative
// eigenvalues remain a warning is issued.
double EigenDecomp_Zeroed(gsl_matrix *G, gsl_matrix *U, gsl_vector *eval,
                          const size_t flag_largematrix) {
  EigenDecomp(G,U,eval,flag_largematrix);
  auto d = 0.0;
  int count_zero_eigenvalues = 0;
  int count_negative_eigenvalues = 0;
  for (size_t i = 0; i < eval->size; i++) {
    // if (std::abs(gsl_vector_get(eval, i)) < EIGEN_MINVALUE)
    if (gsl_vector_get(eval, i) < 1e-10)
      gsl_vector_set(eval, i, 0.0);
    // checks
    if (gsl_vector_get(eval,i) == 0.0)
      count_zero_eigenvalues += 1;
    if (gsl_vector_get(eval,i) < -EIGEN_MINVALUE) // count smaller than -EIGEN_MINVALUE
      count_negative_eigenvalues += 1;
    d += gsl_vector_get(eval, i);
  }
  d /= (double)eval->size;
  if (count_zero_eigenvalues > 1) {
    write(eval,"eigenvalues");
    std::string msg = "Matrix G has ";
    msg += std::to_string(count_zero_eigenvalues);
    msg += " eigenvalues close to zero";
    warning_msg(msg);
  }
  const bool negative_eigen_values = has_negative_values_but_one(eval);
  if (negative_eigen_values) {
    write(eval,"eigenvalues");
    warning_msg("K has more than one negative eigenvalues!");
  }
  return d;
}

double CholeskySolve(gsl_matrix *Omega, gsl_vector *Xty, gsl_vector *OiXty) {
  double logdet_O = 0.0;

  lapack_cholesky_decomp(Omega);
  for (size_t i = 0; i < Omega->size1; ++i) {
    logdet_O += log(gsl_matrix_get(Omega, i, i));
  }
  logdet_O *= 2.0;
  lapack_cholesky_solve(Omega, Xty, OiXty);

  return logdet_O;
}

// LU decomposition.
void LUDecomp(gsl_matrix *LU, gsl_permutation *p, int *signum) {
  // debug_msg("entering");
  enforce_gsl(gsl_linalg_LU_decomp(LU, p, signum));
  return;
}

// LU invert. Returns inverse. Note that GSL does not recommend using
// this function

// These functions compute the inverse of a matrix A from its LU
// decomposition (LU,p), storing the result in the matrix inverse. The
// inverse is computed by solving the system A x = b for each column
// of the identity matrix. It is preferable to avoid direct use of the
// inverse whenever possible, as the linear solver functions can
// obtain the same result more efficiently and reliably (consult any
// introductory textbook on numerical linear algebra for details).
void LUInvert(const gsl_matrix *LU, const gsl_permutation *p, gsl_matrix *ret_inverse) {
  // debug_msg("entering");
  if (is_check_mode())
    LULndet(LU);

  enforce_gsl(gsl_linalg_LU_invert(LU, p, ret_inverse));
}

// LU lndet.

// These functions compute the logarithm of the absolute value of the
// determinant of a matrix A, \ln|\det(A)|, from its LU decomposition,
// LU. This function may be useful if the direct computation of the
// determinant would overflow or underflow.

double LULndet(const gsl_matrix *LU) {
  // debug_msg("entering");

  double res = gsl_linalg_LU_lndet((gsl_matrix *)LU);
  enforce_msg(!is_inf(res), "LU determinant is zero -> LU is not invertable");
  return res;
}

// LU solve.
void LUSolve(const gsl_matrix *LU, const gsl_permutation *p,
             const gsl_vector *b, gsl_vector *x) {
  // debug_msg("entering");
  enforce_gsl(gsl_linalg_LU_solve(LU, p, b, x));
  return;
}

bool lapack_ddot(vector<double> &x, vector<double> &y, double &v) {
  bool flag = false;
  int incx = 1;
  int incy = 1;
  int n = (int)x.size();
  if (x.size() == y.size()) {
    v = ddot_(&n, &x[0], &incx, &y[0], &incy);
    flag = true;
  }

  return flag;
}
