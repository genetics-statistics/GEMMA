/*
    Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright © 2011-2017, Xiang Zhou
    Copyright © 2017, Peter Carbonetto
    Copyright © 2017-2018, Pjotr Prins

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

#include <bitset>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <vector>

// #include "Eigen/Dense"

#include "gsl/gsl_version.h"

#if GSL_MAJOR_VERSION < 2
#pragma message "GSL version " GSL_VERSION
#endif

#include "gsl/gsl_sys.h" // for gsl_isnan, gsl_isinf, gsl_isfinite
#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_eigen.h"

#include "debug.h"
// #include "eigenlib.h"
#include "fastblas.h"
#include "lapack.h"
#include "mathfunc.h"

using namespace std;
// using namespace Eigen;

bool has_nan(const vector<double> v) {
  if (!is_check_mode()) return false;

  for (const auto& e: v) {
    if (is_nan(e))
      return true;
  }
  return false;
}

bool has_nan(const gsl_vector *v) {
  if (!is_check_mode()) return false;

  for (size_t i = 0; i < v->size; ++i)
    if (is_nan(gsl_vector_get(v,i))) return true;
  return false;
}
bool has_inf(const gsl_vector *v) {
  if (!is_check_mode()) return false;

  for (size_t i = 0; i < v->size; ++i) {
    auto value = gsl_vector_get(v,i);
    if (is_inf(value) != 0) return true;
  }
  return false;
}
bool has_nan(const gsl_matrix *m) {
  if (!is_check_mode()) return false;

  for (size_t i = 0; i < m->size1; ++i)
    for (size_t j = 0; j < m->size2; ++j)
      if (is_nan(gsl_matrix_get(m,i,j))) return true;
  return false;
}
bool has_inf(const gsl_matrix *m) {
  if (!is_check_mode()) return false;

  for (size_t i = 0; i < m->size1; ++i)
    for (size_t j = 0; j < m->size2; ++j) {
      auto value = gsl_matrix_get(m,i,j);
      if (is_inf(value) != 0) return true;
    }
  return false;
}

bool is_integer(const std::string & s){
    return std::regex_match(s, std::regex("^[0-9]+$"));
}

bool is_float(const std::string & s){
    return std::regex_match(s, std::regex("^[+-]?([0-9]*[.])?[0-9]+$"));
}

double safe_log(const double d) {
  if (!is_legacy_mode() && (is_check_mode() || is_debug_mode()))
    enforce_msg(d > 0.0, (std::string("Trying to take the log of ") + std::to_string(d)).c_str());
  return log(d);
}

double safe_sqrt(const double d) {
  double d1 = d;
  if (fabs(d < 0.001))
    d1 = fabs(d);
  if (!is_legacy_mode() && (is_check_mode() || is_debug_mode()))
    enforce_msg(d1 >= 0.0, (std::string("Trying to take the sqrt of ") + std::to_string(d)).c_str());
  if (d1 < 0.0 )
    return nan("");
  return sqrt(d1);
}

// calculate variance of a vector
double VectorVar(const gsl_vector *v) {
  double d, m = 0.0, m2 = 0.0;
  for (size_t i = 0; i < v->size; ++i) {
    d = gsl_vector_get(v, i);
    m += d;
    m2 += d * d;
  }
  m /= (double)v->size;
  m2 /= (double)v->size;
  return m2 - m * m;
}

// Center the matrix G.
void CenterMatrix(gsl_matrix *G) {
  double d;
  gsl_vector *w = gsl_vector_safe_alloc(G->size1);
  gsl_vector *Gw = gsl_vector_safe_alloc(G->size1);
  gsl_vector_set_all(w, 1.0);

  //  y := alpha*A*x+ beta*y or Gw = G*w
  gsl_blas_dgemv(CblasNoTrans, 1.0, G, w, 0.0, Gw);

  // int gsl_blas_dsyr2(CBLAS_UPLO_t Uplo, double alpha, const gsl_vector * x, const gsl_vector * y, gsl_matrix * A)
  // compute the symmetric rank-2 update A = \alpha x y^T + \alpha y x^T + A of the symmetric matrix A. Since the matrix A is symmetric only its upper half or lower half need to be stored
  // or G = (UpperTriangle) alpha*Gw*w' + alpha*w*Gw' + G
  gsl_blas_dsyr2(CblasUpper, -1.0 / (double)G->size1, Gw, w, G);
  // compute dot product of vectors w.Gw and store in d
  gsl_blas_ddot(w, Gw, &d);
  // G = (upper) alpha w*w' + G
  gsl_blas_dsyr(CblasUpper, d / ((double)G->size1 * (double)G->size1), w, G);

  // Transpose the matrix
  for (size_t i = 0; i < G->size1; ++i) {
    for (size_t j = 0; j < i; ++j) {
      d = gsl_matrix_get(G, j, i);
      gsl_matrix_set(G, i, j, d);
    }
  }

  gsl_vector_safe_free(w);
  gsl_vector_safe_free(Gw);

  return;
}

// Center the matrix G.
// Only used in vc
void CenterMatrix(gsl_matrix *G, const gsl_vector *w) {
  double d, wtw;
  gsl_vector *Gw = gsl_vector_safe_alloc(G->size1);

  gsl_blas_ddot(w, w, &wtw);
  gsl_blas_dgemv(CblasNoTrans, 1.0, G, w, 0.0, Gw);
  gsl_blas_dsyr2(CblasUpper, -1.0 / wtw, Gw, w, G);
  gsl_blas_ddot(w, Gw, &d);
  gsl_blas_dsyr(CblasUpper, d / (wtw * wtw), w, G);

  for (size_t i = 0; i < G->size1; ++i) {
    for (size_t j = 0; j < i; ++j) {
      d = gsl_matrix_get(G, j, i);
      gsl_matrix_set(G, i, j, d);
    }
  }

  gsl_vector_safe_free(Gw);

  return;
}

// Center the matrix G.
// Only used in vc
void CenterMatrix(gsl_matrix *G, const gsl_matrix *W) {
  gsl_matrix *WtW = gsl_matrix_safe_alloc(W->size2, W->size2);
  gsl_matrix *WtWi = gsl_matrix_safe_alloc(W->size2, W->size2);
  gsl_matrix *WtWiWt = gsl_matrix_safe_alloc(W->size2, G->size1);
  gsl_matrix *GW = gsl_matrix_safe_alloc(G->size1, W->size2);
  gsl_matrix *WtGW = gsl_matrix_safe_alloc(W->size2, W->size2);
  gsl_matrix *Gtmp = gsl_matrix_safe_alloc(G->size1, G->size1);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);

  int sig;
  gsl_permutation *pmt = gsl_permutation_alloc(W->size2);
  LUDecomp(WtW, pmt, &sig);
  LUInvert(WtW, pmt, WtWi);

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, WtWi, W, 0.0, WtWiWt);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, G, W, 0.0, GW);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, GW, WtWiWt, 0.0, Gtmp);

  gsl_matrix_sub(G, Gtmp);
  gsl_matrix_transpose(Gtmp);
  gsl_matrix_sub(G, Gtmp);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, GW, 0.0, WtGW);
  // GW is destroyed.
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, WtWiWt, WtGW, 0.0, GW);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, GW, WtWiWt, 0.0, Gtmp);

  gsl_matrix_add(G, Gtmp);

  gsl_matrix_safe_free(WtW);
  gsl_matrix_safe_free(WtWi);
  gsl_matrix_safe_free(WtWiWt);
  gsl_matrix_safe_free(GW);
  gsl_matrix_safe_free(WtGW);
  gsl_matrix_safe_free(Gtmp);

  return;
}

// "Standardize" the matrix G such that all diagonal elements = 1.
// (only used by vc)
void StandardizeMatrix(gsl_matrix *G) {
  double d = 0.0;
  vector<double> vec_d;

  for (size_t i = 0; i < G->size1; ++i) {
    vec_d.push_back(gsl_matrix_get(G, i, i));
  }
  for (size_t i = 0; i < G->size1; ++i) {
    for (size_t j = i; j < G->size2; ++j) {
      if (j == i) {
        gsl_matrix_set(G, i, j, 1);
      } else {
        d = gsl_matrix_get(G, i, j);
        d /= sqrt(vec_d[i] * vec_d[j]);
        gsl_matrix_set(G, i, j, d);
        gsl_matrix_set(G, j, i, d);
      }
    }
  }

  return;
}

// Scale the matrix G such that the mean diagonal = 1.
double ScaleMatrix(gsl_matrix *G) {
  double d = 0.0;

  // Compute mean of diagonal
  for (size_t i = 0; i < G->size1; ++i) {
    d += gsl_matrix_get(G, i, i);
  }
  d /= (double)G->size1;

  // Scale the matrix using the diagonal mean
  if (d != 0) {
    gsl_matrix_scale(G, 1.0 / d);
  }

  return d;
}

bool isMatrixSymmetric(const gsl_matrix *G) {
  enforce(G->size1 == G->size2);
  auto m = G->data;
  // upper triangle
  auto size = G->size1;
  for(size_t c = 0; c < size; c++) {
    for(size_t r = 0; r < c; r++) {
      int x1 = c, y1 = r, x2 = r, y2 = c;
      auto idx1 = y1*size+x1, idx2 = y2*size+x2;
      // printf("(%d,%d %f - %d,%d %f)",x1,y1,m[idx1],x2,y2,m[idx2]);
      if(m[idx1] != m[idx2]) {
        cout << "Mismatch coordinates (" << c << "," << r << ")" << m[idx1] << ":" << m[idx2] << "!" << endl;
        return false;
      }
    }
  }
  return true;
}

bool isMatrixPositiveDefinite(const gsl_matrix *G) {
  enforce(G->size1 == G->size2);
  auto G2 = gsl_matrix_safe_alloc(G->size1, G->size2);
  enforce_gsl(gsl_matrix_safe_memcpy(G2,G));
  auto handler = gsl_set_error_handler_off();
#if GSL_MAJOR_VERSION >= 2 && GSL_MINOR_VERSION >= 3
  auto s = gsl_linalg_cholesky_decomp1(G2);
#else
  auto s = gsl_linalg_cholesky_decomp(G2);
#endif
  gsl_set_error_handler(handler);
  if (s == GSL_SUCCESS) {
    gsl_matrix_safe_free(G2);
    return true;
  }
  gsl_matrix_free(G2);
  return (false);
}

gsl_vector *getEigenValues(const gsl_matrix *G) {
  enforce(G->size1 == G->size2);
  auto G2 = gsl_matrix_safe_alloc(G->size1, G->size2);
  enforce_gsl(gsl_matrix_safe_memcpy(G2,G));
  auto eworkspace = gsl_eigen_symm_alloc(G->size1);
  enforce(eworkspace);
  gsl_vector *eigenvalues = gsl_vector_safe_alloc(G->size1);
  enforce_gsl(gsl_eigen_symm(G2, eigenvalues, eworkspace));
  gsl_eigen_symm_free(eworkspace);
  gsl_matrix_safe_free(G2);
  return eigenvalues;
}

// Check whether eigen values are larger than *min*
// by default 1E-5.
// Returns success, eigen min, eigen min-but-1, eigen max

tuple<double, double, double> minmax(const gsl_vector *v) {
  auto min  = v->data[0];
  auto min1 = min;
  auto max  = min;
  for (size_t i=1; i<v->size; i++) {
    auto value = std::abs(v->data[i]);
    if (value < min) {
      min1 = min;
      min = value;
    }
    if (value > max)
      max = value;
  }
  return std::make_tuple(min, min1, max);
}

tuple<double, double, double> abs_minmax(const gsl_vector *v) {
  auto min  = std::abs(v->data[0]);
  auto min1 = min;
  auto max  = min;
  for (size_t i=1; i<v->size; i++) {
    auto value = std::abs(v->data[i]);
    if (value < min) {
      min1 = min;
      min = value;
    }
    if (value > max)
      max = value;
  }
  return std::make_tuple(min, min1, max);
}

// Check for negative values. skip_min will leave out
// the lowest value
bool has_negative_values_but_one(const gsl_vector *v) {
  bool one_skipped = false;
  for (size_t i=0; i<v->size; i++) {
    if (v->data[i] < -EIGEN_MINVALUE) {
      if (one_skipped)
        return true;
      one_skipped = true;
    }
  }
  return false;
}

uint count_abs_small_values(const gsl_vector *v, double min) {
  uint count = 0;
  for (size_t i=0; i<v->size; i++) {
    if (std::abs(v->data[i]) < min) {
      count += 1;
    }
  }
  return count;
}

// Check for matrix being ill conditioned by taking the eigen values
// and the ratio of max and min but one (min is expected to be zero).
bool isMatrixIllConditioned(const gsl_vector *eigenvalues, double max_ratio) {
  auto t = abs_minmax(eigenvalues);

#if !defined NDEBUG
  auto absmin = get<0>(t);
#endif

  auto absmin1 = get<1>(t);
  auto absmax = get<2>(t);
  if (absmax/absmin1 > max_ratio) {
    #if !NDEBUG
    cerr << "**** DEBUG: Ratio |eigenmax|/|eigenmin| suggests matrix is ill conditioned for double precision" << endl;
    auto t = minmax(eigenvalues);
    auto min = get<0>(t);
    auto min1 = get<1>(t);
    auto max = get<2>(t);
    cerr << "**** DEBUG: Abs eigenvalues [Min " << absmin << ", " << absmin1 << " ... " << absmax << " Max] Ratio (" << absmax << "/" << absmin1 << ") = " << absmax/absmin1 << endl;
    cerr << "**** DEBUG: Eigenvalues [Min " << min << ", " << min1 << " ... " << max << " Max]" << endl;
    #endif
    return true;
  }
  return false;
}

double sum(const double *m, size_t rows, size_t cols) {
  double sum = 0.0;
  for (size_t i = 0; i<rows*cols; i++)
    sum += m[i];
  return sum;
}

double SumVector(const gsl_vector *v) {
  double sum = 0;
  for (size_t i = 0; i < v->size; i++ ) {
    sum += gsl_vector_get(v, i);
  }
  return( sum );
}

// Center the vector y.
double CenterVector(gsl_vector *y) {
  double d = 0.0;

  for (size_t i = 0; i < y->size; ++i) {
    d += gsl_vector_get(y, i);
  }
  d /= (double)y->size;

  gsl_vector_add_constant(y, -1.0 * d);

  return d;
}

// Center the vector y.
void CenterVector(gsl_vector *y, const gsl_matrix *W) {
  gsl_matrix *WtW = gsl_matrix_safe_alloc(W->size2, W->size2);
  gsl_vector *Wty = gsl_vector_safe_alloc(W->size2);
  gsl_vector *WtWiWty = gsl_vector_safe_alloc(W->size2);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
  gsl_blas_dgemv(CblasTrans, 1.0, W, y, 0.0, Wty);

  int sig;
  gsl_permutation *pmt = gsl_permutation_alloc(W->size2);
  LUDecomp(WtW, pmt, &sig);
  LUSolve(WtW, pmt, Wty, WtWiWty);

  gsl_blas_dgemv(CblasNoTrans, -1.0, W, WtWiWty, 1.0, y);

  gsl_matrix_safe_free(WtW);
  gsl_vector_safe_free(Wty);
  gsl_vector_safe_free(WtWiWty);

  return;
}

// "Standardize" vector y to have mean 0 and y^ty/n=1.
void StandardizeVector(gsl_vector *y) {
  double d = 0.0, m = 0.0, v = 0.0;

  for (size_t i = 0; i < y->size; ++i) {
    d = gsl_vector_get(y, i);
    m += d;
    v += d * d;
  }
  m /= (double)y->size;
  v /= (double)y->size;
  v -= m * m;

  gsl_vector_add_constant(y, -1.0 * m);
  gsl_vector_scale(y, 1.0 / sqrt(v));

  return;
}

// Calculate UtX (U gets transposed)
void CalcUtX(const gsl_matrix *U, gsl_matrix *UtX) {
  gsl_matrix *X = gsl_matrix_safe_alloc(UtX->size1, UtX->size2);
  gsl_matrix_safe_memcpy(X, UtX);
  fast_dgemm("T", "N", 1.0, U, X, 0.0, UtX);
  gsl_matrix_safe_free(X);
}

void CalcUtX(const gsl_matrix *U, const gsl_matrix *X, gsl_matrix *UtX) {
  fast_dgemm("T", "N", 1.0, U, X, 0.0, UtX);
}

void CalcUtX(const gsl_matrix *U, const gsl_vector *x, gsl_vector *Utx) {
  gsl_blas_dgemv(CblasTrans, 1.0, U, x, 0.0, Utx);
}

// Kronecker product.
void Kronecker(const gsl_matrix *K, const gsl_matrix *V, gsl_matrix *H) {
  for (size_t i = 0; i < K->size1; i++) {
    for (size_t j = 0; j < K->size2; j++) {
      gsl_matrix_view H_sub = gsl_matrix_submatrix(
          H, i * V->size1, j * V->size2, V->size1, V->size2);
      gsl_matrix_safe_memcpy(&H_sub.matrix, V);
      gsl_matrix_scale(&H_sub.matrix, gsl_matrix_get(K, i, j));
    }
  }
  return;
}

// Symmetric K matrix.
void KroneckerSym(const gsl_matrix *K, const gsl_matrix *V, gsl_matrix *H) {
  for (size_t i = 0; i < K->size1; i++) {
    for (size_t j = i; j < K->size2; j++) {
      gsl_matrix_view H_sub = gsl_matrix_submatrix(
          H, i * V->size1, j * V->size2, V->size1, V->size2);
      gsl_matrix_safe_memcpy(&H_sub.matrix, V);
      gsl_matrix_scale(&H_sub.matrix, gsl_matrix_get(K, i, j));

      if (i != j) {
        gsl_matrix_view H_sub_sym = gsl_matrix_submatrix(
            H, j * V->size1, i * V->size2, V->size1, V->size2);
        gsl_matrix_safe_memcpy(&H_sub_sym.matrix, &H_sub.matrix);
      }
    }
  }
  return;
}

// This function calculates HWE p value with methods described in
// Wigginton et al. (2005) AJHG; it is based on the code in plink 1.07.
double CalcHWE(const size_t n_hom1, const size_t n_hom2, const size_t n_ab) {
  if ((n_hom1 + n_hom2 + n_ab) == 0) {
    return 1;
  }

  // "AA" is the rare allele.
  int n_aa = n_hom1 < n_hom2 ? n_hom1 : n_hom2;
  int n_bb = n_hom1 < n_hom2 ? n_hom2 : n_hom1;

  int rare_copies = 2 * n_aa + n_ab;
  int genotypes = n_ab + n_bb + n_aa;

  double *het_probs = (double *)malloc((rare_copies + 1) * sizeof(double));
  if (het_probs == NULL)
    cout << "Internal error: SNP-HWE: Unable to allocate array" << endl;

  int i;
  for (i = 0; i <= rare_copies; i++)
    het_probs[i] = 0.0;

  // Start at midpoint.
  // XZ modified to add (long int)
  int mid = ((long int)rare_copies *
             (2 * (long int)genotypes - (long int)rare_copies)) /
            (2 * (long int)genotypes);

  // Check to ensure that midpoint and rare alleles have same
  // parity.
  if ((rare_copies & 1) ^ (mid & 1))
    mid++;

  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;

  het_probs[mid] = 1.0;
  double sum = het_probs[mid];
  for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
    het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets *
                               (curr_hets - 1.0) /
                               (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
    sum += het_probs[curr_hets - 2];

    // Two fewer heterozygotes for next iteration; add one
    // rare, one common homozygote.
    curr_homr++;
    curr_homc++;
  }

  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;
  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
    het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr *
                               curr_homc /
                               ((curr_hets + 2.0) * (curr_hets + 1.0));
    sum += het_probs[curr_hets + 2];

    // Add 2 heterozygotes for next iteration; subtract
    // one rare, one common homozygote.
    curr_homr--;
    curr_homc--;
  }

  for (i = 0; i <= rare_copies; i++)
    het_probs[i] /= sum;

  double p_hwe = 0.0;

  // p-value calculation for p_hwe.
  for (i = 0; i <= rare_copies; i++) {
    if (het_probs[i] > het_probs[n_ab])
      continue;
    p_hwe += het_probs[i];
  }

  p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

  free(het_probs);

  return p_hwe;
}

double UcharToDouble02(const unsigned char c) { return (double)c * 0.01; }

unsigned char Double02ToUchar(const double dosage) {
  return (int)(dosage * 100);
}
