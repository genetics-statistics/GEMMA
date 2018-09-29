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

#include "Eigen/Dense"
// #include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
// #include "gsl/gsl_vector.h"
#include <cmath>
#include <iostream>
#include <vector>
// #include <cblas.h>

using namespace std;
using namespace Eigen;


// On two different clusters, compare eigen vs lapack/gsl:
//
// dgemm, 5x or 0.5x faster or slower than lapack, 5x or 4x faster than gsl
// dgemv, 20x or 4x faster than gsl,
// eigen, 1x or 0.3x slower than lapack
// invert, 20x or 10x faster than lapack
//
void eigenlib_dgemm(const char *TransA, const char *TransB, const double alpha,
                    const gsl_matrix *A, const gsl_matrix *B, const double beta,
                    gsl_matrix *C) {
  Map<Matrix<double, Dynamic, Dynamic, RowMajor>, 0, OuterStride<Dynamic>>
      A_mat(A->data, A->size1, A->size2, OuterStride<Dynamic>(A->tda));
  Map<Matrix<double, Dynamic, Dynamic, RowMajor>, 0, OuterStride<Dynamic>>
      B_mat(B->data, B->size1, B->size2, OuterStride<Dynamic>(B->tda));
  Map<Matrix<double, Dynamic, Dynamic, RowMajor>, 0, OuterStride<Dynamic>>
      C_mat(C->data, C->size1, C->size2, OuterStride<Dynamic>(C->tda));

  if (*TransA == 'N' || *TransA == 'n') {
    if (*TransB == 'N' || *TransB == 'n') {
      C_mat = alpha * A_mat * B_mat + beta * C_mat;
    } else {
      C_mat = alpha * A_mat * B_mat.transpose() + beta * C_mat;
    }
  } else {
    if (*TransB == 'N' || *TransB == 'n') {
      C_mat = alpha * A_mat.transpose() * B_mat + beta * C_mat;
    } else {
      C_mat = alpha * A_mat.transpose() * B_mat.transpose() + beta * C_mat;
    }
  }
}

void eigenlib_dgemv(const char *TransA, const double alpha, const gsl_matrix *A,
                    const gsl_vector *x, const double beta, gsl_vector *y) {
  Map<Matrix<double, Dynamic, Dynamic, RowMajor>, 0, OuterStride<Dynamic>>
      A_mat(A->data, A->size1, A->size2, OuterStride<Dynamic>(A->tda));
  Map<Matrix<double, Dynamic, 1>, 0, InnerStride<Dynamic>> x_vec(
      x->data, x->size, InnerStride<Dynamic>(x->stride));
  Map<Matrix<double, Dynamic, 1>, 0, InnerStride<Dynamic>> y_vec(
      y->data, y->size, InnerStride<Dynamic>(y->stride));

  if (*TransA == 'N' || *TransA == 'n') {
    y_vec = alpha * A_mat * x_vec + beta * y_vec;
  } else {
    y_vec = alpha * A_mat.transpose() * x_vec + beta * y_vec;
  }
}

void eigenlib_invert(gsl_matrix *A) {
  Map<Matrix<double, Dynamic, Dynamic, RowMajor>> A_mat(A->data, A->size1,
                                                        A->size2);
  A_mat = A_mat.inverse();
}

void eigenlib_dsyr(const double alpha, const gsl_vector *b, gsl_matrix *A) {
  Map<Matrix<double, Dynamic, Dynamic, RowMajor>> A_mat(A->data, A->size1,
                                                        A->size2);
  Map<Matrix<double, Dynamic, 1>, 0, OuterStride<Dynamic>> b_vec(
      b->data, b->size, OuterStride<Dynamic>(b->stride));
  A_mat = alpha * b_vec * b_vec.transpose() + A_mat;
}

void eigenlib_eigensymm(const gsl_matrix *G, gsl_matrix *U, gsl_vector *eval) {
  Map<Matrix<double, Dynamic, Dynamic, RowMajor>, 0, OuterStride<Dynamic>>
      G_mat(G->data, G->size1, G->size2, OuterStride<Dynamic>(G->tda));
  Map<Matrix<double, Dynamic, Dynamic, RowMajor>, 0, OuterStride<Dynamic>>
      U_mat(U->data, U->size1, U->size2, OuterStride<Dynamic>(U->tda));
  Map<Matrix<double, Dynamic, 1>, 0, OuterStride<Dynamic>> eval_vec(
      eval->data, eval->size, OuterStride<Dynamic>(eval->stride));

  SelfAdjointEigenSolver<MatrixXd> es(G_mat);
  if (es.info() != Success)
    abort();
  eval_vec = es.eigenvalues();
  U_mat = es.eigenvectors();
}
