#ifndef __FASTOPENBLAS_H__
#define __FASTOPENBLAS_H__

#include <assert.h>
#include <iostream>
#include <cblas.h>   // For OpenBlas
#include "gsl/gsl_matrix.h"

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
                      const size_t ldc);
