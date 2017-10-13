#ifndef __FASTOPENBLAS_H__
#define __FASTOPENBLAS_H__

#include <assert.h>
#include <iostream>
#include <cblas.h>   // For OpenBlas
#include "gsl/gsl_matrix.h"

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
                      OPENBLAS_CONST blasint ldc);

#endif
