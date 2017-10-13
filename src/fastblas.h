#ifndef __FASTBLAS_H__
#define __FASTBLAS_H__

#include <assert.h>
#include <iostream>
#include "gsl/gsl_matrix.h"

gsl_matrix *fast_copy(gsl_matrix *m, const double *mem);

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

void fast_dgemm(const char *TransA, const char *TransB, const double alpha,
                const gsl_matrix *A, const gsl_matrix *B, const double beta,
                gsl_matrix *C);

#endif
