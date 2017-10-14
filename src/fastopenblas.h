#ifndef __FASTOPENBLAS_H__
#define __FASTOPENBLAS_H__

#include <assert.h>
#include <iostream>
#include <cblas.h>   // For OpenBlas
#include "gsl/gsl_matrix.h"

typedef long blasint;

void fast_cblas_dgemm(const enum CBLAS_ORDER Order,
                      const enum CBLAS_TRANSPOSE TransA,
                      const enum CBLAS_TRANSPOSE TransB,
                      const blasint M,
                      const blasint N,
                      const blasint K,
                      const double alpha,
                      const double *A,
                      const blasint lda,
                      const double *B,
                      const blasint ldb,
                      const double beta,
                      double *C,
                      const blasint ldc);
