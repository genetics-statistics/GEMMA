#ifndef __FASTBLAS_H__
#define __FASTBLAS_H__

#include <assert.h>
#include <iostream>
#include "gsl/gsl_matrix.h"

gsl_matrix *fast_copy(gsl_matrix *m, const double *mem);

void fast_dgemm(const char *TransA, const char *TransB, const double alpha,
                const gsl_matrix *A, const gsl_matrix *B, const double beta,
                gsl_matrix *C);

#endif
