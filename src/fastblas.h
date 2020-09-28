/*
    Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright © 2011-2017, Xiang Zhou
    Copyright © 2017, Peter Carbonetto
    Copyright © 2017, Pjotr Prins

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

#ifndef __FASTBLAS_H__
#define __FASTBLAS_H__

#include <assert.h>
#include <iostream>
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_matrix.h"

gsl_matrix *fast_copy(gsl_matrix *m, const double *mem);

extern const char *FastblasTrans;
extern const char *FastblasNoTrans;

void fast_dgemm(const char *TransA, const char *TransB, const double alpha,
                const gsl_matrix *A, const gsl_matrix *B, const double beta,
                gsl_matrix *C);
void fast_eigen_dgemm(const char *TransA, const char *TransB, const double alpha,
                      const gsl_matrix *A, const gsl_matrix *B, const double beta,
                      gsl_matrix *C);

void fast_inverse(gsl_matrix *m);

#endif
