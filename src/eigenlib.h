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

#ifndef __EIGENLIB_H__
#define __EIGENLIB_H__

// #include <vector>

// using namespace std;

void eigenlib_dgemm(const char *TransA, const char *TransB, const double alpha,
                    const gsl_matrix *A, const gsl_matrix *B, const double beta,
                    gsl_matrix *C);
void eigenlib_dgemv(const char *TransA, const double alpha, const gsl_matrix *A,
                    const gsl_vector *x, const double beta, gsl_vector *y);
void eigenlib_invert(gsl_matrix *A);
void eigenlib_dsyr(const double alpha, const gsl_vector *b, gsl_matrix *A);
void eigenlib_eigensymm(const gsl_matrix *G, gsl_matrix *U, gsl_vector *eval);

#endif
