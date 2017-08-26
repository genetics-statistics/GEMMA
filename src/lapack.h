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

#ifndef __LAPACK_H__
#define __LAPACK_H__

#include <vector>

using namespace std;

void lapack_cholesky_decomp(gsl_matrix *A);
void lapack_cholesky_solve(gsl_matrix *A, const gsl_vector *b, gsl_vector *x);
void lapack_dgemm(char *TransA, char *TransB, double alpha, const gsl_matrix *A,
                  const gsl_matrix *B, double beta, gsl_matrix *C);
void lapack_eigen_symmv(gsl_matrix *A, gsl_vector *eval, gsl_matrix *evec,
                        const size_t flag_largematrix);
double EigenDecomp(gsl_matrix *G, gsl_matrix *U, gsl_vector *eval,
                   const size_t flag_largematrix);
double EigenDecomp_Zeroed(gsl_matrix *G, gsl_matrix *U, gsl_vector *eval,
                   const size_t flag_largematrix);
double CholeskySolve(gsl_matrix *Omega, gsl_vector *Xty, gsl_vector *OiXty);
void LUDecomp(gsl_matrix *LU, gsl_permutation *p, int *signum);
void LUInvert(const gsl_matrix *LU, const gsl_permutation *p, gsl_matrix *ret_inverse);
double LULndet(const gsl_matrix *LU);
void LUSolve(const gsl_matrix *LU, const gsl_permutation *p,
             const gsl_vector *b, gsl_vector *x);
bool lapack_ddot(vector<double> &x, vector<double> &y, double &v);

#endif
