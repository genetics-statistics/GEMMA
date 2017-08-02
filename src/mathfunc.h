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
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __MATHFUNC_H__
#define __MATHFUNC_H__

#include "Eigen/Dense"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"

using namespace std;
using namespace Eigen;

double VectorVar(const gsl_vector *v);
void CenterMatrix(gsl_matrix *G);
void CenterMatrix(gsl_matrix *G, const gsl_vector *w);
void CenterMatrix(gsl_matrix *G, const gsl_matrix *W);
void StandardizeMatrix(gsl_matrix *G);
double ScaleMatrix(gsl_matrix *G);
double CenterVector(gsl_vector *y);
void CenterVector(gsl_vector *y, const gsl_matrix *W);
void StandardizeVector(gsl_vector *y);
void CalcUtX(const gsl_matrix *U, gsl_matrix *UtX);
void CalcUtX(const gsl_matrix *U, const gsl_matrix *X, gsl_matrix *UtX);
void CalcUtX(const gsl_matrix *U, const gsl_vector *x, gsl_vector *Utx);
double CalcHWE(const size_t n_hom1, const size_t n_hom2, const size_t n_ab);
void Kronecker(const gsl_matrix *K, const gsl_matrix *V, gsl_matrix *H);
void KroneckerSym(const gsl_matrix *K, const gsl_matrix *V, gsl_matrix *H);

double UcharToDouble02(const unsigned char c);
unsigned char Double02ToUchar(const double dosage);
void uchar_matrix_get_row(const vector<vector<unsigned char>> &X,
                          const size_t i_row, VectorXd &x_row);

#endif
