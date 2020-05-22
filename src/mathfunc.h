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

#ifndef __MATHFUNC_H__
#define __MATHFUNC_H__

// #include "Eigen/Dense"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"

#define CONDITIONED_MAXRATIO 2e+6 // based on http://mathworld.wolfram.com/ConditionNumber.html
#define EIGEN_MINVALUE 1e-10

using namespace std;

inline bool is_nan(double f) {
  return (isnan(f));
}

#define is_inf(d) gsl_isinf(d)
#define is_nan(d) gsl_isnan(d)

bool has_nan(const vector<double> v);
bool has_nan(const gsl_vector *v);
bool has_inf(const gsl_vector *v);
bool has_nan(const gsl_matrix *m);
bool has_inf(const gsl_matrix *m);

bool is_integer(const std::string & s);

double safe_log(const double d);
double safe_sqrt(const double d);

double VectorVar(const gsl_vector *v);
void CenterMatrix(gsl_matrix *G);
void CenterMatrix(gsl_matrix *G, const gsl_vector *w);
void CenterMatrix(gsl_matrix *G, const gsl_matrix *W);
void StandardizeMatrix(gsl_matrix *G);
double ScaleMatrix(gsl_matrix *G);
bool has_negative_values_but_one(const gsl_vector *v);
uint count_abs_small_values(const gsl_vector *v, double min);
bool isMatrixPositiveDefinite(const gsl_matrix *G);
bool isMatrixIllConditioned(const gsl_vector *eigenvalues, double max_ratio=CONDITIONED_MAXRATIO);
bool isMatrixSymmetric(const gsl_matrix *G);
gsl_vector *getEigenValues(const gsl_matrix *G);
double sum(const double *m, size_t rows, size_t cols);
double SumVector(const gsl_vector *v);
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
// void uchar_matrix_get_row(const vector<vector<unsigned char>> &X,
//                          const size_t i_row, Eigen::VectorXd &x_row);

#endif
