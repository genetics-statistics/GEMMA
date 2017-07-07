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

#include <iostream>
#include <fstream>
#include <sstream>

#include <iomanip>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cstring>
#include <algorithm>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_roots.h"
#include "Eigen/Dense"

#include "lapack.h"
#include "param.h"
#include "ldr.h"
#include "lm.h"
#include "mathfunc.h"

using namespace std;
using namespace Eigen;

void LDR::CopyFromParam (PARAM &cPar) {
	a_mode=cPar.a_mode;
	d_pace=cPar.d_pace;

	file_bfile=cPar.file_bfile;
	file_geno=cPar.file_geno;
	file_out=cPar.file_out;
	path_out=cPar.path_out;

	ni_total=cPar.ni_total;
	ns_total=cPar.ns_total;
	ni_test=cPar.ni_test;
	ns_test=cPar.ns_test;
	n_cvt=cPar.n_cvt;

	indicator_idv=cPar.indicator_idv;
	indicator_snp=cPar.indicator_snp;
	snpInfo=cPar.snpInfo;

	return;
}

void LDR::CopyToParam (PARAM &cPar) {
	return;
}

//X is a p by n matrix.
void LDR::VB (const vector<vector<unsigned char> > &Xt,
	      const gsl_matrix *W_gsl, const gsl_vector *y_gsl) {

  // Save gsl_vector and gsl_matrix into Eigen library formats.
  MatrixXd W(W_gsl->size1, W_gsl->size2);
  VectorXd y(y_gsl->size);
  VectorXd x_col(y_gsl->size);

  double d;
  for (size_t i=0; i<W_gsl->size1; i++) {
    d=gsl_vector_get(y_gsl, i);
    y(i)=d;
    for (size_t j=0; j<W_gsl->size2; j++) {
      W(i,j)=gsl_matrix_get(W_gsl, i, j);
    }
  }

  // Initial VB values by lm.
  cout<<indicator_snp[0]<<" "<<indicator_snp[1]<<" "<<indicator_snp[2]<<endl;
  uchar_matrix_get_row (Xt, 0, x_col);

  for (size_t j=0; j<10; j++) {
    cout<<x_col(j)<<endl;
  }

  // Run VB iterations.
  // TO DO.

  // Save results.
  // TO DO.

  return;
}
