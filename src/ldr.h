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

#ifndef __LDR_H__
#define __LDR_H__

#include "param.h"
// #include <gsl/gsl_randist.h>
// #include <gsl/gsl_rng.h>
#include <map>
#include <vector>

using namespace std;

class LDR {

public:
  // IO-related parameters.
  int a_mode;
  size_t d_pace;

  string file_bfile;
  string file_geno;
  string file_out;
  string path_out;

  // Summary statistics.
  size_t ni_total, ns_total; // Total number of individuals & SNPs.
  size_t ni_test, ns_test;   // Number of individuals & SNPs used
                             // for analysis
  size_t n_cvt;              // Number of covariates.

  // Indicator for individuals (phenotypes): 0 missing, 1
  // available for analysis.
  vector<int> indicator_idv;

  // Sequence indicator for SNPs: 0 ignored because of (a) maf,
  // (b) miss, (c) non-poly; 1 available for analysis.
  vector<int> indicator_snp;

  vector<SNPINFO> snpInfo; // Record SNP information.

  // Not included in PARAM.
  // gsl_rng *gsl_r;

  // Main functions.
  void CopyFromParam(PARAM &cPar);
  void CopyToParam(PARAM &cPar);

  void VB(const vector<vector<unsigned char>> &Xt, const gsl_matrix *W_gsl,
          const gsl_vector *y_gsl);
};

#endif
