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

#ifndef __VARCOV_H__
#define __VARCOV_H__

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gemma_io.h"
#include "param.h"

using namespace std;

class VARCOV {

public:
  // IO-related parameters.
  string file_out;
  string path_out;
  string file_geno;
  string file_bfile;
  int d_pace;

  vector<int> indicator_idv;
  vector<int> indicator_snp;

  vector<SNPINFO> snpInfo;

  double time_opt;

  // Class-specific parameters.
  double window_cm;
  size_t window_bp;
  size_t window_ns;

  // Main functions.
  void CopyFromParam(PARAM &cPar);
  void CopyToParam(PARAM &cPar);
  void CalcNB(vector<SNPINFO> &snpInfo_sort);
  void WriteCov(const int flag, const vector<SNPINFO> &snpInfo_sub,
                const vector<vector<double>> &Cov_mat);
  void AnalyzeBimbam();
  void AnalyzePlink();
};

#endif
