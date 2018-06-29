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

#ifndef __MVLMM_H__
#define __MVLMM_H__

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gemma_io.h"
#include "param.h"

using namespace std;

class MVLMM {

public:
  // IO-related parameters.
  int a_mode;    // Analysis mode: 1/2/3/4 for Frequentist tests.
  size_t d_pace; // Display pace.

  string file_bfile;
  string file_geno;
  string file_oxford;
  string file_out;
  string path_out;

  // MVLMM-related parameters.
  double l_min;
  double l_max;
  size_t n_region;
  double logl_remle_H0, logl_mle_H0;
  vector<double> Vg_remle_null, Ve_remle_null, Vg_mle_null, Ve_mle_null;
  vector<double> VVg_remle_null, VVe_remle_null, VVg_mle_null;
  vector<double> VVe_mle_null;
  vector<double> beta_remle_null, se_beta_remle_null, beta_mle_null;
  vector<double> se_beta_mle_null;
  double p_nr;
  size_t em_iter, nr_iter;
  double em_prec, nr_prec;
  size_t crt;

  // Summary statistics.
  size_t ni_total, ni_test; // Number of individuals.
  size_t ns_total, ns_test; // Number of SNPs.
  size_t n_cvt;
  size_t n_ph;
  double time_UtX; // Time spent on optimization iterations.
  double time_opt; // Time spent on optimization iterations.

  // Indicator for individuals (phenotypes): 0 missing, 1
  // available for analysis.
  vector<int> indicator_idv;

  // Sequence indicator for SNPs: 0 ignored because of (a) maf,
  // (b) miss, (c) non-poly; 1 available for analysis.
  vector<int> indicator_snp;

  vector<SNPINFO> snpInfo; // Record SNP information.

  // Not included in PARAM.
  vector<MPHSUMSTAT> sumStat; // Output SNPSummary Data.

  // Main functions
  void CopyFromParam(PARAM &cPar);
  void CopyToParam(PARAM &cPar);
  void AnalyzeBimbam(const gsl_matrix *U, const gsl_vector *eval,
                     const gsl_matrix *UtW, const gsl_matrix *UtY);
  void AnalyzePlink(const gsl_matrix *U, const gsl_vector *eval,
                    const gsl_matrix *UtW, const gsl_matrix *UtY);
  void Analyzebgen(const gsl_matrix *U, const gsl_vector *eval,
                   const gsl_matrix *UtW, const gsl_matrix *UtY);
  void AnalyzeBimbamGXE(const gsl_matrix *U, const gsl_vector *eval,
                        const gsl_matrix *UtW, const gsl_matrix *UtY,
                        const gsl_vector *env);
  void AnalyzePlinkGXE(const gsl_matrix *U, const gsl_vector *eval,
                       const gsl_matrix *UtW, const gsl_matrix *UtY,
                       const gsl_vector *env);
  void WriteFiles();
};

void CalcMvLmmVgVeBeta(const gsl_vector *eval, const gsl_matrix *UtW,
                       const gsl_matrix *UtY, const size_t em_iter,
                       const size_t nr_iter, const double em_prec,
                       const double nr_prec, const double l_min,
                       const double l_max, const size_t n_region,
                       gsl_matrix *V_g, gsl_matrix *V_e, gsl_matrix *B,
                       gsl_matrix *se_B);

#endif
