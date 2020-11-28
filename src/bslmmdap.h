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

#ifndef __BSLMMDAP_H__
#define __BSLMMDAP_H__

#include "param.h"
#include <gsl/gsl_randist.h>
// #include <gsl/gsl_rng.h>
#include <map>
#include <vector>

using namespace std;

class BSLMMDAP {

public:
  // IO-related parameters.
  int a_mode;
  size_t d_pace;

  string file_bfile;
  string file_geno;
  string file_out;
  string path_out;

  // LMM related parameters
  double pve_null;
  double pheno_mean;

  // BSLMM MCMC related parameters
  // long int randseed;
  double trace_G;

  HYPBSLMM cHyp_initial;

  // Summary statistics
  size_t ni_total, ns_total; // Number of total individuals and SNPs.
  size_t ni_test, ns_test;   // Number of individuals and SNPs
                             // used for analysis.

  double h_min, h_max, rho_min, rho_max;
  size_t h_ngrid, rho_ngrid;

  double time_UtZ;
  double time_Omega;    // Time spent on optimization iterations.
  double time_Proposal; // Time spent on constructing the
                        // proposal distribution for gamma
                        // (i.e., lmm or lm analysis).

  // Indicator for individuals (phenotypes): 0 missing, 1
  // available for analysis.
  vector<int> indicator_idv;

  // Sequence indicator for SNPs: 0 ignored because of (a) maf,
  // (b) miss, (c) non-poly; 1 available for analysis.
  vector<int> indicator_snp;

  vector<SNPINFO> snpInfo; // Record SNP information.

  // Main functions.
  void CopyFromParam(PARAM &cPar);
  void CopyToParam(PARAM &cPar);

  void WriteResult(const gsl_matrix *Hyper, const gsl_matrix *BF);
  void WriteResult(const vector<string> &vec_rs, const gsl_matrix *Hyper,
                   const gsl_vector *pip, const gsl_vector *coef);
  double CalcMarginal(const gsl_vector *Uty, const gsl_vector *K_eval,
                      const double sigma_b2, const double tau);
  double CalcMarginal(const gsl_matrix *UtXgamma, const gsl_vector *Uty,
                      const gsl_vector *K_eval, const double sigma_a2,
                      const double sigma_b2, const double tau);
  double CalcPrior(class HYPBSLMM &cHyp);

  void DAP_CalcBF(const gsl_matrix *U, const gsl_matrix *UtX,
                  const gsl_vector *Uty, const gsl_vector *K_eval,
                  const gsl_vector *y);
  void
  DAP_EstimateHyper(const size_t kc, const size_t kd,
                    const vector<string> &vec_rs, const vector<double> &vec_sa2,
                    const vector<double> &vec_sb2, const vector<double> &wab,
                    const vector<vector<vector<double>>> &BF, gsl_matrix *Ac,
                    gsl_matrix_int *Ad, gsl_vector_int *dlevel);
};

void ReadFile_hyb(const string &file_hyp, vector<double> &vec_sa2,
                  vector<double> &vec_sb2, vector<double> &vec_wab);
void ReadFile_bf(const string &file_bf, vector<string> &vec_rs,
                 vector<vector<vector<double>>> &BF);
void ReadFile_cat(const string &file_cat, const vector<string> &vec_rs,
                  gsl_matrix *Ac, gsl_matrix_int *Ad, gsl_vector_int *dlevel,
                  size_t &kc, size_t &kd);

#endif
