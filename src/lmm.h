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

#ifndef __LMM_H__
#define __LMM_H__

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "io.h"
#include "param.h"

using namespace std;

#define LMM_BATCH_SIZE 10000 // used for batch processing

class FUNC_PARAM {

public:
  bool calc_null;
  size_t ni_test;
  size_t n_cvt;
  const gsl_vector *eval;
  const gsl_matrix *Uab;
  const gsl_vector *ab;
  size_t e_mode;
};

class LMM {

public:
  // IO-related parameters
  int a_mode;    // Analysis mode: 1/2/3/4 for Frequentist tests.
  size_t d_pace; // Display pace.

  string file_bfile;
  string file_geno;
  string file_out;
  string path_out;

  string file_gene;

  // LMM related parameters
  double l_min;
  double l_max;
  size_t n_region;
  double l_mle_null;
  double logl_mle_H0;

  // Summary statistics
  size_t ni_total, ni_test; // Number of individuals.
  size_t ns_total, ns_test; // Number of SNPs.
  size_t ng_total, ng_test; // Number of genes.
  size_t n_cvt;
  double time_UtX; // Time spent on optimization iterations.
  double time_opt; // Time spent on optimization iterations.

  // Indicator for individuals (phenotypes): 0 missing, 1
  // available for analysis.
  vector<int> indicator_idv;

  // Sequence indicator for SNPs: 0 ignored because of (a) maf,
  // (b) miss, (c) non-poly; 1 available for analysis.
  vector<int> indicator_snp;

  vector<SNPINFO> snpInfo; // Record SNP information.
  set<string> setGWASnps;  // Record SNP information.

  // Not included in PARAM.
  vector<SUMSTAT> sumStat; // Output SNPSummary Data.

  // Main functions.
  void CopyFromParam(PARAM &cPar);
  void CopyToParam(PARAM &cPar);
  void AnalyzeGene(const gsl_matrix *U, const gsl_vector *eval,
                   const gsl_matrix *UtW, const gsl_vector *Utx,
                   const gsl_matrix *W, const gsl_vector *x);
  void AnalyzePlink(const gsl_matrix *U, const gsl_vector *eval,
                    const gsl_matrix *UtW, const gsl_vector *Uty,
                    const gsl_matrix *W, const gsl_vector *y);
  void AnalyzeBimbam(const gsl_matrix *U, const gsl_vector *eval,
                     const gsl_matrix *UtW, const gsl_vector *Uty,
                     const gsl_matrix *W, const gsl_vector *y,
                     const set<string> gwasnps);
  void AnalyzePlinkGXE(const gsl_matrix *U, const gsl_vector *eval,
                       const gsl_matrix *UtW, const gsl_vector *Uty,
                       const gsl_matrix *W, const gsl_vector *y,
                       const gsl_vector *env);
  void AnalyzeBimbamGXE(const gsl_matrix *U, const gsl_vector *eval,
                        const gsl_matrix *UtW, const gsl_vector *Uty,
                        const gsl_matrix *W, const gsl_vector *y,
                        const gsl_vector *env);
  void WriteFiles();

  void CalcRLWald(const double &lambda, const FUNC_PARAM &params, double &beta,
                  double &se, double &p_wald);
  void CalcRLScore(const double &l, const FUNC_PARAM &params, double &beta,
                   double &se, double &p_score);
};

void MatrixCalcLR(const gsl_matrix *U, const gsl_matrix *UtX,
                  const gsl_vector *Uty, const gsl_vector *K_eval,
                  const double l_min, const double l_max, const size_t n_region,
                  vector<pair<size_t, double>> &pos_loglr);
void CalcLambda(const char func_name, FUNC_PARAM &params, const double l_min,
                const double l_max, const size_t n_region, double &lambda,
                double &logf);
void CalcLambda(const char func_name, const gsl_vector *eval,
                const gsl_matrix *UtW, const gsl_vector *Uty,
                const double l_min, const double l_max, const size_t n_region,
                double &lambda, double &logl_H0);
void CalcPve(const gsl_vector *eval, const gsl_matrix *UtW,
             const gsl_vector *Uty, const double lambda, const double trace_G,
             double &pve, double &pve_se);
void CalcLmmVgVeBeta(const gsl_vector *eval, const gsl_matrix *UtW,
                     const gsl_vector *Uty, const double lambda, double &vg,
                     double &ve, gsl_vector *beta, gsl_vector *se_beta);

#endif
