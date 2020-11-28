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

#ifndef __BSLMM_H__
#define __BSLMM_H__

#include <gsl/gsl_randist.h>
// #include <gsl/gsl_rng.h>
#include <map>
#include <vector>

#include "param.h"

using namespace std;

class BSLMM {

public:
  // IO-related parameters.
  int a_mode;
  size_t d_pace;

  string file_bfile;
  string file_geno;
  string file_out;
  string path_out;

  // LMM-related parameters.
  double l_min;
  double l_max;
  size_t n_region;
  double pve_null;
  double pheno_mean;

  // BSLMM MCMC-related parameters
  double h_min, h_max, h_scale;          // Priors for h.
  double rho_min, rho_max, rho_scale;    // Priors for rho.
  double logp_min, logp_max, logp_scale; // Priors for log(pi).
  size_t s_min, s_max;                   // Min. & max. number of gammas.
  size_t w_step;                         // Number of warm up/burn in
                                         // iterations.
  size_t s_step;                         // Num. sampling iterations.
  size_t r_pace;                         // Record pace.
  size_t w_pace;                         // Write pace.
  size_t n_accept;                       // Number of acceptances.
  size_t n_mh;                           // Number of MH steps per iter.
  double geo_mean;                       // Mean of geometric dist.
  // long int randseed;
  gsl_rng *gsl_r;                        // Track randomizer state
  double trace_G;

  HYPBSLMM cHyp_initial;

  // Summary statistics.
  size_t ni_total, ns_total; // Number of total individuals and SNPs
  size_t ni_test, ns_test;   // Num. individuals & SNPs used in analysis.
  size_t n_cvt;              // Number of covariates.
  double time_UtZ;
  double time_Omega; // Time spent on optimization iterations.

  // Time spent on constructing the proposal distribution for
  // gamma (i.e. lmm or lm analysis).
  double time_Proposal;

  // Indicator for individuals (phenotypes): 0 missing, 1
  // available for analysis.
  vector<int> indicator_idv;

  // Sequence indicator for SNPs: 0 ignored because of (a) maf,
  // (b) miss, (c) non-poly; 1 available for analysis.
  vector<int> indicator_snp;

  // Record SNP information.
  vector<SNPINFO> snpInfo;

  // Not included in PARAM.
  // gsl_rng *gsl_r;
  gsl_ran_discrete_t *gsl_t;
  map<size_t, size_t> mapRank2pos;

  // Main functions.
  void CopyFromParam(PARAM &cPar);
  void CopyToParam(PARAM &cPar);

  void RidgeR(const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty,
              const gsl_vector *eval, const double lambda);

  void MCMC(const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty,
            const gsl_vector *K_eval, const gsl_vector *y);
  void WriteLog();
  void WriteLR();
  void WriteBV(const gsl_vector *bv);
  void WriteParam(vector<pair<double, double>> &beta_g, const gsl_vector *alpha,
                  const size_t w);
  void WriteParam(const gsl_vector *alpha);
  void WriteResult(const int flag, const gsl_matrix *Result_hyp,
                   const gsl_matrix *Result_gamma, const size_t w_col);

  // Subfunctions inside MCMC.
  void CalcPgamma(double *p_gammar);

  double CalcPveLM(const gsl_matrix *UtXgamma, const gsl_vector *Uty,
                   const double sigma_a2);
  void InitialMCMC(const gsl_matrix *UtX, const gsl_vector *Uty,
                   vector<size_t> &rank_old, class HYPBSLMM &cHyp,
                   vector<pair<size_t, double>> &pos_loglr);
  double CalcPosterior(const gsl_vector *Uty, const gsl_vector *K_eval,
                       gsl_vector *Utu, gsl_vector *alpha_prime,
                       class HYPBSLMM &cHyp);
  double CalcPosterior(const gsl_matrix *UtXgamma, const gsl_vector *Uty,
                       const gsl_vector *K_eval, gsl_vector *UtXb,
                       gsl_vector *Utu, gsl_vector *alpha_prime,
                       gsl_vector *beta, class HYPBSLMM &cHyp);
  void CalcCC_PVEnZ(const gsl_matrix *U, const gsl_vector *Utu,
                    gsl_vector *z_hat, class HYPBSLMM &cHyp);
  void CalcCC_PVEnZ(const gsl_matrix *U, const gsl_vector *UtXb,
                    const gsl_vector *Utu, gsl_vector *z_hat,
                    class HYPBSLMM &cHyp);
  double CalcREMLE(const gsl_matrix *Utw, const gsl_vector *Uty,
                   const gsl_vector *K_eval);

  // Calculate the maximum marginal likelihood ratio for each
  // analyzed SNPs with gemma, use it to rank SNPs.
  double CalcLR(const gsl_matrix *U, const gsl_matrix *UtX,
                const gsl_vector *Uty, const gsl_vector *K_eval,
                vector<pair<size_t, double>> &loglr_sort);
  void SampleZ(const gsl_vector *y, const gsl_vector *z_hat, gsl_vector *z);
  double ProposeHnRho(const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new,
                      const size_t &repeat);
  double ProposePi(const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new,
                   const size_t &repeat);
  double ProposeGamma(const vector<size_t> &rank_old, vector<size_t> &rank_new,
                      const double *p_gamma, const class HYPBSLMM &cHyp_old,
                      class HYPBSLMM &cHyp_new, const size_t &repeat);
  void SetXgamma(gsl_matrix *Xgamma, const gsl_matrix *X, vector<size_t> &rank);

  void CalcXtX(const gsl_matrix *X_new, const gsl_vector *y,
               const size_t s_size, gsl_matrix *XtX_new, gsl_vector *Xty_new);
  void SetXgamma(const gsl_matrix *X, const gsl_matrix *X_old,
                 const gsl_matrix *XtX_old, const gsl_vector *Xty_old,
                 const gsl_vector *y, const vector<size_t> &rank_old,
                 const vector<size_t> &rank_new, gsl_matrix *X_new,
                 gsl_matrix *XtX_new, gsl_vector *Xty_new);
  double CalcPosterior(const double yty, class HYPBSLMM &cHyp);
  double CalcPosterior(const gsl_matrix *Xgamma, const gsl_matrix *XtX,
                       const gsl_vector *Xty, const double yty,
                       const size_t s_size, gsl_vector *Xb, gsl_vector *beta,
                       class HYPBSLMM &cHyp);
  void CalcCC_PVEnZ(gsl_vector *z_hat, class HYPBSLMM &cHyp);
  void CalcCC_PVEnZ(const gsl_vector *Xb, gsl_vector *z_hat,
                    class HYPBSLMM &cHyp);
  void MCMC(const gsl_matrix *X, const gsl_vector *y);
};

#endif
