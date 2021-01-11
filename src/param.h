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

#ifndef __PARAM_H__
#define __PARAM_H__

#include "debug.h"
#include "gsl/gsl_matrix.h"
#include <gsl/gsl_rng.h>
#include "gsl/gsl_vector.h"
#include <map>
#include <set>
#include <vector>

#define K_BATCH_SIZE 20000 // #snps used for batched K
#define DEFAULT_PACE 1000  // for display only

using namespace std;

class SNPINFO {
public:
  string chr;
  string rs_number;
  double cM;
  long int base_position;
  string a_minor;
  string a_major;
  size_t n_miss;
  double missingness;
  double maf;
  size_t n_idv;         // Number of non-missing individuals.
  size_t n_nb;          // Number of neighbours on the right hand side.
  size_t file_position; // SNP location in file.
};

// Results for LMM.
class SUMSTAT {
public:
  double beta;         // REML estimator for beta.
  double se;           // SE for beta.
  double lambda_remle; // REML estimator for lambda.
  double lambda_mle;   // MLE estimator for lambda.
  double p_wald;       // p value from a Wald test.
  double p_lrt;        // p value from a likelihood ratio test.
  double p_score;      // p value from a score test.
  double logl_H1;      // log likelihood under the alternative
                       // hypothesis as a measure of goodness of fit,
                       // see https://github.com/genetics-statistics/GEMMA/issues/81
};

// Results for mvLMM.
class MPHSUMSTAT {
public:
  vector<double> v_beta;  // REML estimator for beta.
  double p_wald;          // p value from a Wald test.
  double p_lrt;           // p value from a likelihood ratio test.
  double p_score;         // p value from a score test.
  vector<double> v_Vg;    // Estimator for Vg, right half.
  vector<double> v_Ve;    // Estimator for Ve, right half.
  vector<double> v_Vbeta; // Estimator for Vbeta, right half.
};

// Hyper-parameters for BSLMM.
class HYPBSLMM {
public:
  double h;
  double pve;
  double rho;
  double pge;
  double logp;
  size_t n_gamma;
};

// Header class.
class HEADER {
public:
  size_t rs_col;
  size_t chr_col;
  size_t pos_col;
  size_t cm_col;
  size_t a1_col;
  size_t a0_col;
  size_t z_col;
  size_t beta_col;
  size_t sebeta_col;
  size_t chisq_col;
  size_t p_col;
  size_t n_col;
  size_t nmis_col;
  size_t nobs_col;
  size_t ncase_col;
  size_t ncontrol_col;
  size_t af_col;
  size_t var_col;
  size_t ws_col;
  size_t cor_col;
  size_t coln; // Number of columns.
  set<size_t> catc_col;
  set<size_t> catd_col;
};

class PARAM {
public:
  // IO-related parameters
  // bool mode_check = true;   // run data checks (slower)
  // bool mode_strict = false; // exit on some data checks
  // bool mode_silence;
  // bool mode_debug = false;
  // uint issue; // enable tests for issue on github tracker

  uint a_mode; // Analysis mode, 1/2/3/4 for Frequentist tests
  int k_mode; // Kinship read mode: 1: n by n matrix, 2: id/id/k_value;
  vector<size_t> p_column; // Which phenotype column needs analysis.
  size_t d_pace = DEFAULT_PACE;   // Display pace (-pace switch)

  string file_bfile, file_mbfile;
  string file_geno, file_mgeno;
  string file_pheno;
  string file_anno; // Optional.
  string file_gxe;  // Optional.
  string file_cvt;  // Optional.
  string file_cat, file_mcat;
  string file_catc, file_mcatc;
  string file_var;
  string file_beta;
  string file_cor;
  string file_kin, file_mk;
  string file_ku, file_kd;
  string file_study, file_mstudy;
  string file_ref, file_mref;
  string file_weight, file_wsnp, file_wcat;
  string file_out;
  string file_bf, file_hyp;
  string path_out;

  string file_epm;     // Estimated parameter file.
  string file_ebv;     // Estimated breeding value file.
  string file_log;     // Log file containing mean estimate.
  string file_read;    // File containing total number of reads.
  string file_gene;    // Gene expression file.
  string file_snps;    // File containing analyzed SNPs or genes.
  string file_ksnps;   // File SNPs for computing K
  string file_gwasnps; // File SNPs for computing GWAS

  // QC-related parameters.
  double miss_level;
  double maf_level;
  double hwe_level;
  double r2_level;

  // LMM-related parameters.
  string loco;
  double l_min;
  double l_max;
  size_t n_region;
  double l_mle_null, l_remle_null;
  double logl_mle_H0, logl_remle_H0;
  double pve_null, pve_se_null, pve_total, se_pve_total;
  double vg_remle_null, ve_remle_null, vg_mle_null, ve_mle_null;
  vector<double> Vg_remle_null, Ve_remle_null, Vg_mle_null, Ve_mle_null;
  vector<double> VVg_remle_null, VVe_remle_null, VVg_mle_null;
  vector<double> VVe_mle_null;
  vector<double> beta_remle_null, se_beta_remle_null, beta_mle_null;
  vector<double> se_beta_mle_null;
  double p_nr;
  double em_prec, nr_prec;
  size_t em_iter, nr_iter;
  size_t crt;
  double pheno_mean; // Phenotype mean from BSLMM fitting or prediction.

  // For fitting multiple variance components.
  // The first 3 are of size (n_vc), and the next 2 are of size n_vc+1.
  bool noconstrain;
  vector<double> v_traceG;
  vector<double> v_pve;
  vector<double> v_se_pve;

  vector<double> v_sigma2;
  vector<double> v_se_sigma2;
  vector<double> v_enrich;
  vector<double> v_se_enrich;
  vector<double> v_beta;
  vector<double> v_se_beta;

  // BSLMM/MCMC-related parameters.
  double h_min, h_max, h_scale;          // Priors for h.
  double rho_min, rho_max, rho_scale;    // Priors for rho.
  double logp_min, logp_max, logp_scale; // Priors for log(pi).
  size_t h_ngrid, rho_ngrid;
  size_t s_min, s_max; // Min & max. number of gammas.
  size_t w_step;       // # warm up/burn in iter.
  size_t s_step;       // # sampling iterations.
  size_t r_pace;       // Record pace.
  size_t w_pace;       // Write pace.
  size_t n_accept;     // Number of acceptance.
  size_t n_mh;         // # MH steps in each iter.
  double geo_mean;     // Mean of geometric dist.
  long int randseed;   // holds -seed parameter
  gsl_rng *gsl_r = NULL;      // Track the randomizer
  double trace_G;

  HYPBSLMM cHyp_initial;

  // VARCOV-related parameters.
  double window_cm;
  size_t window_bp;
  size_t window_ns;

  // vc-related parameters.
  size_t n_block;

  // Summary statistics.
  bool error;

  // Number of individuals.
  size_t ni_total, ni_test, ni_cvt, ni_study, ni_ref;
  size_t ni_max = 0; // -nind switch for testing purposes

  // Number of observed and missing phenotypes.
  size_t np_obs, np_miss;

  // Number of SNPs.
  size_t ns_total, ns_test, ns_study, ns_ref;

  size_t ng_total, ng_test;   // Number of genes.
  size_t ni_control, ni_case; // Number of controls and number of cases.
  size_t ni_subsample;        // Number of subsampled individuals.
  size_t n_cvt;               // Number of covariates.
  size_t n_cat;               // Number of continuous categories.
  size_t n_ph;                // Number of phenotypes.
  size_t n_vc;                // Number of variance components
                              // (including the diagonal matrix).
  double time_total;          // Record total time.
  double time_G;              // Time spent on reading files the
                              // second time and calculate K.
  double time_eigen;          // Time spent on eigen-decomposition.
  double time_UtX;            // Time spent on calculating UX and Uy.
  double time_UtZ;            // Time calculating UtZ for probit BSLMM.
  double time_opt;            // Time on optimization iterations/MCMC.
  double time_Omega;          // Time spent on calculating Omega.
  double time_hyp;            // Time sampling hyperparameters in PMM.
  double time_Proposal;       // Time spent on constructing the
                              // proposal distribution (i.e. the
                              // initial LMM or LM analysis).

  // Data.
  // Vector recording all phenotypes (NA replaced with -9).
  vector<vector<double>> pheno;

  // Vector recording all covariates (NA replaced with -9).
  vector<vector<double>> cvt;

  // Vector recording all covariates (NA replaced with -9).
  vector<double> gxe;

  // Vector recording weights for the individuals, which is
  // useful for animal breeding studies.
  vector<double> weight;

  // Matrix recording when a phenotype is missing for an
  // individual; 0 missing, 1 available.
  vector<vector<int>> indicator_pheno;

  // Indicator for individuals (phenotypes): 0 missing, 1
  // available for analysis
  vector<int> indicator_idv;

  // Sequence indicator for SNPs: 0 ignored because of (a) maf,
  // (b) miss, (c) non-poly; 1 available for analysis.
  vector<int> indicator_snp;

  // Sequence indicator for SNPs: 0 ignored because of (a) maf,
  // (b) miss, (c) non-poly; 1 available for analysis.
  vector<vector<int>> mindicator_snp;

  // Indicator for covariates: 0 missing, 1 available for
  // analysis.
  vector<int> indicator_cvt;

  // Indicator for gxe: 0 missing, 1 available for analysis.
  vector<int> indicator_gxe;

  // Indicator for weight: 0 missing, 1 available for analysis.
  vector<int> indicator_weight;

  // Indicator for estimated breeding value file: 0 missing, 1
  // available for analysis.
  vector<int> indicator_bv;

  // Indicator for read file: 0 missing, 1 available for analysis.
  vector<int> indicator_read;
  vector<double> vec_read; // Total number of reads.
  vector<double> vec_bv;   // Breeding values.
  vector<size_t> est_column;

  map<string, int> mapID2num;             // Map small ID to number, 0 to n-1.
  map<string, string> mapRS2chr;          // Map rs# to chromosome location.
  map<string, long int> mapRS2bp;         // Map rs# to base position.
  map<string, double> mapRS2cM;           // Map rs# to cM.
  map<string, double> mapRS2est;          // Map rs# to parameters.
  map<string, size_t> mapRS2cat;          // Map rs# to category number.
  map<string, vector<double>> mapRS2catc; // Map rs# to cont. cat's.
  map<string, double> mapRS2wsnp;         // Map rs# to SNP weights.
  map<string, vector<double>> mapRS2wcat; // Map rs# to SNP cat weights.

  vector<SNPINFO> snpInfo;          // Record SNP information.
  vector<vector<SNPINFO>> msnpInfo; // Record SNP information.
  set<string> setSnps;              // Set of snps for analysis (-snps).
  set<string> setKSnps;             // Set of snps for K (-ksnps and LOCO)
  set<string> setGWASnps;           // Set of snps for GWA (-gwasnps and LOCO)

  // Constructor and destructor
  PARAM();
  ~PARAM();

  // Functions.
  void ReadFiles();
  void CheckParam();
  void CheckData();
  void PrintSummary();
  void ReadGenotypes(gsl_matrix *UtX, gsl_matrix *K, const bool calc_K);
  void ReadGenotypes(vector<vector<unsigned char>> &Xt, gsl_matrix *K,
                     const bool calc_K);
  void CheckCvt();
  void CopyCvt(gsl_matrix *W);
  void CopyA(size_t flag, gsl_matrix *A);
  void CopyGxe(gsl_vector *gxe);
  void CopyWeight(gsl_vector *w);
  void ProcessCvtPhen();
  void CopyCvtPhen(gsl_matrix *W, gsl_vector *y, size_t flag);
  void CopyCvtPhen(gsl_matrix *W, gsl_matrix *Y, size_t flag);
  void CalcKin(gsl_matrix *matrix_kin);
  void CalcS(const map<string, double> &mapRS2wA,
             const map<string, double> &mapRS2wK, const gsl_matrix *W,
             gsl_matrix *A, gsl_matrix *K, gsl_matrix *S, gsl_matrix *Svar,
             gsl_vector *ns);
  void WriteVector(const gsl_vector *q, const gsl_vector *s,
                   const size_t n_total, const string suffix);
  void WriteVar(const string suffix);
  void WriteMatrix(const gsl_matrix *matrix_U, const string suffix);
  void WriteVector(const gsl_vector *vector_D, const string suffix);
  void CopyRead(gsl_vector *log_N);
  void ObtainWeight(const set<string> &setSnps_beta,
                    map<string, double> &mapRS2wK);
  void UpdateWeight(const size_t pve_flag, const map<string, double> &mapRS2wK,
                    const size_t ni_test, const gsl_vector *ns,
                    map<string, double> &mapRS2wA);
  void UpdateSNPnZ(const map<string, double> &mapRS2wA,
                   const map<string, string> &mapRS2A1,
                   const map<string, double> &mapRS2z, gsl_vector *w,
                   gsl_vector *z, vector<size_t> &vec_cat);
  void UpdateSNP(const map<string, double> &mapRS2wA);
};

size_t GetabIndex(const size_t a, const size_t b, const size_t n_cvt);

#endif
