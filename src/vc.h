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

#ifndef __VC_H__
#define __VC_H__

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gemma_io.h"
#include "param.h"

using namespace std;

class VC_PARAM {

public:
  const gsl_matrix *K;
  const gsl_matrix *W;
  const gsl_vector *y;
  gsl_matrix *P;
  gsl_vector *Py;
  gsl_matrix *KPy_mat;
  gsl_matrix *PKPy_mat;
  gsl_matrix *Hessian;
  bool noconstrain;
};

// Methods for variant component (VC) estimation

class VC {

public:
  // IO-related parameters
  size_t a_mode;
  string file_cat;
  string file_beta;
  string file_cor;
  string file_mq;
  string file_ms;

  string file_out;
  string path_out;

  set<string> setSnps;

  size_t ni_total_ref, ns_total_ref, ns_pair;
  size_t ni_total, ns_total, ns_test;
  size_t n_vc;

  double pve_total, se_pve_total;
  vector<double> v_sigma2;
  vector<double> v_se_sigma2;
  vector<double> v_pve;
  vector<double> v_se_pve;
  vector<double> v_traceG;
  vector<double> v_beta;
  vector<double> v_se_beta;

  size_t crt;
  double window_cm, window_bp, window_ns;

  double time_UtX;
  double time_opt;

  // Main functions.
  void CopyFromParam(PARAM &cPar);
  void CopyToParam(PARAM &cPar);
  void WriteFile_qs(const gsl_vector *s_vec, const gsl_vector *q_vec,
                    const gsl_vector *qvar_vec, const gsl_matrix *S_mat,
                    const gsl_matrix *Svar_mat);
  void CalcVChe(const gsl_matrix *K, const gsl_matrix *W, const gsl_vector *y);
  void CalcVCreml(const bool noconstrain, const gsl_matrix *K,
                  const gsl_matrix *W, const gsl_vector *y);
  void CalcVCacl(const gsl_matrix *K, const gsl_matrix *W, const gsl_vector *y);
};

void CalcVCss(const gsl_matrix *Vq, const gsl_matrix *S_mat,
              const gsl_matrix *Svar_mat, const gsl_vector *q_vec,
              const gsl_vector *s_vec, const double df, vector<double> &v_pve,
              vector<double> &v_se_pve, double &pve_total, double &se_pve_total,
              vector<double> &v_sigma2, vector<double> &v_se_sigma2,
              vector<double> &v_enrich, vector<double> &v_se_enrich);

bool BimbamXwz(const string &file_geno, const int display_pace,
               vector<int> &indicator_idv, vector<int> &indicator_snp,
               const vector<size_t> &vec_cat, const gsl_vector *w,
               const gsl_vector *z, size_t ns_test, gsl_matrix *XWz);
bool PlinkXwz(const string &file_bed, const int display_pace,
              vector<int> &indicator_idv, vector<int> &indicator_snp,
              const vector<size_t> &vec_cat, const gsl_vector *w,
              const gsl_vector *z, size_t ns_test, gsl_matrix *XWz);
bool MFILEXwz(const size_t mfile_mode, const string &file_mfile,
              const int display_pace, vector<int> &indicator_idv,
              vector<vector<int>> &mindicator_snp,
              const vector<size_t> &vec_cat, const gsl_vector *w,
              const gsl_vector *z, gsl_matrix *XWz);

bool BimbamXtXwz(const string &file_geno, const int display_pace,
                 vector<int> &indicator_idv, vector<int> &indicator_snp,
                 const gsl_matrix *XWz, size_t ns_test, gsl_matrix *XtXWz);
bool PlinkXtXwz(const string &file_bed, const int display_pace,
                vector<int> &indicator_idv, vector<int> &indicator_snp,
                const gsl_matrix *XWz, size_t ns_test, gsl_matrix *XtXWz);
bool MFILEXtXwz(const size_t mfile_mode, const string &file_mfile,
                const int display_pace, vector<int> &indicator_idv,
                vector<vector<int>> &mindicator_snp, const gsl_matrix *XWz,
                gsl_matrix *XtXWz);

void CalcCIss(const gsl_matrix *Xz, const gsl_matrix *XWz,
              const gsl_matrix *XtXWz, const gsl_matrix *S_mat,
              const gsl_matrix *Svar_mat, const gsl_vector *w,
              const gsl_vector *z, const gsl_vector *s_vec,
              const vector<size_t> &vec_cat, const vector<double> &v_pve,
              vector<double> &v_se_pve, double &pve_total, double &se_pve_total,
              vector<double> &v_sigma2, vector<double> &v_se_sigma2,
              vector<double> &v_enrich, vector<double> &v_se_enrich);

#endif
