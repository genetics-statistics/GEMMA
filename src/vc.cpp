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

#include <fstream>
#include <iostream>
#include <sstream>

#include <bitset>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"

#include "gsl/gsl_cdf.h"
#include "gsl/gsl_min.h"
#include "gsl/gsl_multiroots.h"

// #include "Eigen/Dense"

#include "gzstream.h"
#include "gemma_io.h"
#include "lapack.h"
#include "lmm.h"
#include "mathfunc.h"
#include "param.h"
#include "vc.h"
#include "fastblas.h"

using namespace std;
// using namespace Eigen;

// In this file, X, Y are already transformed (i.e. UtX and UtY).
void VC::CopyFromParam(PARAM &cPar) {
  a_mode = cPar.a_mode;

  file_cat = cPar.file_cat;
  file_beta = cPar.file_beta;
  file_cor = cPar.file_cor;

  setSnps = cPar.setSnps;

  file_out = cPar.file_out;
  path_out = cPar.path_out;

  time_UtX = 0.0;
  time_opt = 0.0;

  v_traceG = cPar.v_traceG;

  ni_total = cPar.ni_total;
  ns_total = cPar.ns_total;
  ns_test = cPar.ns_test;

  crt = cPar.crt;
  window_cm = cPar.window_cm;
  window_bp = cPar.window_bp;
  window_ns = cPar.window_ns;

  n_vc = cPar.n_vc;

  return;
}

void VC::CopyToParam(PARAM &cPar) {
  cPar.time_UtX = time_UtX;
  cPar.time_opt = time_opt;

  cPar.v_pve = v_pve;
  cPar.v_se_pve = v_se_pve;
  cPar.v_sigma2 = v_sigma2;
  cPar.v_se_sigma2 = v_se_sigma2;
  cPar.pve_total = pve_total;
  cPar.se_pve_total = se_pve_total;
  cPar.v_traceG = v_traceG;

  cPar.v_beta = v_beta;
  cPar.v_se_beta = v_se_beta;

  cPar.ni_total = ni_total;
  cPar.ns_total = ns_total;
  cPar.ns_test = ns_test;

  cPar.n_vc = n_vc;

  return;
}

void VC::WriteFile_qs(const gsl_vector *s_vec, const gsl_vector *q_vec,
                      const gsl_vector *qvar_vec, const gsl_matrix *S_mat,
                      const gsl_matrix *Svar_mat) {
  string file_str;
  file_str = path_out + "/" + file_out;
  file_str += ".qvec.txt";

  ofstream outfile_q(file_str.c_str(), ofstream::out);
  if (!outfile_q) {
    cout << "error writing file: " << file_str.c_str() << endl;
    return;
  }

  for (size_t i = 0; i < s_vec->size; i++) {
    outfile_q << gsl_vector_get(s_vec, i) << endl;
  }
  for (size_t i = 0; i < q_vec->size; i++) {
    outfile_q << gsl_vector_get(q_vec, i) << endl;
  }
  for (size_t i = 0; i < qvar_vec->size; i++) {
    outfile_q << gsl_vector_get(qvar_vec, i) << endl;
  }

  outfile_q.clear();
  outfile_q.close();

  file_str = path_out + "/" + file_out;
  file_str += ".smat.txt";

  ofstream outfile_s(file_str.c_str(), ofstream::out);
  if (!outfile_s) {
    cout << "error writing file: " << file_str.c_str() << endl;
    return;
  }

  for (size_t i = 0; i < S_mat->size1; i++) {
    for (size_t j = 0; j < S_mat->size2; j++) {
      outfile_s << gsl_matrix_get(S_mat, i, j) << "\t";
    }
    outfile_s << endl;
  }
  for (size_t i = 0; i < Svar_mat->size1; i++) {
    for (size_t j = 0; j < Svar_mat->size2; j++) {
      outfile_s << gsl_matrix_get(Svar_mat, i, j) << "\t";
    }
    outfile_s << endl;
  }

  outfile_s.clear();
  outfile_s.close();

  return;
}

void UpdateParam(const gsl_vector *log_sigma2, VC_PARAM *p) {
  size_t n1 = (p->K)->size1, n_vc = log_sigma2->size - 1, n_cvt = (p->W)->size2;

  gsl_matrix *K_temp = gsl_matrix_alloc(n1, n1);
  gsl_matrix *HiW = gsl_matrix_alloc(n1, n_cvt);
  gsl_matrix *WtHiW = gsl_matrix_alloc(n_cvt, n_cvt);
  gsl_matrix *WtHiWi = gsl_matrix_alloc(n_cvt, n_cvt);
  gsl_matrix *WtHiWiWtHi = gsl_matrix_alloc(n_cvt, n1);

  double sigma2;

  // Calculate H = \sum_i^{k+1} \sigma_i^2 K_i.
  gsl_matrix_set_zero(p->P);
  for (size_t i = 0; i < n_vc + 1; i++) {
    if (i == n_vc) {
      gsl_matrix_set_identity(K_temp);
    } else {
      gsl_matrix_const_view K_sub =
          gsl_matrix_const_submatrix(p->K, 0, n1 * i, n1, n1);
      gsl_matrix_memcpy(K_temp, &K_sub.matrix);
    }

    // When unconstrained, update on sigma2 instead of log_sigma2.
    if (p->noconstrain) {
      sigma2 = gsl_vector_get(log_sigma2, i);
    } else {
      sigma2 = exp(gsl_vector_get(log_sigma2, i));
    }
    gsl_matrix_scale(K_temp, sigma2);
    gsl_matrix_add(p->P, K_temp);
  }

  // Calculate H^{-1}.
  fast_inverse(p->P);

  fast_dgemm("N", "N", 1.0, p->P, p->W, 0.0, HiW);
  fast_dgemm("T", "N", 1.0, p->W, HiW, 0.0, WtHiW);

  fast_inverse(WtHiW);
  gsl_matrix_memcpy(WtHiWi, WtHiW);

  fast_dgemm("N", "T", 1.0, WtHiWi, HiW, 0.0, WtHiWiWtHi);
  fast_dgemm("N", "N", -1.0, HiW, WtHiWiWtHi, 1.0, p->P);

  // Calculate Py, KPy, PKPy.
  gsl_blas_dgemv(CblasNoTrans, 1.0, p->P, p->y, 0.0, p->Py);

  double d;
  for (size_t i = 0; i < n_vc + 1; i++) {
    gsl_vector_view KPy = gsl_matrix_column(p->KPy_mat, i);
    gsl_vector_view PKPy = gsl_matrix_column(p->PKPy_mat, i);

    if (i == n_vc) {
      gsl_vector_memcpy(&KPy.vector, p->Py);
    } else {
      gsl_matrix_const_view K_sub =
          gsl_matrix_const_submatrix(p->K, 0, n1 * i, n1, n1);

      // Seems to be important to use gsl dgemv here instead of
      // fast_dgemv; otherwise.
      gsl_blas_dgemv(CblasNoTrans, 1.0, &K_sub.matrix, p->Py, 0.0, &KPy.vector);
    }

    gsl_blas_dgemv(CblasNoTrans, 1.0, p->P, &KPy.vector, 0.0, &PKPy.vector);

    // When phenotypes are not normalized well, then some values in
    // the following matrix maybe NaN; change that to 0; this seems to
    // only happen when fast_dgemv was used above.
    for (size_t j = 0; j < p->KPy_mat->size1; j++) {
      d = gsl_matrix_get(p->KPy_mat, j, i);
      if (isnan(d)) {
        gsl_matrix_set(p->KPy_mat, j, i, 0);
        cout << "nan appears in " << i << " " << j << endl;
      }
      d = gsl_matrix_get(p->PKPy_mat, j, i);
      if (isnan(d)) {
        gsl_matrix_set(p->PKPy_mat, j, i, 0);
        cout << "nan appears in " << i << " " << j << endl;
      }
    }
  }

  gsl_matrix_free(K_temp);
  gsl_matrix_free(HiW);
  gsl_matrix_free(WtHiW);
  gsl_matrix_free(WtHiWi);
  gsl_matrix_free(WtHiWiWtHi);

  return;
}

// Below are functions for AI algorithm.
int LogRL_dev1(const gsl_vector *log_sigma2, void *params, gsl_vector *dev1) {
  VC_PARAM *p = (VC_PARAM *)params;

  size_t n1 = (p->K)->size1, n_vc = log_sigma2->size - 1;

  double tr, d;

  // Update parameters.
  UpdateParam(log_sigma2, p);

  // Calculate dev1=-0.5*trace(PK_i)+0.5*yPKPy.
  for (size_t i = 0; i < n_vc + 1; i++) {
    if (i == n_vc) {
      tr = 0;
      for (size_t l = 0; l < n1; l++) {
        tr += gsl_matrix_get(p->P, l, l);
      }
    } else {
      tr = 0;
      for (size_t l = 0; l < n1; l++) {
        gsl_vector_view P_row = gsl_matrix_row(p->P, l);
        gsl_vector_const_view K_col = gsl_matrix_const_column(p->K, n1 * i + l);
        gsl_blas_ddot(&P_row.vector, &K_col.vector, &d);
        tr += d;
      }
    }

    gsl_vector_view KPy_i = gsl_matrix_column(p->KPy_mat, i);
    gsl_blas_ddot(p->Py, &KPy_i.vector, &d);

    if (p->noconstrain) {
      d = (-0.5 * tr + 0.5 * d);
    } else {
      d = (-0.5 * tr + 0.5 * d) * exp(gsl_vector_get(log_sigma2, i));
    }

    gsl_vector_set(dev1, i, d);
  }

  return GSL_SUCCESS;
}

int LogRL_dev2(const gsl_vector *log_sigma2, void *params, gsl_matrix *dev2) {
  VC_PARAM *p = (VC_PARAM *)params;

  size_t n_vc = log_sigma2->size - 1;

  double d, sigma2_i, sigma2_j;

  // Update parameters.
  UpdateParam(log_sigma2, p);

  // Calculate dev2 = 0.5(yPKPKPy).
  for (size_t i = 0; i < n_vc + 1; i++) {
    gsl_vector_view KPy_i = gsl_matrix_column(p->KPy_mat, i);
    if (p->noconstrain) {
      sigma2_i = gsl_vector_get(log_sigma2, i);
    } else {
      sigma2_i = exp(gsl_vector_get(log_sigma2, i));
    }

    for (size_t j = i; j < n_vc + 1; j++) {
      gsl_vector_view PKPy_j = gsl_matrix_column(p->PKPy_mat, j);

      gsl_blas_ddot(&KPy_i.vector, &PKPy_j.vector, &d);
      if (p->noconstrain) {
        sigma2_j = gsl_vector_get(log_sigma2, j);
        d *= -0.5;
      } else {
        sigma2_j = exp(gsl_vector_get(log_sigma2, j));
        d *= -0.5 * sigma2_i * sigma2_j;
      }

      gsl_matrix_set(dev2, i, j, d);
      if (j != i) {
        gsl_matrix_set(dev2, j, i, d);
      }
    }
  }

  gsl_matrix_memcpy(p->Hessian, dev2);
  return GSL_SUCCESS;
}

int LogRL_dev12(const gsl_vector *log_sigma2, void *params, gsl_vector *dev1,
                gsl_matrix *dev2) {
  VC_PARAM *p = (VC_PARAM *)params;

  size_t n1 = (p->K)->size1, n_vc = log_sigma2->size - 1;

  double tr, d, sigma2_i, sigma2_j;

  // Update parameters.
  UpdateParam(log_sigma2, p);

  for (size_t i = 0; i < n_vc + 1; i++) {
    if (i == n_vc) {
      tr = 0;
      for (size_t l = 0; l < n1; l++) {
        tr += gsl_matrix_get(p->P, l, l);
      }
    } else {
      tr = 0;
      for (size_t l = 0; l < n1; l++) {
        gsl_vector_view P_row = gsl_matrix_row(p->P, l);
        gsl_vector_const_view K_col = gsl_matrix_const_column(p->K, n1 * i + l);
        gsl_blas_ddot(&P_row.vector, &K_col.vector, &d);
        tr += d;
      }
    }

    gsl_vector_view KPy_i = gsl_matrix_column(p->KPy_mat, i);
    gsl_blas_ddot(p->Py, &KPy_i.vector, &d);

    if (p->noconstrain) {
      sigma2_i = gsl_vector_get(log_sigma2, i);
      d = (-0.5 * tr + 0.5 * d);
    } else {
      sigma2_i = exp(gsl_vector_get(log_sigma2, i));
      d = (-0.5 * tr + 0.5 * d) * sigma2_i;
    }

    gsl_vector_set(dev1, i, d);

    for (size_t j = i; j < n_vc + 1; j++) {
      gsl_vector_view PKPy_j = gsl_matrix_column(p->PKPy_mat, j);
      gsl_blas_ddot(&KPy_i.vector, &PKPy_j.vector, &d);

      if (p->noconstrain) {
        sigma2_j = gsl_vector_get(log_sigma2, j);
        d *= -0.5;
      } else {
        sigma2_j = exp(gsl_vector_get(log_sigma2, j));
        d *= -0.5 * sigma2_i * sigma2_j;
      }

      gsl_matrix_set(dev2, i, j, d);
      if (j != i) {
        gsl_matrix_set(dev2, j, i, d);
      }
    }
  }

  gsl_matrix_memcpy(p->Hessian, dev2);

  return GSL_SUCCESS;
}

// Read header to determine which column contains which item.
bool ReadHeader_vc(const string &line, HEADER &header) {
  debug_msg("entering");
  string rs_ptr[] = {"rs",   "RS",    "snp",   "SNP",  "snps",
                     "SNPS", "snpid", "SNPID", "rsid", "RSID"};
  set<string> rs_set(rs_ptr, rs_ptr + 10);
  string chr_ptr[] = {"chr", "CHR"};
  set<string> chr_set(chr_ptr, chr_ptr + 2);
  string pos_ptr[] = {
      "ps", "PS", "pos", "POS", "base_position", "BASE_POSITION", "bp", "BP"};
  set<string> pos_set(pos_ptr, pos_ptr + 8);
  string cm_ptr[] = {"cm", "CM"};
  set<string> cm_set(cm_ptr, cm_ptr + 2);
  string a1_ptr[] = {"a1", "A1", "allele1", "ALLELE1"};
  set<string> a1_set(a1_ptr, a1_ptr + 4);
  string a0_ptr[] = {"a0", "A0", "allele0", "ALLELE0"};
  set<string> a0_set(a0_ptr, a0_ptr + 4);

  string z_ptr[] = {"z", "Z", "z_score", "Z_SCORE", "zscore", "ZSCORE"};
  set<string> z_set(z_ptr, z_ptr + 6);
  string beta_ptr[] = {"beta", "BETA", "b", "B"};
  set<string> beta_set(beta_ptr, beta_ptr + 4);
  string sebeta_ptr[] = {"se_beta", "SE_BETA", "se", "SE"};
  set<string> sebeta_set(sebeta_ptr, sebeta_ptr + 4);
  string chisq_ptr[] = {"chisq", "CHISQ", "chisquare", "CHISQUARE"};
  set<string> chisq_set(chisq_ptr, chisq_ptr + 4);
  string p_ptr[] = {"p", "P", "pvalue", "PVALUE", "p-value", "P-VALUE"};
  set<string> p_set(p_ptr, p_ptr + 6);

  string n_ptr[] = {"n", "N", "ntotal", "NTOTAL", "n_total", "N_TOTAL"};
  set<string> n_set(n_ptr, n_ptr + 6);
  string nmis_ptr[] = {"nmis", "NMIS", "n_mis", "N_MIS", "n_miss", "N_MISS"};
  set<string> nmis_set(nmis_ptr, nmis_ptr + 6);
  string nobs_ptr[] = {"nobs", "NOBS", "n_obs", "N_OBS"};
  set<string> nobs_set(nobs_ptr, nobs_ptr + 4);

  string af_ptr[] = {"af",
                     "AF",
                     "maf",
                     "MAF",
                     "f",
                     "F",
                     "allele_freq",
                     "ALLELE_FREQ",
                     "allele_frequency",
                     "ALLELE_FREQUENCY"};
  set<string> af_set(af_ptr, af_ptr + 10);
  string var_ptr[] = {"var", "VAR"};
  set<string> var_set(var_ptr, var_ptr + 2);

  string ws_ptr[] = {"window_size", "WINDOW_SIZE", "ws", "WS"};
  set<string> ws_set(ws_ptr, ws_ptr + 4);
  string cor_ptr[] = {"cor", "COR", "r", "R"};
  set<string> cor_set(cor_ptr, cor_ptr + 4);

  header.rs_col = 0;
  header.chr_col = 0;
  header.pos_col = 0;
  header.a1_col = 0;
  header.a0_col = 0;
  header.z_col = 0;
  header.beta_col = 0;
  header.sebeta_col = 0;
  header.chisq_col = 0;
  header.p_col = 0;
  header.n_col = 0;
  header.nmis_col = 0;
  header.nobs_col = 0;
  header.af_col = 0;
  header.var_col = 0;
  header.ws_col = 0;
  header.cor_col = 0;
  header.coln = 0;

  char *ch_ptr;
  string type;
  size_t n_error = 0;

  ch_ptr = strtok((char *)line.c_str(), " , \t");
  while (ch_ptr != NULL) {
    type = ch_ptr;
    if (rs_set.count(type) != 0) {
      if (header.rs_col == 0) {
        header.rs_col = header.coln + 1;
      } else {
        cout << "error! more than two rs columns in the file." << endl;
        n_error++;
      }
    } else if (chr_set.count(type) != 0) {
      if (header.chr_col == 0) {
        header.chr_col = header.coln + 1;
      } else {
        cout << "error! more than two chr columns in the file." << endl;
        n_error++;
      }
    } else if (pos_set.count(type) != 0) {
      if (header.pos_col == 0) {
        header.pos_col = header.coln + 1;
      } else {
        cout << "error! more than two pos columns in the file." << endl;
        n_error++;
      }
    } else if (cm_set.count(type) != 0) {
      if (header.cm_col == 0) {
        header.cm_col = header.coln + 1;
      } else {
        cout << "error! more than two cm columns in the file." << endl;
        n_error++;
      }
    } else if (a1_set.count(type) != 0) {
      if (header.a1_col == 0) {
        header.a1_col = header.coln + 1;
      } else {
        cout << "error! more than two allele1 columns in the file." << endl;
        n_error++;
      }
    } else if (a0_set.count(type) != 0) {
      if (header.a0_col == 0) {
        header.a0_col = header.coln + 1;
      } else {
        cout << "error! more than two allele0 columns in the file." << endl;
        n_error++;
      }
    } else if (z_set.count(type) != 0) {
      if (header.z_col == 0) {
        header.z_col = header.coln + 1;
      } else {
        cout << "error! more than two z columns in the file." << endl;
        n_error++;
      }
    } else if (beta_set.count(type) != 0) {
      if (header.beta_col == 0) {
        header.beta_col = header.coln + 1;
      } else {
        cout << "error! more than two beta columns in the file." << endl;
        n_error++;
      }
    } else if (sebeta_set.count(type) != 0) {
      if (header.sebeta_col == 0) {
        header.sebeta_col = header.coln + 1;
      } else {
        cout << "error! more than two se_beta columns in the file." << endl;
        n_error++;
      }
    } else if (chisq_set.count(type) != 0) {
      if (header.chisq_col == 0) {
        header.chisq_col = header.coln + 1;
      } else {
        cout << "error! more than two z columns in the file." << endl;
        n_error++;
      }
    } else if (p_set.count(type) != 0) {
      if (header.p_col == 0) {
        header.p_col = header.coln + 1;
      } else {
        cout << "error! more than two p columns in the file." << endl;
        n_error++;
      }
    } else if (n_set.count(type) != 0) {
      if (header.n_col == 0) {
        header.n_col = header.coln + 1;
      } else {
        cout << "error! more than two n_total columns in the file." << endl;
        n_error++;
      }
    } else if (nmis_set.count(type) != 0) {
      if (header.nmis_col == 0) {
        header.nmis_col = header.coln + 1;
      } else {
        cout << "error! more than two n_mis columns in the file." << endl;
        n_error++;
      }
    } else if (nobs_set.count(type) != 0) {
      if (header.nobs_col == 0) {
        header.nobs_col = header.coln + 1;
      } else {
        cout << "error! more than two n_obs columns in the file." << endl;
        n_error++;
      }
    } else if (ws_set.count(type) != 0) {
      if (header.ws_col == 0) {
        header.ws_col = header.coln + 1;
      } else {
        cout << "error! more than two window_size columns in the file." << endl;
        n_error++;
      }
    } else if (af_set.count(type) != 0) {
      if (header.af_col == 0) {
        header.af_col = header.coln + 1;
      } else {
        cout << "error! more than two af columns in the file." << endl;
        n_error++;
      }
    } else if (cor_set.count(type) != 0) {
      if (header.cor_col == 0) {
        header.cor_col = header.coln + 1;
      } else {
        cout << "error! more than two cor columns in the file." << endl;
        n_error++;
      }
    } else {
    }

    ch_ptr = strtok(NULL, " , \t");
    header.coln++;
  }

  if (header.cor_col != 0 && header.cor_col != header.coln) {
    cout << "error! the cor column should be the last column." << endl;
    n_error++;
  }

  if (header.rs_col == 0) {
    if (header.chr_col != 0 && header.pos_col != 0) {
      cout << "missing an rs column. rs id will be replaced by chr:pos" << endl;
    } else {
      cout << "error! missing an rs column." << endl;
      n_error++;
    }
  }

  if (n_error == 0) {
    return true;
  } else {
    return false;
  }
}

// Read cov file the first time, record mapRS2in, mapRS2var (in case
// var is not provided in the z file), store vec_n and vec_rs.
void ReadFile_cor(const string &file_cor, const set<string> &setSnps,
                  vector<string> &vec_rs, vector<size_t> &vec_n,
                  vector<double> &vec_cm, vector<double> &vec_bp,
                  map<string, size_t> &mapRS2in,
                  map<string, double> &mapRS2var) {
  debug_msg("entering");
  vec_rs.clear();
  vec_n.clear();
  mapRS2in.clear();
  mapRS2var.clear();

  igzstream infile(file_cor.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open cov file: " << file_cor << endl;
    return;
  }

  string line;
  char *ch_ptr;

  string rs, chr, a1, a0, pos, cm;
  double af = 0, var_x = 0, d_pos, d_cm;
  size_t n_total = 0, n_mis = 0, n_obs = 0, ni_total = 0;
  size_t ns_test = 0, ns_total = 0;

  HEADER header;

  // Header.
  safeGetline(infile, line).eof();
  ReadHeader_vc(line, header);

  if (header.n_col == 0) {
    if (header.nobs_col == 0 && header.nmis_col == 0) {
      cout << "error! missing sample size in the cor file." << endl;
    } else {
      cout << "total sample size will be replaced by obs/mis sample size."
           << endl;
    }
  }

  while (!safeGetline(infile, line).eof()) {

    // do not read cor values this time; upto col_n-1.
    ch_ptr = strtok_safe((char *)line.c_str(), " , \t");

    n_total = 0;
    n_mis = 0;
    n_obs = 0;
    af = 0;
    var_x = 0;
    d_cm = 0;
    d_pos = 0;
    for (size_t i = 0; i < header.coln - 1; i++) {
      enforce(ch_ptr);
      if (header.rs_col != 0 && header.rs_col == i + 1) {
        rs = ch_ptr;
      }
      if (header.chr_col != 0 && header.chr_col == i + 1) {
        chr = ch_ptr;
      }
      if (header.pos_col != 0 && header.pos_col == i + 1) {
        pos = ch_ptr;
        d_pos = atof(ch_ptr);
      }
      if (header.cm_col != 0 && header.cm_col == i + 1) {
        cm = ch_ptr;
        d_cm = atof(ch_ptr);
      }
      if (header.a1_col != 0 && header.a1_col == i + 1) {
        a1 = ch_ptr;
      }
      if (header.a0_col != 0 && header.a0_col == i + 1) {
        a0 = ch_ptr;
      }

      if (header.n_col != 0 && header.n_col == i + 1) {
        n_total = atoi(ch_ptr);
      }
      if (header.nmis_col != 0 && header.nmis_col == i + 1) {
        n_mis = atoi(ch_ptr);
      }
      if (header.nobs_col != 0 && header.nobs_col == i + 1) {
        n_obs = atoi(ch_ptr);
      }

      if (header.af_col != 0 && header.af_col == i + 1) {
        af = atof(ch_ptr);
      }
      if (header.var_col != 0 && header.var_col == i + 1) {
        var_x = atof(ch_ptr);
      }

      ch_ptr = strtok(NULL, " , \t");
    }

    if (header.rs_col == 0) {
      rs = chr + ":" + pos;
    }

    if (header.n_col == 0) {
      n_total = n_mis + n_obs;
    }

    // Record rs, n.
    vec_rs.push_back(rs);
    vec_n.push_back(n_total);
    if (d_cm > 0) {
      vec_cm.push_back(d_cm);
    } else {
      vec_cm.push_back(d_cm);
    }
    if (d_pos > 0) {
      vec_bp.push_back(d_pos);
    } else {
      vec_bp.push_back(d_pos);
    }

    // Record mapRS2in and mapRS2var.
    if (setSnps.size() == 0 || setSnps.count(rs) != 0) {
      if (mapRS2in.count(rs) == 0) {
        mapRS2in[rs] = 1;

        if (header.var_col != 0) {
          mapRS2var[rs] = var_x;
        } else if (header.af_col != 0) {
          var_x = 2.0 * af * (1.0 - af);
          mapRS2var[rs] = var_x;
        } else {
        }

        ns_test++;

      } else {
        cout << "error! more than one snp has the same id " << rs
             << " in cor file?" << endl;
      }
    }

    // Record max pos.
    ni_total = max(ni_total, n_total);
    ns_total++;
  }

  infile.close();
  infile.clear();

  return;
}

// Read beta file, store mapRS2var if var is provided here, calculate
// q and var_y.
void ReadFile_beta(const bool flag_priorscale, const string &file_beta,
                   const map<string, size_t> &mapRS2cat,
                   map<string, size_t> &mapRS2in,
                   map<string, double> &mapRS2var,
                   map<string, size_t> &mapRS2nsamp, gsl_vector *q_vec,
                   gsl_vector *qvar_vec, gsl_vector *s_vec, size_t &ni_total,
                   size_t &ns_total) {
  debug_msg("entering");
  mapRS2nsamp.clear();

  igzstream infile(file_beta.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open beta file: " << file_beta << endl;
    return;
  }

  string line;
  char *ch_ptr;
  string type;

  string rs, chr, a1, a0, pos, cm;
  double z = 0, beta = 0, se_beta = 0, chisq = 0, pvalue = 0, zsquare = 0,
         af = 0, var_x = 0;
  size_t n_total = 0, n_mis = 0, n_obs = 0;
  size_t ns_test = 0;
  ns_total = 0;
  ni_total = 0;

  vector<double> vec_q, vec_qvar, vec_s;
  for (size_t i = 0; i < q_vec->size; i++) {
    vec_q.push_back(0.0);
    vec_qvar.push_back(0.0);
    vec_s.push_back(0.0);
  }

  // Read header.
  HEADER header;
  safeGetline(infile, line).eof();
  ReadHeader_vc(line, header);

  if (header.n_col == 0) {
    if (header.nobs_col == 0 && header.nmis_col == 0) {
      cout << "error! missing sample size in the beta file." << endl;
    } else {
      cout << "total sample size will be replaced by obs/mis sample size."
           << endl;
    }
  }

  if (header.z_col == 0 && (header.beta_col == 0 || header.sebeta_col == 0) &&
      header.chisq_col == 0 && header.p_col == 0) {
    cout << "error! missing z scores in the beta file." << endl;
  }

  if (header.af_col == 0 && header.var_col == 0 && mapRS2var.size() == 0) {
    cout << "error! missing allele frequency in the beta file." << endl;
  }

  while (!safeGetline(infile, line).eof()) {
    ch_ptr = strtok_safe((char *)line.c_str(), " , \t");

    z = 0;
    beta = 0;
    se_beta = 0;
    chisq = 0;
    pvalue = 0;
    n_total = 0;
    n_mis = 0;
    n_obs = 0;
    af = 0;
    var_x = 0;
    for (size_t i = 0; i < header.coln; i++) {
      enforce(ch_ptr);
      if (header.rs_col != 0 && header.rs_col == i + 1) {
        rs = ch_ptr;
      }
      if (header.chr_col != 0 && header.chr_col == i + 1) {
        chr = ch_ptr;
      }
      if (header.pos_col != 0 && header.pos_col == i + 1) {
        pos = ch_ptr;
      }
      if (header.cm_col != 0 && header.cm_col == i + 1) {
        cm = ch_ptr;
      }
      if (header.a1_col != 0 && header.a1_col == i + 1) {
        a1 = ch_ptr;
      }
      if (header.a0_col != 0 && header.a0_col == i + 1) {
        a0 = ch_ptr;
      }

      if (header.z_col != 0 && header.z_col == i + 1) {
        z = atof(ch_ptr);
      }
      if (header.beta_col != 0 && header.beta_col == i + 1) {
        beta = atof(ch_ptr);
      }
      if (header.sebeta_col != 0 && header.sebeta_col == i + 1) {
        se_beta = atof(ch_ptr);
      }
      if (header.chisq_col != 0 && header.chisq_col == i + 1) {
        chisq = atof(ch_ptr);
      }
      if (header.p_col != 0 && header.p_col == i + 1) {
        pvalue = atof(ch_ptr);
      }

      if (header.n_col != 0 && header.n_col == i + 1) {
        n_total = atoi(ch_ptr);
      }
      if (header.nmis_col != 0 && header.nmis_col == i + 1) {
        n_mis = atoi(ch_ptr);
      }
      if (header.nobs_col != 0 && header.nobs_col == i + 1) {
        n_obs = atoi(ch_ptr);
      }

      if (header.af_col != 0 && header.af_col == i + 1) {
        af = atof(ch_ptr);
      }
      if (header.var_col != 0 && header.var_col == i + 1) {
        var_x = atof(ch_ptr);
      }

      ch_ptr = strtok(NULL, " , \t");
    }

    if (header.rs_col == 0) {
      rs = chr + ":" + pos;
    }

    if (header.n_col == 0) {
      n_total = n_mis + n_obs;
    }

    // Both z values and beta/se_beta have directions, while
    // chisq/pvalue do not.
    if (header.z_col != 0) {
      zsquare = z * z;
    } else if (header.beta_col != 0 && header.sebeta_col != 0) {
      z = beta / se_beta;
      zsquare = z * z;
    } else if (header.chisq_col != 0) {
      zsquare = chisq;
    } else if (header.p_col != 0) {
      zsquare = gsl_cdf_chisq_Qinv(pvalue, 1);
    } else {
      zsquare = 0;
    }

    // If the snp is also present in cor file, then do calculations.
    if ((header.var_col != 0 || header.af_col != 0 ||
         mapRS2var.count(rs) != 0) &&
        mapRS2in.count(rs) != 0 &&
        (mapRS2cat.size() == 0 || mapRS2cat.count(rs) != 0)) {
      if (mapRS2in.at(rs) > 1) {
        cout << "error! more than one snp has the same id " << rs
             << " in beta file?" << endl;
        break;
      }

      if (header.var_col == 0) {
        if (header.af_col != 0) {
          var_x = 2.0 * af * (1.0 - af);
        } else {
          var_x = mapRS2var.at(rs);
        }
      }

      if (flag_priorscale) {
        var_x = 1;
      }

      mapRS2in[rs]++;
      mapRS2var[rs] = var_x;
      mapRS2nsamp[rs] = n_total;

      if (mapRS2cat.size() != 0) {
        vec_q[mapRS2cat.at(rs)] += (zsquare - 1.0) * var_x / (double)n_total;
        vec_s[mapRS2cat.at(rs)] += var_x;
        vec_qvar[mapRS2cat.at(rs)] +=
            var_x * var_x / ((double)n_total * (double)n_total);
      } else {
        vec_q[0] += (zsquare - 1.0) * var_x / (double)n_total;
        vec_s[0] += var_x;
        vec_qvar[0] += var_x * var_x / ((double)n_total * (double)n_total);
      }

      ni_total = max(ni_total, n_total);
      ns_test++;
    }

    ns_total++;
  }

  for (size_t i = 0; i < q_vec->size; i++) {
    gsl_vector_set(q_vec, i, vec_q[i]);
    gsl_vector_set(qvar_vec, i, 2.0 * vec_qvar[i]);
    gsl_vector_set(s_vec, i, vec_s[i]);
  }

  infile.clear();
  infile.close();

  return;
}

// Read covariance file the second time.
// Look for rs, n_mis+n_obs, var, window_size, cov.
// If window_cm/bp/ns is provided, then use these max values to
// calibrate estimates.
void ReadFile_cor(const string &file_cor, const vector<string> &vec_rs,
                  const vector<size_t> &vec_n, const vector<double> &vec_cm,
                  const vector<double> &vec_bp,
                  const map<string, size_t> &mapRS2cat,
                  const map<string, size_t> &mapRS2in,
                  const map<string, double> &mapRS2var,
                  const map<string, size_t> &mapRS2nsamp, const size_t crt,
                  const double &window_cm, const double &window_bp,
                  const double &window_ns, gsl_matrix *S_mat,
                  gsl_matrix *Svar_mat, gsl_vector *qvar_vec, size_t &ni_total,
                  size_t &ns_total, size_t &ns_test, size_t &ns_pair) {
  debug_msg("entering");
  igzstream infile(file_cor.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open cov file: " << file_cor << endl;
    return;
  }

  string line;
  char *ch_ptr;

  string rs1, rs2;
  double d1, d2, d3, cor, var1, var2;
  size_t n_nb, nsamp1, nsamp2, n12, bin_size = 10, bin;

  vector<vector<double>> mat_S, mat_Svar, mat_tmp;
  vector<double> vec_qvar, vec_tmp;
  vector<vector<vector<double>>> mat3d_Sbin;

  for (size_t i = 0; i < S_mat->size1; i++) {
    vec_qvar.push_back(0.0);
  }

  for (size_t i = 0; i < S_mat->size1; i++) {
    mat_S.push_back(vec_qvar);
    mat_Svar.push_back(vec_qvar);
  }

  for (size_t k = 0; k < bin_size; k++) {
    vec_tmp.push_back(0.0);
  }
  for (size_t i = 0; i < S_mat->size1; i++) {
    mat_tmp.push_back(vec_tmp);
  }
  for (size_t i = 0; i < S_mat->size1; i++) {
    mat3d_Sbin.push_back(mat_tmp);
  }

  string rs, chr, a1, a0, type, pos, cm;
  size_t n_total = 0, n_mis = 0, n_obs = 0;
  double d_pos1, d_pos2, d_pos, d_cm1, d_cm2, d_cm;
  ns_test = 0;
  ns_total = 0;
  ns_pair = 0;
  ni_total = 0;

  // Header.
  HEADER header;

  safeGetline(infile, line).eof();
  ReadHeader_vc(line, header);

  while (!safeGetline(infile, line).eof()) {

    // Do not read cor values this time; upto col_n-1.
    d_pos1 = 0;
    d_cm1 = 0;
    ch_ptr = strtok_safe((char *)line.c_str(), " , \t");
    for (size_t i = 0; i < header.coln - 1; i++) {
      enforce(ch_ptr);
      if (header.rs_col != 0 && header.rs_col == i + 1) {
        rs = ch_ptr;
      }
      if (header.chr_col != 0 && header.chr_col == i + 1) {
        chr = ch_ptr;
      }
      if (header.pos_col != 0 && header.pos_col == i + 1) {
        pos = ch_ptr;
        d_pos1 = atof(ch_ptr);
      }
      if (header.cm_col != 0 && header.cm_col == i + 1) {
        cm = ch_ptr;
        d_cm1 = atof(ch_ptr);
      }
      if (header.a1_col != 0 && header.a1_col == i + 1) {
        a1 = ch_ptr;
      }
      if (header.a0_col != 0 && header.a0_col == i + 1) {
        a0 = ch_ptr;
      }

      if (header.n_col != 0 && header.n_col == i + 1) {
        n_total = atoi(ch_ptr);
      }
      if (header.nmis_col != 0 && header.nmis_col == i + 1) {
        n_mis = atoi(ch_ptr);
      }
      if (header.nobs_col != 0 && header.nobs_col == i + 1) {
        n_obs = atoi(ch_ptr);
      }

      ch_ptr = strtok(NULL, " , \t");
    }

    if (header.rs_col == 0) {
      rs = chr + ":" + pos;
    }

    if (header.n_col == 0) {
      n_total = n_mis + n_obs;
    }

    rs1 = rs;

    if ((mapRS2cat.size() == 0 || mapRS2cat.count(rs1) != 0) &&
        mapRS2in.count(rs1) != 0 && mapRS2in.at(rs1) == 2) {
      var1 = mapRS2var.at(rs1);
      nsamp1 = mapRS2nsamp.at(rs1);
      d2 = var1 * var1;

      if (mapRS2cat.size() != 0) {
        mat_S[mapRS2cat.at(rs1)][mapRS2cat.at(rs1)] +=
            (1 - 1.0 / (double)vec_n[ns_total]) * d2;
        mat_Svar[mapRS2cat.at(rs1)][mapRS2cat.at(rs1)] +=
            d2 * d2 / ((double)vec_n[ns_total] * (double)vec_n[ns_total]);
        if (crt == 1) {
          mat3d_Sbin[mapRS2cat.at(rs1)][mapRS2cat.at(rs1)][0] +=
              (1 - 1.0 / (double)vec_n[ns_total]) * d2;
        }
      } else {
        mat_S[0][0] += (1 - 1.0 / (double)vec_n[ns_total]) * d2;
        mat_Svar[0][0] +=
            d2 * d2 / ((double)vec_n[ns_total] * (double)vec_n[ns_total]);
        if (crt == 1) {
          mat3d_Sbin[0][0][0] += (1 - 1.0 / (double)vec_n[ns_total]) * d2;
        }
      }

      n_nb = 0;
      while (ch_ptr != NULL) {
        type = ch_ptr;
        if (type.compare("NA") != 0 && type.compare("na") != 0 &&
            type.compare("nan") != 0 && type.compare("-nan") != 0) {
          cor = atof(ch_ptr);
          rs2 = vec_rs[ns_total + n_nb + 1];
          d_pos2 = vec_bp[ns_total + n_nb + 1];
          d_cm2 = vec_cm[ns_total + n_nb + 1];
          d_pos = abs(d_pos2 - d_pos1);
          d_cm = abs(d_cm2 - d_cm1);

          if ((mapRS2cat.size() == 0 || mapRS2cat.count(rs2) != 0) &&
              mapRS2in.count(rs2) != 0 && mapRS2in.at(rs2) == 2) {
            var2 = mapRS2var.at(rs2);
            nsamp2 = mapRS2nsamp.at(rs2);
            d1 = cor * cor -
                 1.0 / (double)min(vec_n[ns_total], vec_n[ns_total + n_nb + 1]);
            d2 = var1 * var2;
            d3 = cor * cor / ((double)nsamp1 * (double)nsamp2);
            n12 = min(vec_n[ns_total], vec_n[ns_total + n_nb + 1]);

            // Compute bin.
            if (crt == 1) {
              if (window_cm != 0 && d_cm1 != 0 && d_cm2 != 0) {
                bin =
                    min((int)floor(d_cm / window_cm * bin_size), (int)bin_size);
              } else if (window_bp != 0 && d_pos1 != 0 && d_pos2 != 0) {
                bin = min((int)floor(d_pos / window_bp * bin_size),
                          (int)bin_size);
              } else if (window_ns != 0) {
                bin = min((int)floor(((double)n_nb + 1) / window_ns * bin_size),
                          (int)bin_size);
              }
            }

            if (mapRS2cat.size() != 0) {
              if (mapRS2cat.at(rs1) == mapRS2cat.at(rs2)) {
                vec_qvar[mapRS2cat.at(rs1)] += 2 * d3 * d2;
                mat_S[mapRS2cat.at(rs1)][mapRS2cat.at(rs2)] += 2 * d1 * d2;
                mat_Svar[mapRS2cat.at(rs1)][mapRS2cat.at(rs2)] +=
                    2 * d2 * d2 / ((double)n12 * (double)n12);
                if (crt == 1) {
                  mat3d_Sbin[mapRS2cat.at(rs1)][mapRS2cat.at(rs2)][bin] +=
                      2 * d1 * d2;
                }
              } else {
                mat_S[mapRS2cat.at(rs1)][mapRS2cat.at(rs2)] += d1 * d2;
                mat_Svar[mapRS2cat.at(rs1)][mapRS2cat.at(rs2)] +=
                    d2 * d2 / ((double)n12 * (double)n12);
                if (crt == 1) {
                  mat3d_Sbin[mapRS2cat.at(rs1)][mapRS2cat.at(rs2)][bin] +=
                      d1 * d2;
                }
              }
            } else {
              vec_qvar[0] += 2 * d3 * d2;
              mat_S[0][0] += 2 * d1 * d2;
              mat_Svar[0][0] += 2 * d2 * d2 / ((double)n12 * (double)n12);

              if (crt == 1) {
                mat3d_Sbin[0][0][bin] += 2 * d1 * d2;
              }
            }
            ns_pair++;
          }
        }

        ch_ptr = strtok(NULL, " , \t");
        n_nb++;
      }
      ni_total = max(ni_total, n_total);
      ns_test++;
    }

    ns_total++;
  }

  // Use S_bin to fit a rational function y=1/(a+bx)^2, where
  // x=seq(0.5,bin_size-0.5,by=1) and then compute a correlation
  // factor as a percentage.
  double a, b, x, y, n, var_y, var_x, mean_y, mean_x, cov_xy, crt_factor;
  if (crt == 1) {
    for (size_t i = 0; i < S_mat->size1; i++) {
      for (size_t j = i; j < S_mat->size2; j++) {

        // Correct mat_S.
        n = 0;
        var_y = 0;
        var_x = 0;
        mean_y = 0;
        mean_x = 0;
        cov_xy = 0;
        for (size_t k = 0; k < bin_size; k++) {
          if (j == i) {
            y = mat3d_Sbin[i][j][k];
          } else {
            y = mat3d_Sbin[i][j][k] + mat3d_Sbin[j][i][k];
          }
          x = k + 0.5;
          cout << y << ", ";
          if (y > 0) {
            y = 1 / sqrt(y);
            mean_x += x;
            mean_y += y;
            var_x += x * x;
            var_y += y * y;
            cov_xy += x * y;
            n++;
          }
        }
        cout << endl;

        if (n >= 5) {
          mean_x /= n;
          mean_y /= n;
          var_x /= n;
          var_y /= n;
          cov_xy /= n;
          var_x -= mean_x * mean_x;
          var_y -= mean_y * mean_y;
          cov_xy -= mean_x * mean_y;
          b = cov_xy / var_x;
          a = mean_y - b * mean_x;
          crt_factor = a / (b * (bin_size + 0.5)) + 1;
          if (i == j) {
            mat_S[i][j] *= crt_factor;
          } else {
            mat_S[i][j] *= crt_factor;
            mat_S[j][i] *= crt_factor;
          }
          cout << crt_factor << endl;

          // Correct qvar.
          if (i == j) {
            vec_qvar[i] *= crt_factor;
          }
        }
      }
    }
  }

  // Save to gsl_vector and gsl_matrix: qvar_vec, S_mat, Svar_mat.
  for (size_t i = 0; i < S_mat->size1; i++) {
    d1 = gsl_vector_get(qvar_vec, i) + 2 * vec_qvar[i];
    gsl_vector_set(qvar_vec, i, d1);
    for (size_t j = 0; j < S_mat->size2; j++) {
      if (i == j) {
        gsl_matrix_set(S_mat, i, j, mat_S[i][i]);
        gsl_matrix_set(Svar_mat, i, j, 2.0 * mat_Svar[i][i] * ns_test *
                                           ns_test / (2.0 * ns_pair));
      } else {
        gsl_matrix_set(S_mat, i, j, mat_S[i][j] + mat_S[j][i]);
        gsl_matrix_set(Svar_mat, i, j, 2.0 * (mat_Svar[i][j] + mat_Svar[j][i]) *
                                           ns_test * ns_test / (2.0 * ns_pair));
      }
    }
  }

  infile.clear();
  infile.close();

  return;
}

// Use the new method to calculate variance components with summary
// statistics first, use a function CalcS to compute S matrix (where
// the diagonal elements are part of V(q) ), and then use bootstrap to
// compute the variance for S, use a set of genotypes, phenotypes, and
// individual ids, and snp category label.
void CalcVCss(const gsl_matrix *Vq, const gsl_matrix *S_mat,
              const gsl_matrix *Svar_mat, const gsl_vector *q_vec,
              const gsl_vector *s_vec, const double df, vector<double> &v_pve,
              vector<double> &v_se_pve, double &pve_total, double &se_pve_total,
              vector<double> &v_sigma2, vector<double> &v_se_sigma2,
              vector<double> &v_enrich, vector<double> &v_se_enrich) {
  size_t n_vc = S_mat->size1;

  gsl_matrix *Si_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *Var_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *tmp_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *tmp_mat1 = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *VarEnrich_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *qvar_mat = gsl_matrix_alloc(n_vc, n_vc);

  gsl_vector *pve = gsl_vector_alloc(n_vc);
  gsl_vector *pve_plus = gsl_vector_alloc(n_vc + 1);
  gsl_vector *tmp = gsl_vector_alloc(n_vc + 1);
  gsl_vector *sigma2persnp = gsl_vector_alloc(n_vc);
  gsl_vector *enrich = gsl_vector_alloc(n_vc);
  gsl_vector *se_pve = gsl_vector_alloc(n_vc);
  gsl_vector *se_sigma2persnp = gsl_vector_alloc(n_vc);
  gsl_vector *se_enrich = gsl_vector_alloc(n_vc);

  double d;

  // Calculate S^{-1}q.
  gsl_matrix_memcpy(tmp_mat, S_mat);
  int sig;
  gsl_permutation *pmt = gsl_permutation_alloc(n_vc);
  LUDecomp(tmp_mat, pmt, &sig);
  LUInvert(tmp_mat, pmt, Si_mat);

  // Calculate sigma2snp and pve.
  gsl_blas_dgemv(CblasNoTrans, 1.0, Si_mat, q_vec, 0.0, pve);
  gsl_vector_memcpy(sigma2persnp, pve);
  gsl_vector_div(sigma2persnp, s_vec);

  // Get qvar_mat.
  gsl_matrix_memcpy(qvar_mat, Vq);
  gsl_matrix_scale(qvar_mat, 1.0 / (df * df));

  // Calculate variance for these estimates.
  for (size_t i = 0; i < n_vc; i++) {
    for (size_t j = i; j < n_vc; j++) {
      d = gsl_matrix_get(Svar_mat, i, j);
      d *= gsl_vector_get(pve, i) * gsl_vector_get(pve, j);

      d += gsl_matrix_get(qvar_mat, i, j);
      gsl_matrix_set(Var_mat, i, j, d);
      if (i != j) {
        gsl_matrix_set(Var_mat, j, i, d);
      }
    }
  }

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Si_mat, Var_mat, 0.0,
                 tmp_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp_mat, Si_mat, 0.0,
                 Var_mat);

  for (size_t i = 0; i < n_vc; i++) {
    d = sqrt(gsl_matrix_get(Var_mat, i, i));
    gsl_vector_set(se_pve, i, d);
    d /= gsl_vector_get(s_vec, i);
    gsl_vector_set(se_sigma2persnp, i, d);
  }

  // Compute pve_total, se_pve_total.
  pve_total = 0;
  se_pve_total = 0;
  for (size_t i = 0; i < n_vc; i++) {
    pve_total += gsl_vector_get(pve, i);

    for (size_t j = 0; j < n_vc; j++) {
      se_pve_total += gsl_matrix_get(Var_mat, i, j);
    }
  }
  se_pve_total = sqrt(se_pve_total);

  // Compute enrichment and its variance.
  double s_pve = 0, s_snp = 0;
  for (size_t i = 0; i < n_vc; i++) {
    s_pve += gsl_vector_get(pve, i);
    s_snp += gsl_vector_get(s_vec, i);
  }
  gsl_vector_memcpy(enrich, sigma2persnp);
  gsl_vector_scale(enrich, s_snp / s_pve);

  gsl_matrix_set_identity(tmp_mat);

  double d1;
  for (size_t i = 0; i < n_vc; i++) {
    d = gsl_vector_get(pve, i) / s_pve;
    d1 = gsl_vector_get(s_vec, i);
    for (size_t j = 0; j < n_vc; j++) {
      if (i == j) {
        gsl_matrix_set(tmp_mat, i, j, (1 - d) / d1 * s_snp / s_pve);
      } else {
        gsl_matrix_set(tmp_mat, i, j, -1 * d / d1 * s_snp / s_pve);
      }
    }
  }
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp_mat, Var_mat, 0.0,
                 tmp_mat1);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmp_mat1, tmp_mat, 0.0,
                 VarEnrich_mat);

  for (size_t i = 0; i < n_vc; i++) {
    d = sqrt(gsl_matrix_get(VarEnrich_mat, i, i));
    gsl_vector_set(se_enrich, i, d);
  }

  cout << "pve = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << gsl_vector_get(pve, i) << " ";
  }
  cout << endl;

  cout << "se(pve) = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << gsl_vector_get(se_pve, i) << " ";
  }
  cout << endl;

  cout << "sigma2 per snp = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << gsl_vector_get(sigma2persnp, i) << " ";
  }
  cout << endl;

  cout << "se(sigma2 per snp) = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << gsl_vector_get(se_sigma2persnp, i) << " ";
  }
  cout << endl;

  cout << "enrichment = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << gsl_vector_get(enrich, i) << " ";
  }
  cout << endl;

  cout << "se(enrichment) = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << gsl_vector_get(se_enrich, i) << " ";
  }
  cout << endl;

  // Save data.
  v_pve.clear();
  v_se_pve.clear();
  v_sigma2.clear();
  v_se_sigma2.clear();
  v_enrich.clear();
  v_se_enrich.clear();
  for (size_t i = 0; i < n_vc; i++) {
    d = gsl_vector_get(pve, i);
    v_pve.push_back(d);
    d = gsl_vector_get(se_pve, i);
    v_se_pve.push_back(d);

    d = gsl_vector_get(sigma2persnp, i);
    v_sigma2.push_back(d);
    d = gsl_vector_get(se_sigma2persnp, i);
    v_se_sigma2.push_back(d);

    d = gsl_vector_get(enrich, i);
    v_enrich.push_back(d);
    d = gsl_vector_get(se_enrich, i);
    v_se_enrich.push_back(d);
  }

  // Delete matrices.
  gsl_matrix_free(Si_mat);
  gsl_matrix_free(Var_mat);
  gsl_matrix_free(VarEnrich_mat);
  gsl_matrix_free(tmp_mat);
  gsl_matrix_free(tmp_mat1);
  gsl_matrix_free(qvar_mat);

  gsl_vector_free(pve);
  gsl_vector_free(pve_plus);
  gsl_vector_free(tmp);
  gsl_vector_free(sigma2persnp);
  gsl_vector_free(enrich);
  gsl_vector_free(se_pve);
  gsl_vector_free(se_sigma2persnp);
  gsl_vector_free(se_enrich);

  return;
}

// Ks are not scaled.
void VC::CalcVChe(const gsl_matrix *K, const gsl_matrix *W,
                  const gsl_vector *y) {
  size_t n1 = K->size1, n2 = K->size2;
  size_t n_vc = n2 / n1;

  double r = (double)n1 / (double)(n1 - W->size2);
  double var_y, var_y_new;
  double d, tr, s, v;
  vector<double> traceG_new;

  // New matrices/vectors.
  gsl_matrix *K_scale = gsl_matrix_alloc(n1, n2);
  gsl_vector *y_scale = gsl_vector_alloc(n1);
  gsl_matrix *Kry = gsl_matrix_alloc(n1, n_vc);
  gsl_matrix *yKrKKry = gsl_matrix_alloc(n_vc, n_vc * (n_vc + 1));
  gsl_vector *KKry = gsl_vector_alloc(n1);

  // Old matrices/vectors.
  gsl_vector *pve = gsl_vector_alloc(n_vc);
  gsl_vector *se_pve = gsl_vector_alloc(n_vc);
  gsl_vector *q_vec = gsl_vector_alloc(n_vc);
  gsl_matrix *qvar_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *tmp_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *S_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *Si_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *Var_mat = gsl_matrix_alloc(n_vc, n_vc);

  // Center and scale K by W.
  for (size_t i = 0; i < n_vc; i++) {
    gsl_matrix_view Kscale_sub =
        gsl_matrix_submatrix(K_scale, 0, n1 * i, n1, n1);
    gsl_matrix_const_view K_sub =
        gsl_matrix_const_submatrix(K, 0, n1 * i, n1, n1);
    gsl_matrix_memcpy(&Kscale_sub.matrix, &K_sub.matrix);

    CenterMatrix(&Kscale_sub.matrix, W);
    d = ScaleMatrix(&Kscale_sub.matrix);
    traceG_new.push_back(d);
  }

  // Center y by W, and standardize it to have variance 1 (t(y)%*%y/n=1).
  gsl_vector_memcpy(y_scale, y);
  CenterVector(y_scale, W);

  var_y = VectorVar(y);
  var_y_new = VectorVar(y_scale);

  StandardizeVector(y_scale);

  // Compute Kry, which is used for confidence interval; also compute
  // q_vec (*n^2).
  for (size_t i = 0; i < n_vc; i++) {
    gsl_matrix_const_view Kscale_sub =
        gsl_matrix_const_submatrix(K_scale, 0, n1 * i, n1, n1);
    gsl_vector_view Kry_col = gsl_matrix_column(Kry, i);

    gsl_vector_memcpy(&Kry_col.vector, y_scale);
    gsl_blas_dgemv(CblasNoTrans, 1.0, &Kscale_sub.matrix, y_scale, -1.0 * r,
                   &Kry_col.vector);

    gsl_blas_ddot(&Kry_col.vector, y_scale, &d);
    gsl_vector_set(q_vec, i, d);
  }

  // Compute yKrKKry, which is used later for confidence interval.
  for (size_t i = 0; i < n_vc; i++) {
    gsl_vector_const_view Kry_coli = gsl_matrix_const_column(Kry, i);
    for (size_t j = i; j < n_vc; j++) {
      gsl_vector_const_view Kry_colj = gsl_matrix_const_column(Kry, j);
      for (size_t l = 0; l < n_vc; l++) {
        gsl_matrix_const_view Kscale_sub =
            gsl_matrix_const_submatrix(K_scale, 0, n1 * l, n1, n1);
        gsl_blas_dgemv(CblasNoTrans, 1.0, &Kscale_sub.matrix, &Kry_coli.vector,
                       0.0, KKry);
        gsl_blas_ddot(&Kry_colj.vector, KKry, &d);
        gsl_matrix_set(yKrKKry, i, l * n_vc + j, d);
        if (i != j) {
          gsl_matrix_set(yKrKKry, j, l * n_vc + i, d);
        }
      }
      gsl_blas_ddot(&Kry_coli.vector, &Kry_colj.vector, &d);
      gsl_matrix_set(yKrKKry, i, n_vc * n_vc + j, d);
      if (i != j) {
        gsl_matrix_set(yKrKKry, j, n_vc * n_vc + i, d);
      }
    }
  }

  // Compute Sij (*n^2).
  for (size_t i = 0; i < n_vc; i++) {
    for (size_t j = i; j < n_vc; j++) {
      tr = 0;
      for (size_t l = 0; l < n1; l++) {
        gsl_vector_const_view Ki_col =
            gsl_matrix_const_column(K_scale, i * n1 + l);
        gsl_vector_const_view Kj_col =
            gsl_matrix_const_column(K_scale, j * n1 + l);
        gsl_blas_ddot(&Ki_col.vector, &Kj_col.vector, &d);
        tr += d;
      }

      tr = tr - r * (double)n1;
      gsl_matrix_set(S_mat, i, j, tr);
      if (i != j) {
        gsl_matrix_set(S_mat, j, i, tr);
      }
    }
  }

  // Compute S^{-1}q.
  int sig;
  gsl_permutation *pmt = gsl_permutation_alloc(n_vc);
  LUDecomp(S_mat, pmt, &sig);
  LUInvert(S_mat, pmt, Si_mat);

  // Compute pve (on the transformed scale).
  gsl_blas_dgemv(CblasNoTrans, 1.0, Si_mat, q_vec, 0.0, pve);

  // Compute q_var (*n^4).
  gsl_matrix_set_zero(qvar_mat);
  s = 1;
  for (size_t i = 0; i < n_vc; i++) {
    d = gsl_vector_get(pve, i);
    gsl_matrix_view yKrKKry_sub =
        gsl_matrix_submatrix(yKrKKry, 0, i * n_vc, n_vc, n_vc);
    gsl_matrix_memcpy(tmp_mat, &yKrKKry_sub.matrix);
    gsl_matrix_scale(tmp_mat, d);
    gsl_matrix_add(qvar_mat, tmp_mat);
    s -= d;
  }
  gsl_matrix_view yKrKKry_sub =
      gsl_matrix_submatrix(yKrKKry, 0, n_vc * n_vc, n_vc, n_vc);
  gsl_matrix_memcpy(tmp_mat, &yKrKKry_sub.matrix);
  gsl_matrix_scale(tmp_mat, s);
  gsl_matrix_add(qvar_mat, tmp_mat);

  gsl_matrix_scale(qvar_mat, 2.0);

  // Compute S^{-1}var_qS^{-1}.
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Si_mat, qvar_mat, 0.0,
                 tmp_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp_mat, Si_mat, 0.0,
                 Var_mat);

  // Transform pve back to the original scale and save data.
  v_pve.clear();
  v_se_pve.clear();
  v_sigma2.clear();
  v_se_sigma2.clear();

  s = 1.0, v = 0, pve_total = 0, se_pve_total = 0;
  for (size_t i = 0; i < n_vc; i++) {
    d = gsl_vector_get(pve, i);
    v_sigma2.push_back(d * var_y_new / traceG_new[i]);
    v_pve.push_back(d * (var_y_new / traceG_new[i]) * (v_traceG[i] / var_y));
    s -= d;
    pve_total += d * (var_y_new / traceG_new[i]) * (v_traceG[i] / var_y);

    d = sqrt(gsl_matrix_get(Var_mat, i, i));
    v_se_sigma2.push_back(d * var_y_new / traceG_new[i]);
    v_se_pve.push_back(d * (var_y_new / traceG_new[i]) * (v_traceG[i] / var_y));

    for (size_t j = 0; j < n_vc; j++) {
      v += gsl_matrix_get(Var_mat, i, j);
      se_pve_total += gsl_matrix_get(Var_mat, i, j) *
                      (var_y_new / traceG_new[i]) * (v_traceG[i] / var_y) *
                      (var_y_new / traceG_new[j]) * (v_traceG[j] / var_y);
    }
  }
  v_sigma2.push_back(s * r * var_y_new);
  v_se_sigma2.push_back(sqrt(v) * r * var_y_new);
  se_pve_total = sqrt(se_pve_total);

  cout << "sigma2 = ";
  for (size_t i = 0; i < n_vc + 1; i++) {
    cout << v_sigma2[i] << " ";
  }
  cout << endl;

  cout << "se(sigma2) = ";
  for (size_t i = 0; i < n_vc + 1; i++) {
    cout << v_se_sigma2[i] << " ";
  }
  cout << endl;

  cout << "pve = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << v_pve[i] << " ";
  }
  cout << endl;

  cout << "se(pve) = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << v_se_pve[i] << " ";
  }
  cout << endl;

  if (n_vc > 1) {
    cout << "total pve = " << pve_total << endl;
    cout << "se(total pve) = " << se_pve_total << endl;
  }

  gsl_permutation_free(pmt);
  gsl_matrix_free(K_scale);
  gsl_vector_free(y_scale);
  gsl_matrix_free(Kry);
  gsl_matrix_free(yKrKKry);
  gsl_vector_free(KKry);

  // Old matrices/vectors.
  gsl_vector_free(pve);
  gsl_vector_free(se_pve);
  gsl_vector_free(q_vec);
  gsl_matrix_free(qvar_mat);
  gsl_matrix_free(tmp_mat);
  gsl_matrix_free(S_mat);
  gsl_matrix_free(Si_mat);
  gsl_matrix_free(Var_mat);

  return;
}

// REML for log(sigma2) based on the AI algorithm.
void VC::CalcVCreml(bool noconstrain, const gsl_matrix *K, const gsl_matrix *W,
                    const gsl_vector *y) {
  size_t n1 = K->size1, n2 = K->size2;
  size_t n_vc = n2 / n1;
  gsl_vector *log_sigma2 = gsl_vector_alloc(n_vc + 1);
  double d, s;

  // Set up params.
  gsl_matrix *P = gsl_matrix_alloc(n1, n1);
  gsl_vector *Py = gsl_vector_alloc(n1);
  gsl_matrix *KPy_mat = gsl_matrix_alloc(n1, n_vc + 1);
  gsl_matrix *PKPy_mat = gsl_matrix_alloc(n1, n_vc + 1);
  gsl_vector *dev1 = gsl_vector_alloc(n_vc + 1);
  gsl_matrix *dev2 = gsl_matrix_alloc(n_vc + 1, n_vc + 1);
  gsl_matrix *Hessian = gsl_matrix_alloc(n_vc + 1, n_vc + 1);
  VC_PARAM params = {K, W, y, P, Py, KPy_mat, PKPy_mat, Hessian, noconstrain};

  // Initialize sigma2/log_sigma2.
  CalcVChe(K, W, y);

  gsl_blas_ddot(y, y, &s);
  s /= (double)n1;
  for (size_t i = 0; i < n_vc + 1; i++) {
    if (noconstrain) {
      d = v_sigma2[i];
    } else {
      if (v_sigma2[i] <= 0) {
        d = log(0.1);
      } else {
        d = log(v_sigma2[i]);
      }
    }
    gsl_vector_set(log_sigma2, i, d);
  }

  cout << "iteration " << 0 << endl;
  cout << "sigma2 = ";
  for (size_t i = 0; i < n_vc + 1; i++) {
    if (noconstrain) {
      cout << gsl_vector_get(log_sigma2, i) << " ";
    } else {
      cout << exp(gsl_vector_get(log_sigma2, i)) << " ";
    }
  }
  cout << endl;

  // Set up fdf.
  gsl_multiroot_function_fdf FDF;
  FDF.n = n_vc + 1;
  FDF.params = &params;
  FDF.f = &LogRL_dev1;
  FDF.df = &LogRL_dev2;
  FDF.fdf = &LogRL_dev12;

  // Set up solver.
  int status;
  int iter = 0, max_iter = 100;

  const gsl_multiroot_fdfsolver_type *T_fdf;
  gsl_multiroot_fdfsolver *s_fdf;
  T_fdf = gsl_multiroot_fdfsolver_hybridsj;
  s_fdf = gsl_multiroot_fdfsolver_alloc(T_fdf, n_vc + 1);

  gsl_multiroot_fdfsolver_set(s_fdf, &FDF, log_sigma2);

  do {
    iter++;
    status = gsl_multiroot_fdfsolver_iterate(s_fdf);

    if (status)
      break;

    cout << "iteration " << iter << endl;
    cout << "sigma2 = ";
    for (size_t i = 0; i < n_vc + 1; i++) {
      if (noconstrain) {
        cout << gsl_vector_get(s_fdf->x, i) << " ";
      } else {
        cout << exp(gsl_vector_get(s_fdf->x, i)) << " ";
      }
    }
    cout << endl;
    status = gsl_multiroot_test_residual(s_fdf->f, 1e-3);
  } while (status == GSL_CONTINUE && iter < max_iter);

  // Obtain Hessian and Hessian inverse.
  int sig = LogRL_dev12(s_fdf->x, &params, dev1, dev2);

  gsl_permutation *pmt = gsl_permutation_alloc(n_vc + 1);
  LUDecomp(dev2, pmt, &sig);
  LUInvert(dev2, pmt, Hessian);
  gsl_permutation_free(pmt);

  // Save sigma2 and se_sigma2.
  v_sigma2.clear();
  v_se_sigma2.clear();
  for (size_t i = 0; i < n_vc + 1; i++) {
    if (noconstrain) {
      d = gsl_vector_get(s_fdf->x, i);
    } else {
      d = exp(gsl_vector_get(s_fdf->x, i));
    }
    v_sigma2.push_back(d);

    if (noconstrain) {
      d = -1.0 * gsl_matrix_get(Hessian, i, i);
    } else {
      d = -1.0 * d * d * gsl_matrix_get(Hessian, i, i);
    }
    v_se_sigma2.push_back(sqrt(d));
  }

  s = 0;
  for (size_t i = 0; i < n_vc; i++) {
    s += v_traceG[i] * v_sigma2[i];
  }
  s += v_sigma2[n_vc];

  // Compute pve.
  v_pve.clear();
  pve_total = 0;
  for (size_t i = 0; i < n_vc; i++) {
    d = v_traceG[i] * v_sigma2[i] / s;
    v_pve.push_back(d);
    pve_total += d;
  }

  // Compute se_pve; k=n_vc+1: total.
  double d1, d2;
  v_se_pve.clear();
  se_pve_total = 0;
  for (size_t k = 0; k < n_vc + 1; k++) {
    d = 0;
    for (size_t i = 0; i < n_vc + 1; i++) {
      if (noconstrain) {
        d1 = gsl_vector_get(s_fdf->x, i);
        d1 = 1;
      } else {
        d1 = exp(gsl_vector_get(s_fdf->x, i));
      }

      if (k < n_vc) {
        if (i == k) {
          d1 *= v_traceG[k] * (s - v_sigma2[k] * v_traceG[k]) / (s * s);
        } else if (i == n_vc) {
          d1 *= -1 * v_traceG[k] * v_sigma2[k] / (s * s);
        } else {
          d1 *= -1 * v_traceG[i] * v_traceG[k] * v_sigma2[k] / (s * s);
        }
      } else {
        if (i == k) {
          d1 *= -1 * (s - v_sigma2[n_vc]) / (s * s);
        } else {
          d1 *= v_traceG[i] * v_sigma2[n_vc] / (s * s);
        }
      }

      for (size_t j = 0; j < n_vc + 1; j++) {
        if (noconstrain) {
          d2 = gsl_vector_get(s_fdf->x, j);
          d2 = 1;
        } else {
          d2 = exp(gsl_vector_get(s_fdf->x, j));
        }

        if (k < n_vc) {
          if (j == k) {
            d2 *= v_traceG[k] * (s - v_sigma2[k] * v_traceG[k]) / (s * s);
          } else if (j == n_vc) {
            d2 *= -1 * v_traceG[k] * v_sigma2[k] / (s * s);
          } else {
            d2 *= -1 * v_traceG[j] * v_traceG[k] * v_sigma2[k] / (s * s);
          }
        } else {
          if (j == k) {
            d2 *= -1 * (s - v_sigma2[n_vc]) / (s * s);
          } else {
            d2 *= v_traceG[j] * v_sigma2[n_vc] / (s * s);
          }
        }

        d += -1.0 * d1 * d2 * gsl_matrix_get(Hessian, i, j);
      }
    }

    if (k < n_vc) {
      v_se_pve.push_back(sqrt(d));
    } else {
      se_pve_total = sqrt(d);
    }
  }

  gsl_multiroot_fdfsolver_free(s_fdf);

  gsl_vector_free(log_sigma2);
  gsl_matrix_free(P);
  gsl_vector_free(Py);
  gsl_matrix_free(KPy_mat);
  gsl_matrix_free(PKPy_mat);
  gsl_vector_free(dev1);
  gsl_matrix_free(dev2);
  gsl_matrix_free(Hessian);

  return;
}

// Ks are not scaled.
void VC::CalcVCacl(const gsl_matrix *K, const gsl_matrix *W,
                   const gsl_vector *y) {
  size_t n1 = K->size1, n2 = K->size2;
  size_t n_vc = n2 / n1;

  double d, y2_sum;

  // New matrices/vectors.
  gsl_matrix *K_scale = gsl_matrix_alloc(n1, n2);
  gsl_vector *y_scale = gsl_vector_alloc(n1);
  gsl_vector *y2 = gsl_vector_alloc(n1);
  gsl_vector *n1_vec = gsl_vector_alloc(n1);
  gsl_matrix *Ay = gsl_matrix_alloc(n1, n_vc);
  gsl_matrix *K2 = gsl_matrix_alloc(n1, n_vc * n_vc);
  gsl_matrix *K_tmp = gsl_matrix_alloc(n1, n1);
  gsl_matrix *V_mat = gsl_matrix_alloc(n1, n1);

  // Old matrices/vectors.
  gsl_vector *pve = gsl_vector_alloc(n_vc);
  gsl_vector *se_pve = gsl_vector_alloc(n_vc);
  gsl_vector *q_vec = gsl_vector_alloc(n_vc);
  gsl_matrix *S1 = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *S2 = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *S_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *Si_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *J_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *Var_mat = gsl_matrix_alloc(n_vc, n_vc);

  int sig;
  gsl_permutation *pmt = gsl_permutation_alloc(n_vc);

  // Center and scale K by W, and standardize K further so that all
  // diagonal elements are 1
  for (size_t i = 0; i < n_vc; i++) {
    gsl_matrix_view Kscale_sub =
        gsl_matrix_submatrix(K_scale, 0, n1 * i, n1, n1);
    gsl_matrix_const_view K_sub =
        gsl_matrix_const_submatrix(K, 0, n1 * i, n1, n1);
    gsl_matrix_memcpy(&Kscale_sub.matrix, &K_sub.matrix);

    CenterMatrix(&Kscale_sub.matrix, W);
    StandardizeMatrix(&Kscale_sub.matrix);
  }

  // Center y by W, and standardize it to have variance 1 (t(y)%*%y/n=1)
  gsl_vector_memcpy(y_scale, y);
  CenterVector(y_scale, W);

  // Compute y^2 and sum(y^2), which is also the variance of y*n1.
  gsl_vector_memcpy(y2, y_scale);
  gsl_vector_mul(y2, y_scale);

  y2_sum = 0;
  for (size_t i = 0; i < y2->size; i++) {
    y2_sum += gsl_vector_get(y2, i);
  }

  // Compute the n_vc size q vector.
  for (size_t i = 0; i < n_vc; i++) {
    gsl_matrix_const_view Kscale_sub =
        gsl_matrix_const_submatrix(K_scale, 0, n1 * i, n1, n1);

    gsl_blas_dgemv(CblasNoTrans, 1.0, &Kscale_sub.matrix, y_scale, 0.0, n1_vec);

    gsl_blas_ddot(n1_vec, y_scale, &d);
    gsl_vector_set(q_vec, i, d - y2_sum);
  }

  // Compute the n_vc by n_vc S1 and S2 matrix (and eventually
  // S=S1-\tau^{-1}S2).
  for (size_t i = 0; i < n_vc; i++) {
    gsl_matrix_const_view Kscale_sub1 =
        gsl_matrix_const_submatrix(K_scale, 0, n1 * i, n1, n1);

    for (size_t j = i; j < n_vc; j++) {
      gsl_matrix_const_view Kscale_sub2 =
          gsl_matrix_const_submatrix(K_scale, 0, n1 * j, n1, n1);

      gsl_matrix_memcpy(K_tmp, &Kscale_sub1.matrix);
      gsl_matrix_mul_elements(K_tmp, &Kscale_sub2.matrix);

      gsl_vector_set_zero(n1_vec);
      for (size_t t = 0; t < K_tmp->size1; t++) {
        gsl_vector_view Ktmp_col = gsl_matrix_column(K_tmp, t);
        gsl_vector_add(n1_vec, &Ktmp_col.vector);
      }
      gsl_vector_add_constant(n1_vec, -1.0);

      // Compute S1.
      gsl_blas_ddot(n1_vec, y2, &d);
      gsl_matrix_set(S1, i, j, 2 * d);
      if (i != j) {
        gsl_matrix_set(S1, j, i, 2 * d);
      }

      // Compute S2.
      d = 0;
      for (size_t t = 0; t < n1_vec->size; t++) {
        d += gsl_vector_get(n1_vec, t);
      }
      gsl_matrix_set(S2, i, j, d);
      if (i != j) {
        gsl_matrix_set(S2, j, i, d);
      }

      // Save information to compute J.
      gsl_vector_view K2col1 = gsl_matrix_column(K2, n_vc * i + j);
      gsl_vector_view K2col2 = gsl_matrix_column(K2, n_vc * j + i);

      gsl_vector_memcpy(&K2col1.vector, n1_vec);
      if (i != j) {
        gsl_vector_memcpy(&K2col2.vector, n1_vec);
      }
    }
  }

  // Iterate to solve tau and h's.
  size_t it = 0;
  double s = 1;
  double tau_inv = y2_sum / (double)n1 - d / ((double)n1 * ((double)n1 - 1));
  while (abs(s) > 1e-3 && it < 100) {

    // Update tau_inv.
    gsl_blas_ddot(q_vec, pve, &d);
    if (it > 0) {
      s = y2_sum / (double)n1 - d / ((double)n1 * ((double)n1 - 1)) - tau_inv;
    }
    tau_inv = y2_sum / (double)n1 - d / ((double)n1 * ((double)n1 - 1));
    if (it > 0) {
      s /= tau_inv;
    }

    // Update S.
    gsl_matrix_memcpy(S_mat, S2);
    gsl_matrix_scale(S_mat, -1 * tau_inv);
    gsl_matrix_add(S_mat, S1);

    // Update h=S^{-1}q.
    int sig;
    gsl_permutation *pmt = gsl_permutation_alloc(n_vc);
    LUDecomp(S_mat, pmt, &sig);
    LUInvert(S_mat, pmt, Si_mat);
    gsl_blas_dgemv(CblasNoTrans, 1.0, Si_mat, q_vec, 0.0, pve);

    it++;
  }

  // Compute V matrix and A matrix (K_scale is destroyed, so need to
  // compute V first).
  gsl_matrix_set_zero(V_mat);
  for (size_t i = 0; i < n_vc; i++) {
    gsl_matrix_view Kscale_sub =
        gsl_matrix_submatrix(K_scale, 0, n1 * i, n1, n1);

    // Compute V.
    gsl_matrix_memcpy(K_tmp, &Kscale_sub.matrix);
    gsl_matrix_scale(K_tmp, gsl_vector_get(pve, i));
    gsl_matrix_add(V_mat, K_tmp);

    // Compute A; the corresponding Kscale is destroyed.
    gsl_matrix_const_view K2_sub =
        gsl_matrix_const_submatrix(K2, 0, n_vc * i, n1, n_vc);
    gsl_blas_dgemv(CblasNoTrans, 1.0, &K2_sub.matrix, pve, 0.0, n1_vec);

    for (size_t t = 0; t < n1; t++) {
      gsl_matrix_set(K_scale, t, n1 * i + t, gsl_vector_get(n1_vec, t));
    }

    // Compute Ay.
    gsl_vector_view Ay_col = gsl_matrix_column(Ay, i);
    gsl_blas_dgemv(CblasNoTrans, 1.0, &Kscale_sub.matrix, y_scale, 0.0,
                   &Ay_col.vector);
  }
  gsl_matrix_scale(V_mat, tau_inv);

  // Compute J matrix.
  for (size_t i = 0; i < n_vc; i++) {
    gsl_vector_view Ay_col1 = gsl_matrix_column(Ay, i);
    gsl_blas_dgemv(CblasNoTrans, 1.0, V_mat, &Ay_col1.vector, 0.0, n1_vec);

    for (size_t j = i; j < n_vc; j++) {
      gsl_vector_view Ay_col2 = gsl_matrix_column(Ay, j);

      gsl_blas_ddot(&Ay_col2.vector, n1_vec, &d);
      gsl_matrix_set(J_mat, i, j, 2.0 * d);
      if (i != j) {
        gsl_matrix_set(J_mat, j, i, 2.0 * d);
      }
    }
  }

  // Compute H^{-1}JH^{-1} as V(\hat h), where H=S2*tau_inv; this is
  // stored in Var_mat.
  gsl_matrix_memcpy(S_mat, S2);
  gsl_matrix_scale(S_mat, tau_inv);

  LUDecomp(S_mat, pmt, &sig);
  LUInvert(S_mat, pmt, Si_mat);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Si_mat, J_mat, 0.0, S_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, S_mat, Si_mat, 0.0, Var_mat);

  // Compute variance for tau_inv.
  gsl_blas_dgemv(CblasNoTrans, 1.0, V_mat, y_scale, 0.0, n1_vec);
  gsl_blas_ddot(y_scale, n1_vec, &d);
  // auto se_tau_inv = sqrt(2 * d) / (double)n1;  UNUSED

  // Transform pve back to the original scale and save data.
  v_pve.clear();
  v_se_pve.clear();
  v_sigma2.clear();
  v_se_sigma2.clear();

  pve_total = 0, se_pve_total = 0;
  for (size_t i = 0; i < n_vc; i++) {
    d = gsl_vector_get(pve, i);
    pve_total += d;

    v_pve.push_back(d);
    v_sigma2.push_back(d * tau_inv / v_traceG[i]);

    d = sqrt(gsl_matrix_get(Var_mat, i, i));
    v_se_pve.push_back(d);
    v_se_sigma2.push_back(d * tau_inv / v_traceG[i]);

    for (size_t j = 0; j < n_vc; j++) {
      se_pve_total += gsl_matrix_get(Var_mat, i, j);
    }
  }
  v_sigma2.push_back((1 - pve_total) * tau_inv);
  v_se_sigma2.push_back(sqrt(se_pve_total) * tau_inv);
  se_pve_total = sqrt(se_pve_total);

  cout << "sigma2 = ";
  for (size_t i = 0; i < n_vc + 1; i++) {
    cout << v_sigma2[i] << " ";
  }
  cout << endl;

  cout << "se(sigma2) = ";
  for (size_t i = 0; i < n_vc + 1; i++) {
    cout << v_se_sigma2[i] << " ";
  }
  cout << endl;

  cout << "pve = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << v_pve[i] << " ";
  }
  cout << endl;

  cout << "se(pve) = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << v_se_pve[i] << " ";
  }
  cout << endl;

  if (n_vc > 1) {
    cout << "total pve = " << pve_total << endl;
    cout << "se(total pve) = " << se_pve_total << endl;
  }

  gsl_permutation_free(pmt);

  gsl_matrix_free(K_scale);
  gsl_vector_free(y_scale);
  gsl_vector_free(y2);
  gsl_vector_free(n1_vec);
  gsl_matrix_free(Ay);
  gsl_matrix_free(K2);
  gsl_matrix_free(K_tmp);
  gsl_matrix_free(V_mat);

  gsl_vector_free(pve);
  gsl_vector_free(se_pve);
  gsl_vector_free(q_vec);
  gsl_matrix_free(S1);
  gsl_matrix_free(S2);
  gsl_matrix_free(S_mat);
  gsl_matrix_free(Si_mat);
  gsl_matrix_free(J_mat);
  gsl_matrix_free(Var_mat);

  return;
}

// Read bimbam mean genotype file and compute XWz.
bool BimbamXwz(const string &file_geno, const int display_pace,
               vector<int> &indicator_idv, vector<int> &indicator_snp,
               const vector<size_t> &vec_cat, const gsl_vector *w,
               const gsl_vector *z, size_t ns_test, gsl_matrix *XWz) {
  debug_msg("entering");
  igzstream infile(file_geno.c_str(), igzstream::in);
  if (!infile) {
    cout << "error reading genotype file:" << file_geno << endl;
    return false;
  }

  string line;
  char *ch_ptr;

  size_t n_miss;
  double d, geno_mean, geno_var;

  size_t ni_test = XWz->size1;
  gsl_vector *geno = gsl_vector_alloc(ni_test);
  gsl_vector *geno_miss = gsl_vector_alloc(ni_test);
  gsl_vector *wz = gsl_vector_alloc(w->size);
  gsl_vector_memcpy(wz, z);
  gsl_vector_mul(wz, w);

  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    safeGetline(infile, line).eof();
    if (t % display_pace == 0 || t == (indicator_snp.size() - 1)) {
      ProgressBar("Reading SNPs  ", t, indicator_snp.size() - 1);
    }
    if (indicator_snp[t] == 0) {
      continue;
    }

    ch_ptr = strtok_safe((char *)line.c_str(), " , \t");
    ch_ptr = strtok_safe(NULL, " , \t");
    ch_ptr = strtok_safe(NULL, " , \t");

    geno_mean = 0.0;
    n_miss = 0;
    geno_var = 0.0;
    gsl_vector_set_all(geno_miss, 0);

    size_t j = 0;
    for (size_t i = 0; i < indicator_idv.size(); ++i) {
      if (indicator_idv[i] == 0) {
        continue;
      }
      ch_ptr = strtok_safe(NULL, " , \t");
      if (strcmp(ch_ptr, "NA") == 0) {
        gsl_vector_set(geno_miss, i, 0);
        n_miss++;
      } else {
        d = atof(ch_ptr);
        gsl_vector_set(geno, j, d);
        gsl_vector_set(geno_miss, j, 1);
        geno_mean += d;
        geno_var += d * d;
      }
      j++;
    }

    geno_mean /= (double)(ni_test - n_miss);
    geno_var += geno_mean * geno_mean * (double)n_miss;
    geno_var /= (double)ni_test;
    geno_var -= geno_mean * geno_mean;

    for (size_t i = 0; i < ni_test; ++i) {
      if (gsl_vector_get(geno_miss, i) == 0) {
        gsl_vector_set(geno, i, geno_mean);
      }
    }

    gsl_vector_add_constant(geno, -1.0 * geno_mean);

    gsl_vector_view XWz_col = gsl_matrix_column(XWz, vec_cat[ns_test]);
    d = gsl_vector_get(wz, ns_test);
    gsl_blas_daxpy(d / sqrt(geno_var), geno, &XWz_col.vector);

    ns_test++;
  }

  cout << endl;

  gsl_vector_free(geno);
  gsl_vector_free(geno_miss);
  gsl_vector_free(wz);

  infile.close();
  infile.clear();

  return true;
}

// Read PLINK bed file and compute XWz.
bool PlinkXwz(const string &file_bed, const int display_pace,
              vector<int> &indicator_idv, vector<int> &indicator_snp,
              const vector<size_t> &vec_cat, const gsl_vector *w,
              const gsl_vector *z, size_t ns_test, gsl_matrix *XWz) {
  debug_msg("entering");
  ifstream infile(file_bed.c_str(), ios::binary);
  if (!infile) {
    cout << "error reading bed file:" << file_bed << endl;
    return false;
  }

  char ch[1];
  bitset<8> b;

  size_t n_miss, ci_total, ci_test;
  double d, geno_mean, geno_var;

  size_t ni_test = XWz->size1;
  size_t ni_total = indicator_idv.size();
  gsl_vector *geno = gsl_vector_alloc(ni_test);
  gsl_vector *wz = gsl_vector_alloc(w->size);
  gsl_vector_memcpy(wz, z);
  gsl_vector_mul(wz, w);

  int n_bit;

  // Calculate n_bit and c, the number of bit for each snp.
  if (ni_total % 4 == 0) {
    n_bit = ni_total / 4;
  } else {
    n_bit = ni_total / 4 + 1;
  }

  // Print the first three magic numbers.
  for (int i = 0; i < 3; ++i) {
    infile.read(ch, 1);
    b = ch[0];
  }

  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    if (t % display_pace == 0 || t == (indicator_snp.size() - 1)) {
      ProgressBar("Reading SNPs  ", t, indicator_snp.size() - 1);
    }
    if (indicator_snp[t] == 0) {
      continue;
    }

    // n_bit, and 3 is the number of magic numbers.
    infile.seekg(t * n_bit + 3);

    // Read genotypes.
    geno_mean = 0.0;
    n_miss = 0;
    ci_total = 0;
    geno_var = 0.0;
    ci_test = 0;
    for (int i = 0; i < n_bit; ++i) {
      infile.read(ch, 1);
      b = ch[0];

      // Minor allele homozygous: 2.0; major: 0.0.
      for (size_t j = 0; j < 4; ++j) {
        if ((i == (n_bit - 1)) && ci_total == ni_total) {
          break;
        }
        if (indicator_idv[ci_total] == 0) {
          ci_total++;
          continue;
        }

        if (b[2 * j] == 0) {
          if (b[2 * j + 1] == 0) {
            gsl_vector_set(geno, ci_test, 2.0);
            geno_mean += 2.0;
            geno_var += 4.0;
          } else {
            gsl_vector_set(geno, ci_test, 1.0);
            geno_mean += 1.0;
            geno_var += 1.0;
          }
        } else {
          if (b[2 * j + 1] == 1) {
            gsl_vector_set(geno, ci_test, 0.0);
          } else {
            gsl_vector_set(geno, ci_test, -9.0);
            n_miss++;
          }
        }

        ci_test++;
        ci_total++;
      }
    }

    geno_mean /= (double)(ni_test - n_miss);
    geno_var += geno_mean * geno_mean * (double)n_miss;
    geno_var /= (double)ni_test;
    geno_var -= geno_mean * geno_mean;

    for (size_t i = 0; i < ni_test; ++i) {
      d = gsl_vector_get(geno, i);
      if (d == -9.0) {
        gsl_vector_set(geno, i, geno_mean);
      }
    }

    gsl_vector_add_constant(geno, -1.0 * geno_mean);

    gsl_vector_view XWz_col = gsl_matrix_column(XWz, vec_cat[ns_test]);
    d = gsl_vector_get(wz, ns_test);
    gsl_blas_daxpy(d / sqrt(geno_var), geno, &XWz_col.vector);

    ns_test++;
  }
  cout << endl;

  gsl_vector_free(geno);
  gsl_vector_free(wz);

  infile.close();
  infile.clear();

  return true;
}

// Read multiple genotype files and compute XWz.
bool MFILEXwz(const size_t mfile_mode, const string &file_mfile,
              const int display_pace, vector<int> &indicator_idv,
              vector<vector<int>> &mindicator_snp,
              const vector<size_t> &vec_cat, const gsl_vector *w,
              const gsl_vector *z, gsl_matrix *XWz) {
  debug_msg("entering");
  gsl_matrix_set_zero(XWz);

  igzstream infile(file_mfile.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open mfile file: " << file_mfile << endl;
    return false;
  }

  string file_name;
  size_t l = 0, ns_test = 0;

  while (!safeGetline(infile, file_name).eof()) {
    if (mfile_mode == 1) {
      file_name += ".bed";
      PlinkXwz(file_name, display_pace, indicator_idv, mindicator_snp[l],
               vec_cat, w, z, ns_test, XWz);
    } else {
      BimbamXwz(file_name, display_pace, indicator_idv, mindicator_snp[l],
                vec_cat, w, z, ns_test, XWz);
    }

    l++;
  }

  infile.close();
  infile.clear();

  return true;
}

// Read bimbam mean genotype file and compute X_i^TX_jWz.
bool BimbamXtXwz(const string &file_geno, const int display_pace,
                 vector<int> &indicator_idv, vector<int> &indicator_snp,
                 const gsl_matrix *XWz, size_t ns_test, gsl_matrix *XtXWz) {
  debug_msg("entering");
  igzstream infile(file_geno.c_str(), igzstream::in);
  if (!infile) {
    cout << "error reading genotype file:" << file_geno << endl;
    return false;
  }

  string line;
  char *ch_ptr;

  size_t n_miss;
  double d, geno_mean, geno_var;

  size_t ni_test = XWz->size1;
  gsl_vector *geno = gsl_vector_alloc(ni_test);
  gsl_vector *geno_miss = gsl_vector_alloc(ni_test);

  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    safeGetline(infile, line).eof();
    if (t % display_pace == 0 || t == (indicator_snp.size() - 1)) {
      ProgressBar("Reading SNPs  ", t, indicator_snp.size() - 1);
    }
    if (indicator_snp[t] == 0) {
      continue;
    }

    ch_ptr = strtok_safe((char *)line.c_str(), " , \t");
    ch_ptr = strtok_safe(NULL, " , \t");
    ch_ptr = strtok_safe(NULL, " , \t");

    geno_mean = 0.0;
    n_miss = 0;
    geno_var = 0.0;
    gsl_vector_set_all(geno_miss, 0);

    size_t j = 0;
    for (size_t i = 0; i < indicator_idv.size(); ++i) {
      if (indicator_idv[i] == 0) {
        continue;
      }
      ch_ptr = strtok_safe(NULL, " , \t");
      if (strcmp(ch_ptr, "NA") == 0) {
        gsl_vector_set(geno_miss, i, 0);
        n_miss++;
      } else {
        d = atof(ch_ptr);
        gsl_vector_set(geno, j, d);
        gsl_vector_set(geno_miss, j, 1);
        geno_mean += d;
        geno_var += d * d;
      }
      j++;
    }

    geno_mean /= (double)(ni_test - n_miss);
    geno_var += geno_mean * geno_mean * (double)n_miss;
    geno_var /= (double)ni_test;
    geno_var -= geno_mean * geno_mean;

    for (size_t i = 0; i < ni_test; ++i) {
      if (gsl_vector_get(geno_miss, i) == 0) {
        gsl_vector_set(geno, i, geno_mean);
      }
    }

    gsl_vector_add_constant(geno, -1.0 * geno_mean);

    for (size_t i = 0; i < XWz->size2; i++) {
      gsl_vector_const_view XWz_col = gsl_matrix_const_column(XWz, i);
      gsl_blas_ddot(geno, &XWz_col.vector, &d);
      gsl_matrix_set(XtXWz, ns_test, i, d / sqrt(geno_var));
    }

    ns_test++;
  }

  cout << endl;

  gsl_vector_free(geno);
  gsl_vector_free(geno_miss);

  infile.close();
  infile.clear();

  return true;
}

// Read PLINK bed file and compute XWz.
bool PlinkXtXwz(const string &file_bed, const int display_pace,
                vector<int> &indicator_idv, vector<int> &indicator_snp,
                const gsl_matrix *XWz, size_t ns_test, gsl_matrix *XtXWz) {
  debug_msg("entering");
  ifstream infile(file_bed.c_str(), ios::binary);
  if (!infile) {
    cout << "error reading bed file:" << file_bed << endl;
    return false;
  }

  char ch[1];
  bitset<8> b;

  size_t n_miss, ci_total, ci_test;
  double d, geno_mean, geno_var;

  size_t ni_test = XWz->size1;
  size_t ni_total = indicator_idv.size();
  gsl_vector *geno = gsl_vector_alloc(ni_test);

  int n_bit;

  // Calculate n_bit and c, the number of bit for each snp.
  if (ni_total % 4 == 0) {
    n_bit = ni_total / 4;
  } else {
    n_bit = ni_total / 4 + 1;
  }

  // Print the first three magic numbers.
  for (int i = 0; i < 3; ++i) {
    infile.read(ch, 1);
    b = ch[0];
  }

  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    if (t % display_pace == 0 || t == (indicator_snp.size() - 1)) {
      ProgressBar("Reading SNPs  ", t, indicator_snp.size() - 1);
    }
    if (indicator_snp[t] == 0) {
      continue;
    }

    // n_bit, and 3 is the number of magic numbers.
    infile.seekg(t * n_bit + 3);

    // Read genotypes.
    geno_mean = 0.0;
    n_miss = 0;
    ci_total = 0;
    geno_var = 0.0;
    ci_test = 0;
    for (int i = 0; i < n_bit; ++i) {
      infile.read(ch, 1);
      b = ch[0];

      // Minor allele homozygous: 2.0; major: 0.0;
      for (size_t j = 0; j < 4; ++j) {
        if ((i == (n_bit - 1)) && ci_total == ni_total) {
          break;
        }
        if (indicator_idv[ci_total] == 0) {
          ci_total++;
          continue;
        }

        if (b[2 * j] == 0) {
          if (b[2 * j + 1] == 0) {
            gsl_vector_set(geno, ci_test, 2.0);
            geno_mean += 2.0;
            geno_var += 4.0;
          } else {
            gsl_vector_set(geno, ci_test, 1.0);
            geno_mean += 1.0;
            geno_var += 1.0;
          }
        } else {
          if (b[2 * j + 1] == 1) {
            gsl_vector_set(geno, ci_test, 0.0);
          } else {
            gsl_vector_set(geno, ci_test, -9.0);
            n_miss++;
          }
        }

        ci_test++;
        ci_total++;
      }
    }

    geno_mean /= (double)(ni_test - n_miss);
    geno_var += geno_mean * geno_mean * (double)n_miss;
    geno_var /= (double)ni_test;
    geno_var -= geno_mean * geno_mean;

    for (size_t i = 0; i < ni_test; ++i) {
      d = gsl_vector_get(geno, i);
      if (d == -9.0) {
        gsl_vector_set(geno, i, geno_mean);
      }
    }

    gsl_vector_add_constant(geno, -1.0 * geno_mean);

    for (size_t i = 0; i < XWz->size2; i++) {
      gsl_vector_const_view XWz_col = gsl_matrix_const_column(XWz, i);
      gsl_blas_ddot(geno, &XWz_col.vector, &d);
      gsl_matrix_set(XtXWz, ns_test, i, d / sqrt(geno_var));
    }

    ns_test++;
  }
  cout << endl;

  gsl_vector_free(geno);

  infile.close();
  infile.clear();

  return true;
}

// Read multiple genotype files and compute XWz.
bool MFILEXtXwz(const size_t mfile_mode, const string &file_mfile,
                const int display_pace, vector<int> &indicator_idv,
                vector<vector<int>> &mindicator_snp, const gsl_matrix *XWz,
                gsl_matrix *XtXWz) {
  debug_msg("entering");
  gsl_matrix_set_zero(XtXWz);

  igzstream infile(file_mfile.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open mfile file: " << file_mfile << endl;
    return false;
  }

  string file_name;
  size_t l = 0, ns_test = 0;

  while (!safeGetline(infile, file_name).eof()) {
    if (mfile_mode == 1) {
      file_name += ".bed";
      PlinkXtXwz(file_name, display_pace, indicator_idv, mindicator_snp[l], XWz,
                 ns_test, XtXWz);
    } else {
      BimbamXtXwz(file_name, display_pace, indicator_idv, mindicator_snp[l],
                  XWz, ns_test, XtXWz);
    }

    l++;
  }

  infile.close();
  infile.clear();

  return true;
}

// Compute confidence intervals from summary statistics.
void CalcCIss(const gsl_matrix *Xz, const gsl_matrix *XWz,
              const gsl_matrix *XtXWz, const gsl_matrix *S_mat,
              const gsl_matrix *Svar_mat, const gsl_vector *w,
              const gsl_vector *z, const gsl_vector *s_vec,
              const vector<size_t> &vec_cat, const vector<double> &v_pve,
              vector<double> &v_se_pve, double &pve_total, double &se_pve_total,
              vector<double> &v_sigma2, vector<double> &v_se_sigma2,
              vector<double> &v_enrich, vector<double> &v_se_enrich) {
  size_t n_vc = XWz->size2, ns_test = w->size, ni_test = XWz->size1;

  // Set up matrices.
  gsl_vector *w_pve = gsl_vector_alloc(ns_test);
  gsl_vector *wz = gsl_vector_alloc(ns_test);
  gsl_vector *zwz = gsl_vector_alloc(n_vc);
  gsl_vector *zz = gsl_vector_alloc(n_vc);
  gsl_vector *Xz_pve = gsl_vector_alloc(ni_test);
  gsl_vector *WXtXWz = gsl_vector_alloc(ns_test);

  gsl_matrix *Si_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *Var_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *tmp_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *tmp_mat1 = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *VarEnrich_mat = gsl_matrix_alloc(n_vc, n_vc);
  gsl_matrix *qvar_mat = gsl_matrix_alloc(n_vc, n_vc);

  double d, s0, s1, s, s_pve, s_snp;

  // Compute wz and zwz.
  gsl_vector_memcpy(wz, z);
  gsl_vector_mul(wz, w);

  gsl_vector_set_zero(zwz);
  gsl_vector_set_zero(zz);
  for (size_t i = 0; i < w->size; i++) {
    d = gsl_vector_get(wz, i) * gsl_vector_get(z, i);
    d += gsl_vector_get(zwz, vec_cat[i]);
    gsl_vector_set(zwz, vec_cat[i], d);

    d = gsl_vector_get(z, i) * gsl_vector_get(z, i);
    d += gsl_vector_get(zz, vec_cat[i]);
    gsl_vector_set(zz, vec_cat[i], d);
  }

  // Compute wz, ve and Xz_pve.
  gsl_vector_set_zero(Xz_pve);
  s_pve = 0;
  s_snp = 0;
  for (size_t i = 0; i < n_vc; i++) {
    s_pve += v_pve[i];
    s_snp += gsl_vector_get(s_vec, i);

    gsl_vector_const_view Xz_col = gsl_matrix_const_column(Xz, i);
    gsl_blas_daxpy(v_pve[i] / gsl_vector_get(s_vec, i), &Xz_col.vector, Xz_pve);
  }

  // Set up wpve vector.
  for (size_t i = 0; i < w->size; i++) {
    d = v_pve[vec_cat[i]] / gsl_vector_get(s_vec, vec_cat[i]);
    gsl_vector_set(w_pve, i, d);
  }

  // Compute Vq (in qvar_mat).
  s0 = 1 - s_pve;
  for (size_t i = 0; i < n_vc; i++) {
    s0 += gsl_vector_get(zz, i) * v_pve[i] / gsl_vector_get(s_vec, i);
  }

  for (size_t i = 0; i < n_vc; i++) {
    s1 = s0;
    s1 -= gsl_vector_get(zwz, i) * (1 - s_pve) / gsl_vector_get(s_vec, i);

    gsl_vector_const_view XWz_col1 = gsl_matrix_const_column(XWz, i);
    gsl_vector_const_view XtXWz_col1 = gsl_matrix_const_column(XtXWz, i);

    gsl_vector_memcpy(WXtXWz, &XtXWz_col1.vector);
    gsl_vector_mul(WXtXWz, w_pve);

    gsl_blas_ddot(Xz_pve, &XWz_col1.vector, &d);
    s1 -= d / gsl_vector_get(s_vec, i);

    for (size_t j = 0; j < n_vc; j++) {
      s = s1;

      s -= gsl_vector_get(zwz, j) * (1 - s_pve) / gsl_vector_get(s_vec, j);

      gsl_vector_const_view XWz_col2 = gsl_matrix_const_column(XWz, j);
      gsl_vector_const_view XtXWz_col2 = gsl_matrix_const_column(XtXWz, j);

      gsl_blas_ddot(WXtXWz, &XtXWz_col2.vector, &d);
      s += d / (gsl_vector_get(s_vec, i) * gsl_vector_get(s_vec, j));

      gsl_blas_ddot(&XWz_col1.vector, &XWz_col2.vector, &d);
      s += d / (gsl_vector_get(s_vec, i) * gsl_vector_get(s_vec, j)) *
           (1 - s_pve);

      gsl_blas_ddot(Xz_pve, &XWz_col2.vector, &d);
      s -= d / gsl_vector_get(s_vec, j);

      gsl_matrix_set(qvar_mat, i, j, s);
    }
  }

  d = (double)(ni_test - 1);
  gsl_matrix_scale(qvar_mat, 2.0 / (d * d * d));

  // Calculate S^{-1}.
  gsl_matrix_memcpy(tmp_mat, S_mat);
  int sig;
  gsl_permutation *pmt = gsl_permutation_alloc(n_vc);
  LUDecomp(tmp_mat, pmt, &sig);
  LUInvert(tmp_mat, pmt, Si_mat);

  // Calculate variance for the estimates.
  for (size_t i = 0; i < n_vc; i++) {
    for (size_t j = i; j < n_vc; j++) {
      d = gsl_matrix_get(Svar_mat, i, j);
      d *= v_pve[i] * v_pve[j];

      d += gsl_matrix_get(qvar_mat, i, j);
      gsl_matrix_set(Var_mat, i, j, d);
      if (i != j) {
        gsl_matrix_set(Var_mat, j, i, d);
      }
    }
  }

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Si_mat, Var_mat, 0.0,
                 tmp_mat);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp_mat, Si_mat, 0.0,
                 Var_mat);

  // Compute sigma2 per snp, enrich.
  v_sigma2.clear();
  v_enrich.clear();
  for (size_t i = 0; i < n_vc; i++) {
    v_sigma2.push_back(v_pve[i] / gsl_vector_get(s_vec, i));
    v_enrich.push_back(v_pve[i] / gsl_vector_get(s_vec, i) * s_snp / s_pve);
  }

  // Compute se_pve, se_sigma2.
  for (size_t i = 0; i < n_vc; i++) {
    d = sqrt(gsl_matrix_get(Var_mat, i, i));
    v_se_pve.push_back(d);
    v_se_sigma2.push_back(d / gsl_vector_get(s_vec, i));
  }

  // Compute pve_total, se_pve_total.
  pve_total = 0;
  for (size_t i = 0; i < n_vc; i++) {
    pve_total += v_pve[i];
  }

  se_pve_total = 0;
  for (size_t i = 0; i < n_vc; i++) {
    for (size_t j = 0; j < n_vc; j++) {
      se_pve_total += gsl_matrix_get(Var_mat, i, j);
    }
  }
  se_pve_total = sqrt(se_pve_total);

  // Compute se_enrich.
  gsl_matrix_set_identity(tmp_mat);

  double d1;
  for (size_t i = 0; i < n_vc; i++) {
    d = v_pve[i] / s_pve;
    d1 = gsl_vector_get(s_vec, i);
    for (size_t j = 0; j < n_vc; j++) {
      if (i == j) {
        gsl_matrix_set(tmp_mat, i, j, (1 - d) / d1 * s_snp / s_pve);
      } else {
        gsl_matrix_set(tmp_mat, i, j, -1 * d / d1 * s_snp / s_pve);
      }
    }
  }
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp_mat, Var_mat, 0.0,
                 tmp_mat1);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmp_mat1, tmp_mat, 0.0,
                 VarEnrich_mat);

  for (size_t i = 0; i < n_vc; i++) {
    d = sqrt(gsl_matrix_get(VarEnrich_mat, i, i));
    v_se_enrich.push_back(d);
  }

  cout << "pve = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << v_pve[i] << " ";
  }
  cout << endl;

  cout << "se(pve) = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << v_se_pve[i] << " ";
  }
  cout << endl;

  cout << "sigma2 per snp = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << v_sigma2[i] << " ";
  }
  cout << endl;

  cout << "se(sigma2 per snp) = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << v_se_sigma2[i] << " ";
  }
  cout << endl;

  cout << "enrichment = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << v_enrich[i] << " ";
  }
  cout << endl;

  cout << "se(enrichment) = ";
  for (size_t i = 0; i < n_vc; i++) {
    cout << v_se_enrich[i] << " ";
  }
  cout << endl;

  // Delete matrices.
  gsl_matrix_free(Si_mat);
  gsl_matrix_free(Var_mat);
  gsl_matrix_free(VarEnrich_mat);
  gsl_matrix_free(tmp_mat);
  gsl_matrix_free(tmp_mat1);
  gsl_matrix_free(qvar_mat);

  gsl_vector_free(w_pve);
  gsl_vector_free(wz);
  gsl_vector_free(zwz);
  gsl_vector_free(WXtXWz);
  gsl_vector_free(Xz_pve);

  return;
}
