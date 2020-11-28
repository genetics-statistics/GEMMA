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

#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_vector.h"

#include "bslmm.h"
#include "lapack.h"
#include "lm.h"
#include "lmm.h"
#include "mathfunc.h"
#include "param.h"

using namespace std;

void BSLMM::CopyFromParam(PARAM &cPar) {
  a_mode = cPar.a_mode;
  d_pace = cPar.d_pace;

  file_bfile = cPar.file_bfile;
  file_geno = cPar.file_geno;
  file_out = cPar.file_out;
  path_out = cPar.path_out;

  l_min = cPar.h_min;
  l_max = cPar.h_max;
  n_region = cPar.n_region;
  pve_null = cPar.pve_null;
  pheno_mean = cPar.pheno_mean;

  time_UtZ = 0.0;
  time_Omega = 0.0;
  n_accept = 0;

  h_min = cPar.h_min;
  h_max = cPar.h_max;
  h_scale = cPar.h_scale;
  rho_min = cPar.rho_min;
  rho_max = cPar.rho_max;
  rho_scale = cPar.rho_scale;
  logp_min = cPar.logp_min;
  logp_max = cPar.logp_max;
  logp_scale = cPar.logp_scale;

  s_min = cPar.s_min;
  s_max = cPar.s_max;
  w_step = cPar.w_step;
  s_step = cPar.s_step;
  r_pace = cPar.r_pace;
  w_pace = cPar.w_pace;
  n_mh = cPar.n_mh;
  geo_mean = cPar.geo_mean;
  // randseed = cPar.randseed;
  gsl_r = cPar.gsl_r;
  trace_G = cPar.trace_G;

  ni_total = cPar.ni_total;
  ns_total = cPar.ns_total;
  ni_test = cPar.ni_test;
  ns_test = cPar.ns_test;
  n_cvt = cPar.n_cvt;

  indicator_idv = cPar.indicator_idv;
  indicator_snp = cPar.indicator_snp;
  snpInfo = cPar.snpInfo;

  return;
}

void BSLMM::CopyToParam(PARAM &cPar) {
  cPar.time_UtZ = time_UtZ;
  cPar.time_Omega = time_Omega;
  cPar.time_Proposal = time_Proposal;
  cPar.cHyp_initial = cHyp_initial;
  cPar.n_accept = n_accept;
  cPar.pheno_mean = pheno_mean;
  // cPar.randseed = randseed;

  return;
}

void BSLMM::WriteBV(const gsl_vector *bv) {
  string file_str;
  file_str = path_out + "/" + file_out;
  file_str += ".bv.txt";

  ofstream outfile(file_str.c_str(), ofstream::out);
  if (!outfile) {
    cout << "error writing file: " << file_str.c_str() << endl;
    return;
  }

  size_t t = 0;
  for (size_t i = 0; i < ni_total; ++i) {
    if (indicator_idv[i] == 0) {
      outfile << "NA" << endl;
    } else {
      outfile << scientific << setprecision(6) << gsl_vector_get(bv, t) << endl;
      t++;
    }
  }

  outfile.clear();
  outfile.close();
  return;
}

void BSLMM::WriteParam(vector<pair<double, double>> &beta_g,
                       const gsl_vector *alpha, const size_t w) {
  string file_str;
  file_str = path_out + "/" + file_out;
  file_str += ".param.txt";

  ofstream outfile(file_str.c_str(), ofstream::out);
  if (!outfile) {
    cout << "error writing file: " << file_str.c_str() << endl;
    return;
  }

  outfile << "chr"
          << "\t"
          << "rs"
          << "\t"
          << "ps"
          << "\t"
          << "n_miss"
          << "\t"
          << "alpha"
          << "\t"
          << "beta"
          << "\t"
          << "gamma" << endl;

  size_t t = 0;
  for (size_t i = 0; i < ns_total; ++i) {
    if (indicator_snp[i] == 0) {
      continue;
    }

    outfile << snpInfo[i].chr << "\t" << snpInfo[i].rs_number << "\t"
            << snpInfo[i].base_position << "\t" << snpInfo[i].n_miss << "\t";

    outfile << scientific << setprecision(6) << gsl_vector_get(alpha, t)
            << "\t";
    if (beta_g[t].second != 0) {
      outfile << beta_g[t].first / beta_g[t].second << "\t"
              << beta_g[t].second / (double)w << endl;
    } else {
      outfile << 0.0 << "\t" << 0.0 << endl;
    }
    t++;
  }

  outfile.clear();
  outfile.close();
  return;
}

void BSLMM::WriteParam(const gsl_vector *alpha) {
  string file_str;
  file_str = path_out + "/" + file_out;
  file_str += ".param.txt";

  ofstream outfile(file_str.c_str(), ofstream::out);
  if (!outfile) {
    cout << "error writing file: " << file_str.c_str() << endl;
    return;
  }

  outfile << "chr"
          << "\t"
          << "rs"
          << "\t"
          << "ps"
          << "\t"
          << "n_miss"
          << "\t"
          << "alpha"
          << "\t"
          << "beta"
          << "\t"
          << "gamma" << endl;

  size_t t = 0;
  for (size_t i = 0; i < ns_total; ++i) {
    if (indicator_snp[i] == 0) {
      continue;
    }

    outfile << snpInfo[i].chr << "\t" << snpInfo[i].rs_number << "\t"
            << snpInfo[i].base_position << "\t" << snpInfo[i].n_miss << "\t";
    outfile << scientific << setprecision(6) << gsl_vector_get(alpha, t)
            << "\t";
    outfile << 0.0 << "\t" << 0.0 << endl;
    t++;
  }

  outfile.clear();
  outfile.close();
  return;
}

void BSLMM::WriteResult(const int flag, const gsl_matrix *Result_hyp,
                        const gsl_matrix *Result_gamma, const size_t w_col) {
  string file_gamma, file_hyp;
  file_gamma = path_out + "/" + file_out;
  file_gamma += ".gamma.txt";
  file_hyp = path_out + "/" + file_out;
  file_hyp += ".hyp.txt";

  ofstream outfile_gamma, outfile_hyp;

  if (flag == 0) {
    outfile_gamma.open(file_gamma.c_str(), ofstream::out);
    outfile_hyp.open(file_hyp.c_str(), ofstream::out);
    if (!outfile_gamma) {
      cout << "error writing file: " << file_gamma << endl;
      return;
    }
    if (!outfile_hyp) {
      cout << "error writing file: " << file_hyp << endl;
      return;
    }

    outfile_hyp << "h \t pve \t rho \t pge \t pi \t n_gamma" << endl;

    for (size_t i = 0; i < s_max; ++i) {
      outfile_gamma << "s" << i << "\t";
    }
    outfile_gamma << endl;
  } else {
    outfile_gamma.open(file_gamma.c_str(), ofstream::app);
    outfile_hyp.open(file_hyp.c_str(), ofstream::app);
    if (!outfile_gamma) {
      cout << "error writing file: " << file_gamma << endl;
      return;
    }
    if (!outfile_hyp) {
      cout << "error writing file: " << file_hyp << endl;
      return;
    }

    size_t w;
    if (w_col == 0) {
      w = w_pace;
    } else {
      w = w_col;
    }

    for (size_t i = 0; i < w; ++i) {
      outfile_hyp << scientific;
      for (size_t j = 0; j < 4; ++j) {
        outfile_hyp << setprecision(6) << gsl_matrix_get(Result_hyp, i, j)
                    << "\t";
      }
      outfile_hyp << setprecision(6) << exp(gsl_matrix_get(Result_hyp, i, 4))
                  << "\t";
      outfile_hyp << (int)gsl_matrix_get(Result_hyp, i, 5) << "\t";
      outfile_hyp << endl;
    }

    for (size_t i = 0; i < w; ++i) {
      for (size_t j = 0; j < s_max; ++j) {
        outfile_gamma << (int)gsl_matrix_get(Result_gamma, i, j) << "\t";
      }
      outfile_gamma << endl;
    }
  }

  outfile_hyp.close();
  outfile_hyp.clear();
  outfile_gamma.close();
  outfile_gamma.clear();
  return;
}

void BSLMM::CalcPgamma(double *p_gamma) {
  double p, s = 0.0;
  for (size_t i = 0; i < ns_test; ++i) {
    p = 0.7 * gsl_ran_geometric_pdf(i + 1, 1.0 / geo_mean) +
        0.3 / (double)ns_test;
    p_gamma[i] = p;
    s += p;
  }
  for (size_t i = 0; i < ns_test; ++i) {
    p = p_gamma[i];
    p_gamma[i] = p / s;
  }
  return;
}

void BSLMM::SetXgamma(gsl_matrix *Xgamma, const gsl_matrix *X,
                      vector<size_t> &rank) {
  size_t pos;
  for (size_t i = 0; i < rank.size(); ++i) {
    pos = mapRank2pos[rank[i]];
    gsl_vector_view Xgamma_col = gsl_matrix_column(Xgamma, i);
    gsl_vector_const_view X_col = gsl_matrix_const_column(X, pos);
    gsl_vector_memcpy(&Xgamma_col.vector, &X_col.vector);
  }

  return;
}

double BSLMM::CalcPveLM(const gsl_matrix *UtXgamma, const gsl_vector *Uty,
                        const double sigma_a2) {
  double pve, var_y;

  gsl_matrix *Omega = gsl_matrix_alloc(UtXgamma->size2, UtXgamma->size2);
  gsl_vector *Xty = gsl_vector_alloc(UtXgamma->size2);
  gsl_vector *OiXty = gsl_vector_alloc(UtXgamma->size2);

  gsl_matrix_set_identity(Omega);
  gsl_matrix_scale(Omega, 1.0 / sigma_a2);

  lapack_dgemm((char *)"T", (char *)"N", 1.0, UtXgamma, UtXgamma, 1.0, Omega);
  gsl_blas_dgemv(CblasTrans, 1.0, UtXgamma, Uty, 0.0, Xty);

  CholeskySolve(Omega, Xty, OiXty);

  gsl_blas_ddot(Xty, OiXty, &pve);
  gsl_blas_ddot(Uty, Uty, &var_y);

  pve /= var_y;

  gsl_matrix_free(Omega);
  gsl_vector_free(Xty);
  gsl_vector_free(OiXty);

  return pve;
}

void BSLMM::InitialMCMC(const gsl_matrix *UtX, const gsl_vector *Uty,
                        vector<size_t> &rank, class HYPBSLMM &cHyp,
                        vector<pair<size_t, double>> &pos_loglr) {
  double q_genome = gsl_cdf_chisq_Qinv(0.05 / (double)ns_test, 1);

  cHyp.n_gamma = 0;
  for (size_t i = 0; i < pos_loglr.size(); ++i) {
    if (2.0 * pos_loglr[i].second > q_genome) {
      cHyp.n_gamma++;
    }
  }
  if (cHyp.n_gamma < 10) {
    cHyp.n_gamma = 10;
  }

  if (cHyp.n_gamma > s_max) {
    cHyp.n_gamma = s_max;
  }
  if (cHyp.n_gamma < s_min) {
    cHyp.n_gamma = s_min;
  }

  rank.clear();
  for (size_t i = 0; i < cHyp.n_gamma; ++i) {
    rank.push_back(i);
  }

  cHyp.logp = log((double)cHyp.n_gamma / (double)ns_test);
  cHyp.h = pve_null;

  if (cHyp.logp == 0) {
    cHyp.logp = -0.000001;
  }
  if (cHyp.h == 0) {
    cHyp.h = 0.1;
  }

  gsl_matrix *UtXgamma = gsl_matrix_alloc(ni_test, cHyp.n_gamma);
  SetXgamma(UtXgamma, UtX, rank);
  double sigma_a2;
  if (trace_G != 0) {
    sigma_a2 = cHyp.h * 1.0 /
               (trace_G * (1 - cHyp.h) * exp(cHyp.logp) * (double)ns_test);
  } else {
    sigma_a2 = cHyp.h * 1.0 / ((1 - cHyp.h) * exp(cHyp.logp) * (double)ns_test);
  }
  if (sigma_a2 == 0) {
    sigma_a2 = 0.025;
  }
  cHyp.rho = CalcPveLM(UtXgamma, Uty, sigma_a2) / cHyp.h;
  gsl_matrix_free(UtXgamma);

  if (cHyp.rho > 1.0) {
    cHyp.rho = 1.0;
  }

  if (cHyp.h < h_min) {
    cHyp.h = h_min;
  }
  if (cHyp.h > h_max) {
    cHyp.h = h_max;
  }
  if (cHyp.rho < rho_min) {
    cHyp.rho = rho_min;
  }
  if (cHyp.rho > rho_max) {
    cHyp.rho = rho_max;
  }
  if (cHyp.logp < logp_min) {
    cHyp.logp = logp_min;
  }
  if (cHyp.logp > logp_max) {
    cHyp.logp = logp_max;
  }

  cout << "initial value of h = " << cHyp.h << endl;
  cout << "initial value of rho = " << cHyp.rho << endl;
  cout << "initial value of pi = " << exp(cHyp.logp) << endl;
  cout << "initial value of |gamma| = " << cHyp.n_gamma << endl;

  return;
}

double BSLMM::CalcPosterior(const gsl_vector *Uty, const gsl_vector *K_eval,
                            gsl_vector *Utu, gsl_vector *alpha_prime,
                            class HYPBSLMM &cHyp) {
  double sigma_b2 = cHyp.h * (1.0 - cHyp.rho) / (trace_G * (1 - cHyp.h));

  gsl_vector *Utu_rand = gsl_vector_alloc(Uty->size);
  gsl_vector *weight_Hi = gsl_vector_alloc(Uty->size);

  double logpost = 0.0;
  double d, ds, uy, Hi_yy = 0, logdet_H = 0.0;
  for (size_t i = 0; i < ni_test; ++i) {
    d = gsl_vector_get(K_eval, i) * sigma_b2;
    ds = d / (d + 1.0);
    d = 1.0 / (d + 1.0);
    gsl_vector_set(weight_Hi, i, d);

    logdet_H -= log(d);
    uy = gsl_vector_get(Uty, i);
    Hi_yy += d * uy * uy;

    gsl_vector_set(Utu_rand, i, gsl_ran_gaussian(gsl_r, 1) * sqrt(ds));
  }

  // Sample tau.
  double tau = 1.0;
  if (a_mode == 11) {
    tau = gsl_ran_gamma(gsl_r, (double)ni_test / 2.0, 2.0 / Hi_yy);
  }

  // Sample alpha.
  gsl_vector_memcpy(alpha_prime, Uty);
  gsl_vector_mul(alpha_prime, weight_Hi);
  gsl_vector_scale(alpha_prime, sigma_b2);

  // Sample u.
  gsl_vector_memcpy(Utu, alpha_prime);
  gsl_vector_mul(Utu, K_eval);
  if (a_mode == 11) {
    gsl_vector_scale(Utu_rand, sqrt(1.0 / tau));
  }
  gsl_vector_add(Utu, Utu_rand);

  // For quantitative traits, calculate pve and ppe.
  if (a_mode == 11) {
    gsl_blas_ddot(Utu, Utu, &d);
    cHyp.pve = d / (double)ni_test;
    cHyp.pve /= cHyp.pve + 1.0 / tau;
    cHyp.pge = 0.0;
  }

  // Calculate likelihood.
  logpost = -0.5 * logdet_H;
  if (a_mode == 11) {
    logpost -= 0.5 * (double)ni_test * log(Hi_yy);
  } else {
    logpost -= 0.5 * Hi_yy;
  }

  logpost += ((double)cHyp.n_gamma - 1.0) * cHyp.logp +
             ((double)ns_test - (double)cHyp.n_gamma) * log(1 - exp(cHyp.logp));

  gsl_vector_free(Utu_rand);
  gsl_vector_free(weight_Hi);

  return logpost;
}

double BSLMM::CalcPosterior(const gsl_matrix *UtXgamma, const gsl_vector *Uty,
                            const gsl_vector *K_eval, gsl_vector *UtXb,
                            gsl_vector *Utu, gsl_vector *alpha_prime,
                            gsl_vector *beta, class HYPBSLMM &cHyp) {
  clock_t time_start;

  double sigma_a2 = cHyp.h * cHyp.rho /
                    (trace_G * (1 - cHyp.h) * exp(cHyp.logp) * (double)ns_test);
  double sigma_b2 = cHyp.h * (1.0 - cHyp.rho) / (trace_G * (1 - cHyp.h));

  double logpost = 0.0;
  double d, ds, uy, P_yy = 0, logdet_O = 0.0, logdet_H = 0.0;

  gsl_matrix *UtXgamma_eval =
      gsl_matrix_alloc(UtXgamma->size1, UtXgamma->size2);
  gsl_matrix *Omega = gsl_matrix_alloc(UtXgamma->size2, UtXgamma->size2);
  gsl_vector *XtHiy = gsl_vector_alloc(UtXgamma->size2);
  gsl_vector *beta_hat = gsl_vector_alloc(UtXgamma->size2);
  gsl_vector *Utu_rand = gsl_vector_alloc(UtXgamma->size1);
  gsl_vector *weight_Hi = gsl_vector_alloc(UtXgamma->size1);

  gsl_matrix_memcpy(UtXgamma_eval, UtXgamma);

  logdet_H = 0.0;
  P_yy = 0.0;
  for (size_t i = 0; i < ni_test; ++i) {
    gsl_vector_view UtXgamma_row = gsl_matrix_row(UtXgamma_eval, i);
    d = gsl_vector_get(K_eval, i) * sigma_b2;
    ds = d / (d + 1.0);
    d = 1.0 / (d + 1.0);
    gsl_vector_set(weight_Hi, i, d);

    logdet_H -= log(d);
    uy = gsl_vector_get(Uty, i);
    P_yy += d * uy * uy;
    gsl_vector_scale(&UtXgamma_row.vector, d);

    gsl_vector_set(Utu_rand, i, gsl_ran_gaussian(gsl_r, 1) * sqrt(ds));
  }

  // Calculate Omega.
  gsl_matrix_set_identity(Omega);

  time_start = clock();
  lapack_dgemm((char *)"T", (char *)"N", sigma_a2, UtXgamma_eval, UtXgamma, 1.0,
               Omega);
  time_Omega += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

  // Calculate beta_hat.
  gsl_blas_dgemv(CblasTrans, 1.0, UtXgamma_eval, Uty, 0.0, XtHiy);

  logdet_O = CholeskySolve(Omega, XtHiy, beta_hat);

  gsl_vector_scale(beta_hat, sigma_a2);

  gsl_blas_ddot(XtHiy, beta_hat, &d);
  P_yy -= d;

  // Sample tau.
  double tau = 1.0;
  if (a_mode == 11) {
    tau = gsl_ran_gamma(gsl_r, (double)ni_test / 2.0, 2.0 / P_yy);
  }

  // Sample beta.
  for (size_t i = 0; i < beta->size; i++) {
    d = gsl_ran_gaussian(gsl_r, 1);
    gsl_vector_set(beta, i, d);
  }
  gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Omega, beta);

  // This computes inv(L^T(Omega)) %*% beta.
  gsl_vector_scale(beta, sqrt(sigma_a2 / tau));
  gsl_vector_add(beta, beta_hat);
  gsl_blas_dgemv(CblasNoTrans, 1.0, UtXgamma, beta, 0.0, UtXb);

  // Sample alpha.
  gsl_vector_memcpy(alpha_prime, Uty);
  gsl_vector_sub(alpha_prime, UtXb);
  gsl_vector_mul(alpha_prime, weight_Hi);
  gsl_vector_scale(alpha_prime, sigma_b2);

  // Sample u.
  gsl_vector_memcpy(Utu, alpha_prime);
  gsl_vector_mul(Utu, K_eval);

  if (a_mode == 11) {
    gsl_vector_scale(Utu_rand, sqrt(1.0 / tau));
  }
  gsl_vector_add(Utu, Utu_rand);

  // For quantitative traits, calculate pve and pge.
  if (a_mode == 11) {
    gsl_blas_ddot(UtXb, UtXb, &d);
    cHyp.pge = d / (double)ni_test;

    gsl_blas_ddot(Utu, Utu, &d);
    cHyp.pve = cHyp.pge + d / (double)ni_test;

    if (cHyp.pve == 0) {
      cHyp.pge = 0.0;
    } else {
      cHyp.pge /= cHyp.pve;
    }
    cHyp.pve /= cHyp.pve + 1.0 / tau;
  }

  gsl_matrix_free(UtXgamma_eval);
  gsl_matrix_free(Omega);
  gsl_vector_free(XtHiy);
  gsl_vector_free(beta_hat);
  gsl_vector_free(Utu_rand);
  gsl_vector_free(weight_Hi);

  logpost = -0.5 * logdet_H - 0.5 * logdet_O;
  if (a_mode == 11) {
    logpost -= 0.5 * (double)ni_test * log(P_yy);
  } else {
    logpost -= 0.5 * P_yy;
  }
  logpost +=
      ((double)cHyp.n_gamma - 1.0) * cHyp.logp +
      ((double)ns_test - (double)cHyp.n_gamma) * log(1.0 - exp(cHyp.logp));

  return logpost;
}

// Calculate pve and pge, and calculate z_hat for case-control data.
void BSLMM::CalcCC_PVEnZ(const gsl_matrix *U, const gsl_vector *Utu,
                         gsl_vector *z_hat, class HYPBSLMM &cHyp) {
  double d;

  gsl_blas_ddot(Utu, Utu, &d);
  cHyp.pve = d / (double)ni_test;

  gsl_blas_dgemv(CblasNoTrans, 1.0, U, Utu, 0.0, z_hat);

  cHyp.pve /= cHyp.pve + 1.0;
  cHyp.pge = 0.0;

  return;
}

// Calculate pve and pge, and calculate z_hat for case-control data.
void BSLMM::CalcCC_PVEnZ(const gsl_matrix *U, const gsl_vector *UtXb,
                         const gsl_vector *Utu, gsl_vector *z_hat,
                         class HYPBSLMM &cHyp) {
  double d;
  gsl_vector *UtXbU = gsl_vector_alloc(Utu->size);

  gsl_blas_ddot(UtXb, UtXb, &d);
  cHyp.pge = d / (double)ni_test;

  gsl_blas_ddot(Utu, Utu, &d);
  cHyp.pve = cHyp.pge + d / (double)ni_test;

  gsl_vector_memcpy(UtXbU, Utu);
  gsl_vector_add(UtXbU, UtXb);
  gsl_blas_dgemv(CblasNoTrans, 1.0, U, UtXbU, 0.0, z_hat);

  if (cHyp.pve == 0) {
    cHyp.pge = 0.0;
  } else {
    cHyp.pge /= cHyp.pve;
  }

  cHyp.pve /= cHyp.pve + 1.0;

  gsl_vector_free(UtXbU);
  return;
}

void BSLMM::SampleZ(const gsl_vector *y, const gsl_vector *z_hat,
                    gsl_vector *z) {
  double d1, d2, z_rand = 0.0;
  for (size_t i = 0; i < z->size; ++i) {
    d1 = gsl_vector_get(y, i);
    d2 = gsl_vector_get(z_hat, i);

    // y is centered for case control studies.
    if (d1 <= 0.0) {

      // Control, right truncated.
      do {
        z_rand = d2 + gsl_ran_gaussian(gsl_r, 1.0);
      } while (z_rand > 0.0);
    } else {
      do {
        z_rand = d2 + gsl_ran_gaussian(gsl_r, 1.0);
      } while (z_rand < 0.0);
    }

    gsl_vector_set(z, i, z_rand);
  }

  return;
}

double BSLMM::ProposeHnRho(const class HYPBSLMM &cHyp_old,
                           class HYPBSLMM &cHyp_new, const size_t &repeat) {

  double h = cHyp_old.h, rho = cHyp_old.rho;

  double d_h = (h_max - h_min) * h_scale,
         d_rho = (rho_max - rho_min) * rho_scale;

  for (size_t i = 0; i < repeat; ++i) {
    h = h + (gsl_rng_uniform(gsl_r) - 0.5) * d_h;
    if (h < h_min) {
      h = 2 * h_min - h;
    }
    if (h > h_max) {
      h = 2 * h_max - h;
    }

    rho = rho + (gsl_rng_uniform(gsl_r) - 0.5) * d_rho;
    if (rho < rho_min) {
      rho = 2 * rho_min - rho;
    }
    if (rho > rho_max) {
      rho = 2 * rho_max - rho;
    }
  }
  cHyp_new.h = h;
  cHyp_new.rho = rho;
  return 0.0;
}

double BSLMM::ProposePi(const class HYPBSLMM &cHyp_old,
                        class HYPBSLMM &cHyp_new, const size_t &repeat) {
  double logp_old = cHyp_old.logp, logp_new = cHyp_old.logp;
  double log_ratio = 0.0;

  double d_logp = min(0.1, (logp_max - logp_min) * logp_scale);

  for (size_t i = 0; i < repeat; ++i) {
    logp_new = logp_old + (gsl_rng_uniform(gsl_r) - 0.5) * d_logp;
    if (logp_new < logp_min) {
      logp_new = 2 * logp_min - logp_new;
    }
    if (logp_new > logp_max) {
      logp_new = 2 * logp_max - logp_new;
    }
    log_ratio += logp_new - logp_old;
    logp_old = logp_new;
  }
  cHyp_new.logp = logp_new;

  return log_ratio;
}

bool comp_vec(size_t a, size_t b) { return (a < b); }

double BSLMM::ProposeGamma(const vector<size_t> &rank_old,
                           vector<size_t> &rank_new, const double *p_gamma,
                           const class HYPBSLMM &cHyp_old,
                           class HYPBSLMM &cHyp_new, const size_t &repeat) {
  map<size_t, int> mapRank2in;
  size_t r;
  double unif, logp = 0.0;
  int flag_gamma;
  size_t r_add, r_remove, col_id;

  rank_new.clear();
  if (cHyp_old.n_gamma != rank_old.size()) {
    cout << "size wrong" << endl;
  }

  if (cHyp_old.n_gamma != 0) {
    for (size_t i = 0; i < rank_old.size(); ++i) {
      r = rank_old[i];
      rank_new.push_back(r);
      mapRank2in[r] = 1;
    }
  }
  cHyp_new.n_gamma = cHyp_old.n_gamma;

  for (size_t i = 0; i < repeat; ++i) {
    unif = gsl_rng_uniform(gsl_r);

    if (unif < 0.40 && cHyp_new.n_gamma < s_max) {
      flag_gamma = 1;
    } else if (unif >= 0.40 && unif < 0.80 && cHyp_new.n_gamma > s_min) {
      flag_gamma = 2;
    } else if (unif >= 0.80 && cHyp_new.n_gamma > 0 &&
               cHyp_new.n_gamma < ns_test) {
      flag_gamma = 3;
    } else {
      flag_gamma = 4;
    }

    if (flag_gamma == 1) {

      // Add a SNP.
      do {
        r_add = gsl_ran_discrete(gsl_r, gsl_t);
      } while (mapRank2in.count(r_add) != 0);

      double prob_total = 1.0;
      for (size_t i = 0; i < cHyp_new.n_gamma; ++i) {
        r = rank_new[i];
        prob_total -= p_gamma[r];
      }

      mapRank2in[r_add] = 1;
      rank_new.push_back(r_add);
      cHyp_new.n_gamma++;
      logp += -log(p_gamma[r_add] / prob_total) - log((double)cHyp_new.n_gamma);
    } else if (flag_gamma == 2) {

      // Delete a SNP.
      col_id = gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma);
      r_remove = rank_new[col_id];

      double prob_total = 1.0;
      for (size_t i = 0; i < cHyp_new.n_gamma; ++i) {
        r = rank_new[i];
        prob_total -= p_gamma[r];
      }
      prob_total += p_gamma[r_remove];

      mapRank2in.erase(r_remove);
      rank_new.erase(rank_new.begin() + col_id);
      logp +=
          log(p_gamma[r_remove] / prob_total) + log((double)cHyp_new.n_gamma);
      cHyp_new.n_gamma--;
    } else if (flag_gamma == 3) {

      // Switch a SNP.
      col_id = gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma);
      r_remove = rank_new[col_id];

      // Be careful with the proposal.
      do {
        r_add = gsl_ran_discrete(gsl_r, gsl_t);
      } while (mapRank2in.count(r_add) != 0);

      double prob_total = 1.0;
      for (size_t i = 0; i < cHyp_new.n_gamma; ++i) {
        r = rank_new[i];
        prob_total -= p_gamma[r];
      }

      logp += log(p_gamma[r_remove] /
                  (prob_total + p_gamma[r_remove] - p_gamma[r_add]));
      logp -= log(p_gamma[r_add] / prob_total);

      mapRank2in.erase(r_remove);
      mapRank2in[r_add] = 1;
      rank_new.erase(rank_new.begin() + col_id);
      rank_new.push_back(r_add);
    } else {
      logp += 0;
    } // Do not change.
  }

  stable_sort(rank_new.begin(), rank_new.end(), comp_vec);

  mapRank2in.clear();
  return logp;
}

bool comp_lr(pair<size_t, double> a, pair<size_t, double> b) {
  return (a.second > b.second);
}

// If a_mode==13 then Uty==y.
void BSLMM::MCMC(const gsl_matrix *U, const gsl_matrix *UtX,
                 const gsl_vector *Uty, const gsl_vector *K_eval,
                 const gsl_vector *y) {
  clock_t time_start;

  class HYPBSLMM cHyp_old, cHyp_new;

  gsl_matrix *Result_hyp = gsl_matrix_alloc(w_pace, 6);
  gsl_matrix *Result_gamma = gsl_matrix_alloc(w_pace, s_max);

  gsl_vector *alpha_prime = gsl_vector_alloc(ni_test);
  gsl_vector *alpha_new = gsl_vector_alloc(ni_test);
  gsl_vector *alpha_old = gsl_vector_alloc(ni_test);
  gsl_vector *Utu = gsl_vector_alloc(ni_test);
  gsl_vector *Utu_new = gsl_vector_alloc(ni_test);
  gsl_vector *Utu_old = gsl_vector_alloc(ni_test);

  gsl_vector *UtXb_new = gsl_vector_alloc(ni_test);
  gsl_vector *UtXb_old = gsl_vector_alloc(ni_test);

  gsl_vector *z_hat = gsl_vector_alloc(ni_test);
  gsl_vector *z = gsl_vector_alloc(ni_test);
  gsl_vector *Utz = gsl_vector_alloc(ni_test);

  gsl_vector_memcpy(Utz, Uty);

  double logPost_new, logPost_old;
  double logMHratio;
  double mean_z = 0.0;

  gsl_matrix_set_zero(Result_gamma);
  gsl_vector_set_zero(Utu);
  gsl_vector_set_zero(alpha_prime);
  if (a_mode == 13) {
    pheno_mean = 0.0;
  }

  vector<pair<double, double>> beta_g;
  for (size_t i = 0; i < ns_test; i++) {
    beta_g.push_back(make_pair(0.0, 0.0));
  }

  vector<size_t> rank_new, rank_old;
  vector<double> beta_new, beta_old;

  vector<pair<size_t, double>> pos_loglr;

  time_start = clock();
  MatrixCalcLR(U, UtX, Utz, K_eval, l_min, l_max, n_region, pos_loglr);
  time_Proposal = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

  stable_sort(pos_loglr.begin(), pos_loglr.end(), comp_lr);
  for (size_t i = 0; i < ns_test; ++i) {
    mapRank2pos[i] = pos_loglr[i].first;
  }

  // Calculate proposal distribution for gamma (unnormalized),
  // and set up gsl_r and gsl_t.

  double *p_gamma = new double[ns_test];
  CalcPgamma(p_gamma);

  gsl_t = gsl_ran_discrete_preproc(ns_test, p_gamma);

  // Initial parameters.
  InitialMCMC(UtX, Utz, rank_old, cHyp_old, pos_loglr);

  cHyp_initial = cHyp_old;

  if (cHyp_old.n_gamma == 0 || cHyp_old.rho == 0) {
    logPost_old = CalcPosterior(Utz, K_eval, Utu_old, alpha_old, cHyp_old);

    beta_old.clear();
    for (size_t i = 0; i < cHyp_old.n_gamma; ++i) {
      beta_old.push_back(0);
    }
  } else {
    gsl_matrix *UtXgamma = gsl_matrix_alloc(ni_test, cHyp_old.n_gamma);
    gsl_vector *beta = gsl_vector_alloc(cHyp_old.n_gamma);
    SetXgamma(UtXgamma, UtX, rank_old);
    logPost_old = CalcPosterior(UtXgamma, Utz, K_eval, UtXb_old, Utu_old,
                                alpha_old, beta, cHyp_old);

    beta_old.clear();
    for (size_t i = 0; i < beta->size; ++i) {
      beta_old.push_back(gsl_vector_get(beta, i));
    }
    gsl_matrix_free(UtXgamma);
    gsl_vector_free(beta);
  }

  // Calculate centered z_hat, and pve.
  if (a_mode == 13) {
    time_start = clock();
    if (cHyp_old.n_gamma == 0 || cHyp_old.rho == 0) {
      CalcCC_PVEnZ(U, Utu_old, z_hat, cHyp_old);
    } else {
      CalcCC_PVEnZ(U, UtXb_old, Utu_old, z_hat, cHyp_old);
    }
    time_UtZ += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
  }

  // Start MCMC.
  int accept;
  size_t total_step = w_step + s_step;
  size_t w = 0, w_col, pos;
  size_t repeat = 0;

  for (size_t t = 0; t < total_step; ++t) {
    if (t % d_pace == 0 || t == total_step - 1) {
      ProgressBar("Running MCMC ", t, total_step - 1,
                  (double)n_accept / (double)(t * n_mh + 1));
    }

    if (a_mode == 13) {
      SampleZ(y, z_hat, z);
      mean_z = CenterVector(z);

      time_start = clock();
      gsl_blas_dgemv(CblasTrans, 1.0, U, z, 0.0, Utz);
      time_UtZ += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

      // First proposal.
      if (cHyp_old.n_gamma == 0 || cHyp_old.rho == 0) {
        logPost_old = CalcPosterior(Utz, K_eval, Utu_old, alpha_old, cHyp_old);
        beta_old.clear();
        for (size_t i = 0; i < cHyp_old.n_gamma; ++i) {
          beta_old.push_back(0);
        }
      } else {
        gsl_matrix *UtXgamma = gsl_matrix_alloc(ni_test, cHyp_old.n_gamma);
        gsl_vector *beta = gsl_vector_alloc(cHyp_old.n_gamma);
        SetXgamma(UtXgamma, UtX, rank_old);
        logPost_old = CalcPosterior(UtXgamma, Utz, K_eval, UtXb_old, Utu_old,
                                    alpha_old, beta, cHyp_old);

        beta_old.clear();
        for (size_t i = 0; i < beta->size; ++i) {
          beta_old.push_back(gsl_vector_get(beta, i));
        }
        gsl_matrix_free(UtXgamma);
        gsl_vector_free(beta);
      }
    }

    // M-H steps.
    for (size_t i = 0; i < n_mh; ++i) {
      if (gsl_rng_uniform(gsl_r) < 0.33) {
        repeat = 1 + gsl_rng_uniform_int(gsl_r, 20);
      } else {
        repeat = 1;
      }

      logMHratio = 0.0;
      logMHratio += ProposeHnRho(cHyp_old, cHyp_new, repeat);
      logMHratio +=
          ProposeGamma(rank_old, rank_new, p_gamma, cHyp_old, cHyp_new, repeat);
      logMHratio += ProposePi(cHyp_old, cHyp_new, repeat);

      if (cHyp_new.n_gamma == 0 || cHyp_new.rho == 0) {
        logPost_new = CalcPosterior(Utz, K_eval, Utu_new, alpha_new, cHyp_new);
        beta_new.clear();
        for (size_t i = 0; i < cHyp_new.n_gamma; ++i) {
          beta_new.push_back(0);
        }
      } else {
        gsl_matrix *UtXgamma = gsl_matrix_alloc(ni_test, cHyp_new.n_gamma);
        gsl_vector *beta = gsl_vector_alloc(cHyp_new.n_gamma);
        SetXgamma(UtXgamma, UtX, rank_new);
        logPost_new = CalcPosterior(UtXgamma, Utz, K_eval, UtXb_new, Utu_new,
                                    alpha_new, beta, cHyp_new);
        beta_new.clear();
        for (size_t i = 0; i < beta->size; ++i) {
          beta_new.push_back(gsl_vector_get(beta, i));
        }
        gsl_matrix_free(UtXgamma);
        gsl_vector_free(beta);
      }

      logMHratio += logPost_new - logPost_old;

      if (logMHratio > 0 || log(gsl_rng_uniform(gsl_r)) < logMHratio) {
        accept = 1;
        n_accept++;
      } else {
        accept = 0;
      }

      if (accept == 1) {
        logPost_old = logPost_new;
        rank_old.clear();
        beta_old.clear();
        if (rank_new.size() != 0) {
          for (size_t i = 0; i < rank_new.size(); ++i) {
            rank_old.push_back(rank_new[i]);
            beta_old.push_back(beta_new[i]);
          }
        }
        cHyp_old = cHyp_new;
        gsl_vector_memcpy(alpha_old, alpha_new);
        gsl_vector_memcpy(UtXb_old, UtXb_new);
        gsl_vector_memcpy(Utu_old, Utu_new);
      } else {
        cHyp_new = cHyp_old;
      }
    }

    // Calculate z_hat, and pve.
    if (a_mode == 13) {
      time_start = clock();
      if (cHyp_old.n_gamma == 0 || cHyp_old.rho == 0) {
        CalcCC_PVEnZ(U, Utu_old, z_hat, cHyp_old);
      } else {
        CalcCC_PVEnZ(U, UtXb_old, Utu_old, z_hat, cHyp_old);
      }

      // Sample mu and update z_hat.
      gsl_vector_sub(z, z_hat);
      mean_z += CenterVector(z);
      mean_z += gsl_ran_gaussian(gsl_r, sqrt(1.0 / (double)ni_test));
      gsl_vector_add_constant(z_hat, mean_z);

      time_UtZ += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
    }

    // Save data.
    if (t < w_step) {
      continue;
    } else {
      if (t % r_pace == 0) {
        w_col = w % w_pace;
        if (w_col == 0) {
          if (w == 0) {
            WriteResult(0, Result_hyp, Result_gamma, w_col);
          } else {
            WriteResult(1, Result_hyp, Result_gamma, w_col);
            gsl_matrix_set_zero(Result_hyp);
            gsl_matrix_set_zero(Result_gamma);
          }
        }

        gsl_matrix_set(Result_hyp, w_col, 0, cHyp_old.h);
        gsl_matrix_set(Result_hyp, w_col, 1, cHyp_old.pve);
        gsl_matrix_set(Result_hyp, w_col, 2, cHyp_old.rho);
        gsl_matrix_set(Result_hyp, w_col, 3, cHyp_old.pge);
        gsl_matrix_set(Result_hyp, w_col, 4, cHyp_old.logp);
        gsl_matrix_set(Result_hyp, w_col, 5, cHyp_old.n_gamma);

        for (size_t i = 0; i < cHyp_old.n_gamma; ++i) {
          pos = mapRank2pos[rank_old[i]] + 1;

          gsl_matrix_set(Result_gamma, w_col, i, pos);

          beta_g[pos - 1].first += beta_old[i];
          beta_g[pos - 1].second += 1.0;
        }

        gsl_vector_add(alpha_prime, alpha_old);
        gsl_vector_add(Utu, Utu_old);

        if (a_mode == 13) {
          pheno_mean += mean_z;
        }

        w++;
      }
    }
  }
  cout << endl;

  w_col = w % w_pace;
  WriteResult(1, Result_hyp, Result_gamma, w_col);

  gsl_matrix_free(Result_hyp);
  gsl_matrix_free(Result_gamma);

  gsl_vector_free(z_hat);
  gsl_vector_free(z);
  gsl_vector_free(Utz);
  gsl_vector_free(UtXb_new);
  gsl_vector_free(UtXb_old);
  gsl_vector_free(alpha_new);
  gsl_vector_free(alpha_old);
  gsl_vector_free(Utu_new);
  gsl_vector_free(Utu_old);

  gsl_vector_scale(alpha_prime, 1.0 / (double)w);
  gsl_vector_scale(Utu, 1.0 / (double)w);
  if (a_mode == 13) {
    pheno_mean /= (double)w;
  }

  gsl_vector *alpha = gsl_vector_alloc(ns_test);
  gsl_blas_dgemv(CblasTrans, 1.0 / (double)ns_test, UtX, alpha_prime, 0.0,
                 alpha);
  WriteParam(beta_g, alpha, w);
  gsl_vector_free(alpha);

  gsl_blas_dgemv(CblasNoTrans, 1.0, U, Utu, 0.0, alpha_prime);
  WriteBV(alpha_prime);

  gsl_vector_free(alpha_prime);
  gsl_vector_free(Utu);

  delete[] p_gamma;
  beta_g.clear();

  return;
}

void BSLMM::RidgeR(const gsl_matrix *U, const gsl_matrix *UtX,
                   const gsl_vector *Uty, const gsl_vector *eval,
                   const double lambda) {
  gsl_vector *beta = gsl_vector_alloc(UtX->size2);
  gsl_vector *H_eval = gsl_vector_alloc(Uty->size);
  gsl_vector *bv = gsl_vector_alloc(Uty->size);

  gsl_vector_memcpy(H_eval, eval);
  gsl_vector_scale(H_eval, lambda);
  gsl_vector_add_constant(H_eval, 1.0);

  gsl_vector_memcpy(bv, Uty);
  gsl_vector_div(bv, H_eval);

  gsl_blas_dgemv(CblasTrans, lambda / (double)UtX->size2, UtX, bv, 0.0, beta);
  gsl_vector_add_constant(H_eval, -1.0);
  gsl_vector_mul(H_eval, bv);
  gsl_blas_dgemv(CblasNoTrans, 1.0, U, H_eval, 0.0, bv);

  WriteParam(beta);
  WriteBV(bv);

  gsl_vector_free(H_eval);
  gsl_vector_free(beta);
  gsl_vector_free(bv);

  return;
}

// Below fits MCMC for rho=1.
void BSLMM::CalcXtX(const gsl_matrix *X, const gsl_vector *y,
                    const size_t s_size, gsl_matrix *XtX, gsl_vector *Xty) {
  time_t time_start = clock();
  gsl_matrix_const_view X_sub =
      gsl_matrix_const_submatrix(X, 0, 0, X->size1, s_size);
  gsl_matrix_view XtX_sub = gsl_matrix_submatrix(XtX, 0, 0, s_size, s_size);
  gsl_vector_view Xty_sub = gsl_vector_subvector(Xty, 0, s_size);

  lapack_dgemm((char *)"T", (char *)"N", 1.0, &X_sub.matrix, &X_sub.matrix, 0.0,
               &XtX_sub.matrix);
  gsl_blas_dgemv(CblasTrans, 1.0, &X_sub.matrix, y, 0.0, &Xty_sub.vector);

  time_Omega += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

  return;
}

void BSLMM::SetXgamma(const gsl_matrix *X, const gsl_matrix *X_old,
                      const gsl_matrix *XtX_old, const gsl_vector *Xty_old,
                      const gsl_vector *y, const vector<size_t> &rank_old,
                      const vector<size_t> &rank_new, gsl_matrix *X_new,
                      gsl_matrix *XtX_new, gsl_vector *Xty_new) {
  double d;

  // rank_old and rank_new are sorted already inside PorposeGamma
  // calculate vectors rank_remove and rank_add.
  // make sure that v_size is larger than repeat.
  size_t v_size = 20;
  vector<size_t> rank_remove(v_size), rank_add(v_size),
      rank_union(s_max + v_size);
  vector<size_t>::iterator it;

  it = set_difference(rank_old.begin(), rank_old.end(), rank_new.begin(),
                      rank_new.end(), rank_remove.begin());
  rank_remove.resize(it - rank_remove.begin());

  it = set_difference(rank_new.begin(), rank_new.end(), rank_old.begin(),
                      rank_old.end(), rank_add.begin());
  rank_add.resize(it - rank_add.begin());

  it = set_union(rank_new.begin(), rank_new.end(), rank_old.begin(),
                 rank_old.end(), rank_union.begin());
  rank_union.resize(it - rank_union.begin());

  // Map rank_remove and rank_add.
  map<size_t, int> mapRank2in_remove, mapRank2in_add;
  for (size_t i = 0; i < rank_remove.size(); i++) {
    mapRank2in_remove[rank_remove[i]] = 1;
  }
  for (size_t i = 0; i < rank_add.size(); i++) {
    mapRank2in_add[rank_add[i]] = 1;
  }

  // Obtain the subset of matrix/vector.
  gsl_matrix_const_view Xold_sub =
      gsl_matrix_const_submatrix(X_old, 0, 0, X_old->size1, rank_old.size());
  gsl_matrix_const_view XtXold_sub = gsl_matrix_const_submatrix(
      XtX_old, 0, 0, rank_old.size(), rank_old.size());
  gsl_vector_const_view Xtyold_sub =
      gsl_vector_const_subvector(Xty_old, 0, rank_old.size());

  gsl_matrix_view Xnew_sub =
      gsl_matrix_submatrix(X_new, 0, 0, X_new->size1, rank_new.size());
  gsl_matrix_view XtXnew_sub =
      gsl_matrix_submatrix(XtX_new, 0, 0, rank_new.size(), rank_new.size());
  gsl_vector_view Xtynew_sub =
      gsl_vector_subvector(Xty_new, 0, rank_new.size());

  // Get X_new and calculate XtX_new.
  if (rank_remove.size() == 0 && rank_add.size() == 0) {
    gsl_matrix_memcpy(&Xnew_sub.matrix, &Xold_sub.matrix);
    gsl_matrix_memcpy(&XtXnew_sub.matrix, &XtXold_sub.matrix);
    gsl_vector_memcpy(&Xtynew_sub.vector, &Xtyold_sub.vector);
  } else {
    size_t i_old, j_old, i_new, j_new, i_add, j_add, i_flag, j_flag;
    if (rank_add.size() == 0) {
      i_old = 0;
      i_new = 0;
      for (size_t i = 0; i < rank_union.size(); i++) {
        if (mapRank2in_remove.count(rank_old[i_old]) != 0) {
          i_old++;
          continue;
        }

        gsl_vector_view Xnew_col = gsl_matrix_column(X_new, i_new);
        gsl_vector_const_view Xcopy_col = gsl_matrix_const_column(X_old, i_old);
        gsl_vector_memcpy(&Xnew_col.vector, &Xcopy_col.vector);

        d = gsl_vector_get(Xty_old, i_old);
        gsl_vector_set(Xty_new, i_new, d);

        j_old = i_old;
        j_new = i_new;
        for (size_t j = i; j < rank_union.size(); j++) {
          if (mapRank2in_remove.count(rank_old[j_old]) != 0) {
            j_old++;
            continue;
          }

          d = gsl_matrix_get(XtX_old, i_old, j_old);

          gsl_matrix_set(XtX_new, i_new, j_new, d);
          if (i_new != j_new) {
            gsl_matrix_set(XtX_new, j_new, i_new, d);
          }

          j_old++;
          j_new++;
        }
        i_old++;
        i_new++;
      }
    } else {
      gsl_matrix *X_add = gsl_matrix_alloc(X_old->size1, rank_add.size());
      gsl_matrix *XtX_aa = gsl_matrix_alloc(X_add->size2, X_add->size2);
      gsl_matrix *XtX_ao = gsl_matrix_alloc(X_add->size2, X_old->size2);
      gsl_vector *Xty_add = gsl_vector_alloc(X_add->size2);

      // Get X_add.
      SetXgamma(X_add, X, rank_add);

      // Get t(X_add)X_add and t(X_add)X_temp.
      clock_t time_start = clock();

      // Somehow the lapack_dgemm does not work here.
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X_add, X_add, 0.0, XtX_aa);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X_add, X_old, 0.0, XtX_ao);
      gsl_blas_dgemv(CblasTrans, 1.0, X_add, y, 0.0, Xty_add);

      time_Omega += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

      // Save to X_new, XtX_new and Xty_new.
      i_old = 0;
      i_new = 0;
      i_add = 0;
      for (size_t i = 0; i < rank_union.size(); i++) {
        if (mapRank2in_remove.count(rank_old[i_old]) != 0) {
          i_old++;
          continue;
        }
        if (mapRank2in_add.count(rank_new[i_new]) != 0) {
          i_flag = 1;
        } else {
          i_flag = 0;
        }

        gsl_vector_view Xnew_col = gsl_matrix_column(X_new, i_new);
        if (i_flag == 1) {
          gsl_vector_view Xcopy_col = gsl_matrix_column(X_add, i_add);
          gsl_vector_memcpy(&Xnew_col.vector, &Xcopy_col.vector);
        } else {
          gsl_vector_const_view Xcopy_col =
              gsl_matrix_const_column(X_old, i_old);
          gsl_vector_memcpy(&Xnew_col.vector, &Xcopy_col.vector);
        }

        if (i_flag == 1) {
          d = gsl_vector_get(Xty_add, i_add);
        } else {
          d = gsl_vector_get(Xty_old, i_old);
        }
        gsl_vector_set(Xty_new, i_new, d);

        j_old = i_old;
        j_new = i_new;
        j_add = i_add;
        for (size_t j = i; j < rank_union.size(); j++) {
          if (mapRank2in_remove.count(rank_old[j_old]) != 0) {
            j_old++;
            continue;
          }
          if (mapRank2in_add.count(rank_new[j_new]) != 0) {
            j_flag = 1;
          } else {
            j_flag = 0;
          }

          if (i_flag == 1 && j_flag == 1) {
            d = gsl_matrix_get(XtX_aa, i_add, j_add);
          } else if (i_flag == 1) {
            d = gsl_matrix_get(XtX_ao, i_add, j_old);
          } else if (j_flag == 1) {
            d = gsl_matrix_get(XtX_ao, j_add, i_old);
          } else {
            d = gsl_matrix_get(XtX_old, i_old, j_old);
          }

          gsl_matrix_set(XtX_new, i_new, j_new, d);
          if (i_new != j_new) {
            gsl_matrix_set(XtX_new, j_new, i_new, d);
          }

          j_new++;
          if (j_flag == 1) {
            j_add++;
          } else {
            j_old++;
          }
        }
        i_new++;
        if (i_flag == 1) {
          i_add++;
        } else {
          i_old++;
        }
      }

      gsl_matrix_free(X_add);
      gsl_matrix_free(XtX_aa);
      gsl_matrix_free(XtX_ao);
      gsl_vector_free(Xty_add);
    }
  }

  rank_remove.clear();
  rank_add.clear();
  rank_union.clear();
  mapRank2in_remove.clear();
  mapRank2in_add.clear();

  return;
}

double BSLMM::CalcPosterior(const double yty, class HYPBSLMM &cHyp) {
  double logpost = 0.0;

  // For quantitative traits, calculate pve and pge.
  // Pve and pge for case/control data are calculted in CalcCC_PVEnZ.
  if (a_mode == 11) {
    cHyp.pve = 0.0;
    cHyp.pge = 1.0;
  }

  // Calculate likelihood.
  if (a_mode == 11) {
    logpost -= 0.5 * (double)ni_test * log(yty);
  } else {
    logpost -= 0.5 * yty;
  }

  logpost += ((double)cHyp.n_gamma - 1.0) * cHyp.logp +
             ((double)ns_test - (double)cHyp.n_gamma) * log(1 - exp(cHyp.logp));

  return logpost;
}

double BSLMM::CalcPosterior(const gsl_matrix *Xgamma, const gsl_matrix *XtX,
                            const gsl_vector *Xty, const double yty,
                            const size_t s_size, gsl_vector *Xb,
                            gsl_vector *beta, class HYPBSLMM &cHyp) {
  double sigma_a2 = cHyp.h / ((1 - cHyp.h) * exp(cHyp.logp) * (double)ns_test);
  double logpost = 0.0;
  double d, P_yy = yty, logdet_O = 0.0;

  gsl_matrix_const_view Xgamma_sub =
      gsl_matrix_const_submatrix(Xgamma, 0, 0, Xgamma->size1, s_size);
  gsl_matrix_const_view XtX_sub =
      gsl_matrix_const_submatrix(XtX, 0, 0, s_size, s_size);
  gsl_vector_const_view Xty_sub = gsl_vector_const_subvector(Xty, 0, s_size);

  gsl_matrix *Omega = gsl_matrix_alloc(s_size, s_size);
  gsl_matrix *M_temp = gsl_matrix_alloc(s_size, s_size);
  gsl_vector *beta_hat = gsl_vector_alloc(s_size);
  gsl_vector *Xty_temp = gsl_vector_alloc(s_size);

  gsl_vector_memcpy(Xty_temp, &Xty_sub.vector);

  // Calculate Omega.
  gsl_matrix_memcpy(Omega, &XtX_sub.matrix);
  gsl_matrix_scale(Omega, sigma_a2);
  gsl_matrix_set_identity(M_temp);
  gsl_matrix_add(Omega, M_temp);

  // Calculate beta_hat.
  logdet_O = CholeskySolve(Omega, Xty_temp, beta_hat);
  gsl_vector_scale(beta_hat, sigma_a2);

  gsl_blas_ddot(Xty_temp, beta_hat, &d);
  P_yy -= d;

  // Sample tau.
  double tau = 1.0;
  if (a_mode == 11) {
    tau = gsl_ran_gamma(gsl_r, (double)ni_test / 2.0, 2.0 / P_yy);
  }

  // Sample beta.
  for (size_t i = 0; i < s_size; i++) {
    d = gsl_ran_gaussian(gsl_r, 1);
    gsl_vector_set(beta, i, d);
  }
  gsl_vector_view beta_sub = gsl_vector_subvector(beta, 0, s_size);
  gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Omega,
                 &beta_sub.vector);

  // This computes inv(L^T(Omega)) %*% beta.
  gsl_vector_scale(&beta_sub.vector, sqrt(sigma_a2 / tau));
  gsl_vector_add(&beta_sub.vector, beta_hat);
  gsl_blas_dgemv(CblasNoTrans, 1.0, &Xgamma_sub.matrix, &beta_sub.vector, 0.0,
                 Xb);

  // For quantitative traits, calculate pve and pge.
  if (a_mode == 11) {
    gsl_blas_ddot(Xb, Xb, &d);
    cHyp.pve = d / (double)ni_test;
    cHyp.pve /= cHyp.pve + 1.0 / tau;
    cHyp.pge = 1.0;
  }

  logpost = -0.5 * logdet_O;
  if (a_mode == 11) {
    logpost -= 0.5 * (double)ni_test * log(P_yy);
  } else {
    logpost -= 0.5 * P_yy;
  }

  logpost +=
      ((double)cHyp.n_gamma - 1.0) * cHyp.logp +
      ((double)ns_test - (double)cHyp.n_gamma) * log(1.0 - exp(cHyp.logp));

  gsl_matrix_free(Omega);
  gsl_matrix_free(M_temp);
  gsl_vector_free(beta_hat);
  gsl_vector_free(Xty_temp);

  return logpost;
}

// Calculate pve and pge, and calculate z_hat for case-control data.
void BSLMM::CalcCC_PVEnZ(gsl_vector *z_hat, class HYPBSLMM &cHyp) {
  gsl_vector_set_zero(z_hat);
  cHyp.pve = 0.0;
  cHyp.pge = 1.0;
  return;
}

// Calculate pve and pge, and calculate z_hat for case-control data.
void BSLMM::CalcCC_PVEnZ(const gsl_vector *Xb, gsl_vector *z_hat,
                         class HYPBSLMM &cHyp) {
  double d;

  gsl_blas_ddot(Xb, Xb, &d);
  cHyp.pve = d / (double)ni_test;
  cHyp.pve /= cHyp.pve + 1.0;
  cHyp.pge = 1.0;

  gsl_vector_memcpy(z_hat, Xb);

  return;
}

// If a_mode==13, then run probit model.
void BSLMM::MCMC(const gsl_matrix *X, const gsl_vector *y) {
  clock_t time_start;
  double time_set = 0, time_post = 0;

  class HYPBSLMM cHyp_old, cHyp_new;

  gsl_matrix *Result_hyp = gsl_matrix_alloc(w_pace, 6);
  gsl_matrix *Result_gamma = gsl_matrix_alloc(w_pace, s_max);

  gsl_vector *Xb_new = gsl_vector_alloc(ni_test);
  gsl_vector *Xb_old = gsl_vector_alloc(ni_test);
  gsl_vector *z_hat = gsl_vector_alloc(ni_test);
  gsl_vector *z = gsl_vector_alloc(ni_test);

  gsl_matrix *Xgamma_old = gsl_matrix_alloc(ni_test, s_max);
  gsl_matrix *XtX_old = gsl_matrix_alloc(s_max, s_max);
  gsl_vector *Xtz_old = gsl_vector_alloc(s_max);
  gsl_vector *beta_old = gsl_vector_alloc(s_max);

  gsl_matrix *Xgamma_new = gsl_matrix_alloc(ni_test, s_max);
  gsl_matrix *XtX_new = gsl_matrix_alloc(s_max, s_max);
  gsl_vector *Xtz_new = gsl_vector_alloc(s_max);
  gsl_vector *beta_new = gsl_vector_alloc(s_max);

  double ztz = 0.0;
  gsl_vector_memcpy(z, y);

  // For quantitative traits, y is centered already in
  // gemma.cpp, but just in case.
  double mean_z = CenterVector(z);
  gsl_blas_ddot(z, z, &ztz);

  double logPost_new, logPost_old;
  double logMHratio;

  gsl_matrix_set_zero(Result_gamma);
  if (a_mode == 13) {
    pheno_mean = 0.0;
  }

  vector<pair<double, double>> beta_g;
  for (size_t i = 0; i < ns_test; i++) {
    beta_g.push_back(make_pair(0.0, 0.0));
  }

  vector<size_t> rank_new, rank_old;
  vector<pair<size_t, double>> pos_loglr;

  time_start = clock();
  MatrixCalcLmLR(X, z, pos_loglr);
  time_Proposal = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

  stable_sort(pos_loglr.begin(), pos_loglr.end(), comp_lr);
  for (size_t i = 0; i < ns_test; ++i) {
    mapRank2pos[i] = pos_loglr[i].first;
  }

  // Calculate proposal distribution for gamma (unnormalized),
  double *p_gamma = new double[ns_test];
  CalcPgamma(p_gamma);

  gsl_t = gsl_ran_discrete_preproc(ns_test, p_gamma);

  // Initial parameters.
  InitialMCMC(X, z, rank_old, cHyp_old, pos_loglr);

  cHyp_initial = cHyp_old;

  if (cHyp_old.n_gamma == 0) {
    logPost_old = CalcPosterior(ztz, cHyp_old);
  } else {
    SetXgamma(Xgamma_old, X, rank_old);
    CalcXtX(Xgamma_old, z, rank_old.size(), XtX_old, Xtz_old);
    logPost_old = CalcPosterior(Xgamma_old, XtX_old, Xtz_old, ztz,
                                rank_old.size(), Xb_old, beta_old, cHyp_old);
  }

  // Calculate centered z_hat, and pve.
  if (a_mode == 13) {
    if (cHyp_old.n_gamma == 0) {
      CalcCC_PVEnZ(z_hat, cHyp_old);
    } else {
      CalcCC_PVEnZ(Xb_old, z_hat, cHyp_old);
    }
  }

  // Start MCMC.
  int accept;
  size_t total_step = w_step + s_step;
  size_t w = 0, w_col, pos;
  size_t repeat = 0;

  for (size_t t = 0; t < total_step; ++t) {
    if (t % d_pace == 0 || t == total_step - 1) {
      ProgressBar("Running MCMC ", t, total_step - 1,
                  (double)n_accept / (double)(t * n_mh + 1));
    }

    if (a_mode == 13) {
      SampleZ(y, z_hat, z);
      mean_z = CenterVector(z);
      gsl_blas_ddot(z, z, &ztz);

      // First proposal.
      if (cHyp_old.n_gamma == 0) {
        logPost_old = CalcPosterior(ztz, cHyp_old);
      } else {
        gsl_matrix_view Xold_sub =
            gsl_matrix_submatrix(Xgamma_old, 0, 0, ni_test, rank_old.size());
        gsl_vector_view Xtz_sub =
            gsl_vector_subvector(Xtz_old, 0, rank_old.size());
        gsl_blas_dgemv(CblasTrans, 1.0, &Xold_sub.matrix, z, 0.0,
                       &Xtz_sub.vector);
        logPost_old =
            CalcPosterior(Xgamma_old, XtX_old, Xtz_old, ztz, rank_old.size(),
                          Xb_old, beta_old, cHyp_old);
      }
    }

    // M-H steps.
    for (size_t i = 0; i < n_mh; ++i) {
      if (gsl_rng_uniform(gsl_r) < 0.33) {
        repeat = 1 + gsl_rng_uniform_int(gsl_r, 20);
      } else {
        repeat = 1;
      }

      logMHratio = 0.0;
      logMHratio += ProposeHnRho(cHyp_old, cHyp_new, repeat);
      logMHratio +=
          ProposeGamma(rank_old, rank_new, p_gamma, cHyp_old, cHyp_new, repeat);
      logMHratio += ProposePi(cHyp_old, cHyp_new, repeat);

      if (cHyp_new.n_gamma == 0) {
        logPost_new = CalcPosterior(ztz, cHyp_new);
      } else {

        // This makes sure that rank_old.size() ==
        // rank_remove.size() does not happen.
        if (cHyp_new.n_gamma <= 20 || cHyp_old.n_gamma <= 20) {
          time_start = clock();
          SetXgamma(Xgamma_new, X, rank_new);
          CalcXtX(Xgamma_new, z, rank_new.size(), XtX_new, Xtz_new);
          time_set += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
        } else {
          time_start = clock();
          SetXgamma(X, Xgamma_old, XtX_old, Xtz_old, z, rank_old, rank_new,
                    Xgamma_new, XtX_new, Xtz_new);
          time_set += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
        }
        time_start = clock();
        logPost_new =
            CalcPosterior(Xgamma_new, XtX_new, Xtz_new, ztz, rank_new.size(),
                          Xb_new, beta_new, cHyp_new);
        time_post += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
      }
      logMHratio += logPost_new - logPost_old;

      if (logMHratio > 0 || log(gsl_rng_uniform(gsl_r)) < logMHratio) {
        accept = 1;
        n_accept++;
      } else {
        accept = 0;
      }

      if (accept == 1) {
        logPost_old = logPost_new;
        cHyp_old = cHyp_new;
        gsl_vector_memcpy(Xb_old, Xb_new);

        rank_old.clear();
        if (rank_new.size() != 0) {
          for (size_t i = 0; i < rank_new.size(); ++i) {
            rank_old.push_back(rank_new[i]);
          }

          gsl_matrix_view Xold_sub =
              gsl_matrix_submatrix(Xgamma_old, 0, 0, ni_test, rank_new.size());
          gsl_matrix_view XtXold_sub = gsl_matrix_submatrix(
              XtX_old, 0, 0, rank_new.size(), rank_new.size());
          gsl_vector_view Xtzold_sub =
              gsl_vector_subvector(Xtz_old, 0, rank_new.size());
          gsl_vector_view betaold_sub =
              gsl_vector_subvector(beta_old, 0, rank_new.size());

          gsl_matrix_view Xnew_sub =
              gsl_matrix_submatrix(Xgamma_new, 0, 0, ni_test, rank_new.size());
          gsl_matrix_view XtXnew_sub = gsl_matrix_submatrix(
              XtX_new, 0, 0, rank_new.size(), rank_new.size());
          gsl_vector_view Xtznew_sub =
              gsl_vector_subvector(Xtz_new, 0, rank_new.size());
          gsl_vector_view betanew_sub =
              gsl_vector_subvector(beta_new, 0, rank_new.size());

          gsl_matrix_memcpy(&Xold_sub.matrix, &Xnew_sub.matrix);
          gsl_matrix_memcpy(&XtXold_sub.matrix, &XtXnew_sub.matrix);
          gsl_vector_memcpy(&Xtzold_sub.vector, &Xtznew_sub.vector);
          gsl_vector_memcpy(&betaold_sub.vector, &betanew_sub.vector);
        }
      } else {
        cHyp_new = cHyp_old;
      }
    }

    // Calculate z_hat, and pve.
    if (a_mode == 13) {
      if (cHyp_old.n_gamma == 0) {
        CalcCC_PVEnZ(z_hat, cHyp_old);
      } else {
        CalcCC_PVEnZ(Xb_old, z_hat, cHyp_old);
      }

      // Sample mu and update z_hat.
      gsl_vector_sub(z, z_hat);
      mean_z += CenterVector(z);
      mean_z += gsl_ran_gaussian(gsl_r, sqrt(1.0 / (double)ni_test));

      gsl_vector_add_constant(z_hat, mean_z);
    }

    // Save data.
    if (t < w_step) {
      continue;
    } else {
      if (t % r_pace == 0) {
        w_col = w % w_pace;
        if (w_col == 0) {
          if (w == 0) {
            WriteResult(0, Result_hyp, Result_gamma, w_col);
          } else {
            WriteResult(1, Result_hyp, Result_gamma, w_col);
            gsl_matrix_set_zero(Result_hyp);
            gsl_matrix_set_zero(Result_gamma);
          }
        }

        gsl_matrix_set(Result_hyp, w_col, 0, cHyp_old.h);
        gsl_matrix_set(Result_hyp, w_col, 1, cHyp_old.pve);
        gsl_matrix_set(Result_hyp, w_col, 2, cHyp_old.rho);
        gsl_matrix_set(Result_hyp, w_col, 3, cHyp_old.pge);
        gsl_matrix_set(Result_hyp, w_col, 4, cHyp_old.logp);
        gsl_matrix_set(Result_hyp, w_col, 5, cHyp_old.n_gamma);

        for (size_t i = 0; i < cHyp_old.n_gamma; ++i) {
          pos = mapRank2pos[rank_old[i]] + 1;
          gsl_matrix_set(Result_gamma, w_col, i, pos);

          beta_g[pos - 1].first += gsl_vector_get(beta_old, i);
          beta_g[pos - 1].second += 1.0;
        }

        if (a_mode == 13) {
          pheno_mean += mean_z;
        }

        w++;
      }
    }
  }
  cout << endl;

  cout << "time on selecting Xgamma: " << time_set << endl;
  cout << "time on calculating posterior: " << time_post << endl;

  w_col = w % w_pace;
  WriteResult(1, Result_hyp, Result_gamma, w_col);

  gsl_vector *alpha = gsl_vector_alloc(ns_test);
  gsl_vector_set_zero(alpha);
  WriteParam(beta_g, alpha, w);
  gsl_vector_free(alpha);

  gsl_matrix_free(Result_hyp);
  gsl_matrix_free(Result_gamma);

  gsl_vector_free(z_hat);
  gsl_vector_free(z);
  gsl_vector_free(Xb_new);
  gsl_vector_free(Xb_old);

  gsl_matrix_free(Xgamma_old);
  gsl_matrix_free(XtX_old);
  gsl_vector_free(Xtz_old);
  gsl_vector_free(beta_old);

  gsl_matrix_free(Xgamma_new);
  gsl_matrix_free(XtX_new);
  gsl_vector_free(Xtz_new);
  gsl_vector_free(beta_new);

  delete[] p_gamma;
  beta_g.clear();

  return;
}
