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

#include <assert.h>
#include <bitset>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_min.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_vector.h"

#include "fastblas.h"
#include "gzstream.h"
#include "gemma_io.h"
#include "lapack.h"
#include "mathfunc.h"
#include "lmm.h"
#include "mvlmm.h"

using namespace std;

// In this file, X, Y are already transformed (i.e. UtX and UtY).
void MVLMM::CopyFromParam(PARAM &cPar) {
  a_mode = cPar.a_mode;
  d_pace = cPar.d_pace;

  file_bfile = cPar.file_bfile;
  file_geno = cPar.file_geno;
  file_out = cPar.file_out;
  path_out = cPar.path_out;

  l_min = cPar.l_min;
  l_max = cPar.l_max;
  n_region = cPar.n_region;
  p_nr = cPar.p_nr;
  em_iter = cPar.em_iter;
  nr_iter = cPar.nr_iter;
  em_prec = cPar.em_prec;
  nr_prec = cPar.nr_prec;
  crt = cPar.crt;

  Vg_remle_null = cPar.Vg_remle_null;
  Ve_remle_null = cPar.Ve_remle_null;
  Vg_mle_null = cPar.Vg_mle_null;
  Ve_mle_null = cPar.Ve_mle_null;

  time_UtX = 0.0;
  time_opt = 0.0;

  ni_total = cPar.ni_total;
  ns_total = cPar.ns_total;
  ni_test = cPar.ni_test;
  ns_test = cPar.ns_test;
  n_cvt = cPar.n_cvt;

  n_ph = cPar.n_ph;

  indicator_idv = cPar.indicator_idv;
  indicator_snp = cPar.indicator_snp;
  snpInfo = cPar.snpInfo;

  return;
}

void MVLMM::CopyToParam(PARAM &cPar) {
  cPar.time_UtX = time_UtX;
  cPar.time_opt = time_opt;

  cPar.Vg_remle_null = Vg_remle_null;
  cPar.Ve_remle_null = Ve_remle_null;
  cPar.Vg_mle_null = Vg_mle_null;
  cPar.Ve_mle_null = Ve_mle_null;

  cPar.VVg_remle_null = VVg_remle_null;
  cPar.VVe_remle_null = VVe_remle_null;
  cPar.VVg_mle_null = VVg_mle_null;
  cPar.VVe_mle_null = VVe_mle_null;

  cPar.beta_remle_null = beta_remle_null;
  cPar.se_beta_remle_null = se_beta_remle_null;
  cPar.beta_mle_null = beta_mle_null;
  cPar.se_beta_mle_null = se_beta_mle_null;

  cPar.logl_remle_H0 = logl_remle_H0;
  cPar.logl_mle_H0 = logl_mle_H0;
  return;
}

void MVLMM::WriteFiles() {
  string file_str;
  file_str = path_out + "/" + file_out;
  file_str += ".assoc.txt";

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
          << "allele1"
          << "\t"
          << "allele0"
          << "\t"
          << "af"
          << "\t";

  for (size_t i = 0; i < n_ph; i++) {
    outfile << "beta_" << i + 1 << "\t";
  }
  for (size_t i = 0; i < n_ph; i++) {
    for (size_t j = i; j < n_ph; j++) {
      outfile << "Vbeta_" << i + 1 << "_" << j + 1 << "\t";
    }
  }

  if (a_mode == 1) {
    outfile << "p_wald" << endl;
  } else if (a_mode == 2) {
    outfile << "p_lrt" << endl;
  } else if (a_mode == 3) {
    outfile << "p_score" << endl;
  } else if (a_mode == 4) {
    outfile << "p_wald"
            << "\t"
            << "p_lrt"
            << "\t"
            << "p_score" << endl;
  } else {
  }

  size_t t = 0, c = 0;
  for (size_t i = 0; i < snpInfo.size(); ++i) {
    if (indicator_snp[i] == 0) {
      continue;
    }

    outfile << snpInfo[i].chr << "\t" << snpInfo[i].rs_number << "\t"
            << snpInfo[i].base_position << "\t" << snpInfo[i].n_miss << "\t"
            << snpInfo[i].a_minor << "\t" << snpInfo[i].a_major << "\t" << fixed
            << setprecision(3) << snpInfo[i].maf << "\t";

    outfile << scientific << setprecision(6);

    for (size_t i = 0; i < n_ph; i++) {
      outfile << sumStat[t].v_beta[i] << "\t";
    }

    c = 0;
    for (size_t i = 0; i < n_ph; i++) {
      for (size_t j = i; j < n_ph; j++) {
        outfile << sumStat[t].v_Vbeta[c] << "\t";
        c++;
      }
    }

    if (a_mode == 1) {
      outfile << sumStat[t].p_wald << endl;
    } else if (a_mode == 2) {
      outfile << sumStat[t].p_lrt << endl;
    } else if (a_mode == 3) {
      outfile << sumStat[t].p_score << endl;
    } else if (a_mode == 4) {
      outfile << sumStat[t].p_wald << "\t" << sumStat[t].p_lrt << "\t"
              << sumStat[t].p_score << endl;
    } else {
    }

    t++;
  }

  outfile.close();
  outfile.clear();
  return;
}

// Below are functions for EM algorithm.
double EigenProc(const gsl_matrix *V_g, const gsl_matrix *V_e, gsl_vector *D_l,
                 gsl_matrix *UltVeh, gsl_matrix *UltVehi) {
  size_t d_size = V_g->size1;
  double d, logdet_Ve = 0.0;

  // Eigen decomposition of V_e.
  gsl_matrix *Lambda = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *V_e_temp = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *V_e_h = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *V_e_hi = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *VgVehi = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *U_l = gsl_matrix_alloc(d_size, d_size);

  gsl_matrix_memcpy(V_e_temp, V_e);
  EigenDecomp(V_e_temp, U_l, D_l, 0);

  // Calculate V_e_h and V_e_hi.
  gsl_matrix_set_zero(V_e_h);
  gsl_matrix_set_zero(V_e_hi);
  for (size_t i = 0; i < d_size; i++) {
    d = gsl_vector_get(D_l, i);
    if (d <= 0) {
      continue;
    }
    logdet_Ve += safe_log(d);

    gsl_vector_view U_col = gsl_matrix_column(U_l, i);
    d = safe_sqrt(d);
    gsl_blas_dsyr(CblasUpper, d, &U_col.vector, V_e_h);
    d = 1.0 / d;
    gsl_blas_dsyr(CblasUpper, d, &U_col.vector, V_e_hi);
  }

  // Copy the upper part to lower part.
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j < i; j++) {
      gsl_matrix_set(V_e_h, i, j, gsl_matrix_get(V_e_h, j, i));
      gsl_matrix_set(V_e_hi, i, j, gsl_matrix_get(V_e_hi, j, i));
    }
  }

  // Calculate Lambda=V_ehi V_g V_ehi.
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V_g, V_e_hi, 0.0, VgVehi);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V_e_hi, VgVehi, 0.0, Lambda);

  // Eigen decomposition of Lambda.
  EigenDecomp(Lambda, U_l, D_l, 0);

  for (size_t i = 0; i < d_size; i++) {
    d = gsl_vector_get(D_l, i);
    if (d < 0) {
      gsl_vector_set(D_l, i, 0);
    }
  }


  // Calculate UltVeh and UltVehi.
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, U_l, V_e_h, 0.0, UltVeh);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, U_l, V_e_hi, 0.0, UltVehi);

  // free memory
  gsl_matrix_free(Lambda);
  gsl_matrix_free(V_e_temp);
  gsl_matrix_free(V_e_h);
  gsl_matrix_free(V_e_hi);
  gsl_matrix_free(VgVehi);
  gsl_matrix_free(U_l);

  return logdet_Ve;
}

// Qi=(\sum_{k=1}^n x_kx_k^T\otimes(delta_k*Dl+I)^{-1} )^{-1}.
double CalcQi(const gsl_vector *eval, const gsl_vector *D_l,
              const gsl_matrix *X, gsl_matrix *Qi) {
  size_t n_size = eval->size, d_size = D_l->size, dc_size = Qi->size1;
  size_t c_size = dc_size / d_size;

  double delta, dl, d1, d2, d, logdet_Q;

  gsl_matrix *Q = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix_set_zero(Q);

  for (size_t i = 0; i < c_size; i++) {
    for (size_t j = 0; j < c_size; j++) {
      for (size_t l = 0; l < d_size; l++) {
        dl = gsl_vector_get(D_l, l);

        if (j < i) {
          d = gsl_matrix_get(Q, j * d_size + l, i * d_size + l);
        } else {
          d = 0.0;
          for (size_t k = 0; k < n_size; k++) {
            d1 = gsl_matrix_get(X, i, k);
            d2 = gsl_matrix_get(X, j, k);
            delta = gsl_vector_get(eval, k);
            d += d1 * d2 / (dl * delta + 1.0); // @@
          }
        }

        gsl_matrix_set(Q, i * d_size + l, j * d_size + l, d);
      }
    }
  }

  // Calculate LU decomposition of Q, and invert Q and calculate |Q|.
  int sig;
  gsl_permutation *pmt = gsl_permutation_alloc(dc_size);
  LUDecomp(Q, pmt, &sig);
  LUInvert(Q, pmt, Qi);

  logdet_Q = LULndet(Q);

  gsl_matrix_free(Q);
  gsl_permutation_free(pmt);

  return logdet_Q;
}

// xHiy=\sum_{k=1}^n x_k\otimes ((delta_k*Dl+I)^{-1}Ul^TVe^{-1/2}y.
//
// FIXME: mvlmm spends a massive amount of time here
void CalcXHiY(const gsl_vector *eval, const gsl_vector *D_l,
              const gsl_matrix *X, const gsl_matrix *UltVehiY,
              gsl_vector *xHiy) {
  // debug_msg("enter");
  size_t n_size = eval->size, c_size = X->size1, d_size = D_l->size;

  // gsl_vector_set_zero(xHiy);

  double x, delta, dl, y, d;
  for (size_t i = 0; i < d_size; i++) {
    dl = gsl_vector_get(D_l, i);
    for (size_t j = 0; j < c_size; j++) {
      d = 0.0;
      for (size_t k = 0; k < n_size; k++) {
        x = gsl_matrix_get(X, j, k);
        y = gsl_matrix_get(UltVehiY, i, k);
        delta = gsl_vector_get(eval, k);
        d += x * y / (delta * dl + 1.0);
      }
      gsl_vector_set(xHiy, j * d_size + i, d);
    }
  }
  // debug_msg("exit");

  return;
}

// OmegaU=D_l/(delta Dl+I)^{-1}
// OmegaE=delta D_l/(delta Dl+I)^{-1}
void CalcOmega(const gsl_vector *eval, const gsl_vector *D_l,
               gsl_matrix *OmegaU, gsl_matrix *OmegaE) {
  size_t n_size = eval->size, d_size = D_l->size;
  double delta, dl, d_u, d_e;

  for (size_t k = 0; k < n_size; k++) {
    delta = gsl_vector_get(eval, k);
    for (size_t i = 0; i < d_size; i++) {
      dl = gsl_vector_get(D_l, i);

      d_u = dl / (delta * dl + 1.0);  // @@
      d_e = delta * d_u;

      gsl_matrix_set(OmegaU, i, k, d_u);
      gsl_matrix_set(OmegaE, i, k, d_e);
    }
  }

  return;
}

void UpdateU(const gsl_matrix *OmegaE, const gsl_matrix *UltVehiY,
             const gsl_matrix *UltVehiBX, gsl_matrix *UltVehiU) {
  gsl_matrix_memcpy(UltVehiU, UltVehiY);
  gsl_matrix_sub(UltVehiU, UltVehiBX);

  gsl_matrix_mul_elements(UltVehiU, OmegaE);
  return;
}

void UpdateE(const gsl_matrix *UltVehiY, const gsl_matrix *UltVehiBX,
             const gsl_matrix *UltVehiU, gsl_matrix *UltVehiE) {
  gsl_matrix_memcpy(UltVehiE, UltVehiY);
  gsl_matrix_sub(UltVehiE, UltVehiBX);
  gsl_matrix_sub(UltVehiE, UltVehiU);

  return;
}

void UpdateL_B(const gsl_matrix *X, const gsl_matrix *XXti,
               const gsl_matrix *UltVehiY, const gsl_matrix *UltVehiU,
               gsl_matrix *UltVehiBX, gsl_matrix *UltVehiB) {
  size_t c_size = X->size1, d_size = UltVehiY->size1;

  gsl_matrix *YUX = gsl_matrix_alloc(d_size, c_size);

  gsl_matrix_memcpy(UltVehiBX, UltVehiY);
  gsl_matrix_sub(UltVehiBX, UltVehiU);

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, UltVehiBX, X, 0.0, YUX);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, YUX, XXti, 0.0, UltVehiB);

  gsl_matrix_free(YUX);

  return;
}

void UpdateRL_B(const gsl_vector *xHiy, const gsl_matrix *Qi,
                gsl_matrix *UltVehiB) {
  size_t d_size = UltVehiB->size1, c_size = UltVehiB->size2,
         dc_size = Qi->size1;

  gsl_vector *b = gsl_vector_alloc(dc_size);

  // Calculate b=Qiv.
  gsl_blas_dgemv(CblasNoTrans, 1.0, Qi, xHiy, 0.0, b);

  // Copy b to UltVehiB.
  for (size_t i = 0; i < c_size; i++) {
    gsl_vector_view UltVehiB_col = gsl_matrix_column(UltVehiB, i);
    gsl_vector_const_view b_subcol =
        gsl_vector_const_subvector(b, i * d_size, d_size);
    gsl_vector_memcpy(&UltVehiB_col.vector, &b_subcol.vector);
  }

  gsl_vector_free(b);

  return;
}

void UpdateV(const gsl_vector *eval, const gsl_matrix *U, const gsl_matrix *E,
             const gsl_matrix *Sigma_uu, const gsl_matrix *Sigma_ee,
             gsl_matrix *V_g, gsl_matrix *V_e) {
  size_t n_size = eval->size, d_size = U->size1;

  gsl_matrix_set_zero(V_g);
  gsl_matrix_set_zero(V_e);

  double delta;

  // Calculate the first part: UD^{-1}U^T and EE^T.
  for (size_t k = 0; k < n_size; k++) {
    delta = gsl_vector_get(eval, k);
    if (delta == 0) {
      continue;
    }

    gsl_vector_const_view U_col = gsl_matrix_const_column(U, k);
    gsl_blas_dsyr(CblasUpper, 1.0 / delta, &U_col.vector, V_g);
  }

  gsl_blas_dsyrk(CblasUpper, CblasNoTrans, 1.0, E, 0.0, V_e);

  // Copy the upper part to lower part.
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j < i; j++) {
      gsl_matrix_set(V_g, i, j, gsl_matrix_get(V_g, j, i));
      gsl_matrix_set(V_e, i, j, gsl_matrix_get(V_e, j, i));
    }
  }

  // Add Sigma.
  gsl_matrix_add(V_g, Sigma_uu);
  gsl_matrix_add(V_e, Sigma_ee);

  // Scale by 1/n.
  gsl_matrix_scale(V_g, 1.0 / (double)n_size);
  gsl_matrix_scale(V_e, 1.0 / (double)n_size);

  return;
}

void CalcSigma(const char func_name, const gsl_vector *eval,
               const gsl_vector *D_l, const gsl_matrix *X,
               const gsl_matrix *OmegaU, const gsl_matrix *OmegaE,
               const gsl_matrix *UltVeh, const gsl_matrix *Qi,
               gsl_matrix *Sigma_uu, gsl_matrix *Sigma_ee) {
  if (func_name != 'R' && func_name != 'L' && func_name != 'r' &&
      func_name != 'l') {
    cout << "func_name only takes 'R' or 'L': 'R' for log-restricted "
         << "likelihood, 'L' for log-likelihood." << endl;
    return;
  }

  size_t n_size = eval->size, c_size = X->size1;
  size_t d_size = D_l->size, dc_size = Qi->size1;

  gsl_matrix_set_zero(Sigma_uu);
  gsl_matrix_set_zero(Sigma_ee);

  double delta, dl, x, d;

  // Calculate the first diagonal term.
  gsl_vector_view Suu_diag = gsl_matrix_diagonal(Sigma_uu);
  gsl_vector_view See_diag = gsl_matrix_diagonal(Sigma_ee);

  for (size_t k = 0; k < n_size; k++) {
    gsl_vector_const_view OmegaU_col = gsl_matrix_const_column(OmegaU, k);
    gsl_vector_const_view OmegaE_col = gsl_matrix_const_column(OmegaE, k);

    gsl_vector_add(&Suu_diag.vector, &OmegaU_col.vector);
    gsl_vector_add(&See_diag.vector, &OmegaE_col.vector);
  }

  // Calculate the second term for REML.
  if (func_name == 'R' || func_name == 'r') {
    gsl_matrix *M_u = gsl_matrix_alloc(dc_size, d_size);
    gsl_matrix *M_e = gsl_matrix_alloc(dc_size, d_size);
    gsl_matrix *QiM = gsl_matrix_alloc(dc_size, d_size);

    gsl_matrix_set_zero(M_u);
    gsl_matrix_set_zero(M_e);

    for (size_t k = 0; k < n_size; k++) {
      delta = gsl_vector_get(eval, k);

      for (size_t i = 0; i < d_size; i++) {
        dl = gsl_vector_get(D_l, i);
        for (size_t j = 0; j < c_size; j++) {
          x = gsl_matrix_get(X, j, k);
          d = x / (delta * dl + 1.0);
          gsl_matrix_set(M_e, j * d_size + i, i, d);
          gsl_matrix_set(M_u, j * d_size + i, i, d * dl);
        }
      }
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Qi, M_u, 0.0, QiM);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, delta, M_u, QiM, 1.0, Sigma_uu);

      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Qi, M_e, 0.0, QiM);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, M_e, QiM, 1.0, Sigma_ee);
    }

    gsl_matrix_free(M_u);
    gsl_matrix_free(M_e);
    gsl_matrix_free(QiM);
  }

  // Multiply both sides by VehUl.
  gsl_matrix *M = gsl_matrix_alloc(d_size, d_size);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Sigma_uu, UltVeh, 0.0, M);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, UltVeh, M, 0.0, Sigma_uu);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Sigma_ee, UltVeh, 0.0, M);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, UltVeh, M, 0.0, Sigma_ee);

  gsl_matrix_free(M);
  return;
}

// 'R' for restricted likelihood and 'L' for likelihood.
// 'R' update B and 'L' don't.
// only calculate -0.5*\sum_{k=1}^n|H_k|-0.5yPxy.
double MphCalcLogL(const gsl_vector *eval, const gsl_vector *xHiy,
                   const gsl_vector *D_l, const gsl_matrix *UltVehiY,
                   const gsl_matrix *Qi) {
  size_t n_size = eval->size, d_size = D_l->size, dc_size = Qi->size1;
  double logl = 0.0, delta, dl, y, d;

  // Calculate yHiy+log|H_k|.
  for (size_t k = 0; k < n_size; k++) {
    delta = gsl_vector_get(eval, k);
    for (size_t i = 0; i < d_size; i++) {
      y = gsl_matrix_get(UltVehiY, i, k);
      dl = gsl_vector_get(D_l, i);
      d = delta * dl + 1.0;

      logl += y * y / d + safe_log(d);
    }
  }

  // Calculate the rest of yPxy.
  gsl_vector *Qiv = gsl_vector_alloc(dc_size);

  gsl_blas_dgemv(CblasNoTrans, 1.0, Qi, xHiy, 0.0, Qiv);
  gsl_blas_ddot(xHiy, Qiv, &d);

  logl -= d;

  gsl_vector_free(Qiv);

  return -0.5 * logl;
}

// Y is a dxn matrix, X is a cxn matrix, B is a dxc matrix, V_g is a
// dxd matrix, V_e is a dxd matrix, eval is a size n vector
//'R' for restricted likelihood and 'L' for likelihood.
double MphEM(const char func_name, const size_t max_iter, const double max_prec,
             const gsl_vector *eval, const gsl_matrix *X, const gsl_matrix *Y,
             gsl_matrix *U_hat, gsl_matrix *E_hat, gsl_matrix *OmegaU,
             gsl_matrix *OmegaE, gsl_matrix *UltVehiY, gsl_matrix *UltVehiBX,
             gsl_matrix *UltVehiU, gsl_matrix *UltVehiE, gsl_matrix *V_g,
             gsl_matrix *V_e, gsl_matrix *B) {
  if (func_name != 'R' && func_name != 'L' && func_name != 'r' &&
      func_name != 'l') {
    cout << "func_name only takes 'R' or 'L': 'R' for log-restricted "
         << "likelihood, 'L' for log-likelihood." << endl;
    return 0.0;
  }

  size_t n_size = eval->size, c_size = X->size1, d_size = Y->size1;
  size_t dc_size = d_size * c_size;

  gsl_matrix *XXt = gsl_matrix_alloc(c_size, c_size);
  gsl_matrix *XXti = gsl_matrix_alloc(c_size, c_size);
  gsl_vector *D_l = gsl_vector_alloc(d_size);
  gsl_matrix *UltVeh = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *UltVehi = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *UltVehiB = gsl_matrix_alloc(d_size, c_size);
  gsl_matrix *Qi = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix *Sigma_uu = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *Sigma_ee = gsl_matrix_alloc(d_size, d_size);
  gsl_vector *xHiy = gsl_vector_alloc(dc_size);
  gsl_permutation *pmt = gsl_permutation_alloc(c_size);

  double logl_const = 0.0, logl_old = 0.0, logl_new = 0.0;
  double logdet_Q, logdet_Ve;
  int sig;

  // Calculate |XXt| and (XXt)^{-1}.
  gsl_blas_dsyrk(CblasUpper, CblasNoTrans, 1.0, X, 0.0, XXt);
  for (size_t i = 0; i < c_size; ++i) {
    for (size_t j = 0; j < i; ++j) {
      gsl_matrix_set(XXt, i, j, gsl_matrix_get(XXt, j, i));
    }
  }

  LUDecomp(XXt, pmt, &sig);
  LUInvert(XXt, pmt, XXti);

  // Calculate the constant for logl.
  if (func_name == 'R' || func_name == 'r') {
    logl_const =
        -0.5 * (double)(n_size - c_size) * (double)d_size * safe_log(2.0 * M_PI) +
        0.5 * (double)d_size * LULndet(XXt);
  } else {
    logl_const = -0.5 * (double)n_size * (double)d_size * safe_log(2.0 * M_PI);
  }

  // Start EM.
  for (size_t t = 0; t < max_iter; t++) {
    logdet_Ve = EigenProc(V_g, V_e, D_l, UltVeh, UltVehi);

    logdet_Q = CalcQi(eval, D_l, X, Qi);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, UltVehi, Y, 0.0, UltVehiY);
    CalcXHiY(eval, D_l, X, UltVehiY, xHiy);

    // Calculate log likelihood/restricted likelihood value, and
    // terminate if change is small.
    logl_new = logl_const + MphCalcLogL(eval, xHiy, D_l, UltVehiY, Qi) -
               0.5 * (double)n_size * logdet_Ve;
    if (func_name == 'R' || func_name == 'r') {
      logl_new += -0.5 * (logdet_Q - (double)c_size * logdet_Ve);
    }
    if (t != 0 && abs(logl_new - logl_old) < max_prec) {
      break;
    }
    logl_old = logl_new;

    CalcOmega(eval, D_l, OmegaU, OmegaE);

    // Update UltVehiB, UltVehiU.
    if (func_name == 'R' || func_name == 'r') {
      UpdateRL_B(xHiy, Qi, UltVehiB);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, UltVehiB, X, 0.0,
                     UltVehiBX);
    } else if (t == 0) {
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, UltVehi, B, 0.0,
                     UltVehiB);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, UltVehiB, X, 0.0,
                     UltVehiBX);
    }

    UpdateU(OmegaE, UltVehiY, UltVehiBX, UltVehiU);

    if (func_name == 'L' || func_name == 'l') {

      // UltVehiBX is destroyed here.
      UpdateL_B(X, XXti, UltVehiY, UltVehiU, UltVehiBX, UltVehiB);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, UltVehiB, X, 0.0,
                     UltVehiBX);
    }

    UpdateE(UltVehiY, UltVehiBX, UltVehiU, UltVehiE);

    // Calculate U_hat, E_hat and B.
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, UltVeh, UltVehiU, 0.0, U_hat);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, UltVeh, UltVehiE, 0.0, E_hat);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, UltVeh, UltVehiB, 0.0, B);

    // Calculate Sigma_uu and Sigma_ee.
    CalcSigma(func_name, eval, D_l, X, OmegaU, OmegaE, UltVeh, Qi, Sigma_uu,
              Sigma_ee);

    // Update V_g and V_e.
    UpdateV(eval, U_hat, E_hat, Sigma_uu, Sigma_ee, V_g, V_e);
  }

  gsl_matrix_free(XXt);
  gsl_matrix_free(XXti);
  gsl_vector_free(D_l);
  gsl_matrix_free(UltVeh);
  gsl_matrix_free(UltVehi);
  gsl_matrix_free(UltVehiB);
  gsl_matrix_free(Qi);
  gsl_matrix_free(Sigma_uu);
  gsl_matrix_free(Sigma_ee);
  gsl_vector_free(xHiy);
  gsl_permutation_free(pmt);

  return logl_new;
}

// Calculate p-value, beta (d by 1 vector) and V(beta).
double MphCalcP(const gsl_vector *eval, const gsl_vector *x_vec,
                const gsl_matrix *W, const gsl_matrix *Y, const gsl_matrix *V_g,
                const gsl_matrix *V_e, gsl_matrix *UltVehiY, gsl_vector *beta,
                gsl_matrix *Vbeta) {
  size_t n_size = eval->size, c_size = W->size1, d_size = V_g->size1;
  size_t dc_size = d_size * c_size;
  double delta, dl, d, d1, d2, dy, dx, dw; //  logdet_Ve, logdet_Q, p_value;

  gsl_vector *D_l = gsl_vector_alloc(d_size);
  gsl_matrix *UltVeh = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *UltVehi = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *Qi = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix *WHix = gsl_matrix_alloc(dc_size, d_size);
  gsl_matrix *QiWHix = gsl_matrix_alloc(dc_size, d_size);

  gsl_matrix *xPx = gsl_matrix_alloc(d_size, d_size);
  gsl_vector *xPy = gsl_vector_alloc(d_size);
  gsl_vector *WHiy = gsl_vector_alloc(dc_size);

  gsl_matrix_set_zero(xPx);
  gsl_matrix_set_zero(WHix);
  gsl_vector_set_zero(xPy);
  gsl_vector_set_zero(WHiy);

  // Eigen decomposition and calculate log|Ve|.
  // double logdet_Ve = EigenProc(V_g, V_e, D_l, UltVeh, UltVehi);
  EigenProc(V_g, V_e, D_l, UltVeh, UltVehi);

  // Calculate Qi and log|Q|.
  // double logdet_Q = CalcQi(eval, D_l, W, Qi);
  CalcQi(eval, D_l, W, Qi);

  // Calculate UltVehiY.
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, UltVehi, Y, 0.0, UltVehiY);

  // Calculate WHix, WHiy, xHiy, xHix.
  for (size_t i = 0; i < d_size; i++) {
    dl = gsl_vector_get(D_l, i);

    d1 = 0.0;
    d2 = 0.0;
    for (size_t k = 0; k < n_size; k++) {
      delta = gsl_vector_get(eval, k);
      dx = gsl_vector_get(x_vec, k);
      dy = gsl_matrix_get(UltVehiY, i, k);

      d1 += dx * dy / (delta * dl + 1.0);
      d2 += dx * dx / (delta * dl + 1.0);
    }
    gsl_vector_set(xPy, i, d1);
    gsl_matrix_set(xPx, i, i, d2);

    for (size_t j = 0; j < c_size; j++) {
      d1 = 0.0;
      d2 = 0.0;
      for (size_t k = 0; k < n_size; k++) {
        delta = gsl_vector_get(eval, k);
        dx = gsl_vector_get(x_vec, k);
        dw = gsl_matrix_get(W, j, k);
        dy = gsl_matrix_get(UltVehiY, i, k);

        d1 += dx * dw / (delta * dl + 1.0);
        d2 += dy * dw / (delta * dl + 1.0);
      }
      gsl_matrix_set(WHix, j * d_size + i, i, d1);
      gsl_vector_set(WHiy, j * d_size + i, d2);
    }
  }

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Qi, WHix, 0.0, QiWHix);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, -1.0, WHix, QiWHix, 1.0, xPx);
  gsl_blas_dgemv(CblasTrans, -1.0, QiWHix, WHiy, 1.0, xPy);

  // Calculate V(beta) and beta.
  int sig;
  gsl_permutation *pmt = gsl_permutation_alloc(d_size);
  LUDecomp(xPx, pmt, &sig);
  LUSolve(xPx, pmt, xPy, D_l);
  LUInvert(xPx, pmt, Vbeta);

  // Need to multiply UltVehi on both sides or one side.
  gsl_blas_dgemv(CblasTrans, 1.0, UltVeh, D_l, 0.0, beta);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Vbeta, UltVeh, 0.0, xPx);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, UltVeh, xPx, 0.0, Vbeta);

  // Calculate test statistic and p value.
  gsl_blas_ddot(D_l, xPy, &d);

  double p_value = gsl_cdf_chisq_Q(d, (double)d_size);

  gsl_vector_free(D_l);
  gsl_matrix_free(UltVeh);
  gsl_matrix_free(UltVehi);
  gsl_matrix_free(Qi);
  gsl_matrix_free(WHix);
  gsl_matrix_free(QiWHix);

  gsl_matrix_free(xPx);
  gsl_vector_free(xPy);
  gsl_vector_free(WHiy);

  gsl_permutation_free(pmt);

  return p_value;
}

// Calculate B and its standard error (which is a matrix of the same
// dimension as B).
void MphCalcBeta(const gsl_vector *eval, const gsl_matrix *W,
                 const gsl_matrix *Y, const gsl_matrix *V_g,
                 const gsl_matrix *V_e, gsl_matrix *UltVehiY, gsl_matrix *B,
                 gsl_matrix *se_B) {
  size_t n_size = eval->size, c_size = W->size1, d_size = V_g->size1;
  size_t dc_size = d_size * c_size;
  double delta, dl, d, dy, dw; // , logdet_Ve, logdet_Q;

  gsl_vector *D_l = gsl_vector_alloc(d_size);
  gsl_matrix *UltVeh = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *UltVehi = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *Qi = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix *Qi_temp = gsl_matrix_alloc(dc_size, dc_size);
  gsl_vector *WHiy = gsl_vector_alloc(dc_size);
  gsl_vector *QiWHiy = gsl_vector_alloc(dc_size);
  gsl_vector *beta = gsl_vector_alloc(dc_size);
  gsl_matrix *Vbeta = gsl_matrix_alloc(dc_size, dc_size);

  gsl_vector_set_zero(WHiy);

  // Eigen decomposition and calculate log|Ve|.
  // double logdet_Ve = EigenProc(V_g, V_e, D_l, UltVeh, UltVehi);
  EigenProc(V_g, V_e, D_l, UltVeh, UltVehi);

  // Calculate Qi and log|Q|.
  // double logdet_Q = CalcQi(eval, D_l, W, Qi);
  CalcQi(eval, D_l, W, Qi);

  // Calculate UltVehiY.
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, UltVehi, Y, 0.0, UltVehiY);

  // Calculate WHiy.
  for (size_t i = 0; i < d_size; i++) {
    dl = gsl_vector_get(D_l, i);

    for (size_t j = 0; j < c_size; j++) {
      d = 0.0;
      for (size_t k = 0; k < n_size; k++) {
        delta = gsl_vector_get(eval, k);
        dw = gsl_matrix_get(W, j, k);
        dy = gsl_matrix_get(UltVehiY, i, k);

        d += dy * dw / (delta * dl + 1.0);
      }
      gsl_vector_set(WHiy, j * d_size + i, d);
    }
  }

  gsl_blas_dgemv(CblasNoTrans, 1.0, Qi, WHiy, 0.0, QiWHiy);

  // Need to multiply I_c\otimes UltVehi on both sides or one side.
  for (size_t i = 0; i < c_size; i++) {
    gsl_vector_view QiWHiy_sub =
        gsl_vector_subvector(QiWHiy, i * d_size, d_size);
    gsl_vector_view beta_sub = gsl_vector_subvector(beta, i * d_size, d_size);
    gsl_blas_dgemv(CblasTrans, 1.0, UltVeh, &QiWHiy_sub.vector, 0.0,
                   &beta_sub.vector);

    for (size_t j = 0; j < c_size; j++) {
      gsl_matrix_view Qi_sub =
          gsl_matrix_submatrix(Qi, i * d_size, j * d_size, d_size, d_size);
      gsl_matrix_view Qitemp_sub =
          gsl_matrix_submatrix(Qi_temp, i * d_size, j * d_size, d_size, d_size);
      gsl_matrix_view Vbeta_sub =
          gsl_matrix_submatrix(Vbeta, i * d_size, j * d_size, d_size, d_size);

      if (j < i) {
        gsl_matrix_view Vbeta_sym =
            gsl_matrix_submatrix(Vbeta, j * d_size, i * d_size, d_size, d_size);
        gsl_matrix_transpose_memcpy(&Vbeta_sub.matrix, &Vbeta_sym.matrix);
      } else {
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Qi_sub.matrix, UltVeh,
                       0.0, &Qitemp_sub.matrix);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, UltVeh,
                       &Qitemp_sub.matrix, 0.0, &Vbeta_sub.matrix);
      }
    }
  }

  // Copy beta to B, and Vbeta to se_B.
  for (size_t j = 0; j < B->size2; j++) {
    for (size_t i = 0; i < B->size1; i++) {
      gsl_matrix_set(B, i, j, gsl_vector_get(beta, j * d_size + i));
      gsl_matrix_set(se_B, i, j, safe_sqrt(gsl_matrix_get(Vbeta, j * d_size + i,
                                                     j * d_size + i)));
    }
  }

  // Free matrices.
  gsl_vector_free(D_l);
  gsl_matrix_free(UltVeh);
  gsl_matrix_free(UltVehi);
  gsl_matrix_free(Qi);
  gsl_matrix_free(Qi_temp);
  gsl_vector_free(WHiy);
  gsl_vector_free(QiWHiy);
  gsl_vector_free(beta);
  gsl_matrix_free(Vbeta);

  return;
}

// Below are functions for Newton-Raphson's algorithm.

// Calculate all Hi and return logdet_H=\sum_{k=1}^{n}log|H_k|
// and calculate Qi and return logdet_Q
// and calculate yPy.
void CalcHiQi(const gsl_vector *eval, const gsl_matrix *X,
              const gsl_matrix *V_g, const gsl_matrix *V_e, gsl_matrix *Hi_all,
              gsl_matrix *Qi, double &logdet_H, double &logdet_Q) {
  gsl_matrix_set_zero(Hi_all);
  gsl_matrix_set_zero(Qi);
  logdet_H = 0.0;
  logdet_Q = 0.0;

  size_t n_size = eval->size, c_size = X->size1, d_size = V_g->size1;
  double logdet_Ve = 0.0, delta, dl, d;

  gsl_matrix *mat_dd = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *UltVeh = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *UltVehi = gsl_matrix_alloc(d_size, d_size);
  gsl_vector *D_l = gsl_vector_alloc(d_size);

  // Calculate D_l, UltVeh and UltVehi.
  logdet_Ve = EigenProc(V_g, V_e, D_l, UltVeh, UltVehi);

  // Calculate each Hi and log|H_k|.
  logdet_H = (double)n_size * logdet_Ve;
  for (size_t k = 0; k < n_size; k++) {
    delta = gsl_vector_get(eval, k);

    gsl_matrix_memcpy(mat_dd, UltVehi);
    for (size_t i = 0; i < d_size; i++) {
      dl = gsl_vector_get(D_l, i);
      d = delta * dl + 1.0;

      gsl_vector_view mat_row = gsl_matrix_row(mat_dd, i);
      gsl_vector_scale(&mat_row.vector, 1.0 / d); // @@

      logdet_H += safe_log(d);
    }

    gsl_matrix_view Hi_k =
        gsl_matrix_submatrix(Hi_all, 0, k * d_size, d_size, d_size);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, UltVehi, mat_dd, 0.0,
                   &Hi_k.matrix);
  }

  // Calculate Qi, and multiply I\o times UtVeh on both side and
  // calculate logdet_Q, don't forget to substract
  // c_size*logdet_Ve.
  logdet_Q = CalcQi(eval, D_l, X, Qi) - (double)c_size * logdet_Ve;

  for (size_t i = 0; i < c_size; i++) {
    for (size_t j = 0; j < c_size; j++) {
      gsl_matrix_view Qi_sub =
          gsl_matrix_submatrix(Qi, i * d_size, j * d_size, d_size, d_size);
      if (j < i) {
        gsl_matrix_view Qi_sym =
            gsl_matrix_submatrix(Qi, j * d_size, i * d_size, d_size, d_size);
        gsl_matrix_transpose_memcpy(&Qi_sub.matrix, &Qi_sym.matrix);
      } else {
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Qi_sub.matrix, UltVeh,
                       0.0, mat_dd);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, UltVeh, mat_dd, 0.0,
                       &Qi_sub.matrix);
      }
    }
  }

  // Free memory.
  gsl_matrix_free(mat_dd);
  gsl_matrix_free(UltVeh);
  gsl_matrix_free(UltVehi);
  gsl_vector_free(D_l);

  return;
}

// Calculate all Hiy.
void Calc_Hiy_all(const gsl_matrix *Y, const gsl_matrix *Hi_all,
                  gsl_matrix *Hiy_all) {
  gsl_matrix_set_zero(Hiy_all);

  size_t n_size = Y->size2, d_size = Y->size1;

  for (size_t k = 0; k < n_size; k++) {
    gsl_matrix_const_view Hi_k =
        gsl_matrix_const_submatrix(Hi_all, 0, k * d_size, d_size, d_size);
    gsl_vector_const_view y_k = gsl_matrix_const_column(Y, k);
    gsl_vector_view Hiy_k = gsl_matrix_column(Hiy_all, k);

    gsl_blas_dgemv(CblasNoTrans, 1.0, &Hi_k.matrix, &y_k.vector, 0.0,
                   &Hiy_k.vector);
  }

  return;
}

// Calculate all xHi.
void Calc_xHi_all(const gsl_matrix *X, const gsl_matrix *Hi_all,
                  gsl_matrix *xHi_all) {
  gsl_matrix_set_zero(xHi_all);

  size_t n_size = X->size2, c_size = X->size1, d_size = Hi_all->size1;

  double d;

  for (size_t k = 0; k < n_size; k++) {
    gsl_matrix_const_view Hi_k =
        gsl_matrix_const_submatrix(Hi_all, 0, k * d_size, d_size, d_size);

    for (size_t i = 0; i < c_size; i++) {
      d = gsl_matrix_get(X, i, k);
      gsl_matrix_view xHi_sub =
          gsl_matrix_submatrix(xHi_all, i * d_size, k * d_size, d_size, d_size);
      gsl_matrix_memcpy(&xHi_sub.matrix, &Hi_k.matrix);
      gsl_matrix_scale(&xHi_sub.matrix, d);
    }
  }

  return;
}

// Calculate scalar yHiy.
double Calc_yHiy(const gsl_matrix *Y, const gsl_matrix *Hiy_all) {
  double yHiy = 0.0, d;
  size_t n_size = Y->size2;

  for (size_t k = 0; k < n_size; k++) {
    gsl_vector_const_view y_k = gsl_matrix_const_column(Y, k);
    gsl_vector_const_view Hiy_k = gsl_matrix_const_column(Hiy_all, k);

    gsl_blas_ddot(&Hiy_k.vector, &y_k.vector, &d);
    yHiy += d;
  }

  return yHiy;
}

// Calculate the vector xHiy.
void Calc_xHiy(const gsl_matrix *Y, const gsl_matrix *xHi, gsl_vector *xHiy) {
  gsl_vector_set_zero(xHiy);

  size_t n_size = Y->size2, d_size = Y->size1, dc_size = xHi->size1;

  for (size_t k = 0; k < n_size; k++) {
    gsl_matrix_const_view xHi_k =
        gsl_matrix_const_submatrix(xHi, 0, k * d_size, dc_size, d_size);
    gsl_vector_const_view y_k = gsl_matrix_const_column(Y, k);

    gsl_blas_dgemv(CblasNoTrans, 1.0, &xHi_k.matrix, &y_k.vector, 1.0, xHiy);
  }

  return;
}

// 0<=i,j<d_size
size_t GetIndex(const size_t i, const size_t j, const size_t d_size) {
  if (i >= d_size || j >= d_size) {
    cout << "error in GetIndex." << endl;
    return 0;
  }

  size_t s, l;
  if (j < i) {
    s = j;
    l = i;
  } else {
    s = i;
    l = j;
  }

  return (2 * d_size - s + 1) * s / 2 + l - s;
}

void Calc_yHiDHiy(const gsl_vector *eval, const gsl_matrix *Hiy, const size_t i,
                  const size_t j, double &yHiDHiy_g, double &yHiDHiy_e) {
  yHiDHiy_g = 0.0;
  yHiDHiy_e = 0.0;

  size_t n_size = eval->size;

  double delta, d1, d2;

  for (size_t k = 0; k < n_size; k++) {
    delta = gsl_vector_get(eval, k);
    d1 = gsl_matrix_get(Hiy, i, k);
    d2 = gsl_matrix_get(Hiy, j, k);

    if (i == j) {
      yHiDHiy_g += delta * d1 * d2;
      yHiDHiy_e += d1 * d2;
    } else {
      yHiDHiy_g += delta * d1 * d2 * 2.0;
      yHiDHiy_e += d1 * d2 * 2.0;
    }
  }

  return;
}

void Calc_xHiDHiy(const gsl_vector *eval, const gsl_matrix *xHi,
                  const gsl_matrix *Hiy, const size_t i, const size_t j,
                  gsl_vector *xHiDHiy_g, gsl_vector *xHiDHiy_e) {
  gsl_vector_set_zero(xHiDHiy_g);
  gsl_vector_set_zero(xHiDHiy_e);

  size_t n_size = eval->size, d_size = Hiy->size1;

  double delta, d;

  for (size_t k = 0; k < n_size; k++) {
    delta = gsl_vector_get(eval, k);

    gsl_vector_const_view xHi_col_i =
        gsl_matrix_const_column(xHi, k * d_size + i);
    d = gsl_matrix_get(Hiy, j, k);

    gsl_blas_daxpy(d * delta, &xHi_col_i.vector, xHiDHiy_g);
    gsl_blas_daxpy(d, &xHi_col_i.vector, xHiDHiy_e);

    if (i != j) {
      gsl_vector_const_view xHi_col_j =
          gsl_matrix_const_column(xHi, k * d_size + j);
      d = gsl_matrix_get(Hiy, i, k);

      gsl_blas_daxpy(d * delta, &xHi_col_j.vector, xHiDHiy_g);
      gsl_blas_daxpy(d, &xHi_col_j.vector, xHiDHiy_e);
    }
  }

  return;
}

void Calc_xHiDHix(const gsl_vector *eval, const gsl_matrix *xHi, const size_t i,
                  const size_t j, gsl_matrix *xHiDHix_g,
                  gsl_matrix *xHiDHix_e) {
  gsl_matrix_set_zero(xHiDHix_g);
  gsl_matrix_set_zero(xHiDHix_e);

  size_t n_size = eval->size, dc_size = xHi->size1;
  size_t d_size = xHi->size2 / n_size;

  double delta;

  gsl_matrix *mat_dcdc = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix *mat_dcdc_t = gsl_matrix_alloc(dc_size, dc_size);

  for (size_t k = 0; k < n_size; k++) {
    delta = gsl_vector_get(eval, k);

    gsl_vector_const_view xHi_col_i =
        gsl_matrix_const_column(xHi, k * d_size + i);
    gsl_vector_const_view xHi_col_j =
        gsl_matrix_const_column(xHi, k * d_size + j);

    gsl_matrix_set_zero(mat_dcdc);
    gsl_blas_dger(1.0, &xHi_col_i.vector, &xHi_col_j.vector, mat_dcdc);

    gsl_matrix_transpose_memcpy(mat_dcdc_t, mat_dcdc);

    gsl_matrix_add(xHiDHix_e, mat_dcdc);

    gsl_matrix_scale(mat_dcdc, delta);
    gsl_matrix_add(xHiDHix_g, mat_dcdc);

    if (i != j) {
      gsl_matrix_add(xHiDHix_e, mat_dcdc_t);

      gsl_matrix_scale(mat_dcdc_t, delta);
      gsl_matrix_add(xHiDHix_g, mat_dcdc_t);
    }
  }

  gsl_matrix_free(mat_dcdc);
  gsl_matrix_free(mat_dcdc_t);

  return;
}

void Calc_yHiDHiDHiy(const gsl_vector *eval, const gsl_matrix *Hi,
                     const gsl_matrix *Hiy, const size_t i1, const size_t j1,
                     const size_t i2, const size_t j2, double &yHiDHiDHiy_gg,
                     double &yHiDHiDHiy_ee, double &yHiDHiDHiy_ge) {
  yHiDHiDHiy_gg = 0.0;
  yHiDHiDHiy_ee = 0.0;
  yHiDHiDHiy_ge = 0.0;

  size_t n_size = eval->size, d_size = Hiy->size1;

  double delta, d_Hiy_i1, d_Hiy_j1, d_Hiy_i2, d_Hiy_j2;
  double d_Hi_i1i2, d_Hi_i1j2, d_Hi_j1i2, d_Hi_j1j2;

  for (size_t k = 0; k < n_size; k++) {
    delta = gsl_vector_get(eval, k);

    d_Hiy_i1 = gsl_matrix_get(Hiy, i1, k);
    d_Hiy_j1 = gsl_matrix_get(Hiy, j1, k);
    d_Hiy_i2 = gsl_matrix_get(Hiy, i2, k);
    d_Hiy_j2 = gsl_matrix_get(Hiy, j2, k);

    d_Hi_i1i2 = gsl_matrix_get(Hi, i1, k * d_size + i2);
    d_Hi_i1j2 = gsl_matrix_get(Hi, i1, k * d_size + j2);
    d_Hi_j1i2 = gsl_matrix_get(Hi, j1, k * d_size + i2);
    d_Hi_j1j2 = gsl_matrix_get(Hi, j1, k * d_size + j2);

    if (i1 == j1) {
      yHiDHiDHiy_gg += delta * delta * (d_Hiy_i1 * d_Hi_j1i2 * d_Hiy_j2);
      yHiDHiDHiy_ee += (d_Hiy_i1 * d_Hi_j1i2 * d_Hiy_j2);
      yHiDHiDHiy_ge += delta * (d_Hiy_i1 * d_Hi_j1i2 * d_Hiy_j2);

      if (i2 != j2) {
        yHiDHiDHiy_gg += delta * delta * (d_Hiy_i1 * d_Hi_j1j2 * d_Hiy_i2);
        yHiDHiDHiy_ee += (d_Hiy_i1 * d_Hi_j1j2 * d_Hiy_i2);
        yHiDHiDHiy_ge += delta * (d_Hiy_i1 * d_Hi_j1j2 * d_Hiy_i2);
      }
    } else {
      yHiDHiDHiy_gg += delta * delta * (d_Hiy_i1 * d_Hi_j1i2 * d_Hiy_j2 +
                                        d_Hiy_j1 * d_Hi_i1i2 * d_Hiy_j2);
      yHiDHiDHiy_ee +=
          (d_Hiy_i1 * d_Hi_j1i2 * d_Hiy_j2 + d_Hiy_j1 * d_Hi_i1i2 * d_Hiy_j2);
      yHiDHiDHiy_ge += delta * (d_Hiy_i1 * d_Hi_j1i2 * d_Hiy_j2 +
                                d_Hiy_j1 * d_Hi_i1i2 * d_Hiy_j2);

      if (i2 != j2) {
        yHiDHiDHiy_gg += delta * delta * (d_Hiy_i1 * d_Hi_j1j2 * d_Hiy_i2 +
                                          d_Hiy_j1 * d_Hi_i1j2 * d_Hiy_i2);
        yHiDHiDHiy_ee +=
            (d_Hiy_i1 * d_Hi_j1j2 * d_Hiy_i2 + d_Hiy_j1 * d_Hi_i1j2 * d_Hiy_i2);
        yHiDHiDHiy_ge += delta * (d_Hiy_i1 * d_Hi_j1j2 * d_Hiy_i2 +
                                  d_Hiy_j1 * d_Hi_i1j2 * d_Hiy_i2);
      }
    }
  }

  return;
}

void Calc_xHiDHiDHiy(const gsl_vector *eval, const gsl_matrix *Hi,
                     const gsl_matrix *xHi, const gsl_matrix *Hiy,
                     const size_t i1, const size_t j1, const size_t i2,
                     const size_t j2, gsl_vector *xHiDHiDHiy_gg,
                     gsl_vector *xHiDHiDHiy_ee, gsl_vector *xHiDHiDHiy_ge) {
  gsl_vector_set_zero(xHiDHiDHiy_gg);
  gsl_vector_set_zero(xHiDHiDHiy_ee);
  gsl_vector_set_zero(xHiDHiDHiy_ge);

  size_t n_size = eval->size, d_size = Hiy->size1;

  double delta, d_Hiy_i, d_Hiy_j, d_Hi_i1i2, d_Hi_i1j2;
  double d_Hi_j1i2, d_Hi_j1j2;

  for (size_t k = 0; k < n_size; k++) {
    delta = gsl_vector_get(eval, k);

    gsl_vector_const_view xHi_col_i =
        gsl_matrix_const_column(xHi, k * d_size + i1);
    gsl_vector_const_view xHi_col_j =
        gsl_matrix_const_column(xHi, k * d_size + j1);

    d_Hiy_i = gsl_matrix_get(Hiy, i2, k);
    d_Hiy_j = gsl_matrix_get(Hiy, j2, k);

    d_Hi_i1i2 = gsl_matrix_get(Hi, i1, k * d_size + i2);
    d_Hi_i1j2 = gsl_matrix_get(Hi, i1, k * d_size + j2);
    d_Hi_j1i2 = gsl_matrix_get(Hi, j1, k * d_size + i2);
    d_Hi_j1j2 = gsl_matrix_get(Hi, j1, k * d_size + j2);

    if (i1 == j1) {
      gsl_blas_daxpy(delta * delta * d_Hi_j1i2 * d_Hiy_j, &xHi_col_i.vector,
                     xHiDHiDHiy_gg);
      gsl_blas_daxpy(d_Hi_j1i2 * d_Hiy_j, &xHi_col_i.vector, xHiDHiDHiy_ee);
      gsl_blas_daxpy(delta * d_Hi_j1i2 * d_Hiy_j, &xHi_col_i.vector,
                     xHiDHiDHiy_ge);

      if (i2 != j2) {
        gsl_blas_daxpy(delta * delta * d_Hi_j1j2 * d_Hiy_i, &xHi_col_i.vector,
                       xHiDHiDHiy_gg);
        gsl_blas_daxpy(d_Hi_j1j2 * d_Hiy_i, &xHi_col_i.vector, xHiDHiDHiy_ee);
        gsl_blas_daxpy(delta * d_Hi_j1j2 * d_Hiy_i, &xHi_col_i.vector,
                       xHiDHiDHiy_ge);
      }
    } else {
      gsl_blas_daxpy(delta * delta * d_Hi_j1i2 * d_Hiy_j, &xHi_col_i.vector,
                     xHiDHiDHiy_gg);
      gsl_blas_daxpy(d_Hi_j1i2 * d_Hiy_j, &xHi_col_i.vector, xHiDHiDHiy_ee);
      gsl_blas_daxpy(delta * d_Hi_j1i2 * d_Hiy_j, &xHi_col_i.vector,
                     xHiDHiDHiy_ge);

      gsl_blas_daxpy(delta * delta * d_Hi_i1i2 * d_Hiy_j, &xHi_col_j.vector,
                     xHiDHiDHiy_gg);
      gsl_blas_daxpy(d_Hi_i1i2 * d_Hiy_j, &xHi_col_j.vector, xHiDHiDHiy_ee);
      gsl_blas_daxpy(delta * d_Hi_i1i2 * d_Hiy_j, &xHi_col_j.vector,
                     xHiDHiDHiy_ge);

      if (i2 != j2) {
        gsl_blas_daxpy(delta * delta * d_Hi_j1j2 * d_Hiy_i, &xHi_col_i.vector,
                       xHiDHiDHiy_gg);
        gsl_blas_daxpy(d_Hi_j1j2 * d_Hiy_i, &xHi_col_i.vector, xHiDHiDHiy_ee);
        gsl_blas_daxpy(delta * d_Hi_j1j2 * d_Hiy_i, &xHi_col_i.vector,
                       xHiDHiDHiy_ge);

        gsl_blas_daxpy(delta * delta * d_Hi_i1j2 * d_Hiy_i, &xHi_col_j.vector,
                       xHiDHiDHiy_gg);
        gsl_blas_daxpy(d_Hi_i1j2 * d_Hiy_i, &xHi_col_j.vector, xHiDHiDHiy_ee);
        gsl_blas_daxpy(delta * d_Hi_i1j2 * d_Hiy_i, &xHi_col_j.vector,
                       xHiDHiDHiy_ge);
      }
    }
  }

  return;
}

void Calc_xHiDHiDHix(const gsl_vector *eval, const gsl_matrix *Hi,
                     const gsl_matrix *xHi, const size_t i1, const size_t j1,
                     const size_t i2, const size_t j2,
                     gsl_matrix *xHiDHiDHix_gg, gsl_matrix *xHiDHiDHix_ee,
                     gsl_matrix *xHiDHiDHix_ge) {
  gsl_matrix_set_zero(xHiDHiDHix_gg);
  gsl_matrix_set_zero(xHiDHiDHix_ee);
  gsl_matrix_set_zero(xHiDHiDHix_ge);

  size_t n_size = eval->size, d_size = Hi->size1, dc_size = xHi->size1;

  double delta, d_Hi_i1i2, d_Hi_i1j2, d_Hi_j1i2, d_Hi_j1j2;

  gsl_matrix *mat_dcdc = gsl_matrix_alloc(dc_size, dc_size);

  for (size_t k = 0; k < n_size; k++) {
    delta = gsl_vector_get(eval, k);

    gsl_vector_const_view xHi_col_i1 =
        gsl_matrix_const_column(xHi, k * d_size + i1);
    gsl_vector_const_view xHi_col_j1 =
        gsl_matrix_const_column(xHi, k * d_size + j1);
    gsl_vector_const_view xHi_col_i2 =
        gsl_matrix_const_column(xHi, k * d_size + i2);
    gsl_vector_const_view xHi_col_j2 =
        gsl_matrix_const_column(xHi, k * d_size + j2);

    d_Hi_i1i2 = gsl_matrix_get(Hi, i1, k * d_size + i2);
    d_Hi_i1j2 = gsl_matrix_get(Hi, i1, k * d_size + j2);
    d_Hi_j1i2 = gsl_matrix_get(Hi, j1, k * d_size + i2);
    d_Hi_j1j2 = gsl_matrix_get(Hi, j1, k * d_size + j2);

    if (i1 == j1) {
      gsl_matrix_set_zero(mat_dcdc);
      gsl_blas_dger(d_Hi_j1i2, &xHi_col_i1.vector, &xHi_col_j2.vector,
                    mat_dcdc);

      gsl_matrix_add(xHiDHiDHix_ee, mat_dcdc);
      gsl_matrix_scale(mat_dcdc, delta);
      gsl_matrix_add(xHiDHiDHix_ge, mat_dcdc);
      gsl_matrix_scale(mat_dcdc, delta);
      gsl_matrix_add(xHiDHiDHix_gg, mat_dcdc);

      if (i2 != j2) {
        gsl_matrix_set_zero(mat_dcdc);
        gsl_blas_dger(d_Hi_j1j2, &xHi_col_i1.vector, &xHi_col_i2.vector,
                      mat_dcdc);

        gsl_matrix_add(xHiDHiDHix_ee, mat_dcdc);
        gsl_matrix_scale(mat_dcdc, delta);
        gsl_matrix_add(xHiDHiDHix_ge, mat_dcdc);
        gsl_matrix_scale(mat_dcdc, delta);
        gsl_matrix_add(xHiDHiDHix_gg, mat_dcdc);
      }
    } else {
      gsl_matrix_set_zero(mat_dcdc);
      gsl_blas_dger(d_Hi_j1i2, &xHi_col_i1.vector, &xHi_col_j2.vector,
                    mat_dcdc);

      gsl_matrix_add(xHiDHiDHix_ee, mat_dcdc);
      gsl_matrix_scale(mat_dcdc, delta);
      gsl_matrix_add(xHiDHiDHix_ge, mat_dcdc);
      gsl_matrix_scale(mat_dcdc, delta);
      gsl_matrix_add(xHiDHiDHix_gg, mat_dcdc);

      gsl_matrix_set_zero(mat_dcdc);
      gsl_blas_dger(d_Hi_i1i2, &xHi_col_j1.vector, &xHi_col_j2.vector,
                    mat_dcdc);

      gsl_matrix_add(xHiDHiDHix_ee, mat_dcdc);
      gsl_matrix_scale(mat_dcdc, delta);
      gsl_matrix_add(xHiDHiDHix_ge, mat_dcdc);
      gsl_matrix_scale(mat_dcdc, delta);
      gsl_matrix_add(xHiDHiDHix_gg, mat_dcdc);

      if (i2 != j2) {
        gsl_matrix_set_zero(mat_dcdc);
        gsl_blas_dger(d_Hi_j1j2, &xHi_col_i1.vector, &xHi_col_i2.vector,
                      mat_dcdc);

        gsl_matrix_add(xHiDHiDHix_ee, mat_dcdc);
        gsl_matrix_scale(mat_dcdc, delta);
        gsl_matrix_add(xHiDHiDHix_ge, mat_dcdc);
        gsl_matrix_scale(mat_dcdc, delta);
        gsl_matrix_add(xHiDHiDHix_gg, mat_dcdc);

        gsl_matrix_set_zero(mat_dcdc);
        gsl_blas_dger(d_Hi_i1j2, &xHi_col_j1.vector, &xHi_col_i2.vector,
                      mat_dcdc);

        gsl_matrix_add(xHiDHiDHix_ee, mat_dcdc);
        gsl_matrix_scale(mat_dcdc, delta);
        gsl_matrix_add(xHiDHiDHix_ge, mat_dcdc);
        gsl_matrix_scale(mat_dcdc, delta);
        gsl_matrix_add(xHiDHiDHix_gg, mat_dcdc);
      }
    }
  }

  gsl_matrix_free(mat_dcdc);

  return;
}

void Calc_traceHiD(const gsl_vector *eval, const gsl_matrix *Hi, const size_t i,
                   const size_t j, double &tHiD_g, double &tHiD_e) {
  tHiD_g = 0.0;
  tHiD_e = 0.0;

  size_t n_size = eval->size, d_size = Hi->size1;
  double delta, d;

  for (size_t k = 0; k < n_size; k++) {
    delta = gsl_vector_get(eval, k);
    d = gsl_matrix_get(Hi, j, k * d_size + i);

    if (i == j) {
      tHiD_g += delta * d;
      tHiD_e += d;
    } else {
      tHiD_g += delta * d * 2.0;
      tHiD_e += d * 2.0;
    }
  }

  return;
}

void Calc_traceHiDHiD(const gsl_vector *eval, const gsl_matrix *Hi,
                      const size_t i1, const size_t j1, const size_t i2,
                      const size_t j2, double &tHiDHiD_gg, double &tHiDHiD_ee,
                      double &tHiDHiD_ge) {
  tHiDHiD_gg = 0.0;
  tHiDHiD_ee = 0.0;
  tHiDHiD_ge = 0.0;

  size_t n_size = eval->size, d_size = Hi->size1;
  double delta, d_Hi_i1i2, d_Hi_i1j2, d_Hi_j1i2, d_Hi_j1j2;

  for (size_t k = 0; k < n_size; k++) {
    delta = gsl_vector_get(eval, k);

    d_Hi_i1i2 = gsl_matrix_get(Hi, i1, k * d_size + i2);
    d_Hi_i1j2 = gsl_matrix_get(Hi, i1, k * d_size + j2);
    d_Hi_j1i2 = gsl_matrix_get(Hi, j1, k * d_size + i2);
    d_Hi_j1j2 = gsl_matrix_get(Hi, j1, k * d_size + j2);

    if (i1 == j1) {
      tHiDHiD_gg += delta * delta * d_Hi_i1j2 * d_Hi_j1i2;
      tHiDHiD_ee += d_Hi_i1j2 * d_Hi_j1i2;
      tHiDHiD_ge += delta * d_Hi_i1j2 * d_Hi_j1i2;

      if (i2 != j2) {
        tHiDHiD_gg += delta * delta * d_Hi_i1i2 * d_Hi_j1j2;
        tHiDHiD_ee += d_Hi_i1i2 * d_Hi_j1j2;
        tHiDHiD_ge += delta * d_Hi_i1i2 * d_Hi_j1j2;
      }
    } else {
      tHiDHiD_gg +=
          delta * delta * (d_Hi_i1j2 * d_Hi_j1i2 + d_Hi_j1j2 * d_Hi_i1i2);
      tHiDHiD_ee += (d_Hi_i1j2 * d_Hi_j1i2 + d_Hi_j1j2 * d_Hi_i1i2);
      tHiDHiD_ge += delta * (d_Hi_i1j2 * d_Hi_j1i2 + d_Hi_j1j2 * d_Hi_i1i2);

      if (i2 != j2) {
        tHiDHiD_gg +=
            delta * delta * (d_Hi_i1i2 * d_Hi_j1j2 + d_Hi_j1i2 * d_Hi_i1j2);
        tHiDHiD_ee += (d_Hi_i1i2 * d_Hi_j1j2 + d_Hi_j1i2 * d_Hi_i1j2);
        tHiDHiD_ge += delta * (d_Hi_i1i2 * d_Hi_j1j2 + d_Hi_j1i2 * d_Hi_i1j2);
      }
    }
  }

  return;
}

// trace(PD) = trace((Hi-HixQixHi)D)=trace(HiD) - trace(HixQixHiD)
void Calc_tracePD(const gsl_vector *eval, const gsl_matrix *Qi,
                  const gsl_matrix *Hi, const gsl_matrix *xHiDHix_all_g,
                  const gsl_matrix *xHiDHix_all_e, const size_t i,
                  const size_t j, double &tPD_g, double &tPD_e) {
  size_t dc_size = Qi->size1, d_size = Hi->size1;
  size_t v = GetIndex(i, j, d_size);

  double d;

  // Calculate the first part: trace(HiD).
  Calc_traceHiD(eval, Hi, i, j, tPD_g, tPD_e);

  // Calculate the second part: -trace(HixQixHiD).
  for (size_t k = 0; k < dc_size; k++) {
    gsl_vector_const_view Qi_row = gsl_matrix_const_row(Qi, k);
    gsl_vector_const_view xHiDHix_g_col =
        gsl_matrix_const_column(xHiDHix_all_g, v * dc_size + k);
    gsl_vector_const_view xHiDHix_e_col =
        gsl_matrix_const_column(xHiDHix_all_e, v * dc_size + k);

    gsl_blas_ddot(&Qi_row.vector, &xHiDHix_g_col.vector, &d);
    tPD_g -= d;
    gsl_blas_ddot(&Qi_row.vector, &xHiDHix_e_col.vector, &d);
    tPD_e -= d;
  }

  return;
}

// trace(PDPD) = trace((Hi-HixQixHi)D(Hi-HixQixHi)D)
//             = trace(HiDHiD) - trace(HixQixHiDHiD)
//               - trace(HiDHixQixHiD) + trace(HixQixHiDHixQixHiD)
void Calc_tracePDPD(const gsl_vector *eval, const gsl_matrix *Qi,
                    const gsl_matrix *Hi, const gsl_matrix *xHi,
                    const gsl_matrix *QixHiDHix_all_g,
                    const gsl_matrix *QixHiDHix_all_e,
                    const gsl_matrix *xHiDHiDHix_all_gg,
                    const gsl_matrix *xHiDHiDHix_all_ee,
                    const gsl_matrix *xHiDHiDHix_all_ge, const size_t i1,
                    const size_t j1, const size_t i2, const size_t j2,
                    double &tPDPD_gg, double &tPDPD_ee, double &tPDPD_ge) {
  size_t dc_size = Qi->size1, d_size = Hi->size1;
  size_t v_size = d_size * (d_size + 1) / 2;
  size_t v1 = GetIndex(i1, j1, d_size), v2 = GetIndex(i2, j2, d_size);

  double d;

  // Calculate the first part: trace(HiDHiD).
  Calc_traceHiDHiD(eval, Hi, i1, j1, i2, j2, tPDPD_gg, tPDPD_ee, tPDPD_ge);

  // Calculate the second and third parts:
  // -trace(HixQixHiDHiD) - trace(HiDHixQixHiD)
  for (size_t i = 0; i < dc_size; i++) {
    gsl_vector_const_view Qi_row = gsl_matrix_const_row(Qi, i);
    gsl_vector_const_view xHiDHiDHix_gg_col = gsl_matrix_const_column(
        xHiDHiDHix_all_gg, (v1 * v_size + v2) * dc_size + i);
    gsl_vector_const_view xHiDHiDHix_ee_col = gsl_matrix_const_column(
        xHiDHiDHix_all_ee, (v1 * v_size + v2) * dc_size + i);
    gsl_vector_const_view xHiDHiDHix_ge_col = gsl_matrix_const_column(
        xHiDHiDHix_all_ge, (v1 * v_size + v2) * dc_size + i);

    gsl_blas_ddot(&Qi_row.vector, &xHiDHiDHix_gg_col.vector, &d);
    tPDPD_gg -= d * 2.0;
    gsl_blas_ddot(&Qi_row.vector, &xHiDHiDHix_ee_col.vector, &d);
    tPDPD_ee -= d * 2.0;
    gsl_blas_ddot(&Qi_row.vector, &xHiDHiDHix_ge_col.vector, &d);
    tPDPD_ge -= d * 2.0;
  }

  // Calculate the fourth part: trace(HixQixHiDHixQixHiD).
  for (size_t i = 0; i < dc_size; i++) {

    gsl_vector_const_view QixHiDHix_g_fullrow1 =
        gsl_matrix_const_row(QixHiDHix_all_g, i);
    gsl_vector_const_view QixHiDHix_e_fullrow1 =
        gsl_matrix_const_row(QixHiDHix_all_e, i);
    gsl_vector_const_view QixHiDHix_g_row1 = gsl_vector_const_subvector(
        &QixHiDHix_g_fullrow1.vector, v1 * dc_size, dc_size);
    gsl_vector_const_view QixHiDHix_e_row1 = gsl_vector_const_subvector(
        &QixHiDHix_e_fullrow1.vector, v1 * dc_size, dc_size);

    gsl_vector_const_view QixHiDHix_g_col2 =
        gsl_matrix_const_column(QixHiDHix_all_g, v2 * dc_size + i);
    gsl_vector_const_view QixHiDHix_e_col2 =
        gsl_matrix_const_column(QixHiDHix_all_e, v2 * dc_size + i);

    gsl_blas_ddot(&QixHiDHix_g_row1.vector, &QixHiDHix_g_col2.vector, &d);
    tPDPD_gg += d;
    gsl_blas_ddot(&QixHiDHix_e_row1.vector, &QixHiDHix_e_col2.vector, &d);
    tPDPD_ee += d;
    gsl_blas_ddot(&QixHiDHix_g_row1.vector, &QixHiDHix_e_col2.vector, &d);
    tPDPD_ge += d;
  }

  return;
}

// Calculate (xHiDHiy) for every pair (i,j).
void Calc_xHiDHiy_all(const gsl_vector *eval, const gsl_matrix *xHi,
                      const gsl_matrix *Hiy, gsl_matrix *xHiDHiy_all_g,
                      gsl_matrix *xHiDHiy_all_e) {
  gsl_matrix_set_zero(xHiDHiy_all_g);
  gsl_matrix_set_zero(xHiDHiy_all_e);

  size_t d_size = Hiy->size1;
  size_t v;

  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j < d_size; j++) {
      if (j < i) {
        continue;
      }
      v = GetIndex(i, j, d_size);

      gsl_vector_view xHiDHiy_g = gsl_matrix_column(xHiDHiy_all_g, v);
      gsl_vector_view xHiDHiy_e = gsl_matrix_column(xHiDHiy_all_e, v);

      Calc_xHiDHiy(eval, xHi, Hiy, i, j, &xHiDHiy_g.vector, &xHiDHiy_e.vector);
    }
  }
  return;
}

// Calculate (xHiDHix) for every pair (i,j).
void Calc_xHiDHix_all(const gsl_vector *eval, const gsl_matrix *xHi,
                      gsl_matrix *xHiDHix_all_g, gsl_matrix *xHiDHix_all_e) {
  gsl_matrix_set_zero(xHiDHix_all_g);
  gsl_matrix_set_zero(xHiDHix_all_e);

  size_t d_size = xHi->size2 / eval->size, dc_size = xHi->size1;
  size_t v;

  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j < d_size; j++) {
      if (j < i) {
        continue;
      }
      v = GetIndex(i, j, d_size);

      gsl_matrix_view xHiDHix_g =
          gsl_matrix_submatrix(xHiDHix_all_g, 0, v * dc_size, dc_size, dc_size);
      gsl_matrix_view xHiDHix_e =
          gsl_matrix_submatrix(xHiDHix_all_e, 0, v * dc_size, dc_size, dc_size);

      Calc_xHiDHix(eval, xHi, i, j, &xHiDHix_g.matrix, &xHiDHix_e.matrix);
    }
  }
  return;
}

// Calculate (xHiDHiy) for every pair (i,j).
void Calc_xHiDHiDHiy_all(const size_t v_size, const gsl_vector *eval,
                         const gsl_matrix *Hi, const gsl_matrix *xHi,
                         const gsl_matrix *Hiy, gsl_matrix *xHiDHiDHiy_all_gg,
                         gsl_matrix *xHiDHiDHiy_all_ee,
                         gsl_matrix *xHiDHiDHiy_all_ge) {
  gsl_matrix_set_zero(xHiDHiDHiy_all_gg);
  gsl_matrix_set_zero(xHiDHiDHiy_all_ee);
  gsl_matrix_set_zero(xHiDHiDHiy_all_ge);

  size_t d_size = Hiy->size1;
  size_t v1, v2;

  for (size_t i1 = 0; i1 < d_size; i1++) {
    for (size_t j1 = 0; j1 < d_size; j1++) {
      if (j1 < i1) {
        continue;
      }
      v1 = GetIndex(i1, j1, d_size);

      for (size_t i2 = 0; i2 < d_size; i2++) {
        for (size_t j2 = 0; j2 < d_size; j2++) {
          if (j2 < i2) {
            continue;
          }
          v2 = GetIndex(i2, j2, d_size);

          gsl_vector_view xHiDHiDHiy_gg =
              gsl_matrix_column(xHiDHiDHiy_all_gg, v1 * v_size + v2);
          gsl_vector_view xHiDHiDHiy_ee =
              gsl_matrix_column(xHiDHiDHiy_all_ee, v1 * v_size + v2);
          gsl_vector_view xHiDHiDHiy_ge =
              gsl_matrix_column(xHiDHiDHiy_all_ge, v1 * v_size + v2);

          Calc_xHiDHiDHiy(eval, Hi, xHi, Hiy, i1, j1, i2, j2,
                          &xHiDHiDHiy_gg.vector, &xHiDHiDHiy_ee.vector,
                          &xHiDHiDHiy_ge.vector);
        }
      }
    }
  }
  return;
}

// Calculate (xHiDHix) for every pair (i,j).
void Calc_xHiDHiDHix_all(const size_t v_size, const gsl_vector *eval,
                         const gsl_matrix *Hi, const gsl_matrix *xHi,
                         gsl_matrix *xHiDHiDHix_all_gg,
                         gsl_matrix *xHiDHiDHix_all_ee,
                         gsl_matrix *xHiDHiDHix_all_ge) {
  gsl_matrix_set_zero(xHiDHiDHix_all_gg);
  gsl_matrix_set_zero(xHiDHiDHix_all_ee);
  gsl_matrix_set_zero(xHiDHiDHix_all_ge);

  size_t d_size = xHi->size2 / eval->size, dc_size = xHi->size1;
  size_t v1, v2;

  for (size_t i1 = 0; i1 < d_size; i1++) {
    for (size_t j1 = 0; j1 < d_size; j1++) {
      if (j1 < i1) {
        continue;
      }
      v1 = GetIndex(i1, j1, d_size);

      for (size_t i2 = 0; i2 < d_size; i2++) {
        for (size_t j2 = 0; j2 < d_size; j2++) {
          if (j2 < i2) {
            continue;
          }
          v2 = GetIndex(i2, j2, d_size);

          if (v2 < v1) {
            continue;
          }

          gsl_matrix_view xHiDHiDHix_gg1 = gsl_matrix_submatrix(
              xHiDHiDHix_all_gg, 0, (v1 * v_size + v2) * dc_size, dc_size,
              dc_size);
          gsl_matrix_view xHiDHiDHix_ee1 = gsl_matrix_submatrix(
              xHiDHiDHix_all_ee, 0, (v1 * v_size + v2) * dc_size, dc_size,
              dc_size);
          gsl_matrix_view xHiDHiDHix_ge1 = gsl_matrix_submatrix(
              xHiDHiDHix_all_ge, 0, (v1 * v_size + v2) * dc_size, dc_size,
              dc_size);

          Calc_xHiDHiDHix(eval, Hi, xHi, i1, j1, i2, j2, &xHiDHiDHix_gg1.matrix,
                          &xHiDHiDHix_ee1.matrix, &xHiDHiDHix_ge1.matrix);

          if (v2 != v1) {
            gsl_matrix_view xHiDHiDHix_gg2 = gsl_matrix_submatrix(
                xHiDHiDHix_all_gg, 0, (v2 * v_size + v1) * dc_size, dc_size,
                dc_size);
            gsl_matrix_view xHiDHiDHix_ee2 = gsl_matrix_submatrix(
                xHiDHiDHix_all_ee, 0, (v2 * v_size + v1) * dc_size, dc_size,
                dc_size);
            gsl_matrix_view xHiDHiDHix_ge2 = gsl_matrix_submatrix(
                xHiDHiDHix_all_ge, 0, (v2 * v_size + v1) * dc_size, dc_size,
                dc_size);

            gsl_matrix_memcpy(&xHiDHiDHix_gg2.matrix, &xHiDHiDHix_gg1.matrix);
            gsl_matrix_memcpy(&xHiDHiDHix_ee2.matrix, &xHiDHiDHix_ee1.matrix);
            gsl_matrix_memcpy(&xHiDHiDHix_ge2.matrix, &xHiDHiDHix_ge1.matrix);
          }
        }
      }
    }
  }

  return;
}

// Calculate (xHiDHix)Qi(xHiy) for every pair (i,j).
void Calc_xHiDHixQixHiy_all(const gsl_matrix *xHiDHix_all_g,
                            const gsl_matrix *xHiDHix_all_e,
                            const gsl_vector *QixHiy,
                            gsl_matrix *xHiDHixQixHiy_all_g,
                            gsl_matrix *xHiDHixQixHiy_all_e) {
  size_t dc_size = xHiDHix_all_g->size1;
  size_t v_size = xHiDHix_all_g->size2 / dc_size;

  for (size_t i = 0; i < v_size; i++) {
    gsl_matrix_const_view xHiDHix_g = gsl_matrix_const_submatrix(
        xHiDHix_all_g, 0, i * dc_size, dc_size, dc_size);
    gsl_matrix_const_view xHiDHix_e = gsl_matrix_const_submatrix(
        xHiDHix_all_e, 0, i * dc_size, dc_size, dc_size);

    gsl_vector_view xHiDHixQixHiy_g = gsl_matrix_column(xHiDHixQixHiy_all_g, i);
    gsl_vector_view xHiDHixQixHiy_e = gsl_matrix_column(xHiDHixQixHiy_all_e, i);

    gsl_blas_dgemv(CblasNoTrans, 1.0, &xHiDHix_g.matrix, QixHiy, 0.0,
                   &xHiDHixQixHiy_g.vector);
    gsl_blas_dgemv(CblasNoTrans, 1.0, &xHiDHix_e.matrix, QixHiy, 0.0,
                   &xHiDHixQixHiy_e.vector);
  }

  return;
}

// Calculate Qi(xHiDHiy) and Qi(xHiDHix)Qi(xHiy) for each pair of i,j (i<=j).
void Calc_QiVec_all(const gsl_matrix *Qi, const gsl_matrix *vec_all_g,
                    const gsl_matrix *vec_all_e, gsl_matrix *Qivec_all_g,
                    gsl_matrix *Qivec_all_e) {
  for (size_t i = 0; i < vec_all_g->size2; i++) {
    gsl_vector_const_view vec_g = gsl_matrix_const_column(vec_all_g, i);
    gsl_vector_const_view vec_e = gsl_matrix_const_column(vec_all_e, i);

    gsl_vector_view Qivec_g = gsl_matrix_column(Qivec_all_g, i);
    gsl_vector_view Qivec_e = gsl_matrix_column(Qivec_all_e, i);

    gsl_blas_dgemv(CblasNoTrans, 1.0, Qi, &vec_g.vector, 0.0, &Qivec_g.vector);
    gsl_blas_dgemv(CblasNoTrans, 1.0, Qi, &vec_e.vector, 0.0, &Qivec_e.vector);
  }

  return;
}

// Calculate Qi(xHiDHix) for each pair of i,j (i<=j).
void Calc_QiMat_all(const gsl_matrix *Qi, const gsl_matrix *mat_all_g,
                    const gsl_matrix *mat_all_e, gsl_matrix *Qimat_all_g,
                    gsl_matrix *Qimat_all_e) {
  size_t dc_size = Qi->size1;
  size_t v_size = mat_all_g->size2 / mat_all_g->size1;

  for (size_t i = 0; i < v_size; i++) {
    gsl_matrix_const_view mat_g =
        gsl_matrix_const_submatrix(mat_all_g, 0, i * dc_size, dc_size, dc_size);
    gsl_matrix_const_view mat_e =
        gsl_matrix_const_submatrix(mat_all_e, 0, i * dc_size, dc_size, dc_size);

    gsl_matrix_view Qimat_g =
        gsl_matrix_submatrix(Qimat_all_g, 0, i * dc_size, dc_size, dc_size);
    gsl_matrix_view Qimat_e =
        gsl_matrix_submatrix(Qimat_all_e, 0, i * dc_size, dc_size, dc_size);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Qi, &mat_g.matrix, 0.0,
                   &Qimat_g.matrix);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Qi, &mat_e.matrix, 0.0,
                   &Qimat_e.matrix);
  }

  return;
}

// Calculate yPDPy
// yPDPy = y(Hi-HixQixHi)D(Hi-HixQixHi)y
//       = ytHiDHiy - (yHix)Qi(xHiDHiy) - (yHiDHix)Qi(xHiy)
//         + (yHix)Qi(xHiDHix)Qi(xtHiy)
void Calc_yPDPy(const gsl_vector *eval, const gsl_matrix *Hiy,
                const gsl_vector *QixHiy, const gsl_matrix *xHiDHiy_all_g,
                const gsl_matrix *xHiDHiy_all_e,
                const gsl_matrix *xHiDHixQixHiy_all_g,
                const gsl_matrix *xHiDHixQixHiy_all_e, const size_t i,
                const size_t j, double &yPDPy_g, double &yPDPy_e) {
  size_t d_size = Hiy->size1;
  size_t v = GetIndex(i, j, d_size);

  double d;

  // First part: ytHiDHiy.
  Calc_yHiDHiy(eval, Hiy, i, j, yPDPy_g, yPDPy_e);

  // Second and third parts: -(yHix)Qi(xHiDHiy)-(yHiDHix)Qi(xHiy)
  gsl_vector_const_view xHiDHiy_g = gsl_matrix_const_column(xHiDHiy_all_g, v);
  gsl_vector_const_view xHiDHiy_e = gsl_matrix_const_column(xHiDHiy_all_e, v);

  gsl_blas_ddot(QixHiy, &xHiDHiy_g.vector, &d);
  yPDPy_g -= d * 2.0;
  gsl_blas_ddot(QixHiy, &xHiDHiy_e.vector, &d);
  yPDPy_e -= d * 2.0;

  // Fourth part: +(yHix)Qi(xHiDHix)Qi(xHiy).
  gsl_vector_const_view xHiDHixQixHiy_g =
      gsl_matrix_const_column(xHiDHixQixHiy_all_g, v);
  gsl_vector_const_view xHiDHixQixHiy_e =
      gsl_matrix_const_column(xHiDHixQixHiy_all_e, v);

  gsl_blas_ddot(QixHiy, &xHiDHixQixHiy_g.vector, &d);
  yPDPy_g += d;
  gsl_blas_ddot(QixHiy, &xHiDHixQixHiy_e.vector, &d);
  yPDPy_e += d;

  return;
}

// calculate yPDPDPy = y(Hi-HixQixHi)D(Hi-HixQixHi)D(Hi-HixQixHi)y
//           yPDPDPy = yHiDHiDHiy
//                     - (yHix)Qi(xHiDHiDHiy)-(yHiDHiDHix)Qi(xHiy)
//                     - (yHiDHix)Qi(xHiDHiy)
//                     + (yHix)Qi(xHiDHix)Qi(xHiDHiy)
//                     + (yHiDHix)Qi(xHiDHix)Qi(xHiy)
//                     + (yHix)Qi(xHiDHiDHix)Qi(xHiy)
//                     - (yHix)Qi(xHiDHix)Qi(xHiDHix)Qi(xHiy)
void Calc_yPDPDPy(
    const gsl_vector *eval, const gsl_matrix *Hi, const gsl_matrix *xHi,
    const gsl_matrix *Hiy, const gsl_vector *QixHiy,
    const gsl_matrix *xHiDHiy_all_g, const gsl_matrix *xHiDHiy_all_e,
    const gsl_matrix *QixHiDHiy_all_g, const gsl_matrix *QixHiDHiy_all_e,
    const gsl_matrix *xHiDHixQixHiy_all_g,
    const gsl_matrix *xHiDHixQixHiy_all_e,
    const gsl_matrix *QixHiDHixQixHiy_all_g,
    const gsl_matrix *QixHiDHixQixHiy_all_e,
    const gsl_matrix *xHiDHiDHiy_all_gg, const gsl_matrix *xHiDHiDHiy_all_ee,
    const gsl_matrix *xHiDHiDHiy_all_ge, const gsl_matrix *xHiDHiDHix_all_gg,
    const gsl_matrix *xHiDHiDHix_all_ee, const gsl_matrix *xHiDHiDHix_all_ge,
    const size_t i1, const size_t j1, const size_t i2, const size_t j2,
    double &yPDPDPy_gg, double &yPDPDPy_ee, double &yPDPDPy_ge) {
  size_t d_size = Hi->size1, dc_size = xHi->size1;
  size_t v1 = GetIndex(i1, j1, d_size), v2 = GetIndex(i2, j2, d_size);
  size_t v_size = d_size * (d_size + 1) / 2;

  double d;

  gsl_vector *xHiDHiDHixQixHiy = gsl_vector_alloc(dc_size);

  // First part: yHiDHiDHiy.
  Calc_yHiDHiDHiy(eval, Hi, Hiy, i1, j1, i2, j2, yPDPDPy_gg, yPDPDPy_ee,
                  yPDPDPy_ge);

  // Second and third parts:
  // -(yHix)Qi(xHiDHiDHiy) - (yHiDHiDHix)Qi(xHiy).
  gsl_vector_const_view xHiDHiDHiy_gg1 =
      gsl_matrix_const_column(xHiDHiDHiy_all_gg, v1 * v_size + v2);
  gsl_vector_const_view xHiDHiDHiy_ee1 =
      gsl_matrix_const_column(xHiDHiDHiy_all_ee, v1 * v_size + v2);
  gsl_vector_const_view xHiDHiDHiy_ge1 =
      gsl_matrix_const_column(xHiDHiDHiy_all_ge, v1 * v_size + v2);

  gsl_vector_const_view xHiDHiDHiy_gg2 =
      gsl_matrix_const_column(xHiDHiDHiy_all_gg, v2 * v_size + v1);
  gsl_vector_const_view xHiDHiDHiy_ee2 =
      gsl_matrix_const_column(xHiDHiDHiy_all_ee, v2 * v_size + v1);
  gsl_vector_const_view xHiDHiDHiy_ge2 =
      gsl_matrix_const_column(xHiDHiDHiy_all_ge, v2 * v_size + v1);

  gsl_blas_ddot(QixHiy, &xHiDHiDHiy_gg1.vector, &d);
  yPDPDPy_gg -= d;
  gsl_blas_ddot(QixHiy, &xHiDHiDHiy_ee1.vector, &d);
  yPDPDPy_ee -= d;
  gsl_blas_ddot(QixHiy, &xHiDHiDHiy_ge1.vector, &d);
  yPDPDPy_ge -= d;

  gsl_blas_ddot(QixHiy, &xHiDHiDHiy_gg2.vector, &d);
  yPDPDPy_gg -= d;
  gsl_blas_ddot(QixHiy, &xHiDHiDHiy_ee2.vector, &d);
  yPDPDPy_ee -= d;
  gsl_blas_ddot(QixHiy, &xHiDHiDHiy_ge2.vector, &d);
  yPDPDPy_ge -= d;

  // Fourth part: - (yHiDHix)Qi(xHiDHiy).
  gsl_vector_const_view xHiDHiy_g1 = gsl_matrix_const_column(xHiDHiy_all_g, v1);
  gsl_vector_const_view xHiDHiy_e1 = gsl_matrix_const_column(xHiDHiy_all_e, v1);
  gsl_vector_const_view QixHiDHiy_g2 =
      gsl_matrix_const_column(QixHiDHiy_all_g, v2);
  gsl_vector_const_view QixHiDHiy_e2 =
      gsl_matrix_const_column(QixHiDHiy_all_e, v2);

  gsl_blas_ddot(&xHiDHiy_g1.vector, &QixHiDHiy_g2.vector, &d);
  yPDPDPy_gg -= d;
  gsl_blas_ddot(&xHiDHiy_e1.vector, &QixHiDHiy_e2.vector, &d);
  yPDPDPy_ee -= d;
  gsl_blas_ddot(&xHiDHiy_g1.vector, &QixHiDHiy_e2.vector, &d);
  yPDPDPy_ge -= d;

  // Fifth and sixth parts:
  //   + (yHix)Qi(xHiDHix)Qi(xHiDHiy) +
  //   (yHiDHix)Qi(xHiDHix)Qi(xHiy)
  gsl_vector_const_view QixHiDHiy_g1 =
      gsl_matrix_const_column(QixHiDHiy_all_g, v1);
  gsl_vector_const_view QixHiDHiy_e1 =
      gsl_matrix_const_column(QixHiDHiy_all_e, v1);

  gsl_vector_const_view xHiDHixQixHiy_g1 =
      gsl_matrix_const_column(xHiDHixQixHiy_all_g, v1);
  gsl_vector_const_view xHiDHixQixHiy_e1 =
      gsl_matrix_const_column(xHiDHixQixHiy_all_e, v1);
  gsl_vector_const_view xHiDHixQixHiy_g2 =
      gsl_matrix_const_column(xHiDHixQixHiy_all_g, v2);
  gsl_vector_const_view xHiDHixQixHiy_e2 =
      gsl_matrix_const_column(xHiDHixQixHiy_all_e, v2);

  gsl_blas_ddot(&xHiDHixQixHiy_g1.vector, &QixHiDHiy_g2.vector, &d);
  yPDPDPy_gg += d;
  gsl_blas_ddot(&xHiDHixQixHiy_g2.vector, &QixHiDHiy_g1.vector, &d);
  yPDPDPy_gg += d;

  gsl_blas_ddot(&xHiDHixQixHiy_e1.vector, &QixHiDHiy_e2.vector, &d);
  yPDPDPy_ee += d;
  gsl_blas_ddot(&xHiDHixQixHiy_e2.vector, &QixHiDHiy_e1.vector, &d);
  yPDPDPy_ee += d;

  gsl_blas_ddot(&xHiDHixQixHiy_g1.vector, &QixHiDHiy_e2.vector, &d);
  yPDPDPy_ge += d;
  gsl_blas_ddot(&xHiDHixQixHiy_e2.vector, &QixHiDHiy_g1.vector, &d);
  yPDPDPy_ge += d;

  // Seventh part: + (yHix)Qi(xHiDHiDHix)Qi(xHiy)
  gsl_matrix_const_view xHiDHiDHix_gg = gsl_matrix_const_submatrix(
      xHiDHiDHix_all_gg, 0, (v1 * v_size + v2) * dc_size, dc_size, dc_size);
  gsl_matrix_const_view xHiDHiDHix_ee = gsl_matrix_const_submatrix(
      xHiDHiDHix_all_ee, 0, (v1 * v_size + v2) * dc_size, dc_size, dc_size);
  gsl_matrix_const_view xHiDHiDHix_ge = gsl_matrix_const_submatrix(
      xHiDHiDHix_all_ge, 0, (v1 * v_size + v2) * dc_size, dc_size, dc_size);

  gsl_blas_dgemv(CblasNoTrans, 1.0, &xHiDHiDHix_gg.matrix, QixHiy, 0.0,
                 xHiDHiDHixQixHiy);
  gsl_blas_ddot(xHiDHiDHixQixHiy, QixHiy, &d);
  yPDPDPy_gg += d;
  gsl_blas_dgemv(CblasNoTrans, 1.0, &xHiDHiDHix_ee.matrix, QixHiy, 0.0,
                 xHiDHiDHixQixHiy);
  gsl_blas_ddot(xHiDHiDHixQixHiy, QixHiy, &d);
  yPDPDPy_ee += d;
  gsl_blas_dgemv(CblasNoTrans, 1.0, &xHiDHiDHix_ge.matrix, QixHiy, 0.0,
                 xHiDHiDHixQixHiy);
  gsl_blas_ddot(xHiDHiDHixQixHiy, QixHiy, &d);
  yPDPDPy_ge += d;

  // Eighth part: - (yHix)Qi(xHiDHix)Qi(xHiDHix)Qi(xHiy).
  gsl_vector_const_view QixHiDHixQixHiy_g1 =
      gsl_matrix_const_column(QixHiDHixQixHiy_all_g, v1);
  gsl_vector_const_view QixHiDHixQixHiy_e1 =
      gsl_matrix_const_column(QixHiDHixQixHiy_all_e, v1);

  gsl_blas_ddot(&QixHiDHixQixHiy_g1.vector, &xHiDHixQixHiy_g2.vector, &d);
  yPDPDPy_gg -= d;
  gsl_blas_ddot(&QixHiDHixQixHiy_e1.vector, &xHiDHixQixHiy_e2.vector, &d);
  yPDPDPy_ee -= d;
  gsl_blas_ddot(&QixHiDHixQixHiy_g1.vector, &xHiDHixQixHiy_e2.vector, &d);
  yPDPDPy_ge -= d;

  // Free memory.
  gsl_vector_free(xHiDHiDHixQixHiy);

  return;
}

// Calculate Edgeworth correctation factors for small samples notation
// and method follows Thomas J. Rothenberg, Econometirca 1984; 52 (4)
// M=xHiDHix
void CalcCRT(const gsl_matrix *Hessian_inv, const gsl_matrix *Qi,
             const gsl_matrix *QixHiDHix_all_g,
             const gsl_matrix *QixHiDHix_all_e,
             const gsl_matrix *xHiDHiDHix_all_gg,
             const gsl_matrix *xHiDHiDHix_all_ee,
             const gsl_matrix *xHiDHiDHix_all_ge, const size_t d_size,
             double &crt_a, double &crt_b, double &crt_c) {
  crt_a = 0.0;
  crt_b = 0.0;
  crt_c = 0.0;

  size_t dc_size = Qi->size1, v_size = Hessian_inv->size1 / 2;
  size_t c_size = dc_size / d_size;
  double h_gg, h_ge, h_ee, d, B = 0.0, C = 0.0, D = 0.0;
  double trCg1, trCe1, trCg2, trCe2, trB_gg, trB_ge, trB_ee;
  double trCC_gg, trCC_ge, trCC_ee, trD_gg = 0.0, trD_ge = 0.0, trD_ee = 0.0;

  gsl_matrix *QiMQi_g1 = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix *QiMQi_e1 = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix *QiMQi_g2 = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix *QiMQi_e2 = gsl_matrix_alloc(dc_size, dc_size);

  gsl_matrix *QiMQisQisi_g1 = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *QiMQisQisi_e1 = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *QiMQisQisi_g2 = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *QiMQisQisi_e2 = gsl_matrix_alloc(d_size, d_size);

  gsl_matrix *QiMQiMQi_gg = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix *QiMQiMQi_ge = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix *QiMQiMQi_ee = gsl_matrix_alloc(dc_size, dc_size);

  gsl_matrix *QiMMQi_gg = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix *QiMMQi_ge = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix *QiMMQi_ee = gsl_matrix_alloc(dc_size, dc_size);

  gsl_matrix *Qi_si = gsl_matrix_alloc(d_size, d_size);

  gsl_matrix *M_dd = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *M_dcdc = gsl_matrix_alloc(dc_size, dc_size);

  // Invert Qi_sub to Qi_si.
  gsl_matrix *Qi_sub = gsl_matrix_alloc(d_size, d_size);

  gsl_matrix_const_view Qi_s = gsl_matrix_const_submatrix(
      Qi, (c_size - 1) * d_size, (c_size - 1) * d_size, d_size, d_size);

  int sig;
  gsl_permutation *pmt = gsl_permutation_alloc(d_size);

  gsl_matrix_memcpy(Qi_sub, &Qi_s.matrix);
  LUDecomp(Qi_sub, pmt, &sig);
  LUInvert(Qi_sub, pmt, Qi_si);

  gsl_permutation_free(pmt);
  gsl_matrix_free(Qi_sub);

  // Calculate correction factors.
  for (size_t v1 = 0; v1 < v_size; v1++) {

    // Calculate Qi(xHiDHix)Qi, and subpart of it.
    gsl_matrix_const_view QiM_g1 = gsl_matrix_const_submatrix(
        QixHiDHix_all_g, 0, v1 * dc_size, dc_size, dc_size);
    gsl_matrix_const_view QiM_e1 = gsl_matrix_const_submatrix(
        QixHiDHix_all_e, 0, v1 * dc_size, dc_size, dc_size);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiM_g1.matrix, Qi, 0.0,
                   QiMQi_g1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiM_e1.matrix, Qi, 0.0,
                   QiMQi_e1);

    gsl_matrix_view QiMQi_g1_s = gsl_matrix_submatrix(
        QiMQi_g1, (c_size - 1) * d_size, (c_size - 1) * d_size, d_size, d_size);
    gsl_matrix_view QiMQi_e1_s = gsl_matrix_submatrix(
        QiMQi_e1, (c_size - 1) * d_size, (c_size - 1) * d_size, d_size, d_size);

    // Calculate trCg1 and trCe1.
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiMQi_g1_s.matrix, Qi_si,
                   0.0, QiMQisQisi_g1);
    trCg1 = 0.0;
    for (size_t k = 0; k < d_size; k++) {
      trCg1 -= gsl_matrix_get(QiMQisQisi_g1, k, k);
    }

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiMQi_e1_s.matrix, Qi_si,
                   0.0, QiMQisQisi_e1);
    trCe1 = 0.0;
    for (size_t k = 0; k < d_size; k++) {
      trCe1 -= gsl_matrix_get(QiMQisQisi_e1, k, k);
    }

    for (size_t v2 = 0; v2 < v_size; v2++) {
      if (v2 < v1) {
        continue;
      }

      // Calculate Qi(xHiDHix)Qi, and subpart of it.
      gsl_matrix_const_view QiM_g2 = gsl_matrix_const_submatrix(
          QixHiDHix_all_g, 0, v2 * dc_size, dc_size, dc_size);
      gsl_matrix_const_view QiM_e2 = gsl_matrix_const_submatrix(
          QixHiDHix_all_e, 0, v2 * dc_size, dc_size, dc_size);

      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiM_g2.matrix, Qi, 0.0,
                     QiMQi_g2);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiM_e2.matrix, Qi, 0.0,
                     QiMQi_e2);

      gsl_matrix_view QiMQi_g2_s =
          gsl_matrix_submatrix(QiMQi_g2, (c_size - 1) * d_size,
                               (c_size - 1) * d_size, d_size, d_size);
      gsl_matrix_view QiMQi_e2_s =
          gsl_matrix_submatrix(QiMQi_e2, (c_size - 1) * d_size,
                               (c_size - 1) * d_size, d_size, d_size);

      // Calculate trCg2 and trCe2.
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiMQi_g2_s.matrix, Qi_si,
                     0.0, QiMQisQisi_g2);
      trCg2 = 0.0;
      for (size_t k = 0; k < d_size; k++) {
        trCg2 -= gsl_matrix_get(QiMQisQisi_g2, k, k);
      }

      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiMQi_e2_s.matrix, Qi_si,
                     0.0, QiMQisQisi_e2);
      trCe2 = 0.0;
      for (size_t k = 0; k < d_size; k++) {
        trCe2 -= gsl_matrix_get(QiMQisQisi_e2, k, k);
      }

      // Calculate trCC_gg, trCC_ge, trCC_ee.
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, QiMQisQisi_g1,
                     QiMQisQisi_g2, 0.0, M_dd);
      trCC_gg = 0.0;
      for (size_t k = 0; k < d_size; k++) {
        trCC_gg += gsl_matrix_get(M_dd, k, k);
      }

      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, QiMQisQisi_g1,
                     QiMQisQisi_e2, 0.0, M_dd);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, QiMQisQisi_e1,
                     QiMQisQisi_g2, 1.0, M_dd);
      trCC_ge = 0.0;
      for (size_t k = 0; k < d_size; k++) {
        trCC_ge += gsl_matrix_get(M_dd, k, k);
      }

      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, QiMQisQisi_e1,
                     QiMQisQisi_e2, 0.0, M_dd);
      trCC_ee = 0.0;
      for (size_t k = 0; k < d_size; k++) {
        trCC_ee += gsl_matrix_get(M_dd, k, k);
      }

      // Calculate Qi(xHiDHix)Qi(xHiDHix)Qi, and subpart of it.
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiM_g1.matrix, QiMQi_g2,
                     0.0, QiMQiMQi_gg);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiM_g1.matrix, QiMQi_e2,
                     0.0, QiMQiMQi_ge);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiM_e1.matrix, QiMQi_g2,
                     1.0, QiMQiMQi_ge);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiM_e1.matrix, QiMQi_e2,
                     0.0, QiMQiMQi_ee);

      gsl_matrix_view QiMQiMQi_gg_s =
          gsl_matrix_submatrix(QiMQiMQi_gg, (c_size - 1) * d_size,
                               (c_size - 1) * d_size, d_size, d_size);
      gsl_matrix_view QiMQiMQi_ge_s =
          gsl_matrix_submatrix(QiMQiMQi_ge, (c_size - 1) * d_size,
                               (c_size - 1) * d_size, d_size, d_size);
      gsl_matrix_view QiMQiMQi_ee_s =
          gsl_matrix_submatrix(QiMQiMQi_ee, (c_size - 1) * d_size,
                               (c_size - 1) * d_size, d_size, d_size);

      // and part of trB_gg, trB_ge, trB_ee.
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiMQiMQi_gg_s.matrix,
                     Qi_si, 0.0, M_dd);
      trB_gg = 0.0;
      for (size_t k = 0; k < d_size; k++) {
        d = gsl_matrix_get(M_dd, k, k);
        trB_gg -= d;
      }

      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiMQiMQi_ge_s.matrix,
                     Qi_si, 0.0, M_dd);
      trB_ge = 0.0;
      for (size_t k = 0; k < d_size; k++) {
        d = gsl_matrix_get(M_dd, k, k);
        trB_ge -= d;
      }

      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiMQiMQi_ee_s.matrix,
                     Qi_si, 0.0, M_dd);
      trB_ee = 0.0;
      for (size_t k = 0; k < d_size; k++) {
        d = gsl_matrix_get(M_dd, k, k);
        trB_ee -= d;
      }

      // Calculate Qi(xHiDHiDHix)Qi, and subpart of it.
      gsl_matrix_const_view MM_gg = gsl_matrix_const_submatrix(
          xHiDHiDHix_all_gg, 0, (v1 * v_size + v2) * dc_size, dc_size, dc_size);
      gsl_matrix_const_view MM_ge = gsl_matrix_const_submatrix(
          xHiDHiDHix_all_ge, 0, (v1 * v_size + v2) * dc_size, dc_size, dc_size);
      gsl_matrix_const_view MM_ee = gsl_matrix_const_submatrix(
          xHiDHiDHix_all_ee, 0, (v1 * v_size + v2) * dc_size, dc_size, dc_size);

      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Qi, &MM_gg.matrix, 0.0,
                     M_dcdc);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, M_dcdc, Qi, 0.0,
                     QiMMQi_gg);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Qi, &MM_ge.matrix, 0.0,
                     M_dcdc);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, M_dcdc, Qi, 0.0,
                     QiMMQi_ge);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Qi, &MM_ee.matrix, 0.0,
                     M_dcdc);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, M_dcdc, Qi, 0.0,
                     QiMMQi_ee);

      gsl_matrix_view QiMMQi_gg_s =
          gsl_matrix_submatrix(QiMMQi_gg, (c_size - 1) * d_size,
                               (c_size - 1) * d_size, d_size, d_size);
      gsl_matrix_view QiMMQi_ge_s =
          gsl_matrix_submatrix(QiMMQi_ge, (c_size - 1) * d_size,
                               (c_size - 1) * d_size, d_size, d_size);
      gsl_matrix_view QiMMQi_ee_s =
          gsl_matrix_submatrix(QiMMQi_ee, (c_size - 1) * d_size,
                               (c_size - 1) * d_size, d_size, d_size);

      // Calculate the other part of trB_gg, trB_ge, trB_ee.
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiMMQi_gg_s.matrix,
                     Qi_si, 0.0, M_dd);
      for (size_t k = 0; k < d_size; k++) {
        trB_gg += gsl_matrix_get(M_dd, k, k);
      }
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiMMQi_ge_s.matrix,
                     Qi_si, 0.0, M_dd);
      for (size_t k = 0; k < d_size; k++) {
        trB_ge += 2.0 * gsl_matrix_get(M_dd, k, k);
      }
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &QiMMQi_ee_s.matrix,
                     Qi_si, 0.0, M_dd);
      for (size_t k = 0; k < d_size; k++) {
        trB_ee += gsl_matrix_get(M_dd, k, k);
      }

      // Calculate trD_gg, trD_ge, trD_ee.
      trD_gg = 2.0 * trB_gg;
      trD_ge = 2.0 * trB_ge;
      trD_ee = 2.0 * trB_ee;

      // calculate B, C and D
      h_gg = -1.0 * gsl_matrix_get(Hessian_inv, v1, v2);
      h_ge = -1.0 * gsl_matrix_get(Hessian_inv, v1, v2 + v_size);
      h_ee = -1.0 * gsl_matrix_get(Hessian_inv, v1 + v_size, v2 + v_size);

      B += h_gg * trB_gg + h_ge * trB_ge + h_ee * trB_ee;
      C += h_gg * (trCC_gg + 0.5 * trCg1 * trCg2) +
           h_ge * (trCC_ge + 0.5 * trCg1 * trCe2 + 0.5 * trCe1 * trCg2) +
           h_ee * (trCC_ee + 0.5 * trCe1 * trCe2);
      D += h_gg * (trCC_gg + 0.5 * trD_gg) + h_ge * (trCC_ge + 0.5 * trD_ge) +
           h_ee * (trCC_ee + 0.5 * trD_ee);

      if (v1 != v2) {
        B += h_gg * trB_gg + h_ge * trB_ge + h_ee * trB_ee;
        C += h_gg * (trCC_gg + 0.5 * trCg1 * trCg2) +
             h_ge * (trCC_ge + 0.5 * trCg1 * trCe2 + 0.5 * trCe1 * trCg2) +
             h_ee * (trCC_ee + 0.5 * trCe1 * trCe2);
        D += h_gg * (trCC_gg + 0.5 * trD_gg) + h_ge * (trCC_ge + 0.5 * trD_ge) +
             h_ee * (trCC_ee + 0.5 * trD_ee);
      }
    }
  }

  // Calculate a, b, c from B C D.
  crt_a = 2.0 * D - C;
  crt_b = 2.0 * B;
  crt_c = C;

  // Free matrix memory.
  gsl_matrix_free(QiMQi_g1);
  gsl_matrix_free(QiMQi_e1);
  gsl_matrix_free(QiMQi_g2);
  gsl_matrix_free(QiMQi_e2);

  gsl_matrix_free(QiMQisQisi_g1);
  gsl_matrix_free(QiMQisQisi_e1);
  gsl_matrix_free(QiMQisQisi_g2);
  gsl_matrix_free(QiMQisQisi_e2);

  gsl_matrix_free(QiMQiMQi_gg);
  gsl_matrix_free(QiMQiMQi_ge);
  gsl_matrix_free(QiMQiMQi_ee);

  gsl_matrix_free(QiMMQi_gg);
  gsl_matrix_free(QiMMQi_ge);
  gsl_matrix_free(QiMMQi_ee);

  gsl_matrix_free(Qi_si);

  gsl_matrix_free(M_dd);
  gsl_matrix_free(M_dcdc);

  return;
}

// Calculate first-order and second-order derivatives.
void CalcDev(const char func_name, const gsl_vector *eval, const gsl_matrix *Qi,
             const gsl_matrix *Hi, const gsl_matrix *xHi, const gsl_matrix *Hiy,
             const gsl_vector *QixHiy, gsl_vector *gradient,
             gsl_matrix *Hessian_inv, double &crt_a, double &crt_b,
             double &crt_c) {
  if (func_name != 'R' && func_name != 'L' && func_name != 'r' &&
      func_name != 'l') {
    cout << "func_name only takes 'R' or 'L': 'R' for "
         << "log-restricted likelihood, 'L' for log-likelihood." << endl;
    return;
  }

  size_t dc_size = Qi->size1, d_size = Hi->size1;
  size_t c_size = dc_size / d_size;
  size_t v_size = d_size * (d_size + 1) / 2;
  size_t v1, v2;
  double dev1_g, dev1_e, dev2_gg, dev2_ee, dev2_ge;

  gsl_matrix *Hessian = gsl_matrix_alloc(v_size * 2, v_size * 2);

  gsl_matrix *xHiDHiy_all_g = gsl_matrix_alloc(dc_size, v_size);
  gsl_matrix *xHiDHiy_all_e = gsl_matrix_alloc(dc_size, v_size);
  gsl_matrix *xHiDHix_all_g = gsl_matrix_alloc(dc_size, v_size * dc_size);
  gsl_matrix *xHiDHix_all_e = gsl_matrix_alloc(dc_size, v_size * dc_size);
  gsl_matrix *xHiDHixQixHiy_all_g = gsl_matrix_alloc(dc_size, v_size);
  gsl_matrix *xHiDHixQixHiy_all_e = gsl_matrix_alloc(dc_size, v_size);

  gsl_matrix *QixHiDHiy_all_g = gsl_matrix_alloc(dc_size, v_size);
  gsl_matrix *QixHiDHiy_all_e = gsl_matrix_alloc(dc_size, v_size);
  gsl_matrix *QixHiDHix_all_g = gsl_matrix_alloc(dc_size, v_size * dc_size);
  gsl_matrix *QixHiDHix_all_e = gsl_matrix_alloc(dc_size, v_size * dc_size);
  gsl_matrix *QixHiDHixQixHiy_all_g = gsl_matrix_alloc(dc_size, v_size);
  gsl_matrix *QixHiDHixQixHiy_all_e = gsl_matrix_alloc(dc_size, v_size);

  gsl_matrix *xHiDHiDHiy_all_gg = gsl_matrix_alloc(dc_size, v_size * v_size);
  gsl_matrix *xHiDHiDHiy_all_ee = gsl_matrix_alloc(dc_size, v_size * v_size);
  gsl_matrix *xHiDHiDHiy_all_ge = gsl_matrix_alloc(dc_size, v_size * v_size);
  gsl_matrix *xHiDHiDHix_all_gg =
      gsl_matrix_alloc(dc_size, v_size * v_size * dc_size);
  gsl_matrix *xHiDHiDHix_all_ee =
      gsl_matrix_alloc(dc_size, v_size * v_size * dc_size);
  gsl_matrix *xHiDHiDHix_all_ge =
      gsl_matrix_alloc(dc_size, v_size * v_size * dc_size);

  // Calculate xHiDHiy_all, xHiDHix_all and xHiDHixQixHiy_all.
  Calc_xHiDHiy_all(eval, xHi, Hiy, xHiDHiy_all_g, xHiDHiy_all_e);
  Calc_xHiDHix_all(eval, xHi, xHiDHix_all_g, xHiDHix_all_e);
  Calc_xHiDHixQixHiy_all(xHiDHix_all_g, xHiDHix_all_e, QixHiy,
                         xHiDHixQixHiy_all_g, xHiDHixQixHiy_all_e);

  Calc_xHiDHiDHiy_all(v_size, eval, Hi, xHi, Hiy, xHiDHiDHiy_all_gg,
                      xHiDHiDHiy_all_ee, xHiDHiDHiy_all_ge);
  Calc_xHiDHiDHix_all(v_size, eval, Hi, xHi, xHiDHiDHix_all_gg,
                      xHiDHiDHix_all_ee, xHiDHiDHix_all_ge);

  // Calculate QixHiDHiy_all, QixHiDHix_all and QixHiDHixQixHiy_all.
  Calc_QiVec_all(Qi, xHiDHiy_all_g, xHiDHiy_all_e, QixHiDHiy_all_g,
                 QixHiDHiy_all_e);
  Calc_QiVec_all(Qi, xHiDHixQixHiy_all_g, xHiDHixQixHiy_all_e,
                 QixHiDHixQixHiy_all_g, QixHiDHixQixHiy_all_e);
  Calc_QiMat_all(Qi, xHiDHix_all_g, xHiDHix_all_e, QixHiDHix_all_g,
                 QixHiDHix_all_e);

  double tHiD_g, tHiD_e, tPD_g, tPD_e, tHiDHiD_gg, tHiDHiD_ee;
  double tHiDHiD_ge, tPDPD_gg, tPDPD_ee, tPDPD_ge;
  double yPDPy_g, yPDPy_e, yPDPDPy_gg, yPDPDPy_ee, yPDPDPy_ge;

  // Calculate gradient and Hessian for Vg.
  for (size_t i1 = 0; i1 < d_size; i1++) {
    for (size_t j1 = 0; j1 < d_size; j1++) {
      if (j1 < i1) {
        continue;
      }
      v1 = GetIndex(i1, j1, d_size);

      Calc_yPDPy(eval, Hiy, QixHiy, xHiDHiy_all_g, xHiDHiy_all_e,
                 xHiDHixQixHiy_all_g, xHiDHixQixHiy_all_e, i1, j1, yPDPy_g,
                 yPDPy_e);

      if (func_name == 'R' || func_name == 'r') {
        Calc_tracePD(eval, Qi, Hi, xHiDHix_all_g, xHiDHix_all_e, i1, j1, tPD_g,
                     tPD_e);

        dev1_g = -0.5 * tPD_g + 0.5 * yPDPy_g;
        dev1_e = -0.5 * tPD_e + 0.5 * yPDPy_e;
      } else {
        Calc_traceHiD(eval, Hi, i1, j1, tHiD_g, tHiD_e);

        dev1_g = -0.5 * tHiD_g + 0.5 * yPDPy_g;
        dev1_e = -0.5 * tHiD_e + 0.5 * yPDPy_e;
      }

      gsl_vector_set(gradient, v1, dev1_g);
      gsl_vector_set(gradient, v1 + v_size, dev1_e);

      for (size_t i2 = 0; i2 < d_size; i2++) {
        for (size_t j2 = 0; j2 < d_size; j2++) {
          if (j2 < i2) {
            continue;
          }
          v2 = GetIndex(i2, j2, d_size);

          if (v2 < v1) {
            continue;
          }

          Calc_yPDPDPy(eval, Hi, xHi, Hiy, QixHiy, xHiDHiy_all_g, xHiDHiy_all_e,
                       QixHiDHiy_all_g, QixHiDHiy_all_e, xHiDHixQixHiy_all_g,
                       xHiDHixQixHiy_all_e, QixHiDHixQixHiy_all_g,
                       QixHiDHixQixHiy_all_e, xHiDHiDHiy_all_gg,
                       xHiDHiDHiy_all_ee, xHiDHiDHiy_all_ge, xHiDHiDHix_all_gg,
                       xHiDHiDHix_all_ee, xHiDHiDHix_all_ge, i1, j1, i2, j2,
                       yPDPDPy_gg, yPDPDPy_ee, yPDPDPy_ge);

          // AI for REML.
          if (func_name == 'R' || func_name == 'r') {
            Calc_tracePDPD(eval, Qi, Hi, xHi, QixHiDHix_all_g, QixHiDHix_all_e,
                           xHiDHiDHix_all_gg, xHiDHiDHix_all_ee,
                           xHiDHiDHix_all_ge, i1, j1, i2, j2, tPDPD_gg,
                           tPDPD_ee, tPDPD_ge);

            dev2_gg = 0.5 * tPDPD_gg - yPDPDPy_gg;
            dev2_ee = 0.5 * tPDPD_ee - yPDPDPy_ee;
            dev2_ge = 0.5 * tPDPD_ge - yPDPDPy_ge;
          } else {
            Calc_traceHiDHiD(eval, Hi, i1, j1, i2, j2, tHiDHiD_gg, tHiDHiD_ee,
                             tHiDHiD_ge);

            dev2_gg = 0.5 * tHiDHiD_gg - yPDPDPy_gg;
            dev2_ee = 0.5 * tHiDHiD_ee - yPDPDPy_ee;
            dev2_ge = 0.5 * tHiDHiD_ge - yPDPDPy_ge;
          }

          // Set up Hessian.
          gsl_matrix_set(Hessian, v1, v2, dev2_gg);
          gsl_matrix_set(Hessian, v1 + v_size, v2 + v_size, dev2_ee);
          gsl_matrix_set(Hessian, v1, v2 + v_size, dev2_ge);
          gsl_matrix_set(Hessian, v2 + v_size, v1, dev2_ge);

          if (v1 != v2) {
            gsl_matrix_set(Hessian, v2, v1, dev2_gg);
            gsl_matrix_set(Hessian, v2 + v_size, v1 + v_size, dev2_ee);
            gsl_matrix_set(Hessian, v2, v1 + v_size, dev2_ge);
            gsl_matrix_set(Hessian, v1 + v_size, v2, dev2_ge);
          }
        }
      }
    }
  }

  // Invert Hessian.
  int sig;
  gsl_permutation *pmt = gsl_permutation_alloc(v_size * 2);

  LUDecomp(Hessian, pmt, &sig);
  LUInvert(Hessian, pmt, Hessian_inv);

  gsl_permutation_free(pmt);
  gsl_matrix_free(Hessian);

  // Calculate Edgeworth correction factors after inverting
  // Hessian.
  if (c_size > 1) {
    CalcCRT(Hessian_inv, Qi, QixHiDHix_all_g, QixHiDHix_all_e,
            xHiDHiDHix_all_gg, xHiDHiDHix_all_ee, xHiDHiDHix_all_ge, d_size,
            crt_a, crt_b, crt_c);
  } else {
    crt_a = 0.0;
    crt_b = 0.0;
    crt_c = 0.0;
  }

  gsl_matrix_free(xHiDHiy_all_g);
  gsl_matrix_free(xHiDHiy_all_e);
  gsl_matrix_free(xHiDHix_all_g);
  gsl_matrix_free(xHiDHix_all_e);
  gsl_matrix_free(xHiDHixQixHiy_all_g);
  gsl_matrix_free(xHiDHixQixHiy_all_e);

  gsl_matrix_free(QixHiDHiy_all_g);
  gsl_matrix_free(QixHiDHiy_all_e);
  gsl_matrix_free(QixHiDHix_all_g);
  gsl_matrix_free(QixHiDHix_all_e);
  gsl_matrix_free(QixHiDHixQixHiy_all_g);
  gsl_matrix_free(QixHiDHixQixHiy_all_e);

  gsl_matrix_free(xHiDHiDHiy_all_gg);
  gsl_matrix_free(xHiDHiDHiy_all_ee);
  gsl_matrix_free(xHiDHiDHiy_all_ge);
  gsl_matrix_free(xHiDHiDHix_all_gg);
  gsl_matrix_free(xHiDHiDHix_all_ee);
  gsl_matrix_free(xHiDHiDHix_all_ge);

  return;
}

// Update Vg, Ve.
void UpdateVgVe(const gsl_matrix *Hessian_inv, const gsl_vector *gradient,
                const double step_scale, gsl_matrix *V_g, gsl_matrix *V_e) {
  size_t v_size = gradient->size / 2, d_size = V_g->size1;
  size_t v;

  gsl_vector *vec_v = gsl_vector_alloc(v_size * 2);

  double d;

  // Vectorize Vg and Ve.
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j < d_size; j++) {
      if (j < i) {
        continue;
      }
      v = GetIndex(i, j, d_size);

      d = gsl_matrix_get(V_g, i, j);
      gsl_vector_set(vec_v, v, d);

      d = gsl_matrix_get(V_e, i, j);
      gsl_vector_set(vec_v, v + v_size, d);
    }
  }

  gsl_blas_dgemv(CblasNoTrans, -1.0 * step_scale, Hessian_inv, gradient, 1.0,
                 vec_v);

  // Save Vg and Ve.
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j < d_size; j++) {
      if (j < i) {
        continue;
      }
      v = GetIndex(i, j, d_size);

      d = gsl_vector_get(vec_v, v);
      gsl_matrix_set(V_g, i, j, d);
      gsl_matrix_set(V_g, j, i, d);

      d = gsl_vector_get(vec_v, v + v_size);
      gsl_matrix_set(V_e, i, j, d);
      gsl_matrix_set(V_e, j, i, d);
    }
  }

  gsl_vector_free(vec_v);

  return;
}

double MphNR(const char func_name, const size_t max_iter, const double max_prec,
             const gsl_vector *eval, const gsl_matrix *X, const gsl_matrix *Y,
             gsl_matrix *Hi_all, gsl_matrix *xHi_all, gsl_matrix *Hiy_all,
             gsl_matrix *V_g, gsl_matrix *V_e, gsl_matrix *Hessian_inv,
             double &crt_a, double &crt_b, double &crt_c) {
  if (func_name != 'R' && func_name != 'L' && func_name != 'r' &&
      func_name != 'l') {
    cout << "func_name only takes 'R' or 'L': 'R' for log-restricted "
         << "likelihood, 'L' for log-likelihood." << endl;
    return 0.0;
  }
  size_t n_size = eval->size, c_size = X->size1, d_size = Y->size1;
  size_t dc_size = d_size * c_size;
  size_t v_size = d_size * (d_size + 1) / 2;

  double logdet_H, logdet_Q, yPy, logl_const;
  double logl_old = 0.0, logl_new = 0.0, step_scale;
  int sig;
  size_t step_iter, flag_pd;

  gsl_matrix *Vg_save = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *Ve_save = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *V_temp = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *U_temp = gsl_matrix_alloc(d_size, d_size);
  gsl_vector *D_temp = gsl_vector_alloc(d_size);
  gsl_vector *xHiy = gsl_vector_alloc(dc_size);
  gsl_vector *QixHiy = gsl_vector_alloc(dc_size);
  gsl_matrix *Qi = gsl_matrix_alloc(dc_size, dc_size);
  gsl_matrix *XXt = gsl_matrix_alloc(c_size, c_size);

  gsl_vector *gradient = gsl_vector_alloc(v_size * 2);

  // Calculate |XXt| and (XXt)^{-1}.
  gsl_blas_dsyrk(CblasUpper, CblasNoTrans, 1.0, X, 0.0, XXt);
  for (size_t i = 0; i < c_size; ++i) {
    for (size_t j = 0; j < i; ++j) {
      gsl_matrix_set(XXt, i, j, gsl_matrix_get(XXt, j, i));
    }
  }

  gsl_permutation *pmt = gsl_permutation_alloc(c_size);
  LUDecomp(XXt, pmt, &sig);
  gsl_permutation_free(pmt);

  // Calculate the constant for logl.
  if (func_name == 'R' || func_name == 'r') {
    logl_const =
        -0.5 * (double)(n_size - c_size) * (double)d_size * safe_log(2.0 * M_PI) +
        0.5 * (double)d_size * LULndet(XXt);
  } else {
    logl_const = -0.5 * (double)n_size * (double)d_size * safe_log(2.0 * M_PI);
  }

  // Optimization iterations.
  for (size_t t = 0; t < max_iter; t++) {
    gsl_matrix_memcpy(Vg_save, V_g);
    gsl_matrix_memcpy(Ve_save, V_e);

    step_scale = 1.0;
    step_iter = 0;
    do {
      gsl_matrix_memcpy(V_g, Vg_save);
      gsl_matrix_memcpy(V_e, Ve_save);

      // Update Vg, Ve, and invert Hessian.
      if (t != 0) {
        UpdateVgVe(Hessian_inv, gradient, step_scale, V_g, V_e);
      }

      // Check if both Vg and Ve are positive definite.
      flag_pd = 1;
      gsl_matrix_memcpy(V_temp, V_e);
      EigenDecomp(V_temp, U_temp, D_temp, 0);
      for (size_t i = 0; i < d_size; i++) {
        if (gsl_vector_get(D_temp, i) <= 0) {
          flag_pd = 0;
        }
      }
      gsl_matrix_memcpy(V_temp, V_g);
      EigenDecomp(V_temp, U_temp, D_temp, 0);
      for (size_t i = 0; i < d_size; i++) {
        if (gsl_vector_get(D_temp, i) <= 0) {
          flag_pd = 0;
        }
      }

      // If flag_pd==1, continue to calculate quantities
      // and logl.
      if (flag_pd == 1) {
        CalcHiQi(eval, X, V_g, V_e, Hi_all, Qi, logdet_H, logdet_Q);
        Calc_Hiy_all(Y, Hi_all, Hiy_all);
        Calc_xHi_all(X, Hi_all, xHi_all);

        // Calculate QixHiy and yPy.
        Calc_xHiy(Y, xHi_all, xHiy);
        gsl_blas_dgemv(CblasNoTrans, 1.0, Qi, xHiy, 0.0, QixHiy);

        gsl_blas_ddot(QixHiy, xHiy, &yPy);
        yPy = Calc_yHiy(Y, Hiy_all) - yPy;

        // Calculate log likelihood/restricted likelihood value.
        if (func_name == 'R' || func_name == 'r') {
          logl_new = logl_const - 0.5 * logdet_H - 0.5 * logdet_Q - 0.5 * yPy;
        } else {
          logl_new = logl_const - 0.5 * logdet_H - 0.5 * yPy;
        }
      }

      step_scale /= 2.0;
      step_iter++;

    } while (
        (flag_pd == 0 || logl_new < logl_old || logl_new - logl_old > 10) &&
        step_iter < 10 && t != 0);

    // Terminate if change is small.
    if (t != 0) {
      if (logl_new < logl_old || flag_pd == 0) {
        gsl_matrix_memcpy(V_g, Vg_save);
        gsl_matrix_memcpy(V_e, Ve_save);
        break;
      }

      if (logl_new - logl_old < max_prec) {
        break;
      }
    }

    logl_old = logl_new;

    CalcDev(func_name, eval, Qi, Hi_all, xHi_all, Hiy_all, QixHiy, gradient,
            Hessian_inv, crt_a, crt_b, crt_c);
  }

  // Mutiply Hessian_inv with -1.0.
  // Now Hessian_inv is the variance matrix.
  gsl_matrix_scale(Hessian_inv, -1.0);

  gsl_matrix_free(Vg_save);
  gsl_matrix_free(Ve_save);
  gsl_matrix_free(V_temp);
  gsl_matrix_free(U_temp);
  gsl_vector_free(D_temp);
  gsl_vector_free(xHiy);
  gsl_vector_free(QixHiy);

  gsl_matrix_free(Qi);
  gsl_matrix_free(XXt);

  gsl_vector_free(gradient);

  return logl_new;
}

// Initialize Vg, Ve and B.
void MphInitial(const size_t em_iter, const double em_prec,
                const size_t nr_iter, const double nr_prec,
                const gsl_vector *eval, const gsl_matrix *X,
                const gsl_matrix *Y, const double l_min, const double l_max,
                const size_t n_region, gsl_matrix *V_g, gsl_matrix *V_e,
                gsl_matrix *B) {

  debug_msg("MphInitial");
  gsl_matrix_set_zero(V_g);
  gsl_matrix_set_zero(V_e);
  gsl_matrix_set_zero(B);

  size_t n_size = eval->size, c_size = X->size1, d_size = Y->size1;
  double a, b, c;
  double lambda, logl, vg, ve;

  // Initialize the diagonal elements of Vg and Ve using univariate
  // LMM and REML estimates.
  gsl_matrix *Xt = gsl_matrix_alloc(n_size, c_size);
  gsl_vector *beta_temp = gsl_vector_alloc(c_size);
  gsl_vector *se_beta_temp = gsl_vector_alloc(c_size);

  gsl_matrix_transpose_memcpy(Xt, X);

  for (size_t i = 0; i < d_size; i++) {
    gsl_vector_const_view Y_row = gsl_matrix_const_row(Y, i);
    CalcLambda('R', eval, Xt, &Y_row.vector, l_min, l_max, n_region, lambda,
               logl);
    CalcLmmVgVeBeta(eval, Xt, &Y_row.vector, lambda, vg, ve, beta_temp,
                    se_beta_temp);

    gsl_matrix_set(V_g, i, i, vg);
    gsl_matrix_set(V_e, i, i, ve);
  }

  gsl_matrix_free(Xt);
  gsl_vector_free(beta_temp);
  gsl_vector_free(se_beta_temp);

  // If number of phenotypes is above four, then obtain the off
  // diagonal elements with two trait models.
  if (d_size > 4) {

    // First obtain good initial values.
    // Large matrices for EM.
    gsl_matrix *U_hat = gsl_matrix_alloc(2, n_size);
    gsl_matrix *E_hat = gsl_matrix_alloc(2, n_size);
    gsl_matrix *OmegaU = gsl_matrix_alloc(2, n_size);
    gsl_matrix *OmegaE = gsl_matrix_alloc(2, n_size);
    gsl_matrix *UltVehiY = gsl_matrix_alloc(2, n_size);
    gsl_matrix *UltVehiBX = gsl_matrix_alloc(2, n_size);
    gsl_matrix *UltVehiU = gsl_matrix_alloc(2, n_size);
    gsl_matrix *UltVehiE = gsl_matrix_alloc(2, n_size);

    // Large matrices for NR. Each dxd block is H_k^{-1}.
    gsl_matrix *Hi_all = gsl_matrix_alloc(2, 2 * n_size);

    // Each column is H_k^{-1}y_k.
    gsl_matrix *Hiy_all = gsl_matrix_alloc(2, n_size);

    // Each dcxdc block is x_k\otimes H_k^{-1}.
    gsl_matrix *xHi_all = gsl_matrix_alloc(2 * c_size, 2 * n_size);
    gsl_matrix *Hessian = gsl_matrix_alloc(6, 6);

    // 2 by n matrix of Y.
    gsl_matrix *Y_sub = gsl_matrix_alloc(2, n_size);
    gsl_matrix *Vg_sub = gsl_matrix_alloc(2, 2);
    gsl_matrix *Ve_sub = gsl_matrix_alloc(2, 2);
    gsl_matrix *B_sub = gsl_matrix_alloc(2, c_size);

    for (size_t i = 0; i < d_size; i++) {
      gsl_vector_view Y_sub1 = gsl_matrix_row(Y_sub, 0);
      gsl_vector_const_view Y_1 = gsl_matrix_const_row(Y, i);
      gsl_vector_memcpy(&Y_sub1.vector, &Y_1.vector);

      for (size_t j = i + 1; j < d_size; j++) {
        gsl_vector_view Y_sub2 = gsl_matrix_row(Y_sub, 1);
        gsl_vector_const_view Y_2 = gsl_matrix_const_row(Y, j);
        gsl_vector_memcpy(&Y_sub2.vector, &Y_2.vector);

        gsl_matrix_set_zero(Vg_sub);
        gsl_matrix_set_zero(Ve_sub);
        gsl_matrix_set(Vg_sub, 0, 0, gsl_matrix_get(V_g, i, i));
        gsl_matrix_set(Ve_sub, 0, 0, gsl_matrix_get(V_e, i, i));
        gsl_matrix_set(Vg_sub, 1, 1, gsl_matrix_get(V_g, j, j));
        gsl_matrix_set(Ve_sub, 1, 1, gsl_matrix_get(V_e, j, j));

        logl = MphEM('R', em_iter, em_prec, eval, X, Y_sub, U_hat, E_hat,
                     OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE,
                     Vg_sub, Ve_sub, B_sub);
        logl = MphNR('R', nr_iter, nr_prec, eval, X, Y_sub, Hi_all, xHi_all,
                     Hiy_all, Vg_sub, Ve_sub, Hessian, a, b, c);

        gsl_matrix_set(V_g, i, j, gsl_matrix_get(Vg_sub, 0, 1));
        gsl_matrix_set(V_g, j, i, gsl_matrix_get(Vg_sub, 0, 1));

        gsl_matrix_set(V_e, i, j, ve = gsl_matrix_get(Ve_sub, 0, 1));
        gsl_matrix_set(V_e, j, i, ve = gsl_matrix_get(Ve_sub, 0, 1));
      }
    }

    // Free matrices.
    gsl_matrix_free(U_hat);
    gsl_matrix_free(E_hat);
    gsl_matrix_free(OmegaU);
    gsl_matrix_free(OmegaE);
    gsl_matrix_free(UltVehiY);
    gsl_matrix_free(UltVehiBX);
    gsl_matrix_free(UltVehiU);
    gsl_matrix_free(UltVehiE);

    gsl_matrix_free(Hi_all);
    gsl_matrix_free(Hiy_all);
    gsl_matrix_free(xHi_all);
    gsl_matrix_free(Hessian);

    gsl_matrix_free(Y_sub);
    gsl_matrix_free(Vg_sub);
    gsl_matrix_free(Ve_sub);
    gsl_matrix_free(B_sub);
  }

  // Calculate B hat using GSL estimate.
  gsl_matrix *UltVehiY = gsl_matrix_alloc(d_size, n_size);

  gsl_vector *D_l = gsl_vector_alloc(d_size);
  gsl_matrix *UltVeh = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *UltVehi = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *Qi = gsl_matrix_alloc(d_size * c_size, d_size * c_size);
  gsl_vector *XHiy = gsl_vector_alloc(d_size * c_size);
  gsl_vector *beta = gsl_vector_alloc(d_size * c_size);

  gsl_vector_set_zero(XHiy);

  double dl, d, delta, dx, dy;

  // Eigen decomposition and calculate log|Ve|.
  // double logdet_Ve = EigenProc(V_g, V_e, D_l, UltVeh, UltVehi);
  EigenProc(V_g, V_e, D_l, UltVeh, UltVehi);

  // Calculate Qi and log|Q|.
  // double logdet_Q = CalcQi(eval, D_l, X, Qi);
  CalcQi(eval, D_l, X, Qi);

  // Calculate UltVehiY.
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, UltVehi, Y, 0.0, UltVehiY);

  // calculate XHiy
  for (size_t i = 0; i < d_size; i++) {
    dl = gsl_vector_get(D_l, i);

    for (size_t j = 0; j < c_size; j++) {
      d = 0.0;
      for (size_t k = 0; k < n_size; k++) {
        delta = gsl_vector_get(eval, k);
        dx = gsl_matrix_get(X, j, k);
        dy = gsl_matrix_get(UltVehiY, i, k);
        d += dy * dx / (delta * dl + 1.0);
      }
      gsl_vector_set(XHiy, j * d_size + i, d);
    }
  }

  gsl_blas_dgemv(CblasNoTrans, 1.0, Qi, XHiy, 0.0, beta);

  // Multiply beta by UltVeh and save to B.
  for (size_t i = 0; i < c_size; i++) {
    gsl_vector_view B_col = gsl_matrix_column(B, i);
    gsl_vector_view beta_sub = gsl_vector_subvector(beta, i * d_size, d_size);
    gsl_blas_dgemv(CblasTrans, 1.0, UltVeh, &beta_sub.vector, 0.0,
                   &B_col.vector);
  }

  // Free memory.
  gsl_matrix_free(UltVehiY);

  gsl_vector_free(D_l);
  gsl_matrix_free(UltVeh);
  gsl_matrix_free(UltVehi);
  gsl_matrix_free(Qi);
  gsl_vector_free(XHiy);
  gsl_vector_free(beta);

  return;
}

// p-value correction
// mode=1 Wald; mode=2 LRT; mode=3 SCORE;
double PCRT(const size_t mode, const size_t d_size, const double p_value,
            const double crt_a, const double crt_b, const double crt_c) {
  double p_crt = 0.0, chisq_crt = 0.0, q = (double)d_size;
  double chisq = gsl_cdf_chisq_Qinv(p_value, (double)d_size);

  if (mode == 1) {
    double a = crt_c / (2.0 * q * (q + 2.0));
    double b = 1.0 + (crt_a + crt_b) / (2.0 * q);
    chisq_crt = (-1.0 * b + safe_sqrt(b * b + 4.0 * a * chisq)) / (2.0 * a);
  } else if (mode == 2) {
    chisq_crt = chisq / (1.0 + crt_a / (2.0 * q));
  } else {
    chisq_crt = chisq;
  }

  p_crt = gsl_cdf_chisq_Q(chisq_crt, (double)d_size);

  return p_crt;
}

void MVLMM::AnalyzeBimbam(const gsl_matrix *U, const gsl_vector *eval,
                          const gsl_matrix *UtW, const gsl_matrix *UtY) {
  debug_msg("entering");
  igzstream infile(file_geno.c_str(), igzstream::in);
  if (!infile) {
    cout << "error reading genotype file:" << file_geno << endl;
    return;
  }

  clock_t time_start = clock();
  time_UtX = 0;
  time_opt = 0;

  string line;
  char *ch_ptr;

  double logl_H0 = 0.0, logl_H1 = 0.0, p_wald = 0, p_lrt = 0, p_score = 0;
  double crt_a, crt_b, crt_c;
  int n_miss, c_phen;
  double geno, x_mean;
  size_t c = 0;
  size_t n_size = UtY->size1, d_size = UtY->size2, c_size = UtW->size2;

  size_t dc_size = d_size * (c_size + 1), v_size = d_size * (d_size + 1) / 2;

  // Create a large matrix.
  size_t msize = LMM_BATCH_SIZE;
  gsl_matrix *Xlarge = gsl_matrix_alloc(U->size1, msize);
  gsl_matrix *UtXlarge = gsl_matrix_alloc(U->size1, msize);
  gsl_matrix_set_zero(Xlarge);

  // Large matrices for EM.
  gsl_matrix *U_hat = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *E_hat = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *OmegaU = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *OmegaE = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiY = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiBX = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiU = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiE = gsl_matrix_alloc(d_size, n_size);

  // Large matrices for NR.
  // Each dxd block is H_k^{-1}.
  gsl_matrix *Hi_all = gsl_matrix_alloc(d_size, d_size * n_size);

  // Each column is H_k^{-1}y_k.
  gsl_matrix *Hiy_all = gsl_matrix_alloc(d_size, n_size);

  // Each dcxdc block is x_k \otimes H_k^{-1}.
  gsl_matrix *xHi_all = gsl_matrix_alloc(dc_size, d_size * n_size);
  gsl_matrix *Hessian = gsl_matrix_alloc(v_size * 2, v_size * 2);

  gsl_vector *x = gsl_vector_alloc(n_size);
  gsl_vector *x_miss = gsl_vector_alloc(n_size);

  gsl_matrix *Y = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *X = gsl_matrix_alloc(c_size + 1, n_size);
  gsl_matrix *V_g = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *V_e = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *B = gsl_matrix_alloc(d_size, c_size + 1);
  gsl_vector *beta = gsl_vector_alloc(d_size);
  gsl_matrix *Vbeta = gsl_matrix_alloc(d_size, d_size);

  // Null estimates for initial values.
  gsl_matrix *V_g_null = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *V_e_null = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *B_null = gsl_matrix_alloc(d_size, c_size + 1);
  gsl_matrix *se_B_null = gsl_matrix_alloc(d_size, c_size);

  gsl_matrix_view X_sub = gsl_matrix_submatrix(X, 0, 0, c_size, n_size);
  gsl_matrix_view B_sub = gsl_matrix_submatrix(B, 0, 0, d_size, c_size);
  gsl_matrix_view xHi_all_sub =
      gsl_matrix_submatrix(xHi_all, 0, 0, d_size * c_size, d_size * n_size);

  gsl_matrix_transpose_memcpy(Y, UtY);

  gsl_matrix_transpose_memcpy(&X_sub.matrix, UtW);

  gsl_vector_view X_row = gsl_matrix_row(X, c_size);
  gsl_vector_set_zero(&X_row.vector);
  gsl_vector_view B_col = gsl_matrix_column(B, c_size);
  gsl_vector_set_zero(&B_col.vector);

  MphInitial(em_iter, em_prec, nr_iter, nr_prec, eval, &X_sub.matrix, Y, l_min,
             l_max, n_region, V_g, V_e, &B_sub.matrix);
  logl_H0 = MphEM('R', em_iter, em_prec, eval, &X_sub.matrix, Y, U_hat, E_hat,
                  OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g,
                  V_e, &B_sub.matrix);
  logl_H0 = MphNR('R', nr_iter, nr_prec, eval, &X_sub.matrix, Y, Hi_all,
                  &xHi_all_sub.matrix, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b,
                  crt_c);
  MphCalcBeta(eval, &X_sub.matrix, Y, V_g, V_e, UltVehiY, &B_sub.matrix,
              se_B_null);

  c = 0;
  Vg_remle_null.clear();
  Ve_remle_null.clear();
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = i; j < d_size; j++) {
      Vg_remle_null.push_back(gsl_matrix_get(V_g, i, j));
      Ve_remle_null.push_back(gsl_matrix_get(V_e, i, j));
      VVg_remle_null.push_back(gsl_matrix_get(Hessian, c, c));
      VVe_remle_null.push_back(gsl_matrix_get(Hessian, c + v_size, c + v_size));
      c++;
    }
  }
  beta_remle_null.clear();
  se_beta_remle_null.clear();
  for (size_t i = 0; i < se_B_null->size1; i++) {
    for (size_t j = 0; j < se_B_null->size2; j++) {
      beta_remle_null.push_back(gsl_matrix_get(B, i, j));
      se_beta_remle_null.push_back(gsl_matrix_get(se_B_null, i, j));
    }
  }
  logl_remle_H0 = logl_H0;

  cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  cout.precision(4);

  cout << "REMLE estimate for Vg in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_g, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Vg): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c, c)) << "\t";
    }
    cout << endl;
  }
  cout << "REMLE estimate for Ve in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_e, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Ve): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c + v_size, c + v_size)) << "\t";
    }
    cout << endl;
  }
  cout << "REMLE likelihood = " << logl_H0 << endl;

  logl_H0 = MphEM('L', em_iter, em_prec, eval, &X_sub.matrix, Y, U_hat, E_hat,
                  OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g,
                  V_e, &B_sub.matrix);
  logl_H0 = MphNR('L', nr_iter, nr_prec, eval, &X_sub.matrix, Y, Hi_all,
                  &xHi_all_sub.matrix, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b,
                  crt_c);
  MphCalcBeta(eval, &X_sub.matrix, Y, V_g, V_e, UltVehiY, &B_sub.matrix,
              se_B_null);

  c = 0;
  Vg_mle_null.clear();
  Ve_mle_null.clear();
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = i; j < d_size; j++) {
      Vg_mle_null.push_back(gsl_matrix_get(V_g, i, j));
      Ve_mle_null.push_back(gsl_matrix_get(V_e, i, j));
      VVg_mle_null.push_back(gsl_matrix_get(Hessian, c, c));
      VVe_mle_null.push_back(gsl_matrix_get(Hessian, c + v_size, c + v_size));
      c++;
    }
  }
  beta_mle_null.clear();
  se_beta_mle_null.clear();
  for (size_t i = 0; i < se_B_null->size1; i++) {
    for (size_t j = 0; j < se_B_null->size2; j++) {
      beta_mle_null.push_back(gsl_matrix_get(B, i, j));
      se_beta_mle_null.push_back(gsl_matrix_get(se_B_null, i, j));
    }
  }
  logl_mle_H0 = logl_H0;

  cout << "MLE estimate for Vg in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_g, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Vg): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c, c)) << "\t";
    }
    cout << endl;
  }
  cout << "MLE estimate for Ve in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_e, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Ve): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c + v_size, c + v_size)) << "\t";
    }
    cout << endl;
  }
  cout << "MLE likelihood = " << logl_H0 << endl;

  vector<double> v_beta, v_Vg, v_Ve, v_Vbeta;
  for (size_t i = 0; i < d_size; i++) {
    v_beta.push_back(0.0);
  }
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = i; j < d_size; j++) {
      v_Vg.push_back(0.0);
      v_Ve.push_back(0.0);
      v_Vbeta.push_back(0.0);
    }
  }

  gsl_matrix_memcpy(V_g_null, V_g);
  gsl_matrix_memcpy(V_e_null, V_e);
  gsl_matrix_memcpy(B_null, B);

  // Start reading genotypes and analyze.
  size_t csnp = 0, t_last = 0;
  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    if (indicator_snp[t] == 0) {
      continue;
    }
    t_last++;
  }
  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    safeGetline(infile, line).eof();
    if (t % d_pace == 0 || t == (ns_total - 1)) {
      ProgressBar("Reading SNPs", t, ns_total - 1);
    }
    if (indicator_snp[t] == 0) {
      continue;
    }

    ch_ptr = strtok_safe((char *)line.c_str(), " , \t");
    ch_ptr = strtok_safe(NULL, " , \t");
    ch_ptr = strtok_safe(NULL, " , \t");

    x_mean = 0.0;
    c_phen = 0;
    n_miss = 0;
    gsl_vector_set_zero(x_miss);
    for (size_t i = 0; i < ni_total; ++i) {
      ch_ptr = strtok_safe(NULL, " , \t");
      if (indicator_idv[i] == 0) {
        continue;
      }

      if (strcmp(ch_ptr, "NA") == 0) {
        gsl_vector_set(x_miss, c_phen, 0.0);
        n_miss++;
      } else {
        geno = atof(ch_ptr);

        gsl_vector_set(x, c_phen, geno);
        gsl_vector_set(x_miss, c_phen, 1.0);
        x_mean += geno;
      }
      c_phen++;
    }

    x_mean /= (double)(ni_test - n_miss);

    for (size_t i = 0; i < ni_test; ++i) {
      if (gsl_vector_get(x_miss, i) == 0) {
        gsl_vector_set(x, i, x_mean);
      }
      geno = gsl_vector_get(x, i);
    }

    gsl_vector_view Xlarge_col = gsl_matrix_column(Xlarge, csnp % msize);
    gsl_vector_memcpy(&Xlarge_col.vector, x);
    csnp++;

    if (csnp % msize == 0 || csnp == t_last) {
      size_t l = 0;
      if (csnp % msize == 0) {
        l = msize;
      } else {
        l = csnp % msize;
      }

      gsl_matrix_view Xlarge_sub =
          gsl_matrix_submatrix(Xlarge, 0, 0, Xlarge->size1, l);
      gsl_matrix_view UtXlarge_sub =
          gsl_matrix_submatrix(UtXlarge, 0, 0, UtXlarge->size1, l);

      time_start = clock();
      fast_dgemm("T", "N", 1.0, U, &Xlarge_sub.matrix, 0.0,
                 &UtXlarge_sub.matrix);
      time_UtX += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

      gsl_matrix_set_zero(Xlarge);

      for (size_t i = 0; i < l; i++) {
        gsl_vector_view UtXlarge_col = gsl_matrix_column(UtXlarge, i);
        gsl_vector_memcpy(&X_row.vector, &UtXlarge_col.vector);

        // Initial values.
        gsl_matrix_memcpy(V_g, V_g_null);
        gsl_matrix_memcpy(V_e, V_e_null);
        gsl_matrix_memcpy(B, B_null);

        time_start = clock();

        // 3 is before 1.
        if (a_mode == 3 || a_mode == 4) {
          p_score = MphCalcP(eval, &X_row.vector, &X_sub.matrix, Y, V_g_null,
                             V_e_null, UltVehiY, beta, Vbeta);
          if (p_score < p_nr && crt == 1) {
            logl_H1 = MphNR('R', 1, nr_prec * 10, eval, X, Y, Hi_all, xHi_all,
                            Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);
            p_score = PCRT(3, d_size, p_score, crt_a, crt_b, crt_c);
          }
        }

        if (a_mode == 2 || a_mode == 4) {
          logl_H1 = MphEM('L', em_iter / 10, em_prec * 10, eval, X, Y, U_hat,
                          E_hat, OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU,
                          UltVehiE, V_g, V_e, B);

          // Calculate beta and Vbeta.
          p_lrt = MphCalcP(eval, &X_row.vector, &X_sub.matrix, Y, V_g, V_e,
                           UltVehiY, beta, Vbeta);
          p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), (double)d_size);

          if (p_lrt < p_nr) {
            logl_H1 =
                MphNR('L', nr_iter / 10, nr_prec * 10, eval, X, Y, Hi_all,
                      xHi_all, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);

            // Calculate beta and Vbeta.
            p_lrt = MphCalcP(eval, &X_row.vector, &X_sub.matrix, Y, V_g, V_e,
                             UltVehiY, beta, Vbeta);
            p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), (double)d_size);

            if (crt == 1) {
              p_lrt = PCRT(2, d_size, p_lrt, crt_a, crt_b, crt_c);
            }
          }
        }

        if (a_mode == 1 || a_mode == 4) {
          logl_H1 = MphEM('R', em_iter / 10, em_prec * 10, eval, X, Y, U_hat,
                          E_hat, OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU,
                          UltVehiE, V_g, V_e, B);
          p_wald = MphCalcP(eval, &X_row.vector, &X_sub.matrix, Y, V_g, V_e,
                            UltVehiY, beta, Vbeta);

          if (p_wald < p_nr) {
            logl_H1 =
                MphNR('R', nr_iter / 10, nr_prec * 10, eval, X, Y, Hi_all,
                      xHi_all, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);
            p_wald = MphCalcP(eval, &X_row.vector, &X_sub.matrix, Y, V_g, V_e,
                              UltVehiY, beta, Vbeta);

            if (crt == 1) {
              p_wald = PCRT(1, d_size, p_wald, crt_a, crt_b, crt_c);
            }
          }
        }

        time_opt += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

        // Store summary data.
        for (size_t i = 0; i < d_size; i++) {
          v_beta[i] = gsl_vector_get(beta, i);
        }

        c = 0;
        for (size_t i = 0; i < d_size; i++) {
          for (size_t j = i; j < d_size; j++) {
            v_Vg[c] = gsl_matrix_get(V_g, i, j);
            v_Ve[c] = gsl_matrix_get(V_e, i, j);
            v_Vbeta[c] = gsl_matrix_get(Vbeta, i, j);
            c++;
          }
        }

        MPHSUMSTAT SNPs = {v_beta, p_wald, p_lrt, p_score, v_Vg, v_Ve, v_Vbeta};
        sumStat.push_back(SNPs);
      }
    }
  }
  cout << endl;

  infile.close();
  infile.clear();

  gsl_matrix_free(U_hat);
  gsl_matrix_free(E_hat);
  gsl_matrix_free(OmegaU);
  gsl_matrix_free(OmegaE);
  gsl_matrix_free(UltVehiY);
  gsl_matrix_free(UltVehiBX);
  gsl_matrix_free(UltVehiU);
  gsl_matrix_free(UltVehiE);

  gsl_matrix_free(Hi_all);
  gsl_matrix_free(Hiy_all);
  gsl_matrix_free(xHi_all);
  gsl_matrix_free(Hessian);

  gsl_vector_free(x);
  gsl_vector_free(x_miss);

  gsl_matrix_free(Y);
  gsl_matrix_free(X);
  gsl_matrix_free(V_g);
  gsl_matrix_free(V_e);
  gsl_matrix_free(B);
  gsl_vector_free(beta);
  gsl_matrix_free(Vbeta);

  gsl_matrix_free(V_g_null);
  gsl_matrix_free(V_e_null);
  gsl_matrix_free(B_null);
  gsl_matrix_free(se_B_null);

  gsl_matrix_free(Xlarge);
  gsl_matrix_free(UtXlarge);

  return;
}

void MVLMM::AnalyzePlink(const gsl_matrix *U, const gsl_vector *eval,
                         const gsl_matrix *UtW, const gsl_matrix *UtY) {
  debug_msg("entering");
  string file_bed = file_bfile + ".bed";
  ifstream infile(file_bed.c_str(), ios::binary);
  if (!infile) {
    cout << "error reading bed file:" << file_bed << endl;
    return;
  }

  clock_t time_start = clock();
  time_UtX = 0;
  time_opt = 0;

  char ch[1];
  bitset<8> b;

  double logl_H0 = 0.0, logl_H1 = 0.0, p_wald = 0, p_lrt = 0, p_score = 0;
  double crt_a, crt_b, crt_c;
  int n_bit, n_miss, ci_total, ci_test;
  double geno, x_mean;
  size_t c = 0;
  size_t n_size = UtY->size1, d_size = UtY->size2, c_size = UtW->size2;
  size_t dc_size = d_size * (c_size + 1), v_size = d_size * (d_size + 1) / 2;

  // Create a large matrix.
  size_t msize = LMM_BATCH_SIZE;
  gsl_matrix *Xlarge = gsl_matrix_alloc(U->size1, msize);
  gsl_matrix *UtXlarge = gsl_matrix_alloc(U->size1, msize);
  gsl_matrix_set_zero(Xlarge);

  // Large matrices for EM.
  gsl_matrix *U_hat = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *E_hat = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *OmegaU = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *OmegaE = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiY = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiBX = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiU = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiE = gsl_matrix_alloc(d_size, n_size);

  // Large matrices for NR.
  // Each dxd block is H_k^{-1}.
  gsl_matrix *Hi_all = gsl_matrix_alloc(d_size, d_size * n_size);

  // Each column is H_k^{-1}y_k.
  gsl_matrix *Hiy_all = gsl_matrix_alloc(d_size, n_size);

  // Each dcxdc block is x_k\otimes H_k^{-1}.
  gsl_matrix *xHi_all = gsl_matrix_alloc(dc_size, d_size * n_size);

  gsl_matrix *Hessian = gsl_matrix_alloc(v_size * 2, v_size * 2);

  gsl_vector *x = gsl_vector_alloc(n_size);

  gsl_matrix *Y = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *X = gsl_matrix_alloc(c_size + 1, n_size);
  gsl_matrix *V_g = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *V_e = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *B = gsl_matrix_alloc(d_size, c_size + 1);
  gsl_vector *beta = gsl_vector_alloc(d_size);
  gsl_matrix *Vbeta = gsl_matrix_alloc(d_size, d_size);

  // Null estimates for initial values.
  gsl_matrix *V_g_null = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *V_e_null = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *B_null = gsl_matrix_alloc(d_size, c_size + 1);
  gsl_matrix *se_B_null = gsl_matrix_alloc(d_size, c_size);

  gsl_matrix_view X_sub = gsl_matrix_submatrix(X, 0, 0, c_size, n_size);
  gsl_matrix_view B_sub = gsl_matrix_submatrix(B, 0, 0, d_size, c_size);
  gsl_matrix_view xHi_all_sub =
      gsl_matrix_submatrix(xHi_all, 0, 0, d_size * c_size, d_size * n_size);

  gsl_matrix_transpose_memcpy(Y, UtY);
  gsl_matrix_transpose_memcpy(&X_sub.matrix, UtW);

  gsl_vector_view X_row = gsl_matrix_row(X, c_size);
  gsl_vector_set_zero(&X_row.vector);
  gsl_vector_view B_col = gsl_matrix_column(B, c_size);
  gsl_vector_set_zero(&B_col.vector);

  MphInitial(em_iter, em_prec, nr_iter, nr_prec, eval, &X_sub.matrix, Y, l_min,
             l_max, n_region, V_g, V_e, &B_sub.matrix);

  write(eval,"eval4");

  logl_H0 = MphEM('R', em_iter, em_prec, eval, &X_sub.matrix, Y, U_hat, E_hat,
                  OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g,
                  V_e, &B_sub.matrix);
  logl_H0 = MphNR('R', nr_iter, nr_prec, eval, &X_sub.matrix, Y, Hi_all,
                  &xHi_all_sub.matrix, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b,
                  crt_c);
  MphCalcBeta(eval, &X_sub.matrix, Y, V_g, V_e, UltVehiY, &B_sub.matrix,
              se_B_null);

  c = 0;
  Vg_remle_null.clear();
  Ve_remle_null.clear();
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = i; j < d_size; j++) {
      Vg_remle_null.push_back(gsl_matrix_get(V_g, i, j));
      Ve_remle_null.push_back(gsl_matrix_get(V_e, i, j));
      VVg_remle_null.push_back(gsl_matrix_get(Hessian, c, c));
      VVe_remle_null.push_back(gsl_matrix_get(Hessian, c + v_size, c + v_size));
      c++;
    }
  }
  beta_remle_null.clear();
  se_beta_remle_null.clear();
  for (size_t i = 0; i < se_B_null->size1; i++) {
    for (size_t j = 0; j < se_B_null->size2; j++) {
      beta_remle_null.push_back(gsl_matrix_get(B, i, j));
      se_beta_remle_null.push_back(gsl_matrix_get(se_B_null, i, j));
    }
  }
  logl_remle_H0 = logl_H0;

  cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  cout.precision(4);
  cout << "REMLE estimate for Vg in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_g, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Vg): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c, c)) << "\t";
    }
    cout << endl;
  }
  cout << "REMLE estimate for Ve in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_e, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Ve): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      auto v = gsl_matrix_get(Hessian, c + v_size, c + v_size);
      if (is_strict_mode())
        enforce_msg(v >= 0,"se(Ve) is not valid");

      cout << safe_sqrt(v) << "\t";
    }
    cout << endl;
  }
  cout << "REMLE likelihood = " << logl_H0 << endl;

  logl_H0 = MphEM('L', em_iter, em_prec, eval, &X_sub.matrix, Y, U_hat, E_hat,
                  OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g,
                  V_e, &B_sub.matrix);
  logl_H0 = MphNR('L', nr_iter, nr_prec, eval, &X_sub.matrix, Y, Hi_all,
                  &xHi_all_sub.matrix, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b,
                  crt_c);
  MphCalcBeta(eval, &X_sub.matrix, Y, V_g, V_e, UltVehiY, &B_sub.matrix,
              se_B_null);

  c = 0;
  Vg_mle_null.clear();
  Ve_mle_null.clear();
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = i; j < d_size; j++) {
      Vg_mle_null.push_back(gsl_matrix_get(V_g, i, j));
      Ve_mle_null.push_back(gsl_matrix_get(V_e, i, j));
      VVg_mle_null.push_back(gsl_matrix_get(Hessian, c, c));
      VVe_mle_null.push_back(gsl_matrix_get(Hessian, c + v_size, c + v_size));
      c++;
    }
  }
  beta_mle_null.clear();
  se_beta_mle_null.clear();
  for (size_t i = 0; i < se_B_null->size1; i++) {
    for (size_t j = 0; j < se_B_null->size2; j++) {
      beta_mle_null.push_back(gsl_matrix_get(B, i, j));
      se_beta_mle_null.push_back(gsl_matrix_get(se_B_null, i, j));
    }
  }
  logl_mle_H0 = logl_H0;

  cout << "MLE estimate for Vg in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_g, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Vg): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c, c)) << "\t";
    }
    cout << endl;
  }
  cout << "MLE estimate for Ve in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_e, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Ve): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      auto v = gsl_matrix_get(Hessian, c + v_size, c + v_size);
      if (is_strict_mode())
        enforce_msg(v >= 0,"se(Ve) is not valid");
      cout << safe_sqrt(v) << "\t";
    }
    cout << endl;
  }
  cout << "MLE likelihood = " << logl_H0 << endl;

  vector<double> v_beta, v_Vg, v_Ve, v_Vbeta;
  for (size_t i = 0; i < d_size; i++) {
    v_beta.push_back(0.0);
  }
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = i; j < d_size; j++) {
      v_Vg.push_back(0.0);
      v_Ve.push_back(0.0);
      v_Vbeta.push_back(0.0);
    }
  }

  gsl_matrix_memcpy(V_g_null, V_g);
  gsl_matrix_memcpy(V_e_null, V_e);
  gsl_matrix_memcpy(B_null, B);

  // Start reading genotypes and analyze.
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

  size_t csnp = 0, t_last = 0;
  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    if (indicator_snp[t] == 0) {
      continue;
    }
    t_last++;
  }
  for (vector<SNPINFO>::size_type t = 0; t < snpInfo.size(); ++t) {
    if (t % d_pace == 0 || t == snpInfo.size() - 1) {
      ProgressBar("Reading SNPs", t, snpInfo.size() - 1);
    }
    if (indicator_snp[t] == 0) {
      continue;
    }

    // n_bit, and 3 is the number of magic numbers.
    infile.seekg(t * n_bit + 3);

    // read genotypes
    x_mean = 0.0;
    n_miss = 0;
    ci_total = 0;
    ci_test = 0;
    for (int i = 0; i < n_bit; ++i) {
      infile.read(ch, 1);
      b = ch[0];

      // Minor allele homozygous: 2.0; major: 0.0;
      for (size_t j = 0; j < 4; ++j) {
        if ((i == (n_bit - 1)) && ci_total == (int)ni_total) {
          break;
        }
        if (indicator_idv[ci_total] == 0) {
          ci_total++;
          continue;
        }

        if (b[2 * j] == 0) {
          if (b[2 * j + 1] == 0) {
            gsl_vector_set(x, ci_test, 2);
            x_mean += 2.0;
          } else {
            gsl_vector_set(x, ci_test, 1);
            x_mean += 1.0;
          }
        } else {
          if (b[2 * j + 1] == 1) {
            gsl_vector_set(x, ci_test, 0);
          } else {
            gsl_vector_set(x, ci_test, -9);
            n_miss++;
          }
        }

        ci_total++;
        ci_test++;
      }
    }

    x_mean /= (double)(ni_test - n_miss);

    for (size_t i = 0; i < ni_test; ++i) {
      geno = gsl_vector_get(x, i);
      if (geno == -9) {
        gsl_vector_set(x, i, x_mean);
        geno = x_mean;
      }
    }

    gsl_vector_view Xlarge_col = gsl_matrix_column(Xlarge, csnp % msize);
    gsl_vector_memcpy(&Xlarge_col.vector, x);
    csnp++;

    if (csnp % msize == 0 || csnp == t_last) {
      size_t l = 0;
      if (csnp % msize == 0) {
        l = msize;
      } else {
        l = csnp % msize;
      }

      gsl_matrix_view Xlarge_sub =
          gsl_matrix_submatrix(Xlarge, 0, 0, Xlarge->size1, l);
      gsl_matrix_view UtXlarge_sub =
          gsl_matrix_submatrix(UtXlarge, 0, 0, UtXlarge->size1, l);

      time_start = clock();
      fast_dgemm("T", "N", 1.0, U, &Xlarge_sub.matrix, 0.0,
                     &UtXlarge_sub.matrix);
      time_UtX += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

      gsl_matrix_set_zero(Xlarge);

      for (size_t i = 0; i < l; i++) {
        gsl_vector_view UtXlarge_col = gsl_matrix_column(UtXlarge, i);
        gsl_vector_memcpy(&X_row.vector, &UtXlarge_col.vector);

        // Initial values.
        gsl_matrix_memcpy(V_g, V_g_null);
        gsl_matrix_memcpy(V_e, V_e_null);
        gsl_matrix_memcpy(B, B_null);

        time_start = clock();

        // 3 is before 1.
        if (a_mode == 3 || a_mode == 4) {
          p_score = MphCalcP(eval, &X_row.vector, &X_sub.matrix, Y, V_g_null,
                             V_e_null, UltVehiY, beta, Vbeta);

          if (p_score < p_nr && crt == 1) {
            logl_H1 = MphNR('R', 1, nr_prec * 10, eval, X, Y, Hi_all, xHi_all,
                            Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);
            p_score = PCRT(3, d_size, p_score, crt_a, crt_b, crt_c);
          }
        }

        if (a_mode == 2 || a_mode == 4) {
          logl_H1 = MphEM('L', em_iter / 10, em_prec * 10, eval, X, Y, U_hat,
                          E_hat, OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU,
                          UltVehiE, V_g, V_e, B);

          // Calculate beta and Vbeta.
          p_lrt = MphCalcP(eval, &X_row.vector, &X_sub.matrix, Y, V_g, V_e,
                           UltVehiY, beta, Vbeta);
          p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), (double)d_size);

          if (p_lrt < p_nr) {
            logl_H1 =
                MphNR('L', nr_iter / 10, nr_prec * 10, eval, X, Y, Hi_all,
                      xHi_all, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);

            // Calculate beta and Vbeta.
            p_lrt = MphCalcP(eval, &X_row.vector, &X_sub.matrix, Y, V_g, V_e,
                             UltVehiY, beta, Vbeta);
            p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), (double)d_size);
            if (crt == 1) {
              p_lrt = PCRT(2, d_size, p_lrt, crt_a, crt_b, crt_c);
            }
          }
        }

        if (a_mode == 1 || a_mode == 4) {
          logl_H1 = MphEM('R', em_iter / 10, em_prec * 10, eval, X, Y, U_hat,
                          E_hat, OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU,
                          UltVehiE, V_g, V_e, B);
          p_wald = MphCalcP(eval, &X_row.vector, &X_sub.matrix, Y, V_g, V_e,
                            UltVehiY, beta, Vbeta);

          if (p_wald < p_nr) {
            logl_H1 =
                MphNR('R', nr_iter / 10, nr_prec * 10, eval, X, Y, Hi_all,
                      xHi_all, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);
            p_wald = MphCalcP(eval, &X_row.vector, &X_sub.matrix, Y, V_g, V_e,
                              UltVehiY, beta, Vbeta);

            if (crt == 1) {
              p_wald = PCRT(1, d_size, p_wald, crt_a, crt_b, crt_c);
            }
          }
        }

        time_opt += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

        // Store summary data.
        for (size_t i = 0; i < d_size; i++) {
          v_beta[i] = gsl_vector_get(beta, i);
        }

        c = 0;
        for (size_t i = 0; i < d_size; i++) {
          for (size_t j = i; j < d_size; j++) {
            v_Vg[c] = gsl_matrix_get(V_g, i, j);
            v_Ve[c] = gsl_matrix_get(V_e, i, j);
            v_Vbeta[c] = gsl_matrix_get(Vbeta, i, j);
            c++;
          }
        }

        MPHSUMSTAT SNPs = {v_beta, p_wald, p_lrt, p_score, v_Vg, v_Ve, v_Vbeta};
        sumStat.push_back(SNPs);
      }
    }
  }
  cout << endl;

  infile.close();
  infile.clear();

  gsl_matrix_free(U_hat);
  gsl_matrix_free(E_hat);
  gsl_matrix_free(OmegaU);
  gsl_matrix_free(OmegaE);
  gsl_matrix_free(UltVehiY);
  gsl_matrix_free(UltVehiBX);
  gsl_matrix_free(UltVehiU);
  gsl_matrix_free(UltVehiE);

  gsl_matrix_free(Hi_all);
  gsl_matrix_free(Hiy_all);
  gsl_matrix_free(xHi_all);
  gsl_matrix_free(Hessian);

  gsl_vector_free(x);

  gsl_matrix_free(Y);
  gsl_matrix_free(X);
  gsl_matrix_free(V_g);
  gsl_matrix_free(V_e);
  gsl_matrix_free(B);
  gsl_vector_free(beta);
  gsl_matrix_free(Vbeta);

  gsl_matrix_free(V_g_null);
  gsl_matrix_free(V_e_null);
  gsl_matrix_free(B_null);
  gsl_matrix_free(se_B_null);

  gsl_matrix_free(Xlarge);
  gsl_matrix_free(UtXlarge);

  return;
}

// Calculate Vg, Ve, B, se(B) in the null mvLMM model.
// Both B and se_B are d by c matrices.
void CalcMvLmmVgVeBeta(const gsl_vector *eval, const gsl_matrix *UtW,
                       const gsl_matrix *UtY, const size_t em_iter,
                       const size_t nr_iter, const double em_prec,
                       const double nr_prec, const double l_min,
                       const double l_max, const size_t n_region,
                       gsl_matrix *V_g, gsl_matrix *V_e, gsl_matrix *B,
                       gsl_matrix *se_B) {
  size_t n_size = UtY->size1, d_size = UtY->size2, c_size = UtW->size2;
  size_t dc_size = d_size * c_size, v_size = d_size * (d_size + 1) / 2;

  double crt_a, crt_b, crt_c;

  // Large matrices for EM.
  gsl_matrix *U_hat = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *E_hat = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *OmegaU = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *OmegaE = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiY = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiBX = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiU = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiE = gsl_matrix_alloc(d_size, n_size);

  // Large matrices for NR.
  // Each dxd block is H_k^{-1}.
  gsl_matrix *Hi_all = gsl_matrix_alloc(d_size, d_size * n_size);

  // Each column is H_k^{-1}y_k.
  gsl_matrix *Hiy_all = gsl_matrix_alloc(d_size, n_size);

  // Each dcxdc block is x_k\otimes H_k^{-1}.
  gsl_matrix *xHi_all = gsl_matrix_alloc(dc_size, d_size * n_size);
  gsl_matrix *Hessian = gsl_matrix_alloc(v_size * 2, v_size * 2);

  // Transpose matrices.
  gsl_matrix *Y = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *W = gsl_matrix_alloc(c_size, n_size);
  gsl_matrix_transpose_memcpy(Y, UtY);
  gsl_matrix_transpose_memcpy(W, UtW);

  // Initial, EM, NR, and calculate B.
  MphInitial(em_iter, em_prec, nr_iter, nr_prec, eval, W, Y, l_min, l_max,
             n_region, V_g, V_e, B);
  MphEM('R', em_iter, em_prec, eval, W, Y, U_hat, E_hat, OmegaU, OmegaE,
        UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g, V_e, B);
  MphNR('R', nr_iter, nr_prec, eval, W, Y, Hi_all, xHi_all, Hiy_all, V_g,
        V_e, Hessian, crt_a, crt_b, crt_c);
  MphCalcBeta(eval, W, Y, V_g, V_e, UltVehiY, B, se_B);

  // Free matrices.
  gsl_matrix_free(U_hat);
  gsl_matrix_free(E_hat);
  gsl_matrix_free(OmegaU);
  gsl_matrix_free(OmegaE);
  gsl_matrix_free(UltVehiY);
  gsl_matrix_free(UltVehiBX);
  gsl_matrix_free(UltVehiU);
  gsl_matrix_free(UltVehiE);

  gsl_matrix_free(Hi_all);
  gsl_matrix_free(Hiy_all);
  gsl_matrix_free(xHi_all);
  gsl_matrix_free(Hessian);

  gsl_matrix_free(Y);
  gsl_matrix_free(W);

  return;
}

void MVLMM::AnalyzeBimbamGXE(const gsl_matrix *U, const gsl_vector *eval,
                             const gsl_matrix *UtW, const gsl_matrix *UtY,
                             const gsl_vector *env) {
  debug_msg("entering");
  igzstream infile(file_geno.c_str(), igzstream::in);
  if (!infile) {
    cout << "error reading genotype file:" << file_geno << endl;
    return;
  }

  clock_t time_start = clock();
  time_UtX = 0;
  time_opt = 0;

  string line;
  char *ch_ptr;

  double logl_H0 = 0.0, logl_H1 = 0.0, p_wald = 0, p_lrt = 0, p_score = 0;
  double crt_a, crt_b, crt_c;
  int n_miss, c_phen;
  double geno, x_mean;
  size_t c = 0;
  size_t n_size = UtY->size1, d_size = UtY->size2, c_size = UtW->size2 + 2;
  size_t dc_size = d_size * (c_size + 1), v_size = d_size * (d_size + 1) / 2;

  // Large matrices for EM.
  gsl_matrix *U_hat = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *E_hat = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *OmegaU = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *OmegaE = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiY = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiBX = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiU = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiE = gsl_matrix_alloc(d_size, n_size);

  // Large matrices for NR.
  // Each dxd block is H_k^{-1}.
  gsl_matrix *Hi_all = gsl_matrix_alloc(d_size, d_size * n_size);

  // Each column is H_k^{-1}y_k.
  gsl_matrix *Hiy_all = gsl_matrix_alloc(d_size, n_size);

  // Each dcxdc block is x_k\otimes H_k^{-1}.
  gsl_matrix *xHi_all = gsl_matrix_alloc(dc_size, d_size * n_size);
  gsl_matrix *Hessian = gsl_matrix_alloc(v_size * 2, v_size * 2);

  gsl_vector *x = gsl_vector_alloc(n_size);
  gsl_vector *x_miss = gsl_vector_alloc(n_size);

  gsl_matrix *Y = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *X = gsl_matrix_alloc(c_size + 1, n_size);
  gsl_matrix *V_g = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *V_e = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *B = gsl_matrix_alloc(d_size, c_size + 1);
  gsl_vector *beta = gsl_vector_alloc(d_size);
  gsl_matrix *Vbeta = gsl_matrix_alloc(d_size, d_size);

  // Null estimates for initial values; including env but not
  // including x.
  gsl_matrix *V_g_null = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *V_e_null = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *B_null = gsl_matrix_alloc(d_size, c_size + 1);
  gsl_matrix *se_B_null1 = gsl_matrix_alloc(d_size, c_size - 1);
  gsl_matrix *se_B_null2 = gsl_matrix_alloc(d_size, c_size);

  gsl_matrix_view X_sub1 = gsl_matrix_submatrix(X, 0, 0, c_size - 1, n_size);
  gsl_matrix_view B_sub1 = gsl_matrix_submatrix(B, 0, 0, d_size, c_size - 1);
  gsl_matrix_view xHi_all_sub1 = gsl_matrix_submatrix(
      xHi_all, 0, 0, d_size * (c_size - 1), d_size * n_size);

  gsl_matrix_view X_sub2 = gsl_matrix_submatrix(X, 0, 0, c_size, n_size);
  gsl_matrix_view B_sub2 = gsl_matrix_submatrix(B, 0, 0, d_size, c_size);
  gsl_matrix_view xHi_all_sub2 =
      gsl_matrix_submatrix(xHi_all, 0, 0, d_size * c_size, d_size * n_size);

  gsl_matrix_transpose_memcpy(Y, UtY);

  gsl_matrix_view X_sub0 = gsl_matrix_submatrix(X, 0, 0, c_size - 2, n_size);
  gsl_matrix_transpose_memcpy(&X_sub0.matrix, UtW);
  gsl_vector_view X_row0 = gsl_matrix_row(X, c_size - 2);
  gsl_blas_dgemv(CblasTrans, 1.0, U, env, 0.0, &X_row0.vector);

  gsl_vector_view X_row1 = gsl_matrix_row(X, c_size - 1);
  gsl_vector_set_zero(&X_row1.vector);
  gsl_vector_view X_row2 = gsl_matrix_row(X, c_size);
  gsl_vector_set_zero(&X_row2.vector);

  gsl_vector_view B_col1 = gsl_matrix_column(B, c_size - 1);
  gsl_vector_set_zero(&B_col1.vector);
  gsl_vector_view B_col2 = gsl_matrix_column(B, c_size);
  gsl_vector_set_zero(&B_col2.vector);

  MphInitial(em_iter, em_prec, nr_iter, nr_prec, eval, &X_sub1.matrix, Y, l_min,
             l_max, n_region, V_g, V_e, &B_sub1.matrix);
  logl_H0 = MphEM('R', em_iter, em_prec, eval, &X_sub1.matrix, Y, U_hat, E_hat,
                  OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g,
                  V_e, &B_sub1.matrix);
  logl_H0 = MphNR('R', nr_iter, nr_prec, eval, &X_sub1.matrix, Y, Hi_all,
                  &xHi_all_sub1.matrix, Hiy_all, V_g, V_e, Hessian, crt_a,
                  crt_b, crt_c);
  MphCalcBeta(eval, &X_sub1.matrix, Y, V_g, V_e, UltVehiY, &B_sub1.matrix,
              se_B_null1);

  c = 0;
  Vg_remle_null.clear();
  Ve_remle_null.clear();
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = i; j < d_size; j++) {
      Vg_remle_null.push_back(gsl_matrix_get(V_g, i, j));
      Ve_remle_null.push_back(gsl_matrix_get(V_e, i, j));
      VVg_remle_null.push_back(gsl_matrix_get(Hessian, c, c));
      VVe_remle_null.push_back(gsl_matrix_get(Hessian, c + v_size, c + v_size));
      c++;
    }
  }
  beta_remle_null.clear();
  se_beta_remle_null.clear();
  for (size_t i = 0; i < se_B_null1->size1; i++) {
    for (size_t j = 0; j < se_B_null1->size2; j++) {
      beta_remle_null.push_back(gsl_matrix_get(B, i, j));
      se_beta_remle_null.push_back(gsl_matrix_get(se_B_null1, i, j));
    }
  }
  logl_remle_H0 = logl_H0;

  cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  cout.precision(4);

  cout << "REMLE estimate for Vg in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_g, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Vg): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c, c)) << "\t";
    }
    cout << endl;
  }
  cout << "REMLE estimate for Ve in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_e, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Ve): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c + v_size, c + v_size)) << "\t";
    }
    cout << endl;
  }
  cout << "REMLE likelihood = " << logl_H0 << endl;

  logl_H0 = MphEM('L', em_iter, em_prec, eval, &X_sub1.matrix, Y, U_hat, E_hat,
                  OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g,
                  V_e, &B_sub1.matrix);
  logl_H0 = MphNR('L', nr_iter, nr_prec, eval, &X_sub1.matrix, Y, Hi_all,
                  &xHi_all_sub1.matrix, Hiy_all, V_g, V_e, Hessian, crt_a,
                  crt_b, crt_c);
  MphCalcBeta(eval, &X_sub1.matrix, Y, V_g, V_e, UltVehiY, &B_sub1.matrix,
              se_B_null1);

  c = 0;
  Vg_mle_null.clear();
  Ve_mle_null.clear();
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = i; j < d_size; j++) {
      Vg_mle_null.push_back(gsl_matrix_get(V_g, i, j));
      Ve_mle_null.push_back(gsl_matrix_get(V_e, i, j));
      VVg_mle_null.push_back(gsl_matrix_get(Hessian, c, c));
      VVe_mle_null.push_back(gsl_matrix_get(Hessian, c + v_size, c + v_size));
      c++;
    }
  }
  beta_mle_null.clear();
  se_beta_mle_null.clear();
  for (size_t i = 0; i < se_B_null1->size1; i++) {
    for (size_t j = 0; j < se_B_null1->size2; j++) {
      beta_mle_null.push_back(gsl_matrix_get(B, i, j));
      se_beta_mle_null.push_back(gsl_matrix_get(se_B_null1, i, j));
    }
  }
  logl_mle_H0 = logl_H0;

  cout << "MLE estimate for Vg in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_g, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Vg): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c, c)) << "\t";
    }
    cout << endl;
  }
  cout << "MLE estimate for Ve in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_e, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Ve): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c + v_size, c + v_size)) << "\t";
    }
    cout << endl;
  }
  cout << "MLE likelihood = " << logl_H0 << endl;

  vector<double> v_beta, v_Vg, v_Ve, v_Vbeta;
  for (size_t i = 0; i < d_size; i++) {
    v_beta.push_back(0.0);
  }
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = i; j < d_size; j++) {
      v_Vg.push_back(0.0);
      v_Ve.push_back(0.0);
      v_Vbeta.push_back(0.0);
    }
  }

  gsl_matrix_memcpy(V_g_null, V_g);
  gsl_matrix_memcpy(V_e_null, V_e);
  gsl_matrix_memcpy(B_null, B);

  // Start reading genotypes and analyze.
  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    safeGetline(infile, line).eof();
    if (t % d_pace == 0 || t == (ns_total - 1)) {
      ProgressBar("Reading SNPs", t, ns_total - 1);
    }
    if (indicator_snp[t] == 0) {
      continue;
    }

    ch_ptr = strtok_safe((char *)line.c_str(), " , \t");
    ch_ptr = strtok_safe(NULL, " , \t");
    ch_ptr = strtok_safe(NULL, " , \t");

    x_mean = 0.0;
    c_phen = 0;
    n_miss = 0;
    gsl_vector_set_zero(x_miss);
    for (size_t i = 0; i < ni_total; ++i) {
      ch_ptr = strtok_safe(NULL, " , \t");
      if (indicator_idv[i] == 0) {
        continue;
      }

      if (strcmp(ch_ptr, "NA") == 0) {
        gsl_vector_set(x_miss, c_phen, 0.0);
        n_miss++;
      } else {
        geno = atof(ch_ptr);

        gsl_vector_set(x, c_phen, geno);
        gsl_vector_set(x_miss, c_phen, 1.0);
        x_mean += geno;
      }
      c_phen++;
    }

    x_mean /= (double)(ni_test - n_miss);

    for (size_t i = 0; i < ni_test; ++i) {
      if (gsl_vector_get(x_miss, i) == 0) {
        gsl_vector_set(x, i, x_mean);
      }
      geno = gsl_vector_get(x, i);
      if (x_mean > 1) {
        gsl_vector_set(x, i, 2 - geno);
      }
    }

    // Calculate statistics.
    time_start = clock();
    gsl_blas_dgemv(CblasTrans, 1.0, U, x, 0.0, &X_row1.vector);
    gsl_vector_mul(x, env);
    gsl_blas_dgemv(CblasTrans, 1.0, U, x, 0.0, &X_row2.vector);
    time_UtX += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    // initial values
    gsl_matrix_memcpy(V_g, V_g_null);
    gsl_matrix_memcpy(V_e, V_e_null);
    gsl_matrix_memcpy(B, B_null);

    if (a_mode == 2 || a_mode == 3 || a_mode == 4) {
      if (a_mode == 3 || a_mode == 4) {
        logl_H0 = MphEM('R', em_iter / 10, em_prec * 10, eval, &X_sub2.matrix,
                        Y, U_hat, E_hat, OmegaU, OmegaE, UltVehiY, UltVehiBX,
                        UltVehiU, UltVehiE, V_g, V_e, &B_sub2.matrix);
        logl_H0 = MphNR('R', nr_iter / 10, nr_prec * 10, eval, &X_sub2.matrix,
                        Y, Hi_all, &xHi_all_sub2.matrix, Hiy_all, V_g, V_e,
                        Hessian, crt_a, crt_b, crt_c);
        MphCalcBeta(eval, &X_sub2.matrix, Y, V_g, V_e, UltVehiY, &B_sub2.matrix,
                    se_B_null2);
      }

      if (a_mode == 2 || a_mode == 4) {
        logl_H0 = MphEM('L', em_iter / 10, em_prec * 10, eval, &X_sub2.matrix,
                        Y, U_hat, E_hat, OmegaU, OmegaE, UltVehiY, UltVehiBX,
                        UltVehiU, UltVehiE, V_g, V_e, &B_sub2.matrix);
        logl_H0 = MphNR('L', nr_iter / 10, nr_prec * 10, eval, &X_sub2.matrix,
                        Y, Hi_all, &xHi_all_sub2.matrix, Hiy_all, V_g, V_e,
                        Hessian, crt_a, crt_b, crt_c);
        MphCalcBeta(eval, &X_sub2.matrix, Y, V_g, V_e, UltVehiY, &B_sub2.matrix,
                    se_B_null2);
      }
    }

    time_start = clock();

    // 3 is before 1.
    if (a_mode == 3 || a_mode == 4) {
      p_score = MphCalcP(eval, &X_row2.vector, &X_sub2.matrix, Y, V_g_null,
                         V_e_null, UltVehiY, beta, Vbeta);
      if (p_score < p_nr && crt == 1) {
        logl_H1 = MphNR('R', 1, nr_prec * 10, eval, X, Y, Hi_all, xHi_all,
                        Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);
        p_score = PCRT(3, d_size, p_score, crt_a, crt_b, crt_c);
      }
    }

    if (a_mode == 2 || a_mode == 4) {
      logl_H1 = MphEM('L', em_iter / 10, em_prec * 10, eval, X, Y, U_hat, E_hat,
                      OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE,
                      V_g, V_e, B);

      // Calculate beta and Vbeta.
      p_lrt = MphCalcP(eval, &X_row2.vector, &X_sub2.matrix, Y, V_g, V_e,
                       UltVehiY, beta, Vbeta);
      p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), (double)d_size);

      if (p_lrt < p_nr) {
        logl_H1 =
            MphNR('L', nr_iter / 10, nr_prec * 10, eval, X, Y, Hi_all, xHi_all,
                  Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);

        // Calculate beta and Vbeta.
        p_lrt = MphCalcP(eval, &X_row2.vector, &X_sub2.matrix, Y, V_g, V_e,
                         UltVehiY, beta, Vbeta);
        p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), (double)d_size);

        if (crt == 1) {
          p_lrt = PCRT(2, d_size, p_lrt, crt_a, crt_b, crt_c);
        }
      }
    }

    if (a_mode == 1 || a_mode == 4) {
      logl_H1 = MphEM('R', em_iter / 10, em_prec * 10, eval, X, Y, U_hat, E_hat,
                      OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE,
                      V_g, V_e, B);
      p_wald = MphCalcP(eval, &X_row2.vector, &X_sub2.matrix, Y, V_g, V_e,
                        UltVehiY, beta, Vbeta);

      if (p_wald < p_nr) {
        logl_H1 =
            MphNR('R', nr_iter / 10, nr_prec * 10, eval, X, Y, Hi_all, xHi_all,
                  Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);
        p_wald = MphCalcP(eval, &X_row2.vector, &X_sub2.matrix, Y, V_g, V_e,
                          UltVehiY, beta, Vbeta);

        if (crt == 1) {
          p_wald = PCRT(1, d_size, p_wald, crt_a, crt_b, crt_c);
        }
      }
    }

    if (x_mean > 1) {
      gsl_vector_scale(beta, -1.0);
    }

    time_opt += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    // Store summary data.
    for (size_t i = 0; i < d_size; i++) {
      v_beta[i] = gsl_vector_get(beta, i);
    }

    c = 0;
    for (size_t i = 0; i < d_size; i++) {
      for (size_t j = i; j < d_size; j++) {
        v_Vg[c] = gsl_matrix_get(V_g, i, j);
        v_Ve[c] = gsl_matrix_get(V_e, i, j);
        v_Vbeta[c] = gsl_matrix_get(Vbeta, i, j);
        c++;
      }
    }

    MPHSUMSTAT SNPs = {v_beta, p_wald, p_lrt, p_score, v_Vg, v_Ve, v_Vbeta};
    sumStat.push_back(SNPs);
  }
  cout << endl;

  infile.close();
  infile.clear();

  gsl_matrix_free(U_hat);
  gsl_matrix_free(E_hat);
  gsl_matrix_free(OmegaU);
  gsl_matrix_free(OmegaE);
  gsl_matrix_free(UltVehiY);
  gsl_matrix_free(UltVehiBX);
  gsl_matrix_free(UltVehiU);
  gsl_matrix_free(UltVehiE);

  gsl_matrix_free(Hi_all);
  gsl_matrix_free(Hiy_all);
  gsl_matrix_free(xHi_all);
  gsl_matrix_free(Hessian);

  gsl_vector_free(x);
  gsl_vector_free(x_miss);

  gsl_matrix_free(Y);
  gsl_matrix_free(X);
  gsl_matrix_free(V_g);
  gsl_matrix_free(V_e);
  gsl_matrix_free(B);
  gsl_vector_free(beta);
  gsl_matrix_free(Vbeta);

  gsl_matrix_free(V_g_null);
  gsl_matrix_free(V_e_null);
  gsl_matrix_free(B_null);
  gsl_matrix_free(se_B_null1);
  gsl_matrix_free(se_B_null2);

  return;
}

void MVLMM::AnalyzePlinkGXE(const gsl_matrix *U, const gsl_vector *eval,
                            const gsl_matrix *UtW, const gsl_matrix *UtY,
                            const gsl_vector *env) {
  debug_msg("entering");
  string file_bed = file_bfile + ".bed";
  ifstream infile(file_bed.c_str(), ios::binary);
  if (!infile) {
    cout << "error reading bed file:" << file_bed << endl;
    return;
  }

  clock_t time_start = clock();
  time_UtX = 0;
  time_opt = 0;

  char ch[1];
  bitset<8> b;

  double logl_H0 = 0.0, logl_H1 = 0.0, p_wald = 0, p_lrt = 0, p_score = 0;
  double crt_a, crt_b, crt_c;
  int n_bit, n_miss, ci_total, ci_test;
  double geno, x_mean;
  size_t c = 0;
  size_t n_size = UtY->size1, d_size = UtY->size2, c_size = UtW->size2 + 2;
  size_t dc_size = d_size * (c_size + 1), v_size = d_size * (d_size + 1) / 2;

  // Large matrices for EM.
  gsl_matrix *U_hat = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *E_hat = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *OmegaU = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *OmegaE = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiY = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiBX = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiU = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *UltVehiE = gsl_matrix_alloc(d_size, n_size);

  // Large matrices for NR.
  // Each dxd block is H_k^{-1}.
  gsl_matrix *Hi_all = gsl_matrix_alloc(d_size, d_size * n_size);

  // Each column is H_k^{-1}y_k
  gsl_matrix *Hiy_all = gsl_matrix_alloc(d_size, n_size);

  // Each dcxdc block is x_k\otimes H_k^{-1}.
  gsl_matrix *xHi_all = gsl_matrix_alloc(dc_size, d_size * n_size);
  gsl_matrix *Hessian = gsl_matrix_alloc(v_size * 2, v_size * 2);

  gsl_vector *x = gsl_vector_alloc(n_size);

  gsl_matrix *Y = gsl_matrix_alloc(d_size, n_size);
  gsl_matrix *X = gsl_matrix_alloc(c_size + 1, n_size);
  gsl_matrix *V_g = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *V_e = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *B = gsl_matrix_alloc(d_size, c_size + 1);
  gsl_vector *beta = gsl_vector_alloc(d_size);
  gsl_matrix *Vbeta = gsl_matrix_alloc(d_size, d_size);

  // Null estimates for initial values.
  gsl_matrix *V_g_null = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *V_e_null = gsl_matrix_alloc(d_size, d_size);
  gsl_matrix *B_null = gsl_matrix_alloc(d_size, c_size + 1);
  gsl_matrix *se_B_null1 = gsl_matrix_alloc(d_size, c_size - 1);
  gsl_matrix *se_B_null2 = gsl_matrix_alloc(d_size, c_size);

  gsl_matrix_view X_sub1 = gsl_matrix_submatrix(X, 0, 0, c_size - 1, n_size);
  gsl_matrix_view B_sub1 = gsl_matrix_submatrix(B, 0, 0, d_size, c_size - 1);
  gsl_matrix_view xHi_all_sub1 = gsl_matrix_submatrix(
      xHi_all, 0, 0, d_size * (c_size - 1), d_size * n_size);

  gsl_matrix_view X_sub2 = gsl_matrix_submatrix(X, 0, 0, c_size, n_size);
  gsl_matrix_view B_sub2 = gsl_matrix_submatrix(B, 0, 0, d_size, c_size);
  gsl_matrix_view xHi_all_sub2 =
      gsl_matrix_submatrix(xHi_all, 0, 0, d_size * c_size, d_size * n_size);

  gsl_matrix_transpose_memcpy(Y, UtY);

  gsl_matrix_view X_sub0 = gsl_matrix_submatrix(X, 0, 0, c_size - 2, n_size);
  gsl_matrix_transpose_memcpy(&X_sub0.matrix, UtW);
  gsl_vector_view X_row0 = gsl_matrix_row(X, c_size - 2);
  gsl_blas_dgemv(CblasTrans, 1.0, U, env, 0.0, &X_row0.vector);

  gsl_vector_view X_row1 = gsl_matrix_row(X, c_size - 1);
  gsl_vector_set_zero(&X_row1.vector);
  gsl_vector_view X_row2 = gsl_matrix_row(X, c_size);
  gsl_vector_set_zero(&X_row2.vector);

  gsl_vector_view B_col1 = gsl_matrix_column(B, c_size - 1);
  gsl_vector_set_zero(&B_col1.vector);
  gsl_vector_view B_col2 = gsl_matrix_column(B, c_size);
  gsl_vector_set_zero(&B_col2.vector);

  MphInitial(em_iter, em_prec, nr_iter, nr_prec, eval, &X_sub1.matrix, Y, l_min,
             l_max, n_region, V_g, V_e, &B_sub1.matrix);

  logl_H0 = MphEM('R', em_iter, em_prec, eval, &X_sub1.matrix, Y, U_hat, E_hat,
                  OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g,
                  V_e, &B_sub1.matrix);
  logl_H0 = MphNR('R', nr_iter, nr_prec, eval, &X_sub1.matrix, Y, Hi_all,
                  &xHi_all_sub1.matrix, Hiy_all, V_g, V_e, Hessian, crt_a,
                  crt_b, crt_c);
  MphCalcBeta(eval, &X_sub1.matrix, Y, V_g, V_e, UltVehiY, &B_sub1.matrix,
              se_B_null1);

  c = 0;
  Vg_remle_null.clear();
  Ve_remle_null.clear();
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = i; j < d_size; j++) {
      Vg_remle_null.push_back(gsl_matrix_get(V_g, i, j));
      Ve_remle_null.push_back(gsl_matrix_get(V_e, i, j));
      VVg_remle_null.push_back(gsl_matrix_get(Hessian, c, c));
      VVe_remle_null.push_back(gsl_matrix_get(Hessian, c + v_size, c + v_size));
      c++;
    }
  }
  beta_remle_null.clear();
  se_beta_remle_null.clear();
  for (size_t i = 0; i < se_B_null1->size1; i++) {
    for (size_t j = 0; j < se_B_null1->size2; j++) {
      beta_remle_null.push_back(gsl_matrix_get(B, i, j));
      se_beta_remle_null.push_back(gsl_matrix_get(se_B_null1, i, j));
    }
  }
  logl_remle_H0 = logl_H0;

  cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  cout.precision(4);
  cout << "REMLE estimate for Vg in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_g, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Vg): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c, c)) << "\t";
    }
    cout << endl;
  }
  cout << "REMLE estimate for Ve in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_e, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Ve): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c + v_size, c + v_size)) << "\t";
    }
    cout << endl;
  }
  cout << "REMLE likelihood = " << logl_H0 << endl;

  logl_H0 = MphEM('L', em_iter, em_prec, eval, &X_sub1.matrix, Y, U_hat, E_hat,
                  OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g,
                  V_e, &B_sub1.matrix);
  logl_H0 = MphNR('L', nr_iter, nr_prec, eval, &X_sub1.matrix, Y, Hi_all,
                  &xHi_all_sub1.matrix, Hiy_all, V_g, V_e, Hessian, crt_a,
                  crt_b, crt_c);
  MphCalcBeta(eval, &X_sub1.matrix, Y, V_g, V_e, UltVehiY, &B_sub1.matrix,
              se_B_null1);

  c = 0;
  Vg_mle_null.clear();
  Ve_mle_null.clear();
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = i; j < d_size; j++) {
      Vg_mle_null.push_back(gsl_matrix_get(V_g, i, j));
      Ve_mle_null.push_back(gsl_matrix_get(V_e, i, j));
      VVg_mle_null.push_back(gsl_matrix_get(Hessian, c, c));
      VVe_mle_null.push_back(gsl_matrix_get(Hessian, c + v_size, c + v_size));
      c++;
    }
  }
  beta_mle_null.clear();
  se_beta_mle_null.clear();
  for (size_t i = 0; i < se_B_null1->size1; i++) {
    for (size_t j = 0; j < se_B_null1->size2; j++) {
      beta_mle_null.push_back(gsl_matrix_get(B, i, j));
      se_beta_mle_null.push_back(gsl_matrix_get(se_B_null1, i, j));
    }
  }
  logl_mle_H0 = logl_H0;

  cout << "MLE estimate for Vg in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_g, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Vg): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c, c)) << "\t";
    }
    cout << endl;
  }
  cout << "MLE estimate for Ve in the null model: " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << gsl_matrix_get(V_e, i, j) << "\t";
    }
    cout << endl;
  }
  cout << "se(Ve): " << endl;
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = 0; j <= i; j++) {
      c = GetIndex(i, j, d_size);
      cout << safe_sqrt(gsl_matrix_get(Hessian, c + v_size, c + v_size)) << "\t";
    }
    cout << endl;
  }
  cout << "MLE likelihood = " << logl_H0 << endl;

  vector<double> v_beta, v_Vg, v_Ve, v_Vbeta;
  for (size_t i = 0; i < d_size; i++) {
    v_beta.push_back(0.0);
  }
  for (size_t i = 0; i < d_size; i++) {
    for (size_t j = i; j < d_size; j++) {
      v_Vg.push_back(0.0);
      v_Ve.push_back(0.0);
      v_Vbeta.push_back(0.0);
    }
  }

  gsl_matrix_memcpy(V_g_null, V_g);
  gsl_matrix_memcpy(V_e_null, V_e);
  gsl_matrix_memcpy(B_null, B);

  // Start reading genotypes and analyze.
  // Calculate n_bit and c, the number of bit for each SNP.
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

  for (vector<SNPINFO>::size_type t = 0; t < snpInfo.size(); ++t) {
    if (t % d_pace == 0 || t == snpInfo.size() - 1) {
      ProgressBar("Reading SNPs", t, snpInfo.size() - 1);
    }
    if (indicator_snp[t] == 0) {
      continue;
    }

    // n_bit, and 3 is the number of magic numbers.
    infile.seekg(t * n_bit + 3);

    // Read genotypes.
    x_mean = 0.0;
    n_miss = 0;
    ci_total = 0;
    ci_test = 0;
    for (int i = 0; i < n_bit; ++i) {
      infile.read(ch, 1);
      b = ch[0];

      // Minor allele homozygous: 2.0; major: 0.0.
      for (size_t j = 0; j < 4; ++j) {

        if ((i == (n_bit - 1)) && ci_total == (int)ni_total) {
          break;
        }
        if (indicator_idv[ci_total] == 0) {
          ci_total++;
          continue;
        }

        if (b[2 * j] == 0) {
          if (b[2 * j + 1] == 0) {
            gsl_vector_set(x, ci_test, 2);
            x_mean += 2.0;
          } else {
            gsl_vector_set(x, ci_test, 1);
            x_mean += 1.0;
          }
        } else {
          if (b[2 * j + 1] == 1) {
            gsl_vector_set(x, ci_test, 0);
          } else {
            gsl_vector_set(x, ci_test, -9);
            n_miss++;
          }
        }

        ci_total++;
        ci_test++;
      }
    }

    x_mean /= (double)(ni_test - n_miss);

    for (size_t i = 0; i < ni_test; ++i) {
      geno = gsl_vector_get(x, i);
      if (geno == -9) {
        gsl_vector_set(x, i, x_mean);
        geno = x_mean;
      }
      if (x_mean > 1) {
        gsl_vector_set(x, i, 2 - geno);
      }
    }

    // Calculate statistics.
    time_start = clock();
    gsl_blas_dgemv(CblasTrans, 1.0, U, x, 0.0, &X_row1.vector);
    gsl_vector_mul(x, env);
    gsl_blas_dgemv(CblasTrans, 1.0, U, x, 0.0, &X_row2.vector);
    time_UtX += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    // Initial values.
    gsl_matrix_memcpy(V_g, V_g_null);
    gsl_matrix_memcpy(V_e, V_e_null);
    gsl_matrix_memcpy(B, B_null);

    if (a_mode == 2 || a_mode == 3 || a_mode == 4) {
      if (a_mode == 3 || a_mode == 4) {
        logl_H0 = MphEM('R', em_iter / 10, em_prec * 10, eval, &X_sub2.matrix,
                        Y, U_hat, E_hat, OmegaU, OmegaE, UltVehiY, UltVehiBX,
                        UltVehiU, UltVehiE, V_g, V_e, &B_sub2.matrix);
        logl_H0 = MphNR('R', nr_iter / 10, nr_prec * 10, eval, &X_sub2.matrix,
                        Y, Hi_all, &xHi_all_sub2.matrix, Hiy_all, V_g, V_e,
                        Hessian, crt_a, crt_b, crt_c);
        MphCalcBeta(eval, &X_sub2.matrix, Y, V_g, V_e, UltVehiY, &B_sub2.matrix,
                    se_B_null2);
      }

      if (a_mode == 2 || a_mode == 4) {
        logl_H0 = MphEM('L', em_iter / 10, em_prec * 10, eval, &X_sub2.matrix,
                        Y, U_hat, E_hat, OmegaU, OmegaE, UltVehiY, UltVehiBX,
                        UltVehiU, UltVehiE, V_g, V_e, &B_sub2.matrix);
        logl_H0 = MphNR('L', nr_iter / 10, nr_prec * 10, eval, &X_sub2.matrix,
                        Y, Hi_all, &xHi_all_sub2.matrix, Hiy_all, V_g, V_e,
                        Hessian, crt_a, crt_b, crt_c);
        MphCalcBeta(eval, &X_sub2.matrix, Y, V_g, V_e, UltVehiY, &B_sub2.matrix,
                    se_B_null2);
      }
    }

    time_start = clock();

    // 3 is before 1.
    if (a_mode == 3 || a_mode == 4) {
      p_score = MphCalcP(eval, &X_row2.vector, &X_sub2.matrix, Y, V_g_null,
                         V_e_null, UltVehiY, beta, Vbeta);

      if (p_score < p_nr && crt == 1) {
        logl_H1 = MphNR('R', 1, nr_prec * 10, eval, X, Y, Hi_all, xHi_all,
                        Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);
        p_score = PCRT(3, d_size, p_score, crt_a, crt_b, crt_c);
      }
    }

    if (a_mode == 2 || a_mode == 4) {
      logl_H1 = MphEM('L', em_iter / 10, em_prec * 10, eval, X, Y, U_hat, E_hat,
                      OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE,
                      V_g, V_e, B);

      // Calculate beta and Vbeta.
      p_lrt = MphCalcP(eval, &X_row2.vector, &X_sub2.matrix, Y, V_g, V_e,
                       UltVehiY, beta, Vbeta);
      p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), (double)d_size);

      if (p_lrt < p_nr) {
        logl_H1 =
            MphNR('L', nr_iter / 10, nr_prec * 10, eval, X, Y, Hi_all, xHi_all,
                  Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);

        // Calculate beta and Vbeta.
        p_lrt = MphCalcP(eval, &X_row2.vector, &X_sub2.matrix, Y, V_g, V_e,
                         UltVehiY, beta, Vbeta);
        p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), (double)d_size);
        if (crt == 1) {
          p_lrt = PCRT(2, d_size, p_lrt, crt_a, crt_b, crt_c);
        }
      }
    }

    if (a_mode == 1 || a_mode == 4) {
      logl_H1 = MphEM('R', em_iter / 10, em_prec * 10, eval, X, Y, U_hat, E_hat,
                      OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE,
                      V_g, V_e, B);
      p_wald = MphCalcP(eval, &X_row2.vector, &X_sub2.matrix, Y, V_g, V_e,
                        UltVehiY, beta, Vbeta);

      if (p_wald < p_nr) {
        logl_H1 =
            MphNR('R', nr_iter / 10, nr_prec * 10, eval, X, Y, Hi_all, xHi_all,
                  Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);
        p_wald = MphCalcP(eval, &X_row2.vector, &X_sub2.matrix, Y, V_g, V_e,
                          UltVehiY, beta, Vbeta);

        if (crt == 1) {
          p_wald = PCRT(1, d_size, p_wald, crt_a, crt_b, crt_c);
        }
      }
    }

    if (x_mean > 1) {
      gsl_vector_scale(beta, -1.0);
    }

    time_opt += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    // Store summary data.
    for (size_t i = 0; i < d_size; i++) {
      v_beta[i] = gsl_vector_get(beta, i);
    }

    c = 0;
    for (size_t i = 0; i < d_size; i++) {
      for (size_t j = i; j < d_size; j++) {
        v_Vg[c] = gsl_matrix_get(V_g, i, j);
        v_Ve[c] = gsl_matrix_get(V_e, i, j);
        v_Vbeta[c] = gsl_matrix_get(Vbeta, i, j);
        c++;
      }
    }

    MPHSUMSTAT SNPs = {v_beta, p_wald, p_lrt, p_score, v_Vg, v_Ve, v_Vbeta};
    sumStat.push_back(SNPs);
  }
  cout << endl;

  infile.close();
  infile.clear();

  gsl_matrix_free(U_hat);
  gsl_matrix_free(E_hat);
  gsl_matrix_free(OmegaU);
  gsl_matrix_free(OmegaE);
  gsl_matrix_free(UltVehiY);
  gsl_matrix_free(UltVehiBX);
  gsl_matrix_free(UltVehiU);
  gsl_matrix_free(UltVehiE);

  gsl_matrix_free(Hi_all);
  gsl_matrix_free(Hiy_all);
  gsl_matrix_free(xHi_all);
  gsl_matrix_free(Hessian);

  gsl_vector_free(x);

  gsl_matrix_free(Y);
  gsl_matrix_free(X);
  gsl_matrix_free(V_g);
  gsl_matrix_free(V_e);
  gsl_matrix_free(B);
  gsl_vector_free(beta);
  gsl_matrix_free(Vbeta);

  gsl_matrix_free(V_g_null);
  gsl_matrix_free(V_e_null);
  gsl_matrix_free(B_null);
  gsl_matrix_free(se_B_null1);
  gsl_matrix_free(se_B_null2);

  return;
}
