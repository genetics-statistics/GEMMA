/*
    Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright © 2011-2017, Xiang Zhou
    Copyright © 2017, Peter Carbonetto
    Copyright © 2017-2018 Pjotr Prins

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
#include <regex>
#include <stdio.h>
#include <stdlib.h>

#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_min.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_vector.h"

#include "gzstream.h"
#include "gemma.h"
#include "gemma_io.h"
#include "fastblas.h"
#include "lapack.h"
#include "lmm.h"
#include "mathfunc.h"

using namespace std;

void LMM::CopyFromParam(PARAM &cPar) {
  a_mode = cPar.a_mode;
  d_pace = cPar.d_pace;

  file_bfile = cPar.file_bfile;
  file_geno = cPar.file_geno;
  file_out = cPar.file_out;
  path_out = cPar.path_out;
  file_gene = cPar.file_gene;

  l_min = cPar.l_min;
  l_max = cPar.l_max;
  n_region = cPar.n_region;
  l_mle_null = cPar.l_mle_null;
  logl_mle_H0 = cPar.logl_mle_H0;

  time_UtX = 0.0;
  time_opt = 0.0;

  ni_total = cPar.ni_total;
  ns_total = cPar.ns_total;
  ni_test = cPar.ni_test;
  ns_test = cPar.ns_test;
  n_cvt = cPar.n_cvt;

  ng_total = cPar.ng_total;
  ng_test = 0;

  indicator_idv = cPar.indicator_idv;
  indicator_snp = cPar.indicator_snp;
  snpInfo = cPar.snpInfo;
  setGWASnps = cPar.setGWASnps;

  return;
}

void LMM::CopyToParam(PARAM &cPar) {
  cPar.time_UtX = time_UtX;
  cPar.time_opt = time_opt;

  cPar.ng_test = ng_test;

  return;
}

void LMM::WriteFiles() {
  string file_str;
  debug_msg("LMM::WriteFiles");

  file_str = path_out + "/" + file_out;
  file_str += ".assoc.txt";

  ofstream outfile(file_str.c_str(), ofstream::out);
  if (!outfile) {
    cout << "error writing file: " << file_str.c_str() << endl;
    return;
  }

  auto common_header = [&] () {
    if (a_mode != M_LMM2) {
      outfile << "beta" << "\t";
      outfile << "se" << "\t";
    }

    if (a_mode != M_LMM3 && a_mode != M_LMM9)
      outfile << "logl_H1" << "\t";

    switch(a_mode) {
    case M_LMM1:
      outfile << "l_remle" << "\t"
              << "p_wald" << endl;
      break;
    case M_LMM2:
    case M_LMM9:
      outfile << "l_mle" << "\t"
              << "p_lrt" << endl;
      break;
    case M_LMM3:
      outfile << "p_score" << endl;
      break;
    case M_LMM4:
      outfile << "l_remle" << "\t"
              << "l_mle" << "\t"
              << "p_wald" << "\t"
              << "p_lrt" << "\t"
              << "p_score" << endl;
      break;
    }
  };

  auto sumstats = [&] (SUMSTAT st) {
    outfile << scientific << setprecision(6);

    if (a_mode != M_LMM2) {
      outfile << st.beta << "\t";
      outfile << st.se << "\t";
    }

    if (a_mode != M_LMM3 && a_mode != M_LMM9)
      outfile << st.logl_H1 << "\t";

    switch(a_mode) {
    case M_LMM1:
      outfile << st.lambda_remle << "\t"
              << st.p_wald << endl;
      break;
    case M_LMM2:
    case M_LMM9:
      outfile << st.lambda_mle << "\t"
              << st.p_lrt << endl;
      break;
    case M_LMM3:
      outfile << st.p_score << endl;
      break;
    case M_LMM4:
      outfile << st.lambda_remle << "\t"
              << st.lambda_mle << "\t"
              << st.p_wald << "\t"
              << st.p_lrt << "\t"
              << st.p_score << endl;
      break;
    }
  };



  if (!file_gene.empty()) {
    outfile << "geneID" << "\t";

    common_header();

    for (vector<SUMSTAT>::size_type t = 0; t < sumStat.size(); ++t) {
      outfile << snpInfo[t].rs_number << "\t";
      sumstats(sumStat[t]);
    }
  } else {
    bool process_gwasnps = setGWASnps.size();

    outfile << "chr" << "\t"
            << "rs" << "\t"
            << "ps" << "\t"
            << "n_miss" << "\t"
            << "allele1" << "\t"
            << "allele0" << "\t"
            << "af" << "\t";

    common_header();

    size_t t = 0;
    for (size_t i = 0; i < snpInfo.size(); ++i) {
      if (indicator_snp[i] == 0)
        continue;
      auto snp = snpInfo[i].rs_number;
      if (process_gwasnps && setGWASnps.count(snp) == 0)
        continue;
      // cout << t << endl;
      outfile << snpInfo[i].chr << "\t" << snpInfo[i].rs_number << "\t"
              << snpInfo[i].base_position << "\t" << snpInfo[i].n_miss << "\t"
              << snpInfo[i].a_minor << "\t" << snpInfo[i].a_major << "\t"
              << fixed << setprecision(3) << snpInfo[i].maf << "\t";

      sumstats(sumStat[t]);
      t++;
    }
  }

  outfile.close();
  outfile.clear();
  return;
}

/*
   As explained in
   https://github.com/genetics-statistics/GEMMA/issues/94 CalcPab
   returns the Pab matrix. As described For pab, it stores all
   variables in the form of v_a P_p v_b. (Similarly, ppab stores all
   v_a P_p P_p v_b, while pppab stores all v_a P_p P_p P_p v_b. These
   quantities are defined according to the page 6 of this
   supplementary information
   http://xzlab.org/papers/2012_Zhou&Stephens_NG_SI.pdf).

   In the code, p, a, b are indexes: when p=n_cvt+1, P_p is P_x as in
   that supplementary information; when a=n_cvt+1, v_a=x; and when
   a=n_cvt+2, v_a=y.

   e_mode determines which model the algorithm is fitting: when
   e_mode==1, it computes all the above quantities for the alternative
   model (with the random effects term); when e_mode==0, it compute
   these quantities for the null model (without the random effects
   term).  Note that e==0 is only used here.

   These quantities were computed based on the initial GEMMA paper,
   and the goal is to finally compute y P_x y, y P_x P_xy,
   y P_x P_x P_x y and the few trace forms (section 3.1.4 on page 5 of
   the supplementary information). Sometimes I was wondering if we
   should compute all these final quantities directly, instead of
   through these complicated recursions. Direct computation may only
   make computation a little slower, but will make the code much
   easier to follow and easier to modify

   a typical call sends n_cvt a vector
   Hi_eval, a vector ab and a matrix Uab.

  CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);

(gdb) p Uab->size1
$1 = 247
(gdb) p n_cvt
$2 = 1
(gdb) p e_mode
$3 = 0
(gdb) p Uab->size2
$4 = 6
(gdb) p Hi_eval->size
$5 = 247
(gdb) p ab->size
$6 = 6
(gdb) p Pab->size1
$7 = 3
(gdb) p Pab->size2
$8 = 6

  Hi_eval [0..ind] x Uab [ind, n_index] x ab [n_index]

Iterating through a dataset Hi_eval differs and Uab (last row)
 */

void CalcPab(const size_t n_cvt, const size_t e_mode, const gsl_vector *Hi_eval,
             const gsl_matrix *Uab, const gsl_vector *unused, gsl_matrix *Pab) {

#if !defined NDEBUG
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2; // result size
  auto ni_test = Uab->size1; // inds
  assert(Uab->size1 == Hi_eval->size);
  assert(Uab->size2 == n_index);

  assert(Pab->size1 == n_cvt+2);
  assert(Pab->size2 == n_index);
  assert(Hi_eval->size == ni_test);
  // assert(ab->size == n_index);
#endif // DEBUG

  // compute Hi_eval (inds)  * Uab  (inds x n_index) * ab (n_index) and return in Pab (cvt x n_index).

  double p_ab = 0.0;

  write("CalcPab");
  write(Hi_eval,"Hi_eval");
  write(Uab,"Uab");
  // write(ab,"ab");
  if (is_check_mode())
    assert(!has_nan(Hi_eval));
  assert(!has_nan(Uab));
  // assert(!has_nan(ab));

  for (size_t p = 0; p <= n_cvt + 1; ++p) {        // p walks rows phenotypes + covariates
    for (size_t a = p + 1; a <= n_cvt + 2; ++a) {  // a walks cols in p+1..rest
      for (size_t b = a; b <= n_cvt + 2; ++b) {    // b in a..rest
        size_t index_ab = GetabIndex(a, b, n_cvt); // index in top half matrix, see above
        if (p == 0) { // fills row 0 for each a,b using dot product of Hi_eval . Uab(a)
          // cout << "p is 0 " << index_ab; // walk row 0
          gsl_vector_const_view Uab_col = gsl_matrix_const_column(Uab, index_ab); // get the column
          gsl_blas_ddot(Hi_eval, &Uab_col.vector, &p_ab); // dot product with H_eval
          if (e_mode != 0) { // if not null model (defunct right now)
            if (! is_legacy_mode()) assert(false); // disabling to see when it is used; allow with legacy mode
            p_ab = gsl_vector_get(unused, index_ab) - p_ab; // was ab
          }
          // cout << p << "r, index_ab " << index_ab << ":" << p_ab << endl;
          gsl_matrix_set(Pab, 0, index_ab, p_ab);
          write(Pab,"Pab int");
        } else {
          // walk the rest of the upper triangle of the matrix (row 1..n). Cols jump with 2 at a time
          // cout << "a" << a << "b" << b << "p" << p << "n_cvt" << n_cvt << endl;
          write(Pab,"Pab int");
          size_t index_aw = GetabIndex(a, p, n_cvt);
          size_t index_bw = GetabIndex(b, p, n_cvt);
          size_t index_ww = GetabIndex(p, p, n_cvt);

          // auto rows = Pab->size1; // n_cvt+2
          double ps_ab = gsl_matrix_safe_get(Pab, p - 1, index_ab);
          double ps_aw = gsl_matrix_safe_get(Pab, p - 1, index_aw);
          double ps_bw = gsl_matrix_safe_get(Pab, p - 1, index_bw);
          double ps_ww = gsl_matrix_safe_get(Pab, p - 1, index_ww);

          // cout << "unsafe " << p-1 << "," << index_ww << ":" << gsl_matrix_get(Pab,p-1,index_ww) << endl;
          // if (is_check_mode() || is_debug_mode()) assert(ps_ww != 0.0);
          if (ps_ww != 0)
            p_ab = ps_ab - ps_aw * ps_bw / ps_ww;
          else
            p_ab = ps_ab;

          // cout << "set " << p << "r, index_ab " << index_ab << "c: " << p_ab << endl;
          gsl_matrix_set(Pab, p, index_ab, p_ab);
        }
      }
    }
  }
  write(Pab,"Pab");
  // if (is_strict_mode() && (has_nan(Uab) || has_nan(Pab) || has_nan(Hi_eval)))
  //   exit(2);
  return;
}

void CalcPPab(const size_t n_cvt, const size_t e_mode,
              const gsl_vector *HiHi_eval, const gsl_matrix *Uab,
              const gsl_vector *unused_ab, const gsl_matrix *Pab, gsl_matrix *PPab) {
  size_t index_ab, index_aw, index_bw, index_ww;
  double p2_ab;
  double ps2_ab, ps_aw, ps_bw, ps_ww, ps2_aw, ps2_bw, ps2_ww;

  write("CalcPPab");
  write(HiHi_eval,"Hi_eval");
  write(Uab,"Uab");
  // write(ab,"ab");

  for (size_t p = 0; p <= n_cvt + 1; ++p) {
    for (size_t a = p + 1; a <= n_cvt + 2; ++a) {
      for (size_t b = a; b <= n_cvt + 2; ++b) {
        index_ab = GetabIndex(a, b, n_cvt);
        if (p == 0) {
          gsl_vector_const_view Uab_col =
              gsl_matrix_const_column(Uab, index_ab);
          gsl_blas_ddot(HiHi_eval, &Uab_col.vector, &p2_ab);
          if (e_mode != 0) {
            assert(false);
            p2_ab = p2_ab - gsl_vector_get(unused_ab, index_ab) +
                    2.0 * gsl_matrix_safe_get(Pab, 0, index_ab);
          }
          gsl_matrix_set(PPab, 0, index_ab, p2_ab);
        } else {
          index_aw = GetabIndex(a, p, n_cvt);
          index_bw = GetabIndex(b, p, n_cvt);
          index_ww = GetabIndex(p, p, n_cvt);

          ps2_ab = gsl_matrix_safe_get(PPab, p - 1, index_ab);
          ps_aw = gsl_matrix_safe_get(Pab, p - 1, index_aw);
          ps_bw = gsl_matrix_safe_get(Pab, p - 1, index_bw);
          ps_ww = gsl_matrix_safe_get(Pab, p - 1, index_ww);
          ps2_aw = gsl_matrix_safe_get(PPab, p - 1, index_aw);
          ps2_bw = gsl_matrix_safe_get(PPab, p - 1, index_bw);
          ps2_ww = gsl_matrix_safe_get(PPab, p - 1, index_ww);

          // if (is_check_mode() || is_debug_mode()) assert(ps_ww != 0.0);
          if (ps_ww != 0) {
            p2_ab = ps2_ab + ps_aw * ps_bw * ps2_ww / (ps_ww * ps_ww);
            p2_ab -= (ps_aw * ps2_bw + ps_bw * ps2_aw) / ps_ww;
          }
          else {
            p2_ab = ps2_ab;
          }
          gsl_matrix_set(PPab, p, index_ab, p2_ab);

        }
      }
    }
  }
  write(PPab,"PPab");
  // if (is_strict_mode() && (has_nan(Uab) || has_nan(PPab) || has_nan(HiHi_eval)))
  //   exit(2);
  return;
}

void CalcPPPab(const size_t n_cvt, const size_t e_mode,
               const gsl_vector *HiHiHi_eval, const gsl_matrix *Uab,
               const gsl_vector *unused_ab, const gsl_matrix *Pab,
               const gsl_matrix *PPab, gsl_matrix *PPPab) {
  size_t index_ab, index_aw, index_bw, index_ww;
  double p3_ab;
  double ps3_ab, ps_aw, ps_bw, ps_ww, ps2_aw, ps2_bw, ps2_ww, ps3_aw, ps3_bw,
      ps3_ww;
  write("CalcPPPab");
  write(HiHiHi_eval,"HiHiHi_eval");
  write(Uab,"Uab");

  for (size_t p = 0; p <= n_cvt + 1; ++p) {
    for (size_t a = p + 1; a <= n_cvt + 2; ++a) {
      for (size_t b = a; b <= n_cvt + 2; ++b) {
        index_ab = GetabIndex(a, b, n_cvt);
        if (p == 0) {
          gsl_vector_const_view Uab_col =
              gsl_matrix_const_column(Uab, index_ab);
          gsl_blas_ddot(HiHiHi_eval, &Uab_col.vector, &p3_ab);
          if (e_mode != 0) {
            assert(false);
            p3_ab = gsl_vector_get(unused_ab, index_ab) - p3_ab +
                    3.0 * gsl_matrix_get(PPab, 0, index_ab) -
                    3.0 * gsl_matrix_get(Pab, 0, index_ab);
          }
          gsl_matrix_set(PPPab, 0, index_ab, p3_ab);
        } else {
          index_aw = GetabIndex(a, p, n_cvt);
          index_bw = GetabIndex(b, p, n_cvt);
          index_ww = GetabIndex(p, p, n_cvt);

          ps3_ab = gsl_matrix_safe_get(PPPab, p - 1, index_ab);
          ps_aw = gsl_matrix_safe_get(Pab, p - 1, index_aw);
          ps_bw = gsl_matrix_safe_get(Pab, p - 1, index_bw);
          ps_ww = gsl_matrix_safe_get(Pab, p - 1, index_ww);
          ps2_aw = gsl_matrix_safe_get(PPab, p - 1, index_aw);
          ps2_bw = gsl_matrix_safe_get(PPab, p - 1, index_bw);
          ps2_ww = gsl_matrix_safe_get(PPab, p - 1, index_ww);
          ps3_aw = gsl_matrix_safe_get(PPPab, p - 1, index_aw);
          ps3_bw = gsl_matrix_safe_get(PPPab, p - 1, index_bw);
          ps3_ww = gsl_matrix_safe_get(PPPab, p - 1, index_ww);

          // if (is_check_mode() || is_debug_mode()) assert(ps_ww != 0.0);
          if (ps_ww != 0) {
            p3_ab = ps3_ab -
              ps_aw * ps_bw * ps2_ww * ps2_ww / (ps_ww * ps_ww * ps_ww);
            p3_ab -= (ps_aw * ps3_bw + ps_bw * ps3_aw + ps2_aw * ps2_bw) / ps_ww;
            p3_ab += (ps_aw * ps2_bw * ps2_ww + ps_bw * ps2_aw * ps2_ww +
                      ps_aw * ps_bw * ps3_ww) /
              (ps_ww * ps_ww);
          }
          else {
            p3_ab = ps3_ab;
          }
          gsl_matrix_set(PPPab, p, index_ab, p3_ab);
        }
      }
    }
  }
  write(PPPab,"PPPab");
  // if (is_strict_mode() && (has_nan(Uab) || has_nan(PPPab) || has_nan(HiHiHi_eval)))
  //  exit(2);
  return;
}

double LogL_f(double l, void *params) {
  FUNC_PARAM *p = (FUNC_PARAM *)params;
  size_t n_cvt = p->n_cvt;
  size_t ni_test = p->ni_test;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  size_t nc_total;
  if (p->calc_null == true) {
    nc_total = n_cvt;
  } else {
    nc_total = n_cvt + 1;
  }

  double f = 0.0, logdet_h = 0.0, d;
  size_t index_yy;

  gsl_matrix *Pab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_vector *Hi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *v_temp = gsl_vector_safe_alloc((p->eval)->size);

  gsl_vector_safe_memcpy(v_temp, p->eval);
  gsl_vector_scale(v_temp, l);
  if (p->e_mode == 0) {
    gsl_vector_set_all(Hi_eval, 1.0);
  } else {
    gsl_vector_safe_memcpy(Hi_eval, v_temp);
  }
  gsl_vector_add_constant(v_temp, 1.0);
  gsl_vector_div(Hi_eval, v_temp);

  for (size_t i = 0; i < (p->eval)->size; ++i) {
    d = gsl_vector_get(v_temp, i);
    logdet_h += safe_log(fabs(d));
  }

  CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);

  double c =
      0.5 * (double)ni_test * (safe_log((double)ni_test) - safe_log(2 * M_PI) - 1.0);

  index_yy = GetabIndex(n_cvt + 2, n_cvt + 2, n_cvt);
  double P_yy = gsl_matrix_safe_get(Pab, nc_total, index_yy);

  if (is_check_mode() || is_debug_mode()) {
    // cerr << "P_yy is" << P_yy << endl;
    assert(!is_nan(P_yy));
    assert(P_yy > 0);
  }
  f = c - 0.5 * logdet_h - 0.5 * (double)ni_test * safe_log(P_yy);
  if (is_check_mode() || is_debug_mode()) {
    assert(!is_nan(f));
  }
  gsl_matrix_free(Pab); // FIXME
  gsl_vector_safe_free(Hi_eval);
  gsl_vector_safe_free(v_temp);
  return f;
}

double LogL_dev1(double l, void *params) {
  FUNC_PARAM *p = (FUNC_PARAM *)params;
  size_t n_cvt = p->n_cvt;
  size_t ni_test = p->ni_test;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2; // represents top half of covariate matrix

  size_t nc_total;
  if (p->calc_null == true) {
    nc_total = n_cvt;
  } else {
    nc_total = n_cvt + 1;
  }

  double dev1 = 0.0, trace_Hi = 0.0;
  size_t index_yy;

  gsl_matrix *Pab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_matrix *PPab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_vector *Hi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *HiHi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *v_temp = gsl_vector_safe_alloc((p->eval)->size);

  gsl_vector_safe_memcpy(v_temp, p->eval);
  gsl_vector_scale(v_temp, l);
  if (p->e_mode == 0) {
    gsl_vector_set_all(Hi_eval, 1.0);
  } else {
    gsl_vector_safe_memcpy(Hi_eval, v_temp);
  }
  gsl_vector_add_constant(v_temp, 1.0);
  gsl_vector_div(Hi_eval, v_temp);

  gsl_vector_safe_memcpy(HiHi_eval, Hi_eval);
  gsl_vector_mul(HiHi_eval, Hi_eval);

  gsl_vector_set_all(v_temp, 1.0);
  gsl_blas_ddot(Hi_eval, v_temp, &trace_Hi);

  if (p->e_mode != 0) {
    trace_Hi = (double)ni_test - trace_Hi;
  }

  /*
(gdb) p Uab->size1
$1 = 247
(gdb) p n_cvt
$2 = 1
(gdb) p e_mode
$3 = 0
(gdb) p Uab->size2
$4 = 6
(gdb) p Hi_eval->size
$5 = 247
(gdb) p ab->size
$6 = 6
(gdb) p Pab->size1
$7 = 3
(gdb) p Pab->size2
$8 = 6
  */

#if !defined NDEBUG
  auto Uab = p->Uab;
  auto ab = p->ab;
  assert(n_index == (n_cvt + 2 + 1) * (n_cvt + 2) / 2);
  assert(Uab->size1 == ni_test);
  assert(Uab->size2 == n_index); // n_cvt == 1 -> n_index == 6?

  assert(Pab->size1 == n_cvt+2);
  assert(Pab->size2 == n_index);

  assert(ab->size == n_index);

  assert(p->e_mode == 0);
  assert(Hi_eval->size == ni_test);
#endif // DEBUG

  CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);
  CalcPPab(n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);

  double trace_HiK = ((double)ni_test - trace_Hi) / l;

  index_yy = GetabIndex(n_cvt + 2, n_cvt + 2, n_cvt);

  double P_yy = gsl_matrix_safe_get(Pab, nc_total, index_yy);
  double PP_yy = gsl_matrix_safe_get(PPab, nc_total, index_yy);
  double yPKPy = (P_yy - PP_yy) / l;
  dev1 = -0.5 * trace_HiK + 0.5 * (double)ni_test * yPKPy / P_yy;

  gsl_matrix_free(Pab);   // FIXME: may contain NaN
  gsl_matrix_free(PPab);  // FIXME: may contain NaN
  gsl_vector_safe_free(Hi_eval);
  gsl_vector_safe_free(HiHi_eval);
  gsl_vector_safe_free(v_temp);

  return dev1;
}

double LogL_dev2(double l, void *params) {
  FUNC_PARAM *p = (FUNC_PARAM *)params;
  size_t n_cvt = p->n_cvt;
  size_t ni_test = p->ni_test;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  size_t nc_total;
  if (p->calc_null == true) {
    nc_total = n_cvt;
  } else {
    nc_total = n_cvt + 1;
  }

  double dev2 = 0.0, trace_Hi = 0.0, trace_HiHi = 0.0;
  size_t index_yy;

  gsl_matrix *Pab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_matrix *PPab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_matrix *PPPab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_vector *Hi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *HiHi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *HiHiHi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *v_temp = gsl_vector_safe_alloc((p->eval)->size);

  gsl_vector_safe_memcpy(v_temp, p->eval);
  gsl_vector_scale(v_temp, l);
  if (p->e_mode == 0) {
    gsl_vector_set_all(Hi_eval, 1.0);
  } else {
    gsl_vector_safe_memcpy(Hi_eval, v_temp);
  }
  gsl_vector_add_constant(v_temp, 1.0);
  gsl_vector_div(Hi_eval, v_temp);

  gsl_vector_safe_memcpy(HiHi_eval, Hi_eval);
  gsl_vector_mul(HiHi_eval, Hi_eval);
  gsl_vector_safe_memcpy(HiHiHi_eval, HiHi_eval);
  gsl_vector_mul(HiHiHi_eval, Hi_eval);

  gsl_vector_set_all(v_temp, 1.0);
  gsl_blas_ddot(Hi_eval, v_temp, &trace_Hi);
  gsl_blas_ddot(HiHi_eval, v_temp, &trace_HiHi);

  if (p->e_mode != 0) {
    trace_Hi = (double)ni_test - trace_Hi;
    trace_HiHi = 2 * trace_Hi + trace_HiHi - (double)ni_test;
  }

  CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);
  CalcPPab(n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);
  CalcPPPab(n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, Pab, PPab, PPPab);

  double trace_HiKHiK = ((double)ni_test + trace_HiHi - 2 * trace_Hi) / (l * l);

  index_yy = GetabIndex(n_cvt + 2, n_cvt + 2, n_cvt);
  double P_yy = gsl_matrix_safe_get(Pab, nc_total, index_yy);
  double PP_yy = gsl_matrix_safe_get(PPab, nc_total, index_yy);
  double PPP_yy = gsl_matrix_safe_get(PPPab, nc_total, index_yy);

  double yPKPy = (P_yy - PP_yy) / l;
  double yPKPKPy = (P_yy + PPP_yy - 2.0 * PP_yy) / (l * l);

  dev2 = 0.5 * trace_HiKHiK -
         0.5 * (double)ni_test * (2.0 * yPKPKPy * P_yy - yPKPy * yPKPy) /
             (P_yy * P_yy);

  gsl_matrix_free(Pab);  // FIXME
  gsl_matrix_free(PPab);
  gsl_matrix_free(PPPab);
  gsl_vector_safe_free(Hi_eval);
  gsl_vector_safe_free(HiHi_eval);
  gsl_vector_safe_free(HiHiHi_eval);
  gsl_vector_safe_free(v_temp);

  return dev2;
}

void LogL_dev12(double l, void *params, double *dev1, double *dev2) {
  FUNC_PARAM *p = (FUNC_PARAM *)params;
  size_t n_cvt = p->n_cvt;
  size_t ni_test = p->ni_test;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  size_t nc_total;
  if (p->calc_null == true) {
    nc_total = n_cvt;
  } else {
    nc_total = n_cvt + 1;
  }

  double trace_Hi = 0.0, trace_HiHi = 0.0;
  size_t index_yy;

  gsl_matrix *Pab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_matrix *PPab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_matrix *PPPab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_vector *Hi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *HiHi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *HiHiHi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *v_temp = gsl_vector_safe_alloc((p->eval)->size);

  gsl_vector_safe_memcpy(v_temp, p->eval);
  gsl_vector_scale(v_temp, l);
  if (p->e_mode == 0) {
    gsl_vector_set_all(Hi_eval, 1.0);
  } else {
    gsl_vector_safe_memcpy(Hi_eval, v_temp);
  }
  gsl_vector_add_constant(v_temp, 1.0);
  gsl_vector_div(Hi_eval, v_temp);

  gsl_vector_safe_memcpy(HiHi_eval, Hi_eval);
  gsl_vector_mul(HiHi_eval, Hi_eval);
  gsl_vector_safe_memcpy(HiHiHi_eval, HiHi_eval);
  gsl_vector_mul(HiHiHi_eval, Hi_eval);

  gsl_vector_set_all(v_temp, 1.0);
  gsl_blas_ddot(Hi_eval, v_temp, &trace_Hi);
  gsl_blas_ddot(HiHi_eval, v_temp, &trace_HiHi);

  if (p->e_mode != 0) {
    trace_Hi = (double)ni_test - trace_Hi;
    trace_HiHi = 2 * trace_Hi + trace_HiHi - (double)ni_test;
  }

  CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);
  CalcPPab(n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);
  CalcPPPab(n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, Pab, PPab, PPPab);

  double trace_HiK = ((double)ni_test - trace_Hi) / l;
  double trace_HiKHiK = ((double)ni_test + trace_HiHi - 2 * trace_Hi) / (l * l);

  index_yy = GetabIndex(n_cvt + 2, n_cvt + 2, n_cvt);

  double P_yy = gsl_matrix_safe_get(Pab, nc_total, index_yy);
  double PP_yy = gsl_matrix_safe_get(PPab, nc_total, index_yy);
  double PPP_yy = gsl_matrix_safe_get(PPPab, nc_total, index_yy);

  double yPKPy = (P_yy - PP_yy) / l;
  double yPKPKPy = (P_yy + PPP_yy - 2.0 * PP_yy) / (l * l);

  *dev1 = -0.5 * trace_HiK + 0.5 * (double)ni_test * yPKPy / P_yy;
  *dev2 = 0.5 * trace_HiKHiK -
          0.5 * (double)ni_test * (2.0 * yPKPKPy * P_yy - yPKPy * yPKPy) /
              (P_yy * P_yy);

  gsl_matrix_free(Pab);   // FIXME: may contain NaN
  gsl_matrix_free(PPab);  // FIXME: may contain NaN
  gsl_matrix_free(PPPab); // FIXME: may contain NaN
  gsl_vector_safe_free(Hi_eval);
  gsl_vector_safe_free(HiHi_eval);
  gsl_vector_safe_free(HiHiHi_eval);
  gsl_vector_safe_free(v_temp);

  return;
}

double LogRL_f(double l, void *params) {
  FUNC_PARAM *p = (FUNC_PARAM *)params;
  size_t n_cvt = p->n_cvt;
  size_t ni_test = p->ni_test;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  double df;
  size_t nc_total;
  if (p->calc_null == true) {
    nc_total = n_cvt;
    df = (double)ni_test - (double)n_cvt;
  } else {
    nc_total = n_cvt + 1;
    df = (double)ni_test - (double)n_cvt - 1.0;
  }

  double f = 0.0, logdet_h = 0.0, logdet_hiw = 0.0, d;
  size_t index_ww;

  gsl_matrix *Pab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_matrix *Iab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_vector *Hi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *v_temp = gsl_vector_safe_alloc((p->eval)->size);

  gsl_vector_safe_memcpy(v_temp, p->eval);
  gsl_vector_scale(v_temp, l);
  if (p->e_mode == 0) {
    gsl_vector_set_all(Hi_eval, 1.0);
  } else {
    gsl_vector_safe_memcpy(Hi_eval, v_temp);
  }
  gsl_vector_add_constant(v_temp, 1.0);
  gsl_vector_div(Hi_eval, v_temp);

  for (size_t i = 0; i < (p->eval)->size; ++i) {
    d = gsl_vector_get(v_temp, i);
    logdet_h += safe_log(fabs(d));
  }

  CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);
  gsl_vector_set_all(v_temp, 1.0);
  CalcPab(n_cvt, p->e_mode, v_temp, p->Uab, p->ab, Iab);

  // Calculate |WHiW|-|WW|.
  logdet_hiw = 0.0;
  for (size_t i = 0; i < nc_total; ++i) {
    index_ww = GetabIndex(i + 1, i + 1, n_cvt);
    d = gsl_matrix_safe_get(Pab, i, index_ww);
    logdet_hiw += safe_log(d);
    d = gsl_matrix_safe_get(Iab, i, index_ww);
    logdet_hiw -= safe_log(d);
  }
  index_ww = GetabIndex(n_cvt + 2, n_cvt + 2, n_cvt);
  double P_yy = gsl_matrix_safe_get(Pab, nc_total, index_ww);

  double c = 0.5 * df * (safe_log(df) - safe_log(2 * M_PI) - 1.0);
  f = c - 0.5 * logdet_h - 0.5 * logdet_hiw - 0.5 * df * safe_log(P_yy);

  gsl_matrix_free(Pab);
  gsl_matrix_free(Iab); // contains NaN
  gsl_vector_safe_free(Hi_eval);
  gsl_vector_safe_free(v_temp);
  return f;
}

double LogRL_dev1(double l, void *params) {
  FUNC_PARAM *p = (FUNC_PARAM *)params;
  size_t n_cvt = p->n_cvt;
  size_t ni_test = p->ni_test;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  double df;
  size_t nc_total;
  if (p->calc_null == true) {
    nc_total = n_cvt;
    df = (double)ni_test - (double)n_cvt;
  } else {
    nc_total = n_cvt + 1;
    df = (double)ni_test - (double)n_cvt - 1.0;
  }

  double dev1 = 0.0, trace_Hi = 0.0;
  size_t index_ww;

  gsl_matrix *Pab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_matrix *PPab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_vector *Hi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *HiHi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *v_temp = gsl_vector_safe_alloc((p->eval)->size);

  // write(p->eval, "p->eval");
  gsl_vector_safe_memcpy(v_temp, p->eval); // initialize with eval
  gsl_vector_scale(v_temp, l);
  if (p->e_mode == 0) {
    gsl_vector_set_all(Hi_eval, 1.0);
  } else {
    gsl_vector_safe_memcpy(Hi_eval, v_temp);
  }
  gsl_vector_add_constant(v_temp, 1.0);
  gsl_vector_div(Hi_eval, v_temp);

  gsl_vector_safe_memcpy(HiHi_eval, Hi_eval);
  gsl_vector_mul(HiHi_eval, Hi_eval);

  gsl_vector_set_all(v_temp, 1.0);
  gsl_blas_ddot(Hi_eval, v_temp, &trace_Hi);

  if (p->e_mode != 0) {
    trace_Hi = (double)ni_test - trace_Hi;
  }

  write(p->eval, "p->eval2");
  write(p->ab, "p->ab");
  CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);
  CalcPPab(n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);

  // Calculate tracePK and trace PKPK.
  double trace_P = trace_Hi;
  double ps_ww, ps2_ww;
  for (size_t i = 0; i < nc_total; ++i) {
    index_ww = GetabIndex(i + 1, i + 1, n_cvt);
    ps_ww = gsl_matrix_safe_get(Pab, i, index_ww);
    ps2_ww = gsl_matrix_safe_get(PPab, i, index_ww);
    trace_P -= ps2_ww / ps_ww;
  }
  double trace_PK = (df - trace_P) / l;

  // Calculate yPKPy, yPKPKPy.
  index_ww = GetabIndex(n_cvt + 2, n_cvt + 2, n_cvt);
  double P_yy = gsl_matrix_safe_get(Pab, nc_total, index_ww);
  double PP_yy = gsl_matrix_safe_get(PPab, nc_total, index_ww);
  double yPKPy = (P_yy - PP_yy) / l;

  dev1 = -0.5 * trace_PK + 0.5 * df * yPKPy / P_yy;

  gsl_matrix_free(Pab);  // FIXME: may contain NaN
  gsl_matrix_free(PPab); // FIXME: may contain NaN
  gsl_vector_safe_free(Hi_eval);
  gsl_vector_safe_free(HiHi_eval);
  gsl_vector_safe_free(v_temp);

  return dev1;
}

double LogRL_dev2(double l, void *params) {
  FUNC_PARAM *p = (FUNC_PARAM *)params;
  size_t n_cvt = p->n_cvt;
  size_t ni_test = p->ni_test;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  double df;
  size_t nc_total;
  if (p->calc_null == true) {
    nc_total = n_cvt;
    df = (double)ni_test - (double)n_cvt;
  } else {
    nc_total = n_cvt + 1;
    df = (double)ni_test - (double)n_cvt - 1.0;
  }

  double dev2 = 0.0, trace_Hi = 0.0, trace_HiHi = 0.0;
  size_t index_ww;

  gsl_matrix *Pab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_matrix *PPab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_matrix *PPPab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_vector *Hi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *HiHi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *HiHiHi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *v_temp = gsl_vector_safe_alloc((p->eval)->size);

  gsl_vector_safe_memcpy(v_temp, p->eval);
  gsl_vector_scale(v_temp, l);
  if (p->e_mode == 0) {
    gsl_vector_set_all(Hi_eval, 1.0);
  } else {
    gsl_vector_safe_memcpy(Hi_eval, v_temp);
  }
  gsl_vector_add_constant(v_temp, 1.0);
  gsl_vector_div(Hi_eval, v_temp);

  gsl_vector_safe_memcpy(HiHi_eval, Hi_eval);
  gsl_vector_mul(HiHi_eval, Hi_eval);
  gsl_vector_safe_memcpy(HiHiHi_eval, HiHi_eval);
  gsl_vector_mul(HiHiHi_eval, Hi_eval);

  gsl_vector_set_all(v_temp, 1.0);
  gsl_blas_ddot(Hi_eval, v_temp, &trace_Hi);
  gsl_blas_ddot(HiHi_eval, v_temp, &trace_HiHi);

  if (p->e_mode != 0) {
    trace_Hi = (double)ni_test - trace_Hi;
    trace_HiHi = 2 * trace_Hi + trace_HiHi - (double)ni_test;
  }

  CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);
  CalcPPab(n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);
  CalcPPPab(n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, Pab, PPab, PPPab);

  // Calculate tracePK and trace PKPK.
  double trace_P = trace_Hi, trace_PP = trace_HiHi;
  double ps_ww, ps2_ww, ps3_ww;
  for (size_t i = 0; i < nc_total; ++i) {
    index_ww = GetabIndex(i + 1, i + 1, n_cvt);
    ps_ww = gsl_matrix_safe_get(Pab, i, index_ww);
    ps2_ww = gsl_matrix_safe_get(PPab, i, index_ww);
    ps3_ww = gsl_matrix_safe_get(PPPab, i, index_ww);
    trace_P -= ps2_ww / ps_ww;
    trace_PP += ps2_ww * ps2_ww / (ps_ww * ps_ww) - 2.0 * ps3_ww / ps_ww;
  }
  double trace_PKPK = (df + trace_PP - 2.0 * trace_P) / (l * l);

  // Calculate yPKPy, yPKPKPy.
  index_ww = GetabIndex(n_cvt + 2, n_cvt + 2, n_cvt);
  double P_yy = gsl_matrix_safe_get(Pab, nc_total, index_ww);
  double PP_yy = gsl_matrix_safe_get(PPab, nc_total, index_ww);
  double PPP_yy = gsl_matrix_safe_get(PPPab, nc_total, index_ww);
  double yPKPy = (P_yy - PP_yy) / l;
  double yPKPKPy = (P_yy + PPP_yy - 2.0 * PP_yy) / (l * l);

  dev2 = 0.5 * trace_PKPK -
         0.5 * df * (2.0 * yPKPKPy * P_yy - yPKPy * yPKPy) / (P_yy * P_yy);

  gsl_matrix_free(Pab);  // FIXME
  gsl_matrix_free(PPab);
  gsl_matrix_free(PPPab);
  gsl_vector_safe_free(Hi_eval);
  gsl_vector_safe_free(HiHi_eval);
  gsl_vector_safe_free(HiHiHi_eval);
  gsl_vector_safe_free(v_temp);

  return dev2;
}

void LogRL_dev12(double l, void *params, double *dev1, double *dev2) {
  FUNC_PARAM *p = (FUNC_PARAM *)params;
  size_t n_cvt = p->n_cvt;
  size_t ni_test = p->ni_test;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  double df;
  size_t nc_total;
  if (p->calc_null == true) {
    nc_total = n_cvt;
    df = (double)ni_test - (double)n_cvt;
  } else {
    nc_total = n_cvt + 1;
    df = (double)ni_test - (double)n_cvt - 1.0;
  }

  double trace_Hi = 0.0, trace_HiHi = 0.0;
  size_t index_ww;

  gsl_matrix *Pab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_matrix *PPab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_matrix *PPPab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_vector *Hi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *HiHi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *HiHiHi_eval = gsl_vector_safe_alloc((p->eval)->size);
  gsl_vector *v_temp = gsl_vector_safe_alloc((p->eval)->size);

  gsl_vector_safe_memcpy(v_temp, p->eval);
  gsl_vector_scale(v_temp, l);
  if (p->e_mode == 0) {
    gsl_vector_set_all(Hi_eval, 1.0);
  } else {
    gsl_vector_safe_memcpy(Hi_eval, v_temp);
  }
  gsl_vector_add_constant(v_temp, 1.0);
  gsl_vector_div(Hi_eval, v_temp);

  gsl_vector_safe_memcpy(HiHi_eval, Hi_eval);
  gsl_vector_mul(HiHi_eval, Hi_eval);
  gsl_vector_safe_memcpy(HiHiHi_eval, HiHi_eval);
  gsl_vector_mul(HiHiHi_eval, Hi_eval);

  gsl_vector_set_all(v_temp, 1.0);
  gsl_blas_ddot(Hi_eval, v_temp, &trace_Hi);
  gsl_blas_ddot(HiHi_eval, v_temp, &trace_HiHi);

  if (p->e_mode != 0) {
    trace_Hi = (double)ni_test - trace_Hi;
    trace_HiHi = 2 * trace_Hi + trace_HiHi - (double)ni_test;
  }

  CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);
  CalcPPab(n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);
  CalcPPPab(n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, Pab, PPab, PPPab);

  // Calculate tracePK and trace PKPK.
  double trace_P = trace_Hi, trace_PP = trace_HiHi;
  double ps_ww, ps2_ww, ps3_ww;
  for (size_t i = 0; i < nc_total; ++i) {
    index_ww = GetabIndex(i + 1, i + 1, n_cvt);
    ps_ww = gsl_matrix_safe_get(Pab, i, index_ww);
    ps2_ww = gsl_matrix_safe_get(PPab, i, index_ww);
    ps3_ww = gsl_matrix_safe_get(PPPab, i, index_ww);
    trace_P -= ps2_ww / ps_ww;
    trace_PP += ps2_ww * ps2_ww / (ps_ww * ps_ww) - 2.0 * ps3_ww / ps_ww;
  }
  double trace_PK = (df - trace_P) / l;
  double trace_PKPK = (df + trace_PP - 2.0 * trace_P) / (l * l);

  // Calculate yPKPy, yPKPKPy.
  index_ww = GetabIndex(n_cvt + 2, n_cvt + 2, n_cvt);
  double P_yy = gsl_matrix_safe_get(Pab, nc_total, index_ww);
  double PP_yy = gsl_matrix_safe_get(PPab, nc_total, index_ww);
  double PPP_yy = gsl_matrix_safe_get(PPPab, nc_total, index_ww);
  double yPKPy = (P_yy - PP_yy) / l;
  double yPKPKPy = (P_yy + PPP_yy - 2.0 * PP_yy) / (l * l);

  *dev1 = -0.5 * trace_PK + 0.5 * df * yPKPy / P_yy;
  *dev2 = 0.5 * trace_PKPK -
          0.5 * df * (2.0 * yPKPKPy * P_yy - yPKPy * yPKPy) / (P_yy * P_yy);

  gsl_matrix_free(Pab);  // FIXME
  gsl_matrix_free(PPab);
  gsl_matrix_free(PPPab);
  gsl_vector_safe_free(Hi_eval);
  gsl_vector_safe_free(HiHi_eval);
  gsl_vector_safe_free(HiHiHi_eval);
  gsl_vector_safe_free(v_temp);

  return;
}

void LMM::CalcRLWald(const double l, const FUNC_PARAM &params, double &beta,
                     double &se, double &p_wald) {
  size_t n_cvt = params.n_cvt;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  int df = (int)ni_test - (int)n_cvt - 1;

  gsl_matrix *Pab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_vector *Hi_eval = gsl_vector_safe_alloc(params.eval->size);
  gsl_vector *v_temp = gsl_vector_safe_alloc(params.eval->size);

  gsl_vector_safe_memcpy(v_temp, params.eval);
  gsl_vector_scale(v_temp, l);
  if (params.e_mode == 0) {
    gsl_vector_set_all(Hi_eval, 1.0);
  } else {
    gsl_vector_safe_memcpy(Hi_eval, v_temp);
  }
  gsl_vector_add_constant(v_temp, 1.0);
  gsl_vector_div(Hi_eval, v_temp);

  CalcPab(n_cvt, params.e_mode, Hi_eval, params.Uab, params.ab, Pab);

  size_t index_yy = GetabIndex(n_cvt + 2, n_cvt + 2, n_cvt);
  size_t index_xx = GetabIndex(n_cvt + 1, n_cvt + 1, n_cvt);
  size_t index_xy = GetabIndex(n_cvt + 2, n_cvt + 1, n_cvt);
  double P_yy = gsl_matrix_safe_get(Pab, n_cvt, index_yy);
  double P_xx = gsl_matrix_safe_get(Pab, n_cvt, index_xx);
  double P_xy = gsl_matrix_safe_get(Pab, n_cvt, index_xy);
  double Px_yy = gsl_matrix_safe_get(Pab, n_cvt + 1, index_yy);

  beta = P_xy / P_xx;
  double tau = (double)df / Px_yy;
  se = safe_sqrt(1.0 / (tau * P_xx));
  p_wald = gsl_cdf_fdist_Q((P_yy - Px_yy) * tau, 1.0, df);

  gsl_matrix_free(Pab);
  gsl_vector_safe_free(Hi_eval);
  gsl_vector_safe_free(v_temp);
  return;
}


void LMM::CalcRLScore(const double l, const FUNC_PARAM &params, double &beta,
                      double &se, double &p_score) {
  size_t n_cvt = params.n_cvt;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  int df = (int)ni_test - (int)n_cvt - 1;

  gsl_matrix *Pab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_vector *Hi_eval = gsl_vector_safe_alloc(params.eval->size);
  gsl_vector *v_temp = gsl_vector_safe_alloc(params.eval->size);

  gsl_vector_safe_memcpy(v_temp, params.eval);
  gsl_vector_scale(v_temp, l);
  if (params.e_mode == 0) {
    gsl_vector_set_all(Hi_eval, 1.0);
  } else {
    gsl_vector_safe_memcpy(Hi_eval, v_temp);
  }
  gsl_vector_add_constant(v_temp, 1.0);
  gsl_vector_div(Hi_eval, v_temp);

  CalcPab(n_cvt, params.e_mode, Hi_eval, params.Uab, params.ab, Pab);

  size_t index_yy = GetabIndex(n_cvt + 2, n_cvt + 2, n_cvt);
  size_t index_xx = GetabIndex(n_cvt + 1, n_cvt + 1, n_cvt);
  size_t index_xy = GetabIndex(n_cvt + 2, n_cvt + 1, n_cvt);
  double P_yy = gsl_matrix_safe_get(Pab, n_cvt, index_yy);
  double P_xx = gsl_matrix_safe_get(Pab, n_cvt, index_xx);
  double P_xy = gsl_matrix_safe_get(Pab, n_cvt, index_xy);
  double Px_yy = gsl_matrix_safe_get(Pab, n_cvt + 1, index_yy);

  beta = P_xy / P_xx;
  double tau = (double)df / Px_yy;
  se = safe_sqrt(1.0 / (tau * P_xx));

  p_score =
      gsl_cdf_fdist_Q((double)ni_test * P_xy * P_xy / (P_yy * P_xx), 1.0, df);

  gsl_matrix_free(Pab);
  gsl_vector_safe_free(Hi_eval);
  gsl_vector_safe_free(v_temp);
}

void CalcUab(const gsl_matrix *UtW, const gsl_vector *Uty, gsl_matrix *Uab) {
  size_t index_ab;
  size_t n_cvt = UtW->size2;

  // debug_msg("entering");

  gsl_vector *u_a = gsl_vector_safe_alloc(Uty->size);

  for (size_t a = 1; a <= n_cvt + 2; ++a) { // walk columns of pheno+cvt
    if (a == n_cvt + 1) {
      continue;
    }

    if (a == n_cvt + 2) {
      gsl_vector_safe_memcpy(u_a, Uty); // last column is phenotype
    } else {
      gsl_vector_const_view UtW_col = gsl_matrix_const_column(UtW, a - 1);
      gsl_vector_safe_memcpy(u_a, &UtW_col.vector);
    }

    for (size_t b = a; b >= 1; --b) { // back fill other columns
      if (b == n_cvt + 1) {
        continue;
      }

      index_ab = GetabIndex(a, b, n_cvt);
      gsl_vector_view Uab_col = gsl_matrix_column(Uab, index_ab);

      if (b == n_cvt + 2) {
        gsl_vector_safe_memcpy(&Uab_col.vector, Uty);
      } else {
        gsl_vector_const_view UtW_col = gsl_matrix_const_column(UtW, b - 1);
        gsl_vector_safe_memcpy(&Uab_col.vector, &UtW_col.vector);
      }

      gsl_vector_mul(&Uab_col.vector, u_a);
    }
    // cout << "a" << a << endl;
    write(Uab,"Uab iteration");
  }

  gsl_vector_safe_free(u_a);
  return;
}

void CalcUab(const gsl_matrix *UtW, const gsl_vector *Uty,
             const gsl_vector *Utx, gsl_matrix *Uab) {
  size_t index_ab;
  size_t n_cvt = UtW->size2;

  for (size_t b = 1; b <= n_cvt + 2; ++b) {
    index_ab = GetabIndex(n_cvt + 1, b, n_cvt);
    gsl_vector_view Uab_col = gsl_matrix_column(Uab, index_ab);

    if (b == n_cvt + 2) {
      gsl_vector_safe_memcpy(&Uab_col.vector, Uty);
    } else if (b == n_cvt + 1) {
      gsl_vector_safe_memcpy(&Uab_col.vector, Utx);
    } else {
      gsl_vector_const_view UtW_col = gsl_matrix_const_column(UtW, b - 1);
      gsl_vector_safe_memcpy(&Uab_col.vector, &UtW_col.vector);
    }

    gsl_vector_mul(&Uab_col.vector, Utx);
  }

  return;
}

void Calcab(const gsl_matrix *W, const gsl_vector *y, gsl_vector *ab) {
  write(W,"W");
  write(y,"y");

  gsl_vector_set_zero(ab); // not sure, but emulates v95 behaviour

  size_t n_cvt = W->size2;

  gsl_vector *v_a = gsl_vector_safe_alloc(y->size);
  gsl_vector *v_b = gsl_vector_safe_alloc(y->size);

  double d;

  for (size_t a = 1; a <= n_cvt + 2; ++a) {
    if (a == n_cvt + 1) {
      continue;
    }

    if (a == n_cvt + 2) {
      gsl_vector_safe_memcpy(v_a, y);
    } else {
      gsl_vector_const_view W_col = gsl_matrix_const_column(W, a - 1);
      gsl_vector_safe_memcpy(v_a, &W_col.vector);
    }
    write(v_a,"v_a");

    for (size_t b = a; b >= 1; --b) {
      if (b == n_cvt + 1) {
        continue;
      }

      auto index_ab = GetabIndex(a, b, n_cvt);

      if (b == n_cvt + 2) {
        gsl_vector_safe_memcpy(v_b, y);
      } else {
        gsl_vector_const_view W_col = gsl_matrix_const_column(W, b - 1);
        gsl_vector_safe_memcpy(v_b, &W_col.vector);
      }

      write(v_b,"v_b");
      gsl_blas_ddot(v_a, v_b, &d);
      gsl_vector_set(ab, index_ab, d);
    }
  }

  write(ab,"ab");

  gsl_vector_safe_free(v_a);
  gsl_vector_safe_free(v_b);
  return;
}

void Calcab(const gsl_matrix *W, const gsl_vector *y, const gsl_vector *x,
            gsl_vector *ab) {
  size_t index_ab;
  size_t n_cvt = W->size2;

  gsl_vector_set_zero(ab); // not sure, but emulates v95 behaviour

  double d;
  gsl_vector *v_b = gsl_vector_safe_alloc(y->size);

  for (size_t b = 1; b <= n_cvt + 2; ++b) {
    index_ab = GetabIndex(n_cvt + 1, b, n_cvt);

    if (b == n_cvt + 2) {
      gsl_vector_safe_memcpy(v_b, y);
    } else if (b == n_cvt + 1) {
      gsl_vector_safe_memcpy(v_b, x);
    } else {
      gsl_vector_const_view W_col = gsl_matrix_const_column(W, b - 1);
      gsl_vector_safe_memcpy(v_b, &W_col.vector);
    }

    gsl_blas_ddot(x, v_b, &d);
    gsl_vector_set(ab, index_ab, d);
  }

  gsl_vector_safe_free(v_b);
  return;
}

void LMM::AnalyzeGene(const gsl_matrix *U, const gsl_vector *eval,
                      const gsl_matrix *UtW, const gsl_vector *Utx,
                      const gsl_matrix *W, const gsl_vector *x) {
  debug_msg(file_gene);
  igzstream infile(file_gene.c_str(), igzstream::in);
  if (!infile) {
    cout << "error reading gene expression file:" << file_gene << endl;
    return;
  }

  clock_t time_start = clock();

  string line;
  char *ch_ptr;

  double lambda_mle = 0, lambda_remle = 0, beta = 0, se = 0, p_wald = 0;
  double p_lrt = 0, p_score = 0;
  double logl_H1 = 0.0, logl_H0 = 0.0, l_H0;
  int c_phen;
  string rs; // Gene id.
  double d;

  // Calculate basic quantities.
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  gsl_vector *y = gsl_vector_safe_alloc(U->size1);
  gsl_vector *Uty = gsl_vector_safe_alloc(U->size2);
  gsl_matrix *Uab = gsl_matrix_safe_alloc(U->size2, n_index);
  gsl_vector *ab = gsl_vector_safe_alloc(n_index);

  // Header.
  getline(infile, line);

  for (size_t t = 0; t < ng_total; t++) {
    safeGetline(infile, line).eof();
    if (t % d_pace == 0 || t == ng_total - 1) {
      ProgressBar("Performing Analysis", t, ng_total - 1);
    }
    ch_ptr = strtok_safe2((char *)line.c_str(), " , \t",file_gene.c_str());
    rs = ch_ptr;

    c_phen = 0;
    for (size_t i = 0; i < indicator_idv.size(); ++i) {
      ch_ptr = strtok_safe2(NULL, " , \t",file_gene.c_str());
      if (indicator_idv[i] == 0) {
        continue;
      }

      d = atof(ch_ptr);
      gsl_vector_set(y, c_phen, d);

      c_phen++;
    }

    time_start = clock();
    gsl_blas_dgemv(CblasTrans, 1.0, U, y, 0.0, Uty);
    time_UtX += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    // Calculate null.
    time_start = clock();

    gsl_matrix_set_zero(Uab);

    CalcUab(UtW, Uty, Uab);
    FUNC_PARAM param0 = {false, ni_test, n_cvt, eval, Uab, ab, 0};

    if (a_mode == M_LMM2 || a_mode == M_LMM3 || a_mode == M_LMM4) {
      CalcLambda('L', param0, l_min, l_max, n_region, l_H0, logl_H0);
    }

    // Calculate alternative.
    CalcUab(UtW, Uty, Utx, Uab);
    FUNC_PARAM param1 = {false, ni_test, n_cvt, eval, Uab, ab, 0};

    // 3 is before 1.
    if (a_mode == M_LMM3 || a_mode == M_LMM4) {
      CalcRLScore(l_H0, param1, beta, se, p_score);
    }

    if (a_mode == M_LMM1 || a_mode == M_LMM4) {
      CalcLambda('R', param1, l_min, l_max, n_region, lambda_remle, logl_H1);
      CalcRLWald(lambda_remle, param1, beta, se, p_wald);
    }

    if (a_mode == M_LMM2 || a_mode == M_LMM4) {
      CalcLambda('L', param1, l_min, l_max, n_region, lambda_mle, logl_H1);
      p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), 1);
    }

    time_opt += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    // Store summary data.
    SUMSTAT SNPs = {beta, se, lambda_remle, lambda_mle, p_wald, p_lrt, p_score, logl_H1};
    sumStat.push_back(SNPs);
  }
  cout << endl;

  gsl_vector_safe_free(y);
  gsl_vector_safe_free(Uty);
  gsl_matrix_safe_free(Uab);
  gsl_vector_free(ab); // unused

  infile.close();
  infile.clear();

  return;
}


void LMM::Analyze(std::function< SnpNameValues(size_t) >& fetch_snp,
                  const gsl_matrix *U, const gsl_vector *eval,
                  const gsl_matrix *UtW, const gsl_vector *Uty,
                  const gsl_matrix *W, const gsl_vector *y,
                  const set<string> gwasnps) {
  clock_t time_start = clock();

  // Subset/LOCO support
  bool process_gwasnps = gwasnps.size();
  if (process_gwasnps)
    debug_msg("Analyze subset of SNPs (LOCO)");

  // Calculate basic quantities.
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  const size_t inds = U->size1;
  enforce(inds == ni_test);
  gsl_vector *x = gsl_vector_safe_alloc(inds); // #inds
  gsl_vector *x_miss = gsl_vector_safe_alloc(inds);
  gsl_vector *Utx = gsl_vector_safe_alloc(U->size2);
  gsl_matrix *Uab = gsl_matrix_safe_alloc(U->size2, n_index);
  gsl_vector *ab = gsl_vector_safe_alloc(n_index);

  // Create a large matrix with LMM_BATCH_SIZE columns for batched processing
  // const size_t msize=(process_gwasnps ? 1 : LMM_BATCH_SIZE);
  const size_t msize = LMM_BATCH_SIZE;
  gsl_matrix *Xlarge = gsl_matrix_safe_alloc(inds, msize);
  gsl_matrix *UtXlarge = gsl_matrix_safe_alloc(inds, msize);
  enforce_msg(Xlarge && UtXlarge, "Xlarge memory check"); // just to be sure
  enforce(Xlarge->size1 == inds);
  gsl_matrix_set_zero(Xlarge);
  gsl_matrix_set_zero(Uab);
  CalcUab(UtW, Uty, Uab);

  // start reading genotypes and analyze
  size_t c = 0;

  auto batch_compute = [&](size_t l) { // using a C++ closure
    // Compute SNPs in batch, note the computations are independent per SNP
    // debug_msg("enter batch_compute");
    gsl_matrix_view Xlarge_sub = gsl_matrix_submatrix(Xlarge, 0, 0, inds, l);
    gsl_matrix_view UtXlarge_sub =
        gsl_matrix_submatrix(UtXlarge, 0, 0, inds, l);

    time_start = clock();
    fast_dgemm("T", "N", 1.0, U, &Xlarge_sub.matrix, 0.0,
                   &UtXlarge_sub.matrix);
    time_UtX += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    gsl_matrix_set_zero(Xlarge);
    for (size_t i = 0; i < l; i++) {
      // for every batch...
      gsl_vector_view UtXlarge_col = gsl_matrix_column(UtXlarge, i);
      gsl_vector_safe_memcpy(Utx, &UtXlarge_col.vector);

      CalcUab(UtW, Uty, Utx, Uab);

      time_start = clock();
      FUNC_PARAM param1 = {false, ni_test, n_cvt, eval, Uab, ab, 0};

      double lambda_mle = 0.0, lambda_remle = 0.0, beta = 0.0, se = 0.0, p_wald = 0.0;
      double p_lrt = 0.0, p_score = 0.0;
      double logl_H1 = 0.0;

      // 3 is before 1.
      if (a_mode == M_LMM3 || a_mode == M_LMM4 || a_mode == M_LMM9 ) {
        CalcRLScore(l_mle_null, param1, beta, se, p_score);
      }

      if (a_mode == M_LMM1 || a_mode == M_LMM4) {
        // for univariate a_mode is 1
        CalcLambda('R', param1, l_min, l_max, n_region, lambda_remle, logl_H1);
        CalcRLWald(lambda_remle, param1, beta, se, p_wald);
      }

      if (a_mode == M_LMM2 || a_mode == M_LMM9 || a_mode == M_LMM4) {
        CalcLambda('L', param1, l_min, l_max, n_region, lambda_mle, logl_H1);
        p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_mle_H0), 1);
      }

      time_opt += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

      // Store summary data.
      SUMSTAT SNPs = {beta,   se,    lambda_remle, lambda_mle,
                      p_wald, p_lrt, p_score, logl_H1};
      sumStat.push_back(SNPs);
    }
    // debug_msg("exit batch_compute");
  };

  const auto num_snps = indicator_snp.size();
  enforce_msg(num_snps > 0,"Zero SNPs to process - data corrupt?");
  if (num_snps < 50) {
    cerr << num_snps << " SNPs" << endl;
    warning_msg("very few SNPs processed");
  }
  const size_t progress_step = (num_snps/50>d_pace ? num_snps/50 : d_pace);

  for (size_t t = 0; t < num_snps; ++t) {
    if (t % progress_step == 0 || t == (num_snps - 1)) {
      ProgressBar("Reading SNPs", t, num_snps - 1);
    }
    if (indicator_snp[t] == 0)
      continue;

    auto tup = fetch_snp(t);
    auto snp = get<0>(tup);
    auto gs = get<1>(tup);

    // check whether SNP is included in gwasnps (used by LOCO)
    if (process_gwasnps && gwasnps.count(snp) == 0)
      continue;

    // drop missing idv and plug mean values for missing geno
    double x_total = 0.0; // sum genotype values to compute x_mean
    uint pos = 0;         // position in target vector
    uint n_miss = 0;
    gsl_vector_set_zero(x_miss);
    for (size_t i = 0; i < ni_total; ++i) {
      // get the genotypes per individual and compute stats per SNP
      if (indicator_idv[i] == 0) // skip individual
        continue;

      double geno = gs[i];
      if (isnan(geno)) {
        gsl_vector_set(x_miss, pos, 1.0);
        n_miss++;
      } else {
        gsl_vector_set(x, pos, geno);
        x_total += geno;
      }
      pos++;
    }
    enforce(pos == ni_test);

    const double x_mean = x_total/(double)(ni_test - n_miss);

    // plug x_mean back into missing values
    for (size_t i = 0; i < ni_test; ++i) {
      if (gsl_vector_get(x_miss, i) == 1.0) {
        gsl_vector_set(x, i, x_mean);
      }
    }

    /* this is what below GxE does
    for (size_t i = 0; i < ni_test; ++i) {
      auto geno = gsl_vector_get(x, i);
      if (std::isnan(geno)) {
        gsl_vector_set(x, i, x_mean);
        geno = x_mean;
      }
      if (x_mean > 1.0) {
        gsl_vector_set(x, i, 2 - geno);
      }
    }
    */
    enforce(x->size == ni_test);

    // copy genotype values for SNP into Xlarge cache
    gsl_vector_view Xlarge_col = gsl_matrix_column(Xlarge, c % msize);
    gsl_vector_safe_memcpy(&Xlarge_col.vector, x);
    c++; // count SNPs going in

    if (c % msize == 0) {
      batch_compute(msize);
    }
  }

  batch_compute(c % msize);
  ProgressBar("Reading SNPs", num_snps - 1, num_snps - 1);
  // cout << "Counted SNPs " << c << " sumStat " << sumStat.size() << endl;
  cout << endl;

  gsl_vector_safe_free(x);
  gsl_vector_safe_free(x_miss);
  gsl_vector_safe_free(Utx);
  gsl_matrix_safe_free(Uab);
  gsl_vector_free(ab); // unused

  gsl_matrix_safe_free(Xlarge);
  gsl_matrix_safe_free(UtXlarge);

}

void LMM::AnalyzeBimbam(const gsl_matrix *U, const gsl_vector *eval,
                        const gsl_matrix *UtW, const gsl_vector *Uty,
                        const gsl_matrix *W, const gsl_vector *y,
                        const set<string> gwasnps) {
  debug_msg(file_geno);
  auto infilen = file_geno.c_str();

  igzstream infile(infilen, igzstream::in);
  enforce_msg(infile, "error reading genotype file");
  size_t prev_line = 0;

  std::vector <double> gs;
  gs.resize(ni_total);

  // fetch_snp is a callback function for every SNP row
  std::function<SnpNameValues(size_t)>  fetch_snp = [&](size_t num) {
    string line;
    while (prev_line <= num) {
      // also read SNPs that were skipped
      safeGetline(infile, line);
      prev_line++;
    }
    char *ch_ptr = strtok_safe2((char *)line.c_str(), " , \t",infilen);
    // enforce_msg(ch_ptr, "Parsing BIMBAM genofile"); // ch_ptr should not be NULL

    auto snp = string(ch_ptr);
    ch_ptr = strtok_safe2(NULL, " , \t",infilen); // skip column
    ch_ptr = strtok_safe2(NULL, " , \t",infilen); // skip column

    gs.assign (ni_total,nan("")); // wipe values

    for (size_t i = 0; i < ni_total; ++i) {
      ch_ptr = strtok_safe2(NULL, " , \t",infilen);
      if (strcmp(ch_ptr, "NA") != 0) {
        gs[i] = atof(ch_ptr);
        if (is_strict_mode() && gs[i] == 0.0)
          enforce_is_float(std::string(ch_ptr)); // only allow for NA and valid numbers
      }
    }
    return std::make_tuple(snp,gs);
  };

  LMM::Analyze(fetch_snp,U,eval,UtW,Uty,W,y,gwasnps);

  infile.close();
  infile.clear();
}

#include "eigenlib.h"

void LMM::AnalyzePlink(const gsl_matrix *U, const gsl_vector *eval,
                       const gsl_matrix *UtW, const gsl_vector *Uty,
                       const gsl_matrix *W, const gsl_vector *y,
                       const set<string> gwasnps) {
  string file_bed = file_bfile + ".bed";
  debug_msg(file_bed);

  ifstream infile(file_bed.c_str(), ios::binary);
  enforce_msg(infile,"error reading genotype (.bed) file");

  clock_t time_start = clock();

  char ch[1];
  bitset<8> b;

  double lambda_mle = 0, lambda_remle = 0, beta = 0, se = 0, p_wald = 0;
  double p_lrt = 0, p_score = 0;
  double logl_H1 = 0.0;
  int n_bit, n_miss, ci_total, ci_test;
  double geno, x_mean;

  // Calculate basic quantities.
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  gsl_vector *x = gsl_vector_alloc(U->size1);
  gsl_vector *Utx = gsl_vector_alloc(U->size2);
  gsl_matrix *Uab = gsl_matrix_alloc(U->size2, n_index);
  gsl_vector *ab = gsl_vector_alloc(n_index);

  // Create a large matrix.
  size_t msize = LMM_BATCH_SIZE;
  gsl_matrix *Xlarge = gsl_matrix_alloc(U->size1, msize);
  gsl_matrix *UtXlarge = gsl_matrix_alloc(U->size1, msize);
  gsl_matrix_set_zero(Xlarge);

  gsl_matrix_set_zero(Uab);
  CalcUab(UtW, Uty, Uab);

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

  size_t c = 0, t_last = 0;
  for (size_t t = 0; t < snpInfo.size(); ++t) {
    if (indicator_snp[t] == 0)
      continue;
    t_last++;
  }
  for (vector<SNPINFO>::size_type t = 0; t < snpInfo.size(); ++t) {
    if (t % d_pace == 0 || t == snpInfo.size() - 1) {
      ProgressBar("Reading SNPs  ", t, snpInfo.size() - 1);
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
    }

    gsl_vector_view Xlarge_col = gsl_matrix_column(Xlarge, c % msize);
    gsl_vector_memcpy(&Xlarge_col.vector, x);
    c++;

    if (c % msize == 0 || c == t_last) {
      size_t l = 0;
      if (c % msize == 0) {
        l = msize;
      } else {
        l = c % msize;
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
        gsl_vector_memcpy(Utx, &UtXlarge_col.vector);

        CalcUab(UtW, Uty, Utx, Uab);

        time_start = clock();
        FUNC_PARAM param1 = {false, ni_test, n_cvt, eval, Uab, ab, 0};

        // 3 is before 1, for beta.
        if (a_mode == M_LMM3 || a_mode == M_LMM4) {
          CalcRLScore(l_mle_null, param1, beta, se, p_score);
        }

        if (a_mode == M_LMM1 || a_mode == M_LMM4) {
          CalcLambda('R', param1, l_min, l_max, n_region, lambda_remle,
                     logl_H1);
          CalcRLWald(lambda_remle, param1, beta, se, p_wald);
        }

        if (a_mode == M_LMM2 || a_mode == M_LMM4) {
          CalcLambda('L', param1, l_min, l_max, n_region, lambda_mle, logl_H1);
          p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_mle_H0), 1);
        }

        time_opt += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

        // Store summary data.
        SUMSTAT SNPs = {beta,   se,    lambda_remle, lambda_mle,
                        p_wald, p_lrt, p_score, logl_H1};
        sumStat.push_back(SNPs);
      }
    }
  }
  cout << endl;

  gsl_vector_free(x);
  gsl_vector_free(Utx);
  gsl_matrix_free(Uab);
  gsl_vector_free(ab);

  gsl_matrix_free(Xlarge);
  gsl_matrix_free(UtXlarge);

  infile.close();
  infile.clear();
}

void MatrixCalcLR(const gsl_matrix *U, const gsl_matrix *UtX,
                  const gsl_vector *Uty, const gsl_vector *K_eval,
                  const double l_min, const double l_max, const size_t n_region,
                  vector<pair<size_t, double>> &pos_loglr) {
  double logl_H0, logl_H1, log_lr, lambda0, lambda1;

  gsl_vector *w = gsl_vector_safe_alloc(Uty->size);
  gsl_matrix *Utw = gsl_matrix_safe_alloc(Uty->size, 1);
  gsl_matrix *Uab = gsl_matrix_safe_alloc(Uty->size, 6);
  gsl_vector *ab = gsl_vector_safe_alloc(6);

  gsl_vector_set_zero(ab);
  gsl_vector_set_all(w, 1.0);
  gsl_vector_view Utw_col = gsl_matrix_column(Utw, 0);
  gsl_blas_dgemv(CblasTrans, 1.0, U, w, 0.0, &Utw_col.vector);

  CalcUab(Utw, Uty, Uab);
  FUNC_PARAM param0 = {true, Uty->size, 1, K_eval, Uab, ab, 0};

  CalcLambda('L', param0, l_min, l_max, n_region, lambda0, logl_H0);

  for (size_t i = 0; i < UtX->size2; ++i) {
    gsl_vector_const_view UtX_col = gsl_matrix_const_column(UtX, i);
    CalcUab(Utw, Uty, &UtX_col.vector, Uab);
    FUNC_PARAM param1 = {false, UtX->size1, 1, K_eval, Uab, ab, 0};

    CalcLambda('L', param1, l_min, l_max, n_region, lambda1, logl_H1);
    log_lr = logl_H1 - logl_H0;

    pos_loglr.push_back(make_pair(i, log_lr));
  }

  gsl_vector_safe_free(w);
  gsl_matrix_safe_free(Utw);
  gsl_matrix_safe_free(Uab);
  gsl_vector_free(ab); // unused

  return;
}

void CalcLambda(const char func_name, FUNC_PARAM &params, const double l_min,
                const double l_max, const size_t n_region, double &lambda,
                double &logf) {
  if (func_name != 'R' && func_name != 'L' && func_name != 'r' &&
      func_name != 'l') {
    cout << "func_name only takes 'R' or 'L': 'R' for "
         << "log-restricted likelihood, 'L' for log-likelihood." << endl;
    return;
  }

  vector<pair<double, double>> lambda_lh;

  // Evaluate first-order derivates in different intervals.
  double lambda_l, lambda_h,
      lambda_interval = safe_log(l_max / l_min) / (double)n_region;
  double dev1_l, dev1_h, logf_l, logf_h;

  for (size_t i = 0; i < n_region; ++i) {
    lambda_l = l_min * exp(lambda_interval * i);
    lambda_h = l_min * exp(lambda_interval * (i + 1.0));

    if (func_name == 'R' || func_name == 'r') { // log-restricted likelihood
      dev1_l = LogRL_dev1(lambda_l, &params);
      dev1_h = LogRL_dev1(lambda_h, &params);
    } else {
      dev1_l = LogL_dev1(lambda_l, &params);
      dev1_h = LogL_dev1(lambda_h, &params);
    }

    if (dev1_l * dev1_h <= 0) {
      lambda_lh.push_back(make_pair(lambda_l, lambda_h));
    }
  }

  // If derivates do not change signs in any interval.
  if (lambda_lh.empty()) {
    if (func_name == 'R' || func_name == 'r') {
      logf_l = LogRL_f(l_min, &params);
      logf_h = LogRL_f(l_max, &params);
    } else {
      logf_l = LogL_f(l_min, &params);
      logf_h = LogL_f(l_max, &params);
    }

    if (logf_l >= logf_h) {
      lambda = l_min;
      logf = logf_l;
    } else {
      lambda = l_max;
      logf = logf_h;
    }
  } else {

    // If derivates change signs.
    int status;
    int iter = 0, max_iter = 100;
    double l, l_temp;

    gsl_function F;
    gsl_function_fdf FDF;

    F.params = &params;
    FDF.params = &params;

    if (func_name == 'R' || func_name == 'r') {
      F.function = &LogRL_dev1;
      FDF.f = &LogRL_dev1;
      FDF.df = &LogRL_dev2;
      FDF.fdf = &LogRL_dev12;
    } else {
      F.function = &LogL_dev1;
      FDF.f = &LogL_dev1;
      FDF.df = &LogL_dev2;
      FDF.fdf = &LogL_dev12;
    }

    const gsl_root_fsolver_type *T_f;
    gsl_root_fsolver *s_f;
    T_f = gsl_root_fsolver_brent;
    s_f = gsl_root_fsolver_alloc(T_f);

    const gsl_root_fdfsolver_type *T_fdf;
    gsl_root_fdfsolver *s_fdf;
    T_fdf = gsl_root_fdfsolver_newton;
    s_fdf = gsl_root_fdfsolver_alloc(T_fdf);

    for (vector<double>::size_type i = 0; i < lambda_lh.size(); ++i) {
      lambda_l = lambda_lh[i].first;
      lambda_h = lambda_lh[i].second;
      gsl_root_fsolver_set(s_f, &F, lambda_l, lambda_h);

      do {
        iter++;
        status = gsl_root_fsolver_iterate(s_f);
        l = gsl_root_fsolver_root(s_f);
        lambda_l = gsl_root_fsolver_x_lower(s_f);
        lambda_h = gsl_root_fsolver_x_upper(s_f);
        status = gsl_root_test_interval(lambda_l, lambda_h, 0, 1e-1);
      } while (status == GSL_CONTINUE && iter < max_iter);

      iter = 0;

      gsl_root_fdfsolver_set(s_fdf, &FDF, l);

      do {
        iter++;
        status = gsl_root_fdfsolver_iterate(s_fdf);
        l_temp = l;
        l = gsl_root_fdfsolver_root(s_fdf);
        status = gsl_root_test_delta(l, l_temp, 0, 1e-5);
      } while (status == GSL_CONTINUE && iter < max_iter && l > l_min &&
               l < l_max);

      l = l_temp;
      if (l < l_min) {
        l = l_min;
      }
      if (l > l_max) {
        l = l_max;
      }
      if (func_name == 'R' || func_name == 'r') {
        logf_l = LogRL_f(l, &params);
      } else {
        logf_l = LogL_f(l, &params);
      }

      if (i == 0) {
        logf = logf_l;
        lambda = l;
      } else if (logf < logf_l) {
        logf = logf_l;
        lambda = l;
      } else {
      }
    }
    gsl_root_fsolver_free(s_f);
    gsl_root_fdfsolver_free(s_fdf);

    if (func_name == 'R' || func_name == 'r') {
      logf_l = LogRL_f(l_min, &params);
      logf_h = LogRL_f(l_max, &params);
    } else {
      logf_l = LogL_f(l_min, &params);
      logf_h = LogL_f(l_max, &params);
    }

    if (logf_l > logf) {
      lambda = l_min;
      logf = logf_l;
    }
    if (logf_h > logf) {
      lambda = l_max;
      logf = logf_h;
    }
  }

  return;
}

// Calculate lambda in the null model.
void CalcLambda(const char func_name, const gsl_vector *eval,
                const gsl_matrix *UtW, const gsl_vector *Uty,
                const double l_min, const double l_max, const size_t n_region,
                double &lambda, double &logl_H0) {

  write(eval,"eval6");

  if (func_name != 'R' && func_name != 'L' && func_name != 'r' &&
      func_name != 'l') {
    cout << "func_name only takes 'R' or 'L': 'R' for "
         << "log-restricted likelihood, 'L' for log-likelihood." << endl;
    return;
  }

  size_t n_cvt = UtW->size2, ni_test = UtW->size1;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  // cout << "n_cvt " << n_cvt << ", ni_test " << ni_test << ", n_index " << n_index << endl;

  gsl_matrix *Uab = gsl_matrix_safe_alloc(ni_test, n_index);
  gsl_vector *ab = gsl_vector_safe_alloc(n_index);

  gsl_matrix_set_zero(Uab);
  write(UtW,"UtW");
  write(UtW,"Uty");
  CalcUab(UtW, Uty, Uab);
  write(Uab,"Uab");
  Calcab(UtW, Uty, ab);

  FUNC_PARAM param0 = {true, ni_test, n_cvt, eval, Uab, ab, 0};

  CalcLambda(func_name, param0, l_min, l_max, n_region, lambda, logl_H0);

  gsl_matrix_safe_free(Uab);
  gsl_vector_free(ab); // unused

  return;
}

// Obtain REMLE estimate for PVE using lambda_remle.
void CalcPve(const gsl_vector *eval, const gsl_matrix *UtW,
             const gsl_vector *Uty, const double lambda, const double trace_G,
             double &pve, double &pve_se) {
  size_t n_cvt = UtW->size2, ni_test = UtW->size1;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  gsl_matrix *Uab = gsl_matrix_safe_alloc(ni_test, n_index);
  gsl_vector *ab = gsl_vector_safe_alloc(n_index);

  gsl_matrix_set_zero(Uab);
  CalcUab(UtW, Uty, Uab);

  FUNC_PARAM param0 = {true, ni_test, n_cvt, eval, Uab, ab, 0};

  double se = safe_sqrt(-1.0 / LogRL_dev2(lambda, &param0));

  pve = trace_G * lambda / (trace_G * lambda + 1.0);
  pve_se = trace_G / ((trace_G * lambda + 1.0) * (trace_G * lambda + 1.0)) * se;

  gsl_matrix_safe_free(Uab);
  gsl_vector_free(ab); // unused
  return;
}

// Obtain REML estimate for Vg and Ve using lambda_remle.
// Obtain beta and se(beta) for coefficients.
// ab is not used when e_mode==0.
void CalcLmmVgVeBeta(const gsl_vector *eval, const gsl_matrix *UtW,
                     const gsl_vector *Uty, const double lambda, double &vg,
                     double &ve, gsl_vector *beta, gsl_vector *se_beta) {
  size_t n_cvt = UtW->size2, ni_test = UtW->size1;
  size_t n_index = (n_cvt + 2 + 1) * (n_cvt + 2) / 2;

  gsl_matrix *Uab = gsl_matrix_safe_alloc(ni_test, n_index);
  gsl_vector *ab = gsl_vector_safe_alloc(n_index);
  gsl_matrix *Pab = gsl_matrix_safe_alloc(n_cvt + 2, n_index);
  gsl_vector *Hi_eval = gsl_vector_safe_alloc(eval->size);
  gsl_vector *v_temp = gsl_vector_safe_alloc(eval->size);
  gsl_matrix *HiW = gsl_matrix_safe_alloc(eval->size, UtW->size2);
  gsl_matrix *WHiW = gsl_matrix_safe_alloc(UtW->size2, UtW->size2);
  gsl_vector *WHiy = gsl_vector_safe_alloc(UtW->size2);
  gsl_matrix *Vbeta = gsl_matrix_safe_alloc(UtW->size2, UtW->size2);

  gsl_matrix_set_zero(Uab);
  CalcUab(UtW, Uty, Uab);

  gsl_vector_safe_memcpy(v_temp, eval);
  gsl_vector_scale(v_temp, lambda);
  gsl_vector_set_all(Hi_eval, 1.0);
  gsl_vector_add_constant(v_temp, 1.0);
  gsl_vector_div(Hi_eval, v_temp);

  // Calculate beta.
  gsl_matrix_safe_memcpy(HiW, UtW);
  for (size_t i = 0; i < UtW->size2; i++) {
    gsl_vector_view HiW_col = gsl_matrix_column(HiW, i);
    gsl_vector_mul(&HiW_col.vector, Hi_eval);
  }
  fast_dgemm("T", "N", 1.0, HiW, UtW, 0.0, WHiW);
  gsl_blas_dgemv(CblasTrans, 1.0, HiW, Uty, 0.0, WHiy);

  int sig;
  gsl_permutation *pmt = gsl_permutation_alloc(UtW->size2);
  LUDecomp(WHiW, pmt, &sig);
  LUSolve(WHiW, pmt, WHiy, beta);
  LUInvert(WHiW, pmt, Vbeta);

  // Calculate vg and ve.
  CalcPab(n_cvt, 0, Hi_eval, Uab, ab, Pab);

  size_t index_yy = GetabIndex(n_cvt + 2, n_cvt + 2, n_cvt);
  double P_yy = gsl_matrix_safe_get(Pab, n_cvt, index_yy);

  ve = P_yy / (double)(ni_test - n_cvt);
  vg = ve * lambda;

  // With ve, calculate se(beta).
  gsl_matrix_scale(Vbeta, ve);

  // Obtain se_beta.
  for (size_t i = 0; i < Vbeta->size1; i++) {
    gsl_vector_set(se_beta, i, safe_sqrt(gsl_matrix_get(Vbeta, i, i)));
  }

  gsl_matrix_safe_free(Uab);
  gsl_matrix_free(Pab);
  gsl_vector_free(ab); // ab is unused
  gsl_vector_safe_free(Hi_eval);
  gsl_vector_safe_free(v_temp);
  gsl_matrix_safe_free(HiW);
  gsl_matrix_safe_free(WHiW);
  gsl_vector_safe_free(WHiy);
  gsl_matrix_safe_free(Vbeta);

  gsl_permutation_free(pmt);
  return;
}

void LMM::AnalyzeBimbamGXE(const gsl_matrix *U, const gsl_vector *eval,
                           const gsl_matrix *UtW, const gsl_vector *Uty,
                           const gsl_matrix *W, const gsl_vector *y,
                           const gsl_vector *env) {
  debug_msg("entering");
  auto infilen = file_gene.c_str();
  igzstream infile(infilen, igzstream::in);
  if (!infile) {
    cout << "error reading genotype file:" << file_geno << endl;
    return;
  }

  clock_t time_start = clock();

  string line;
  char *ch_ptr;

  double lambda_mle = 0, lambda_remle = 0, beta = 0, se = 0, p_wald = 0;
  double p_lrt = 0, p_score = 0;
  double logl_H1 = 0.0, logl_H0 = 0.0;
  int n_miss, c_phen;
  double geno, x_mean;

  // Calculate basic quantities.
  size_t n_index = (n_cvt + 2 + 2 + 1) * (n_cvt + 2 + 2) / 2;

  gsl_vector *x = gsl_vector_safe_alloc(U->size1);
  gsl_vector *x_miss = gsl_vector_safe_alloc(U->size1);
  gsl_vector *Utx = gsl_vector_safe_alloc(U->size2);
  gsl_matrix *Uab = gsl_matrix_safe_alloc(U->size2, n_index);
  gsl_vector *ab = gsl_vector_safe_alloc(n_index);

  gsl_matrix *UtW_expand = gsl_matrix_safe_alloc(U->size1, UtW->size2 + 2);
  gsl_matrix_view UtW_expand_mat =
      gsl_matrix_submatrix(UtW_expand, 0, 0, U->size1, UtW->size2);
  gsl_matrix_safe_memcpy(&UtW_expand_mat.matrix, UtW);
  gsl_vector_view UtW_expand_env = gsl_matrix_column(UtW_expand, UtW->size2);
  gsl_blas_dgemv(CblasTrans, 1.0, U, env, 0.0, &UtW_expand_env.vector);
  gsl_vector_view UtW_expand_x = gsl_matrix_column(UtW_expand, UtW->size2 + 1);

  // Start reading genotypes and analyze.
  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    safeGetline(infile, line).eof();
    if (t % d_pace == 0 || t == (ns_total - 1)) {
      ProgressBar("Reading SNPs", t, ns_total - 1);
    }
    if (indicator_snp[t] == 0) {
      continue;
    }

    ch_ptr = strtok_safe2((char *)line.c_str(), " , \t",infilen);
    ch_ptr = strtok_safe2(NULL, " , \t",infilen);
    ch_ptr = strtok_safe2(NULL, " , \t",infilen);

    x_mean = 0.0;
    c_phen = 0;
    n_miss = 0;
    gsl_vector_set_zero(x_miss);
    for (size_t i = 0; i < ni_total; ++i) {
      ch_ptr = strtok_safe2(NULL, " , \t",infilen);
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
    gsl_blas_dgemv(CblasTrans, 1.0, U, x, 0.0, &UtW_expand_x.vector);
    gsl_vector_mul(x, env);
    gsl_blas_dgemv(CblasTrans, 1.0, U, x, 0.0, Utx);
    time_UtX += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    gsl_matrix_set_zero(Uab);
    CalcUab(UtW_expand, Uty, Uab);

    if (a_mode == 2 || a_mode == 4) {
      FUNC_PARAM param0 = {true, ni_test, n_cvt + 2, eval, Uab, ab, 0};
      CalcLambda('L', param0, l_min, l_max, n_region, lambda_mle, logl_H0);
    }

    CalcUab(UtW_expand, Uty, Utx, Uab);

    time_start = clock();
    FUNC_PARAM param1 = {false, ni_test, n_cvt + 2, eval, Uab, ab, 0};

    // 3 is before 1.
    if (a_mode == 3 || a_mode == 4) {
      CalcRLScore(l_mle_null, param1, beta, se, p_score);
    }

    if (a_mode == 1 || a_mode == 4) {
      CalcLambda('R', param1, l_min, l_max, n_region, lambda_remle, logl_H1);
      CalcRLWald(lambda_remle, param1, beta, se, p_wald);
    }

    if (a_mode == 2 || a_mode == 4) {
      CalcLambda('L', param1, l_min, l_max, n_region, lambda_mle, logl_H1);
      p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), 1);
    }

    if (x_mean > 1) {
      beta *= -1;
    }

    time_opt += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    // Store summary data.
    SUMSTAT SNPs = {beta, se, lambda_remle, lambda_mle, p_wald, p_lrt, p_score, logl_H1};
    sumStat.push_back(SNPs);
  }
  cout << endl;

  gsl_vector_safe_free(x);
  gsl_vector_safe_free(x_miss);
  gsl_vector_safe_free(Utx);
  gsl_matrix_safe_free(Uab);
  gsl_vector_free(ab); // unused

  gsl_matrix_safe_free(UtW_expand);

  infile.close();
  infile.clear();

  return;
}

void LMM::AnalyzePlinkGXE(const gsl_matrix *U, const gsl_vector *eval,
                          const gsl_matrix *UtW, const gsl_vector *Uty,
                          const gsl_matrix *W, const gsl_vector *y,
                          const gsl_vector *env) {
  string file_bed = file_bfile + ".bed";
  debug_msg(file_bed);
  ifstream infile(file_bed.c_str(), ios::binary);
  if (!infile) {
    cout << "error reading bed file:" << file_bed << endl;
    return;
  }

  clock_t time_start = clock();

  char ch[1];
  bitset<8> b;

  double lambda_mle = 0, lambda_remle = 0, beta = 0, se = 0, p_wald = 0;
  double p_lrt = 0, p_score = 0;
  double logl_H1 = 0.0, logl_H0 = 0.0;
  int n_bit, n_miss, ci_total, ci_test;
  double geno, x_mean;

  // Calculate basic quantities.
  size_t n_index = (n_cvt + 2 + 2 + 1) * (n_cvt + 2 + 2) / 2;

  gsl_vector *x = gsl_vector_safe_alloc(U->size1);
  gsl_vector *Utx = gsl_vector_safe_alloc(U->size2);
  gsl_matrix *Uab = gsl_matrix_safe_alloc(U->size2, n_index);
  gsl_vector *ab = gsl_vector_safe_alloc(n_index);

  gsl_matrix *UtW_expand = gsl_matrix_safe_alloc(U->size1, UtW->size2 + 2);
  gsl_matrix_view UtW_expand_mat =
      gsl_matrix_submatrix(UtW_expand, 0, 0, U->size1, UtW->size2);
  gsl_matrix_safe_memcpy(&UtW_expand_mat.matrix, UtW);
  gsl_vector_view UtW_expand_env = gsl_matrix_column(UtW_expand, UtW->size2);
  gsl_blas_dgemv(CblasTrans, 1.0, U, env, 0.0, &UtW_expand_env.vector);
  gsl_vector_view UtW_expand_x = gsl_matrix_column(UtW_expand, UtW->size2 + 1);

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

    // n_bit, and 3 is the number of magic numbers
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
    gsl_blas_dgemv(CblasTrans, 1.0, U, x, 0.0, &UtW_expand_x.vector);
    gsl_vector_mul(x, env);
    gsl_blas_dgemv(CblasTrans, 1.0, U, x, 0.0, Utx);
    time_UtX += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    gsl_matrix_set_zero(Uab);
    CalcUab(UtW_expand, Uty, Uab);

    if (a_mode == 2 || a_mode == 4) {
      FUNC_PARAM param0 = {true, ni_test, n_cvt + 2, eval, Uab, ab, 0};
      CalcLambda('L', param0, l_min, l_max, n_region, lambda_mle, logl_H0);
    }

    CalcUab(UtW_expand, Uty, Utx, Uab);

    time_start = clock();
    FUNC_PARAM param1 = {false, ni_test, n_cvt + 2, eval, Uab, ab, 0};

    // 3 is before 1, for beta.
    if (a_mode == 3 || a_mode == 4) {
      CalcRLScore(l_mle_null, param1, beta, se, p_score);
    }

    if (a_mode == 1 || a_mode == 4) {
      CalcLambda('R', param1, l_min, l_max, n_region, lambda_remle, logl_H1);
      CalcRLWald(lambda_remle, param1, beta, se, p_wald);
    }

    if (a_mode == 2 || a_mode == 4) {
      CalcLambda('L', param1, l_min, l_max, n_region, lambda_mle, logl_H1);
      p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), 1);
    }

    if (x_mean > 1) {
      beta *= -1;
    }

    time_opt += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    // Store summary data.
    SUMSTAT SNPs = {beta, se, lambda_remle, lambda_mle, p_wald, p_lrt, p_score, logl_H1};
    sumStat.push_back(SNPs);
  }
  cout << endl;

  gsl_vector_safe_free(x);
  gsl_vector_safe_free(Utx);
  gsl_matrix_safe_free(Uab);
  gsl_vector_free(ab); // unused

  gsl_matrix_safe_free(UtW_expand);

  infile.close();
  infile.clear();

  return;
}
