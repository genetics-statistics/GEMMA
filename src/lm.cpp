/*
 Genome-wide Efficient Mixed Model Association (GEMMA)
 Copyright (C) 2011-2017 Xiang Zhou

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"

#include "gsl/gsl_cdf.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_min.h"
#include "gsl/gsl_roots.h"

// #include "eigenlib.h"
#include "gzstream.h"
#include "lapack.h"
#include "lm.h"

using namespace std;

void LM::CopyFromParam(PARAM &cPar) {
  a_mode = cPar.a_mode;
  d_pace = cPar.d_pace;

  file_bfile = cPar.file_bfile;
  file_geno = cPar.file_geno;
  file_out = cPar.file_out;
  path_out = cPar.path_out;
  file_gene = cPar.file_gene;

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

  return;
}

void LM::CopyToParam(PARAM &cPar) {
  cPar.time_opt = time_opt;
  cPar.ng_test = ng_test;
  return;
}

void LM::WriteFiles() {
  string file_str;
  file_str = path_out + "/" + file_out;
  file_str += ".assoc.txt";

  ofstream outfile(file_str.c_str(), ofstream::out);
  if (!outfile) {
    cout << "error writing file: " << file_str.c_str() << endl;
    return;
  }

  if (!file_gene.empty()) {
    outfile << "geneID"
            << "\t";

    if (a_mode == 51) {
      outfile << "beta"
              << "\t"
              << "se"
              << "\t"
              << "p_wald" << endl;
    } else if (a_mode == 52) {
      outfile << "p_lrt" << endl;
    } else if (a_mode == 53) {
      outfile << "beta"
              << "\t"
              << "se"
              << "\t"
              << "p_score" << endl;
    } else if (a_mode == 54) {
      outfile << "beta"
              << "\t"
              << "se"
              << "\t"
              << "p_wald"
              << "\t"
              << "p_lrt"
              << "\t"
              << "p_score" << endl;
    } else {
    }

    for (vector<SUMSTAT>::size_type t = 0; t < sumStat.size(); ++t) {
      outfile << snpInfo[t].rs_number << "\t";

      if (a_mode == 51) {
        outfile << scientific << setprecision(6) << sumStat[t].beta << "\t"
                << sumStat[t].se << "\t" << sumStat[t].p_wald << endl;
      } else if (a_mode == 52) {
        outfile << scientific << setprecision(6) << "\t" << sumStat[t].p_lrt
                << endl;
      } else if (a_mode == 53) {
        outfile << scientific << setprecision(6) << sumStat[t].beta << "\t"
                << sumStat[t].se << "\t" << sumStat[t].p_score << endl;
      } else if (a_mode == 54) {
        outfile << scientific << setprecision(6) << sumStat[t].beta << "\t"
                << sumStat[t].se << "\t" << sumStat[t].p_wald << "\t"
                << sumStat[t].p_lrt << "\t" << sumStat[t].p_score << endl;
      } else {
      }
    }
  } else {
    outfile << "chr"
            << "\t"
            << "rs"
            << "\t"
            << "ps"
            << "\t"
            << "n_mis"
            << "\t"
            << "n_obs"
            << "\t"
            << "allele1"
            << "\t"
            << "allele0"
            << "\t"
            << "af"
            << "\t";

    if (a_mode == 51) {
      outfile << "beta"
              << "\t"
              << "se"
              << "\t"
              << "p_wald" << endl;
    } else if (a_mode == 52) {
      outfile << "p_lrt" << endl;
    } else if (a_mode == 53) {
      outfile << "beta"
              << "\t"
              << "se"
              << "\t"
              << "p_score" << endl;
    } else if (a_mode == 54) {
      outfile << "beta"
              << "\t"
              << "se"
              << "\t"
              << "p_wald"
              << "\t"
              << "p_lrt"
              << "\t"
              << "p_score" << endl;
    } else {
    }

    size_t t = 0;
    for (size_t i = 0; i < snpInfo.size(); ++i) {
      if (indicator_snp[i] == 0) {
        continue;
      }

      outfile << snpInfo[i].chr << "\t" << snpInfo[i].rs_number << "\t"
              << snpInfo[i].base_position << "\t" << snpInfo[i].n_miss << "\t"
              << ni_test - snpInfo[i].n_miss << "\t" << snpInfo[i].a_minor
              << "\t" << snpInfo[i].a_major << "\t" << fixed << setprecision(3)
              << snpInfo[i].maf << "\t";

      if (a_mode == 51) {
        outfile << scientific << setprecision(6) << sumStat[t].beta << "\t"
                << sumStat[t].se << "\t" << sumStat[t].p_wald << endl;
      } else if (a_mode == 52) {
        outfile << scientific << setprecision(6) << sumStat[t].p_lrt << endl;
      } else if (a_mode == 53) {
        outfile << scientific << setprecision(6) << sumStat[t].beta << "\t"
                << sumStat[t].se << "\t" << sumStat[t].p_score << endl;
      } else if (a_mode == 54) {
        outfile << scientific << setprecision(6) << sumStat[t].beta << "\t"
                << sumStat[t].se << "\t" << sumStat[t].p_wald << "\t"
                << sumStat[t].p_lrt << "\t" << sumStat[t].p_score << endl;
      } else {
      }
      t++;
    }
  }

  outfile.close();
  outfile.clear();
  return;
}

void CalcvPv(const gsl_matrix *WtWi, const gsl_vector *Wty,
             const gsl_vector *Wtx, const gsl_vector *y, const gsl_vector *x,
             double &xPwy, double &xPwx) {
  size_t c_size = Wty->size;
  double d;

  gsl_vector *WtWiWtx = gsl_vector_alloc(c_size);

  gsl_blas_ddot(x, x, &xPwx);
  gsl_blas_ddot(x, y, &xPwy);
  gsl_blas_dgemv(CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);

  gsl_blas_ddot(WtWiWtx, Wtx, &d);
  xPwx -= d;

  gsl_blas_ddot(WtWiWtx, Wty, &d);
  xPwy -= d;

  gsl_vector_free(WtWiWtx);

  return;
}

void CalcvPv(const gsl_matrix *WtWi, const gsl_vector *Wty, const gsl_vector *y,
             double &yPwy) {
  size_t c_size = Wty->size;
  double d;

  gsl_vector *WtWiWty = gsl_vector_alloc(c_size);

  gsl_blas_ddot(y, y, &yPwy);
  gsl_blas_dgemv(CblasNoTrans, 1.0, WtWi, Wty, 0.0, WtWiWty);

  gsl_blas_ddot(WtWiWty, Wty, &d);
  yPwy -= d;

  gsl_vector_free(WtWiWty);

  return;
}

// Calculate p-values and beta/se in a linear model.
void LmCalcP(const size_t test_mode, const double yPwy, const double xPwy,
             const double xPwx, const double df, const size_t n_size,
             double &beta, double &se, double &p_wald, double &p_lrt,
             double &p_score) {
  double yPxy = yPwy - xPwy * xPwy / xPwx;
  double se_wald, se_score;

  beta = xPwy / xPwx;
  se_wald = sqrt(yPxy / (df * xPwx));
  se_score = sqrt(yPwy / ((double)n_size * xPwx));

  p_wald = gsl_cdf_fdist_Q(beta * beta / (se_wald * se_wald), 1.0, df);
  p_score = gsl_cdf_fdist_Q(beta * beta / (se_score * se_score), 1.0, df);
  p_lrt = gsl_cdf_chisq_Q((double)n_size * (log(yPwy) - log(yPxy)), 1);

  if (test_mode == 3) {
    se = se_score;
  } else {
    se = se_wald;
  }

  return;
}

void LM::AnalyzeGene(const gsl_matrix *W, const gsl_vector *x) {
  debug_msg("entering");
  ifstream infile(file_gene.c_str(), ifstream::in);
  if (!infile) {
    cout << "error reading gene expression file:" << file_gene << endl;
    return;
  }

  clock_t time_start = clock();

  string line;
  char *ch_ptr;

  double beta = 0, se = 0, p_wald = 0, p_lrt = 0, p_score = 0;
  int c_phen;
  string rs; // Gene id.
  double d;

  // Calculate some basic quantities.
  double yPwy, xPwy, xPwx;
  double df = (double)W->size1 - (double)W->size2 - 1.0;

  gsl_vector *y = gsl_vector_alloc(W->size1);

  gsl_matrix *WtW = gsl_matrix_alloc(W->size2, W->size2);
  gsl_matrix *WtWi = gsl_matrix_alloc(W->size2, W->size2);
  gsl_vector *Wty = gsl_vector_alloc(W->size2);
  gsl_vector *Wtx = gsl_vector_alloc(W->size2);
  gsl_permutation *pmt = gsl_permutation_alloc(W->size2);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
  int sig;
  LUDecomp(WtW, pmt, &sig);
  LUInvert(WtW, pmt, WtWi);

  gsl_blas_dgemv(CblasTrans, 1.0, W, x, 0.0, Wtx);
  CalcvPv(WtWi, Wtx, x, xPwx);

  // Header.
  getline(infile, line);

  for (size_t t = 0; t < ng_total; t++) {
    getline(infile, line);
    if (t % d_pace == 0 || t == ng_total - 1) {
      ProgressBar("Performing Analysis", t, ng_total - 1);
    }
    ch_ptr = strtok_safe((char *)line.c_str(), " , \t");
    rs = ch_ptr;

    c_phen = 0;
    for (size_t i = 0; i < indicator_idv.size(); ++i) {
      ch_ptr = strtok_safe(NULL, " , \t");
      if (indicator_idv[i] == 0) {
        continue;
      }

      d = atof(ch_ptr);
      gsl_vector_set(y, c_phen, d);

      c_phen++;
    }

    // Calculate statistics.
    time_start = clock();

    gsl_blas_dgemv(CblasTrans, 1.0, W, y, 0.0, Wty);
    CalcvPv(WtWi, Wtx, Wty, x, y, xPwy, yPwy);
    LmCalcP(a_mode - 50, yPwy, xPwy, xPwx, df, W->size1, beta, se, p_wald,
            p_lrt, p_score);

    time_opt += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    // Store summary data.
    SUMSTAT SNPs = {beta, se, 0.0, 0.0, p_wald, p_lrt, p_score, -0.0 };
    sumStat.push_back(SNPs);
  }
  cout << endl;

  gsl_vector_free(y);

  gsl_matrix_free(WtW);
  gsl_matrix_free(WtWi);
  gsl_vector_free(Wty);
  gsl_vector_free(Wtx);
  gsl_permutation_free(pmt);

  infile.close();
  infile.clear();

  return;
}

void LM::AnalyzeBimbam(const gsl_matrix *W, const gsl_vector *y) {
  debug_msg("entering");
  igzstream infile(file_geno.c_str(), igzstream::in);
  if (!infile) {
    cout << "error reading genotype file:" << file_geno << endl;
    return;
  }

  clock_t time_start = clock();

  string line;
  char *ch_ptr;

  double beta = 0, se = 0, p_wald = 0, p_lrt = 0, p_score = 0;
  int n_miss, c_phen;
  double geno, x_mean;

  // Calculate some basic quantities.
  double yPwy, xPwy, xPwx;
  double df = (double)W->size1 - (double)W->size2 - 1.0;

  gsl_vector *x = gsl_vector_alloc(W->size1);
  gsl_vector *x_miss = gsl_vector_alloc(W->size1);

  gsl_matrix *WtW = gsl_matrix_alloc(W->size2, W->size2);
  gsl_matrix *WtWi = gsl_matrix_alloc(W->size2, W->size2);
  gsl_vector *Wty = gsl_vector_alloc(W->size2);
  gsl_vector *Wtx = gsl_vector_alloc(W->size2);
  gsl_permutation *pmt = gsl_permutation_alloc(W->size2);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
  int sig;
  LUDecomp(WtW, pmt, &sig);
  LUInvert(WtW, pmt, WtWi);

  gsl_blas_dgemv(CblasTrans, 1.0, W, y, 0.0, Wty);
  CalcvPv(WtWi, Wty, y, yPwy);

  // Start reading genotypes and analyze.
  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    getline(infile, line);
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

    // Calculate statistics.
    time_start = clock();

    gsl_blas_dgemv(CblasTrans, 1.0, W, x, 0.0, Wtx);
    CalcvPv(WtWi, Wty, Wtx, y, x, xPwy, xPwx);
    LmCalcP(a_mode - 50, yPwy, xPwy, xPwx, df, W->size1, beta, se, p_wald,
            p_lrt, p_score);

    time_opt += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    // Store summary data.
    SUMSTAT SNPs = {beta, se, 0.0, 0.0, p_wald, p_lrt, p_score, -0.0};
    sumStat.push_back(SNPs);
  }
  cout << endl;

  gsl_vector_free(x);
  gsl_vector_free(x_miss);

  gsl_matrix_free(WtW);
  gsl_matrix_free(WtWi);
  gsl_vector_free(Wty);
  gsl_vector_free(Wtx);
  gsl_permutation_free(pmt);

  infile.close();
  infile.clear();

  return;
}

void LM::AnalyzePlink(const gsl_matrix *W, const gsl_vector *y) {
  debug_msg("entering");
  string file_bed = file_bfile + ".bed";
  ifstream infile(file_bed.c_str(), ios::binary);
  if (!infile) {
    cout << "error reading bed file:" << file_bed << endl;
    return;
  }

  clock_t time_start = clock();

  char ch[1];
  bitset<8> b;

  double beta = 0, se = 0, p_wald = 0, p_lrt = 0, p_score = 0;
  int n_bit, n_miss, ci_total, ci_test;
  double geno, x_mean;

  // Calculate some basic quantities.
  double yPwy, xPwy, xPwx;
  double df = (double)W->size1 - (double)W->size2 - 1.0;

  gsl_vector *x = gsl_vector_alloc(W->size1);

  gsl_matrix *WtW = gsl_matrix_alloc(W->size2, W->size2);
  gsl_matrix *WtWi = gsl_matrix_alloc(W->size2, W->size2);
  gsl_vector *Wty = gsl_vector_alloc(W->size2);
  gsl_vector *Wtx = gsl_vector_alloc(W->size2);
  gsl_permutation *pmt = gsl_permutation_alloc(W->size2);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
  int sig;
  LUDecomp(WtW, pmt, &sig);
  LUInvert(WtW, pmt, WtWi);

  gsl_blas_dgemv(CblasTrans, 1.0, W, y, 0.0, Wty);
  CalcvPv(WtWi, Wty, y, yPwy);

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

    // Calculate statistics.
    time_start = clock();

    gsl_blas_dgemv(CblasTrans, 1.0, W, x, 0.0, Wtx);
    CalcvPv(WtWi, Wty, Wtx, y, x, xPwy, xPwx);
    LmCalcP(a_mode - 50, yPwy, xPwy, xPwx, df, W->size1, beta, se, p_wald,
            p_lrt, p_score);

    // store summary data
    SUMSTAT SNPs = {beta, se, 0.0, 0.0, p_wald, p_lrt, p_score, -0.0};
    sumStat.push_back(SNPs);

    time_opt += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
  }
  cout << endl;

  gsl_vector_free(x);

  gsl_matrix_free(WtW);
  gsl_matrix_free(WtWi);
  gsl_vector_free(Wty);
  gsl_vector_free(Wtx);
  gsl_permutation_free(pmt);

  infile.close();
  infile.clear();

  return;
}

// Make sure that both y and X are centered already.
void MatrixCalcLmLR(const gsl_matrix *X, const gsl_vector *y,
                    vector<pair<size_t, double>> &pos_loglr) {
  double yty, xty, xtx, log_lr;
  gsl_blas_ddot(y, y, &yty);

  for (size_t i = 0; i < X->size2; ++i) {
    gsl_vector_const_view X_col = gsl_matrix_const_column(X, i);
    gsl_blas_ddot(&X_col.vector, &X_col.vector, &xtx);
    gsl_blas_ddot(&X_col.vector, y, &xty);

    log_lr = 0.5 * (double)y->size * (log(yty) - log(yty - xty * xty / xtx));
    pos_loglr.push_back(make_pair(i, log_lr));
  }

  return;
}
