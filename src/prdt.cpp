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

#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include <bitset>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "gzstream.h"
#include "gemma_io.h"
#include "lapack.h"
#include "mathfunc.h"
#include "prdt.h"

using namespace std;

void PRDT::CopyFromParam(PARAM &cPar) {
  a_mode = cPar.a_mode;
  d_pace = cPar.d_pace;

  file_bfile = cPar.file_bfile;
  file_geno = cPar.file_geno;
  file_out = cPar.file_out;
  path_out = cPar.path_out;

  indicator_pheno = cPar.indicator_pheno;
  indicator_cvt = cPar.indicator_cvt;
  indicator_idv = cPar.indicator_idv;

  snpInfo = cPar.snpInfo;
  mapRS2est = cPar.mapRS2est;

  time_eigen = 0;

  n_ph = cPar.n_ph;
  np_obs = cPar.np_obs;
  np_miss = cPar.np_miss;
  ns_total = cPar.ns_total;
  ns_test = 0;

  return;
}

void PRDT::CopyToParam(PARAM &cPar) {
  cPar.ns_test = ns_test;
  cPar.time_eigen = time_eigen;

  return;
}

void PRDT::WriteFiles(gsl_vector *y_prdt) {
  string file_str;
  file_str = path_out + "/" + file_out;
  file_str += ".";
  file_str += "prdt";
  file_str += ".txt";

  ofstream outfile(file_str.c_str(), ofstream::out);
  if (!outfile) {
    cout << "error writing file: " << file_str.c_str() << endl;
    return;
  }

  size_t ci_test = 0;
  for (size_t i = 0; i < indicator_idv.size(); i++) {
    if (indicator_idv[i] == 1) {
      outfile << "NA" << endl;
    } else {
      outfile << gsl_vector_get(y_prdt, ci_test) << endl;
      ci_test++;
    }
  }

  outfile.close();
  outfile.clear();
  return;
}

void PRDT::WriteFiles(gsl_matrix *Y_full) {
  string file_str;
  file_str = path_out + "/" + file_out;
  file_str += ".prdt.txt";

  ofstream outfile(file_str.c_str(), ofstream::out);
  if (!outfile) {
    cout << "error writing file: " << file_str.c_str() << endl;
    return;
  }

  size_t ci_test = 0;
  for (size_t i = 0; i < indicator_cvt.size(); i++) {
    if (indicator_cvt[i] == 0) {
      outfile << "NA" << endl;
    } else {
      for (size_t j = 0; j < Y_full->size2; j++) {
        outfile << gsl_matrix_get(Y_full, ci_test, j) << "\t";
      }
      outfile << endl;
      ci_test++;
    }
  }

  outfile.close();
  outfile.clear();
  return;
}

void PRDT::AddBV(gsl_matrix *G, const gsl_vector *u_hat, gsl_vector *y_prdt) {
  size_t ni_test = u_hat->size, ni_total = G->size1;

  gsl_matrix *Goo = gsl_matrix_alloc(ni_test, ni_test);
  gsl_matrix *Gfo = gsl_matrix_alloc(ni_total - ni_test, ni_test);
  gsl_matrix *U = gsl_matrix_alloc(ni_test, ni_test);
  gsl_vector *eval = gsl_vector_alloc(ni_test);
  gsl_vector *Utu = gsl_vector_alloc(ni_test);
  gsl_vector *w = gsl_vector_alloc(ni_total);
  gsl_permutation *pmt = gsl_permutation_alloc(ni_test);

  // center matrix G based on indicator_idv
  for (size_t i = 0; i < ni_total; i++) {
    gsl_vector_set(w, i, indicator_idv[i]);
  }
  CenterMatrix(G, w);

  // obtain Koo and Kfo
  size_t o_i = 0, o_j = 0;
  double d;
  for (size_t i = 0; i < indicator_idv.size(); i++) {
    o_j = 0;
    for (size_t j = 0; j < indicator_idv.size(); j++) {
      d = gsl_matrix_get(G, i, j);
      if (indicator_idv[i] == 1 && indicator_idv[j] == 1) {
        gsl_matrix_set(Goo, o_i, o_j, d);
      }
      if (indicator_idv[i] == 0 && indicator_idv[j] == 1) {
        gsl_matrix_set(Gfo, i - o_i, o_j, d);
      }
      if (indicator_idv[j] == 1) {
        o_j++;
      }
    }
    if (indicator_idv[i] == 1) {
      o_i++;
    }
  }

  // matrix operations to get u_prdt
  cout << "Start Eigen-Decomposition..." << endl;
  clock_t time_start = clock();
  EigenDecomp(Goo, U, eval, 0);
  for (size_t i = 0; i < eval->size; i++) {
    if (gsl_vector_get(eval, i) < 1e-10) {
      gsl_vector_set(eval, i, 0);
    }
  }

  time_eigen = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

  gsl_blas_dgemv(CblasTrans, 1.0, U, u_hat, 0.0, Utu);
  for (size_t i = 0; i < eval->size; i++) {
    d = gsl_vector_get(eval, i);
    if (d != 0) {
      d = gsl_vector_get(Utu, i) / d;
      gsl_vector_set(Utu, i, d);
    }
  }
  gsl_blas_dgemv(CblasNoTrans, 1.0, U, Utu, 0.0, eval);
  gsl_blas_dgemv(CblasNoTrans, 1.0, Gfo, eval, 1.0, y_prdt);

  // Free matrices.
  gsl_matrix_free(Goo);
  gsl_matrix_free(Gfo);
  gsl_matrix_free(U);
  gsl_vector_free(eval);
  gsl_vector_free(Utu);
  gsl_vector_free(w);
  gsl_permutation_free(pmt);

  return;
}

void PRDT::AnalyzeBimbam(gsl_vector *y_prdt) {
  debug_msg("entering");
  igzstream infile(file_geno.c_str(), igzstream::in);
  if (!infile) {
    cout << "error reading genotype file:" << file_geno << endl;
    return;
  }

  string line;
  char *ch_ptr;
  string rs;

  size_t n_miss, n_train_nomiss, c_phen;
  double geno, x_mean, x_train_mean, effect_size;

  gsl_vector *x = gsl_vector_alloc(y_prdt->size);
  gsl_vector *x_miss = gsl_vector_alloc(y_prdt->size);

  ns_test = 0;

  // Start reading genotypes and analyze.
  for (size_t t = 0; t < ns_total; ++t) {
    safeGetline(infile, line).eof();
    if (t % d_pace == 0 || t == (ns_total - 1)) {
      ProgressBar("Reading SNPs  ", t, ns_total - 1);
    }

    ch_ptr = strtok((char *)line.c_str(), " , \t");
    rs = ch_ptr;
    ch_ptr = strtok(NULL, " , \t");
    ch_ptr = strtok(NULL, " , \t");

    if (mapRS2est.count(rs) == 0) {
      continue;
    } else {
      effect_size = mapRS2est[rs];
    }

    x_mean = 0.0;
    c_phen = 0;
    n_miss = 0;
    x_train_mean = 0;
    n_train_nomiss = 0;

    gsl_vector_set_zero(x_miss);

    for (size_t i = 0; i < indicator_idv.size(); ++i) {
      ch_ptr = strtok(NULL, " , \t");
      if (indicator_idv[i] == 1) {
        if (strcmp(ch_ptr, "NA") != 0) {
          geno = atof(ch_ptr);
          x_train_mean += geno;
          n_train_nomiss++;
        }
      } else {
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
    }

    if (x->size == n_miss) {
      cout << "snp " << rs << " has missing genotype for all "
           << "individuals and will be ignored." << endl;
      continue;
    }

    x_mean /= (double)(x->size - n_miss);
    x_train_mean /= (double)(n_train_nomiss);

    for (size_t i = 0; i < x->size; ++i) {
      geno = gsl_vector_get(x, i);
      if (gsl_vector_get(x_miss, i) == 0) {
        gsl_vector_set(x, i, x_mean - x_train_mean);
      } else {
        gsl_vector_set(x, i, geno - x_train_mean);
      }
    }

    gsl_vector_scale(x, effect_size);
    gsl_vector_add(y_prdt, x);

    ns_test++;
  }
  cout << endl;

  gsl_vector_free(x);
  gsl_vector_free(x_miss);

  infile.close();
  infile.clear();

  return;
}

void PRDT::AnalyzePlink(gsl_vector *y_prdt) {
  debug_msg("entering");
  string file_bed = file_bfile + ".bed";
  ifstream infile(file_bed.c_str(), ios::binary);
  if (!infile) {
    cout << "error reading bed file:" << file_bed << endl;
    return;
  }

  char ch[1];
  bitset<8> b;
  string rs;

  size_t n_bit, n_miss, ci_total, ci_test, n_train_nomiss;
  double geno, x_mean, x_train_mean, effect_size;

  gsl_vector *x = gsl_vector_alloc(y_prdt->size);

  // Calculate n_bit and c, the number of bit for each SNP.
  if (indicator_idv.size() % 4 == 0) {
    n_bit = indicator_idv.size() / 4;
  } else {
    n_bit = indicator_idv.size() / 4 + 1;
  }

  // Print the first 3 magic numbers.
  for (size_t i = 0; i < 3; ++i) {
    infile.read(ch, 1);
    b = ch[0];
  }

  ns_test = 0;

  for (vector<SNPINFO>::size_type t = 0; t < snpInfo.size(); ++t) {
    if (t % d_pace == 0 || t == snpInfo.size() - 1) {
      ProgressBar("Reading SNPs  ", t, snpInfo.size() - 1);
    }

    rs = snpInfo[t].rs_number;

    if (mapRS2est.count(rs) == 0) {
      continue;
    } else {
      effect_size = mapRS2est[rs];
    }

    // n_bit, and 3 is the number of magic numbers.
    infile.seekg(t * n_bit + 3);

    // Read genotypes.
    x_mean = 0.0;
    n_miss = 0;
    ci_total = 0;
    ci_test = 0;
    x_train_mean = 0;
    n_train_nomiss = 0;
    for (size_t i = 0; i < n_bit; ++i) {
      infile.read(ch, 1);
      b = ch[0];

      // Minor allele homozygous: 2.0; major: 0.0.
      for (size_t j = 0; j < 4; ++j) {
        if ((i == (n_bit - 1)) && ci_total == indicator_idv.size()) {
          break;
        }
        if (indicator_idv[ci_total] == 1) {
          if (b[2 * j] == 0) {
            if (b[2 * j + 1] == 0) {
              x_train_mean += 2.0;
              n_train_nomiss++;
            } else {
              x_train_mean += 1.0;
              n_train_nomiss++;
            }
          } else {
            if (b[2 * j + 1] == 1) {
              n_train_nomiss++;
            } else {
            }
          }
        } else {
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
          ci_test++;
        }
        ci_total++;
      }
    }

    if (x->size == n_miss) {
      cout << "snp " << rs << " has missing genotype for all "
           << "individuals and will be ignored." << endl;
      continue;
    }

    x_mean /= (double)(x->size - n_miss);
    x_train_mean /= (double)(n_train_nomiss);

    for (size_t i = 0; i < x->size; ++i) {
      geno = gsl_vector_get(x, i);
      if (geno == -9) {
        gsl_vector_set(x, i, x_mean - x_train_mean);
      } else {
        gsl_vector_set(x, i, geno - x_train_mean);
      }
    }

    gsl_vector_scale(x, effect_size);
    gsl_vector_add(y_prdt, x);

    ns_test++;
  }
  cout << endl;

  gsl_vector_free(x);

  infile.close();
  infile.clear();

  return;
}

// Predict missing phenotypes using ridge regression.
// Y_hat contains fixed effects
void PRDT::MvnormPrdt(const gsl_matrix *Y_hat, const gsl_matrix *H,
                      gsl_matrix *Y_full) {
  gsl_vector *y_obs = gsl_vector_alloc(np_obs);
  gsl_vector *y_miss = gsl_vector_alloc(np_miss);
  gsl_matrix *H_oo = gsl_matrix_alloc(np_obs, np_obs);
  gsl_matrix *H_mo = gsl_matrix_alloc(np_miss, np_obs);
  gsl_vector *Hiy = gsl_vector_alloc(np_obs);

  size_t c_obs1 = 0, c_obs2 = 0, c_miss1 = 0, c_miss2 = 0;

  // Obtain H_oo, H_mo.
  c_obs1 = 0;
  c_miss1 = 0;
  for (vector<int>::size_type i1 = 0; i1 < indicator_pheno.size(); ++i1) {
    if (indicator_cvt[i1] == 0) {
      continue;
    }
    for (vector<int>::size_type j1 = 0; j1 < n_ph; ++j1) {

      c_obs2 = 0;
      c_miss2 = 0;
      for (vector<int>::size_type i2 = 0; i2 < indicator_pheno.size(); ++i2) {
        if (indicator_cvt[i2] == 0) {
          continue;
        }
        for (vector<int>::size_type j2 = 0; j2 < n_ph; j2++) {

          if (indicator_pheno[i2][j2] == 1) {
            if (indicator_pheno[i1][j1] == 1) {
              gsl_matrix_set(
                  H_oo, c_obs1, c_obs2,
                  gsl_matrix_get(H, c_obs1 + c_miss1, c_obs2 + c_miss2));
            } else {
              gsl_matrix_set(
                  H_mo, c_miss1, c_obs2,
                  gsl_matrix_get(H, c_obs1 + c_miss1, c_obs2 + c_miss2));
            }
            c_obs2++;
          } else {
            c_miss2++;
          }
        }
      }

      if (indicator_pheno[i1][j1] == 1) {
        c_obs1++;
      } else {
        c_miss1++;
      }
    }
  }

  // Do LU decomposition of H_oo.
  int sig;
  gsl_permutation *pmt = gsl_permutation_alloc(np_obs);
  LUDecomp(H_oo, pmt, &sig);

  // Obtain y_obs=y_full-y_hat.
  // Add the fixed effects part to y_miss: y_miss=y_hat.
  c_obs1 = 0;
  c_miss1 = 0;
  for (vector<int>::size_type i = 0; i < indicator_pheno.size(); ++i) {
    if (indicator_cvt[i] == 0) {
      continue;
    }

    for (vector<int>::size_type j = 0; j < n_ph; ++j) {
      if (indicator_pheno[i][j] == 1) {
        gsl_vector_set(y_obs, c_obs1, gsl_matrix_get(Y_full, i, j) -
                                          gsl_matrix_get(Y_hat, i, j));
        c_obs1++;
      } else {
        gsl_vector_set(y_miss, c_miss1, gsl_matrix_get(Y_hat, i, j));
        c_miss1++;
      }
    }
  }

  LUSolve(H_oo, pmt, y_obs, Hiy);

  gsl_blas_dgemv(CblasNoTrans, 1.0, H_mo, Hiy, 1.0, y_miss);

  // Put back predicted y_miss to Y_full.
  c_miss1 = 0;
  for (vector<int>::size_type i = 0; i < indicator_pheno.size(); ++i) {
    if (indicator_cvt[i] == 0) {
      continue;
    }

    for (vector<int>::size_type j = 0; j < n_ph; ++j) {
      if (indicator_pheno[i][j] == 0) {
        gsl_matrix_set(Y_full, i, j, gsl_vector_get(y_miss, c_miss1));
        c_miss1++;
      }
    }
  }

  // Free matrices.
  gsl_vector_free(y_obs);
  gsl_vector_free(y_miss);
  gsl_matrix_free(H_oo);
  gsl_matrix_free(H_mo);
  gsl_vector_free(Hiy);

  return;
}
