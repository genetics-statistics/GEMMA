/*
    Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright © 2011-2017, Xiang Zhou
    Copyright © 2017, Peter Carbonetto
    Copyright © 2017-2020, Pjotr Prins

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

#include <assert.h>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"

#include "debug.h"
// #include "eigenlib.h"
#include "fastblas.h"
#include "gzstream.h"
#include "gemma_io.h"
#include "lapack.h"
#include "mathfunc.h"

using namespace std;

// Print progress bar.
void ProgressBar(string str, double p, double total, double ratio) {
  assert(p<=total);
  assert(p>=0);
  if (total <= 0.0) return;
  const double progress = (100.0 * p / total);
  const uint barsize = (int)(progress / 2.0); // characters
  // cout << barsize << endl;
  // cout << str << " ";
  // cout << p << "/" << total << endl;
  assert(barsize < 101); // corrupted data somehow
  if (barsize > 0) {
    cout << std::string(barsize,'=');
  }
  cout << std::string(50-barsize,' ');
  cout << setprecision(0) << fixed << " " << progress << "%";
  if (ratio != -1.0)
    cout << setprecision(2) << "    " << ratio;
  cout << "\r" << flush;
}


bool isBlankLine(char const *line) {
  for (char const *cp = line; *cp; ++cp) {
    if (!isspace(*cp))
      return false;
  }
  return true;
}

typedef vector<const char *> token_list;

token_list tokenize_whitespace(const string line, uint num, const char *infilen) {
  token_list v;
  v.reserve(num);
  auto token = strtok_safe2((char *)line.c_str(), " ,\t",infilen);
  while (token) {
    const char *token2 = strndup(token,256);
    v.push_back(token2);
    // cout << token << ",";
    token = strtok(NULL, " ,\t");
  }
  return v;
}

// Faster version of safeGetline because of less copying and char
// iteration. Note behaviour differs somewhat when it comes to eol. I
// think the version had to deal with files from differing platforms.
// You can still run that with -legacy switch.
inline istream &safe_get_line(istream &is, string &t) {
  if (is_legacy_mode()) return safeGetline(is,t);

  std::getline(is,t);
  // if(!is.fail);
  return is;
}

bool isBlankLine(std::string const &line) { return isBlankLine(line.c_str()); }

// In case files are ended with "\r" or "\r\n". safeGetline fetches
// lines from a stream and returns them in t. It returns stream so it
// can be tested for eof. This function is a bottleneck in legacy gemma
// and can be replaced with safe_get_line.
std::istream &safeGetline(std::istream &is, std::string &t) {
  t.clear();

  // The characters in the stream are read one-by-one using a
  // std::streambuf. That is faster than reading them one-by-one
  // using the std::istream. Code that uses streambuf this way must
  // be guarded by a sentry object. The sentry object performs
  // various tasks, such as thread synchronization and updating the
  // stream state.
  std::istream::sentry se(is, true);
  std::streambuf *sb = is.rdbuf();

  for (;;) {
    int c = sb->sbumpc();
    switch (c) {
    case '\n':
      return is;
    case '\r':
      if (sb->sgetc() == '\n')
        sb->sbumpc();
      return is;
    case EOF:

      // Also handle the case when the last line has no line
      // ending.
      if (t.empty())
        is.setstate(std::ios::eofbit);
      return is;
    default:
      t += (char)c;
    }
  }
}

// Read SNP file. A single column of SNP names.
bool ReadFile_snps(const string file_snps, set<string> &setSnps) {
  debug_msg("entered");
  setSnps.clear();

  igzstream infile(file_snps.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open snps file: " << file_snps << endl;
    return false;
  }

  string line;
  char *ch_ptr;

  while (getline(infile, line)) {
    ch_ptr = strtok((char *)line.c_str(), " ,\t");
    enforce_msg(ch_ptr,"Problem reading SNP file");
    setSnps.insert(ch_ptr);
  }

  infile.close();
  infile.clear();

  return true;
}

// Read SNP file using a header. The header determines how the
// values for each row are parsed. A valid header can be, for example,
// RS POS CHR
bool ReadFile_snps_header(const string &file_snps, set<string> &setSnps) {
  debug_msg("entered");
  setSnps.clear();

  igzstream infile(file_snps.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open snps file: " << file_snps << endl;
    return false;
  }

  string line, rs, chr, pos;
  char *ch_ptr;

  // Read header.
  HEADER header;
  safeGetline(infile, line).eof();
  ReadHeader_io(line, header);

  if (header.rs_col == 0 && (header.chr_col == 0 || header.pos_col == 0)) {
    cout << "missing rs id in the header" << endl;
  }

  while (!safeGetline(infile, line).eof()) {
    if (isBlankLine(line)) {
      continue;
    }
    ch_ptr = strtok((char *)line.c_str(), " ,\t");
    enforce_msg(ch_ptr,"Problem reading SNP header");

    for (size_t i = 0; i < header.coln; i++) {
      enforce_msg(ch_ptr,"Problem reading SNP file");
      if (header.rs_col != 0 && header.rs_col == i + 1) {
        rs = ch_ptr;
      }
      if (header.chr_col != 0 && header.chr_col == i + 1) {
        chr = ch_ptr;
      }
      if (header.pos_col != 0 && header.pos_col == i + 1) {
        pos = ch_ptr;
      }

      ch_ptr = strtok(NULL, " ,\t");
    }

    if (header.rs_col == 0) {
      rs = chr + ":" + pos;
    }

    setSnps.insert(rs);
  }

  infile.close();
  infile.clear();

  return true;
}

// Read log file.
bool ReadFile_log(const string &file_log, double &pheno_mean) {
  debug_msg("ReadFile_log");
  ifstream infile(file_log.c_str(), ifstream::in);
  if (!infile) {
    cout << "error! fail to open log file: " << file_log << endl;
    return false;
  }

  string line;
  char *ch_ptr;
  size_t flag = 0;

  auto infilen = file_log.c_str();
  while (getline(infile, line)) {
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",infilen);
    ch_ptr = strtok(NULL, " ,\t");

    if (ch_ptr != NULL && strcmp(ch_ptr, "estimated") == 0) {
      ch_ptr = strtok(NULL, " ,\t");
      if (ch_ptr != NULL && strcmp(ch_ptr, "mean") == 0) {
        ch_ptr = strtok(NULL, " ,\t");
        if (ch_ptr != NULL && strcmp(ch_ptr, "=") == 0) {
          ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
          pheno_mean = atof(ch_ptr);
          flag = 1;
        }
      }
    }

    if (flag == 1) {
      break;
    }
  }

  infile.close();
  infile.clear();

  return true;
}

// Read bimbam annotation file which consists of rows of SNP, POS and CHR
bool ReadFile_anno(const string &file_anno, map<string, string> &mapRS2chr,
                   map<string, long int> &mapRS2bp,
                   map<string, double> &mapRS2cM) {
  debug_msg("ReadFile_anno");
  mapRS2chr.clear();
  mapRS2bp.clear();

  ifstream infile(file_anno.c_str(), ifstream::in);
  if (!infile) {
    cout << "error opening annotation file: " << file_anno << endl;
    return false;
  }

  string line;

  while (!safeGetline(infile, line).eof()) {
    const char *ch_ptr = strtok((char *)line.c_str(), " ,\t");
    enforce_str(ch_ptr, line + " Bad RS format");
    const string rs = ch_ptr;
    enforce_str(rs != "", line + " Bad RS format");

    ch_ptr = strtok(NULL, " ,\t");
    enforce_str(ch_ptr, line + " Bad format");
    long b_pos;
    if (strcmp(ch_ptr, "NA") == 0) {
      b_pos = -9;
    } else {
      b_pos = atol(ch_ptr);
    }
    enforce_str(b_pos,line + " Bad pos format (is zero)");

    string chr;
    ch_ptr = strtok(NULL, " ,\t");
    if (ch_ptr == NULL || strcmp(ch_ptr, "NA") == 0) {
      chr = "-9";
    } else {
      chr = ch_ptr;
      enforce_str(chr != "", line + " Bad chr format");
    }

    double cM;
    ch_ptr = strtok(NULL, " ,\t");
    if (ch_ptr == NULL || strcmp(ch_ptr, "NA") == 0) {
      cM = -9;
    } else {
      cM = atof(ch_ptr);
      enforce_str(b_pos, line + "Bad cM format (is zero)");
    }

    mapRS2chr[rs] = chr;
    mapRS2bp[rs] = b_pos;
    mapRS2cM[rs] = cM;
  }

  infile.close();
  infile.clear();

  return true;
}

// Read 1 column of phenotype.
bool ReadFile_column(const string &file_pheno, vector<int> &indicator_idv,
                     vector<double> &pheno, const int &p_column) {
  debug_msg("entered");
  indicator_idv.clear();
  pheno.clear();

  igzstream infile(file_pheno.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open phenotype file: " << file_pheno << endl;
    return false;
  }

  string line;
  char *ch_ptr;

  string id;
  double p;
  auto infilen = file_pheno.c_str();
  while (!safeGetline(infile, line).eof()) {
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",infilen);
    for (int i = 0; i < (p_column - 1); ++i) {
      ch_ptr = strtok(NULL, " ,\t");
    }
    enforce_msg(ch_ptr,"Problem reading PHENO column");
    if (strcmp(ch_ptr, "NA") == 0) {
      indicator_idv.push_back(0);
      pheno.push_back(-9);
    } else {
      // Pheno is different from pimass2.
      p = atof(ch_ptr);
      indicator_idv.push_back(1);
      pheno.push_back(p);
    }
  }

  infile.close();
  infile.clear();

  return true;
}

// Read bimbam phenotype file, p_column=1, 2,...
bool ReadFile_pheno(const string &file_pheno,
                    vector<vector<int>> &indicator_pheno,
                    vector<vector<double>> &pheno,
                    const vector<size_t> &p_column) {
  debug_msg("entered");
  indicator_pheno.clear();
  pheno.clear();

  igzstream infile(file_pheno.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open phenotype file: " << file_pheno << endl;
    return false;
  }

  string line;
  char *ch_ptr;

  string id;
  double p;

  vector<double> pheno_row;
  vector<int> ind_pheno_row;

  size_t p_max = *max_element(p_column.begin(), p_column.end());
  map<size_t, size_t> mapP2c;
  for (size_t i = 0; i < p_column.size(); i++) {
    mapP2c[p_column[i]] = i;
    pheno_row.push_back(-9);
    ind_pheno_row.push_back(0);
  }

  while (!safeGetline(infile, line).eof()) {
    ch_ptr = strtok((char *)line.c_str(), " ,\t");
    size_t i = 0;
    while (i < p_max) {
      enforce_msg(ch_ptr,"Number of phenotypes in pheno file do not match phenotypes in geno file");
      if (mapP2c.count(i + 1) != 0) {
        if (strcmp(ch_ptr, "NA") == 0) {
          ind_pheno_row[mapP2c[i + 1]] = 0;
          pheno_row[mapP2c[i + 1]] = -9;
        } else {
          p = atof(ch_ptr);
          ind_pheno_row[mapP2c[i + 1]] = 1;
          pheno_row[mapP2c[i + 1]] = p;
        }
      }
      i++;
      ch_ptr = strtok(NULL, " ,\t");
    }

    indicator_pheno.push_back(ind_pheno_row);
    pheno.push_back(pheno_row);
  }

  infile.close();
  infile.clear();

  return true;
}

bool ReadFile_cvt(const string &file_cvt, vector<int> &indicator_cvt,
                  vector<vector<double>> &cvt, size_t &n_cvt) {
  debug_msg("entered");
  indicator_cvt.clear();

  ifstream infile(file_cvt.c_str(), ifstream::in);
  if (!infile) {
    cout << "error! fail to open covariates file: " << file_cvt << endl;
    return false;
  }

  string line;
  char *ch_ptr;
  double d;

  int flag_na = 0;

  while (!safeGetline(infile, line).eof()) {
    vector<double> v_d;
    flag_na = 0;
    ch_ptr = strtok((char *)line.c_str(), " ,\t");
    while (ch_ptr != NULL) {
      if (strcmp(ch_ptr, "NA") == 0) {
        flag_na = 1;
        d = -9;
      } else {
        d = atof(ch_ptr);
      }

      v_d.push_back(d);
      ch_ptr = strtok(NULL, " ,\t");
    }
    if (flag_na == 0) {
      indicator_cvt.push_back(1);
    } else {
      indicator_cvt.push_back(0);
    }
    cvt.push_back(v_d);
  }

  if (indicator_cvt.empty()) {
    n_cvt = 0;
  } else {
    flag_na = 0;
    for (vector<int>::size_type i = 0; i < indicator_cvt.size(); ++i) {
      if (indicator_cvt[i] == 0) {
        continue;
      }

      if (flag_na == 0) {
        flag_na = 1;
        n_cvt = cvt[i].size();
      }
      if (flag_na != 0 && n_cvt != cvt[i].size()) {
        cout << "error! number of covariates in row " << i
             << " do not match other rows." << endl;
        return false;
      }
    }
  }

  infile.close();
  infile.clear();

  return true;
}

// Read .bim file.
bool ReadFile_bim(const string &file_bim, vector<SNPINFO> &snpInfo) {
  debug_msg("entered");
  snpInfo.clear();

  ifstream infile(file_bim.c_str(), ifstream::in);
  if (!infile) {
    cout << "error opening .bim file: " << file_bim << endl;
    return false;
  }

  string line;
  char *ch_ptr;

  string rs;
  long int b_pos;
  string chr;
  double cM;
  string major;
  string minor;

  auto infilen = file_bim.c_str();
  while (getline(infile, line)) {
    ch_ptr = strtok_safe2((char *)line.c_str(), " \t",infilen);
    chr = ch_ptr;
    ch_ptr = strtok_safe2(NULL, " \t",infilen);
    rs = ch_ptr;
    ch_ptr = strtok_safe2(NULL, " \t",infilen);
    cM = atof(ch_ptr);
    ch_ptr = strtok_safe2(NULL, " \t",infilen);
    b_pos = atol(ch_ptr);
    ch_ptr = strtok_safe2(NULL, " \t",infilen);
    minor = ch_ptr;
    ch_ptr = strtok_safe2(NULL, " \t",infilen);
    major = ch_ptr;

    SNPINFO sInfo = {chr, rs, cM, b_pos, minor, major, 0, -9, -9, 0, 0, 0};
    snpInfo.push_back(sInfo);
  }

  infile.close();
  infile.clear();
  return true;
}

// Read .fam file (ignored with -p phenotypes switch)
bool ReadFile_fam(const string &file_fam, vector<vector<int>> &indicator_pheno,
                  vector<vector<double>> &pheno, map<string, int> &mapID2num,
                  const vector<size_t> &p_column) {
  debug_msg("entered");
  indicator_pheno.clear();
  pheno.clear();
  mapID2num.clear();

  igzstream infile(file_fam.c_str(), igzstream::in);
  if (!infile) {
    cout << "error opening .fam file: " << file_fam << endl;
    return false;
  }

  string line;
  char *ch_ptr;

  string id;
  int c = 0;
  double p;

  vector<double> pheno_row;
  vector<int> ind_pheno_row;

  size_t p_max = *max_element(p_column.begin(), p_column.end());
  map<size_t, size_t> mapP2c;
  for (size_t i = 0; i < p_column.size(); i++) {
    mapP2c[p_column[i]] = i;
    pheno_row.push_back(-9);
    ind_pheno_row.push_back(0);
  }

  auto infilen = file_fam.c_str();
  while (!safeGetline(infile, line).eof()) {
    ch_ptr = strtok_safe2((char *)line.c_str(), " \t",infilen);
    ch_ptr = strtok_safe2(NULL, " \t",infilen);
    id = ch_ptr;
    ch_ptr = strtok_safe2(NULL, " \t",infilen);
    ch_ptr = strtok_safe2(NULL, " \t",infilen);
    ch_ptr = strtok_safe2(NULL, " \t",infilen);
    ch_ptr = strtok_safe2(NULL, " \t",infilen);

    size_t i = 0;
    while (i < p_max) {
      if (mapP2c.count(i + 1) != 0) {
        enforce_msg(ch_ptr,"Problem reading FAM file (phenotypes do not match geno file)");

        if (strcmp(ch_ptr, "NA") == 0) {
          ind_pheno_row[mapP2c[i + 1]] = 0;
          pheno_row[mapP2c[i + 1]] = -9;
        } else {
          p = atof(ch_ptr);

          if (p == -9) {
            ind_pheno_row[mapP2c[i + 1]] = 0;
            pheno_row[mapP2c[i + 1]] = -9;
          } else {
            ind_pheno_row[mapP2c[i + 1]] = 1;
            pheno_row[mapP2c[i + 1]] = p;
          }
        }
      }
      i++;
      ch_ptr = strtok(NULL, " ,\t");
    }

    indicator_pheno.push_back(ind_pheno_row);
    pheno.push_back(pheno_row);

    mapID2num[id] = c;
    c++;
  }

  infile.close();
  infile.clear();
  return true;
}

// Read bimbam mean genotype file, the first time, to obtain #SNPs for
// analysis (ns_test) and total #SNP (ns_total).
bool ReadFile_geno(const string &file_geno, const set<string> &setSnps,
                   const gsl_matrix *W, vector<int> &indicator_idv,
                   vector<int> &indicator_snp, const double &maf_level,
                   const double &miss_level, const double &hwe_level,
                   const double &r2_level, map<string, string> &mapRS2chr,
                   map<string, long int> &mapRS2bp,
                   map<string, double> &mapRS2cM, vector<SNPINFO> &snpInfo,
                   size_t &ns_test) {
  debug_msg("entered");
  indicator_snp.clear();
  snpInfo.clear();

  igzstream infile(file_geno.c_str(), igzstream::in);
  if (!infile) {
    cout << "error reading genotype file:" << file_geno << endl;
    return false;
  }

  gsl_vector *genotype = gsl_vector_safe_alloc(W->size1);
  gsl_vector *genotype_miss = gsl_vector_safe_alloc(W->size1);
  gsl_matrix *WtW = gsl_matrix_safe_alloc(W->size2, W->size2);
  gsl_matrix *WtWi = gsl_matrix_safe_alloc(W->size2, W->size2);
  gsl_vector *Wtx = gsl_vector_safe_alloc(W->size2);
  gsl_vector *WtWiWtx = gsl_vector_safe_alloc(W->size2);
  gsl_permutation *pmt = gsl_permutation_alloc(W->size2);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
  int sig;
  LUDecomp(WtW, pmt, &sig);

  LUInvert(WtW, pmt, WtWi); // @@

  double v_x, v_w;
  int c_idv = 0;

  string line;
  char *ch_ptr;

  string rs;
  long int b_pos;
  string chr;
  string major;
  string minor;
  double cM;
  size_t file_pos;

  double maf, geno, geno_old;
  size_t n_miss;
  size_t n_0, n_1, n_2;
  int flag_poly;

  int ni_total = indicator_idv.size();
  int ni_test = 0;
  for (int i = 0; i < ni_total; ++i) {
    ni_test += indicator_idv[i];
  }
  ns_test = 0;

  file_pos = 0;
  auto count_warnings = 0;
  auto infilen = file_geno.c_str();
  while (!safe_get_line(infile, line).eof()) {
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",infilen);
    rs = ch_ptr;
    ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
    minor = ch_ptr;
    ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
    major = ch_ptr;

    if (setSnps.size() != 0 && setSnps.count(rs) == 0) {
      // if SNP in geno but not in -snps we add an missing value
      SNPINFO sInfo = {"-9", rs, -9, -9, minor, major,
                       0,    -9, -9, 0,  0,     file_pos};
      snpInfo.push_back(sInfo);
      indicator_snp.push_back(0);

      file_pos++;
      continue;
    }

    if (mapRS2bp.count(rs) == 0) {
      if (is_debug_mode() && count_warnings++ < 10) {
        std::string msg = "Can't figure out position for ";
        msg += rs;
        debug_msg(msg);
        if (count_warnings == 10)
          debug_msg("Skipping similar warnings");
      }
      chr = "-9";
      b_pos = -9;
      cM = -9;
    } else {
      b_pos = mapRS2bp[rs];
      chr = mapRS2chr[rs];
      cM = mapRS2cM[rs];
    }

    // Start on a new marker/SNP
    maf = 0;
    n_miss = 0;
    flag_poly = 0;
    geno_old = -9;
    n_0 = 0;
    n_1 = 0;
    n_2 = 0;
    c_idv = 0;
    gsl_vector_set_zero(genotype_miss);
    auto infilen = file_geno.c_str();
    for (int i = 0; i < ni_total; ++i) {
      ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
      if (indicator_idv[i] == 0)
        continue;

      enforce_msg(ch_ptr,"Problem reading geno file (not enough genotypes in line)");
      if (strcmp(ch_ptr, "NA") == 0) {
        gsl_vector_set(genotype_miss, c_idv, 1);
        n_miss++;
        c_idv++;
        continue;
      }

      geno = atof(ch_ptr);
      if (geno >= 0 && geno <= 0.5) {
        n_0++;
      }
      if (geno > 0.5 && geno < 1.5) {
        n_1++;
      }
      if (geno >= 1.5 && geno <= 2.0) {
        n_2++;
      }

      gsl_vector_set(genotype, c_idv, geno);

      // going through genotypes with 0.0 < geno < 2.0
      if (flag_poly == 0) { // first init in marker
        geno_old = geno;    // set geno_old (double) to previous genotype
        flag_poly = 2;      // initialized state
      }
      if (flag_poly == 2 && geno != geno_old) {
        flag_poly = 1;      // genotypes differ
      }

      maf += geno;

      c_idv++;
    }
    maf /= 2.0 * (double)(ni_test - n_miss);

    SNPINFO sInfo = {chr,    rs,
                     cM,     b_pos,
                     minor,  major,
                     n_miss, (double)n_miss / (double)ni_test,
                     maf,    ni_test - n_miss,
                     0,      file_pos};
    snpInfo.push_back(sInfo);
    file_pos++;

    // -miss flag
    if ((double)n_miss / (double)ni_test > miss_level) {
      indicator_snp.push_back(0);
      continue;
    }

    // -maf flag
    if ((maf < maf_level || maf > (1.0 - maf_level)) && maf_level != -1) {
      indicator_snp.push_back(0);
      continue;
    }

    // remove genotype lines that are identical to the one read before
    if (flag_poly != 1) {
      indicator_snp.push_back(0);
      continue;
    }

    // -hwe flag
    if (hwe_level != 0 && maf_level != -1) {
      if (CalcHWE(n_0, n_2, n_1) < hwe_level) {
        indicator_snp.push_back(0);
        continue;
      }
    }


    // -r2 flag
    for (size_t i = 0; i < genotype->size; ++i) {
      if (gsl_vector_get(genotype_miss, i) == 1) {
        geno = maf * 2.0;
        gsl_vector_set(genotype, i, geno);
      }
    }

    gsl_blas_dgemv(CblasTrans, 1.0, W, genotype, 0.0, Wtx);
    gsl_blas_dgemv(CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
    gsl_blas_ddot(genotype, genotype, &v_x);
    gsl_blas_ddot(Wtx, WtWiWtx, &v_w);

    // Filter SNP if it is correlated with covariates W, unless W has
    // only one column, of 1s (-r2 flag)
    if (W->size2 != 1 && v_w / v_x > r2_level) {
      indicator_snp.push_back(0);
      continue;
    }

    indicator_snp.push_back(1);
    ns_test++;
  }

  gsl_vector_free(genotype);
  gsl_vector_free(genotype_miss);
  gsl_matrix_free(WtW);
  gsl_matrix_free(WtWi);
  gsl_vector_free(Wtx);
  gsl_vector_free(WtWiWtx);
  gsl_permutation_free(pmt);

  infile.close();
  infile.clear();

  return true;
}

// Read bed file, the first time.
bool ReadFile_bed(const string &file_bed, const set<string> &setSnps,
                  const gsl_matrix *W, vector<int> &indicator_idv,
                  vector<int> &indicator_snp, vector<SNPINFO> &snpInfo,
                  const double &maf_level, const double &miss_level,
                  const double &hwe_level, const double &r2_level,
                  size_t &ns_test) {
  debug_msg("entered");
  indicator_snp.clear();
  size_t ns_total = snpInfo.size();

  ifstream infile(file_bed.c_str(), ios::binary);
  if (!infile) {
    cout << "error reading bed file:" << file_bed << endl;
    return false;
  }

  gsl_vector *genotype = gsl_vector_safe_alloc(W->size1);
  gsl_vector *genotype_miss = gsl_vector_safe_alloc(W->size1);
  gsl_matrix *WtW = gsl_matrix_safe_alloc(W->size2, W->size2);
  gsl_matrix *WtWi = gsl_matrix_safe_alloc(W->size2, W->size2);
  gsl_vector *Wtx = gsl_vector_safe_alloc(W->size2);
  gsl_vector *WtWiWtx = gsl_vector_safe_alloc(W->size2);
  gsl_permutation *pmt = gsl_permutation_alloc(W->size2);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
  int sig;
  LUDecomp(WtW, pmt, &sig);
  LUInvert(WtW, pmt, WtWi);

  double v_x, v_w, geno;
  size_t c_idv = 0;

  char ch[1];
  bitset<8> b;

  size_t ni_total = indicator_idv.size();
  size_t ni_test = 0;
  for (size_t i = 0; i < ni_total; ++i) {
    ni_test += indicator_idv[i];
  }
  ns_test = 0;

  // Calculate n_bit and c, the number of bit for each snp.
  size_t n_bit;
  if (ni_total % 4 == 0) {
    n_bit = ni_total / 4;
  } else {
    n_bit = ni_total / 4 + 1;
  }

  // Ignore the first three magic numbers.
  for (int i = 0; i < 3; ++i) {
    infile.read(ch, 1);
    b = ch[0];
  }

  double maf;
  size_t n_miss;
  size_t n_0, n_1, n_2, c;

  // Start reading snps and doing association test.
  for (size_t t = 0; t < ns_total; ++t) {

    // n_bit, and 3 is the number of magic numbers.
    infile.seekg(t * n_bit + 3);

    if (setSnps.size() != 0 && setSnps.count(snpInfo[t].rs_number) == 0) {
      snpInfo[t].n_miss = -9;
      snpInfo[t].missingness = -9;
      snpInfo[t].maf = -9;
      snpInfo[t].file_position = t;
      indicator_snp.push_back(0);
      continue;
    }

    // Read genotypes.
    c = 0;
    maf = 0.0;
    n_miss = 0;
    n_0 = 0;
    n_1 = 0;
    n_2 = 0;
    c_idv = 0;
    gsl_vector_set_zero(genotype_miss);
    for (size_t i = 0; i < n_bit; ++i) {
      infile.read(ch, 1);
      b = ch[0];

      // Minor allele homozygous: 2.0; major: 0.0;
      for (size_t j = 0; j < 4; ++j) {
        if ((i == (n_bit - 1)) && c == ni_total) {
          break;
        }
        if (indicator_idv[c] == 0) {
          c++;
          continue;
        }
        c++;

        if (b[2 * j] == 0) {
          if (b[2 * j + 1] == 0) {
            gsl_vector_set(genotype, c_idv, 2.0);
            maf += 2.0;
            n_2++;
          } else {
            gsl_vector_set(genotype, c_idv, 1.0);
            maf += 1.0;
            n_1++;
          }
        } else {
          if (b[2 * j + 1] == 1) {
            gsl_vector_set(genotype, c_idv, 0.0);
            maf += 0.0;
            n_0++;
          } else {
            gsl_vector_set(genotype_miss, c_idv, 1);
            n_miss++;
          }
        }
        c_idv++;
      }
    }
    maf /= 2.0 * (double)(ni_test - n_miss);

    snpInfo[t].n_miss = n_miss;
    snpInfo[t].missingness = (double)n_miss / (double)ni_test;
    snpInfo[t].maf = maf;
    snpInfo[t].n_idv = ni_test - n_miss;
    snpInfo[t].n_nb = 0;
    snpInfo[t].file_position = t;

    if ((double)n_miss / (double)ni_test > miss_level) {
      indicator_snp.push_back(0);
      continue;
    }

    if ((maf < maf_level || maf > (1.0 - maf_level)) && maf_level != -1) {
      indicator_snp.push_back(0);
      continue;
    }

    if ((n_0 + n_1) == 0 || (n_1 + n_2) == 0 || (n_2 + n_0) == 0) {
      indicator_snp.push_back(0);
      continue;
    }

    if (hwe_level != 0 && maf_level != -1) {
      if (CalcHWE(n_0, n_2, n_1) < hwe_level) {
        indicator_snp.push_back(0);
        continue;
      }
    }

    // Filter SNP if it is correlated with W unless W has
    // only one column, of 1s.
    for (size_t i = 0; i < genotype->size; ++i) {
      if (gsl_vector_get(genotype_miss, i) == 1) {
        geno = maf * 2.0;
        gsl_vector_set(genotype, i, geno);
      }
    }

    gsl_blas_dgemv(CblasTrans, 1.0, W, genotype, 0.0, Wtx);
    gsl_blas_dgemv(CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
    gsl_blas_ddot(genotype, genotype, &v_x);
    gsl_blas_ddot(Wtx, WtWiWtx, &v_w);

    if (W->size2 != 1 && v_w / v_x > r2_level) {
      indicator_snp.push_back(0);
      continue;
    }

    indicator_snp.push_back(1);
    ns_test++;
  }

  gsl_vector_free(genotype);
  gsl_vector_free(genotype_miss);
  gsl_matrix_free(WtW);
  gsl_matrix_free(WtWi);
  gsl_vector_free(Wtx);
  gsl_vector_free(WtWiWtx);
  gsl_permutation_free(pmt);

  infile.close();
  infile.clear();

  return true;
}

// Read the genotype for one SNP; remember to read empty lines.
// Geno stores original genotypes without centering.
// Missing values are replaced by mean.
bool Bimbam_ReadOneSNP(const size_t inc, const vector<int> &indicator_idv,
                       igzstream &infile, gsl_vector *geno, double &geno_mean) {
  debug_msg("entered");
  size_t ni_total = indicator_idv.size();

  string line;
  char *ch_ptr;
  bool flag = false;

  for (size_t i = 0; i < inc; i++) {
    safeGetline(infile, line).eof();
  }

  if (!safeGetline(infile, line).eof()) {
    ch_ptr = strtok_safe((char *)line.c_str(), " ,\t");
    ch_ptr = strtok_safe(NULL, " ,\t");
    ch_ptr = strtok_safe(NULL, " ,\t");

    geno_mean = 0.0;
    double d;
    size_t c_idv = 0;
    vector<size_t> geno_miss;

    for (size_t i = 0; i < ni_total; ++i) {
      ch_ptr = strtok_safe(NULL, " ,\t");
      if (indicator_idv[i] == 0) {
        continue;
      }

      if (strcmp(ch_ptr, "NA") == 0) {
        geno_miss.push_back(c_idv);
      } else {
        d = atof(ch_ptr);
        gsl_vector_set(geno, c_idv, d);
        geno_mean += d;
      }
      c_idv++;
    }

    geno_mean /= (double)(c_idv - geno_miss.size());

    for (size_t i = 0; i < geno_miss.size(); ++i) {
      gsl_vector_set(geno, geno_miss[i], geno_mean);
    }
    flag = true;
  }

  return flag;
}

// For PLINK, store SNPs as double too.
void Plink_ReadOneSNP(const int pos, const vector<int> &indicator_idv,
                      ifstream &infile, gsl_vector *geno, double &geno_mean) {
  debug_msg("entered");
  size_t ni_total = indicator_idv.size(), n_bit;
  if (ni_total % 4 == 0) {
    n_bit = ni_total / 4;
  } else {
    n_bit = ni_total / 4 + 1;
  }

  // n_bit, and 3 is the number of magic numbers.
  infile.seekg(pos * n_bit + 3);

  // Read genotypes.
  char ch[1];
  bitset<8> b;

  geno_mean = 0.0;
  size_t c = 0, c_idv = 0;
  vector<size_t> geno_miss;

  for (size_t i = 0; i < n_bit; ++i) {
    infile.read(ch, 1);
    b = ch[0];

    // Minor allele homozygous: 2.0; major: 0.0.
    for (size_t j = 0; j < 4; ++j) {
      if ((i == (n_bit - 1)) && c == ni_total) {
        break;
      }
      if (indicator_idv[c] == 0) {
        c++;
        continue;
      }
      c++;

      if (b[2 * j] == 0) {
        if (b[2 * j + 1] == 0) {
          gsl_vector_set(geno, c_idv, 2);
          geno_mean += 2.0;
        } else {
          gsl_vector_set(geno, c_idv, 1);
          geno_mean += 1.0;
        }
      } else {
        if (b[2 * j + 1] == 1) {
          gsl_vector_set(geno, c_idv, 0);
          geno_mean += 0.0;
        } else {
          geno_miss.push_back(c_idv);
        }
      }

      c_idv++;
    }
  }

  geno_mean /= (double)(c_idv - geno_miss.size());

  for (size_t i = 0; i < geno_miss.size(); ++i) {
    gsl_vector_set(geno, geno_miss[i], geno_mean);
  }

  return;
}

void ReadFile_kin(const string &file_kin, vector<int> &indicator_idv,
                  map<string, int> &mapID2num, const size_t k_mode, bool &error,
                  gsl_matrix *G) {
  debug_msg("entered");
  igzstream infile(file_kin.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open kinship file: " << file_kin << endl;
    error = true;
    return;
  }

  size_t ni_total = indicator_idv.size();

  gsl_matrix_set_zero(G);

  string line;
  char *ch_ptr;
  double d;

  if (k_mode == 1) {
    size_t i_test = 0, i_total = 0, j_test = 0, j_total = 0;
    while (getline(infile, line)) {
      if (i_total == ni_total) {
        fail_msg("number of rows in the kinship file is larger than the number of phenotypes");
      }

      if (indicator_idv[i_total] == 0) {
        i_total++;
        continue;
      }

      j_total = 0;
      j_test = 0;
      ch_ptr = strtok((char *)line.c_str(), " ,\t");
      while (ch_ptr != NULL) {
        if (j_total == ni_total) {
          fail_msg(string("number of columns in the kinship file is larger than the number of individuals for row = ")+to_string(i_total));
        }

        d = atof(ch_ptr);
        if (indicator_idv[j_total] == 1) {
          gsl_matrix_set(G, i_test, j_test, d);
          j_test++;
        }
        j_total++;

        ch_ptr = strtok(NULL, " ,\t");
      }
      if (j_total != ni_total) {
        string msg = "number of columns in the kinship file does not match the number of individuals for row = " + to_string( i_total );
        fail_msg(msg);
      }
      i_total++;
      i_test++;
    }
    if (i_total != ni_total) {
      fail_msg("number of rows in the kinship file does not match the number of individuals.");
    }
  } else {
    map<size_t, size_t> mapID2ID;
    size_t c = 0;
    for (size_t i = 0; i < indicator_idv.size(); i++) {
      if (indicator_idv[i] == 1) {
        mapID2ID[i] = c;
        c++;
      }
    }

    string id1, id2;
    double Cov_d;
    size_t n_id1, n_id2;

    auto infilen=file_kin.c_str();
    while (getline(infile, line)) {
      ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",infilen);
      id1 = ch_ptr;
      ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
      id2 = ch_ptr;
      ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
      d = atof(ch_ptr);
      if (mapID2num.count(id1) == 0 || mapID2num.count(id2) == 0) {
        continue;
      }
      if (indicator_idv[mapID2num[id1]] == 0 ||
          indicator_idv[mapID2num[id2]] == 0) {
        continue;
      }

      n_id1 = mapID2ID[mapID2num[id1]];
      n_id2 = mapID2ID[mapID2num[id2]];

      Cov_d = gsl_matrix_get(G, n_id1, n_id2);
      if (Cov_d != 0 && Cov_d != d) {
        cerr << "error! redundant and unequal terms in the "
             << "kinship file, for id1 = " << id1 << " and id2 = " << id2
             << endl;
        fail_msg("");
      } else {
        gsl_matrix_set(G, n_id1, n_id2, d);
        gsl_matrix_set(G, n_id2, n_id1, d);
      }
    }
  }

  infile.close();
  infile.clear();

  return;
}

void ReadFile_mk(const string &file_mk, vector<int> &indicator_idv,
                 map<string, int> &mapID2num, const size_t k_mode, bool &error,
                 gsl_matrix *G) {
  debug_msg("entered");
  igzstream infile(file_mk.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open file: " << file_mk << endl;
    error = true;
    return;
  }

  string file_kin, line;

  size_t i = 0;
  while (getline(infile, line)) {
    file_kin = line.c_str();
    gsl_matrix_view G_sub =
        gsl_matrix_submatrix(G, 0, i * G->size1, G->size1, G->size1);
    ReadFile_kin(file_kin, indicator_idv, mapID2num, k_mode, error,
                 &G_sub.matrix);
    i++;
  }

  infile.close();
  infile.clear();
}

void ReadFile_eigenU(const string &file_ku, bool &error, gsl_matrix *U) {
  debug_msg("entered");
  igzstream infile(file_ku.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open the U file: " << file_ku << endl;
    error = true;
    return;
  }

  size_t n_row = U->size1, n_col = U->size2, i_row = 0, i_col = 0;

  gsl_matrix_set_zero(U);

  string line;
  char *ch_ptr;
  double d;

  while (getline(infile, line)) {
    if (i_row == n_row) {
      cout << "error! number of rows in the U file is larger "
           << "than expected." << endl;
      error = true;
    }

    i_col = 0;
    ch_ptr = strtok((char *)line.c_str(), " ,\t");
    while (ch_ptr != NULL) {
      if (i_col == n_col) {
        cout << "error! number of columns in the U file "
             << "is larger than expected, for row = " << i_row << endl;
        error = true;
      }

      d = atof(ch_ptr);
      gsl_matrix_set(U, i_row, i_col, d);
      i_col++;

      ch_ptr = strtok(NULL, " ,\t");
    }

    i_row++;
  }

  infile.close();
  infile.clear();

  return;
}

void ReadFile_eigenD(const string &file_kd, bool &error, gsl_vector *eval) {
  debug_msg("entered");
  igzstream infile(file_kd.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open the D file: " << file_kd << endl;
    error = true;
    return;
  }

  size_t n_row = eval->size, i_row = 0;

  gsl_vector_set_zero(eval);

  string line;
  char *ch_ptr;
  double d;

  while (getline(infile, line)) {
    if (i_row == n_row) {
      cout << "error! number of rows in the D file is larger "
           << "than expected." << endl;
      error = true;
    }

    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",file_kd.c_str());
    d = atof(ch_ptr);

    ch_ptr = strtok(NULL, " ,\t");
    if (ch_ptr != NULL) {
      cout << "error! number of columns in the D file is larger "
           << "than expected, for row = " << i_row << endl;
      error = true;
    }

    gsl_vector_set(eval, i_row, d);

    i_row++;
  }

  infile.close();
  infile.clear();

  return;
}

// Read bimbam mean genotype file and calculate kinship matrix.
bool BimbamKin(const string file_geno, const set<string> ksnps,
               vector<int> &indicator_snp, const int k_mode,
               const int display_pace, gsl_matrix *matrix_kin,
               const bool test_nind) {
  debug_msg("entered");
  auto infilen = file_geno.c_str();
  igzstream infile(infilen, igzstream::in);
  enforce_msg(infilen, "error reading genotype file");

  size_t n_miss;
  double geno_mean, geno_var;

  // setKSnp and/or LOCO support
  bool process_ksnps = ksnps.size();

  size_t ni_total = matrix_kin->size1;
  gsl_vector *geno = gsl_vector_safe_alloc(ni_total);
  gsl_vector *geno_miss = gsl_vector_safe_alloc(ni_total);

  // Xlarge contains inds x markers
  const size_t msize = K_BATCH_SIZE;
  gsl_matrix *Xlarge = gsl_matrix_safe_alloc(ni_total, msize);
  enforce_msg(Xlarge, "allocate Xlarge");

  gsl_matrix_set_zero(Xlarge);
  write(matrix_kin,"K before");

  // For every SNP read the genotype per individual
  size_t ns_test = 0;
  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    string line;
    safeGetline(infile, line).eof();
    // if (t % display_pace == 0 || t == (indicator_snp.size() - 1)) {
    if (t % display_pace == 0) {
      ProgressBar("Reading SNPs", t, indicator_snp.size() - 1);
    }
    if (indicator_snp[t] == 0)
      continue;

    // Using regular expressions is slow:
    // std::regex_token_iterator<std::string::iterator> rend;
    // regex split_on("[,[:blank:]]+");
    // regex_token_iterator<string::iterator> tokens(line.begin(), line.end(),
    //                                               split_on, -1);
    auto tokens = tokenize_whitespace(line,ni_total+3,infilen);
    if (test_nind) {
      // ascertain the number of genotype fields match
      uint token_num = tokens.size();
      if (token_num != ni_total+3) {
        cerr << line << endl;
        cerr << token_num-3 << " != " << ni_total << endl;
        warning_msg("Columns in geno file do not match # individuals in phenotypes");
      }
      enforce_msg(token_num <= ni_total + 3,"not enough genotype fields for marker");
    }

    auto token_i = tokens.begin();
    // const char *snp = *token_i; // first field
    const char *snp = tokens[0]; // first field
    // cout << snp << "!";
    // check whether SNP is included in ksnps (used by LOCO)
    if (process_ksnps && ksnps.count(snp) == 0)
      continue;

    token_i++; // skip snp name
    token_i++; // skip nucleotide field
    token_i++; // skip nucleotide field

    // calc SNP stats
    geno_mean = 0.0;
    n_miss = 0;
    geno_var = 0.0;
    gsl_vector_set_all(geno_miss, 0);
    for (size_t i = 0; i < ni_total; ++i) {
      enforce_str(token_i != tokens.end(), line + " number of fields");
      auto field = *token_i;
      auto sfield = std::string(field);
      // cout << i << ":" << sfield << "," << endl;
      if (strncmp(field,"NA",2)==0) { // missing value
        gsl_vector_set(geno_miss, i, 0);
        n_miss++;
      } else {
        double d = atof(field);
        if (is_strict_mode() && d == 0.0)
          enforce_is_float(std::string(field));  // rule out non NA and non-float fields
        gsl_vector_set(geno, i, d);
        gsl_vector_set(geno_miss, i, 1);
        geno_mean += d;
        geno_var += d * d;
      }
      token_i++;
    }

    geno_mean /= (double)(ni_total - n_miss);
    geno_var += geno_mean * geno_mean * (double)n_miss;
    geno_var /= (double)ni_total;
    geno_var -= geno_mean * geno_mean;

    if (ns_test<1) {
      write(geno,"geno raw");
      write(geno_mean,"geno mean");
    }

    // impute missing values by plugging in the mean
    for (size_t i = 0; i < ni_total; ++i) {
      if (gsl_vector_get(geno_miss, i) == 0) {
        gsl_vector_set(geno, i, geno_mean);
      }
    }

    if (ns_test<1) write(geno,"geno imputed");

    // subtract the mean (centering genotype values)
    gsl_vector_add_constant(geno, -1.0 * geno_mean);
    if (ns_test<1) write(geno,"geno mean");

    // scale the genotypes
    if (k_mode == 2 && geno_var != 0) { // some confusion here, -gk 2
                                        // flag does this
      gsl_vector_scale(geno, 1.0 / sqrt(geno_var));
    }

    if (ns_test<1) {
      write(geno_var,"geno var");
      write(geno,"geno z-scored");
    }

    // set the SNP column ns_test to copy geno into the compute matrix Xlarge
    gsl_vector_view Xlarge_col = gsl_matrix_column(Xlarge, ns_test % msize);
    enforce_gsl(gsl_vector_memcpy(&Xlarge_col.vector, geno));

    ns_test++;

    // Every msize rows batch compute kinship matrix and return
    // by adding to matrix_kin
    if (ns_test % msize == 0) {
      fast_eigen_dgemm("N", "T", 1.0, Xlarge, Xlarge, 1.0, matrix_kin);
      gsl_matrix_set_zero(Xlarge);
      write(matrix_kin,"K updated");
    }
  }
  if (ns_test % msize != 0) { // compute last batch
    fast_eigen_dgemm("N", "T", 1.0, Xlarge, Xlarge, 1.0, matrix_kin);
    write(matrix_kin,"K updated");
  }
  ProgressBar("Reading SNPs", 100, 100);
  cout << endl;

  write(matrix_kin,"K prescaled");

  // scale the kinship matrix (ns_test is used SNPs minus 1)
  write(ns_test,"ns_test scale");
  enforce_gsl(gsl_matrix_scale(matrix_kin, 1.0 / (double)ns_test));

  write(matrix_kin,"K scaled");
  // and transpose
  // FIXME: the following is not so slow but appears to generate an
  // identical matrix

  /*
  for (size_t i = 0; i < ni_total; ++i) {
    for (size_t j = 0; j < i; ++j) {
      double d = gsl_matrix_get(matrix_kin, j, i);
      gsl_matrix_set(matrix_kin, i, j, d);
    }
  }
  write(matrix_kin,"K rotated");
  */
  // GSL is faster - and there are even faster methods
  // enforce_gsl(gsl_matrix_transpose(matrix_kin));

  gsl_vector_free(geno);
  gsl_vector_free(geno_miss);
  gsl_matrix_free(Xlarge);

  infile.close();
  infile.clear();

  return true;
}

bool PlinkKin(const string &file_bed, vector<int> &indicator_snp,
              const int k_mode, const int display_pace,
              gsl_matrix *matrix_kin) {
  debug_msg("entered");
  ifstream infile(file_bed.c_str(), ios::binary);
  if (!infile) {
    cout << "error reading bed file:" << file_bed << endl;
    return false;
  }

  char ch[1];
  bitset<8> b;

  size_t n_miss, ci_total;
  double d, geno_mean, geno_var;

  size_t ni_total = matrix_kin->size1;
  gsl_vector *geno = gsl_vector_safe_alloc(ni_total);

  size_t ns_test = 0;
  int n_bit;

  // Create a large matrix.
  const size_t msize = K_BATCH_SIZE;
  gsl_matrix *Xlarge = gsl_matrix_safe_alloc(ni_total, msize);
  gsl_matrix_set_zero(Xlarge);

  // Calculate n_bit and c, the number of bit for each snp.
  if (ni_total % 4 == 0) {
    n_bit = ni_total / 4;
  } else {
    n_bit = ni_total / 4 + 1;
  }

  // print the first three magic numbers
  for (int i = 0; i < 3; ++i) {
    infile.read(ch, 1);
    b = ch[0];
  }

  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    if (t % display_pace == 0 || t == (indicator_snp.size() - 1)) {
      ProgressBar("Reading SNPs", t, indicator_snp.size() - 1);
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
    for (int i = 0; i < n_bit; ++i) {
      infile.read(ch, 1);
      b = ch[0];

      // Minor allele homozygous: 2.0; major: 0.0.
      for (size_t j = 0; j < 4; ++j) {
        if ((i == (n_bit - 1)) && ci_total == ni_total) {
          break;
        }

        if (b[2 * j] == 0) {
          if (b[2 * j + 1] == 0) {
            gsl_vector_set(geno, ci_total, 2.0);
            geno_mean += 2.0;
            geno_var += 4.0;
          } else {
            gsl_vector_set(geno, ci_total, 1.0);
            geno_mean += 1.0;
            geno_var += 1.0;
          }
        } else {
          if (b[2 * j + 1] == 1) {
            gsl_vector_set(geno, ci_total, 0.0);
          } else {
            gsl_vector_set(geno, ci_total, -9.0);
            n_miss++;
          }
        }

        ci_total++;
      }
    }

    geno_mean /= (double)(ni_total - n_miss);
    geno_var += geno_mean * geno_mean * (double)n_miss;
    geno_var /= (double)ni_total;
    geno_var -= geno_mean * geno_mean;

    for (size_t i = 0; i < ni_total; ++i) {
      d = gsl_vector_get(geno, i);
      if (d == -9.0) {
        gsl_vector_set(geno, i, geno_mean);
      }
    }

    gsl_vector_add_constant(geno, -1.0 * geno_mean);

    if (k_mode == 2 && geno_var != 0) { // k_mode is a_mode-20
      gsl_vector_scale(geno, 1.0 / sqrt(geno_var));
    }
    gsl_vector_view Xlarge_col = gsl_matrix_column(Xlarge, ns_test % msize);
    gsl_vector_memcpy(&Xlarge_col.vector, geno);

    ns_test++;

    if (ns_test % msize == 0) {
      fast_eigen_dgemm("N", "T", 1.0, Xlarge, Xlarge, 1.0, matrix_kin);
      gsl_matrix_set_zero(Xlarge);
    }
  }

  if (ns_test % msize != 0) {
    fast_eigen_dgemm("N", "T", 1.0, Xlarge, Xlarge, 1.0, matrix_kin);
  }

  cout << endl;

  gsl_matrix_scale(matrix_kin, 1.0 / (double)ns_test);

  for (size_t i = 0; i < ni_total; ++i) {
    for (size_t j = 0; j < i; ++j) {
      d = gsl_matrix_get(matrix_kin, j, i);
      gsl_matrix_set(matrix_kin, i, j, d);
    }
  }

  gsl_vector_free(geno);
  gsl_matrix_free(Xlarge);

  infile.close();
  infile.clear();

  return true;
}

// Read bimbam mean genotype file, the second time, recode "mean"
// genotype and calculate K.
bool ReadFile_geno(const string file_geno, vector<int> &indicator_idv,
                   vector<int> &indicator_snp, gsl_matrix *UtX, gsl_matrix *K,
                   const bool calc_K) {
  debug_msg("entered");
  igzstream infile(file_geno.c_str(), igzstream::in);
  if (!infile) {
    cout << "error reading genotype file:" << file_geno << endl;
    return false;
  }

  string line;
  char *ch_ptr;

  if (calc_K == true) {
    gsl_matrix_set_zero(K);
  }

  gsl_vector *genotype = gsl_vector_safe_alloc(UtX->size1);
  gsl_vector *genotype_miss = gsl_vector_safe_alloc(UtX->size1);
  double geno, geno_mean;
  size_t n_miss;

  int ni_total = (int)indicator_idv.size();
  int ns_total = (int)indicator_snp.size();
  int ni_test = UtX->size1;
  int ns_test = UtX->size2;

  int c_idv = 0, c_snp = 0;

  auto infilen = file_geno.c_str();
  for (int i = 0; i < ns_total; ++i) {
    safeGetline(infile, line).eof();
    if (indicator_snp[i] == 0) {
      continue;
    }

    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",infilen);
    ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
    ch_ptr = strtok_safe2(NULL, " ,\t",infilen);

    c_idv = 0;
    geno_mean = 0;
    n_miss = 0;
    gsl_vector_set_zero(genotype_miss);
    for (int j = 0; j < ni_total; ++j) {
      ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
      if (indicator_idv[j] == 0) {
        continue;
      }

      if (strcmp(ch_ptr, "NA") == 0) {
        gsl_vector_set(genotype_miss, c_idv, 1);
        n_miss++;
      } else {
        geno = atof(ch_ptr);
        gsl_vector_set(genotype, c_idv, geno);
        geno_mean += geno;
      }
      c_idv++;
    }

    geno_mean /= (double)(ni_test - n_miss);

    for (size_t i = 0; i < genotype->size; ++i) {
      if (gsl_vector_get(genotype_miss, i) == 1) {
        geno = 0;
      } else {
        geno = gsl_vector_get(genotype, i);
        geno -= geno_mean;
      }

      gsl_vector_set(genotype, i, geno);
      gsl_matrix_set(UtX, i, c_snp, geno);
    }

    if (calc_K == true) {
      gsl_blas_dsyr(CblasUpper, 1.0, genotype, K);
    }

    c_snp++;
  }

  if (calc_K == true) {
    gsl_matrix_scale(K, 1.0 / (double)ns_test);

    for (size_t i = 0; i < genotype->size; ++i) {
      for (size_t j = 0; j < i; ++j) {
        geno = gsl_matrix_get(K, j, i);
        gsl_matrix_set(K, i, j, geno);
      }
    }
  }

  gsl_vector_free(genotype);
  gsl_vector_free(genotype_miss);

  infile.clear();
  infile.close();

  return true;
}

// Compact version of the above function, using uchar instead of
// gsl_matrix.
bool ReadFile_geno(const string &file_geno, vector<int> &indicator_idv,
                   vector<int> &indicator_snp,
                   vector<vector<unsigned char>> &Xt, gsl_matrix *K,
                   const bool calc_K, const size_t ni_test,
                   const size_t ns_test) {
  debug_msg("entered");
  igzstream infile(file_geno.c_str(), igzstream::in);
  if (!infile) {
    cout << "error reading genotype file:" << file_geno << endl;
    return false;
  }

  Xt.clear();
  vector<unsigned char> Xt_row;
  for (size_t i = 0; i < ni_test; i++) {
    Xt_row.push_back(0);
  }

  string line;
  char *ch_ptr;

  if (calc_K == true) {
    gsl_matrix_set_zero(K);
  }

  gsl_vector *genotype = gsl_vector_safe_alloc(ni_test);
  gsl_vector *genotype_miss = gsl_vector_safe_alloc(ni_test);
  double geno, geno_mean;
  size_t n_miss;

  size_t ni_total = indicator_idv.size();
  size_t ns_total = indicator_snp.size();

  size_t c_idv = 0, c_snp = 0;

  auto infilen = file_geno.c_str();
  for (size_t i = 0; i < ns_total; ++i) {
    safeGetline(infile, line).eof();
    if (indicator_snp[i] == 0) {
      continue;
    }

    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",infilen);
    ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
    ch_ptr = strtok_safe2(NULL, " ,\t",infilen);

    c_idv = 0;
    geno_mean = 0;
    n_miss = 0;
    gsl_vector_set_zero(genotype_miss);
    for (uint j = 0; j < ni_total; ++j) {
      ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
      if (indicator_idv[j] == 0) {
        continue;
      }

      if (strcmp(ch_ptr, "NA") == 0) {
        gsl_vector_set(genotype_miss, c_idv, 1);
        n_miss++;
      } else {
        geno = atof(ch_ptr);
        gsl_vector_set(genotype, c_idv, geno);
        geno_mean += geno;
      }
      c_idv++;
    }

    geno_mean /= (double)(ni_test - n_miss);

    for (size_t j = 0; j < genotype->size; ++j) {
      if (gsl_vector_get(genotype_miss, j) == 1) {
        geno = geno_mean;
      } else {
        geno = gsl_vector_get(genotype, j);
      }

      Xt_row[j] = Double02ToUchar(geno);
      gsl_vector_set(genotype, j, (geno - geno_mean));
    }
    Xt.push_back(Xt_row);

    if (calc_K == true) {
      gsl_blas_dsyr(CblasUpper, 1.0, genotype, K);
    }

    c_snp++;
  }

  if (calc_K == true) {
    gsl_matrix_scale(K, 1.0 / (double)ns_test);

    for (size_t i = 0; i < genotype->size; ++i) {
      for (size_t j = 0; j < i; ++j) {
        geno = gsl_matrix_get(K, j, i);
        gsl_matrix_set(K, i, j, geno);
      }
    }
  }

  gsl_vector_free(genotype);
  gsl_vector_free(genotype_miss);

  infile.clear();
  infile.close();

  return true;
}

// Read bimbam mean genotype file, the second time, recode "mean"
// genotype and calculate K.
bool ReadFile_bed(const string &file_bed, vector<int> &indicator_idv,
                  vector<int> &indicator_snp, gsl_matrix *UtX, gsl_matrix *K,
                  const bool calc_K) {
  debug_msg("entered");
  ifstream infile(file_bed.c_str(), ios::binary);
  if (!infile) {
    cout << "error reading bed file:" << file_bed << endl;
    return false;
  }

  char ch[1];
  bitset<8> b;

  size_t ni_total = indicator_idv.size();
  size_t ns_total = indicator_snp.size();
  size_t ni_test = UtX->size1;
  size_t ns_test = UtX->size2;
  int n_bit;

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

  if (calc_K == true) {
    gsl_matrix_set_zero(K);
  }

  gsl_vector *genotype = gsl_vector_safe_alloc(UtX->size1);

  double geno, geno_mean;
  size_t n_miss;
  size_t c_idv = 0, c_snp = 0, c = 0;

  // Start reading snps and doing association test.
  for (size_t t = 0; t < ns_total; ++t) {
    if (indicator_snp[t] == 0) {
      continue;
    }

    // n_bit, and 3 is the number of magic numbers.
    infile.seekg(t * n_bit + 3);

    // Read genotypes.
    c_idv = 0;
    geno_mean = 0.0;
    n_miss = 0;
    c = 0;
    for (int i = 0; i < n_bit; ++i) {
      infile.read(ch, 1);
      b = ch[0];

      // Minor allele homozygous: 2.0; major: 0.0.
      for (size_t j = 0; j < 4; ++j) {
        if ((i == (n_bit - 1)) && c == ni_total) {
          break;
        }
        if (indicator_idv[c] == 0) {
          c++;
          continue;
        }
        c++;

        if (b[2 * j] == 0) {
          if (b[2 * j + 1] == 0) {
            gsl_vector_set(genotype, c_idv, 2.0);
            geno_mean += 2.0;
          } else {
            gsl_vector_set(genotype, c_idv, 1.0);
            geno_mean += 1.0;
          }
        } else {
          if (b[2 * j + 1] == 1) {
            gsl_vector_set(genotype, c_idv, 0.0);
            geno_mean += 0.0;
          } else {
            gsl_vector_set(genotype, c_idv, -9.0);
            n_miss++;
          }
        }
        c_idv++;
      }
    }

    geno_mean /= (double)(ni_test - n_miss);

    for (size_t i = 0; i < genotype->size; ++i) {
      geno = gsl_vector_get(genotype, i);
      if (geno == -9) {
        geno = 0;
      } else {
        geno -= geno_mean;
      }

      gsl_vector_set(genotype, i, geno);
      gsl_matrix_set(UtX, i, c_snp, geno);
    }

    if (calc_K == true) {
      gsl_blas_dsyr(CblasUpper, 1.0, genotype, K);
    }

    c_snp++;
  }

  if (calc_K == true) {
    gsl_matrix_scale(K, 1.0 / (double)ns_test);

    for (size_t i = 0; i < genotype->size; ++i) {
      for (size_t j = 0; j < i; ++j) {
        geno = gsl_matrix_get(K, j, i);
        gsl_matrix_set(K, i, j, geno);
      }
    }
  }

  gsl_vector_free(genotype);
  infile.clear();
  infile.close();

  return true;
}

// Compact version of the above function, using uchar instead of gsl_matrix.
bool ReadFile_bed(const string &file_bed, vector<int> &indicator_idv,
                  vector<int> &indicator_snp, vector<vector<unsigned char>> &Xt,
                  gsl_matrix *K, const bool calc_K, const size_t ni_test,
                  const size_t ns_test) {
  debug_msg("entered");
  ifstream infile(file_bed.c_str(), ios::binary);
  if (!infile) {
    cout << "error reading bed file:" << file_bed << endl;
    return false;
  }

  Xt.clear();
  vector<unsigned char> Xt_row;
  for (size_t i = 0; i < ni_test; i++) {
    Xt_row.push_back(0);
  }

  char ch[1];
  bitset<8> b;

  size_t ni_total = indicator_idv.size();
  size_t ns_total = indicator_snp.size();
  int n_bit;

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

  if (calc_K == true) {
    gsl_matrix_set_zero(K);
  }

  gsl_vector *genotype = gsl_vector_safe_alloc(ni_test);

  double geno, geno_mean;
  size_t n_miss;
  size_t c_idv = 0, c_snp = 0, c = 0;

  // Start reading SNPs and doing association test.
  for (size_t t = 0; t < ns_total; ++t) {
    if (indicator_snp[t] == 0) {
      continue;
    }

    // n_bit, and 3 is the number of magic numbers.
    infile.seekg(t * n_bit + 3);

    // Read genotypes.
    c_idv = 0;
    geno_mean = 0.0;
    n_miss = 0;
    c = 0;
    for (int i = 0; i < n_bit; ++i) {
      infile.read(ch, 1);
      b = ch[0];

      // Minor allele homozygous: 2.0; major: 0.0.
      for (size_t j = 0; j < 4; ++j) {
        if ((i == (n_bit - 1)) && c == ni_total) {
          break;
        }
        if (indicator_idv[c] == 0) {
          c++;
          continue;
        }
        c++;

        if (b[2 * j] == 0) {
          if (b[2 * j + 1] == 0) {
            gsl_vector_set(genotype, c_idv, 2.0);
            geno_mean += 2.0;
          } else {
            gsl_vector_set(genotype, c_idv, 1.0);
            geno_mean += 1.0;
          }
        } else {
          if (b[2 * j + 1] == 1) {
            gsl_vector_set(genotype, c_idv, 0.0);
            geno_mean += 0.0;
          } else {
            gsl_vector_set(genotype, c_idv, -9.0);
            n_miss++;
          }
        }
        c_idv++;
      }
    }

    geno_mean /= (double)(ni_test - n_miss);

    for (size_t i = 0; i < genotype->size; ++i) {
      geno = gsl_vector_get(genotype, i);
      if (geno == -9) {
        geno = geno_mean;
      }

      Xt_row[i] = Double02ToUchar(geno);

      geno -= geno_mean;

      gsl_vector_set(genotype, i, geno);
    }
    Xt.push_back(Xt_row);

    if (calc_K == true) {
      gsl_blas_dsyr(CblasUpper, 1.0, genotype, K);
    }

    c_snp++;
  }

  if (calc_K == true) {
    gsl_matrix_scale(K, 1.0 / (double)ns_test);

    for (size_t i = 0; i < genotype->size; ++i) {
      for (size_t j = 0; j < i; ++j) {
        geno = gsl_matrix_get(K, j, i);
        gsl_matrix_set(K, i, j, geno);
      }
    }
  }

  gsl_vector_free(genotype);
  infile.clear();
  infile.close();

  return true;
}

bool ReadFile_est(const string &file_est, const vector<size_t> &est_column,
                  map<string, double> &mapRS2est) {
  debug_msg("entered");
  mapRS2est.clear();

  ifstream infile(file_est.c_str(), ifstream::in);
  if (!infile) {
    cout << "error opening estimated parameter file: " << file_est << endl;
    return false;
  }

  string line;
  char *ch_ptr;

  string rs;
  double alpha, beta, gamma, d;

  // Header.
  getline(infile, line);

  size_t n = *max_element(est_column.begin(), est_column.end());

  auto infilen = file_est.c_str();
  while (getline(infile, line)) {
    ch_ptr = strtok_safe2((char *)line.c_str(), " \t",infilen);

    alpha = 0.0;
    beta = 0.0;
    gamma = 1.0;
    for (size_t i = 0; i < n + 1; ++i) {
      if (i == est_column[0] - 1) {
        enforce(ch_ptr);
        rs = ch_ptr;
      }
      if (i == est_column[1] - 1) {
        enforce(ch_ptr);
        alpha = atof(ch_ptr);
      }
      if (i == est_column[2] - 1) {
        enforce(ch_ptr);
        beta = atof(ch_ptr);
      }
      if (i == est_column[3] - 1) {
        enforce(ch_ptr);
        gamma = atof(ch_ptr);
      }
      if (i < n) {
        ch_ptr = strtok(NULL, " \t");
      }
    }

    d = alpha + beta * gamma;

    if (mapRS2est.count(rs) == 0) {
      mapRS2est[rs] = d;
    } else {
      cout << "the same SNP occurs more than once in estimated "
           << "parameter file: " << rs << endl;
      return false;
    }
  }

  infile.clear();
  infile.close();
  return true;
}

bool CountFileLines(const string &file_input, size_t &n_lines) {
  debug_msg("entered");
  igzstream infile(file_input.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open file: " << file_input << endl;
    return false;
  }

  n_lines = count(istreambuf_iterator<char>(infile),
                  istreambuf_iterator<char>(), '\n');
  infile.seekg(0, ios::beg);

  return true;
}

// Read gene expression file.
bool ReadFile_gene(const string &file_gene, vector<double> &vec_read,
                   vector<SNPINFO> &snpInfo, size_t &ng_total) {
  debug_msg("entered");
  vec_read.clear();
  ng_total = 0;

  igzstream infile(file_gene.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open gene expression file: " << file_gene << endl;
    return false;
  }

  string line;
  char *ch_ptr;
  string rs;

  size_t n_idv = 0, t = 0;

  // Header.
  getline(infile, line);

  while (getline(infile, line)) {
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",file_gene.c_str());
    rs = ch_ptr;

    ch_ptr = strtok(NULL, " ,\t");

    t = 0;
    while (ch_ptr != NULL) {
      if (ng_total == 0) {
        vec_read.push_back(0);
        t++;
        n_idv++;
      } else {
        vec_read[t] += atof(ch_ptr);
        t++;
      }

      ch_ptr = strtok(NULL, " ,\t");
    }

    if (t != n_idv) {
      cout << "error! number of columns doesn't match in row: " << ng_total
           << endl;
      return false;
    }

    SNPINFO sInfo = {"-9", rs, -9, -9, "-9", "-9", 0, -9, -9, 0, 0, 0};
    snpInfo.push_back(sInfo);

    ng_total++;
  }

  infile.close();
  infile.clear();

  return true;
}

// Read header to determine which column contains which item.
bool ReadHeader_io(const string &line, HEADER &header) {
  debug_msg("entered");
  string rs_ptr[] = {"rs",    "RS",    "snp",  "SNP",  "snps",      "SNPS",
                     "snpid", "SNPID", "rsid", "RSID", "MarkerName"};
  set<string> rs_set(rs_ptr, rs_ptr + 11); // create a set of 11 items
  string chr_ptr[] = {"chr", "CHR"};
  set<string> chr_set(chr_ptr, chr_ptr + 2);
  string pos_ptr[] = {
      "ps", "PS", "pos", "POS", "base_position", "BASE_POSITION", "bp", "BP"};
  set<string> pos_set(pos_ptr, pos_ptr + 8);
  string cm_ptr[] = {"cm", "CM"};
  set<string> cm_set(cm_ptr, cm_ptr + 2);
  string a1_ptr[] = {"a1", "A1", "allele1", "ALLELE1", "Allele1", "INC_ALLELE"};
  set<string> a1_set(a1_ptr, a1_ptr + 5);
  string a0_ptr[] = {"a0", "A0",      "allele0", "ALLELE0", "Allele0",   "a2",
                     "A2", "allele2", "ALLELE2", "Allele2", "DEC_ALLELE"};
  set<string> a0_set(a0_ptr, a0_ptr + 10);

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
  string ncase_ptr[] = {"ncase", "NCASE", "n_case", "N_CASE"};
  set<string> ncase_set(ncase_ptr, ncase_ptr + 4);
  string ncontrol_ptr[] = {"ncontrol", "NCONTROL", "n_control", "N_CONTROL"};
  set<string> ncontrol_set(ncontrol_ptr, ncontrol_ptr + 4);

  string af_ptr[] = {"af",
                     "AF",
                     "maf",
                     "MAF",
                     "f",
                     "F",
                     "allele_freq",
                     "ALLELE_FREQ",
                     "allele_frequency",
                     "ALLELE_FREQUENCY",
                     "Freq.Allele1.HapMapCEU",
                     "FreqAllele1HapMapCEU",
                     "Freq1.Hapmap"};
  set<string> af_set(af_ptr, af_ptr + 13);
  string var_ptr[] = {"var", "VAR"};
  set<string> var_set(var_ptr, var_ptr + 2);

  string ws_ptr[] = {"window_size", "WINDOW_SIZE", "ws", "WS"};
  set<string> ws_set(ws_ptr, ws_ptr + 4);
  string cor_ptr[] = {"cor", "COR", "r", "R"};
  set<string> cor_set(cor_ptr, cor_ptr + 4);

  header.rs_col = 0;
  header.chr_col = 0;
  header.pos_col = 0;
  header.cm_col = 0;
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
  header.ncase_col = 0;
  header.ncontrol_col = 0;
  header.af_col = 0;
  header.var_col = 0;
  header.ws_col = 0;
  header.cor_col = 0;
  header.coln = 0;

  char *ch_ptr;
  string type;
  size_t n_error = 0;

  ch_ptr = strtok((char *)line.c_str(), " ,\t");
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
    } else if (ncase_set.count(type) != 0) {
      if (header.ncase_col == 0) {
        header.ncase_col = header.coln + 1;
      } else {
        cout << "error! more than two n_case columns in the file." << endl;
        n_error++;
      }
    } else if (ncontrol_set.count(type) != 0) {
      if (header.ncontrol_col == 0) {
        header.ncontrol_col = header.coln + 1;
      } else {
        cout << "error! more than two n_control columns in the file." << endl;
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
      string str = ch_ptr;
      string cat = str.substr(str.size() - 2, 2);

      if (cat == "_c" || cat == "_C") {

        // continuous
        header.catc_col.insert(header.coln + 1);
      } else {

        // discrete
        header.catd_col.insert(header.coln + 1);
      }
    }

    ch_ptr = strtok(NULL, " ,\t");
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

// Read category file, record mapRS2 in the category file does not
// contain a null category so if a snp has 0 for all categories, then
// it is not included in the analysis.
bool ReadFile_cat(const string &file_cat, map<string, size_t> &mapRS2cat,
                  size_t &n_vc) {
  debug_msg("entered");
  mapRS2cat.clear();

  igzstream infile(file_cat.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open category file: " << file_cat << endl;
    return false;
  }

  string line;
  char *ch_ptr;

  string rs, chr, a1, a0, pos, cm;
  size_t i_cat;

  // Read header.
  HEADER header;
  safeGetline(infile, line).eof();
  ReadHeader_io(line, header);

  // Use the header to count the number of categories.
  n_vc = header.coln;
  if (header.rs_col != 0) {
    n_vc--;
  }
  if (header.chr_col != 0) {
    n_vc--;
  }
  if (header.pos_col != 0) {
    n_vc--;
  }
  if (header.cm_col != 0) {
    n_vc--;
  }
  if (header.a1_col != 0) {
    n_vc--;
  }
  if (header.a0_col != 0) {
    n_vc--;
  }

  // Read the following lines to record mapRS2cat.
  while (!safeGetline(infile, line).eof()) {
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",file_cat.c_str());

    i_cat = 0;
    for (size_t i = 0; i < header.coln; i++) {
      enforce(ch_ptr);
      if (header.rs_col != 0 && header.rs_col == i + 1) {
        rs = ch_ptr;
      } else if (header.chr_col != 0 && header.chr_col == i + 1) {
        chr = ch_ptr;
      } else if (header.pos_col != 0 && header.pos_col == i + 1) {
        pos = ch_ptr;
      } else if (header.cm_col != 0 && header.cm_col == i + 1) {
        cm = ch_ptr;
      } else if (header.a1_col != 0 && header.a1_col == i + 1) {
        a1 = ch_ptr;
      } else if (header.a0_col != 0 && header.a0_col == i + 1) {
        a0 = ch_ptr;
      } else if (atoi(ch_ptr) == 1 || atoi(ch_ptr) == 0) {
        if (i_cat == 0) {
          if (header.rs_col == 0) {
            rs = chr + ":" + pos;
          }
        }

        if (atoi(ch_ptr) == 1 && mapRS2cat.count(rs) == 0) {
          mapRS2cat[rs] = i_cat;
        }
        i_cat++;
      } else {
      }

      ch_ptr = strtok(NULL, " ,\t");
    }
  }

  infile.clear();
  infile.close();

  return true;
}

bool ReadFile_mcat(const string &file_mcat, map<string, size_t> &mapRS2cat,
                   size_t &n_vc) {
  debug_msg("entered");
  mapRS2cat.clear();

  igzstream infile(file_mcat.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open mcategory file: " << file_mcat << endl;
    return false;
  }

  string file_name;
  map<string, size_t> mapRS2cat_tmp;
  size_t n_vc_tmp, t = 0;

  while (!safeGetline(infile, file_name).eof()) {
    mapRS2cat_tmp.clear();
    ReadFile_cat(file_name, mapRS2cat_tmp, n_vc_tmp);
    mapRS2cat.insert(mapRS2cat_tmp.begin(), mapRS2cat_tmp.end());
    if (t == 0) {
      n_vc = n_vc_tmp;
    } else {
      n_vc = max(n_vc, n_vc_tmp);
    }
    t++;
  }

  return true;
}

// Read bimbam mean genotype file and calculate kinship matrix; this
// time, the kinship matrix is not centered, and can contain multiple
// K matrix.
bool BimbamKinUncentered(const string &file_geno, const set<string> ksnps,
                         const int display_pace,
                         const vector<int> &indicator_idv,
                         const vector<int> &indicator_snp,
                         const map<string, double> &mapRS2weight,
                         const map<string, size_t> &mapRS2cat,
                         const vector<SNPINFO> &snpInfo, const gsl_matrix *W,
                         gsl_matrix *matrix_kin, gsl_vector *vector_ns) {
  debug_msg("entered");
  debug_msg("BimbamKinUncentered");
  igzstream infile(file_geno.c_str(), igzstream::in);
  if (!infile) {
    cout << "error reading genotype file:" << file_geno << endl;
    return false;
  }

  string line;
  char *ch_ptr;

  size_t n_miss;
  double d, geno_mean, geno_var;

  size_t ni_test = matrix_kin->size1;
  gsl_vector *geno = gsl_vector_safe_alloc(ni_test);
  gsl_vector *geno_miss = gsl_vector_safe_alloc(ni_test);

  gsl_vector *Wtx = gsl_vector_safe_alloc(W->size2);
  gsl_matrix *WtW = gsl_matrix_safe_alloc(W->size2, W->size2);
  gsl_matrix *WtWi = gsl_matrix_safe_alloc(W->size2, W->size2);
  gsl_vector *WtWiWtx = gsl_vector_safe_alloc(W->size2);
  gsl_permutation *pmt = gsl_permutation_alloc(W->size2);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
  int sig;
  LUDecomp(WtW, pmt, &sig);
  LUInvert(WtW, pmt, WtWi);

  size_t n_vc = matrix_kin->size2 / ni_test, i_vc;
  string rs;
  vector<size_t> ns_vec;
  for (size_t i = 0; i < n_vc; i++) {
    ns_vec.push_back(0);
  }

  // Create a large matrix.
  const size_t msize = K_BATCH_SIZE;
  gsl_matrix *Xlarge = gsl_matrix_safe_alloc(ni_test, msize * n_vc);
  gsl_matrix_set_zero(Xlarge);

  size_t ns_test = 0;
  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    safeGetline(infile, line).eof();
    if (t % display_pace == 0 || t == (indicator_snp.size() - 1)) {
      ProgressBar("Reading SNPs", t, indicator_snp.size() - 1);
    }
    if (indicator_snp[t] == 0)
      continue;

    auto infilen = file_geno.c_str();
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",infilen);
    ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
    ch_ptr = strtok_safe2(NULL, " ,\t",infilen);

    rs = snpInfo[t].rs_number; // This line is new.

    geno_mean = 0.0;
    n_miss = 0;
    geno_var = 0.0;
    gsl_vector_set_all(geno_miss, 0);

    size_t j = 0;
    for (size_t i = 0; i < indicator_idv.size(); ++i) {
      if (indicator_idv[i] == 0) {
        continue;
      }
      ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
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

    gsl_blas_dgemv(CblasTrans, 1.0, W, geno, 0.0, Wtx);
    gsl_blas_dgemv(CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
    gsl_blas_dgemv(CblasNoTrans, -1.0, W, WtWiWtx, 1.0, geno);
    gsl_blas_ddot(geno, geno, &geno_var);
    geno_var /= (double)ni_test;

    if (geno_var != 0 &&
        (mapRS2weight.size() == 0 || mapRS2weight.count(rs) != 0)) {
      if (mapRS2weight.size() == 0) {
        d = 1.0 / geno_var;
      } else {
        d = mapRS2weight.at(rs) / geno_var;
      }

      gsl_vector_scale(geno, sqrt(d));
      if (n_vc == 1 || mapRS2cat.size() == 0) {
        gsl_vector_view Xlarge_col =
            gsl_matrix_column(Xlarge, ns_vec[0] % msize);
        gsl_vector_memcpy(&Xlarge_col.vector, geno);
        ns_vec[0]++;

        if (ns_vec[0] % msize == 0) {
          fast_eigen_dgemm("N", "T", 1.0, Xlarge, Xlarge, 1.0, matrix_kin);
          gsl_matrix_set_zero(Xlarge);
        }
      } else if (mapRS2cat.count(rs) != 0) {
        i_vc = mapRS2cat.at(rs);

        gsl_vector_view Xlarge_col =
            gsl_matrix_column(Xlarge, msize * i_vc + ns_vec[i_vc] % msize);
        gsl_vector_memcpy(&Xlarge_col.vector, geno);

        ns_vec[i_vc]++;

        if (ns_vec[i_vc] % msize == 0) {
          gsl_matrix_view X_sub =
              gsl_matrix_submatrix(Xlarge, 0, msize * i_vc, ni_test, msize);
          gsl_matrix_view kin_sub = gsl_matrix_submatrix(
              matrix_kin, 0, ni_test * i_vc, ni_test, ni_test);
          fast_eigen_dgemm("N", "T", 1.0, &X_sub.matrix, &X_sub.matrix, 1.0,
                         &kin_sub.matrix);

          gsl_matrix_set_zero(&X_sub.matrix);
        }
      }
    }
    ns_test++;
  }

  for (size_t i_vc = 0; i_vc < n_vc; i_vc++) {
    if (ns_vec[i_vc] % msize != 0) {
      gsl_matrix_view X_sub =
          gsl_matrix_submatrix(Xlarge, 0, msize * i_vc, ni_test, msize);
      gsl_matrix_view kin_sub =
          gsl_matrix_submatrix(matrix_kin, 0, ni_test * i_vc, ni_test, ni_test);
      fast_eigen_dgemm("N", "T", 1.0, &X_sub.matrix, &X_sub.matrix, 1.0,
                     &kin_sub.matrix);
    }
  }

  cout << endl;

  for (size_t t = 0; t < n_vc; t++) {
    gsl_vector_set(vector_ns, t, ns_vec[t]);

    for (size_t i = 0; i < ni_test; ++i) {
      for (size_t j = 0; j <= i; ++j) {
        d = gsl_matrix_get(matrix_kin, j, i + ni_test * t);
        d /= (double)ns_vec[t];
        gsl_matrix_set(matrix_kin, i, j + ni_test * t, d);
        gsl_matrix_set(matrix_kin, j, i + ni_test * t, d);
      }
    }
  }

  gsl_vector_free(geno);
  gsl_vector_free(geno_miss);

  gsl_vector_free(Wtx);
  gsl_matrix_free(WtW);
  gsl_matrix_free(WtWi);
  gsl_vector_free(WtWiWtx);
  gsl_permutation_free(pmt);

  gsl_matrix_free(Xlarge);

  infile.close();
  infile.clear();

  return true;
}

bool PlinkKin(const string &file_bed, const int display_pace,
              const vector<int> &indicator_idv,
              const vector<int> &indicator_snp,
              const map<string, double> &mapRS2weight,
              const map<string, size_t> &mapRS2cat,
              const vector<SNPINFO> &snpInfo, const gsl_matrix *W,
              gsl_matrix *matrix_kin, gsl_vector *vector_ns) {
  debug_msg("entered");
  ifstream infile(file_bed.c_str(), ios::binary);
  if (!infile) {
    cout << "error reading bed file:" << file_bed << endl;
    return false;
  }

  char ch[1];
  bitset<8> b;

  size_t n_miss, ci_total, ci_test;
  double d, geno_mean, geno_var;

  size_t ni_test = matrix_kin->size1;
  size_t ni_total = indicator_idv.size();
  gsl_vector *geno = gsl_vector_safe_alloc(ni_test);

  gsl_vector *Wtx = gsl_vector_safe_alloc(W->size2);
  gsl_matrix *WtW = gsl_matrix_safe_alloc(W->size2, W->size2);
  gsl_matrix *WtWi = gsl_matrix_safe_alloc(W->size2, W->size2);
  gsl_vector *WtWiWtx = gsl_vector_safe_alloc(W->size2);
  gsl_permutation *pmt = gsl_permutation_alloc(W->size2);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
  int sig;
  LUDecomp(WtW, pmt, &sig);
  LUInvert(WtW, pmt, WtWi);

  size_t ns_test = 0;
  int n_bit;

  size_t n_vc = matrix_kin->size2 / ni_test, i_vc;
  string rs;
  vector<size_t> ns_vec;
  for (size_t i = 0; i < n_vc; i++) {
    ns_vec.push_back(0);
  }

  // Create a large matrix.
  const size_t msize = K_BATCH_SIZE;
  gsl_matrix *Xlarge = gsl_matrix_safe_alloc(ni_test, msize * n_vc);
  gsl_matrix_set_zero(Xlarge);

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

  for (size_t t = 0; t < indicator_snp.size(); ++t) {
    if (t % display_pace == 0 || t == (indicator_snp.size() - 1)) {
      ProgressBar("Reading SNPs", t, indicator_snp.size() - 1);
    }
    if (indicator_snp[t] == 0) {
      continue;
    }

    // n_bit, and 3 is the number of magic numbers
    infile.seekg(t * n_bit + 3);

    rs = snpInfo[t].rs_number; // This line is new.

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

    gsl_blas_dgemv(CblasTrans, 1.0, W, geno, 0.0, Wtx);
    gsl_blas_dgemv(CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
    gsl_blas_dgemv(CblasNoTrans, -1.0, W, WtWiWtx, 1.0, geno);
    gsl_blas_ddot(geno, geno, &geno_var);
    geno_var /= (double)ni_test;

    if (geno_var != 0 &&
        (mapRS2weight.size() == 0 || mapRS2weight.count(rs) != 0)) {
      if (mapRS2weight.size() == 0) {
        d = 1.0 / geno_var;
      } else {
        d = mapRS2weight.at(rs) / geno_var;
      }

      gsl_vector_scale(geno, sqrt(d));
      if (n_vc == 1 || mapRS2cat.size() == 0) {
        gsl_vector_view Xlarge_col =
            gsl_matrix_column(Xlarge, ns_vec[0] % msize);
        gsl_vector_memcpy(&Xlarge_col.vector, geno);
        ns_vec[0]++;

        if (ns_vec[0] % msize == 0) {
          fast_eigen_dgemm("N", "T", 1.0, Xlarge, Xlarge, 1.0, matrix_kin);
          gsl_matrix_set_zero(Xlarge);
        }
      } else if (mapRS2cat.count(rs) != 0) {
        i_vc = mapRS2cat.at(rs);

        gsl_vector_view Xlarge_col =
            gsl_matrix_column(Xlarge, msize * i_vc + ns_vec[i_vc] % msize);
        gsl_vector_memcpy(&Xlarge_col.vector, geno);

        ns_vec[i_vc]++;

        if (ns_vec[i_vc] % msize == 0) {
          gsl_matrix_view X_sub =
              gsl_matrix_submatrix(Xlarge, 0, msize * i_vc, ni_test, msize);
          gsl_matrix_view kin_sub = gsl_matrix_submatrix(
              matrix_kin, 0, ni_test * i_vc, ni_test, ni_test);
          fast_eigen_dgemm("N", "T", 1.0, &X_sub.matrix, &X_sub.matrix, 1.0,
                         &kin_sub.matrix);

          gsl_matrix_set_zero(&X_sub.matrix);
        }
      }
    }
    ns_test++;
  }

  for (size_t i_vc = 0; i_vc < n_vc; i_vc++) {
    if (ns_vec[i_vc] % msize != 0) {
      gsl_matrix_view X_sub =
          gsl_matrix_submatrix(Xlarge, 0, msize * i_vc, ni_test, msize);
      gsl_matrix_view kin_sub =
          gsl_matrix_submatrix(matrix_kin, 0, ni_test * i_vc, ni_test, ni_test);
      fast_eigen_dgemm("N", "T", 1.0, &X_sub.matrix, &X_sub.matrix, 1.0,
                     &kin_sub.matrix);
    }
  }

  cout << endl;

  for (size_t t = 0; t < n_vc; t++) {
    gsl_vector_set(vector_ns, t, ns_vec[t]);

    for (size_t i = 0; i < ni_test; ++i) {
      for (size_t j = 0; j <= i; ++j) {
        d = gsl_matrix_get(matrix_kin, j, i + ni_test * t);
        d /= (double)ns_vec[t];
        gsl_matrix_set(matrix_kin, i, j + ni_test * t, d);
        gsl_matrix_set(matrix_kin, j, i + ni_test * t, d);
      }
    }
  }

  gsl_vector_free(geno);

  gsl_vector_free(Wtx);
  gsl_matrix_free(WtW);
  gsl_matrix_free(WtWi);
  gsl_vector_free(WtWiWtx);
  gsl_permutation_free(pmt);

  gsl_matrix_free(Xlarge);

  infile.close();
  infile.clear();

  return true;
}

bool MFILEKin(const size_t mfile_mode, const string &file_mfile,
              const set<string> setKSnps, const int display_pace,
              const vector<int> &indicator_idv,
              const vector<vector<int>> &mindicator_snp,
              const map<string, double> &mapRS2weight,
              const map<string, size_t> &mapRS2cat,
              const vector<vector<SNPINFO>> &msnpInfo, const gsl_matrix *W,
              gsl_matrix *matrix_kin, gsl_vector *vector_ns) {
  debug_msg("entered");
  size_t n_vc = vector_ns->size, ni_test = matrix_kin->size1;
  gsl_matrix_set_zero(matrix_kin);
  gsl_vector_set_zero(vector_ns);

  igzstream infile(file_mfile.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open mfile file: " << file_mfile << endl;
    return false;
  }

  string file_name;

  gsl_matrix *kin_tmp = gsl_matrix_safe_alloc(matrix_kin->size1, matrix_kin->size2);
  gsl_vector *ns_tmp = gsl_vector_safe_alloc(vector_ns->size);

  size_t l = 0;
  double d;
  while (!safeGetline(infile, file_name).eof()) {
    gsl_matrix_set_zero(kin_tmp);
    gsl_vector_set_zero(ns_tmp);

    if (mfile_mode == 1) {
      file_name += ".bed";
      PlinkKin(file_name, display_pace, indicator_idv, mindicator_snp[l],
               mapRS2weight, mapRS2cat, msnpInfo[l], W, kin_tmp, ns_tmp);
    } else {
      BimbamKinUncentered(file_name, setKSnps, display_pace, indicator_idv,
                          mindicator_snp[l], mapRS2weight, mapRS2cat,
                          msnpInfo[l], W, kin_tmp, ns_tmp);
    }

    // Add ns.
    gsl_vector_add(vector_ns, ns_tmp);

    // Add kin.
    for (size_t t = 0; t < n_vc; t++) {
      for (size_t i = 0; i < ni_test; ++i) {
        for (size_t j = 0; j <= i; ++j) {
          d = gsl_matrix_get(matrix_kin, j, i + ni_test * t) +
              gsl_matrix_get(kin_tmp, j, i + ni_test * t) *
                  gsl_vector_get(ns_tmp, t);

          gsl_matrix_set(matrix_kin, i, j + ni_test * t, d);
          gsl_matrix_set(matrix_kin, j, i + ni_test * t, d);
        }
      }
    }
    l++;
  }

  // Renormalize kin.
  for (size_t t = 0; t < n_vc; t++) {
    for (size_t i = 0; i < ni_test; ++i) {
      for (size_t j = 0; j <= i; ++j) {
        d = gsl_matrix_get(matrix_kin, j, i + ni_test * t) /
            gsl_vector_get(vector_ns, t);

        gsl_matrix_set(matrix_kin, i, j + ni_test * t, d);
        gsl_matrix_set(matrix_kin, j, i + ni_test * t, d);
      }
    }
  }
  cout << endl;

  infile.close();
  infile.clear();

  gsl_matrix_free(kin_tmp);
  gsl_vector_free(ns_tmp);

  return true;
}

// Read var file, store mapRS2wsnp.
bool ReadFile_wsnp(const string &file_wsnp, map<string, double> &mapRS2weight) {
  debug_msg("entered");
  mapRS2weight.clear();

  igzstream infile(file_wsnp.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open snp weight file: " << file_wsnp << endl;
    return false;
  }

  char *ch_ptr;
  string line, rs;
  double weight;

  auto infilen = file_wsnp.c_str();
  while (!safeGetline(infile, line).eof()) {
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",infilen);
    rs = ch_ptr;
    ch_ptr = strtok_safe2(NULL, " ,\t",infilen);
    weight = atof(ch_ptr);
    mapRS2weight[rs] = weight;
  }

  return true;
}

bool ReadFile_wsnp(const string &file_wcat, const size_t n_vc,
                   map<string, vector<double>> &mapRS2wvector) {
  debug_msg("entered");
  mapRS2wvector.clear();

  igzstream infile(file_wcat.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open snp weight file: " << file_wcat << endl;
    return false;
  }

  char *ch_ptr;
  vector<double> weight;
  for (size_t i = 0; i < n_vc; i++) {
    weight.push_back(0.0);
  }

  string line, rs, chr, a1, a0, pos, cm;

  // Read header.
  HEADER header;
  safeGetline(infile, line).eof();
  ReadHeader_io(line, header);

  while (!safeGetline(infile, line).eof()) {
    if (isBlankLine(line)) {
      continue;
    }
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",file_wcat.c_str());

    size_t t = 0;
    for (size_t i = 0; i < header.coln; i++) {
      enforce(ch_ptr);
      if (header.rs_col != 0 && header.rs_col == i + 1) {
        rs = ch_ptr;
      } else if (header.chr_col != 0 && header.chr_col == i + 1) {
        chr = ch_ptr;
      } else if (header.pos_col != 0 && header.pos_col == i + 1) {
        pos = ch_ptr;
      } else if (header.cm_col != 0 && header.cm_col == i + 1) {
        cm = ch_ptr;
      } else if (header.a1_col != 0 && header.a1_col == i + 1) {
        a1 = ch_ptr;
      } else if (header.a0_col != 0 && header.a0_col == i + 1) {
        a0 = ch_ptr;
      } else {
        weight[t] = atof(ch_ptr);
        t++;
        if (t > n_vc) {
          cout << "error! Number of columns in the wcat file does not "
               << "match that of cat file.";
          return false;
        }
      }

      ch_ptr = strtok(NULL, " ,\t");
    }

    if (t != n_vc) {
      cout << "error! Number of columns in the wcat file does not "
           << "match that of cat file.";
      return false;
    }

    if (header.rs_col == 0) {
      rs = chr + ":" + pos;
    }

    mapRS2wvector[rs] = weight;
  }

  return true;
}

// Read the beta file, save snp z scores in to z2_score, and save
// category into indicator_snp based on mapRS2var and set, and
// indicator_snp record the category number (from 1 to n_vc), and
// provide var if maf/var is not provided in the beta file notice that
// indicator_snp contains ns_test snps, instead of ns_total snps read
// the beta file for the second time, compute q, and Vq based on block
// jacknife use the mapRS2var to select snps (and to ), calculate q do
// a block-wise jacknife, and compute Vq
void ReadFile_beta(const string &file_beta,
                   const map<string, size_t> &mapRS2cat,
                   const map<string, double> &mapRS2wA, vector<size_t> &vec_cat,
                   vector<size_t> &vec_ni, vector<double> &vec_weight,
                   vector<double> &vec_z2, size_t &ni_total, size_t &ns_total,
                   size_t &ns_test) {
  debug_msg("entered");
  vec_cat.clear();
  vec_ni.clear();
  vec_weight.clear();
  vec_z2.clear();
  ni_total = 0;
  ns_total = 0;
  ns_test = 0;

  igzstream infile(file_beta.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open beta file: " << file_beta << endl;
    return;
  }

  string line;
  char *ch_ptr;
  string type;

  string rs, chr, a1, a0, pos, cm;
  double z = 0, beta = 0, se_beta = 0, pvalue = 0, zsquare = 0; // af = 0;
  size_t n_total = 0, n_mis = 0, n_obs = 0, n_case = 0, n_control = 0;

  // Read header.
  HEADER header;
  safeGetline(infile, line).eof();
  ReadHeader_io(line, header);

  if (header.n_col == 0) {
    if ((header.nobs_col == 0 && header.nmis_col == 0) &&
        (header.ncase_col == 0 && header.ncontrol_col == 0)) {
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

  while (!safeGetline(infile, line).eof()) {
    if (isBlankLine(line)) {
      continue;
    }
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",file_beta.c_str());

    z = 0;
    beta = 0;
    se_beta = 0;
    auto chisq = 0.0;
    pvalue = 0;
    n_total = 0;
    n_mis = 0;
    n_obs = 0;
    n_case = 0;
    n_control = 0;
    // af = 0;
    // auto var_x = 0.0;
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
      if (header.ncase_col != 0 && header.ncase_col == i + 1) {
        n_case = atoi(ch_ptr);
      }
      if (header.ncontrol_col != 0 && header.ncontrol_col == i + 1) {
        n_control = atoi(ch_ptr);
      }
      // if (header.af_col != 0 && header.af_col == i + 1) {
      //   af = atof(ch_ptr);
      // }
      // if (header.var_col != 0 && header.var_col == i + 1) {
      //   var_x = atof(ch_ptr);
      // }

      ch_ptr = strtok(NULL, " ,\t");
    }

    if (header.rs_col == 0) {
      rs = chr + ":" + pos;
    }

    if (header.n_col == 0) {
      if (header.nmis_col != 0 && header.nobs_col != 0) {
        n_total = n_mis + n_obs;
      } else {
        n_total = n_case + n_control;
      }
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

    // Obtain var_x.
    // if (header.var_col == 0 && header.af_col != 0) {
    //   var_x = 2.0 * af * (1.0 - af);
    // }

    // If the SNP is also present in cor file, then do calculations.
    if ((mapRS2wA.size() == 0 || mapRS2wA.count(rs) != 0) &&
        (mapRS2cat.size() == 0 || mapRS2cat.count(rs) != 0) && zsquare != 0) {
      if (mapRS2cat.size() != 0) {
        vec_cat.push_back(mapRS2cat.at(rs));
      } else {
        vec_cat.push_back(0);
      }
      vec_ni.push_back(n_total);
      if (mapRS2wA.size() == 0) {
        vec_weight.push_back(1);
      } else {
        vec_weight.push_back(mapRS2wA.at(rs));
      }
      vec_z2.push_back(zsquare);

      ni_total = max(ni_total, n_total);
      ns_test++;
    }

    ns_total++;
  }

  infile.clear();
  infile.close();

  return;
}

void ReadFile_beta(const string &file_beta, const map<string, double> &mapRS2wA,
                   map<string, string> &mapRS2A1,
                   map<string, double> &mapRS2z) {
  debug_msg("entered");
  mapRS2A1.clear();
  mapRS2z.clear();

  igzstream infile(file_beta.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open beta file: " << file_beta << endl;
    return;
  }

  string line;
  char *ch_ptr;
  string type;

  string rs, chr, a1, a0, pos, cm;
  double z = 0, beta = 0, se_beta = 0; // pvalue = 0, chisq=0, af = 0 , var_x = 0;
  size_t n_total = 0, n_mis = 0, n_obs = 0, n_case = 0, n_control = 0;
  size_t ni_total = 0, ns_total = 0, ns_test = 0;

  // Read header.
  HEADER header;
  safeGetline(infile, line).eof();
  ReadHeader_io(line, header);

  if (header.n_col == 0) {
    if ((header.nobs_col == 0 && header.nmis_col == 0) &&
        (header.ncase_col == 0 && header.ncontrol_col == 0)) {
      cout << "error! missing sample size in the beta file." << endl;
    } else {
      cout << "total sample size will be replaced by obs/mis sample size."
           << endl;
    }
  }

  if (header.z_col == 0 && (header.beta_col == 0 || header.sebeta_col == 0)) {
    cout << "error! missing z scores in the beta file." << endl;
  }

  while (!safeGetline(infile, line).eof()) {
    if (isBlankLine(line)) {
      continue;
    }
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",file_beta.c_str());

    z = 0;
    beta = 0;
    se_beta = 0;
    // chisq = 0;
    // pvalue = 0;
    n_total = 0;
    n_mis = 0;
    n_obs = 0;
    n_case = 0;
    n_control = 0;
    // af = 0;
    // double var_x = 0;
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
      // if (header.chisq_col != 0 && header.chisq_col == i + 1) {
      //   chisq = atof(ch_ptr);
      // }
      // if (header.p_col != 0 && header.p_col == i + 1) {
      //   pvalue = atof(ch_ptr);
      // }

      if (header.n_col != 0 && header.n_col == i + 1) {
        n_total = atoi(ch_ptr);
      }
      if (header.nmis_col != 0 && header.nmis_col == i + 1) {
        n_mis = atoi(ch_ptr);
      }
      if (header.nobs_col != 0 && header.nobs_col == i + 1) {
        n_obs = atoi(ch_ptr);
      }
      if (header.ncase_col != 0 && header.ncase_col == i + 1) {
        n_case = atoi(ch_ptr);
      }
      if (header.ncontrol_col != 0 && header.ncontrol_col == i + 1) {
        n_control = atoi(ch_ptr);
      }

      // if (header.af_col != 0 && header.af_col == i + 1) {
      //   af = atof(ch_ptr);
      // }

      // if (header.var_col != 0 && header.var_col == i + 1) {
      //   var_x = atof(ch_ptr);
      // }

      ch_ptr = strtok(NULL, " ,\t");
    }

    if (header.rs_col == 0) {
      rs = chr + ":" + pos;
    }

    if (header.n_col == 0) {
      if (header.nmis_col != 0 && header.nobs_col != 0) {
        n_total = n_mis + n_obs;
      } else {
        n_total = n_case + n_control;
      }
    }

    // Both z values and beta/se_beta have directions, while
    // chisq/pvalue do not.
    if (header.z_col != 0) {
      z = z;
    } else if (header.beta_col != 0 && header.sebeta_col != 0) {
      z = beta / se_beta;
    } else {
      z = 0;
    }

    // If the snp is also present in cor file, then do calculations.
    if ((mapRS2wA.size() == 0 || mapRS2wA.count(rs) != 0)) {
      mapRS2z[rs] = z;
      mapRS2A1[rs] = a1;

      ni_total = max(ni_total, n_total);
      ns_test++;
    }

    ns_total++;
  }

  infile.clear();
  infile.close();

  return;
}

void Calcq(const size_t n_block, const vector<size_t> &vec_cat,
           const vector<size_t> &vec_ni, const vector<double> &vec_weight,
           const vector<double> &vec_z2, gsl_matrix *Vq, gsl_vector *q,
           gsl_vector *s) {
  debug_msg("entered");
  gsl_matrix_set_zero(Vq);
  gsl_vector_set_zero(q);
  gsl_vector_set_zero(s);

  size_t cat, n_total;
  double w, zsquare;

  vector<double> vec_q, vec_s, n_snps;
  for (size_t i = 0; i < q->size; i++) {
    vec_q.push_back(0.0);
    vec_s.push_back(0.0);
    n_snps.push_back(0.0);
  }

  vector<vector<double>> mat_q, mat_s;
  for (size_t i = 0; i < n_block; i++) {
    mat_q.push_back(vec_q);
    mat_s.push_back(vec_s);
  }

  // Compute q and s.
  for (size_t i = 0; i < vec_cat.size(); i++) {

    // Extract quantities.
    cat = vec_cat[i];
    n_total = vec_ni[i];
    w = vec_weight[i];
    zsquare = vec_z2[i];

    // Compute q and s.
    vec_q[cat] += (zsquare - 1.0) * w / (double)n_total;
    vec_s[cat] += w;
    n_snps[cat]++;
  }

  // Update q; vec_q is used again for computing Vq below.
  for (size_t i = 0; i < q->size; i++) {
    if (vec_s[i] != 0) {
      gsl_vector_set(q, i, vec_q[i] / vec_s[i]);
    }
    gsl_vector_set(s, i, vec_s[i]);
  }

  // Compute Vq; divide SNPs in each category into evenly distributed
  // blocks.
  size_t t = 0, b = 0, n_snp = 0;
  double d, m, n;
  for (size_t l = 0; l < q->size; l++) {
    n_snp = floor(n_snps[l] / n_block);
    t = 0;
    b = 0;
    if (n_snp == 0) {
      continue;
    }

    // Initiate everything to zero.
    for (size_t i = 0; i < n_block; i++) {
      for (size_t j = 0; j < q->size; j++) {
        mat_q[i][j] = 0;
        mat_s[i][j] = 0;
      }
    }

    // Record values.
    for (size_t i = 0; i < vec_cat.size(); i++) {

      // Extract quantities.
      cat = vec_cat[i];
      n_total = vec_ni[i];
      w = vec_weight[i];
      zsquare = vec_z2[i];

      // Save quantities for computing Vq (which is not divided by
      // n_total).
      mat_q[b][cat] += (zsquare - 1.0) * w;
      mat_s[b][cat] += w;

      if (cat == l) {
        if (b < n_block - 1) {
          if (t < n_snp - 1) {
            t++;
          } else {
            b++;
            t = 0;
          }
        } else {
          t++;
        }
      }
    }

    // Center mat_q.
    for (size_t i = 0; i < q->size; i++) {
      m = 0;
      n = 0;
      for (size_t k = 0; k < n_block; k++) {
        if (mat_s[k][i] != 0 && vec_s[i] != mat_s[k][i]) {
          d = (vec_q[i] - mat_q[k][i]) / (vec_s[i] - mat_s[k][i]);
          mat_q[k][i] = d;
          m += d;
          n++;
        }
      }
      if (n != 0) {
        m /= n;
      }

      for (size_t k = 0; k < n_block; k++) {
        if (mat_q[k][i] != 0) {
          mat_q[k][i] -= m;
        }
      }
    }

    // Compute Vq for l'th row and l'th column only.
    for (size_t i = 0; i < q->size; i++) {
      d = 0;
      n = 0;
      for (size_t k = 0; k < n_block; k++) {
        if (mat_q[k][l] != 0 && mat_q[k][i] != 0) {
          d += mat_q[k][l] * mat_q[k][i];
          n++;
        }
      }
      if (n != 0) {
        d /= n;
        d *= n - 1;
      }
      d += gsl_matrix_get(Vq, i, l);
      gsl_matrix_set(Vq, i, l, d);
      if (i != l) {
        gsl_matrix_set(Vq, l, i, d);
      }
    }
  }

  // divide the off diagonal elements of Vq by 2
  for (size_t i = 0; i < q->size; i++) {
    for (size_t j = i; j < q->size; j++) {
      if (i == j) {
        continue;
      }
      d = gsl_matrix_get(Vq, i, j);
      gsl_matrix_set(Vq, i, j, d / 2);
      gsl_matrix_set(Vq, j, i, d / 2);
    }
  }

  return;
}

// Read vector file.
void ReadFile_vector(const string &file_vec, gsl_vector *vec) {
  debug_msg("entered");
  igzstream infile(file_vec.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open vector file: " << file_vec << endl;
    return;
  }

  string line;
  char *ch_ptr;

  for (size_t i = 0; i < vec->size; i++) {
    safeGetline(infile, line).eof();
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",file_vec.c_str());
    gsl_vector_set(vec, i, atof(ch_ptr));
  }

  infile.clear();
  infile.close();

  return;
}

void ReadFile_matrix(const string &file_mat, gsl_matrix *mat) {
  debug_msg("entered");
  igzstream infile(file_mat.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open matrix file: " << file_mat << endl;
    return;
  }

  string line;
  char *ch_ptr;

  for (size_t i = 0; i < mat->size1; i++) {
    safeGetline(infile, line).eof();
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",file_mat.c_str());
    for (size_t j = 0; j < mat->size2; j++) {
      enforce(ch_ptr);
      gsl_matrix_set(mat, i, j, atof(ch_ptr));
      ch_ptr = strtok(NULL, " ,\t");
    }
  }

  infile.clear();
  infile.close();

  return;
}

void ReadFile_matrix(const string &file_mat, gsl_matrix *mat1,
                     gsl_matrix *mat2) {
  debug_msg("entered");
  igzstream infile(file_mat.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open matrix file: " << file_mat << endl;
    return;
  }

  string line;
  char *ch_ptr;

  for (size_t i = 0; i < mat1->size1; i++) {
    safeGetline(infile, line).eof();
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",file_mat.c_str());
    for (size_t j = 0; j < mat1->size2; j++) {
      enforce(ch_ptr);
      gsl_matrix_set(mat1, i, j, atof(ch_ptr));
      ch_ptr = strtok(NULL, " ,\t");
    }
  }

  for (size_t i = 0; i < mat2->size1; i++) {
    safeGetline(infile, line).eof();
    ch_ptr = strtok_safe2((char *)line.c_str(), " ,\t",file_mat.c_str());
    for (size_t j = 0; j < mat2->size2; j++) {
      enforce(ch_ptr);
      gsl_matrix_set(mat2, i, j, atof(ch_ptr));
      ch_ptr = strtok(NULL, " ,\t");
    }
  }

  infile.clear();
  infile.close();

  return;
}

// Read study file.
void ReadFile_study(const string &file_study, gsl_matrix *Vq_mat,
                    gsl_vector *q_vec, gsl_vector *s_vec, size_t &ni) {
  debug_msg("entered");
  string Vqfile = file_study + ".Vq.txt";
  string sfile = file_study + ".size.txt";
  string qfile = file_study + ".q.txt";

  gsl_vector *s = gsl_vector_safe_alloc(s_vec->size + 1);

  ReadFile_matrix(Vqfile, Vq_mat);
  ReadFile_vector(sfile, s);
  ReadFile_vector(qfile, q_vec);

  double d;
  for (size_t i = 0; i < s_vec->size; i++) {
    d = gsl_vector_get(s, i);
    gsl_vector_set(s_vec, i, d);
  }
  ni = gsl_vector_get(s, s_vec->size);

  gsl_vector_free(s);

  return;
}

// Read reference file.
void ReadFile_ref(const string &file_ref, gsl_matrix *S_mat,
                  gsl_matrix *Svar_mat, gsl_vector *s_vec, size_t &ni) {
  debug_msg("entered");
  string sfile = file_ref + ".size.txt";
  string Sfile = file_ref + ".S.txt";

  gsl_vector *s = gsl_vector_safe_alloc(s_vec->size + 1);

  ReadFile_vector(sfile, s);
  ReadFile_matrix(Sfile, S_mat, Svar_mat);

  double d;
  for (size_t i = 0; i < s_vec->size; i++) {
    d = gsl_vector_get(s, i);
    gsl_vector_set(s_vec, i, d);
  }
  ni = gsl_vector_get(s, s_vec->size);

  gsl_vector_free(s);

  return;
}

// Read mstudy file.
void ReadFile_mstudy(const string &file_mstudy, gsl_matrix *Vq_mat,
                     gsl_vector *q_vec, gsl_vector *s_vec, size_t &ni) {
  debug_msg("entered");
  gsl_matrix_set_zero(Vq_mat);
  gsl_vector_set_zero(q_vec);
  gsl_vector_set_zero(s_vec);
  ni = 0;

  gsl_matrix *Vq_sub = gsl_matrix_safe_alloc(Vq_mat->size1, Vq_mat->size2);
  gsl_vector *q_sub = gsl_vector_safe_alloc(q_vec->size);
  gsl_vector *s = gsl_vector_safe_alloc(s_vec->size + 1);

  igzstream infile(file_mstudy.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open mstudy file: " << file_mstudy << endl;
    return;
  }

  string file_name;
  double d1, d2, d;

  while (!safeGetline(infile, file_name).eof()) {
    string Vqfile = file_name + ".Vq.txt";
    string sfile = file_name + ".size.txt";
    string qfile = file_name + ".q.txt";

    ReadFile_matrix(Vqfile, Vq_sub);
    ReadFile_vector(sfile, s);
    ReadFile_vector(qfile, q_sub);

    ni = max(ni, (size_t)gsl_vector_get(s, s_vec->size));

    for (size_t i = 0; i < s_vec->size; i++) {
      d1 = gsl_vector_get(s, i);
      if (d1 == 0) {
        continue;
      }

      d = gsl_vector_get(q_vec, i) + gsl_vector_get(q_sub, i) * d1;
      gsl_vector_set(q_vec, i, d);

      d = gsl_vector_get(s_vec, i) + d1;
      gsl_vector_set(s_vec, i, d);

      for (size_t j = i; j < s_vec->size; j++) {
        d2 = gsl_vector_get(s, j);
        if (d2 == 0) {
          continue;
        }

        d = gsl_matrix_get(Vq_mat, i, j) +
            gsl_matrix_get(Vq_sub, i, j) * d1 * d2;
        gsl_matrix_set(Vq_mat, i, j, d);
        if (i != j) {
          gsl_matrix_set(Vq_mat, j, i, d);
        }
      }
    }
  }

  for (size_t i = 0; i < s_vec->size; i++) {
    d1 = gsl_vector_get(s_vec, i);
    if (d1 == 0) {
      continue;
    }

    d = gsl_vector_get(q_vec, i);
    gsl_vector_set(q_vec, i, d / d1);

    for (size_t j = i; j < s_vec->size; j++) {
      d2 = gsl_vector_get(s_vec, j);
      if (d2 == 0) {
        continue;
      }

      d = gsl_matrix_get(Vq_mat, i, j) / (d1 * d2);
      gsl_matrix_set(Vq_mat, i, j, d);
      if (i != j) {
        gsl_matrix_set(Vq_mat, j, i, d);
      }
    }
  }

  gsl_matrix_free(Vq_sub);
  gsl_vector_free(q_sub);
  gsl_vector_free(s);

  return;
}

// Read reference file.
void ReadFile_mref(const string &file_mref, gsl_matrix *S_mat,
                   gsl_matrix *Svar_mat, gsl_vector *s_vec, size_t &ni) {
  debug_msg("entered");
  gsl_matrix_set_zero(S_mat);
  gsl_matrix_set_zero(Svar_mat);
  gsl_vector_set_zero(s_vec);
  ni = 0;

  gsl_matrix *S_sub = gsl_matrix_safe_alloc(S_mat->size1, S_mat->size2);
  gsl_matrix *Svar_sub = gsl_matrix_safe_alloc(Svar_mat->size1, Svar_mat->size2);
  gsl_vector *s = gsl_vector_safe_alloc(s_vec->size + 1);

  igzstream infile(file_mref.c_str(), igzstream::in);
  if (!infile) {
    cout << "error! fail to open mref file: " << file_mref << endl;
    return;
  }

  string file_name;
  double d1, d2, d;

  while (!safeGetline(infile, file_name).eof()) {
    string sfile = file_name + ".size.txt";
    string Sfile = file_name + ".S.txt";

    ReadFile_vector(sfile, s);
    ReadFile_matrix(Sfile, S_sub, Svar_sub);

    // Update s_vec and ni.
    for (size_t i = 0; i < s_vec->size; i++) {
      d = gsl_vector_get(s, i) + gsl_vector_get(s_vec, i);
      gsl_vector_set(s_vec, i, d);
    }
    ni = max(ni, (size_t)gsl_vector_get(s, s_vec->size));

    // Update S and Svar from each file.
    for (size_t i = 0; i < S_mat->size1; i++) {
      d1 = gsl_vector_get(s, i);
      for (size_t j = 0; j < S_mat->size2; j++) {
        d2 = gsl_vector_get(s, j);

        d = gsl_matrix_get(S_sub, i, j) * d1 * d2;
        gsl_matrix_set(S_sub, i, j, d);
        d = gsl_matrix_get(Svar_sub, i, j) * d1 * d2 * d1 * d2;
        gsl_matrix_set(Svar_sub, i, j, d);
      }
    }

    gsl_matrix_add(S_mat, S_sub);
    gsl_matrix_add(Svar_mat, Svar_sub);
  }

  // Final: update S and Svar.
  for (size_t i = 0; i < S_mat->size1; i++) {
    d1 = gsl_vector_get(s_vec, i);
    if (d1 == 0) {
      continue;
    }
    for (size_t j = i; j < S_mat->size2; j++) {
      d2 = gsl_vector_get(s_vec, j);
      if (d2 == 0) {
        continue;
      }

      d = gsl_matrix_get(S_mat, i, j) / (d1 * d2);
      gsl_matrix_set(S_mat, i, j, d);
      if (i != j) {
        gsl_matrix_set(S_mat, j, i, d);
      }

      d = gsl_matrix_get(Svar_mat, i, j) / (d1 * d2 * d1 * d2);
      gsl_matrix_set(Svar_mat, i, j, d);
      if (i != j) {
        gsl_matrix_set(Svar_mat, j, i, d);
      }
    }
  }

  // Free matrices.
  gsl_matrix_free(S_sub);
  gsl_matrix_free(Svar_sub);
  gsl_vector_free(s);

  return;
}
