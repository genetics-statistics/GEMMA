/*
	Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright (C) 2011  Xiang Zhou

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

#ifndef __IO_H__
#define __IO_H__


#include <vector>
#include <map>
#include <algorithm>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

#include "gzstream.h"

#ifdef FORCE_FLOAT
#include "param_float.h"
#else
#include "param.h"
#endif

using namespace std;




void ProgressBar (string str, double p, double total);
void ProgressBar (string str, double p, double total, double ratio);
std::istream& safeGetline(std::istream& is, std::string& t);

bool ReadFile_snps (const string &file_snps, set<string> &setSnps);
bool ReadFile_log (const string &file_log, double &pheno_mean);

bool ReadFile_bim (const string &file_bim, vector<SNPINFO> &snpInfo);
bool ReadFile_fam (const string &file_fam, vector<vector<int> > &indicator_pheno, vector<vector<double> > &pheno, map<string, int> &mapID2num, const vector<size_t> &p_column);

bool ReadFile_cvt (const string &file_cvt, vector<int> &indicator_cvt, vector<vector<double> > &cvt, size_t &n_cvt);
bool ReadFile_anno (const string &file_bim, map<string, string> &mapRS2chr, map<string, long int> &mapRS2bp, map<string, double> &mapRS2cM);
bool ReadFile_pheno (const string &file_pheno, vector<vector<int> > &indicator_pheno, vector<vector<double> > &pheno, const vector<size_t> &p_column);
bool ReadFile_column (const string &file_pheno, vector<int> &indicator_idv, vector<double> &pheno, const int &p_column);

bool ReadFile_geno (const string &file_geno, const set<string> &setSnps, const gsl_matrix *W, vector<int> &indicator_idv, vector<int> &indicator_snp, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, map<string, string> &mapRS2chr, map<string, long int> &mapRS2bp, map<string, double> &mapRS2cM, vector<SNPINFO> &snpInfo, size_t &ns_test);
bool ReadFile_bed (const string &file_bed, const set<string> &setSnps, const gsl_matrix *W, vector<int> &indicator_idv, vector<int> &indicator_snp, vector<SNPINFO> &snpInfo, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, size_t &ns_test);
bool Bimbam_ReadOneSNP (const size_t inc, const vector<int> &indicator_idv, igzstream &infile, gsl_vector *geno, double &geno_mean);
void Plink_ReadOneSNP (const int pos, const vector<int> &indicator_idv, ifstream &infile, gsl_vector *geno, double &geno_mean);

void ReadFile_kin (const string &file_kin, vector<int> &indicator_idv, map<string, int> &mapID2num, const size_t k_mode, bool &error, gsl_matrix *G);
void ReadFile_mk (const string &file_mk, vector<int> &indicator_idv, map<string, int> &mapID2num, const size_t k_mode, bool &error, gsl_matrix *G);
void ReadFile_eigenU (const string &file_u, bool &error, gsl_matrix *U);
void ReadFile_eigenD (const string &file_d, bool &error, gsl_vector *eval);

bool BimbamKin (const string &file_geno, vector<int> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin);
bool PlinkKin (const string &file_bed, vector<int> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin);

bool ReadFile_geno (const string &file_geno, vector<int> &indicator_idv, vector<int> &indicator_snp, gsl_matrix *UtX, gsl_matrix *K, const bool calc_K);
bool ReadFile_bed (const string &file_bed, vector<int> &indicator_idv, vector<int> &indicator_snp, gsl_matrix *UtX, gsl_matrix *K, const bool calc_K);
bool ReadFile_geno (const string &file_geno, vector<int> &indicator_idv, vector<int> &indicator_snp, vector<vector<unsigned char> > &Xt, gsl_matrix *K, const bool calc_K, const size_t ni_test, const size_t ns_test);
bool ReadFile_bed (const string &file_bed, vector<int> &indicator_idv, vector<int> &indicator_snp, vector<vector<unsigned char> > &Xt, gsl_matrix *K, const bool calc_K, const size_t ni_test, const size_t ns_test);

bool ReadFile_est (const string &file_est, const vector<size_t> &est_column, map<string, double> &mapRS2est);

bool CountFileLines (const string &file_input, size_t &n_lines);

bool ReadFile_gene (const string &file_gene, vector<double> &vec_read, vector<SNPINFO> &snpInfo, size_t &ng_total);

bool ReadHeader (const string &line, HEADER &header);
bool ReadFile_cat (const string &file_cat, map<string, size_t> &mapRS2cat, size_t &n_vc);

bool BimbamKin (const string &file_geno, vector<int> &indicator_idv, vector<int> &indicator_snp, const int k_mode, const int display_pace, const map<string, size_t> &mapRS2cat, map<string, double> &mapRS2var, vector<SNPINFO> &snpInfo, gsl_matrix *matrix_kin);
bool PlinkKin (const string &file_bed, vector<int> &indicator_idv, vector<int> &indicator_snp, const int k_mode, const int display_pace, const map<string, size_t> &mapRS2cat, map<string, double> &mapRS2var, vector<SNPINFO> &snpInfo, gsl_matrix *matrix_kin);

bool ReadFile_var (const string &file_var, map<string, double> &mapRS2var);
void ReadFile_beta (const string &file_beta, const int k_mode, const map<string, size_t> &mapRS2cat, const map<string, double> &mapRS2var, gsl_vector *q, gsl_vector *s, size_t &ni_total, size_t &ns_total, size_t &ns_test);


void ReadFile_s (const string &file_s, gsl_matrix *S, gsl_matrix *Svar);
void ReadFile_ms (const string &file_ms, gsl_matrix *S, gsl_matrix *Svar);
void ReadFile_v (const string &file_v, gsl_matrix *V);
void ReadFile_mv (const string &file_mq, gsl_matrix *V);
void ReadFile_q (const string &file_s, gsl_vector *q_vec, gsl_vector *s_vec, double &df);
void ReadFile_mq (const string &file_mq, gsl_vector *q_vec, gsl_vector *s_vec, double &df);

// WJA added
bool bgenKin (const string &file_geno, vector<int> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin);
bool ReadFile_bgen(const string &file_bgen, const set<string> &setSnps, const gsl_matrix *W, vector<int> &indicator_idv, vector<int> &indicator_snp, vector<SNPINFO> &snpInfo, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, size_t &ns_test);
bool ReadFile_sample(const string &file_sample, vector<vector<int> > &indicator_pheno, vector<vector<double> > &pheno, const vector<size_t> &p_column, vector<int> &indicator_cvt, vector<vector<double> > &cvt, size_t &n_cvt);


#endif







