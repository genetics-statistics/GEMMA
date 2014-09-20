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

#ifndef __LM_H__                
#define __LM_H__

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"


#ifdef FORCE_FLOAT
#include "param_float.h"
#include "io_float.h"
#else
#include "param.h"
#include "io.h"
#endif

using namespace std;


class LM {
	
public:
	// IO related parameters
	int a_mode;				//analysis mode, 50+1/2/3/4 for Frequentist tests
	size_t d_pace;		//display pace
	
	string file_bfile;
	string file_geno;
	string file_out;
	
	string file_gene;
	
	// Summary statistics
	size_t ni_total, ni_test;	//number of individuals
	size_t ns_total, ns_test;	//number of snps
	size_t ng_total, ng_test;	//number of genes
	size_t n_cvt;
	double time_opt;		//time spent
	
	vector<int> indicator_idv;				//indicator for individuals (phenotypes), 0 missing, 1 available for analysis
	vector<int> indicator_snp;				//sequence indicator for SNPs: 0 ignored because of (a) maf, (b) miss, (c) non-poly; 1 available for analysis
	
	vector<SNPINFO> snpInfo;		//record SNP information
	
	// Not included in PARAM
	vector<SUMSTAT> sumStat;		//Output SNPSummary Data
	
	// Main functions
	void CopyFromParam (PARAM &cPar);
	void CopyToParam (PARAM &cPar);
	void AnalyzeGene (const gsl_matrix *W, const gsl_vector *x);
	void AnalyzePlink (const gsl_matrix *W, const gsl_vector *y);
	void AnalyzeBimbam (const gsl_matrix *W, const gsl_vector *y);
	void WriteFiles ();
};
void MatrixCalcLmLR (const gsl_matrix *X, const gsl_vector *y, vector<pair<size_t, double> > &pos_loglr);
#endif
