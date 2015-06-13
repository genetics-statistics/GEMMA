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

#ifndef __PARAM_H__                
#define __PARAM_H__

#include <vector>
#include <map>
#include <set>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

using namespace std;



class SNPINFO {
public:
	string chr;
	string rs_number;
	double cM;
	long int base_position;
	string a_minor;
	string a_major;
	size_t n_miss;
	double missingness;
	double maf;	
};

//results for lmm
class SUMSTAT {
public:
	double beta;				//REML estimator for beta
	double se;				//SE for beta  
	double lambda_remle;		//REML estimator for lambda
	double lambda_mle;		//MLE estimator for lambda
	double p_wald;			//p value from a Wald test
	double p_lrt;				//p value from a likelihood ratio test
	double p_score;			//p value from a score test
};

//results for mvlmm
class MPHSUMSTAT {
public:
	vector<double> v_beta;	//REML estimator for beta
	double p_wald;			//p value from a Wald test
	double p_lrt;				//p value from a likelihood ratio test
	double p_score;			//p value from a score test
	vector<double> v_Vg;	//estimator for Vg, right half
	vector<double> v_Ve;	//estimator for Ve, right half
	vector<double> v_Vbeta;	//estimator for Vbeta, right half
};


//hyper-parameters for bslmm
class HYPBSLMM {
public:
	double h;
	double pve;
	double rho;
	double pge;
	double logp;
	
	size_t n_gamma;
};




class PARAM {
public:	
	// IO related parameters
	bool mode_silence;
	int a_mode;				//analysis mode, 1/2/3/4 for Frequentist tests
	int k_mode;				//kinship read mode: 1: n by n matrix, 2: id/id/k_value; 		
	vector<size_t> p_column;			//which phenotype column needs analysis
	size_t d_pace;		//display pace
	
	string file_bfile;
	string file_geno;
	string file_pheno;
	string file_anno;		//optional
	string file_cvt;		//optional
	string file_kin;
	string file_ku, file_kd;
	string file_mk;
	string file_out;
	string path_out;
	
	string file_epm;		//estimated parameter file
	string file_ebv;		//estimated breeding value file
	string file_log;		//log file containing mean estimate
	
	string file_read;		//file containing total number of reads
	string file_gene;		//gene expression file
	
	string file_snps;		//file containing analyzed snps or genes

	// WJA Added	
	string file_bgenfile;	
	
	
	// QC related parameters	
	double miss_level;
	double maf_level;	
	double hwe_level;
	double r2_level;
	
	// LMM related parameters
	double l_min;
	double l_max;
	size_t n_region;
	double l_mle_null, l_remle_null;
	double logl_mle_H0, logl_remle_H0;
	double pve_null, pve_se_null;
	double vg_remle_null, ve_remle_null, vg_mle_null, ve_mle_null;
	vector<double> Vg_remle_null, Ve_remle_null, Vg_mle_null, Ve_mle_null;
	vector<double> VVg_remle_null, VVe_remle_null, VVg_mle_null, VVe_mle_null;
	vector<double> beta_remle_null, se_beta_remle_null, beta_mle_null, se_beta_mle_null;
	double p_nr;	
	double em_prec, nr_prec;
	size_t em_iter, nr_iter;
	size_t crt;
	double pheno_mean;		//phenotype mean from bslmm fitting or for prediction

	//for fitting multiple variance components
	//the first three are of size n_vc, and the next two are of size n_vc+1
	vector<double> v_traceG;
	vector<double> v_pve;
	vector<double> v_se_pve;

	vector<double> v_sigma2;
	vector<double> v_se_sigma2;	
	vector<double> v_beta;
	vector<double> v_se_beta;	
	
	// BSLMM MCMC related parameters
	double h_min, h_max, h_scale;			//priors for h
	double rho_min, rho_max, rho_scale;		//priors for rho
	double logp_min, logp_max, logp_scale;		//priors for log(pi)
	size_t s_min, s_max;			//minimum and maximum number of gammas
	size_t w_step;					//number of warm up/burn in iterations
	size_t s_step;					//number of sampling iterations
	size_t r_pace;					//record pace
	size_t w_pace;					//write pace
	size_t n_accept;				//number of acceptance
	size_t n_mh;					//number of MH steps within each iteration
	double geo_mean;				//mean of the geometric distribution
	long int randseed;
	double trace_G;

	HYPBSLMM cHyp_initial;
		
	// Summary statistics
	bool error;
	size_t ni_total, ni_test, ni_cvt;	//number of individuals
	size_t np_obs, np_miss;		//number of observed and missing phenotypes
	size_t ns_total, ns_test;	//number of snps
	size_t ng_total, ng_test;	//number of genes
	size_t ni_control, ni_case;	//number of controls and number of cases
	size_t n_cvt;			//number of covariates
	size_t n_ph;			//number of phenotypes
	size_t n_vc;			//number of variance components (including the diagonal matrix)
	double time_total;		//record total time
	double time_G;			//time spent on reading files the second time and calculate K
	double time_eigen;		//time spent on eigen-decomposition
	double time_UtX;		//time spent on calculating UX and Uy
	double time_UtZ;		//time spent on calculating UtZ, for probit BSLMM
	double time_opt;		//time spent on optimization iterations/or mcmc
	double time_Omega;		//time spent on calculating Omega
	double time_hyp;		//time spent on sampling hyper-parameters, in PMM
	double time_Proposal;  //time spend on constructing the proposal distribution (i.e. the initial lmm or lm analysis)

	// Data
	vector<vector<double> > pheno;			//a vector record all phenotypes, NA replaced with -9
	vector<vector<double> > cvt;			//a vector record all covariates, NA replaced with -9	
	vector<vector<int> > indicator_pheno;			//a matrix record when a phenotype is missing for an individual; 0 missing, 1 available
	vector<int> indicator_idv;				//indicator for individuals (phenotypes), 0 missing, 1 available for analysis
	vector<int> indicator_snp;				//sequence indicator for SNPs: 0 ignored because of (a) maf, (b) miss, (c) non-poly; 1 available for analysis
	vector<int> indicator_cvt;				//indicator for covariates, 0 missing, 1 available for analysis
	
	vector<int> indicator_bv;				//indicator for estimated breeding value file, 0 missing, 1 available for analysis
	vector<int> indicator_read;				//indicator for read file, 0 missing, 1 available for analysis
	vector<double> vec_read;				//total number of reads
	vector<double> vec_bv;					//breeding values
	vector<size_t> est_column;
	
	map<string, int> mapID2num;		//map small ID number to number, from 0 to n-1
	map<string, string> mapRS2chr;		//map rs# to chromosome location
	map<string, long int> mapRS2bp;		//map rs# to base position
	map<string, double> mapRS2cM;		//map rs# to cM
	map<string, double> mapRS2est;			//map rs# to parameters
	
	vector<SNPINFO> snpInfo;		//record SNP information
	set<string> setSnps;			//a set of snps for analysis
	
	//constructor
	PARAM();
	
	//functions
	void ReadFiles ();		
	void CheckParam (); 
	void CheckData ();	
	void PrintSummary ();
	void ReadGenotypes (gsl_matrix *UtX, gsl_matrix *K, const bool calc_K);	
	void CheckCvt ();
	void CopyCvt (gsl_matrix *W);
	void ProcessCvtPhen();
	void CopyCvtPhen (gsl_matrix *W, gsl_vector *y, size_t flag);
	void CopyCvtPhen (gsl_matrix *W, gsl_matrix *Y, size_t flag);
	void CalcKin (gsl_matrix *matrix_kin);
	void WriteMatrix (const gsl_matrix *matrix_U, const string suffix);
	void WriteVector (const gsl_vector *vector_D, const string suffix);
	void CopyRead (gsl_vector *log_N);
};


#endif

