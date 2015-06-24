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

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <ctime>
#include <cmath>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_cdf.h"

#include "lapack.h"  //for functions EigenDecomp

#ifdef FORCE_FLOAT
#include "io_float.h"   //for function ReadFile_kin
#include "gemma_float.h"
#include "vc_float.h"
#include "lm_float.h"  //for LM class
#include "bslmm_float.h"  //for BSLMM class
#include "lmm_float.h"  //for LMM class, and functions CalcLambda, CalcPve, CalcVgVe
#include "mvlmm_float.h"  //for MVLMM class
#include "prdt_float.h"	//for PRDT class
#include "mathfunc_float.h"	//for a few functions
#else
#include "io.h"
#include "gemma.h"
#include "vc.h"
#include "lm.h"
#include "bslmm.h"
#include "lmm.h"
#include "mvlmm.h"
#include "prdt.h"
#include "mathfunc.h"
#endif


using namespace std;



GEMMA::GEMMA(void):	
version("0.95alpha"), date("08/08/2014"), year("2011")
{}

void GEMMA::PrintHeader (void)
{
	cout<<endl;
	cout<<"*********************************************************"<<endl;
	cout<<"  Genome-wide Efficient Mixed Model Association (GEMMA) "<<endl;
	cout<<"  Version "<<version<<", "<<date<<"                              "<<endl;
	cout<<"  Visit                                                 "<<endl;
	cout<<"     http://stephenslab.uchicago.edu/software.html      "<<endl;
	cout<<"     http://home.uchicago.edu/~xz7/software.html        "<<endl;
	cout<<"  For Possible Updates                                  "<<endl;
	cout<<"  (C) "<<year<<" Xiang Zhou                                   "<<endl;
	cout<<"  GNU General Public License                            "<<endl;
	cout<<"  For Help, Type ./gemma -h                             "<<endl;
	cout<<"*********************************************************"<<endl;
	cout<<endl;
	
	return;
}


void GEMMA::PrintLicense (void)
{
	cout<<endl;
	cout<<"The Software Is Distributed Under GNU General Public License, But May Also Require The Following Notifications."<<endl;
	cout<<endl;
	
	cout<<"Including Lapack Routines In The Software May Require The Following Notification:"<<endl;
	cout<<"Copyright (c) 1992-2010 The University of Tennessee and The University of Tennessee Research Foundation.  All rights reserved."<<endl;
	cout<<"Copyright (c) 2000-2010 The University of California Berkeley. All rights reserved."<<endl;
	cout<<"Copyright (c) 2006-2010 The University of Colorado Denver.  All rights reserved."<<endl;	
	cout<<endl;
	
	cout<<"$COPYRIGHT$"<<endl;
	cout<<"Additional copyrights may follow"<<endl;
	cout<<"$HEADER$"<<endl;
	cout<<"Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:"<<endl;
	cout<<"- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer."<<endl;
	cout<<"- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer listed in this license in the documentation and/or other materials provided with the distribution."<<endl;
	cout<<"- Neither the name of the copyright holders nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission."<<endl;
	cout<<"The copyright holders provide no reassurances that the source code provided does not infringe any patent, copyright, or any other "
		<<"intellectual property rights of third parties.  The copyright holders disclaim any liability to any recipient for claims brought against "
		<<"recipient by any third party for infringement of that parties intellectual property rights. "<<endl;
	cout<<"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT "
		<<"LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT "
		<<"OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT "
		<<"LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY "
		<<"THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE "
		<<"OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."<<endl;
	cout<<endl;
	
	
	
	return;
}



void GEMMA::PrintHelp(size_t option)
{
	if (option==0) {
		cout<<endl; 
		cout<<" GEMMA version "<<version<<", released on "<<date<<endl;
		cout<<" implemented by Xiang Zhou"<<endl; 
		cout<<endl;
		cout<<" type ./gemma -h [num] for detailed helps"<<endl;
		cout<<" options: " << endl;
		cout<<" 1: quick guide"<<endl;
		cout<<" 2: file I/O related"<<endl;
		cout<<" 3: SNP QC"<<endl;
		cout<<" 4: calculate relatedness matrix"<<endl;
		cout<<" 5: perform eigen decomposition"<<endl;
		cout<<" 6: perform variance component estiamtion"<<endl;
		cout<<" 7: fit a linear model"<<endl;
		cout<<" 8: fit a linear mixed model"<<endl;
		cout<<" 9: fit a multivariate linear mixed model"<<endl;
		cout<<" 10: fit a Bayesian sparse linear mixed model"<<endl;
		cout<<" 11: obtain predicted values"<<endl;
		cout<<" 12: note"<<endl;
		cout<<endl;
	}	
	
	if (option==1) {
		cout<<" QUICK GUIDE" << endl;
		cout<<" to generate a relatedness matrix: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -gk [num] -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -gk [num] -o [prefix]"<<endl;
		cout<<" to perform eigen decomposition of the relatedness matrix: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -k [filename] -eigen -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -k [filename] -eigen -o [prefix]"<<endl;
		cout<<" to estimate variance components: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -k [filename] -vc -o [prefix]"<<endl;
		cout<<"         ./gemma -p [filename] -k [filename] -vc -o [prefix]"<<endl;
		cout<<"         ./gemma -bfile [prefix] -mk [filename] -vc -o [prefix]"<<endl;
		cout<<"         ./gemma -p [filename] -mk [filename] -vc -o [prefix]"<<endl;
		cout<<" to fit a linear mixed model: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -k [filename] -lmm [num] -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -a [filename] -k [filename] -lmm [num] -o [prefix]"<<endl;	
		cout<<" to fit a multivariate linear mixed model: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -k [filename] -lmm [num] -n [num1] [num2] -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -a [filename] -k [filename] -lmm [num] -n [num1] [num2] -o [prefix]"<<endl;	
		cout<<" to fit a Bayesian sparse linear mixed model: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -bslmm [num] -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -a [filename] -bslmm [num] -o [prefix]"<<endl;
		cout<<" to obtain predicted values: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -epm [filename] -emu [filename] -ebv [filename] -k [filename] -predict [num] -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -epm [filename] -emu [filename] -ebv [filename] -k [filename] -predict [num] -o [prefix]"<<endl;
		cout<<endl;
	}
	
	if (option==2) {
		cout<<" FILE I/O RELATED OPTIONS" << endl;
		cout<<" -bfile    [prefix]       "<<" specify input PLINK binary ped file prefix."<<endl;	
		cout<<"          requires: *.fam, *.bim and *.bed files"<<endl;
		// WJA added
		cout<<" -bgenfile    [prefix]       "<<" specify input Oxford genotype bgen file prefix."<<endl;	
		cout<<"          requires: *.bggen, *.sample files"<<endl;	
	
		cout<<"          missing value: -9"<<endl;
		cout<<" -g        [filename]     "<<" specify input BIMBAM mean genotype file name"<<endl;
		cout<<"          format: rs#1, allele0, allele1, genotype for individual 1, genotype for individual 2, ..."<<endl;	
		cout<<"                  rs#2, allele0, allele1, genotype for individual 1, genotype for individual 2, ..."<<endl;	
		cout<<"                  ..."<<endl;	
		cout<<"          missing value: NA"<<endl;	
		cout<<" -p        [filename]     "<<" specify input BIMBAM phenotype file name"<<endl;
		cout<<"          format: phenotype for individual 1"<<endl;	
		cout<<"                  phenotype for individual 2"<<endl;	
		cout<<"                  ..."<<endl;
		cout<<"          missing value: NA"<<endl;	
		cout<<" -a        [filename]     "<<" specify input BIMBAM SNP annotation file name (optional)"<<endl;	
		cout<<"          format: rs#1, base_position, chr_number"<<endl;	
		cout<<"                  rs#2, base_position, chr_number"<<endl;	
		cout<<"                  ..."<<endl;
		cout<<" -k        [filename]     "<<" specify input kinship/relatedness matrix file name"<<endl;	
		cout<<" -mk       [filename]     "<<" specify input file which contains a list of kinship/relatedness matrices"<<endl;	
		cout<<" -u        [filename]     "<<" specify input file containing the eigen vectors of the kinship/relatedness matrix"<<endl;	
		cout<<" -d        [filename]     "<<" specify input file containing the eigen values of the kinship/relatedness matrix"<<endl;	
		cout<<" -c        [filename]     "<<" specify input covariates file name (optional)"<<endl;	
		cout<<"          format: covariate 1 for individual 1, ... , covariate c for individual 1"<<endl;	
		cout<<"                  covariate 1 for individual 2, ... , covariate c for individual 2"<<endl;	
		cout<<"                  ..."<<endl;
		cout<<"          missing value: NA"<<endl;	
		cout<<"          note: the intercept (a column of 1s) may need to be included"<<endl;
		cout<<" -epm      [filename]     "<<" specify input estimated parameter file name"<<endl;
		cout<<" -en [n1] [n2] [n3] [n4]  "<<" specify values for the input estimated parameter file (with a header)"<<endl;
		cout<<"          options: n1: rs column number"<<endl;
		cout<<"                   n2: estimated alpha column number (0 to ignore)"<<endl;
		cout<<"                   n3: estimated beta column number (0 to ignore)"<<endl;
		cout<<"                   n4: estimated gamma column number (0 to ignore)"<<endl;
		cout<<"          default: 2 4 5 6 if -ebv is not specified; 2 0 5 6 if -ebv is specified"<<endl;
		cout<<" -ebv      [filename]     "<<" specify input estimated random effect (breeding value) file name"<<endl;
		cout<<"          format: value for individual 1"<<endl;	
		cout<<"                  value for individual 2"<<endl;	
		cout<<"                  ..."<<endl;
		cout<<"          missing value: NA"<<endl;	
		cout<<" -emu      [filename]     "<<" specify input log file name containing estimated mean"<<endl;
		cout<<" -mu       [num]          "<<" specify input estimated mean value"<<endl;
		cout<<" -gene     [filename]     "<<" specify input gene expression file name"<<endl;
		cout<<"          format: header"<<endl;	
		cout<<"                  gene1, count for individual 1, count for individual 2, ..."<<endl;	
		cout<<"                  gene2, count for individual 1, count for individual 2, ..."<<endl;	
		cout<<"                  ..."<<endl;
		cout<<"          missing value: not allowed"<<endl;	
		cout<<" -r        [filename]     "<<" specify input total read count file name"<<endl;
		cout<<"          format: total read count for individual 1"<<endl;	
		cout<<"                  total read count for individual 2"<<endl;	
		cout<<"                  ..."<<endl;
		cout<<"          missing value: NA"<<endl;	
		cout<<" -snps     [filename]     "<<" specify input snps file name to only analyze a certain set of snps"<<endl;
		cout<<"          format: rs#1"<<endl;	
		cout<<"                  rs#2"<<endl;	
		cout<<"                  ..."<<endl;
		cout<<"          missing value: NA"<<endl;	
		cout<<" -silence                 "<<" silent terminal display"<<endl;
		cout<<" -km       [num]          "<<" specify input kinship/relatedness file type (default 1)."<<endl;
		cout<<"          options: 1: \"n by n matrix\" format"<<endl;
		cout<<"                   2: \"id  id  value\" format"<<endl;
		cout<<" -n        [num]          "<<" specify phenotype column in the phenotype/*.fam file (optional; default 1)"<<endl;	
		cout<<" -pace     [num]          "<<" specify terminal display update pace (default 100000 SNPs or 100000 iterations)."<<endl;
		cout<<" -outdir   [path]         "<<" specify output directory path (default \"./output/\")"<<endl; 
		cout<<" -o        [prefix]       "<<" specify output file prefix (default \"result\")"<<endl;  
		cout<<"          output: prefix.cXX.txt or prefix.sXX.txt from kinship/relatedness matrix estimation"<<endl;	
		cout<<"          output: prefix.assoc.txt and prefix.log.txt form association tests"<<endl;	
		cout<<endl;
	}
	
	if (option==3) {
		cout<<" SNP QC OPTIONS" << endl;
		cout<<" -miss     [num]          "<<" specify missingness threshold (default 0.05)" << endl; 
		cout<<" -maf      [num]          "<<" specify minor allele frequency threshold (default 0.01)" << endl; 
		cout<<" -hwe      [num]          "<<" specify HWE test p value threshold (default 0; no test)" << endl; 
		cout<<" -r2       [num]          "<<" specify r-squared threshold (default 0.9999)" << endl; 
		cout<<" -notsnp                  "<<" minor allele frequency cutoff is not used" << endl; 
		cout<<endl;
	}
	
	if (option==4) {
		cout<<" RELATEDNESS MATRIX CALCULATION OPTIONS" << endl;
		cout<<" -gk       [num]          "<<" specify which type of kinship/relatedness matrix to generate (default 1)" << endl; 
		cout<<"          options: 1: centered XX^T/p"<<endl;
		cout<<"                   2: standardized XX^T/p"<<endl;
		cout<<"          note: non-polymorphic SNPs are excluded "<<endl;
		cout<<endl;
	}
	
	if (option==5) {
		cout<<" EIGEN-DECOMPOSITION OPTIONS" << endl;
		cout<<" -eigen                   "<<" specify to perform eigen decomposition of the loaded relatedness matrix" << endl; 
		cout<<endl;
	}

	if (option==6) {
		cout<<" VARIANCE COMPONENT ESTIMATION OPTIONS" << endl;
		cout<<" -vc                      "<<" specify to perform variance component estimation for the loaded relatedness matrix/matrices" << endl; 
		cout<<endl;
	}
	
	if (option==7) {
		cout<<" LINEAR MODEL OPTIONS" << endl;		
		cout<<" -lm       [num]         "<<" specify analysis options (default 1)."<<endl;
		cout<<"          options: 1: Wald test"<<endl;
		cout<<"                   2: Likelihood ratio test"<<endl;
		cout<<"                   3: Score test"<<endl;
		cout<<"                   4: 1-3"<<endl;
		cout<<endl;
	}
	
	if (option==8) {
		cout<<" LINEAR MIXED MODEL OPTIONS" << endl;		
		cout<<" -lmm      [num]         "<<" specify analysis options (default 1)."<<endl;
		cout<<"          options: 1: Wald test"<<endl;		
		cout<<"                   2: Likelihood ratio test"<<endl;
		cout<<"                   3: Score test"<<endl;
		cout<<"                   4: 1-3"<<endl;
		cout<<"                   5: Parameter estimation in the null model only"<<endl;
		cout<<" -lmin     [num]          "<<" specify minimal value for lambda (default 1e-5)" << endl; 
		cout<<" -lmax     [num]          "<<" specify maximum value for lambda (default 1e+5)" << endl; 
		cout<<" -region   [num]          "<<" specify the number of regions used to evaluate lambda (default 10)" << endl; 
		cout<<endl;
	}
	
	if (option==9) {
		cout<<" MULTIVARIATE LINEAR MIXED MODEL OPTIONS" << endl;
		cout<<" -pnr				     "<<" specify the pvalue threshold to use the Newton-Raphson's method (default 0.001)"<<endl;
		cout<<" -emi				     "<<" specify the maximum number of iterations for the PX-EM method in the null (default 10000)"<<endl;
		cout<<" -nri				     "<<" specify the maximum number of iterations for the Newton-Raphson's method in the null (default 100)"<<endl;
		cout<<" -emp				     "<<" specify the precision for the PX-EM method in the null (default 0.0001)"<<endl;
		cout<<" -nrp				     "<<" specify the precision for the Newton-Raphson's method in the null (default 0.0001)"<<endl;
		cout<<" -crt				     "<<" specify to output corrected pvalues for these pvalues that are below the -pnr threshold"<<endl;
		cout<<endl;
	}
	
	if (option==10) {
		cout<<" MULTI-LOCUS ANALYSIS OPTIONS" << endl;
		cout<<" -bslmm	  [num]			 "<<" specify analysis options (default 1)."<<endl;
		cout<<"          options: 1: BSLMM"<<endl;	
		cout<<"                   2: standard ridge regression/GBLUP (no mcmc)"<<endl;	
		cout<<"                   3: probit BSLMM (requires 0/1 phenotypes)"<<endl;			
		
		cout<<"   MCMC OPTIONS" << endl;
		cout<<"   Prior" << endl;	
		cout<<" -hmin     [num]          "<<" specify minimum value for h (default 0)" << endl; 
		cout<<" -hmax     [num]          "<<" specify maximum value for h (default 1)" << endl; 
		cout<<" -rmin     [num]          "<<" specify minimum value for rho (default 0)" << endl; 
		cout<<" -rmax     [num]          "<<" specify maximum value for rho (default 1)" << endl; 
		cout<<" -pmin     [num]          "<<" specify minimum value for log10(pi) (default log10(1/p), where p is the number of analyzed SNPs )" << endl; 
		cout<<" -pmax     [num]          "<<" specify maximum value for log10(pi) (default log10(1) )" << endl; 	
		cout<<" -smin     [num]          "<<" specify minimum value for |gamma| (default 0)" << endl; 
		cout<<" -smax     [num]          "<<" specify maximum value for |gamma| (default 300)" << endl; 
		
		cout<<"   Proposal" << endl;
		cout<<" -gmean    [num]          "<<" specify the mean for the geometric distribution (default: 2000)" << endl; 
		cout<<" -hscale   [num]          "<<" specify the step size scale for the proposal distribution of h (value between 0 and 1, default min(10/sqrt(n),1) )" << endl; 
		cout<<" -rscale   [num]          "<<" specify the step size scale for the proposal distribution of rho (value between 0 and 1, default min(10/sqrt(n),1) )" << endl; 
		cout<<" -pscale   [num]          "<<" specify the step size scale for the proposal distribution of log10(pi) (value between 0 and 1, default min(5/sqrt(n),1) )" << endl; 
		
		cout<<"   Others" << endl;
		cout<<" -w        [num]          "<<" specify burn-in steps (default 100,000)" << endl; 
		cout<<" -s        [num]          "<<" specify sampling steps (default 1,000,000)" << endl; 
		cout<<" -rpace    [num]          "<<" specify recording pace, record one state in every [num] steps (default 10)" << endl; 	
		cout<<" -wpace    [num]          "<<" specify writing pace, write values down in every [num] recorded steps (default 1000)" << endl; 	
		cout<<" -seed     [num]          "<<" specify random seed (a random seed is generated by default)" << endl; 	
		cout<<" -mh       [num]          "<<" specify number of MH steps in each iteration (default 10)" << endl; 
		cout<<"          requires: 0/1 phenotypes and -bslmm 3 option"<<endl;	
		cout<<endl;
	}
	
	if (option==11) {
		cout<<" PREDICTION OPTIONS" << endl;
		cout<<" -predict  [num]			 "<<" specify prediction options (default 1)."<<endl;
		cout<<"          options: 1: predict for individuals with missing phenotypes"<<endl;	
		cout<<"                   2: predict for individuals with missing phenotypes, and convert the predicted values to probability scale. Use only for files fitted with -bslmm 3 option"<<endl;	
		cout<<endl;
	}
	
	if (option==12) {
		cout<<" NOTE"<<endl;
		cout<<" 1. Only individuals with non-missing phenotoypes and covariates will be analyzed."<<endl;
		cout<<" 2. Missing genotoypes will be repalced with the mean genotype of that SNP."<<endl;
		cout<<" 3. For lmm analysis, memory should be large enough to hold the relatedness matrix and to perform eigen decomposition."<<endl;
		cout<<" 4. For multivariate lmm analysis, use a large -pnr for each snp will increase computation time dramatically."<<endl;
		cout<<" 5. For bslmm analysis, in addition to 3, memory should be large enough to hold the whole genotype matrix."<<endl;
		cout<<endl;
	}
	
	return;
}



void GEMMA::Assign(int argc, char ** argv, PARAM &cPar)
{
	string str;
	
	for(int i = 1; i < argc; i++) {		
		if (strcmp(argv[i], "-bfile")==0 || strcmp(argv[i], "--bfile")==0 || strcmp(argv[i], "-b")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_bfile=str;
		}
		// WJA added	
		else if (strcmp(argv[i], "-bgenfile")==0 || strcmp(argv[i], "--bgenfile")==0 || strcmp(argv[i], "-bgen")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_bgenfile=str;
		}

		else if (strcmp(argv[i], "-silence")==0) {
			cPar.mode_silence=true;
		}
		else if (strcmp(argv[i], "-g")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_geno=str;
		}
		else if (strcmp(argv[i], "-p")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_pheno=str;
		}
		else if (strcmp(argv[i], "-a")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_anno=str;
		}
		else if (strcmp(argv[i], "-k")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_kin=str;
		}
		else if (strcmp(argv[i], "-mk")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_mk=str;
		}
		else if (strcmp(argv[i], "-u")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_ku=str;
		}
		else if (strcmp(argv[i], "-d")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_kd=str;
		}
		else if (strcmp(argv[i], "-c")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_cvt=str;
		}
		else if (strcmp(argv[i], "-epm")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_epm=str;
		}
		else if (strcmp(argv[i], "-en")==0) {			
			while (argv[i+1] != NULL && argv[i+1][0] != '-') {
				++i;
				str.clear();
				str.assign(argv[i]);
				cPar.est_column.push_back(atoi(str.c_str()));
			}
		}
		else if (strcmp(argv[i], "-ebv")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_ebv=str;
		}
		else if (strcmp(argv[i], "-emu")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_log=str;
		}
		else if (strcmp(argv[i], "-mu")==0) {
			if(argv[i+1] == NULL) {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.pheno_mean=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-gene")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_gene=str;
		}
		else if (strcmp(argv[i], "-r")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_read=str;
		}
		else if (strcmp(argv[i], "-snps")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_snps=str;
		}
		else if (strcmp(argv[i], "-km")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.k_mode=atoi(str.c_str());
		}		
		else if (strcmp(argv[i], "-n")==0) {
			(cPar.p_column).clear();
			while (argv[i+1] != NULL && argv[i+1][0] != '-') {
				++i;
				str.clear();
				str.assign(argv[i]);
				(cPar.p_column).push_back(atoi(str.c_str()));
			}
		}
		else if (strcmp(argv[i], "-pace")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.d_pace=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-outdir")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.path_out=str;
		}
		else if (strcmp(argv[i], "-o")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_out=str;
		}		
		else if (strcmp(argv[i], "-miss")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.miss_level=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-maf")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			if (cPar.maf_level!=-1) {cPar.maf_level=atof(str.c_str());}
		}
		else if (strcmp(argv[i], "-hwe")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.hwe_level=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-r2")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.r2_level=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-notsnp")==0) {
			cPar.maf_level=-1;
		}
		else if (strcmp(argv[i], "-gk")==0) {
			if (cPar.a_mode!=0) {cPar.error=true; cout<<"error! only one of -gk -eigen -vc -lm -lmm -bslmm -predict options is allowed."<<endl; break;}
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.a_mode=21; continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.a_mode=20+atoi(str.c_str());
		}	
		else if (strcmp(argv[i], "-eigen")==0) {
			if (cPar.a_mode!=0) {cPar.error=true; cout<<"error! only one of -gk -eigen -vc -lm -lmm -bslmm -predict options is allowed."<<endl; break;}
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.a_mode=31; continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.a_mode=30+atoi(str.c_str());
		}	
		else if (strcmp(argv[i], "-vc")==0) {
			if (cPar.a_mode!=0) {cPar.error=true; cout<<"error! only one of -gk -eigen -vc -lm -lmm -bslmm -predict options is allowed."<<endl; break;}
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.a_mode=61; continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.a_mode=60+atoi(str.c_str());
		}	
		else if (strcmp(argv[i], "-lm")==0) {
			if (cPar.a_mode!=0) {cPar.error=true; cout<<"error! only one of -gk -eigen -vc -lm -lmm -bslmm -predict options is allowed."<<endl; break;}
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.a_mode=51; continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.a_mode=50+atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-fa")==0 || strcmp(argv[i], "-lmm")==0) {
			if (cPar.a_mode!=0) {cPar.error=true; cout<<"error! only one of -gk -eigen -vc -lm -lmm -bslmm -predict options is allowed."<<endl; break;}
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.a_mode=1; continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.a_mode=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-lmin")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.l_min=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-lmax")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.l_max=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-region")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.n_region=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-pnr")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.p_nr=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-emi")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.em_iter=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-nri")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.nr_iter=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-emp")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.em_prec=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-nrp")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.nr_prec=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-crt")==0) {
			cPar.crt=1;
		}
		else if (strcmp(argv[i], "-bslmm")==0) {
			if (cPar.a_mode!=0) {cPar.error=true; cout<<"error! only one of -gk -eigen -vc -lm -lmm -bslmm -predict options is allowed."<<endl; break;}
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.a_mode=11; continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.a_mode=10+atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-hmin")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.h_min=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-hmax")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.h_max=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-rmin")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.rho_min=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-rmax")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.rho_max=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-pmin")==0) {
			if(argv[i+1] == NULL) {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.logp_min=atof(str.c_str())*log(10.0);
		}
		else if (strcmp(argv[i], "-pmax")==0) {
			if(argv[i+1] == NULL) {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.logp_max=atof(str.c_str())*log(10.0);
		}
		else if (strcmp(argv[i], "-smin")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.s_min=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-smax")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.s_max=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-gmean")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.geo_mean=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-hscale")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.h_scale=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-rscale")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.rho_scale=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-pscale")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.logp_scale=atof(str.c_str())*log(10.0);
		}
		else if (strcmp(argv[i], "-w")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.w_step=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-s")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.s_step=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-rpace")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.r_pace=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-wpace")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.w_pace=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-seed")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.randseed=atol(str.c_str());
		}
		else if (strcmp(argv[i], "-mh")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.n_mh=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-predict")==0) {
			if (cPar.a_mode!=0) {cPar.error=true; cout<<"error! only one of -gk -eigen -vc -lm -lmm -bslmm -predict options is allowed."<<endl; break;}
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.a_mode=41; continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.a_mode=40+atoi(str.c_str());
		}
		else {cout<<"error! unrecognized option: "<<argv[i]<<endl; cPar.error=true; continue;}
	}
	
	//change prediction mode to 43, if the epm file is not provided
	if (cPar.a_mode==41 && cPar.file_epm.empty()) {cPar.a_mode=43;}
	
	return;
}



void GEMMA::BatchRun (PARAM &cPar) 
{
	clock_t time_begin, time_start;
	time_begin=clock();

	//Read Files
	cout<<"Reading Files ... "<<endl;
	cPar.ReadFiles();
	if (cPar.error==true) {cout<<"error! fail to read files. "<<endl; return;}

	cPar.CheckData();

	if (cPar.error==true) {cout<<"error! fail to check data. "<<endl; return;}

	//Prediction for bslmm	
	if (cPar.a_mode==41 || cPar.a_mode==42) {
		gsl_vector *y_prdt;
		
		y_prdt=gsl_vector_alloc (cPar.ni_total-cPar.ni_test);

		//set to zero
		gsl_vector_set_zero (y_prdt);
		
		PRDT cPRDT;
		cPRDT.CopyFromParam(cPar);
		
		//add breeding value if needed
		if (!cPar.file_kin.empty() && !cPar.file_ebv.empty()) {
			cout<<"Adding Breeding Values ... "<<endl;
			
			gsl_matrix *G=gsl_matrix_alloc (cPar.ni_total, cPar.ni_total);
			gsl_vector *u_hat=gsl_vector_alloc (cPar.ni_test);
			
			//read kinship matrix and set u_hat
			vector<int> indicator_all;
			size_t c_bv=0;
			for (size_t i=0; i<cPar.indicator_idv.size(); i++) {
				indicator_all.push_back(1);
				if (cPar.indicator_bv[i]==1) {gsl_vector_set(u_hat, c_bv, cPar.vec_bv[i]); c_bv++;}
			}
			
			ReadFile_kin (cPar.file_kin, indicator_all, cPar.mapID2num, cPar.k_mode, cPar.error, G);
			if (cPar.error==true) {cout<<"error! fail to read kinship/relatedness file. "<<endl; return;}
			
			//read u			
			cPRDT.AddBV(G, u_hat, y_prdt);					
			
			gsl_matrix_free(G);
			gsl_vector_free(u_hat);
		}

		//add beta
		if (!cPar.file_bfile.empty()) {
			cPRDT.AnalyzePlink (y_prdt);
		}
		else {
			cPRDT.AnalyzeBimbam (y_prdt);
		}
		
		//add mu
		gsl_vector_add_constant(y_prdt, cPar.pheno_mean);
		
		//convert y to probability if needed
		if (cPar.a_mode==42) {
			double d;
			for (size_t i=0; i<y_prdt->size; i++) {
				d=gsl_vector_get(y_prdt, i);
				d=gsl_cdf_gaussian_P(d, 1.0);
				gsl_vector_set(y_prdt, i, d);
			}
		}
			
			
		cPRDT.CopyToParam(cPar);
		
		cPRDT.WriteFiles(y_prdt);
		
		gsl_vector_free(y_prdt);
	}
	
	
	//Prediction with kinship matrix only; for one or more phenotypes
	if (cPar.a_mode==43) {
		//first, use individuals with full phenotypes to obtain estimates of Vg and Ve		
		gsl_matrix *Y=gsl_matrix_alloc (cPar.ni_test, cPar.n_ph);
		gsl_matrix *W=gsl_matrix_alloc (Y->size1, cPar.n_cvt);		
		gsl_matrix *G=gsl_matrix_alloc (Y->size1, Y->size1);
		gsl_matrix *U=gsl_matrix_alloc (Y->size1, Y->size1); 
		gsl_matrix *UtW=gsl_matrix_alloc (Y->size1, W->size2);
		gsl_matrix *UtY=gsl_matrix_alloc (Y->size1, Y->size2);
		gsl_vector *eval=gsl_vector_alloc (Y->size1);
		
		gsl_matrix *Y_full=gsl_matrix_alloc (cPar.ni_cvt, cPar.n_ph);
		gsl_matrix *W_full=gsl_matrix_alloc (Y_full->size1, cPar.n_cvt);
		//set covariates matrix W and phenotype matrix Y
		//an intercept should be included in W, 
		cPar.CopyCvtPhen (W, Y, 0);
		cPar.CopyCvtPhen (W_full, Y_full, 1);
				
		gsl_matrix *Y_hat=gsl_matrix_alloc (Y_full->size1, cPar.n_ph);		
		gsl_matrix *G_full=gsl_matrix_alloc (Y_full->size1, Y_full->size1);		
		gsl_matrix *H_full=gsl_matrix_alloc (Y_full->size1*Y_hat->size2, Y_full->size1*Y_hat->size2);
				
		//read relatedness matrix G, and matrix G_full
		ReadFile_kin (cPar.file_kin, cPar.indicator_idv, cPar.mapID2num, cPar.k_mode, cPar.error, G);
		if (cPar.error==true) {cout<<"error! fail to read kinship/relatedness file. "<<endl; return;}
		ReadFile_kin (cPar.file_kin, cPar.indicator_cvt, cPar.mapID2num, cPar.k_mode, cPar.error, G_full);
		if (cPar.error==true) {cout<<"error! fail to read kinship/relatedness file. "<<endl; return;}
				
		//center matrix G
		CenterMatrix (G);
		CenterMatrix (G_full);
		
		//eigen-decomposition and calculate trace_G
		cout<<"Start Eigen-Decomposition..."<<endl;
		time_start=clock();	
		cPar.trace_G=EigenDecomp (G, U, eval, 0);
		cPar.trace_G=0.0;
		for (size_t i=0; i<eval->size; i++) {
			if (gsl_vector_get (eval, i)<1e-10) {gsl_vector_set (eval, i, 0);}
			cPar.trace_G+=gsl_vector_get (eval, i);
		}
		cPar.trace_G/=(double)eval->size;
		cPar.time_eigen=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);	
		
		//calculate UtW and Uty
		CalcUtX (U, W, UtW);
		CalcUtX (U, Y, UtY);

		//calculate variance component and beta estimates
		//and then obtain predicted values
		if (cPar.n_ph==1) {
			gsl_vector *beta=gsl_vector_alloc (W->size2);
			gsl_vector *se_beta=gsl_vector_alloc (W->size2);
			
			double lambda, logl, vg, ve;
			gsl_vector_view UtY_col=gsl_matrix_column (UtY, 0);

			//obtain estimates
			CalcLambda ('R', eval, UtW, &UtY_col.vector, cPar.l_min, cPar.l_max, cPar.n_region, lambda, logl);
			CalcLmmVgVeBeta (eval, UtW, &UtY_col.vector, lambda, vg, ve, beta, se_beta);

			cout<<"REMLE estimate for vg in the null model = "<<vg<<endl;
			cout<<"REMLE estimate for ve in the null model = "<<ve<<endl;
			cPar.vg_remle_null=vg; cPar.ve_remle_null=ve;
			
			//obtain Y_hat from fixed effects
			gsl_vector_view Yhat_col=gsl_matrix_column (Y_hat, 0);			
			gsl_blas_dgemv (CblasNoTrans, 1.0, W_full, beta, 0.0, &Yhat_col.vector);
			
			//obtain H
			gsl_matrix_set_identity (H_full);
			gsl_matrix_scale (H_full, ve);
			gsl_matrix_scale (G_full, vg);
			gsl_matrix_add (H_full, G_full);
			
			//free matrices			
			gsl_vector_free(beta);
			gsl_vector_free(se_beta);
		} else {			
			gsl_matrix *Vg=gsl_matrix_alloc (cPar.n_ph, cPar.n_ph);
			gsl_matrix *Ve=gsl_matrix_alloc (cPar.n_ph, cPar.n_ph);
			gsl_matrix *B=gsl_matrix_alloc (cPar.n_ph, W->size2);
			gsl_matrix *se_B=gsl_matrix_alloc (cPar.n_ph, W->size2);
			
			//obtain estimates
			CalcMvLmmVgVeBeta (eval, UtW, UtY, cPar.em_iter, cPar.nr_iter, cPar.em_prec, cPar.nr_prec, cPar.l_min, cPar.l_max, cPar.n_region, Vg, Ve, B, se_B);
			
			cout<<"REMLE estimate for Vg in the null model: "<<endl;
			for (size_t i=0; i<Vg->size1; i++) {
				for (size_t j=0; j<=i; j++) {
					cout<<gsl_matrix_get(Vg, i, j)<<"\t";
				}
				cout<<endl;
			}
			cout<<"REMLE estimate for Ve in the null model: "<<endl;
			for (size_t i=0; i<Ve->size1; i++) {
				for (size_t j=0; j<=i; j++) {
					cout<<gsl_matrix_get(Ve, i, j)<<"\t";
				}
				cout<<endl;
			}
			cPar.Vg_remle_null.clear();
			cPar.Ve_remle_null.clear();
			for (size_t i=0; i<Vg->size1; i++) {
				for (size_t j=i; j<Vg->size2; j++) {
					cPar.Vg_remle_null.push_back(gsl_matrix_get (Vg, i, j) );
					cPar.Ve_remle_null.push_back(gsl_matrix_get (Ve, i, j) );
				}
			}
			
			//obtain Y_hat from fixed effects
			gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, W_full, B, 0.0, Y_hat);
			
			//obtain H
			KroneckerSym(G_full, Vg, H_full);
			for (size_t i=0; i<G_full->size1; i++) {
				gsl_matrix_view H_sub=gsl_matrix_submatrix (H_full, i*Ve->size1, i*Ve->size2, Ve->size1, Ve->size2);
				gsl_matrix_add (&H_sub.matrix, Ve);
			}
			
			//free matrices					
			gsl_matrix_free (Vg);
			gsl_matrix_free (Ve);
			gsl_matrix_free (B);
			gsl_matrix_free (se_B);
		}
					
		PRDT cPRDT;
		
		cPRDT.CopyFromParam(cPar);
		
		cout<<"Predicting Missing Phentypes ... "<<endl;
		time_start=clock();	
		cPRDT.MvnormPrdt(Y_hat, H_full, Y_full);
		cPar.time_opt=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);	

		cPRDT.WriteFiles(Y_full);
		
		gsl_matrix_free(Y);
		gsl_matrix_free(W);		
		gsl_matrix_free(G);
		gsl_matrix_free(U); 
		gsl_matrix_free(UtW);
		gsl_matrix_free(UtY);
		gsl_vector_free(eval);
		
		gsl_matrix_free(Y_full);
		gsl_matrix_free(Y_hat);
		gsl_matrix_free(W_full);
		gsl_matrix_free(G_full);		
		gsl_matrix_free(H_full);
	}
	
	
	//Generate Kinship matrix
	if (cPar.a_mode==21 || cPar.a_mode==22) {  
		cout<<"Calculating Relatedness Matrix ... "<<endl;
		
		gsl_matrix *G=gsl_matrix_alloc (cPar.ni_total, cPar.ni_total);
		
		time_start=clock();
		cPar.CalcKin (G);
		cPar.time_G=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		if (cPar.error==true) {cout<<"error! fail to calculate relatedness matrix. "<<endl; return;}
		
		if (cPar.a_mode==21) {
			cPar.WriteMatrix (G, "cXX");
		} else {
			cPar.WriteMatrix (G, "sXX");
		}
		
		gsl_matrix_free (G);
	}
	
	
	//LM
	if (cPar.a_mode==51 || cPar.a_mode==52 || cPar.a_mode==53 || cPar.a_mode==54) {  //Fit LM
		gsl_matrix *Y=gsl_matrix_alloc (cPar.ni_test, cPar.n_ph);
		gsl_matrix *W=gsl_matrix_alloc (Y->size1, cPar.n_cvt);	
		
		//set covariates matrix W and phenotype matrix Y		
		//an intercept should be included in W, 
		cPar.CopyCvtPhen (W, Y, 0);
		
		//Fit LM or mvLM
		if (cPar.n_ph==1) {			
			LM cLm;
			cLm.CopyFromParam(cPar);
			
			gsl_vector_view Y_col=gsl_matrix_column (Y, 0);
			
			if (!cPar.file_gene.empty()) {		
				cLm.AnalyzeGene (W, &Y_col.vector); //y is the predictor, not the phenotype
			} else if (!cPar.file_bfile.empty()) {
				cLm.AnalyzePlink (W, &Y_col.vector);
			} else {
				cLm.AnalyzeBimbam (W, &Y_col.vector);
			}
			
			cLm.WriteFiles();
			cLm.CopyToParam(cPar);
		}
		/*
		else {			 
			MVLM cMvlm;
			cMvlm.CopyFromParam(cPar);			
			
			if (!cPar.file_bfile.empty()) {
				cMvlm.AnalyzePlink (W, Y);
			} else {
				cMvlm.AnalyzeBimbam (W, Y);
			}
			
			cMvlm.WriteFiles();
			cMvlm.CopyToParam(cPar);
		}
		*/
		//release all matrices and vectors
		gsl_matrix_free (Y);
		gsl_matrix_free (W);
	} 


	//VC estimation with one or multiple kinship matrices
	//REML approach only
	//if file_kin or file_ku/kd is provided, then a_mode is changed to 5 already, in param.cpp
	//for one phenotype only; 
	if (cPar.a_mode==61) {
		gsl_matrix *Y=gsl_matrix_alloc (cPar.ni_test, cPar.n_ph);
		gsl_matrix *W=gsl_matrix_alloc (Y->size1, cPar.n_cvt);
		gsl_matrix *G=gsl_matrix_alloc (Y->size1, Y->size1*cPar.n_vc );

		//set covariates matrix W and phenotype matrix Y		
		//an intercept should be included in W, 
		cPar.CopyCvtPhen (W, Y, 0);

		//read kinship matrices
		if (!(cPar.file_mk).empty()) {
		  ReadFile_mk (cPar.file_mk, cPar.indicator_idv, cPar.mapID2num, cPar.k_mode, cPar.error, G);
		  if (cPar.error==true) {cout<<"error! fail to read kinship/relatedness file. "<<endl; return;}
	
		  //center matrix G, and obtain v_traceG
		  double d=0;
		  (cPar.v_traceG).clear();
		  for (size_t i=0; i<cPar.n_vc; i++) {
		    gsl_matrix_view G_sub=gsl_matrix_submatrix (G, 0, i*G->size1, G->size1, G->size1);
		    CenterMatrix (&G_sub.matrix);
		    d=0;
		    for (size_t j=0; j<G->size1; j++) {
		      d+=gsl_matrix_get (&G_sub.matrix, j, j);
		    }
		    d/=(double)G->size1;
		    (cPar.v_traceG).push_back(d);
		  }
		} else if (!(cPar.file_kin).empty()) {
			ReadFile_kin (cPar.file_kin, cPar.indicator_idv, cPar.mapID2num, cPar.k_mode, cPar.error, G);
			if (cPar.error==true) {cout<<"error! fail to read kinship/relatedness file. "<<endl; return;}
						
			//center matrix G
			CenterMatrix (G);

			(cPar.v_traceG).clear();
			double d=0;
			for (size_t j=0; j<G->size1; j++) {
			  d+=gsl_matrix_get (G, j, j);
			}
			d/=(double)G->size1;
			(cPar.v_traceG).push_back(d);
		}
			/*
			//eigen-decomposition and calculate trace_G
			cout<<"Start Eigen-Decomposition..."<<endl;
			time_start=clock();	
	
			if (cPar.a_mode==31) {
				cPar.trace_G=EigenDecomp (G, U, eval, 1);
			} else {
				cPar.trace_G=EigenDecomp (G, U, eval, 0);
			}

			cPar.trace_G=0.0;
			for (size_t i=0; i<eval->size; i++) {
				if (gsl_vector_get (eval, i)<1e-10) {gsl_vector_set (eval, i, 0);}
				cPar.trace_G+=gsl_vector_get (eval, i);
			}
			cPar.trace_G/=(double)eval->size;

			cPar.time_eigen=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);	
		} else {
			ReadFile_eigenU (cPar.file_ku, cPar.error, U);
			if (cPar.error==true) {cout<<"error! fail to read the U file. "<<endl; return;}
			
			ReadFile_eigenD (cPar.file_kd, cPar.error, eval);			
			if (cPar.error==true) {cout<<"error! fail to read the D file. "<<endl; return;}
			
			cPar.trace_G=0.0;
			for (size_t i=0; i<eval->size; i++) {
				if (gsl_vector_get(eval, i)<1e-10) {gsl_vector_set(eval, i, 0);}
			  	cPar.trace_G+=gsl_vector_get(eval, i);
			}
			cPar.trace_G/=(double)eval->size;
		}
		*/
		//fit multiple variance components
		if (cPar.n_ph==1) {
		  //		  if (cPar.n_vc==1) {
		    /*
		    //calculate UtW and Uty	
		    CalcUtX (U, W, UtW);
		    CalcUtX (U, Y, UtY);

		    gsl_vector_view beta=gsl_matrix_row (B, 0);
		    gsl_vector_view se_beta=gsl_matrix_row (se_B, 0);
		    gsl_vector_view UtY_col=gsl_matrix_column (UtY, 0);

		    CalcLambda ('L', eval, UtW, &UtY_col.vector, cPar.l_min, cPar.l_max, cPar.n_region, cPar.l_mle_null, cPar.logl_mle_H0);
		    CalcLmmVgVeBeta (eval, UtW, &UtY_col.vector, cPar.l_mle_null, cPar.vg_mle_null, cPar.ve_mle_null, &beta.vector, &se_beta.vector);

		    cPar.beta_mle_null.clear();
		    cPar.se_beta_mle_null.clear();
		    for (size_t i=0; i<B->size2; i++) {
		      cPar.beta_mle_null.push_back(gsl_matrix_get(B, 0, i) );
		      cPar.se_beta_mle_null.push_back(gsl_matrix_get(se_B, 0, i) );
		    }

		    CalcLambda ('R', eval, UtW, &UtY_col.vector, cPar.l_min, cPar.l_max, cPar.n_region, cPar.l_remle_null, cPar.logl_remle_H0);
		    CalcLmmVgVeBeta (eval, UtW, &UtY_col.vector, cPar.l_remle_null, cPar.vg_remle_null, cPar.ve_remle_null, &beta.vector, &se_beta.vector);
		    cPar.beta_remle_null.clear();
		    cPar.se_beta_remle_null.clear();
		    for (size_t i=0; i<B->size2; i++) {
		      cPar.beta_remle_null.push_back(gsl_matrix_get(B, 0, i) );
		      cPar.se_beta_remle_null.push_back(gsl_matrix_get(se_B, 0, i) );
		    }
				
		    CalcPve (eval, UtW, &UtY_col.vector, cPar.l_remle_null, cPar.trace_G, cPar.pve_null, cPar.pve_se_null);
		    cPar.PrintSummary();
				
		    //calculate and output residuals
		    if (cPar.a_mode==5) {
		      gsl_vector *Utu_hat=gsl_vector_alloc (Y->size1);
		      gsl_vector *Ute_hat=gsl_vector_alloc (Y->size1);
		      gsl_vector *u_hat=gsl_vector_alloc (Y->size1);
		      gsl_vector *e_hat=gsl_vector_alloc (Y->size1);
		      gsl_vector *y_hat=gsl_vector_alloc (Y->size1);
					
		      //obtain Utu and Ute
		      gsl_vector_memcpy (y_hat, &UtY_col.vector);
		      gsl_blas_dgemv (CblasNoTrans, -1.0, UtW, &beta.vector, 1.0, y_hat);
		      
		      double d, u, e;
		      for (size_t i=0; i<eval->size; i++) {
			d=gsl_vector_get (eval, i);
			u=cPar.l_remle_null*d/(cPar.l_remle_null*d+1.0)*gsl_vector_get(y_hat, i);
			e=1.0/(cPar.l_remle_null*d+1.0)*gsl_vector_get(y_hat, i);
			gsl_vector_set (Utu_hat, i, u);
			gsl_vector_set (Ute_hat, i, e);
		      }
					
		      //obtain u and e
		      gsl_blas_dgemv (CblasNoTrans, 1.0, U, Utu_hat, 0.0, u_hat);
		      gsl_blas_dgemv (CblasNoTrans, 1.0, U, Ute_hat, 0.0, e_hat);
		      
		      //output residuals					
		      cPar.WriteVector(u_hat, "residU");
		      cPar.WriteVector(e_hat, "residE");
		      
		      gsl_vector_free(u_hat);
		      gsl_vector_free(e_hat);
		      gsl_vector_free(y_hat);
		    }	
*/	
		  //		  } else {
		    gsl_vector_view Y_col=gsl_matrix_column (Y, 0);
		    VC cVc;
		    cVc.CopyFromParam(cPar); 
		    cVc.CalcVCreml (G, W, &Y_col.vector);			
		    cVc.CopyToParam(cPar);

		    //obtain pve from sigma2
		    //obtain se_pve from se_sigma2
		    
		    //}
		} 

		
	}
	
	
	//LMM or mvLMM or Eigen-Decomposition
	if (cPar.a_mode==1 || cPar.a_mode==2 || cPar.a_mode==3 || cPar.a_mode==4 || cPar.a_mode==5 || cPar.a_mode==31) {  //Fit LMM or mvLMM or eigen
		gsl_matrix *Y=gsl_matrix_alloc (cPar.ni_test, cPar.n_ph);
		gsl_matrix *W=gsl_matrix_alloc (Y->size1, cPar.n_cvt);
		gsl_matrix *B=gsl_matrix_alloc (Y->size2, W->size2);	//B is a d by c matrix
		gsl_matrix *se_B=gsl_matrix_alloc (Y->size2, W->size2);
		gsl_matrix *G=gsl_matrix_alloc (Y->size1, Y->size1);
		gsl_matrix *U=gsl_matrix_alloc (Y->size1, Y->size1); 
		gsl_matrix *UtW=gsl_matrix_alloc (Y->size1, W->size2);
		gsl_matrix *UtY=gsl_matrix_alloc (Y->size1, Y->size2);
		gsl_vector *eval=gsl_vector_alloc (Y->size1);
				
		//set covariates matrix W and phenotype matrix Y		
		//an intercept should be included in W, 
		cPar.CopyCvtPhen (W, Y, 0);
		
		//read relatedness matrix G	
		if (!(cPar.file_kin).empty()) {
			ReadFile_kin (cPar.file_kin, cPar.indicator_idv, cPar.mapID2num, cPar.k_mode, cPar.error, G);
			if (cPar.error==true) {cout<<"error! fail to read kinship/relatedness file. "<<endl; return;}
						
			//center matrix G
			CenterMatrix (G);
			
			//eigen-decomposition and calculate trace_G
			cout<<"Start Eigen-Decomposition..."<<endl;
			time_start=clock();	
	
			if (cPar.a_mode==31) {
				cPar.trace_G=EigenDecomp (G, U, eval, 1);
			} else {
				cPar.trace_G=EigenDecomp (G, U, eval, 0);
			}

			cPar.trace_G=0.0;
			for (size_t i=0; i<eval->size; i++) {
				if (gsl_vector_get (eval, i)<1e-10) {gsl_vector_set (eval, i, 0);}
				cPar.trace_G+=gsl_vector_get (eval, i);
			}
			cPar.trace_G/=(double)eval->size;

			cPar.time_eigen=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);	
		} else {
			ReadFile_eigenU (cPar.file_ku, cPar.error, U);
			if (cPar.error==true) {cout<<"error! fail to read the U file. "<<endl; return;}
			
			ReadFile_eigenD (cPar.file_kd, cPar.error, eval);			
			if (cPar.error==true) {cout<<"error! fail to read the D file. "<<endl; return;}
			
			cPar.trace_G=0.0;
			for (size_t i=0; i<eval->size; i++) {
				if (gsl_vector_get(eval, i)<1e-10) {gsl_vector_set(eval, i, 0);}
			  	cPar.trace_G+=gsl_vector_get(eval, i);
			}
			cPar.trace_G/=(double)eval->size;
		}
		
		if (cPar.a_mode==31) {
			cPar.WriteMatrix(U, "eigenU");
			cPar.WriteVector(eval, "eigenD");
		} else {
			//calculate UtW and Uty	
			CalcUtX (U, W, UtW);
			CalcUtX (U, Y, UtY);			

			//calculate REMLE/MLE estimate and pve for univariate model
			if (cPar.n_ph==1) {
				gsl_vector_view beta=gsl_matrix_row (B, 0);
				gsl_vector_view se_beta=gsl_matrix_row (se_B, 0);
				gsl_vector_view UtY_col=gsl_matrix_column (UtY, 0);

				CalcLambda ('L', eval, UtW, &UtY_col.vector, cPar.l_min, cPar.l_max, cPar.n_region, cPar.l_mle_null, cPar.logl_mle_H0);
				CalcLmmVgVeBeta (eval, UtW, &UtY_col.vector, cPar.l_mle_null, cPar.vg_mle_null, cPar.ve_mle_null, &beta.vector, &se_beta.vector);

				cPar.beta_mle_null.clear();
				cPar.se_beta_mle_null.clear();
				for (size_t i=0; i<B->size2; i++) {
					cPar.beta_mle_null.push_back(gsl_matrix_get(B, 0, i) );
					cPar.se_beta_mle_null.push_back(gsl_matrix_get(se_B, 0, i) );
				}

				CalcLambda ('R', eval, UtW, &UtY_col.vector, cPar.l_min, cPar.l_max, cPar.n_region, cPar.l_remle_null, cPar.logl_remle_H0);
				CalcLmmVgVeBeta (eval, UtW, &UtY_col.vector, cPar.l_remle_null, cPar.vg_remle_null, cPar.ve_remle_null, &beta.vector, &se_beta.vector);
				cPar.beta_remle_null.clear();
				cPar.se_beta_remle_null.clear();
				for (size_t i=0; i<B->size2; i++) {
					cPar.beta_remle_null.push_back(gsl_matrix_get(B, 0, i) );
					cPar.se_beta_remle_null.push_back(gsl_matrix_get(se_B, 0, i) );
				}
				
				CalcPve (eval, UtW, &UtY_col.vector, cPar.l_remle_null, cPar.trace_G, cPar.pve_null, cPar.pve_se_null);
				cPar.PrintSummary();
				
				//calculate and output residuals
				if (cPar.a_mode==5) {
					gsl_vector *Utu_hat=gsl_vector_alloc (Y->size1);
					gsl_vector *Ute_hat=gsl_vector_alloc (Y->size1);
					gsl_vector *u_hat=gsl_vector_alloc (Y->size1);
					gsl_vector *e_hat=gsl_vector_alloc (Y->size1);
					gsl_vector *y_hat=gsl_vector_alloc (Y->size1);
					
					//obtain Utu and Ute
					gsl_vector_memcpy (y_hat, &UtY_col.vector);
					gsl_blas_dgemv (CblasNoTrans, -1.0, UtW, &beta.vector, 1.0, y_hat);
					
					double d, u, e;
					for (size_t i=0; i<eval->size; i++) {
						d=gsl_vector_get (eval, i);
						u=cPar.l_remle_null*d/(cPar.l_remle_null*d+1.0)*gsl_vector_get(y_hat, i);
						e=1.0/(cPar.l_remle_null*d+1.0)*gsl_vector_get(y_hat, i);
						gsl_vector_set (Utu_hat, i, u);
						gsl_vector_set (Ute_hat, i, e);
					}
					
					//obtain u and e
					gsl_blas_dgemv (CblasNoTrans, 1.0, U, Utu_hat, 0.0, u_hat);
					gsl_blas_dgemv (CblasNoTrans, 1.0, U, Ute_hat, 0.0, e_hat);
					
					//output residuals					
					cPar.WriteVector(u_hat, "residU");
					cPar.WriteVector(e_hat, "residE");
					
					gsl_vector_free(u_hat);
					gsl_vector_free(e_hat);
					gsl_vector_free(y_hat);
				}							
			} 
			
			//Fit LMM or mvLMM
			if (cPar.a_mode==1 || cPar.a_mode==2 || cPar.a_mode==3 || cPar.a_mode==4) {
				if (cPar.n_ph==1) {			
					LMM cLmm;
					cLmm.CopyFromParam(cPar);
					
					gsl_vector_view Y_col=gsl_matrix_column (Y, 0);
					gsl_vector_view UtY_col=gsl_matrix_column (UtY, 0);
					
					if (!cPar.file_gene.empty()) 		
						cLmm.AnalyzeGene (U, eval, UtW, &UtY_col.vector, W, &Y_col.vector); //y is the predictor, not the phenotype
					
					else if (!cPar.file_bfile.empty()) 
						cLmm.AnalyzePlink (U, eval, UtW, &UtY_col.vector, W, &Y_col.vector);
							
						
						 // WJA added
				       	else if(!cPar.file_bgenfile.empty()) 
						cLmm.Analyzebgen (U, eval, UtW, &UtY_col.vector, W, &Y_col.vector);	
													 		
				   	else 
						cLmm.AnalyzeBimbam (U, eval, UtW, &UtY_col.vector, W, &Y_col.vector);
											
					cLmm.WriteFiles();
					cLmm.CopyToParam(cPar);
				} else {			 
					MVLMM cMvlmm;
					cMvlmm.CopyFromParam(cPar);			
					
					if (!cPar.file_bfile.empty()) {
						cMvlmm.AnalyzePlink (U, eval, UtW, UtY);
					} else {
						cMvlmm.AnalyzeBimbam (U, eval, UtW, UtY);
					}
					
					cMvlmm.WriteFiles();
					cMvlmm.CopyToParam(cPar);
				}
			}
		}
		
				
		//release all matrices and vectors
		gsl_matrix_free (Y);
		gsl_matrix_free (W);
		gsl_matrix_free(B);
		gsl_matrix_free(se_B);
		gsl_matrix_free (G);	
		gsl_matrix_free (U);
		gsl_matrix_free (UtW);
		gsl_matrix_free (UtY);
		gsl_vector_free (eval);
	} 
	
	
	//BSLMM
	if (cPar.a_mode==11 || cPar.a_mode==12 || cPar.a_mode==13) {
		gsl_vector *y=gsl_vector_alloc (cPar.ni_test);
		gsl_matrix *W=gsl_matrix_alloc (y->size, cPar.n_cvt);	
		gsl_matrix *G=gsl_matrix_alloc (y->size, y->size);
		gsl_matrix *UtX=gsl_matrix_alloc (y->size, cPar.ns_test);	
		
		//set covariates matrix W and phenotype vector y		
		//an intercept should be included in W, 
		cPar.CopyCvtPhen (W, y, 0);
		
		//center y, even for case/control data
		cPar.pheno_mean=CenterVector(y);

		//run bslmm if rho==1
		if (cPar.rho_min==1 && cPar.rho_max==1) {
		  //read genotypes X (not UtX)
		  cPar.ReadGenotypes (UtX, G, false);

		  //perform BSLMM analysis
		  BSLMM cBslmm;
		  cBslmm.CopyFromParam(cPar);
		  time_start=clock();	
		  cBslmm.MCMC(UtX, y);
		  cPar.time_opt=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		  cBslmm.CopyToParam(cPar);
		  //else, if rho!=1
		} else {
		gsl_matrix *U=gsl_matrix_alloc (y->size, y->size); 
		gsl_vector *eval=gsl_vector_alloc (y->size);
		gsl_matrix *UtW=gsl_matrix_alloc (y->size, W->size2);
		gsl_vector *Uty=gsl_vector_alloc (y->size);

		
		//read relatedness matrix G		
		if (!(cPar.file_kin).empty()) {		
			cPar.ReadGenotypes (UtX, G, false);
			
			//read relatedness matrix G
			ReadFile_kin (cPar.file_kin, cPar.indicator_idv, cPar.mapID2num, cPar.k_mode, cPar.error, G);
			if (cPar.error==true) {cout<<"error! fail to read kinship/relatedness file. "<<endl; return;}
			
			//center matrix G
			CenterMatrix (G);
		} else {
			cPar.ReadGenotypes (UtX, G, true);
		}
		
		//eigen-decomposition and calculate trace_G
		cout<<"Start Eigen-Decomposition..."<<endl;
		time_start=clock();
		cPar.trace_G=EigenDecomp (G, U, eval, 0);
		cPar.trace_G=0.0;
		for (size_t i=0; i<eval->size; i++) {
			if (gsl_vector_get (eval, i)<1e-10) {gsl_vector_set (eval, i, 0);}
			cPar.trace_G+=gsl_vector_get (eval, i);
		}
		cPar.trace_G/=(double)eval->size;
		cPar.time_eigen=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);			
		
		//calculate UtW and Uty		
		CalcUtX (U, W, UtW);
		CalcUtX (U, y, Uty);
		
		//calculate REMLE/MLE estimate and pve
		CalcLambda ('L', eval, UtW, Uty, cPar.l_min, cPar.l_max, cPar.n_region, cPar.l_mle_null, cPar.logl_mle_H0);
		CalcLambda ('R', eval, UtW, Uty, cPar.l_min, cPar.l_max, cPar.n_region, cPar.l_remle_null, cPar.logl_remle_H0);
		CalcPve (eval, UtW, Uty, cPar.l_remle_null, cPar.trace_G, cPar.pve_null, cPar.pve_se_null);
		
		cPar.PrintSummary();
				
		//Creat and calcualte UtX, use a large memory
		cout<<"Calculating UtX..."<<endl;
		time_start=clock();							
		CalcUtX (U, UtX);
		cPar.time_UtX=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		//perform BSLMM analysis
		BSLMM cBslmm;
		cBslmm.CopyFromParam(cPar);
		time_start=clock();	
		if (cPar.a_mode==12) {  //ridge regression				
			cBslmm.RidgeR(U, UtX, Uty, eval, cPar.l_remle_null);
		} else {	//Run MCMC
			cBslmm.MCMC(U, UtX, Uty, eval, y);
		}
		cPar.time_opt=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		cBslmm.CopyToParam(cPar);
		
		//release all matrices and vectors
		gsl_matrix_free (G);	
		gsl_matrix_free (U);
		gsl_matrix_free (UtW);
		gsl_vector_free (eval);
		gsl_vector_free (Uty);

		}
		gsl_matrix_free (W);
		gsl_vector_free (y);
		gsl_matrix_free (UtX);
	} 
	
	
		
	cPar.time_total=(clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0);
	
	return;
}




void GEMMA::WriteLog (int argc, char ** argv, PARAM &cPar) 
{
	string file_str;
	file_str=cPar.path_out+"/"+cPar.file_out;
	file_str+=".log.txt";
	
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing log file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"##"<<endl;
	outfile<<"## GEMMA Version = "<<version<<endl;
	
	outfile<<"##"<<endl;
	outfile<<"## Command Line Input = ";
	for(int i = 1; i < argc; i++) {	
		outfile<<argv[i]<<" ";
	}
	outfile<<endl;

	outfile<<"##"<<endl;
	time_t  rawtime; 
	time(&rawtime);
	tm *ptm = localtime (&rawtime);

	outfile<<"## Date = "<<asctime(ptm)<<endl;
	  //ptm->tm_year<<":"<<ptm->tm_month<<":"<<ptm->tm_day":"<<ptm->tm_hour<<":"<<ptm->tm_min<<endl;
	
	outfile<<"##"<<endl;
	outfile<<"## Summary Statistics:"<<endl;
	outfile<<"## number of total individuals = "<<cPar.ni_total<<endl;	
	if (cPar.a_mode==43) {
		outfile<<"## number of analyzed individuals = "<<cPar.ni_cvt<<endl;
		outfile<<"## number of individuals with full phenotypes = "<<cPar.ni_test<<endl;
	} else {
		outfile<<"## number of analyzed individuals = "<<cPar.ni_test<<endl;
	}
	outfile<<"## number of covariates = "<<cPar.n_cvt<<endl;
	outfile<<"## number of phenotypes = "<<cPar.n_ph<<endl;
	if (cPar.a_mode==43) {
		outfile<<"## number of observed data = "<<cPar.np_obs<<endl;
		outfile<<"## number of missing data = "<<cPar.np_miss<<endl;
	}
	if (cPar.a_mode==61) {
		outfile<<"## number of variance components = "<<cPar.n_vc<<endl;
	}
		
	if (!(cPar.file_gene).empty()) {
		outfile<<"## number of total genes = "<<cPar.ng_total<<endl;
		outfile<<"## number of analyzed genes = "<<cPar.ng_test<<endl;		
	} else if (cPar.file_epm.empty()) {	
		outfile<<"## number of total SNPs = "<<cPar.ns_total<<endl;	
		outfile<<"## number of analyzed SNPs = "<<cPar.ns_test<<endl;
	} else {
		outfile<<"## number of analyzed SNPs = "<<cPar.ns_test<<endl;
	}
	
	if (cPar.a_mode==13) {
		outfile<<"## number of cases = "<<cPar.ni_case<<endl;
		outfile<<"## number of controls = "<<cPar.ni_control<<endl;
	}


	if (cPar.a_mode==61) {
	  //	        outfile<<"## REMLE log-likelihood in the null model = "<<cPar.logl_remle_H0<<endl;
		if (cPar.n_ph==1) {
		  outfile<<"## pve estimate in the null model = ";
		  for (size_t i=0; i<cPar.v_pve.size(); i++) {
		    outfile<<"  "<<cPar.v_pve[i];
		  }
		  outfile<<endl;

		  outfile<<"## se(pve) in the null model = ";
		  for (size_t i=0; i<cPar.v_se_pve.size(); i++) {
		    outfile<<"  "<<cPar.v_se_pve[i];
		  }
		  outfile<<endl;

		  outfile<<"## sigma2 estimate in the null model = ";
		  for (size_t i=0; i<cPar.v_sigma2.size(); i++) {
		    outfile<<"  "<<cPar.v_sigma2[i];
		  }
		  outfile<<endl;

		  outfile<<"## se(sigma2) in the null model = ";
		  for (size_t i=0; i<cPar.v_se_sigma2.size(); i++) {
		    outfile<<"  "<<cPar.v_se_sigma2[i];
		  }
		  outfile<<endl;
		  /*
			outfile<<"## beta estimate in the null model = ";
			for (size_t i=0; i<cPar.beta_remle_null.size(); i++) {
				outfile<<"  "<<cPar.beta_remle_null[i];
			}
			outfile<<endl;
			outfile<<"## se(beta) = ";
			for (size_t i=0; i<cPar.se_beta_remle_null.size(); i++) {
				outfile<<"  "<<cPar.se_beta_remle_null[i];
			}
			outfile<<endl;
		  */
		}
	}
	
	if (cPar.a_mode==1 || cPar.a_mode==2 || cPar.a_mode==3 || cPar.a_mode==4 || cPar.a_mode==5 || cPar.a_mode==11 || cPar.a_mode==12 || cPar.a_mode==13) {
		outfile<<"## REMLE log-likelihood in the null model = "<<cPar.logl_remle_H0<<endl;
		outfile<<"## MLE log-likelihood in the null model = "<<cPar.logl_mle_H0<<endl;
		if (cPar.n_ph==1) {
			//outfile<<"## lambda REMLE estimate in the null (linear mixed) model = "<<cPar.l_remle_null<<endl;
			//outfile<<"## lambda MLE estimate in the null (linear mixed) model = "<<cPar.l_mle_null<<endl;	
			outfile<<"## pve estimate in the null model = "<<cPar.pve_null<<endl;
			outfile<<"## se(pve) in the null model = "<<cPar.pve_se_null<<endl;	
			outfile<<"## vg estimate in the null model = "<<cPar.vg_remle_null<<endl;
			outfile<<"## ve estimate in the null model = "<<cPar.ve_remle_null<<endl;	
			outfile<<"## beta estimate in the null model = ";
			for (size_t i=0; i<cPar.beta_remle_null.size(); i++) {
				outfile<<"  "<<cPar.beta_remle_null[i];
			}
			outfile<<endl;
			outfile<<"## se(beta) = ";
			for (size_t i=0; i<cPar.se_beta_remle_null.size(); i++) {
				outfile<<"  "<<cPar.se_beta_remle_null[i];
			}
			outfile<<endl;
			
		} else {
			size_t c;
			outfile<<"## REMLE estimate for Vg in the null model: "<<endl;			
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<cPar.Vg_remle_null[c]<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## se(Vg): "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<sqrt(cPar.VVg_remle_null[c])<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## REMLE estimate for Ve in the null model: "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<cPar.Ve_remle_null[c]<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## se(Ve): "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<sqrt(cPar.VVe_remle_null[c])<<"\t";
				}
				outfile<<endl;
			}
			
			outfile<<"## MLE estimate for Vg in the null model: "<<endl;
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<cPar.n_ph; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<cPar.Vg_mle_null[c]<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## se(Vg): "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<sqrt(cPar.VVg_mle_null[c])<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## MLE estimate for Ve in the null model: "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<cPar.n_ph; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<cPar.Ve_mle_null[c]<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## se(Ve): "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<sqrt(cPar.VVe_mle_null[c])<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## estimate for B (d by c) in the null model (columns correspond to the covariates provided in the file): "<<endl;
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<cPar.n_cvt; j++) {
					c=i*cPar.n_cvt+j;
					outfile<<cPar.beta_remle_null[c]<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## se(B): "<<endl;
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<cPar.n_cvt; j++) {
					c=i*cPar.n_cvt+j;
					outfile<<cPar.se_beta_remle_null[c]<<"\t";
				}
				outfile<<endl;
			}
		}
	}
	
	/*
	if (cPar.a_mode==1 || cPar.a_mode==2 || cPar.a_mode==3 || cPar.a_mode==4 || cPar.a_mode==11 || cPar.a_mode==12 || cPar.a_mode==13) {
		if (cPar.n_ph==1) {
			outfile<<"## REMLE vg estimate in the null model = "<<cPar.vg_remle_null<<endl;
			outfile<<"## REMLE ve estimate in the null model = "<<cPar.ve_remle_null<<endl;	
		} else {
			size_t c;
			outfile<<"## REMLE estimate for Vg in the null model: "<<endl;			
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<cPar.Vg_remle_null[c]<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## REMLE estimate for Ve in the null model: "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<cPar.Ve_remle_null[c]<<"\t";
				}
				outfile<<endl;
			}
		}
	}
	 */
	
	
	if (cPar.a_mode==11 || cPar.a_mode==12 || cPar.a_mode==13) {
		outfile<<"## estimated mean = "<<cPar.pheno_mean<<endl;
	}
	
	if (cPar.a_mode==11 || cPar.a_mode==13) {	
		outfile<<"##"<<endl;
		outfile<<"## MCMC related:"<<endl;	
		outfile<<"## initial value of h = "<<cPar.cHyp_initial.h<<endl;
		outfile<<"## initial value of rho = "<<cPar.cHyp_initial.rho<<endl;
		outfile<<"## initial value of pi = "<<exp(cPar.cHyp_initial.logp)<<endl;
		outfile<<"## initial value of |gamma| = "<<cPar.cHyp_initial.n_gamma<<endl;
		outfile<<"## random seed = "<<cPar.randseed<<endl;
		outfile<<"## acceptance ratio = "<<(double)cPar.n_accept/(double)((cPar.w_step+cPar.s_step)*cPar.n_mh)<<endl;
	}
	
	outfile<<"##"<<endl;
	outfile<<"## Computation Time:"<<endl;
	outfile<<"## total computation time = "<<cPar.time_total<<" min "<<endl;
	outfile<<"## computation time break down: "<<endl;
	if (cPar.a_mode==21 || cPar.a_mode==22 || cPar.a_mode==11 || cPar.a_mode==13) {
		outfile<<"##      time on calculating relatedness matrix = "<<cPar.time_G<<" min "<<endl;
	}
	if (cPar.a_mode==31) {
		outfile<<"##      time on eigen-decomposition = "<<cPar.time_eigen<<" min "<<endl;
	}
	if (cPar.a_mode==1 || cPar.a_mode==2 || cPar.a_mode==3 || cPar.a_mode==4 || cPar.a_mode==5 || cPar.a_mode==11 || cPar.a_mode==12 || cPar.a_mode==13) {
		outfile<<"##      time on eigen-decomposition = "<<cPar.time_eigen<<" min "<<endl;
		outfile<<"##      time on calculating UtX = "<<cPar.time_UtX<<" min "<<endl;		
	}
	if ((cPar.a_mode>=1 && cPar.a_mode<=4) || (cPar.a_mode>=51 && cPar.a_mode<=54) ) {
		outfile<<"##      time on optimization = "<<cPar.time_opt<<" min "<<endl;
	}
	if (cPar.a_mode==11 || cPar.a_mode==13) {
		outfile<<"##      time on proposal = "<<cPar.time_Proposal<<" min "<<endl;
		outfile<<"##      time on mcmc = "<<cPar.time_opt<<" min "<<endl;
		outfile<<"##      time on Omega = "<<cPar.time_Omega<<" min "<<endl;
	}
	if (cPar.a_mode==41 || cPar.a_mode==42) {
		outfile<<"##      time on eigen-decomposition = "<<cPar.time_eigen<<" min "<<endl;
	}
	if (cPar.a_mode==43) {
		outfile<<"##      time on eigen-decomposition = "<<cPar.time_eigen<<" min "<<endl;
		outfile<<"##      time on predicting phenotypes = "<<cPar.time_opt<<" min "<<endl;
	}
	outfile<<"##"<<endl;
	
	outfile.close();
	outfile.clear();
	return;
}


