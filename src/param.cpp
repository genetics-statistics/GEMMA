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
#include <cmath>
#include <algorithm>

#include "gsl/gsl_randist.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"

#include "eigenlib.h"
#include "mathfunc.h"

#ifdef FORCE_FLOAT
#include "param_float.h"
#include "io_float.h"
#else
#include "param.h"
#include "io.h"
#endif

using namespace std;





PARAM::PARAM(void):
mode_silence (false), a_mode (0), k_mode(1), d_pace (100000),
file_out("result"), path_out("./output/"),
miss_level(0.05), maf_level(0.01), hwe_level(0), r2_level(0.9999),
l_min(1e-5), l_max(1e5), n_region(10),p_nr(0.001),em_prec(0.0001),nr_prec(0.0001),em_iter(10000),nr_iter(100),crt(0),
pheno_mean(0), noconstrain (false),
h_min(-1), h_max(-1),	h_scale(-1),
rho_min(0.0), rho_max(1.0),	rho_scale(-1),
logp_min(0.0), logp_max(0.0), logp_scale(-1),
s_min(0), s_max(300),
w_step(100000),	s_step(1000000),
r_pace(10), w_pace(1000),
n_accept(0),
n_mh(10),
geo_mean(2000.0),
randseed(-1),
window_cm(0), window_bp(0), window_ns(0),
error(false),
ni_subsample(0), n_cvt(1), n_vc(1),
time_total(0.0), time_G(0.0), time_eigen(0.0), time_UtX(0.0), time_UtZ(0.0), time_opt(0.0), time_Omega(0.0)
{}


//read files
//obtain ns_total, ng_total, ns_test, ni_test
void PARAM::ReadFiles (void)
{
	string file_str;


	if (!file_cat.empty()) {
	  if (ReadFile_cat (file_cat, mapRS2cat, n_vc)==false) {error=true;}
	}

	if (!file_var.empty()) {
	  if (ReadFile_var (file_var, mapRS2var)==false) {error=true;}
	}

	if (!file_mk.empty()) {
	  if (CountFileLines (file_mk, n_vc)==false) {error=true;}
	}

	if (!file_snps.empty()) {
		if (ReadFile_snps (file_snps, setSnps)==false) {error=true;}
	} else {
		setSnps.clear();
	}

	//for prediction
	if (!file_epm.empty()) {
		if (ReadFile_est (file_epm, est_column, mapRS2est)==false) {error=true;}

		if (!file_bfile.empty()) {
			file_str=file_bfile+".bim";
			if (ReadFile_bim (file_str, snpInfo)==false) {error=true;}

			file_str=file_bfile+".fam";
			if (ReadFile_fam (file_str, indicator_pheno, pheno, mapID2num, p_column)==false) {error=true;}
		}

		if (!file_geno.empty()) {
			if (ReadFile_pheno (file_pheno, indicator_pheno, pheno, p_column)==false) {error=true;}

			if (CountFileLines (file_geno, ns_total)==false) {error=true;}
		}

		if (!file_ebv.empty() ) {
			if (ReadFile_column (file_ebv, indicator_bv, vec_bv, 1)==false) {error=true;}
		}

		if (!file_log.empty() ) {
			if (ReadFile_log (file_log, pheno_mean)==false) {error=true;}
		}

		//convert indicator_pheno to indicator_idv
		int k=1;
		for (size_t i=0; i<indicator_pheno.size(); i++) {
			k=1;
			for (size_t j=0; j<indicator_pheno[i].size(); j++) {
				if (indicator_pheno[i][j]==0) {k=0;}
			}
			indicator_idv.push_back(k);
		}

		ns_test=0;

		return;
	}

	//read covariates before the genotype files
	if (!file_cvt.empty() ) {
		if (ReadFile_cvt (file_cvt, indicator_cvt, cvt, n_cvt)==false) {error=true;}

		if ((indicator_cvt).size()==0) {
			n_cvt=1;
		}
	} else {
		n_cvt=1;
	}

	if (!file_gxe.empty() ) {
	  if (ReadFile_column (file_gxe, indicator_gxe, gxe, 1)==false) {error=true;}
	}
	if (!file_weight.empty() ) {
	  if (ReadFile_column (file_weight, indicator_weight, weight, 1)==false) {error=true;}
	}

	//read genotype and phenotype file for plink format
	if (!file_bfile.empty()) {
		file_str=file_bfile+".bim";
		if (ReadFile_bim (file_str, snpInfo)==false) {error=true;}

		file_str=file_bfile+".fam";
		if (ReadFile_fam (file_str, indicator_pheno, pheno, mapID2num, p_column)==false) {error=true;}

		//post-process covariates and phenotypes, obtain ni_test, save all useful covariates
		ProcessCvtPhen();

		//obtain covariate matrix
		gsl_matrix *W=gsl_matrix_alloc (ni_test, n_cvt);
		CopyCvt (W);

		file_str=file_bfile+".bed";
		if (ReadFile_bed (file_str, setSnps, W, indicator_idv, indicator_snp, snpInfo, maf_level, miss_level, hwe_level, r2_level, ns_test)==false) {error=true;}

		gsl_matrix_free(W);

		ns_total=indicator_snp.size();
	}

	//read genotype and phenotype file for bimbam format
	if (!file_geno.empty()) {
		//annotation file before genotype file
		if (!file_anno.empty() ) {
			if (ReadFile_anno (file_anno, mapRS2chr, mapRS2bp, mapRS2cM)==false) {error=true;}
		}

		//phenotype file before genotype file
		if (ReadFile_pheno (file_pheno, indicator_pheno, pheno, p_column)==false) {error=true;}

		//post-process covariates and phenotypes, obtain ni_test, save all useful covariates
		ProcessCvtPhen();

		//obtain covariate matrix
		gsl_matrix *W=gsl_matrix_alloc (ni_test, n_cvt);
		CopyCvt (W);

		if (ReadFile_geno (file_geno, setSnps, W, indicator_idv, indicator_snp, maf_level, miss_level, hwe_level, r2_level, mapRS2chr, mapRS2bp, mapRS2cM, snpInfo, ns_test)==false) {error=true;}

		gsl_matrix_free(W);

		ns_total=indicator_snp.size();
	}

	if (!file_gene.empty()) {
		if (ReadFile_pheno (file_pheno, indicator_pheno, pheno, p_column)==false) {error=true;}

		//convert indicator_pheno to indicator_idv
		int k=1;
		for (size_t i=0; i<indicator_pheno.size(); i++) {
			k=1;
			for (size_t j=0; j<indicator_pheno[i].size(); j++) {
				if (indicator_pheno[i][j]==0) {k=0;}
			}
			indicator_idv.push_back(k);
		}

        //post-process covariates and phenotypes, obtain ni_test, save all useful covariates
        ProcessCvtPhen();

        //obtain covariate matrix
        gsl_matrix *W=gsl_matrix_alloc (ni_test, n_cvt);
        CopyCvt (W);

		if (ReadFile_gene (file_gene, vec_read, snpInfo, ng_total)==false) {error=true;}
	}


	//read is after gene file
	if (!file_read.empty() ) {
		if (ReadFile_column (file_read, indicator_read, vec_read, 1)==false) {error=true;}

		ni_test=0;
		for (vector<int>::size_type i=0; i<(indicator_idv).size(); ++i) {
			indicator_idv[i]*=indicator_read[i];
			ni_test+=indicator_idv[i];
		}

		if (ni_test==0) {
			error=true;
			cout<<"error! number of analyzed individuals equals 0. "<<endl;
			return;
		}
	}

	//for ridge prediction, read phenotype only
	if (file_geno.empty() && file_gene.empty() && !file_pheno.empty()) {
		if (ReadFile_pheno (file_pheno, indicator_pheno, pheno, p_column)==false) {error=true;}

		//post-process covariates and phenotypes, obtain ni_test, save all useful covariates
		ProcessCvtPhen();
	}

	return;
}






void PARAM::CheckParam (void)
{
	struct stat fileInfo;
	string str;

	//check parameters
	if (k_mode!=1 && k_mode!=2) {cout<<"error! unknown kinship/relatedness input mode: "<<k_mode<<endl; error=true;}
	if (a_mode!=1 && a_mode!=2 && a_mode!=3 && a_mode!=4 && a_mode!=5 && a_mode!=11 && a_mode!=12 && a_mode!=13 && a_mode!=14 && a_mode!=21 && a_mode!=22 && a_mode!=25 && a_mode!=26 && a_mode!=27 && a_mode!=28 && a_mode!=31 && a_mode!=41 && a_mode!=42 && a_mode!=43 && a_mode!=51 && a_mode!=52 && a_mode!=53 && a_mode!=54 && a_mode!=61 && a_mode!=62 && a_mode!=71)
	{cout<<"error! unknown analysis mode: "<<a_mode<<". make sure -gk or -eigen or -lmm or -bslmm -predict or -calccov is sepcified correctly."<<endl; error=true;}
	if (miss_level>1) {cout<<"error! missing level needs to be between 0 and 1. current value = "<<miss_level<<endl; error=true;}
	if (maf_level>0.5) {cout<<"error! maf level needs to be between 0 and 0.5. current value = "<<maf_level<<endl; error=true;}
	if (hwe_level>1) {cout<<"error! hwe level needs to be between 0 and 1. current value = "<<hwe_level<<endl; error=true;}
	if (r2_level>1) {cout<<"error! r2 level needs to be between 0 and 1. current value = "<<r2_level<<endl; error=true;}

	if (l_max<l_min) {cout<<"error! maximum lambda value must be larger than the minimal value. current values = "<<l_max<<" and "<<l_min<<endl; error=true;}
	if (h_max<h_min) {cout<<"error! maximum h value must be larger than the minimal value. current values = "<<h_max<<" and "<<h_min<<endl; error=true;}
	if (s_max<s_min) {cout<<"error! maximum s value must be larger than the minimal value. current values = "<<s_max<<" and "<<s_min<<endl; error=true;}
	if (rho_max<rho_min) {cout<<"error! maximum rho value must be larger than the minimal value. current values = "<<rho_max<<" and "<<rho_min<<endl; error=true;}
	if (logp_max<logp_min) {cout<<"error! maximum logp value must be larger than the minimal value. current values = "<<logp_max/log(10)<<" and "<<logp_min/log(10)<<endl; error=true;}

	if (h_max>1) {cout<<"error! h values must be bewtween 0 and 1. current values = "<<h_max<<" and "<<h_min<<endl; error=true;}
	if (rho_max>1) {cout<<"error! rho values must be between 0 and 1. current values = "<<rho_max<<" and "<<rho_min<<endl; error=true;}
	if (logp_max>0) {cout<<"error! maximum logp value must be smaller than 0. current values = "<<logp_max/log(10)<<" and "<<logp_min/log(10)<<endl; error=true;}
	if (l_max<l_min) {cout<<"error! maximum lambda value must be larger than the minimal value. current values = "<<l_max<<" and "<<l_min<<endl; error=true;}

	if (h_scale>1.0) {cout<<"error! hscale value must be between 0 and 1. current value = "<<h_scale<<endl; error=true;}
	if (rho_scale>1.0) {cout<<"error! rscale value must be between 0 and 1. current value = "<<rho_scale<<endl; error=true;}
	if (logp_scale>1.0) {cout<<"error! pscale value must be between 0 and 1. current value = "<<logp_scale<<endl; error=true;}

	if (rho_max==1 && rho_min==1 && a_mode==12) {cout<<"error! ridge regression does not support a rho parameter. current values = "<<rho_max<<" and "<<rho_min<<endl; error=true;}

	if (window_cm<0) {cout<<"error! windowcm values must be non-negative. current values = "<<window_cm<<endl; error=true;}

	if (window_cm==0 && window_bp==0 && window_ns==0) {
	  window_bp=1000000;
	}

	//check p_column, and (no need to) sort p_column into ascending order
	if (p_column.size()==0) {
		p_column.push_back(1);
	} else {
		for (size_t i=0; i<p_column.size(); i++) {
			for (size_t j=0; j<i; j++) {
				if (p_column[i]==p_column[j]) {cout<<"error! identical phenotype columns: "<<p_column[i]<<endl; error=true;}
			}
		}
	}

	//sort (p_column.begin(), p_column.end() );
	n_ph=p_column.size();



	//only lmm option (and one prediction option) can deal with multiple phenotypes
	//and no gene expression files
	if (n_ph>1 && a_mode!=1 && a_mode!=2 && a_mode!=3 && a_mode!=4 && a_mode!=43) {
		cout<<"error! the current analysis mode "<<a_mode<<" can not deal with multiple phenotypes."<<endl; error=true;
	}
	if (n_ph>1 && !file_gene.empty() ) {
		cout<<"error! multiple phenotype analysis option not allowed with gene expression files. "<<endl; error=true;
	}

	if (p_nr>1) {
		cout<<"error! pnr value must be between 0 and 1. current value = "<<p_nr<<endl; error=true;
	}

	//check est_column
	if (est_column.size()==0) {
		if (file_ebv.empty()) {
			est_column.push_back(2);
			est_column.push_back(5);
			est_column.push_back(6);
			est_column.push_back(7);
		} else {
			est_column.push_back(2);
			est_column.push_back(0);
			est_column.push_back(6);
			est_column.push_back(7);
		}
	}

	if (est_column.size()!=4) {cout<<"error! -en not followed by four numbers. current number = "<<est_column.size()<<endl; error=true;}
	if (est_column[0]==0) {cout<<"error! -en rs column can not be zero. current number = "<<est_column.size()<<endl; error=true;}

	//check if files are compatible with each other, and if files exist
	if (!file_bfile.empty()) {
		str=file_bfile+".bim";
		if (stat(str.c_str(),&fileInfo)==-1) {cout<<"error! fail to open .bim file: "<<str<<endl; error=true;}
		str=file_bfile+".bed";
		if (stat(str.c_str(),&fileInfo)==-1) {cout<<"error! fail to open .bed file: "<<str<<endl; error=true;}
		str=file_bfile+".fam";
		if (stat(str.c_str(),&fileInfo)==-1) {cout<<"error! fail to open .fam file: "<<str<<endl; error=true;}
	}

	if ((!file_geno.empty() || !file_gene.empty()) ) {
		str=file_pheno;
		if (stat(str.c_str(),&fileInfo)==-1) {cout<<"error! fail to open phenotype file: "<<str<<endl; error=true;}
	}

	str=file_geno;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open mean genotype file: "<<str<<endl; error=true;}

	str=file_gene;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open gene expression file: "<<str<<endl; error=true;}

	str=file_cat;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open category file: "<<str<<endl; error=true;}

	str=file_var;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open category file: "<<str<<endl; error=true;}

	str=file_beta;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open beta file: "<<str<<endl; error=true;}

	str=file_cor;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open correlation file: "<<str<<endl; error=true;}

	str=file_q;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open q file: "<<str<<endl; error=true;}

	str=file_s;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open s file: "<<str<<endl; error=true;}

	str=file_v;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open v file: "<<str<<endl; error=true;}

	str=file_mq;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open mq file: "<<str<<endl; error=true;}

	str=file_ms;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open ms file: "<<str<<endl; error=true;}

	str=file_mv;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open mv file: "<<str<<endl; error=true;}

	size_t flag=0;
	if (!file_bfile.empty()) {flag++;}
	if (!file_geno.empty()) {flag++;}
	if (!file_gene.empty()) {flag++;}

	if (flag!=1 && a_mode!=27 && a_mode!=28 && a_mode!=43 && a_mode!=5 && a_mode!=61 && a_mode!=62) {
		cout<<"error! either plink binary files, or bimbam mean genotype files, or gene expression files are required."<<endl; error=true;
	}

	if (file_pheno.empty() && (a_mode==43 || a_mode==5) ) {
		cout<<"error! phenotype file is required."<<endl; error=true;
	}

	if (a_mode==61 || a_mode==62) {
	  if (!file_pheno.empty()) {
	    if (file_kin.empty() && (file_ku.empty()||file_kd.empty()) && file_mk.empty() ) {
	      cout<<"error! missing relatedness file. "<<endl;  error=true;
	    }
	  } else if (!file_cor.empty()) {
	    if (file_beta.empty() ) {
	      cout<<"error! missing cor file."<<endl; error=true;
	    }
	  } else {
	    if ( (file_mq.empty() || file_ms.empty() || file_mv.empty() ) && (file_q.empty() || file_s.empty() || file_v.empty() )  ) {
	      cout<<"error! either phenotype/kinship files or ms/mq/mv s/q/v files are required."<<endl; error=true;
	    }
	  }
	}



	if (!file_epm.empty() && file_bfile.empty() && file_geno.empty() ) {cout<<"error! estimated parameter file also requires genotype file."<<endl; error=true;}
	if (!file_ebv.empty() && file_kin.empty()) {cout<<"error! estimated breeding value file also requires relatedness file."<<endl; error=true;}

	if (!file_log.empty() && pheno_mean!=0) {cout<<"error! either log file or mu value can be provide."<<endl; error=true;}

	str=file_snps;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open snps file: "<<str<<endl; error=true;}

	str=file_log;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open log file: "<<str<<endl; error=true;}

	str=file_anno;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open annotation file: "<<str<<endl; error=true;}

	str=file_kin;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open relatedness matrix file: "<<str<<endl; error=true;}

	str=file_mk;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open relatedness matrix file: "<<str<<endl; error=true;}

	str=file_cvt;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open covariates file: "<<str<<endl; error=true;}

	str=file_gxe;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open environmental covariate file: "<<str<<endl; error=true;}

	str=file_weight;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open the residual weight file: "<<str<<endl; error=true;}

	str=file_epm;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open estimated parameter file: "<<str<<endl; error=true;}

	str=file_ebv;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open estimated breeding value file: "<<str<<endl; error=true;}

	str=file_read;
	if (!str.empty() && stat(str.c_str(),&fileInfo)==-1 ) {cout<<"error! fail to open total read file: "<<str<<endl; error=true;}

	//check if files are compatible with analysis mode
	if (k_mode==2 && !file_geno.empty() ) {cout<<"error! use \"-km 1\" when using bimbam mean genotype file. "<<endl; error=true;}

	if ((a_mode==1 || a_mode==2 || a_mode==3 || a_mode==4 || a_mode==5 || a_mode==31) && (file_kin.empty() && (file_ku.empty()||file_kd.empty())) )  {cout<<"error! missing relatedness file. "<<endl;  error=true;}

	if ((a_mode==43) && file_kin.empty())  {cout<<"error! missing relatedness file. -predict option requires -k option to provide a relatedness file."<<endl;  error=true;}

	if ((a_mode==11 || a_mode==12 || a_mode==13) && !file_cvt.empty() ) {cout<<"error! -bslmm option does not support covariates files."<<endl; error=true;}

	if (a_mode==41 || a_mode==42) {
		if (!file_cvt.empty() ) {cout<<"error! -predict option does not support covariates files."<<endl; error=true;}
		if (file_epm.empty() ) {cout<<"error! -predict option requires estimated parameter files."<<endl; error=true;}
	}

	if (file_beta.empty() && (a_mode==27 || a_mode==28) ) {
		cout<<"error! beta effects file is required."<<endl; error=true;
	}

	return;
}





void PARAM::CheckData (void) {
	if ((file_cvt).empty() || (indicator_cvt).size()==0) {
		n_cvt=1;
	}
	if ( (indicator_cvt).size()!=0 && (indicator_cvt).size()!=(indicator_idv).size()) {
		error=true;
		cout<<"error! number of rows in the covariates file do not match the number of individuals. "<<endl;
		return;
	}
	if ( (indicator_gxe).size()!=0 && (indicator_gxe).size()!=(indicator_idv).size()) {
		error=true;
		cout<<"error! number of rows in the gxe file do not match the number of individuals. "<<endl;
		return;
	}
	if ( (indicator_weight).size()!=0 && (indicator_weight).size()!=(indicator_idv).size()) {
		error=true;
		cout<<"error! number of rows in the weight file do not match the number of individuals. "<<endl;
		return;
	}

	if ( (indicator_read).size()!=0 && (indicator_read).size()!=(indicator_idv).size()) {
		error=true;
		cout<<"error! number of rows in the total read file do not match the number of individuals. "<<endl;
		return;
	}

	//calculate ni_total and ni_test, and set indicator_idv to 0 whenever indicator_cvt=0
	//and calculate np_obs and np_miss
	ni_total=(indicator_idv).size();

	ni_test=0;
	for (vector<int>::size_type i=0; i<(indicator_idv).size(); ++i) {
		if (indicator_idv[i]==0) {continue;}
		ni_test++;
	}

	ni_cvt=0;
	for (size_t i=0; i<indicator_cvt.size(); i++) {
		if (indicator_cvt[i]==0) {continue;}
		ni_cvt++;
	}

	np_obs=0; np_miss=0;
	for (size_t i=0; i<indicator_pheno.size(); i++) {
		if (indicator_cvt.size()!=0) {
			if (indicator_cvt[i]==0) {continue;}
		}

		if (indicator_gxe.size()!=0) {
			if (indicator_gxe[i]==0) {continue;}
		}

		if (indicator_weight.size()!=0) {
			if (indicator_weight[i]==0) {continue;}
		}

		for (size_t j=0; j<indicator_pheno[i].size(); j++) {
			if (indicator_pheno[i][j]==0) {
				np_miss++;
			} else {
				np_obs++;
			}
		}
	}

	/*
	if ((indicator_cvt).size()!=0) {
		ni_test=0;
		for (vector<int>::size_type i=0; i<(indicator_idv).size(); ++i) {
			indicator_idv[i]*=indicator_cvt[i];
			ni_test+=indicator_idv[i];
		}
	}

	if ((indicator_read).size()!=0) {
		ni_test=0;
		for (vector<int>::size_type i=0; i<(indicator_idv).size(); ++i) {
			indicator_idv[i]*=indicator_read[i];
			ni_test+=indicator_idv[i];
		}
	}
	*/
	if (ni_test==0 && file_cor.empty() && file_mq.empty() && file_q.empty() && file_beta.empty() ) {
		error=true;
		cout<<"error! number of analyzed individuals equals 0. "<<endl;
		return;
	}

	if (a_mode==43) {
		if (ni_cvt==ni_test) {
			error=true;
			cout<<"error! no individual has missing phenotypes."<<endl;
			return;
		}
		if ((np_obs+np_miss)!=(ni_cvt*n_ph)) {
			error=true;
			//cout<<ni_cvt<<"\t"<<ni_test<<"\t"<<ni_total<<"\t"<<np_obs<<"\t"<<np_miss<<"\t"<<indicator_cvt.size()<<endl;
			cout<<"error! number of phenotypes do not match the summation of missing and observed phenotypes."<<endl;
			return;
		}
	}

	//output some information
	if (file_cor.empty() && file_mq.empty() && file_q.empty() ) {
	  cout<<"## number of total individuals = "<<ni_total<<endl;
	  if (a_mode==43) {
	    cout<<"## number of analyzed individuals = "<<ni_cvt<<endl;
	    cout<<"## number of individuals with full phenotypes = "<<ni_test<<endl;
	  } else {
	    cout<<"## number of analyzed individuals = "<<ni_test<<endl;
	  }
	  cout<<"## number of covariates = "<<n_cvt<<endl;
	  cout<<"## number of phenotypes = "<<n_ph<<endl;
	  if (a_mode==43) {
	    cout<<"## number of observed data = "<<np_obs<<endl;
	    cout<<"## number of missing data = "<<np_miss<<endl;
	  }
	  if (!file_gene.empty()) {
	    cout<<"## number of total genes = "<<ng_total<<endl;
	  } else if (file_epm.empty() && a_mode!=43 && a_mode!=5) {
	    cout<<"## number of total SNPs = "<<ns_total<<endl;
	    cout<<"## number of analyzed SNPs = "<<ns_test<<endl;
	  } else {}
	}

	//set d_pace to 1000 for gene expression
	if (!file_gene.empty() && d_pace==100000) {
		d_pace=1000;
	}

	//for case-control studies, count #cases and #controls
	int flag_cc=0;
	if (a_mode==13) {
		ni_case=0;
		ni_control=0;
		for (size_t i=0; i<indicator_idv.size(); i++) {
			if (indicator_idv[i]==0) {continue;}

			if (pheno[i][0]==0) {ni_control++;}
			else if (pheno[i][0]==1) {ni_case++;}
			else {flag_cc=1;}
		}
		cout<<"## number of cases = "<<ni_case<<endl;
		cout<<"## number of controls = "<<ni_control<<endl;
	}

	if (flag_cc==1) {cout<<"Unexpected non-binary phenotypes for case/control analysis. Use default (BSLMM) analysis instead."<<endl; a_mode=11;}

	//set parameters for BSLMM
	//and check for predict
	if (a_mode==11 || a_mode==12 || a_mode==13) {
		if (a_mode==11) {n_mh=1;}
		if (logp_min==0) {logp_min=-1.0*log((double)ns_test);}

		if (h_scale==-1) {h_scale=min(1.0, 10.0/sqrt((double)ni_test) );}
		if (rho_scale==-1) {rho_scale=min(1.0, 10.0/sqrt((double)ni_test) );}
		if (logp_scale==-1) {logp_scale=min(1.0, 5.0/sqrt((double)ni_test) );}

		if (h_min==-1) {h_min=0.0;}
		if (h_max==-1) {h_max=1.0;}

		if (s_max>ns_test) {s_max=ns_test; cout<<"s_max is re-set to the number of analyzed SNPs."<<endl;}
		if (s_max<s_min) {cout<<"error! maximum s value must be larger than the minimal value. current values = "<<s_max<<" and "<<s_min<<endl; error=true;}
	} else if (a_mode==41 || a_mode==42) {
		if (indicator_bv.size()!=0) {
			if (indicator_idv.size()!=indicator_bv.size()) {
				cout<<"error! number of rows in the phenotype file does not match that in the estimated breeding value file: "<<indicator_idv.size()<<"\t"<<indicator_bv.size()<<endl;
				error=true;
			} else {
				size_t flag_bv=0;
				for (size_t i=0; i<(indicator_bv).size(); ++i) {
					if (indicator_idv[i]!=indicator_bv[i]) {flag_bv++;}
				}
				if (flag_bv!=0) {
					cout<<"error! individuals with missing value in the phenotype file does not match that in the estimated breeding value file: "<<flag_bv<<endl;
					error=true;
				}
			}
		}
	}

	//file_mk needs to contain more than one line
	if (n_vc==1 && !file_mk.empty()) {cout<<"error! -mk file should contain more than one line."<<endl; error=true;}

	return;
}


void PARAM::PrintSummary ()
{
	if (n_ph==1) {
		cout<<"pve estimate ="<<pve_null<<endl;
		cout<<"se(pve) ="<<pve_se_null<<endl;
	} else {

	}
	return;
}



void PARAM::ReadGenotypes (gsl_matrix *UtX, gsl_matrix *K, const bool calc_K) {
	string file_str;

	if (!file_bfile.empty()) {
		file_str=file_bfile+".bed";
		if (ReadFile_bed (file_str, indicator_idv, indicator_snp, UtX, K, calc_K)==false) {error=true;}
	}
	else {
		if (ReadFile_geno (file_geno, indicator_idv, indicator_snp, UtX, K, calc_K)==false) {error=true;}
	}

	return;
}


void PARAM::ReadGenotypes (vector<vector<unsigned char> > &Xt, gsl_matrix *K, const bool calc_K) {
	string file_str;

	if (!file_bfile.empty()) {
		file_str=file_bfile+".bed";
		if (ReadFile_bed (file_str, indicator_idv, indicator_snp, Xt, K, calc_K, ni_test, ns_test)==false) {error=true;}
	} else {
		if (ReadFile_geno (file_geno, indicator_idv, indicator_snp, Xt, K, calc_K, ni_test, ns_test)==false) {error=true;}
	}

	return;
}




void PARAM::CalcKin (gsl_matrix *matrix_kin)  {
	string file_str;

	gsl_matrix_set_zero (matrix_kin);

	if (!file_bfile.empty() ) {
		file_str=file_bfile+".bed";
		if (PlinkKin (file_str, indicator_snp, a_mode-20, d_pace, matrix_kin)==false) {error=true;}
	}
	else {
		file_str=file_geno;
		if (BimbamKin (file_str, indicator_snp, a_mode-20, d_pace, matrix_kin)==false) {error=true;}
	}

	return;
}



//from an existing n by nd G matrix, compute the d by d S matrix
void compKtoS (const gsl_matrix *G, gsl_matrix *S) {
  size_t n_vc=S->size1, ni_test=G->size1;
  double di, dj, tr_KiKj, sum_Ki, sum_Kj, s_Ki, s_Kj, s_KiKj, si, sj, d;

  for (size_t i=0; i<n_vc; i++) {
    for (size_t j=i; j<n_vc; j++) {
      tr_KiKj=0; sum_Ki=0; sum_Kj=0; s_KiKj=0; si=0; sj=0;
      for (size_t l=0; l<ni_test; l++) {
	s_Ki=0; s_Kj=0;
	for (size_t k=0; k<ni_test; k++) {
	  di=gsl_matrix_get(G, l, k+ni_test*i);
	  dj=gsl_matrix_get(G, l, k+ni_test*j);
	  s_Ki+=di; s_Kj+=dj;

	  tr_KiKj+=di*dj; sum_Ki+=di; sum_Kj+=dj;
	  if (l==k) {si+=di; sj+=dj;}
	}
	s_KiKj+=s_Ki*s_Kj;
      }

      sum_Ki/=(double)ni_test;
      sum_Kj/=(double)ni_test;
      s_KiKj/=(double)ni_test;
      si-=sum_Ki;
      sj-=sum_Kj;
      d=tr_KiKj-2*s_KiKj+sum_Ki*sum_Kj;
      d=d/(si*sj)-1/(double)(ni_test-1);

      gsl_matrix_set (S, i, j, d);
      if (i!=j) {gsl_matrix_set (S, j, i, d);}
    }
  }
  //cout<<tr_KiKj<<" "<<s_KiKj<<" "<<sum_Ki<<" "<<sum_Kj<<" "<<si<<" "<<sj<<" "<<d*1000000<<endl;
  return;
}



//copied from lmm.cpp; is used in the following function compKtoQ
//map a number 1-(n_cvt+2) to an index between 0 and [(n_c+2)^2+(n_c+2)]/2-1
size_t GetabIndex (const size_t a, const size_t b, const size_t n_cvt) {
	if (a>n_cvt+2 || b>n_cvt+2 || a<=0 || b<=0) {cout<<"error in GetabIndex."<<endl; return 0;}
	size_t index;
	size_t l, h;
	if (b>a) {l=a; h=b;} else {l=b; h=a;}

	size_t n=n_cvt+2;
	index=(2*n-l+2)*(l-1)/2+h-l;

	return index;
}

//from an existing n by nd (centered) G matrix, compute the d+1 by d*(d+1) Q matrix
//where inside i'th d+1 by d+1 matrix, each element is tr(KiKjKiKl)-r*tr(KjKiKl)-r*tr(KlKiKj)+r^2*tr(KjKl), where r=n/(n-1)
void compKtoQ (const gsl_matrix *G, gsl_matrix *Q) {
  size_t n_vc=G->size2/G->size1, ni_test=G->size1;

  gsl_matrix *KiKj=gsl_matrix_alloc(ni_test, n_vc*(n_vc+1)/2*ni_test);
  gsl_vector *trKiKjKi=gsl_vector_alloc ( n_vc*n_vc );
  gsl_vector *trKiKj=gsl_vector_alloc( n_vc*(n_vc+1)/2 );
  gsl_vector *trKi=gsl_vector_alloc(n_vc);

  double d, tr, r=(double)ni_test/(double)(ni_test-1);
  size_t t, t_ij, t_il, t_jl, t_ii;

  //compute KiKj for all pairs of i and j (including the identity matrix)
  t=0;
  for (size_t i=0; i<n_vc; i++) {
    gsl_matrix_const_view Ki=gsl_matrix_const_submatrix(G, 0, i*ni_test, ni_test, ni_test);
    for (size_t j=i; j<n_vc; j++) {
      gsl_matrix_const_view Kj=gsl_matrix_const_submatrix(G, 0, j*ni_test, ni_test, ni_test);
      gsl_matrix_view KiKj_sub=gsl_matrix_submatrix (KiKj, 0, t*ni_test, ni_test, ni_test);
      eigenlib_dgemm ("N", "N", 1.0, &Ki.matrix, &Kj.matrix, 0.0, &KiKj_sub.matrix);
      t++;
    }
  }
  /*
  for (size_t i=0; i<5; i++) {
    for (size_t j=0; j<5; j++) {
      cout<<gsl_matrix_get (G, i, j)<<" ";
    }
    cout<<endl;
  }
  */

  //compute trKi, trKiKj
  t=0;
  for (size_t i=0; i<n_vc; i++) {
    for (size_t j=i; j<n_vc; j++) {
      tr=0;
      for (size_t k=0; k<ni_test; k++) {
	tr+=gsl_matrix_get (KiKj, k, t*ni_test+k);
      }
      gsl_vector_set (trKiKj, t, tr);

      t++;
    }

    tr=0;
    for (size_t k=0; k<ni_test; k++) {
      tr+=gsl_matrix_get (G, k, i*ni_test+k);
    }
    gsl_vector_set (trKi, i, tr);
  }

  //compute trKiKjKi (it is not symmetric w.r.t. i and j)
  for (size_t i=0; i<n_vc; i++) {
    for (size_t j=0; j<n_vc; j++) {
      tr=0;
      t=GetabIndex (i+1, j+1, n_vc-2);
      for (size_t k=0; k<ni_test; k++) {
	gsl_vector_const_view KiKj_row=gsl_matrix_const_subrow (KiKj, k, t*ni_test, ni_test);
	gsl_vector_const_view KiKj_col=gsl_matrix_const_column (KiKj, t*ni_test+k);

	gsl_vector_const_view Ki_col=gsl_matrix_const_column (G, i*ni_test+k);

	if (i<=j) {
	  gsl_blas_ddot (&KiKj_row.vector, &Ki_col.vector, &d);
	  tr+=d;
	} else {
	  gsl_blas_ddot (&KiKj_col.vector, &Ki_col.vector, &d);
	  tr+=d;
	}
      }
      gsl_vector_set (trKiKjKi, i*n_vc+j, tr);
    }
  }

  //compute Q
  for (size_t i=0; i<n_vc; i++) {
    for (size_t j=0; j<n_vc+1; j++) {
      for (size_t l=j; l<n_vc+1; l++) {
	if (j!=n_vc && l!=n_vc) {
	  t_ij=GetabIndex (i+1, j+1, n_vc-2);
	  t_il=GetabIndex (i+1, l+1, n_vc-2);
	  t_jl=GetabIndex (j+1, l+1, n_vc-2);

	  //cout<<ni_test<<" "<<r<<t_ij<<" "<<t_il<<" "<<t_jl<<" "<<endl;
	  tr=0;
	  for (size_t k=0; k<ni_test; k++) {
	    gsl_vector_const_view KiKj_row=gsl_matrix_const_subrow (KiKj, k, t_ij*ni_test, ni_test);
	    gsl_vector_const_view KiKj_col=gsl_matrix_const_column (KiKj, t_ij*ni_test+k);
	    gsl_vector_const_view KiKl_row=gsl_matrix_const_subrow (KiKj, k, t_il*ni_test, ni_test);
	    gsl_vector_const_view KiKl_col=gsl_matrix_const_column (KiKj, t_il*ni_test+k);

	    gsl_vector_const_view Kj_row=gsl_matrix_const_subrow (G, k, j*ni_test, ni_test);
	    gsl_vector_const_view Kl_row=gsl_matrix_const_subrow (G, k, l*ni_test, ni_test);

	    if (i<=j && i<=l) {
	      gsl_blas_ddot (&KiKj_row.vector, &KiKl_col.vector, &d);
	      tr+=d;
	      gsl_blas_ddot (&Kj_row.vector, &KiKl_col.vector, &d);
	      tr-=r*d;
	      gsl_blas_ddot (&Kl_row.vector, &KiKj_col.vector, &d);
	      tr-=r*d;
	    } else if (i<=j && i>l) {
	      gsl_blas_ddot (&KiKj_row.vector, &KiKl_row.vector, &d);
	      tr+=d;
	      gsl_blas_ddot (&Kj_row.vector, &KiKl_row.vector, &d);
	      tr-=r*d;
	      gsl_blas_ddot (&Kl_row.vector, &KiKj_col.vector, &d);
	      tr-=r*d;
	    } else if (i>j && i<=l) {
	      gsl_blas_ddot (&KiKj_col.vector, &KiKl_col.vector, &d);
	      tr+=d;
	      gsl_blas_ddot (&Kj_row.vector, &KiKl_col.vector, &d);
	      tr-=r*d;
	      gsl_blas_ddot (&Kl_row.vector, &KiKj_row.vector, &d);
	      tr-=r*d;
	    } else {
	      gsl_blas_ddot (&KiKj_col.vector, &KiKl_row.vector, &d);
	      tr+=d;
	      gsl_blas_ddot (&Kj_row.vector, &KiKl_row.vector, &d);
	      tr-=r*d;
	      gsl_blas_ddot (&Kl_row.vector, &KiKj_row.vector, &d);
	      tr-=r*d;
	    }
	  }

	  tr+=r*r*gsl_vector_get (trKiKj, t_jl);
	} else if (j!=n_vc && l==n_vc) {
	  t_ij=GetabIndex (i+1, j+1, n_vc-2);
	  tr=gsl_vector_get (trKiKjKi, i*n_vc+j)-2*r*gsl_vector_get (trKiKj, t_ij)+r*r*gsl_vector_get (trKi, j);
	} else if (j==n_vc && l==n_vc) {
	  t_ii=GetabIndex (i+1, i+1, n_vc-2);
	  tr=gsl_vector_get (trKiKj, t_ii)-2*r*gsl_vector_get (trKi, i)+r*r*(double)(ni_test-1);
	}

	gsl_matrix_set (Q, j, i*(n_vc+1)+l, tr);
	if (l!=j) {gsl_matrix_set (Q, l, i*(n_vc+1)+j, tr);}
      }
    }
  }

  gsl_matrix_scale (Q, 1.0/pow((double)ni_test, 2) );

  gsl_matrix_free(KiKj);
  gsl_vector_free(trKiKjKi);
  gsl_vector_free(trKiKj);
  gsl_vector_free(trKi);

  return;
}



//perform Jacknife sampling for variance of S
void JacknifeGtoS (const gsl_matrix *G, gsl_matrix *S, gsl_matrix *Svar) {
  size_t n_vc=Svar->size1, ni_test=G->size1;
  vector<vector<vector<double> > > tr_KiKj, s_KiKj;
  vector<vector<double> > sum_Ki, s_Ki, si;
  vector<double> vec_tmp;
  double di, dj, d, m, v;

  //initialize and set all elements to zero
  for (size_t i=0; i<ni_test; i++) {
    vec_tmp.push_back(0);
  }

  for (size_t i=0; i<n_vc; i++) {
    sum_Ki.push_back(vec_tmp);
    s_Ki.push_back(vec_tmp);
    si.push_back(vec_tmp);
  }

  for (size_t i=0; i<n_vc; i++) {
    tr_KiKj.push_back(sum_Ki);
    s_KiKj.push_back(sum_Ki);
  }

  //run jacknife
  for (size_t i=0; i<n_vc; i++) {
    for (size_t l=0; l<ni_test; l++) {
      for (size_t k=0; k<ni_test; k++) {
	di=gsl_matrix_get(G, l, k+ni_test*i);

	for (size_t t=0; t<ni_test; t++) {
	  if (t==l || t==k) {continue;}
	  sum_Ki[i][t]+=di;
	  if (l==k) {si[i][t]+=di;}
	}
	s_Ki[i][l]+=di;
      }
    }

    for (size_t t=0; t<ni_test; t++) {
      sum_Ki[i][t]/=(double)(ni_test-1);
    }
  }

  for (size_t i=0; i<n_vc; i++) {
    for (size_t j=i; j<n_vc; j++) {
      for (size_t l=0; l<ni_test; l++) {
	for (size_t k=0; k<ni_test; k++) {
	  di=gsl_matrix_get(G, l, k+ni_test*i);
	  dj=gsl_matrix_get(G, l, k+ni_test*j);
	  d=di*dj;

	  for (size_t t=0; t<ni_test; t++) {
	    if (t==l || t==k) {continue;}
	    tr_KiKj[i][j][t]+=d;
          }
	}

	for (size_t t=0; t<ni_test; t++) {
	  if (t==l) {continue;}
	  di=gsl_matrix_get(G, l, t+ni_test*i);
	  dj=gsl_matrix_get(G, l, t+ni_test*j);

	  s_KiKj[i][j][t]+=(s_Ki[i][l]-di)*(s_Ki[j][l]-dj);
	}
      }

      for (size_t t=0; t<ni_test; t++) {
	s_KiKj[i][j][t]/=(double)(ni_test-1);
      }

      m=0; v=0;
      for (size_t t=0; t<ni_test; t++) {
	d=tr_KiKj[i][j][t]-2*s_KiKj[i][j][t]+sum_Ki[i][t]*sum_Ki[j][t];
	d/=(si[i][t]-sum_Ki[i][t])*(si[j][t]-sum_Ki[j][t]);
	d-=1/(double)(ni_test-2);

	m+=d; v+=d*d;
      }
      m/=(double)ni_test;
      v/=(double)ni_test;
      v-=m*m;
      v*=(double)(ni_test-1);

      gsl_matrix_set (Svar, i, j, v);
      d=gsl_matrix_get (S, i, j);
      d=(double)ni_test*d-(double)(ni_test-1)*m;
      gsl_matrix_set (S, i, j, d);
      if (i!=j) {gsl_matrix_set (Svar, j, i, v); gsl_matrix_set (S, j, i, d);}
    }
  }

  return;
}



//compute the d by d S matrix with its d by d variance matrix of Svar, and the d+1 by d(d+1) matrix of Q for V(q)
void PARAM::CalcS (gsl_matrix *S, gsl_matrix *Svar, gsl_matrix *Q)  {
  string file_str;

  gsl_matrix_set_zero (S);
  gsl_matrix_set_zero (Svar);
  gsl_matrix_set_zero (Q);

  //compute the kinship matrix G for multiple categories; these matrices are not centered, for convienence of Jacknife sampling
  gsl_matrix *G=gsl_matrix_alloc (ni_test, n_vc*ni_test);
  gsl_matrix_set_zero (G);

  if (!file_bfile.empty() ) {
    file_str=file_bfile+".bed";
    if (PlinkKin (file_str, indicator_idv, indicator_snp, a_mode-24, d_pace, mapRS2cat, mapRS2var, snpInfo, G)==false) {error=true;}
  } else {
    file_str=file_geno;
    if (BimbamKin (file_str, indicator_idv, indicator_snp, a_mode-24, d_pace, mapRS2cat, mapRS2var, snpInfo, G)==false) {error=true;}
  }

  //center and scale every kinship matrix inside G
  double d;
  for (size_t i=0; i<n_vc; i++) {
    gsl_matrix_view K=gsl_matrix_submatrix(G, 0, i*ni_test, ni_test, ni_test);
    CenterMatrix(&K.matrix);
    d=ScaleMatrix(&K.matrix);
  }

  //based on G, compute S
  compKtoS (G, S);

  //based on G, compute a matrix Q that can be used to calculate the variance of q
  compKtoQ (G, Q);

  /*
  //set up random environment
  gsl_rng_env_setup();
  gsl_rng *gsl_r;
  const gsl_rng_type * gslType;
  gslType = gsl_rng_default;
  if (randseed<0) {
    time_t rawtime;
    time (&rawtime);
    tm * ptm = gmtime (&rawtime);

    randseed = (unsigned) (ptm->tm_hour%24*3600+ptm->tm_min*60+ptm->tm_sec);
  }
  gsl_r = gsl_rng_alloc(gslType);
  gsl_rng_set(gsl_r, randseed);

  //bootstrap: in each iteration, sample individuals and compute S_pmt
  size_t n_pmt=100;
  vector<size_t> idv_order, idv_remove;
  for (size_t i=0; i<ni_test; i++) {
    idv_order.push_back(i);
  }
  for (size_t i=0; i<n_pmt; i++) {
    idv_remove.push_back(0);
  }
  gsl_ran_choose (gsl_r, static_cast<void*>(&idv_remove[0]), n_pmt, static_cast<void*>(&idv_order[0]), ni_test, sizeof(size_t));

  gsl_matrix *S_pmt=gsl_matrix_alloc(n_vc, n_vc*n_pmt);
  for (size_t i=0; i<n_pmt; i++) {
    gsl_matrix_view S_sub=gsl_matrix_submatrix (S_pmt, 0, n_vc*i, n_vc, n_vc);
    compKtoS (G, idv_remove[i], &S_sub.matrix);
  }

  //based on S_pmt, compute Svar
  double m, v, d;
  for (size_t i=0; i<n_vc; i++) {
    for (size_t j=i; j<n_vc; j++) {
      m=0; v=0;
      for (size_t t=0; t<n_pmt; t++) {
	d=gsl_matrix_get(S_pmt, i, j);
	m+=d; v+=d*d;
      }
      m/=(double)n_pmt; v/=(double)n_pmt;
      v=v-m*m;
      gsl_matrix_set(Svar, i, j, v);
      if (i!=j) {gsl_matrix_set(Svar, j, i, v);}
    }
  }
  */

  //compute Svar and update S with Jacknife
  JacknifeGtoS (G, S, Svar);

  gsl_matrix_free(G);
  return;
}



void PARAM::WriteVector (const gsl_vector *q, const gsl_vector *s, const size_t n_total, const string suffix)
{
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".";
	file_str+=suffix;
	file_str+=".txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

	outfile.precision(10);

	for (size_t i=0; i<q->size; ++i) {
	  outfile<<gsl_vector_get (q, i)<<endl;
	}

	for (size_t i=0; i<s->size; ++i) {
	  outfile<<gsl_vector_get (s, i)<<endl;
	}

	outfile<<n_total<<endl;

	outfile.close();
	outfile.clear();
	return;
}



void PARAM::WriteVar (const string suffix)
{
  string file_str, rs;
	file_str=path_out+"/"+file_out;
	file_str+=".";
	file_str+=suffix;
	file_str+=".txt.gz";

	ogzstream outfile (file_str.c_str(), ogzstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

	outfile.precision(10);

	for (size_t i=0; i<indicator_snp.size(); i++) {
	  if (indicator_snp[i]==0) {continue;}
	  rs=snpInfo[i].rs_number;
	  if (mapRS2var.count(rs)!=0) {
	    outfile<<rs<<"\t"<<mapRS2var.at(rs)<<endl;
	  }
	}

	outfile.close();
	outfile.clear();
	return;
}




void PARAM::WriteMatrix (const gsl_matrix *matrix_U, const string suffix)
{
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".";
	file_str+=suffix;
	file_str+=".txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

	outfile.precision(10);

	for (size_t i=0; i<matrix_U->size1; ++i) {
		for (size_t j=0; j<matrix_U->size2; ++j) {
			outfile<<gsl_matrix_get (matrix_U, i, j)<<"\t";
		}
		outfile<<endl;
	}

	outfile.close();
	outfile.clear();
	return;
}


void PARAM::WriteVector (const gsl_vector *vector_D, const string suffix)
{
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".";
	file_str+=suffix;
	file_str+=".txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

	outfile.precision(10);

	for (size_t i=0; i<vector_D->size; ++i) {
		outfile<<gsl_vector_get (vector_D, i)<<endl;
	}

	outfile.close();
	outfile.clear();
	return;
}


void PARAM::CheckCvt ()
{
	if (indicator_cvt.size()==0) {return;}

	size_t ci_test=0;

	gsl_matrix *W=gsl_matrix_alloc (ni_test, n_cvt);

	for (vector<int>::size_type i=0; i<indicator_idv.size(); ++i) {
		if (indicator_idv[i]==0 || indicator_cvt[i]==0) {continue;}
		for (size_t j=0; j<n_cvt; ++j) {
			gsl_matrix_set (W, ci_test, j, (cvt)[i][j]);
		}
		ci_test++;
	}

	size_t flag_ipt=0;
	double v_min, v_max;
	set<size_t> set_remove;

	//check if any columns is an intercept
	for (size_t i=0; i<W->size2; i++) {
		gsl_vector_view w_col=gsl_matrix_column (W, i);
		gsl_vector_minmax (&w_col.vector, &v_min, &v_max);
		if (v_min==v_max) {flag_ipt=1; set_remove.insert (i);}
	}

	//add an intecept term if needed
	if (n_cvt==set_remove.size()) {
		indicator_cvt.clear();
		n_cvt=1;
	} else if (flag_ipt==0) {
		cout<<"no intecept term is found in the cvt file. a column of 1s is added."<<endl;
		for (vector<int>::size_type i=0; i<indicator_idv.size(); ++i) {
			if (indicator_idv[i]==0 || indicator_cvt[i]==0) {continue;}
			cvt[i].push_back(1.0);
		}

		n_cvt++;
	} else {}

	gsl_matrix_free(W);

	return;
}


//post-process phentoypes, covariates
void PARAM::ProcessCvtPhen ()
{
	//convert indicator_pheno to indicator_idv
	int k=1;
	indicator_idv.clear();
	for (size_t i=0; i<indicator_pheno.size(); i++) {
		k=1;
		for (size_t j=0; j<indicator_pheno[i].size(); j++) {
			if (indicator_pheno[i][j]==0) {k=0;}
		}
		indicator_idv.push_back(k);
	}

	//remove individuals with missing covariates
	if ((indicator_cvt).size()!=0) {
		for (vector<int>::size_type i=0; i<(indicator_idv).size(); ++i) {
			indicator_idv[i]*=indicator_cvt[i];
		}
	}

	//remove individuals with missing gxe variables
	if ((indicator_gxe).size()!=0) {
		for (vector<int>::size_type i=0; i<(indicator_idv).size(); ++i) {
			indicator_idv[i]*=indicator_gxe[i];
		}
	}

	//remove individuals with missing residual weights
	if ((indicator_weight).size()!=0) {
		for (vector<int>::size_type i=0; i<(indicator_idv).size(); ++i) {
			indicator_idv[i]*=indicator_weight[i];
		}
	}

	//obtain ni_test
	ni_test=0;
	for (vector<int>::size_type i=0; i<(indicator_idv).size(); ++i) {
	    if (indicator_idv[i]==0) {continue;}
		ni_test++;
	}



	//if subsample number is set, perform a random sub-sampling to determine the subsampled ids
	if (ni_subsample!=0) {
	  if (ni_test<ni_subsample) {
	    cout<<"error! number of subsamples is less than number of analyzed individuals. "<<endl;
	  } else {
	    //set up random environment
	    gsl_rng_env_setup();
	    gsl_rng *gsl_r;
	    const gsl_rng_type * gslType;
	    gslType = gsl_rng_default;
	    if (randseed<0) {
	      time_t rawtime;
	      time (&rawtime);
	      tm * ptm = gmtime (&rawtime);

	      randseed = (unsigned) (ptm->tm_hour%24*3600+ptm->tm_min*60+ptm->tm_sec);
	    }
	    gsl_r = gsl_rng_alloc(gslType);
	    gsl_rng_set(gsl_r, randseed);

	    //from ni_test, sub-sample ni_subsample
	    vector<size_t> a, b;
	    for (size_t i=0; i<ni_subsample; i++) {
              a.push_back(0);
	    }
	    for (size_t i=0; i<ni_test; i++) {
	      b.push_back(i);
	    }

	    gsl_ran_choose (gsl_r, static_cast<void*>(&a[0]), ni_subsample, static_cast<void*>(&b[0]), ni_test, sizeof (size_t) );

	    //re-set indicator_idv and ni_test
	    int j=0;
	    for (vector<int>::size_type i=0; i<(indicator_idv).size(); ++i) {
	      if (indicator_idv[i]==0) {continue;}
	      if(find(a.begin(), a.end(), j) == a.end()) {
		indicator_idv[i]=0;
	      }
	      j++;
	    }
	    ni_test=ni_subsample;
	  }
	}

	//check ni_test
	if (ni_test==0) {
		error=true;
		cout<<"error! number of analyzed individuals equals 0. "<<endl;
		return;
	}

	//check covariates to see if they are correlated with each other, and to see if the intercept term is included
	//after getting ni_test
	//add or remove covariates
	if (indicator_cvt.size()!=0) {
		CheckCvt();
	} else {
		vector<double> cvt_row;
		cvt_row.push_back(1);

		for (vector<int>::size_type i=0; i<(indicator_idv).size(); ++i) {
			indicator_cvt.push_back(1);

			cvt.push_back(cvt_row);
		}
	}

	return;
}




void PARAM::CopyCvt (gsl_matrix *W)
{
	size_t ci_test=0;

	for (vector<int>::size_type i=0; i<indicator_idv.size(); ++i) {
		if (indicator_idv[i]==0 || indicator_cvt[i]==0) {continue;}
		for (size_t j=0; j<n_cvt; ++j) {
			gsl_matrix_set (W, ci_test, j, (cvt)[i][j]);
		}
		ci_test++;
	}

	return;
}


void PARAM::CopyGxe (gsl_vector *env)
{
	size_t ci_test=0;

	for (vector<int>::size_type i=0; i<indicator_idv.size(); ++i) {
		if (indicator_idv[i]==0 || indicator_gxe[i]==0) {continue;}
		gsl_vector_set (env, ci_test, gxe[i]);
		ci_test++;
	}

	return;
}

void PARAM::CopyWeight (gsl_vector *w)
{
	size_t ci_test=0;

	for (vector<int>::size_type i=0; i<indicator_idv.size(); ++i) {
		if (indicator_idv[i]==0 || indicator_weight[i]==0) {continue;}
		gsl_vector_set (w, ci_test, weight[i]);
		ci_test++;
	}

	return;
}


//if flag=0, then use indicator_idv to load W and Y
//else, use indicator_cvt to load them
void PARAM::CopyCvtPhen (gsl_matrix *W, gsl_vector *y, size_t flag)
{
	size_t ci_test=0;

	for (vector<int>::size_type i=0; i<indicator_idv.size(); ++i) {
		if (flag==0) {
			if (indicator_idv[i]==0) {continue;}
		} else {
			if (indicator_cvt[i]==0) {continue;}
		}

		gsl_vector_set (y, ci_test, (pheno)[i][0]);

		for (size_t j=0; j<n_cvt; ++j) {
			gsl_matrix_set (W, ci_test, j, (cvt)[i][j]);
		}
		ci_test++;
	}

	return;
}

//if flag=0, then use indicator_idv to load W and Y
//else, use indicator_cvt to load them
void PARAM::CopyCvtPhen (gsl_matrix *W, gsl_matrix *Y, size_t flag)
{
	size_t ci_test=0;

	for (vector<int>::size_type i=0; i<indicator_idv.size(); ++i) {
		if (flag==0) {
			if (indicator_idv[i]==0) {continue;}
		} else {
			if (indicator_cvt[i]==0) {continue;}
		}

        for (size_t j=0; j<n_ph; ++j) {
			gsl_matrix_set (Y, ci_test, j, (pheno)[i][j]);
		}
		for (size_t j=0; j<n_cvt; ++j) {
			gsl_matrix_set (W, ci_test, j, (cvt)[i][j]);
		}

		ci_test++;
	}

	return;
}





void PARAM::CopyRead (gsl_vector *log_N)
{
	size_t ci_test=0;

	for (vector<int>::size_type i=0; i<indicator_idv.size(); ++i) {
		if (indicator_idv[i]==0) {continue;}
		gsl_vector_set (log_N, ci_test, log(vec_read[i]) );
		ci_test++;
	}

	return;
}



