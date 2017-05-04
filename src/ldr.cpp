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
#include <sstream>

#include <iomanip>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h> 
#include <ctime>
#include <cstring>
#include <algorithm>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_roots.h"
#include "Eigen/Dense"



#include "lapack.h"

#ifdef FORCE_FLOAT
#include "param_float.h"
#include "ldr_float.h"
#include "lm_float.h"
#include "mathfunc_float.h"  //for function CenterVector
#else
#include "param.h"
#include "ldr.h"
#include "lm.h"
#include "mathfunc.h"
#endif

using namespace std;
using namespace Eigen;



void LDR::CopyFromParam (PARAM &cPar) 
{
	a_mode=cPar.a_mode;
	d_pace=cPar.d_pace;
	
	file_bfile=cPar.file_bfile;
	file_geno=cPar.file_geno;
	file_out=cPar.file_out;
	path_out=cPar.path_out;

	ni_total=cPar.ni_total;
	ns_total=cPar.ns_total;
	ni_test=cPar.ni_test;
	ns_test=cPar.ns_test;
	n_cvt=cPar.n_cvt;
	
	indicator_idv=cPar.indicator_idv;
	indicator_snp=cPar.indicator_snp;
	snpInfo=cPar.snpInfo;
	
	return;
}


void LDR::CopyToParam (PARAM &cPar) 
{
	//cPar.pheno_mean=pheno_mean;
	//cPar.randseed=randseed;
	
	return;
}


/*
void BSLMM::WriteBV (const gsl_vector *bv) 
{
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".bv.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	size_t t=0;
	for (size_t i=0; i<ni_total; ++i) {
		if (indicator_idv[i]==0) {
			outfile<<"NA"<<endl;
		}		
		else {
			outfile<<scientific<<setprecision(6)<<gsl_vector_get(bv, t)<<endl;
			t++;
		}
	}		
	
	outfile.clear();	
	outfile.close();	
	return;
}




void BSLMM::WriteParam (vector<pair<double, double> > &beta_g, const gsl_vector *alpha, const size_t w) 
{
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".param.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"chr"<<"\t"<<"rs"<<"\t"
			<<"ps"<<"\t"<<"n_miss"<<"\t"<<"alpha"<<"\t"
			<<"beta"<<"\t"<<"gamma"<<endl;
	
	size_t t=0;
	for (size_t i=0; i<ns_total; ++i) {
		if (indicator_snp[i]==0) {continue;}		
		
		outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"
		<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t";	
				
		outfile<<scientific<<setprecision(6)<<gsl_vector_get(alpha, t)<<"\t";
		if (beta_g[t].second!=0) {
			outfile<<beta_g[t].first/beta_g[t].second<<"\t"<<beta_g[t].second/(double)w<<endl;
		}
		else {
			outfile<<0.0<<"\t"<<0.0<<endl;
		}
		t++;
	}		
	
	outfile.clear();	
	outfile.close();	
	return;
}


void BSLMM::WriteParam (const gsl_vector *alpha) 
{
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".param.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"chr"<<"\t"<<"rs"<<"\t"
			<<"ps"<<"\t"<<"n_miss"<<"\t"<<"alpha"<<"\t"
			<<"beta"<<"\t"<<"gamma"<<endl;
	
	size_t t=0;
	for (size_t i=0; i<ns_total; ++i) {
		if (indicator_snp[i]==0) {continue;}		

		outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"
				<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t";				
		outfile<<scientific<<setprecision(6)<<gsl_vector_get(alpha, t)<<"\t";
		outfile<<0.0<<"\t"<<0.0<<endl;
		t++;
	}		
	
	outfile.clear();	
	outfile.close();	
	return;
}


void BSLMM::WriteResult (const int flag, const gsl_matrix *Result_hyp, const gsl_matrix *Result_gamma, const size_t w_col) 
{
	string file_gamma, file_hyp;
	file_gamma=path_out+"/"+file_out;
	file_gamma+=".gamma.txt";
	file_hyp=path_out+"/"+file_out;
	file_hyp+=".hyp.txt";

	ofstream outfile_gamma, outfile_hyp;
		
	if (flag==0) {
		outfile_gamma.open (file_gamma.c_str(), ofstream::out);
		outfile_hyp.open (file_hyp.c_str(), ofstream::out);
		if (!outfile_gamma) {cout<<"error writing file: "<<file_gamma<<endl; return;}
		if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}
		
		outfile_hyp<<"h \t pve \t rho \t pge \t pi \t n_gamma"<<endl;
		
		for (size_t i=0; i<s_max; ++i) {
			outfile_gamma<<"s"<<i<<"\t";
		}
		outfile_gamma<<endl;
	}
	else {
		outfile_gamma.open (file_gamma.c_str(), ofstream::app);
		outfile_hyp.open (file_hyp.c_str(), ofstream::app);
		if (!outfile_gamma) {cout<<"error writing file: "<<file_gamma<<endl; return;}
		if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}
		
		size_t w;
		if (w_col==0) {w=w_pace;}
		else {w=w_col;}
		
		for (size_t i=0; i<w; ++i) {
			outfile_hyp<<scientific;
			for (size_t j=0; j<4; ++j) {
				outfile_hyp<<setprecision(6)<<gsl_matrix_get (Result_hyp, i, j)<<"\t";
			}
			outfile_hyp<<setprecision(6)<<exp(gsl_matrix_get (Result_hyp, i, 4))<<"\t";
			outfile_hyp<<(int)gsl_matrix_get (Result_hyp, i, 5)<<"\t";
			outfile_hyp<<endl;
		}
		
		for (size_t i=0; i<w; ++i) {
			for (size_t j=0; j<s_max; ++j) {
				outfile_gamma<<(int)gsl_matrix_get (Result_gamma, i, j)<<"\t";
			}
			outfile_gamma<<endl;
		}
		
	}
	
	outfile_hyp.close();
	outfile_hyp.clear();
	outfile_gamma.close();
	outfile_gamma.clear();	
	return;
}

*/

//X is a p by n matrix
void LDR::VB (const vector<vector<unsigned char> > &Xt, const gsl_matrix *W_gsl, const gsl_vector *y_gsl)
{ 
  //save gsl_vector and gsl_matrix into eigen library formats
  MatrixXd W(W_gsl->size1, W_gsl->size2);
  VectorXd y(y_gsl->size);
  VectorXd x_col(y_gsl->size);

  double d;
  for (size_t i=0; i<W_gsl->size1; i++) {
    d=gsl_vector_get(y_gsl, i);
    y(i)=d;
    for (size_t j=0; j<W_gsl->size2; j++) {
      W(i,j)=gsl_matrix_get(W_gsl, i, j);
    }
  }

  //initial VB values by lm
  cout<<indicator_snp[0]<<" "<<indicator_snp[1]<<" "<<indicator_snp[2]<<endl;
  uchar_matrix_get_row (Xt, 0, x_col);

  for (size_t j=0; j<10; j++) {
    cout<<x_col(j)<<endl;
  }


  //run VB iterations



  //save results

  return;
}
