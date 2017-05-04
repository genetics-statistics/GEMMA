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



#include "logistic.h"
#include "lapack.h"

#ifdef FORCE_FLOAT
#include "param_float.h"
#include "bslmmdap_float.h"
#include "lmm_float.h"  //for class FUNC_PARAM and MatrixCalcLR
#include "lm_float.h"
#include "mathfunc_float.h"  //for function CenterVector
#else
#include "param.h"
#include "bslmmdap.h"
#include "lmm.h"
#include "lm.h"
#include "mathfunc.h"
#endif

using namespace std;




void BSLMMDAP::CopyFromParam (PARAM &cPar)
{
	file_out=cPar.file_out;
	path_out=cPar.path_out;

	time_UtZ=0.0;
	time_Omega=0.0;

	h_min=cPar.h_min;
	h_max=cPar.h_max;
	h_ngrid=cPar.h_ngrid;
	rho_min=cPar.rho_min;
	rho_max=cPar.rho_max;
	rho_ngrid=cPar.rho_ngrid;

	if (h_min<=0) {h_min=0.01;}
	if (h_max>=1) {h_max=0.99;}
	if (rho_min<=0) {rho_min=0.01;}
	if (rho_max>=1) {rho_max=0.99;}

	trace_G=cPar.trace_G;

	ni_total=cPar.ni_total;
	ns_total=cPar.ns_total;
	ni_test=cPar.ni_test;
	ns_test=cPar.ns_test;

	indicator_idv=cPar.indicator_idv;
	indicator_snp=cPar.indicator_snp;
	snpInfo=cPar.snpInfo;

	return;
}


void BSLMMDAP::CopyToParam (PARAM &cPar)
{
	cPar.time_UtZ=time_UtZ;
	cPar.time_Omega=time_Omega;

	return;
}



//read hyp file
void ReadFile_hyb (const string &file_hyp, vector<double> &vec_sa2, vector<double> &vec_sb2, vector<double> &vec_wab)
{
  vec_sa2.clear(); vec_sb2.clear(); vec_wab.clear();

  igzstream infile (file_hyp.c_str(), igzstream::in);
  if (!infile) {cout<<"error! fail to open hyp file: "<<file_hyp<<endl; return;}

  string line;
  char *ch_ptr;

  getline(infile, line);

  while (!safeGetline(infile, line).eof()) {
    ch_ptr=strtok ((char *)line.c_str(), " , \t");
    ch_ptr=strtok (NULL, " , \t");

    ch_ptr=strtok (NULL, " , \t");
    vec_sa2.push_back(atof(ch_ptr));

    ch_ptr=strtok (NULL, " , \t");
    vec_sb2.push_back(atof(ch_ptr));

    ch_ptr=strtok (NULL, " , \t");
    vec_wab.push_back(atof(ch_ptr));
  }

  infile.close();
  infile.clear();

  return;
}


//read bf file
void ReadFile_bf (const string &file_bf, vector<string> &vec_rs, vector<vector<vector<double> > > &BF)
{
  BF.clear(); vec_rs.clear();

  igzstream infile (file_bf.c_str(), igzstream::in);
  if (!infile) {cout<<"error! fail to open bf file: "<<file_bf<<endl; return;}

  string line, rs, block;
  vector<double> vec_bf;
  vector<vector<double> > mat_bf;
  char *ch_ptr;

  size_t bf_size, flag_block;

  getline(infile, line);

  size_t t=0;
  while (!safeGetline(infile, line).eof()) {
    flag_block=0;

    ch_ptr=strtok ((char *)line.c_str(), " , \t");
    rs=ch_ptr;
    vec_rs.push_back(rs);

    ch_ptr=strtok (NULL, " , \t");
    if (t==0) {
      block=ch_ptr;
    } else {
      if (strcmp(ch_ptr, block.c_str() )!=0) {
	flag_block=1;
	block=ch_ptr;
      }
    }

    ch_ptr=strtok (NULL, " , \t");
    while (ch_ptr!=NULL) {
      vec_bf.push_back(atof(ch_ptr));
      ch_ptr=strtok (NULL, " , \t");
    }

    if (t==0) {
      bf_size=vec_bf.size();
    } else {
      if (bf_size!=vec_bf.size()) {cout<<"error! unequal row size in bf file."<<endl;}
    }

    if (flag_block==0) {
      mat_bf.push_back(vec_bf);
    } else {
      BF.push_back(mat_bf);
      mat_bf.clear();
    }
    vec_bf.clear();

    t++;
  }

  infile.close();
  infile.clear();

  return;
}


//read category files
//read both continuous and discrete category file, record mapRS2catc
void ReadFile_cat (const string &file_cat, const vector<string> &vec_rs, gsl_matrix *Ac, gsl_matrix_int *Ad, gsl_vector_int *dlevel, size_t &kc, size_t &kd)
{
  igzstream infile (file_cat.c_str(), igzstream::in);
  if (!infile) {cout<<"error! fail to open category file: "<<file_cat<<endl; return;}

  string line;
  char *ch_ptr;

  string rs, chr, a1, a0, pos, cm;

  //read header
  HEADER header;
  !safeGetline(infile, line).eof();
  ReadHeader (line, header);

  //use the header to determine the number of categories
  kc=header.catc_col.size(); kd=header.catd_col.size();

  //set up storage and mapper
  map<string, vector<double> > mapRS2catc;
  map<string, vector<int> > mapRS2catd;
  vector<double> catc;
  vector<int> catd;

  //read the following lines to record mapRS2cat
  while (!safeGetline(infile, line).eof()) {
    ch_ptr=strtok ((char *)line.c_str(), " , \t");

    if (header.rs_col==0) {
      rs=chr+":"+pos;
    }

    catc.clear(); catd.clear();

    for (size_t i=0; i<header.coln; i++) {
      if (header.rs_col!=0 && header.rs_col==i+1) {
	rs=ch_ptr;
      } else if (header.chr_col!=0 && header.chr_col==i+1) {
	chr=ch_ptr;
      } else if (header.pos_col!=0 && header.pos_col==i+1) {
	pos=ch_ptr;
      } else if (header.cm_col!=0 && header.cm_col==i+1) {
	cm=ch_ptr;
      } else if (header.a1_col!=0 && header.a1_col==i+1) {
	a1=ch_ptr;
      } else if (header.a0_col!=0 && header.a0_col==i+1) {
	a0=ch_ptr;
      } else if (header.catc_col.size()!=0 && header.catc_col.count(i+1)!=0 ) {
	catc.push_back(atof(ch_ptr));
      } else if (header.catd_col.size()!=0 && header.catd_col.count(i+1)!=0 ) {
	catd.push_back(atoi(ch_ptr));
      } else {}

      ch_ptr=strtok (NULL, " , \t");
    }

    if (mapRS2catc.count(rs)==0 && kc>0) {mapRS2catc[rs]=catc;}
    if (mapRS2catd.count(rs)==0 && kd>0) {mapRS2catd[rs]=catd;}
  }

  //load into Ad and Ac
  if (kc>0) {
    Ac=gsl_matrix_alloc(vec_rs.size(), kc);
    for (size_t i=0; i<vec_rs.size(); i++) {
      if (mapRS2catc.count(vec_rs[i])!=0) {
	for (size_t j=0; j<kc; j++) {
	  gsl_matrix_set(Ac, i, j, mapRS2catc[vec_rs[i]][j]);
	}
      } else {
	for (size_t j=0; j<kc; j++) {
	  gsl_matrix_set(Ac, i, j, 0);
	}
      }
    }
  }

  if (kd>0) {
    Ad=gsl_matrix_int_alloc(vec_rs.size(), kd);

    for (size_t i=0; i<vec_rs.size(); i++) {
      if (mapRS2catd.count(vec_rs[i])!=0) {
	for (size_t j=0; j<kd; j++) {
	  gsl_matrix_int_set(Ad, i, j, mapRS2catd[vec_rs[i]][j]);
	}
      } else {
	for (size_t j=0; j<kd; j++) {
	  gsl_matrix_int_set(Ad, i, j, 0);
	}
      }
    }

    dlevel=gsl_vector_int_alloc(kd);
    map<int, int> rcd;
    int val;
    for (size_t j=0; j<kd; j++) {
      rcd.clear();
      for (size_t i=0; i<Ad->size1; i++) {
	val = gsl_matrix_int_get(Ad, i, j);
	rcd[val] = 1;
      }
      gsl_vector_int_set (dlevel, j, rcd.size());
    }
  }

  infile.clear();
  infile.close();

  return;
}








void BSLMMDAP::WriteResult (const gsl_matrix *Hyper, const gsl_matrix *BF)
{
  string file_bf, file_hyp;
	file_bf=path_out+"/"+file_out;
	file_bf+=".bf.txt";
	file_hyp=path_out+"/"+file_out;
	file_hyp+=".hyp.txt";

	ofstream outfile_bf, outfile_hyp;

	outfile_bf.open (file_bf.c_str(), ofstream::out);
	outfile_hyp.open (file_hyp.c_str(), ofstream::out);

	if (!outfile_bf) {cout<<"error writing file: "<<file_bf<<endl; return;}
	if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}

	outfile_hyp<<"h"<<"\t"<<"rho"<<"\t"<<"sa2"<<"\t"<<"sb2"<<"\t"<<"weight"<<endl;
	outfile_hyp<<scientific;
	for (size_t i=0; i<Hyper->size1; i++) {
	  for (size_t j=0; j<Hyper->size2; j++) {
	    outfile_hyp<<setprecision(6)<<gsl_matrix_get (Hyper, i, j)<<"\t";
	  }
	  outfile_hyp<<endl;
	}

	outfile_bf<<"chr"<<"\t"<<"rs"<<"\t"<<"ps"<<"\t"<<"n_miss";
	for (size_t i=0; i<BF->size2; i++) {
	  outfile_bf<<"\t"<<"BF"<<i+1;
	}
	outfile_bf<<endl;

	size_t t=0;
	for (size_t i=0; i<ns_total; ++i) {
	  if (indicator_snp[i]==0) {continue;}

	  outfile_bf<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"
		    <<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss;

	  outfile_bf<<scientific;
	  for (size_t j=0; j<BF->size2; j++) {
	    outfile_bf<<"\t"<<setprecision(6)<<gsl_matrix_get (BF, t, j);
	  }
	  outfile_bf<<endl;

	  t++;
	}

	outfile_hyp.close();
	outfile_hyp.clear();
	outfile_bf.close();
	outfile_bf.clear();
	return;
}



void BSLMMDAP::WriteResult (const vector<string> &vec_rs, const gsl_matrix *Hyper, const gsl_vector *pip, const gsl_vector *coef)
{
  string file_gamma, file_hyp, file_coef;
	file_gamma=path_out+"/"+file_out;
	file_gamma+=".gamma.txt";
	file_hyp=path_out+"/"+file_out;
	file_hyp+=".hyp.txt";
	file_coef=path_out+"/"+file_out;
	file_coef+=".coef.txt";

	ofstream outfile_gamma, outfile_hyp, outfile_coef;

	outfile_gamma.open (file_gamma.c_str(), ofstream::out);
	outfile_hyp.open (file_hyp.c_str(), ofstream::out);
	outfile_coef.open (file_coef.c_str(), ofstream::out);

	if (!outfile_gamma) {cout<<"error writing file: "<<file_gamma<<endl; return;}
	if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}
	if (!outfile_coef) {cout<<"error writing file: "<<file_coef<<endl; return;}

	outfile_hyp<<"h"<<"\t"<<"rho"<<"\t"<<"sa2"<<"\t"<<"sb2"<<"\t"<<"weight"<<endl;
	outfile_hyp<<scientific;
	for (size_t i=0; i<Hyper->size1; i++) {
	  for (size_t j=0; j<Hyper->size2; j++) {
	    outfile_hyp<<setprecision(6)<<gsl_matrix_get (Hyper, i, j)<<"\t";
	  }
	  outfile_hyp<<endl;
	}


	outfile_gamma<<"rs"<<"\t"<<"gamma"<<endl;
	for (size_t i=0; i<vec_rs.size(); ++i) {
	  outfile_gamma<<vec_rs[i]<<"\t"<<scientific<<setprecision(6)<<gsl_vector_get(pip, i)<<endl;
	}

	outfile_coef<<"coef"<<endl;
	outfile_coef<<scientific;
	for (size_t i=0; i<coef->size; i++) {
	  outfile_coef<<setprecision(6)<<gsl_vector_get (coef, i)<<endl;
	}

	outfile_coef.close();
	outfile_coef.clear();
	outfile_hyp.close();
	outfile_hyp.clear();
	outfile_gamma.close();
	outfile_gamma.clear();
	return;
}




/*
void BSLMMDAP::SetXgamma (gsl_matrix *Xgamma, const gsl_matrix *X, vector<size_t> &rank)
{
	size_t pos;
	for (size_t i=0; i<rank.size(); ++i) {
		pos=mapRank2pos[rank[i]];
		gsl_vector_view Xgamma_col=gsl_matrix_column (Xgamma, i);
		gsl_vector_const_view X_col=gsl_matrix_const_column (X, pos);
		gsl_vector_memcpy (&Xgamma_col.vector, &X_col.vector);
	}

	return;
}
*/

double BSLMMDAP::CalcMarginal (const gsl_vector *Uty, const gsl_vector *K_eval, const double sigma_b2, const double tau)
{
	gsl_vector *weight_Hi=gsl_vector_alloc (Uty->size);

	double logm=0.0;
	double d, uy, Hi_yy=0, logdet_H=0.0;
	for (size_t i=0; i<ni_test; ++i) {
		d=gsl_vector_get (K_eval, i)*sigma_b2;
		d=1.0/(d+1.0);
		gsl_vector_set (weight_Hi, i, d);

		logdet_H-=log(d);
		uy=gsl_vector_get (Uty, i);
		Hi_yy+=d*uy*uy;
	}

	//calculate likelihood
	logm=-0.5*logdet_H-0.5*tau*Hi_yy+0.5*log(tau)*(double)ni_test;

	gsl_vector_free (weight_Hi);

	return logm;
}


double BSLMMDAP::CalcMarginal (const gsl_matrix *UtXgamma, const gsl_vector *Uty, const gsl_vector *K_eval, const double sigma_a2, const double sigma_b2, const double tau)
{
  clock_t  time_start;
	double logm=0.0;
	double d, uy, P_yy=0, logdet_O=0.0, logdet_H=0.0;

	gsl_matrix *UtXgamma_eval=gsl_matrix_alloc (UtXgamma->size1, UtXgamma->size2);
	gsl_matrix *Omega=gsl_matrix_alloc (UtXgamma->size2, UtXgamma->size2);
	gsl_vector *XtHiy=gsl_vector_alloc (UtXgamma->size2);
	gsl_vector *beta_hat=gsl_vector_alloc (UtXgamma->size2);
	gsl_vector *weight_Hi=gsl_vector_alloc (UtXgamma->size1);

	gsl_matrix_memcpy (UtXgamma_eval, UtXgamma);

	logdet_H=0.0; P_yy=0.0;
	for (size_t i=0; i<ni_test; ++i) {
		gsl_vector_view UtXgamma_row=gsl_matrix_row (UtXgamma_eval, i);
		d=gsl_vector_get (K_eval, i)*sigma_b2;
		d=1.0/(d+1.0);
		gsl_vector_set (weight_Hi, i, d);

		logdet_H-=log(d);
		uy=gsl_vector_get (Uty, i);
		P_yy+=d*uy*uy;
		gsl_vector_scale (&UtXgamma_row.vector, d);
	}

	//calculate Omega
	gsl_matrix_set_identity (Omega);

	time_start=clock();
#ifdef WITH_LAPACK
	lapack_dgemm ((char *)"T", (char *)"N", sigma_a2, UtXgamma_eval, UtXgamma, 1.0, Omega);
#else
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, sigma_a2, UtXgamma_eval, UtXgamma, 1.0, Omega);
#endif
	time_Omega+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);

	//calculate beta_hat
	gsl_blas_dgemv (CblasTrans, 1.0, UtXgamma_eval, Uty, 0.0, XtHiy);

	logdet_O=CholeskySolve(Omega, XtHiy, beta_hat);

	gsl_vector_scale (beta_hat, sigma_a2);

	gsl_blas_ddot (XtHiy, beta_hat, &d);
	P_yy-=d;

	gsl_matrix_free (UtXgamma_eval);
	gsl_matrix_free (Omega);
	gsl_vector_free (XtHiy);
	gsl_vector_free (beta_hat);
	gsl_vector_free (weight_Hi);

	logm=-0.5*logdet_H-0.5*logdet_O-0.5*tau*P_yy+0.5*log(tau)*(double)ni_test;

	return logm;
}


double BSLMMDAP::CalcPrior (class HYPBSLMM &cHyp) {
  double logprior=0;
  logprior=((double)cHyp.n_gamma-1.0)*cHyp.logp+((double)ns_test-(double)cHyp.n_gamma)*log(1.0-exp(cHyp.logp));
  return logprior;
}


//where A is the ni_test by n_cat matrix of annotations
void BSLMMDAP::DAP_CalcBF (const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty, const gsl_vector *K_eval, const gsl_vector *y) {
	clock_t time_start;

	//set up BF
	double tau, h, rho, sigma_a2, sigma_b2, d;
	size_t ns_causal=10;
	size_t n_grid=h_ngrid*rho_ngrid;
	vector<double> vec_sa2, vec_sb2, logm_null;

	gsl_matrix *BF=gsl_matrix_alloc(ns_test, n_grid);
	gsl_matrix *Xgamma=gsl_matrix_alloc(ni_test, 1);
	gsl_matrix *Hyper=gsl_matrix_alloc(n_grid, 5);

	//compute tau by using yty
	gsl_blas_ddot (Uty, Uty, &tau);
	tau=(double)ni_test/tau;

	//set up grid values for sigma_a2 and sigma_b2 based on an approximately even grid for h and rho, and a fixed number of causals
	size_t ij=0;
	for (size_t i=0; i<h_ngrid; i++) {
	  h=h_min+(h_max-h_min)*(double)i/((double)h_ngrid-1);
	  for (size_t j=0; j<rho_ngrid; j++) {
	    rho=rho_min+(rho_max-rho_min)*(double)j/((double)rho_ngrid-1);

	    sigma_a2=h*rho/((1-h)*(double)ns_causal);
	    sigma_b2=h*(1.0-rho)/(trace_G*(1-h));

	    vec_sa2.push_back(sigma_a2);
	    vec_sb2.push_back(sigma_b2);
	    logm_null.push_back(CalcMarginal (Uty, K_eval, 0.0, tau));

	    gsl_matrix_set (Hyper, ij, 0, h);
	    gsl_matrix_set (Hyper, ij, 1, rho);
	    gsl_matrix_set (Hyper, ij, 2, sigma_a2);
	    gsl_matrix_set (Hyper, ij, 3, sigma_b2);
	    gsl_matrix_set (Hyper, ij, 4, 1/(double)n_grid);
	    ij++;
	  }
	}

	//compute BF factors
	time_start=clock();
	cout<<"Calculating BF..."<<endl;
	for (size_t t=0; t<ns_test; t++) {
	  gsl_vector_view Xgamma_col=gsl_matrix_column (Xgamma, 0);
	  gsl_vector_const_view X_col=gsl_matrix_const_column (UtX, t);
	  gsl_vector_memcpy (&Xgamma_col.vector, &X_col.vector);

	  for (size_t ij=0; ij<n_grid; ij++) {
	    sigma_a2=vec_sa2[ij];
	    sigma_b2=vec_sb2[ij];

	    d=CalcMarginal (Xgamma, Uty, K_eval, sigma_a2, sigma_b2, tau);
	    d-=logm_null[ij];
	    d=exp(d);

	    gsl_matrix_set(BF, t, ij, d);
	  }
	}
	time_Proposal=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);

	//save results
	WriteResult (Hyper, BF);

	//free matrices and vectors
	gsl_matrix_free(BF);
	gsl_matrix_free(Xgamma);
	gsl_matrix_free(Hyper);
	return;
}





void single_ct_regression(const gsl_matrix_int *Xd, const gsl_vector_int *dlevel, const gsl_vector *pip_vec, gsl_vector *coef, gsl_vector *prior_vec) {

  map<int,double> sum_pip;
  map<int,double> sum;

  int levels = gsl_vector_int_get(dlevel,0);

  for(int i=0;i<levels;i++){
    sum_pip[i] = sum[i] = 0;
  }

  for(int i=0;i<Xd->size1;i++){
    int cat = gsl_matrix_int_get(Xd,i,0);
    sum_pip[cat] += gsl_vector_get(pip_vec,i);
    sum[cat] += 1;
  }

  for(int i=0;i<Xd->size1;i++){
    int cat = gsl_matrix_int_get(Xd,i,0);
    gsl_vector_set(prior_vec,i,sum_pip[cat]/sum[cat]);
  }

  //double baseline=0;
  for(int i=0;i<levels;i++){
    double new_prior = sum_pip[i]/sum[i];
    //gsl_vector_set(coef, i, log(new_prior/(1-new_prior))-baseline);
    //if(i==0){
    //baseline = log(new_prior/(1-new_prior));
    //}
    gsl_vector_set(coef, i, log(new_prior/(1-new_prior)) );
  }

  return;
}




//where A is the ni_test by n_cat matrix of annotations
void BSLMMDAP::DAP_EstimateHyper (const size_t kc, const size_t kd, const vector<string> &vec_rs, const vector<double> &vec_sa2, const vector<double> &vec_sb2, const vector<double> &wab, const vector<vector<vector<double> > > &BF, gsl_matrix *Ac, gsl_matrix_int *Ad, gsl_vector_int *dlevel) {
	clock_t time_start;

	//set up BF
	double h, rho, sigma_a2, sigma_b2, d, s, logm, logm_save;
	size_t t1, t2;
	size_t n_grid=wab.size(), ns_test=vec_rs.size();

	gsl_vector *prior_vec=gsl_vector_alloc(ns_test);
	gsl_matrix *Hyper=gsl_matrix_alloc(n_grid, 5);
	gsl_vector *pip=gsl_vector_alloc(ns_test);
	gsl_vector *coef=gsl_vector_alloc(kc+kd+1);

	//perform the EM algorithm
	vector<double> vec_wab, vec_wab_new;

	//initial values
	for (size_t t=0; t<ns_test; t++) {
	  gsl_vector_set (prior_vec, t, (double)BF.size()/(double)ns_test);
	}
	for (size_t ij=0; ij<n_grid; ij++) {
	  vec_wab.push_back(wab[ij]);
	  vec_wab_new.push_back(wab[ij]);
	}

	//EM iteration
	size_t it=0;
	double dif=1;
	while (it<100 && dif>1e-3) {
	  //update E_gamma
	  t1=0, t2=0;
	  for (size_t b=0; b<BF.size(); b++) {
	    s=1;
	    for (size_t m=0; m<BF[b].size(); m++) {
	      d=0;
	      for (size_t ij=0; ij<n_grid; ij++) {
		d+=vec_wab_new[ij]*BF[b][m][ij];
	      }
	      d*=gsl_vector_get(prior_vec, t1)/(1-gsl_vector_get(prior_vec, t1));

	      gsl_vector_set(pip, t1, d);
	      s+=d;
	      t1++;
	    }

	    for (size_t m=0; m<BF[b].size(); m++) {
	      d=gsl_vector_get(pip, t2)/s;
	      gsl_vector_set(pip, t2, d);
	      t2++;
	    }
	  }

	  //update E_wab
	  s=0;
	  for (size_t ij=0; ij<n_grid; ij++) {
	    vec_wab_new[ij]=0;

	    t1=0;
	    for (size_t b=0; b<BF.size(); b++) {
	      d=1;
	      for (size_t m=0; m<BF[b].size(); m++) {
		d+=gsl_vector_get(prior_vec, t1)/(1-gsl_vector_get(prior_vec, t1))*vec_wab[ij]*BF[b][m][ij];
		t1++;
	      }
	      vec_wab_new[ij]+=log(d);
	    }

	    s=max(s, vec_wab_new[ij]);
	  }

	  d=0;
	  for (size_t ij=0; ij<n_grid; ij++) {
	    vec_wab_new[ij]=exp(vec_wab_new[ij]-s);
	    d+=vec_wab_new[ij];
	  }

	  for (size_t ij=0; ij<n_grid; ij++) {
	    vec_wab_new[ij]/=d;
	    //	    vec_wab[ij]=vec_wab_new[ij];
	  }

	  //update coef, and pi
	  if(kc==0 && kd==0){//no annotation
	    s=0;
	    for (size_t t=0; t<pip->size; t++) {
	      s+=gsl_vector_get(pip, t);
	    }
	    s=s/(double)pip->size;
	    for (size_t t=0; t<pip->size; t++) {
	      gsl_vector_set(prior_vec, t, s);
	    }

	    gsl_vector_set (coef, 0, log(s/(1-s)));
	  } else if(kc==0 && kd!=0){//only discrete annotations
	    if(kd == 1){
	      single_ct_regression(Ad, dlevel, pip, coef, prior_vec);
	    }else{
	      logistic_cat_fit(coef, Ad, dlevel, pip, 0, 0);
	      logistic_cat_pred(coef, Ad, dlevel, prior_vec);
	    }
	  } else if (kc!=0 && kd==0) {//only continuous annotations
	    logistic_cont_fit(coef, Ac, pip, 0, 0);
	    logistic_cont_pred(coef, Ac, prior_vec);
	  } else if (kc!=0 && kd!=0) {//both continuous and categorical annotations
	    logistic_mixed_fit(coef, Ad, dlevel, Ac, pip, 0, 0);
	    logistic_mixed_pred(coef, Ad, dlevel, Ac, prior_vec);
	  }

	  //compute marginal likelihood
	  logm=0;

	  t1=0;
	  for (size_t b=0; b<BF.size(); b++) {
	    d=1; s=0;
	    for (size_t m=0; m<BF[b].size(); m++) {
	      s+=log(1-gsl_vector_get(prior_vec, t1));
	      for (size_t ij=0; ij<n_grid; ij++) {
		d+=gsl_vector_get(prior_vec, t1)/(1-gsl_vector_get(prior_vec, t1))*vec_wab[ij]*BF[b][m][ij];
	      }
	    }
	    logm+=log(d)+s;
	    t1++;
	  }

	  if (it>0) {
	    dif=logm-logm_save;
	  }
	  logm_save=logm;
	  it++;

	  cout<<"iteration = "<<it<<"; marginal likelihood = "<<logm<<endl;
	}

	//update h and rho that correspond to w_ab
	for (size_t ij=0; ij<n_grid; ij++) {
	  sigma_a2=vec_sa2[ij];
	  sigma_b2=vec_sb2[ij];

	  d=exp(gsl_vector_get(coef, coef->size-1))/(1+exp(gsl_vector_get(coef, coef->size-1)));
	  h=(d*(double)ns_test*sigma_a2+1*sigma_b2)/(1+d*(double)ns_test*sigma_a2+1*sigma_b2);
	  rho=d*(double)ns_test*sigma_a2/(d*(double)ns_test*sigma_a2+1*sigma_b2);

	  gsl_matrix_set (Hyper, ij, 0, h);
	  gsl_matrix_set (Hyper, ij, 1, rho);
	  gsl_matrix_set (Hyper, ij, 2, sigma_a2);
	  gsl_matrix_set (Hyper, ij, 3, sigma_b2);
	  gsl_matrix_set (Hyper, ij, 4, vec_wab_new[ij]);
	}

	//obtain beta and alpha parameters


	//save results
	WriteResult (vec_rs, Hyper, pip, coef);

	//free matrices and vectors
	gsl_vector_free(prior_vec);
	gsl_matrix_free(Hyper);
	gsl_vector_free(pip);
	gsl_vector_free(coef);
	return;
}

/*
//readin the estimated hyper-parameters and perform fine mapping for each region
void BSLMM::DAP_FineMapping (const gsl_matrix *U, const gsl_matrix *UtX, const gsl_matrix *A, const gsl_vector *Uty, const gsl_vector *K_eval, const gsl_vector *y, gsl_matrix *Hyper, gsl_vector *alpha, gsl_vector *pip) {
	clock_t time_start;

	//two priority sets: S_1 contains all candidate causal SNPs; S_2 contains the prioritized combintion of them
	//two marginal probability sets: P_1 contains marginals for S_1; P_2 contains marginals for S_2;
	set<size_t> S1set, S2set;
	vector<size_t> S1vec;
	vector<set<size_t> > S2vec;
	vector<double> P1, P2;

	//calculate P0 (null) and P1 (for every SNP)



	//loop through the number of combinations
	for (size_t s=0; s<p; s++) {
	  //if (s==0), set up S_1: compute marginal of the null model, then compute P_1, then compute BF_1 and use them to select S_1; compute C_1



	  //if (s==1), set up S_2: compute pair-wise P_2 and use them to select S_2; compute C_2

	  //otherwise, match each combination of S_2 with each SNP from S_1, select into S_3; and replace S_2 with S_3; compute C_s


	  //stop when the stopping critieria are reached (if S_2 is empty; if t; if kappa); add the residual component R

	for (size_t t=0; t<total_step; ++t) {
		if (t%d_pace==0 || t==total_step-1) {ProgressBar ("Running MCMC ", t, total_step-1, (double)n_accept/(double)(t*n_mh+1));}
//		if (t>10) {break;}

		if (a_mode==13) {
			SampleZ (y, z_hat, z);
			mean_z=CenterVector (z);

			time_start=clock();
			gsl_blas_dgemv (CblasTrans, 1.0, U, z, 0.0, Utz);
			time_UtZ+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);

			//First proposal
			if (cHyp_old.n_gamma==0 || cHyp_old.rho==0) {
				logPost_old=CalcPosterior(Utz, K_eval, Utu_old, alpha_old, cHyp_old);
				beta_old.clear();
				for (size_t i=0; i<cHyp_old.n_gamma; ++i) {
				  beta_old.push_back(0);
				}
			}
			else {
				gsl_matrix *UtXgamma=gsl_matrix_alloc (ni_test, cHyp_old.n_gamma);
				gsl_vector *beta=gsl_vector_alloc (cHyp_old.n_gamma);
				SetXgamma (UtXgamma, UtX, rank_old);
				logPost_old=CalcPosterior(UtXgamma, Utz, K_eval, UtXb_old, Utu_old, alpha_old, beta, cHyp_old);

				beta_old.clear();
				for (size_t i=0; i<beta->size; ++i) {
					beta_old.push_back(gsl_vector_get(beta, i));
				}
				gsl_matrix_free (UtXgamma);
				gsl_vector_free (beta);
			}
		}


	delete [] p_gamma;
	beta_g.clear();

	return;
}

*/






/*
//below fits MCMC for rho=1
void BSLMM::CalcXtX (const gsl_matrix *X, const gsl_vector *y, const size_t s_size, gsl_matrix *XtX, gsl_vector *Xty)
{
  time_t time_start=clock();
  gsl_matrix_const_view X_sub=gsl_matrix_const_submatrix(X, 0, 0, X->size1, s_size);
  gsl_matrix_view XtX_sub=gsl_matrix_submatrix(XtX, 0, 0, s_size, s_size);
  gsl_vector_view Xty_sub=gsl_vector_subvector(Xty, 0, s_size);

#ifdef WITH_LAPACK
  lapack_dgemm ((char *)"T", (char *)"N", 1.0, &X_sub.matrix, &X_sub.matrix, 0.0, &XtX_sub.matrix);
#else
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, &X_sub.matrix, &X_sub.matrix, 0.0, &XtX_sub.matrix);
#endif
  gsl_blas_dgemv(CblasTrans, 1.0, &X_sub.matrix, y, 0.0, &Xty_sub.vector);

  time_Omega+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);

  return;
}



double BSLMM::CalcPosterior (const double yty, class HYPBSLMM &cHyp)
{
	double logpost=0.0;

	//for quantitative traits, calculate pve and pge
	//pve and pge for case/control data are calculted in CalcCC_PVEnZ
	if (a_mode==11) {
		cHyp.pve=0.0;
		cHyp.pge=1.0;
	}

	//calculate likelihood
	if (a_mode==11) {logpost-=0.5*(double)ni_test*log(yty);}
	else {logpost-=0.5*yty;}

	logpost+=((double)cHyp.n_gamma-1.0)*cHyp.logp+((double)ns_test-(double)cHyp.n_gamma)*log(1-exp(cHyp.logp));

	return logpost;
}


double BSLMM::CalcPosterior (const gsl_matrix *Xgamma, const gsl_matrix *XtX, const gsl_vector *Xty, const double yty, const size_t s_size, gsl_vector *Xb, gsl_vector *beta, class HYPBSLMM &cHyp)
{
	double sigma_a2=cHyp.h/( (1-cHyp.h)*exp(cHyp.logp)*(double)ns_test);
	double logpost=0.0;
	double d, P_yy=yty, logdet_O=0.0;

	gsl_matrix_const_view Xgamma_sub=gsl_matrix_const_submatrix (Xgamma, 0, 0, Xgamma->size1, s_size);
	gsl_matrix_const_view XtX_sub=gsl_matrix_const_submatrix (XtX, 0, 0, s_size, s_size);
	gsl_vector_const_view Xty_sub=gsl_vector_const_subvector (Xty, 0, s_size);

	gsl_matrix *Omega=gsl_matrix_alloc (s_size, s_size);
	gsl_matrix *M_temp=gsl_matrix_alloc (s_size, s_size);
	gsl_vector *beta_hat=gsl_vector_alloc (s_size);
	gsl_vector *Xty_temp=gsl_vector_alloc (s_size);

	gsl_vector_memcpy (Xty_temp, &Xty_sub.vector);

	//calculate Omega
	gsl_matrix_memcpy (Omega, &XtX_sub.matrix);
	gsl_matrix_scale (Omega, sigma_a2);
	gsl_matrix_set_identity (M_temp);
	gsl_matrix_add (Omega, M_temp);

	//calculate beta_hat
	logdet_O=CholeskySolve(Omega, Xty_temp, beta_hat);
	gsl_vector_scale (beta_hat, sigma_a2);

	gsl_blas_ddot (Xty_temp, beta_hat, &d);
	P_yy-=d;

	//sample tau
	double tau=1.0;
	if (a_mode==11) {tau =gsl_ran_gamma (gsl_r, (double)ni_test/2.0,  2.0/P_yy); }

	//sample beta
	for (size_t i=0; i<s_size; i++)
	{
		d=gsl_ran_gaussian(gsl_r, 1);
		gsl_vector_set(beta, i, d);
	}
	gsl_vector_view beta_sub=gsl_vector_subvector(beta, 0, s_size);
	gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Omega, &beta_sub.vector);

	//it compuates inv(L^T(Omega)) %*% beta;
	gsl_vector_scale(&beta_sub.vector, sqrt(sigma_a2/tau));
	gsl_vector_add(&beta_sub.vector, beta_hat);
	gsl_blas_dgemv (CblasNoTrans, 1.0, &Xgamma_sub.matrix, &beta_sub.vector, 0.0, Xb);

	//for quantitative traits, calculate pve and pge
	if (a_mode==11) {
		gsl_blas_ddot (Xb, Xb, &d);
		cHyp.pve=d/(double)ni_test;
		cHyp.pve/=cHyp.pve+1.0/tau;
		cHyp.pge=1.0;
	}

	logpost=-0.5*logdet_O;
	if (a_mode==11) {logpost-=0.5*(double)ni_test*log(P_yy);}
	else {logpost-=0.5*P_yy;}

	logpost+=((double)cHyp.n_gamma-1.0)*cHyp.logp+((double)ns_test-(double)cHyp.n_gamma)*log(1.0-exp(cHyp.logp));

	gsl_matrix_free (Omega);
	gsl_matrix_free (M_temp);
	gsl_vector_free (beta_hat);
	gsl_vector_free (Xty_temp);

	return logpost;
}
*/


