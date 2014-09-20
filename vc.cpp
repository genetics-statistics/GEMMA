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
#include <bitset>
#include <cstring>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"

#include "gsl/gsl_cdf.h"
#include "gsl/gsl_multiroots.h"
#include "gsl/gsl_min.h"

#include "io.h"
#include "lapack.h"
#include "gzstream.h"

#ifdef FORCE_FLOAT
#include "lmm_float.h"
#include "vc_float.h"
#else
#include "lmm.h"
#include "vc.h"
#endif



using namespace std;


//in this file, X, Y are already transformed (i.e. UtX and UtY)


void VC::CopyFromParam (PARAM &cPar) 
{	
	file_out=cPar.file_out;
	
	//	v_sigma2=cPar.v_sigma2;
	
	time_UtX=0.0;
	time_opt=0.0;

	v_traceG=cPar.v_traceG;
	
	return;
}


void VC::CopyToParam (PARAM &cPar) 
{
	cPar.time_UtX=time_UtX;
	cPar.time_opt=time_opt;	
		
	cPar.v_sigma2=v_sigma2;
	cPar.v_se_sigma2=v_se_sigma2;
	cPar.v_pve=v_pve;
	cPar.v_se_pve=v_se_pve;
	cPar.v_traceG=v_traceG;
	
	cPar.v_beta=v_beta;
	cPar.v_se_beta=v_se_beta;
	
	return;
}



void UpdateParam (const gsl_vector *log_sigma2, VC_PARAM *p)
{
  size_t n1=(p->K)->size1, n_vc=log_sigma2->size-1, n_cvt=(p->W)->size2;
 
  gsl_matrix *K_temp=gsl_matrix_alloc(n1, n1);
  gsl_matrix *HiW=gsl_matrix_alloc(n1, n_cvt);
  gsl_matrix *WtHiW=gsl_matrix_alloc(n_cvt, n_cvt);
  gsl_matrix *WtHiWi=gsl_matrix_alloc(n_cvt, n_cvt);
  gsl_matrix *WtHiWiWtHi=gsl_matrix_alloc(n_cvt, n1);

  double sigma2;  
  //calculate H=\sum_i^{k+1} \sigma_i^2 K_i
  gsl_matrix_set_zero (p->P);
  for (size_t i=0; i<n_vc+1; i++) {
    if (i==n_vc) {
      gsl_matrix_set_identity (K_temp);      
    } else {
      gsl_matrix_const_view K_sub=gsl_matrix_const_submatrix (p->K, 0, n1*i, n1, n1);
      gsl_matrix_memcpy (K_temp, &K_sub.matrix);
    }

    sigma2=exp(gsl_vector_get (log_sigma2, i) );
    gsl_matrix_scale(K_temp, sigma2);
    gsl_matrix_add (p->P, K_temp);
  }

  //calculate H^{-1}
  int sig;
  gsl_permutation * pmt1=gsl_permutation_alloc (n1);
  LUDecomp (p->P, pmt1, &sig);	
  LUInvert (p->P, pmt1, K_temp);
  gsl_permutation_free(pmt1);

  gsl_matrix_memcpy (p->P, K_temp);

  //calculate P=H^{-1}-H^{-1}W(W^TH^{-1}W)^{-1}W^TH^{-1}
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, p->P, p->W, 0.0, HiW);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, p->W, HiW, 0.0, WtHiW);

  gsl_permutation * pmt2=gsl_permutation_alloc (n_cvt);
  LUDecomp (WtHiW, pmt2, &sig);	
  LUInvert (WtHiW, pmt2, WtHiWi);
  gsl_permutation_free(pmt2);

  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, WtHiWi, HiW, 0.0, WtHiWiWtHi);  
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, -1.0, HiW, WtHiWiWtHi, 1.0, p->P);
  
  //calculate Py, KPy, PKPy
  gsl_blas_dgemv(CblasNoTrans, 1.0, p->P, p->y, 0.0, p->Py);    

  for (size_t i=0; i<n_vc+1; i++) {
    gsl_vector_view KPy=gsl_matrix_column (p->KPy_mat, i);
    gsl_vector_view PKPy=gsl_matrix_column (p->PKPy_mat, i);

    if (i==n_vc) {
      gsl_vector_memcpy (&KPy.vector, p->Py);
    } else {
      gsl_matrix_const_view K_sub=gsl_matrix_const_submatrix (p->K, 0, n1*i, n1, n1);      
      gsl_blas_dgemv(CblasNoTrans, 1.0, &K_sub.matrix, p->Py, 0.0, &KPy.vector);
    }
    
    gsl_blas_dgemv(CblasNoTrans, 1.0, p->P, &KPy.vector, 0.0, &PKPy.vector);
  }

  gsl_matrix_free (K_temp);
  gsl_matrix_free (HiW);
  gsl_matrix_free (WtHiW);
  gsl_matrix_free (WtHiWi);
  gsl_matrix_free (WtHiWiWtHi);

  return;
}


//below are functions for AI algorithm
int LogRL_dev1 (const gsl_vector *log_sigma2, void *params, gsl_vector *dev1)
{
  VC_PARAM *p=(VC_PARAM *) params;

  size_t n1=(p->K)->size1, n_vc=log_sigma2->size-1;
  
  double tr, d;

  //update parameters
  UpdateParam (log_sigma2, p);

  //calculate dev1=-0.5*trace(PK_i)+0.5*yPKPy
  for (size_t i=0; i<n_vc+1; i++) {
    if (i==n_vc) {
      tr=0;
      for (size_t l=0; l<n1; l++) {
	tr+=gsl_matrix_get (p->P, l, l);
      }
    } else {
      tr=0;
      for (size_t l=0; l<n1; l++) {
	gsl_vector_view P_row=gsl_matrix_row (p->P, l);
	gsl_vector_const_view K_col=gsl_matrix_const_column (p->K, n1*i+l);
	gsl_blas_ddot(&P_row.vector, &K_col.vector, &d);
	tr+=d;
      }
    }

    gsl_vector_view KPy_i=gsl_matrix_column (p->KPy_mat, i);
    gsl_blas_ddot(p->Py, &KPy_i.vector, &d);

    d=(-0.5*tr+0.5*d)*exp(gsl_vector_get(log_sigma2, i));
    
    gsl_vector_set(dev1, i, d);
  }

  return GSL_SUCCESS;
}



int LogRL_dev2 (const gsl_vector *log_sigma2, void *params, gsl_matrix *dev2)
{
  VC_PARAM *p=(VC_PARAM *) params;

  size_t n_vc=log_sigma2->size-1;
  
  double d, sigma2_i, sigma2_j;

  //update parameters
  UpdateParam (log_sigma2, p);
  
  //calculate dev2=0.5(yPKPKPy)
  for (size_t i=0; i<n_vc+1; i++) {
    gsl_vector_view KPy_i=gsl_matrix_column (p->KPy_mat, i);
    sigma2_i=exp(gsl_vector_get(log_sigma2, i));

    for (size_t j=i; j<n_vc+1; j++) {
      gsl_vector_view PKPy_j=gsl_matrix_column (p->PKPy_mat, j);

      gsl_blas_ddot(&KPy_i.vector, &PKPy_j.vector, &d);
      sigma2_j=exp(gsl_vector_get(log_sigma2, j));

      d*=-0.5*sigma2_i*sigma2_j;

      gsl_matrix_set(dev2, i, j, d);
      if (j!=i) {gsl_matrix_set(dev2, j, i, d);}
    }   
  }

  gsl_matrix_memcpy (p->Hessian, dev2);

  return GSL_SUCCESS;
}



int LogRL_dev12 (const gsl_vector *log_sigma2, void *params, gsl_vector *dev1, gsl_matrix *dev2)
{
  VC_PARAM *p=(VC_PARAM *) params;

  size_t n1=(p->K)->size1, n_vc=log_sigma2->size-1;
  
  double tr, d, sigma2_i, sigma2_j;

  //update parameters
  UpdateParam (log_sigma2, p);

  //calculate dev1=-0.5*trace(PK_i)+0.5*yPKPy
  //calculate dev2=0.5(yPKPKPy)
  for (size_t i=0; i<n_vc+1; i++) {
    if (i==n_vc) {
      tr=0;
      for (size_t l=0; l<n1; l++) {
	tr+=gsl_matrix_get (p->P, l, l);
      }
    } else {
      tr=0;
      for (size_t l=0; l<n1; l++) {
	gsl_vector_view P_row=gsl_matrix_row (p->P, l);
	gsl_vector_const_view K_col=gsl_matrix_const_column (p->K, n1*i+l);
	gsl_blas_ddot(&P_row.vector, &K_col.vector, &d);
	tr+=d;
      }
    }

    gsl_vector_view KPy_i=gsl_matrix_column (p->KPy_mat, i);
    gsl_blas_ddot(p->Py, &KPy_i.vector, &d);

    sigma2_i=exp(gsl_vector_get(log_sigma2, i));
    d=(-0.5*tr+0.5*d)*sigma2_i;
 
    gsl_vector_set(dev1, i, d);
      
    for (size_t j=i; j<n_vc+1; j++) {
      gsl_vector_view PKPy_j=gsl_matrix_column (p->PKPy_mat, j);
      gsl_blas_ddot(&KPy_i.vector, &PKPy_j.vector, &d);

      sigma2_j=exp(gsl_vector_get(log_sigma2, j));
      d*=-0.5*sigma2_i*sigma2_j;

      gsl_matrix_set(dev2, i, j, d);
      if (j!=i) {gsl_matrix_set(dev2, j, i, d);}
    }   

  }

  gsl_matrix_memcpy (p->Hessian, dev2);

  return GSL_SUCCESS;
}




void VC::CalcVCreml (const gsl_matrix *K, const gsl_matrix *W, const gsl_vector *y)
{
  size_t n1=K->size1, n2=K->size2;
  size_t n_vc=n2/n1;
  gsl_vector *log_sigma2=gsl_vector_alloc (n_vc+1);
  double d, s;

  //set up params
  gsl_matrix *P=gsl_matrix_alloc (n1, n1);
  gsl_vector *Py=gsl_vector_alloc (n1);
  gsl_matrix *KPy_mat=gsl_matrix_alloc (n1, n_vc+1);
  gsl_matrix *PKPy_mat=gsl_matrix_alloc (n1, n_vc+1);
  gsl_vector *dev1=gsl_vector_alloc (n_vc+1);
  gsl_matrix *dev2=gsl_matrix_alloc (n_vc+1, n_vc+1);
  gsl_matrix *Hessian=gsl_matrix_alloc (n_vc+1, n_vc+1);
  VC_PARAM params={K, W, y, P, Py, KPy_mat, PKPy_mat, Hessian};

  //initialize sigma2/log_sigma2
  gsl_blas_ddot (y, y, &s);
  s/=(double)n1;
  for (size_t i=0; i<n_vc+1; i++) {
    if (i==n_vc) {
      d=s/((double)n_vc+1.0);
    } else {
      d=s/( ((double)n_vc+1.0)*v_traceG[i]);
    }

    gsl_vector_set (log_sigma2, i, d);
  }
  //  gsl_vector_set (log_sigma2, 0, 0.38);
  //  gsl_vector_set (log_sigma2, 1, -1.08);

  cout<<"iteration "<<0<<endl;
  cout<<"sigma2 = ";
  for (size_t i=0; i<n_vc+1; i++) {
    cout<<exp(gsl_vector_get(log_sigma2, i))<<" ";
  }
  cout<<endl;

  //set up fdf
  gsl_multiroot_function_fdf FDF;
  FDF.n=n_vc+1;
  FDF.params=&params;
  FDF.f=&LogRL_dev1;
  FDF.df=&LogRL_dev2;
  FDF.fdf=&LogRL_dev12;
  
  //set up solver 	
  int status;
  int iter=0, max_iter=100;

  const gsl_multiroot_fdfsolver_type *T_fdf;
  gsl_multiroot_fdfsolver *s_fdf;
  T_fdf=gsl_multiroot_fdfsolver_hybridsj;
  s_fdf=gsl_multiroot_fdfsolver_alloc (T_fdf, n_vc+1);	

  gsl_multiroot_fdfsolver_set (s_fdf, &FDF, log_sigma2);

  do {
    iter++;
    status=gsl_multiroot_fdfsolver_iterate (s_fdf);

    if (status) break;

    cout<<"iteration "<<iter<<endl;
    cout<<"sigma2 = ";
    for (size_t i=0; i<n_vc+1; i++) {
      cout<<exp(gsl_vector_get(s_fdf->x, i))<<" ";
    }
    cout<<endl;
    cout<<"derivatives = ";
    for (size_t i=0; i<n_vc+1; i++) {
      cout<<gsl_vector_get(s_fdf->f, i)<<" ";
    }
    cout<<endl;

    status=gsl_multiroot_test_residual (s_fdf->f, 1e-3);		
  }
  while (status==GSL_CONTINUE && iter<max_iter); 

  //obtain Hessian inverse
  int sig=LogRL_dev12 (s_fdf->f, &params, dev1, dev2);

  gsl_permutation * pmt=gsl_permutation_alloc (n_vc+1);
  LUDecomp (dev2, pmt, &sig);	
  LUInvert (dev2, pmt, Hessian);
  gsl_permutation_free(pmt);

  //save data
  v_sigma2.clear(); 
  for (size_t i=0; i<n_vc+1; i++) {
    d=exp(gsl_vector_get(s_fdf->x, i));
    v_sigma2.push_back(d);
  }

  v_se_sigma2.clear();
  for (size_t i=0; i<n_vc+1; i++) {
    d=-1.0*v_sigma2[i]*v_sigma2[i]*gsl_matrix_get(Hessian, i, i);
    v_se_sigma2.push_back(sqrt(d));
  }

  s=0;
  for (size_t i=0; i<n_vc; i++) {
    s+=v_traceG[i]*v_sigma2[i];
  }
  s+=v_sigma2[n_vc];
  
  v_pve.clear();
  for (size_t i=0; i<n_vc; i++) {
    d=v_traceG[i]*v_sigma2[i]/s;
    v_pve.push_back(d);
  }

  v_se_pve.clear();
  for (size_t i=0; i<n_vc; i++) {
    d=v_traceG[i]*(s-v_sigma2[i]*v_traceG[i])/(s*s)*v_se_sigma2[i]*v_se_sigma2[i];
    v_se_pve.push_back(sqrt(d) );
  }
  
  gsl_multiroot_fdfsolver_free(s_fdf);	

  gsl_vector_free(log_sigma2);
  gsl_matrix_free(P);
  gsl_vector_free(Py);
  gsl_matrix_free(KPy_mat);
  gsl_matrix_free(PKPy_mat);
  gsl_vector_free(dev1);
  gsl_matrix_free(dev2);
  gsl_matrix_free(Hessian);

  return;
}


	



