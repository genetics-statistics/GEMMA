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
#include "gsl/gsl_roots.h"
#include "gsl/gsl_min.h"
#include "gsl/gsl_integration.h"

#include "io.h"
#include "lapack.h"
#include "gzstream.h"

#ifdef FORCE_FLOAT
#include "lmm_float.h"
#else
#include "lmm.h"
#endif


using namespace std;





void LMM::CopyFromParam (PARAM &cPar) 
{
	a_mode=cPar.a_mode;
	d_pace=cPar.d_pace;
	
	file_bfile=cPar.file_bfile;
	file_geno=cPar.file_geno;
	file_out=cPar.file_out;
	path_out=cPar.path_out;
	file_gene=cPar.file_gene;
	// WJA added
	file_bgenfile=cPar.file_bgenfile;
	
	l_min=cPar.l_min;
	l_max=cPar.l_max;
	n_region=cPar.n_region;	
	l_mle_null=cPar.l_mle_null;
	logl_mle_H0=cPar.logl_mle_H0;
	
	time_UtX=0.0;
	time_opt=0.0;
	
	ni_total=cPar.ni_total;
	ns_total=cPar.ns_total;
	ni_test=cPar.ni_test;
	ns_test=cPar.ns_test;
	n_cvt=cPar.n_cvt;
	
	ng_total=cPar.ng_total;
	ng_test=0;
	
	indicator_idv=cPar.indicator_idv;	
	indicator_snp=cPar.indicator_snp;	
	snpInfo=cPar.snpInfo;
	
	return;
}


void LMM::CopyToParam (PARAM &cPar) 
{
	cPar.time_UtX=time_UtX;
	cPar.time_opt=time_opt;	
	
	cPar.ng_test=ng_test;
	
	return;
}



void LMM::WriteFiles () 
{
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".assoc.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

	if (!file_gene.empty()) {
		outfile<<"geneID"<<"\t";
		
		if (a_mode==1) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"l_remle"<<"\t"<<"p_wald"<<endl;
		} else if (a_mode==2) {
			outfile<<"l_mle"<<"\t"<<"p_lrt"<<endl;
		} else if (a_mode==3) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"p_score"<<endl;
		} else if (a_mode==4) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"l_remle"<<"\t"<<"l_mle"<<"\t"<<"p_wald"<<"\t"<<"p_lrt"<<"\t"<<"p_score"<<endl;
		} else {}
				
		for (vector<SUMSTAT>::size_type t=0; t<sumStat.size(); ++t) {	
			outfile<<snpInfo[t].rs_number<<"\t";
			
			if (a_mode==1) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].lambda_remle<<"\t"<<sumStat[t].p_wald <<endl;
			} else if (a_mode==2) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].lambda_mle<<"\t"<<sumStat[t].p_lrt<<endl;
			} else if (a_mode==3) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].p_score<<endl;
			} else if (a_mode==4) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].lambda_remle<<"\t"<<sumStat[t].lambda_mle<<"\t"<<sumStat[t].p_wald <<"\t"<<sumStat[t].p_lrt<<"\t"<<sumStat[t].p_score<<endl;
			} else {}
		}	
	}  else {
		outfile<<"chr"<<"\t"<<"rs"<<"\t"<<"ps"<<"\t"<<"n_miss"<<"\t"<<"allele1"<<"\t"<<"allele0"<<"\t"<<"af"<<"\t";
		
		if (a_mode==1) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"l_remle"<<"\t"<<"p_wald"<<endl;
		} else if (a_mode==2) {
			outfile<<"l_mle"<<"\t"<<"p_lrt"<<endl;
		} else if (a_mode==3) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"p_score"<<endl;
		} else if (a_mode==4) {
			outfile<<"beta"<<"\t"<<"se"<<"\t"<<"l_remle"<<"\t"<<"l_mle"<<"\t"<<"p_wald"<<"\t"<<"p_lrt"<<"\t"<<"p_score"<<endl;
		} else {}
		
		size_t t=0;
		for (size_t i=0; i<snpInfo.size(); ++i) {
			if (indicator_snp[i]==0) {continue;}
			
			outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t"<<snpInfo[i].a_minor<<"\t"<<snpInfo[i].a_major<<"\t"<<fixed<<setprecision(3)<<snpInfo[i].maf<<"\t";
			
			if (a_mode==1) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].lambda_remle<<"\t"<<sumStat[t].p_wald <<endl;
			} else if (a_mode==2) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].lambda_mle<<"\t"<<sumStat[t].p_lrt<<endl;
			} else if (a_mode==3) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].p_score<<endl;
			} else if (a_mode==4) {
				outfile<<scientific<<setprecision(6)<<sumStat[t].beta<<"\t"<<sumStat[t].se<<"\t"<<sumStat[t].lambda_remle<<"\t"<<sumStat[t].lambda_mle<<"\t"<<sumStat[t].p_wald <<"\t"<<sumStat[t].p_lrt<<"\t"<<sumStat[t].p_score<<endl;
			} else {}
			t++;
		}
	}
	
		
	outfile.close();
	outfile.clear();
	return;
}











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


void CalcPab (const size_t n_cvt, const size_t e_mode, const gsl_vector *Hi_eval, const gsl_matrix *Uab, const gsl_vector *ab, gsl_matrix *Pab)
{
	size_t index_ab, index_aw, index_bw, index_ww;
	double p_ab;
	double ps_ab, ps_aw, ps_bw, ps_ww;
	
	for (size_t p=0; p<=n_cvt+1; ++p) {
		for (size_t a=p+1; a<=n_cvt+2; ++a) {
			for (size_t b=a; b<=n_cvt+2; ++b) {
				index_ab=GetabIndex (a, b, n_cvt);
				if (p==0) {			
					gsl_vector_const_view Uab_col=gsl_matrix_const_column (Uab, index_ab);
					gsl_blas_ddot (Hi_eval, &Uab_col.vector, &p_ab);
					if (e_mode!=0) {p_ab=gsl_vector_get (ab, index_ab)-p_ab;}
					gsl_matrix_set (Pab, 0, index_ab, p_ab);
				}
				else {
					index_aw=GetabIndex (a, p, n_cvt);
					index_bw=GetabIndex (b, p, n_cvt);
					index_ww=GetabIndex (p, p, n_cvt);
					
					ps_ab=gsl_matrix_get (Pab, p-1, index_ab);
					ps_aw=gsl_matrix_get (Pab, p-1, index_aw);
					ps_bw=gsl_matrix_get (Pab, p-1, index_bw);
					ps_ww=gsl_matrix_get (Pab, p-1, index_ww);
					
					p_ab=ps_ab-ps_aw*ps_bw/ps_ww;
					gsl_matrix_set (Pab, p, index_ab, p_ab);
				}
			}
		}
	}
	return;
}


void CalcPPab (const size_t n_cvt, const size_t e_mode, const gsl_vector *HiHi_eval, const gsl_matrix *Uab, const gsl_vector *ab, const gsl_matrix *Pab, gsl_matrix *PPab)
{
	size_t index_ab, index_aw, index_bw, index_ww;
	double p2_ab;
	double ps2_ab, ps_aw, ps_bw, ps_ww, ps2_aw, ps2_bw, ps2_ww;
	
	for (size_t p=0; p<=n_cvt+1; ++p) {
		for (size_t a=p+1; a<=n_cvt+2; ++a) {
			for (size_t b=a; b<=n_cvt+2; ++b) {
				index_ab=GetabIndex (a, b, n_cvt);
				if (p==0) {					
					gsl_vector_const_view Uab_col=gsl_matrix_const_column (Uab, index_ab);
					gsl_blas_ddot (HiHi_eval, &Uab_col.vector, &p2_ab);
					if (e_mode!=0) {p2_ab=p2_ab-gsl_vector_get (ab, index_ab)+2.0*gsl_matrix_get (Pab, 0, index_ab);}
					gsl_matrix_set (PPab, 0, index_ab, p2_ab);
				}
				else {
					index_aw=GetabIndex (a, p, n_cvt);
					index_bw=GetabIndex (b, p, n_cvt);
					index_ww=GetabIndex (p, p, n_cvt);
					
					ps2_ab=gsl_matrix_get (PPab, p-1, index_ab);
					ps_aw=gsl_matrix_get (Pab, p-1, index_aw);
					ps_bw=gsl_matrix_get (Pab, p-1, index_bw);
					ps_ww=gsl_matrix_get (Pab, p-1, index_ww);
					ps2_aw=gsl_matrix_get (PPab, p-1, index_aw);
					ps2_bw=gsl_matrix_get (PPab, p-1, index_bw);
					ps2_ww=gsl_matrix_get (PPab, p-1, index_ww);
					
					p2_ab=ps2_ab+ps_aw*ps_bw*ps2_ww/(ps_ww*ps_ww);
					p2_ab-=(ps_aw*ps2_bw+ps_bw*ps2_aw)/ps_ww;
					gsl_matrix_set (PPab, p, index_ab, p2_ab);
					
				}
			}
		}
	}
	return;
}


void CalcPPPab (const size_t n_cvt, const size_t e_mode, const gsl_vector *HiHiHi_eval, const gsl_matrix *Uab, const gsl_vector *ab, const gsl_matrix *Pab, const gsl_matrix *PPab, gsl_matrix *PPPab)
{
	size_t index_ab, index_aw, index_bw, index_ww;
	double p3_ab;
	double ps3_ab, ps_aw, ps_bw, ps_ww, ps2_aw, ps2_bw, ps2_ww, ps3_aw, ps3_bw, ps3_ww;
	
	for (size_t p=0; p<=n_cvt+1; ++p) {
		for (size_t a=p+1; a<=n_cvt+2; ++a) {
			for (size_t b=a; b<=n_cvt+2; ++b) {
				index_ab=GetabIndex (a, b, n_cvt);
				if (p==0) {					
					gsl_vector_const_view Uab_col=gsl_matrix_const_column (Uab, index_ab);
					gsl_blas_ddot (HiHiHi_eval, &Uab_col.vector, &p3_ab);
					if (e_mode!=0) {p3_ab=gsl_vector_get (ab, index_ab)-p3_ab+3.0*gsl_matrix_get (PPab, 0, index_ab)-3.0*gsl_matrix_get (Pab, 0, index_ab);}
					gsl_matrix_set (PPPab, 0, index_ab, p3_ab);
				}
				else {
					index_aw=GetabIndex (a, p, n_cvt);
					index_bw=GetabIndex (b, p, n_cvt);
					index_ww=GetabIndex (p, p, n_cvt);
					
					ps3_ab=gsl_matrix_get (PPPab, p-1, index_ab);
					ps_aw=gsl_matrix_get (Pab, p-1, index_aw);
					ps_bw=gsl_matrix_get (Pab, p-1, index_bw);
					ps_ww=gsl_matrix_get (Pab, p-1, index_ww);
					ps2_aw=gsl_matrix_get (PPab, p-1, index_aw);
					ps2_bw=gsl_matrix_get (PPab, p-1, index_bw);
					ps2_ww=gsl_matrix_get (PPab, p-1, index_ww);
					ps3_aw=gsl_matrix_get (PPPab, p-1, index_aw);
					ps3_bw=gsl_matrix_get (PPPab, p-1, index_bw);
					ps3_ww=gsl_matrix_get (PPPab, p-1, index_ww);
					
					p3_ab=ps3_ab-ps_aw*ps_bw*ps2_ww*ps2_ww/(ps_ww*ps_ww*ps_ww);
					p3_ab-=(ps_aw*ps3_bw+ps_bw*ps3_aw+ps2_aw*ps2_bw)/ps_ww;
					p3_ab+=(ps_aw*ps2_bw*ps2_ww+ps_bw*ps2_aw*ps2_ww+ps_aw*ps_bw*ps3_ww)/(ps_ww*ps_ww);
					
					gsl_matrix_set (PPPab, p, index_ab, p3_ab);
				}
			}
		}
	}
	return;
}



double LogL_f (double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt;} else {nc_total=n_cvt+1;}
	
	double f=0.0, logdet_h=0.0, d;
	size_t index_yy;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
				
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);	
	
	for (size_t i=0; i<(p->eval)->size; ++i) {
		d=gsl_vector_get (v_temp, i);
		logdet_h+=log(fabs(d));
	}	
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	
	double c=0.5*(double)ni_test*(log((double)ni_test)-log(2*M_PI)-1.0);
	
	index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);	
	double P_yy=gsl_matrix_get (Pab, nc_total, index_yy);
	f=c-0.5*logdet_h-0.5*(double)ni_test*log(P_yy);
	
	gsl_matrix_free (Pab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (v_temp);
	return f;
}

 
 



double LogL_dev1 (double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;	
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt;} else {nc_total=n_cvt+1;}
	
	double dev1=0.0, trace_Hi=0.0;
	size_t index_yy;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
	
	gsl_vector_memcpy (HiHi_eval, Hi_eval);
	gsl_vector_mul (HiHi_eval, Hi_eval);	
	
	gsl_vector_set_all (v_temp, 1.0);
	gsl_blas_ddot (Hi_eval, v_temp, &trace_Hi);
		
	if (p->e_mode!=0) {trace_Hi=(double)ni_test-trace_Hi;}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	CalcPPab (n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);	
	
	double trace_HiK=((double)ni_test-trace_Hi)/l;	
	
	index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);
	
	double P_yy=gsl_matrix_get (Pab, nc_total, index_yy);
	double PP_yy=gsl_matrix_get (PPab, nc_total, index_yy);
	double yPKPy=(P_yy-PP_yy)/l;	
	dev1=-0.5*trace_HiK+0.5*(double)ni_test*yPKPy/P_yy;
			
	gsl_matrix_free (Pab);
	gsl_matrix_free (PPab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (HiHi_eval);
	gsl_vector_free (v_temp);	
	
	return dev1;
}
	
	


double LogL_dev2 (double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;	
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt;} else {nc_total=n_cvt+1;}
	
	double dev2=0.0, trace_Hi=0.0, trace_HiHi=0.0;
	size_t index_yy;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
		
	gsl_vector_memcpy (HiHi_eval, Hi_eval);
	gsl_vector_mul (HiHi_eval, Hi_eval);	
	gsl_vector_memcpy (HiHiHi_eval, HiHi_eval);
	gsl_vector_mul (HiHiHi_eval, Hi_eval);
	
	gsl_vector_set_all (v_temp, 1.0);
	gsl_blas_ddot (Hi_eval, v_temp, &trace_Hi);
	gsl_blas_ddot (HiHi_eval, v_temp, &trace_HiHi);
	
	if (p->e_mode!=0) {		
		trace_Hi=(double)ni_test-trace_Hi;
		trace_HiHi=2*trace_Hi+trace_HiHi-(double)ni_test;
	}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	CalcPPab (n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);	
	CalcPPPab (n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, Pab, PPab, PPPab);	
	
	double trace_HiKHiK=((double)ni_test+trace_HiHi-2*trace_Hi)/(l*l);
	
	index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);
	double P_yy=gsl_matrix_get (Pab, nc_total, index_yy);
	double PP_yy=gsl_matrix_get (PPab, nc_total, index_yy);
	double PPP_yy=gsl_matrix_get (PPPab, nc_total, index_yy);		
		
	double yPKPy=(P_yy-PP_yy)/l;
	double yPKPKPy=(P_yy+PPP_yy-2.0*PP_yy)/(l*l);
		
	dev2=0.5*trace_HiKHiK-0.5*(double)ni_test*(2.0*yPKPKPy*P_yy-yPKPy*yPKPy)/(P_yy*P_yy);
		
	gsl_matrix_free (Pab);
	gsl_matrix_free (PPab);
	gsl_matrix_free (PPPab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (HiHi_eval);
	gsl_vector_free (HiHiHi_eval);
	gsl_vector_free (v_temp);	
	
	return dev2;
}
	
	
	
	
	
void LogL_dev12 (double l, void *params, double *dev1, double *dev2)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt;} else {nc_total=n_cvt+1;}
	
	double trace_Hi=0.0, trace_HiHi=0.0;
	size_t index_yy;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
		
	gsl_vector_memcpy (HiHi_eval, Hi_eval);
	gsl_vector_mul (HiHi_eval, Hi_eval);	
	gsl_vector_memcpy (HiHiHi_eval, HiHi_eval);
	gsl_vector_mul (HiHiHi_eval, Hi_eval);
	
	gsl_vector_set_all (v_temp, 1.0);
	gsl_blas_ddot (Hi_eval, v_temp, &trace_Hi);
	gsl_blas_ddot (HiHi_eval, v_temp, &trace_HiHi);
	
	if (p->e_mode!=0) {	
		trace_Hi=(double)ni_test-trace_Hi;
		trace_HiHi=2*trace_Hi+trace_HiHi-(double)ni_test;
	}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	CalcPPab (n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);	
	CalcPPPab (n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, Pab, PPab, PPPab);	
	
	double trace_HiK=((double)ni_test-trace_Hi)/l;
	double trace_HiKHiK=((double)ni_test+trace_HiHi-2*trace_Hi)/(l*l);
	
	index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);
	
	double P_yy=gsl_matrix_get (Pab, nc_total, index_yy);
	double PP_yy=gsl_matrix_get (PPab, nc_total, index_yy);
	double PPP_yy=gsl_matrix_get (PPPab, nc_total, index_yy);		
		
	double yPKPy=(P_yy-PP_yy)/l;	
	double yPKPKPy=(P_yy+PPP_yy-2.0*PP_yy)/(l*l);
		
	*dev1=-0.5*trace_HiK+0.5*(double)ni_test*yPKPy/P_yy;
	*dev2=0.5*trace_HiKHiK-0.5*(double)ni_test*(2.0*yPKPKPy*P_yy-yPKPy*yPKPy)/(P_yy*P_yy);
			
	gsl_matrix_free (Pab);
	gsl_matrix_free (PPab);
	gsl_matrix_free (PPPab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (HiHi_eval);
	gsl_vector_free (HiHiHi_eval);
	gsl_vector_free (v_temp);	
	
	return;
}



double LogRL_f (double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;	
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	double df;
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt; df=(double)ni_test-(double)n_cvt; }
	else {nc_total=n_cvt+1; df=(double)ni_test-(double)n_cvt-1.0;}
	
	double f=0.0, logdet_h=0.0, logdet_hiw=0.0, d;
	size_t index_ww;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *Iab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);	
	
	for (size_t i=0; i<(p->eval)->size; ++i) {
		d=gsl_vector_get (v_temp, i);
		logdet_h+=log(fabs(d));
	}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	gsl_vector_set_all (v_temp, 1.0);
	CalcPab (n_cvt, p->e_mode, v_temp, p->Uab, p->ab, Iab);	
	
	//calculate |WHiW|-|WW|
	logdet_hiw=0.0;
	for (size_t i=0; i<nc_total; ++i) {
		index_ww=GetabIndex (i+1, i+1, n_cvt);
		d=gsl_matrix_get (Pab, i, index_ww);
		logdet_hiw+=log(d);
		d=gsl_matrix_get (Iab, i, index_ww);
		logdet_hiw-=log(d);
	}
	index_ww=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);	
	double P_yy=gsl_matrix_get (Pab, nc_total, index_ww);
	
	double c=0.5*df*(log(df)-log(2*M_PI)-1.0);		
	f=c-0.5*logdet_h-0.5*logdet_hiw-0.5*df*log(P_yy);
		
	gsl_matrix_free (Pab);
	gsl_matrix_free (Iab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (v_temp);
	return f;
}



double LogRL_dev1 (double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;	
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	double df;
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt; df=(double)ni_test-(double)n_cvt; }
	else {nc_total=n_cvt+1; df=(double)ni_test-(double)n_cvt-1.0;}
	
	double dev1=0.0, trace_Hi=0.0;
	size_t index_ww;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
	
	gsl_vector_memcpy (HiHi_eval, Hi_eval);
	gsl_vector_mul (HiHi_eval, Hi_eval);	
	
	gsl_vector_set_all (v_temp, 1.0);
	gsl_blas_ddot (Hi_eval, v_temp, &trace_Hi);
	
	if (p->e_mode!=0) {	
		trace_Hi=(double)ni_test-trace_Hi;
	}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	CalcPPab (n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);	
	
	//calculate tracePK and trace PKPK
	double trace_P=trace_Hi;
	double ps_ww, ps2_ww;
	for (size_t i=0; i<nc_total; ++i) {
		index_ww=GetabIndex (i+1, i+1, n_cvt);
		ps_ww=gsl_matrix_get (Pab, i, index_ww);
		ps2_ww=gsl_matrix_get (PPab, i, index_ww);
		trace_P-=ps2_ww/ps_ww;
	}
	double trace_PK=(df-trace_P)/l;
	
	//calculate yPKPy, yPKPKPy
	index_ww=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);
	double P_yy=gsl_matrix_get (Pab, nc_total, index_ww);
	double PP_yy=gsl_matrix_get (PPab, nc_total, index_ww);		
	double yPKPy=(P_yy-PP_yy)/l;	
	
	dev1=-0.5*trace_PK+0.5*df*yPKPy/P_yy;	
			
	gsl_matrix_free (Pab);
	gsl_matrix_free (PPab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (HiHi_eval);
	gsl_vector_free (v_temp);	
	
	return dev1;
}




double LogRL_dev2 (double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;	
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	double df;
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt; df=(double)ni_test-(double)n_cvt; }
	else {nc_total=n_cvt+1; df=(double)ni_test-(double)n_cvt-1.0;}
	
	double dev2=0.0, trace_Hi=0.0, trace_HiHi=0.0;
	size_t index_ww;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
		
	gsl_vector_memcpy (HiHi_eval, Hi_eval);
	gsl_vector_mul (HiHi_eval, Hi_eval);	
	gsl_vector_memcpy (HiHiHi_eval, HiHi_eval);
	gsl_vector_mul (HiHiHi_eval, Hi_eval);
	
	gsl_vector_set_all (v_temp, 1.0);
	gsl_blas_ddot (Hi_eval, v_temp, &trace_Hi);
	gsl_blas_ddot (HiHi_eval, v_temp, &trace_HiHi);
	
	if (p->e_mode!=0) {	
		trace_Hi=(double)ni_test-trace_Hi;
		trace_HiHi=2*trace_Hi+trace_HiHi-(double)ni_test;
	}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	CalcPPab (n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);	
	CalcPPPab (n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, Pab, PPab, PPPab);	
	
	//calculate tracePK and trace PKPK
	double trace_P=trace_Hi, trace_PP=trace_HiHi;
	double ps_ww, ps2_ww, ps3_ww;
	for (size_t i=0; i<nc_total; ++i) {
		index_ww=GetabIndex (i+1, i+1, n_cvt);
		ps_ww=gsl_matrix_get (Pab, i, index_ww);
		ps2_ww=gsl_matrix_get (PPab, i, index_ww);
		ps3_ww=gsl_matrix_get (PPPab, i, index_ww);
		trace_P-=ps2_ww/ps_ww;
		trace_PP+=ps2_ww*ps2_ww/(ps_ww*ps_ww)-2.0*ps3_ww/ps_ww;
	}
	double trace_PKPK=(df+trace_PP-2.0*trace_P)/(l*l);
	
	//calculate yPKPy, yPKPKPy
	index_ww=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);
	double P_yy=gsl_matrix_get (Pab, nc_total, index_ww);
	double PP_yy=gsl_matrix_get (PPab, nc_total, index_ww);
	double PPP_yy=gsl_matrix_get (PPPab, nc_total, index_ww);				
	double yPKPy=(P_yy-PP_yy)/l;	
	double yPKPKPy=(P_yy+PPP_yy-2.0*PP_yy)/(l*l);
	
	dev2=0.5*trace_PKPK-0.5*df*(2.0*yPKPKPy*P_yy-yPKPy*yPKPy)/(P_yy*P_yy);
	
	gsl_matrix_free (Pab);
	gsl_matrix_free (PPab);
	gsl_matrix_free (PPPab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (HiHi_eval);
	gsl_vector_free (HiHiHi_eval);
	gsl_vector_free (v_temp);	
	
	return dev2;
}
	



void LogRL_dev12 (double l, void *params, double *dev1, double *dev2)
{
	FUNC_PARAM *p=(FUNC_PARAM *) params;	
	size_t n_cvt=p->n_cvt;
	size_t ni_test=p->ni_test;	
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	double df;
	size_t nc_total;
	if (p->calc_null==true) {nc_total=n_cvt; df=(double)ni_test-(double)n_cvt; }
	else {nc_total=n_cvt+1; df=(double)ni_test-(double)n_cvt-1.0;}
	
	double trace_Hi=0.0, trace_HiHi=0.0;
	size_t index_ww;
	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_matrix *PPPab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *HiHiHi_eval=gsl_vector_alloc((p->eval)->size);
	gsl_vector *v_temp=gsl_vector_alloc((p->eval)->size);
	
	gsl_vector_memcpy (v_temp, p->eval);
	gsl_vector_scale (v_temp, l);
	if (p->e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
		
	gsl_vector_memcpy (HiHi_eval, Hi_eval);
	gsl_vector_mul (HiHi_eval, Hi_eval);	
	gsl_vector_memcpy (HiHiHi_eval, HiHi_eval);
	gsl_vector_mul (HiHiHi_eval, Hi_eval);
	
	gsl_vector_set_all (v_temp, 1.0);
	gsl_blas_ddot (Hi_eval, v_temp, &trace_Hi);
	gsl_blas_ddot (HiHi_eval, v_temp, &trace_HiHi);
	
	if (p->e_mode!=0) {	
		trace_Hi=(double)ni_test-trace_Hi;
		trace_HiHi=2*trace_Hi+trace_HiHi-(double)ni_test;
	}
	
	CalcPab (n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, Pab);	
	CalcPPab (n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, Pab, PPab);	
	CalcPPPab (n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, Pab, PPab, PPPab);	
		
	//calculate tracePK and trace PKPK
	double trace_P=trace_Hi, trace_PP=trace_HiHi;
	double ps_ww, ps2_ww, ps3_ww;
	for (size_t i=0; i<nc_total; ++i) {
		index_ww=GetabIndex (i+1, i+1, n_cvt);
		ps_ww=gsl_matrix_get (Pab, i, index_ww);
		ps2_ww=gsl_matrix_get (PPab, i, index_ww);
		ps3_ww=gsl_matrix_get (PPPab, i, index_ww);
		trace_P-=ps2_ww/ps_ww;
		trace_PP+=ps2_ww*ps2_ww/(ps_ww*ps_ww)-2.0*ps3_ww/ps_ww;
	}
	double trace_PK=(df-trace_P)/l;
	double trace_PKPK=(df+trace_PP-2.0*trace_P)/(l*l);
	
	//calculate yPKPy, yPKPKPy
	index_ww=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);
	double P_yy=gsl_matrix_get (Pab, nc_total, index_ww);
	double PP_yy=gsl_matrix_get (PPab, nc_total, index_ww);
	double PPP_yy=gsl_matrix_get (PPPab, nc_total, index_ww);				
	double yPKPy=(P_yy-PP_yy)/l;	
	double yPKPKPy=(P_yy+PPP_yy-2.0*PP_yy)/(l*l);
	
	*dev1=-0.5*trace_PK+0.5*df*yPKPy/P_yy;
	*dev2=0.5*trace_PKPK-0.5*df*(2.0*yPKPKPy*P_yy-yPKPy*yPKPy)/(P_yy*P_yy);
	
	gsl_matrix_free (Pab);
	gsl_matrix_free (PPab);
	gsl_matrix_free (PPPab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (HiHi_eval);
	gsl_vector_free (HiHiHi_eval);
	gsl_vector_free (v_temp);	
	
	return ;
}
	







void LMM::CalcRLWald (const double &l, const FUNC_PARAM &params, double &beta, double &se, double &p_wald)
{
	size_t n_cvt=params.n_cvt;
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	int df=(int)ni_test-(int)n_cvt-1;
			
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc(params.eval->size);
	gsl_vector *v_temp=gsl_vector_alloc(params.eval->size);
	
	gsl_vector_memcpy (v_temp, params.eval);
	gsl_vector_scale (v_temp, l);
	if (params.e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);	
	
	CalcPab (n_cvt, params.e_mode, Hi_eval, params.Uab, params.ab, Pab);	
	
	size_t index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);	
	size_t index_xx=GetabIndex (n_cvt+1, n_cvt+1, n_cvt);
	size_t index_xy=GetabIndex (n_cvt+2, n_cvt+1, n_cvt);
	double P_yy=gsl_matrix_get (Pab, n_cvt, index_yy);
	double P_xx=gsl_matrix_get (Pab, n_cvt, index_xx);
	double P_xy=gsl_matrix_get (Pab, n_cvt, index_xy);	
	double Px_yy=gsl_matrix_get (Pab, n_cvt+1, index_yy);	
	
	beta=P_xy/P_xx;
	double tau=(double)df/Px_yy;
	se=sqrt(1.0/(tau*P_xx));	
	p_wald=gsl_cdf_fdist_Q ((P_yy-Px_yy)*tau, 1.0, df);	
//	p_wald=gsl_cdf_chisq_Q ((P_yy-Px_yy)*tau, 1);	
	
	gsl_matrix_free (Pab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (v_temp);
	return ;
}


void LMM::CalcRLScore (const double &l, const FUNC_PARAM &params, double &beta, double &se, double &p_score)
{
	size_t n_cvt=params.n_cvt;
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	int df=(int)ni_test-(int)n_cvt-1;
			
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc(params.eval->size);
	gsl_vector *v_temp=gsl_vector_alloc(params.eval->size);
	
	gsl_vector_memcpy (v_temp, params.eval);
	gsl_vector_scale (v_temp, l);
	if (params.e_mode==0) {gsl_vector_set_all (Hi_eval, 1.0);} else {gsl_vector_memcpy (Hi_eval, v_temp);}
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);	
	
	CalcPab (n_cvt, params.e_mode, Hi_eval, params.Uab, params.ab, Pab);	
	
	size_t index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);	
	size_t index_xx=GetabIndex (n_cvt+1, n_cvt+1, n_cvt);
	size_t index_xy=GetabIndex (n_cvt+2, n_cvt+1, n_cvt);
	double P_yy=gsl_matrix_get (Pab, n_cvt, index_yy);
	double P_xx=gsl_matrix_get (Pab, n_cvt, index_xx);
	double P_xy=gsl_matrix_get (Pab, n_cvt, index_xy);	
	double Px_yy=gsl_matrix_get (Pab, n_cvt+1, index_yy);	
	
	beta=P_xy/P_xx;
	double tau=(double)df/Px_yy;
	se=sqrt(1.0/(tau*P_xx));	
	
	p_score=gsl_cdf_fdist_Q ((double)ni_test*P_xy*P_xy/(P_yy*P_xx), 1.0, df);
//	p_score=gsl_cdf_chisq_Q ((double)ni_test*P_xy*P_xy/(P_yy*P_xx), 1);	
	
	gsl_matrix_free (Pab);
	gsl_vector_free (Hi_eval);
	gsl_vector_free (v_temp);
	return ;
}








void CalcUab (const gsl_matrix *UtW, const gsl_vector *Uty, gsl_matrix *Uab) 
{
	size_t index_ab;
	size_t n_cvt=UtW->size2;
	
	gsl_vector *u_a=gsl_vector_alloc (Uty->size);
	
	for (size_t a=1; a<=n_cvt+2; ++a) {
		if (a==n_cvt+1) {continue;}
		
		if (a==n_cvt+2) {gsl_vector_memcpy (u_a, Uty);}
		else {
			gsl_vector_const_view UtW_col=gsl_matrix_const_column (UtW, a-1);
			gsl_vector_memcpy (u_a, &UtW_col.vector);
		}
		
		for (size_t b=a; b>=1; --b) {		
			if (b==n_cvt+1) {continue;}
			
			index_ab=GetabIndex (a, b, n_cvt);
			gsl_vector_view Uab_col=gsl_matrix_column (Uab, index_ab);
			
			if (b==n_cvt+2) {gsl_vector_memcpy (&Uab_col.vector, Uty);}
			else {
				gsl_vector_const_view UtW_col=gsl_matrix_const_column (UtW, b-1);
				gsl_vector_memcpy (&Uab_col.vector, &UtW_col.vector);
			}			
			
			gsl_vector_mul(&Uab_col.vector, u_a);
		}
	}
	
	gsl_vector_free (u_a);
	return;
}


void CalcUab (const gsl_matrix *UtW, const gsl_vector *Uty, const gsl_vector *Utx, gsl_matrix *Uab) 
{	
	size_t index_ab;
	size_t n_cvt=UtW->size2;
	
	for (size_t b=1; b<=n_cvt+2; ++b) {			
		index_ab=GetabIndex (n_cvt+1, b, n_cvt);
		gsl_vector_view Uab_col=gsl_matrix_column (Uab, index_ab);
		
		if (b==n_cvt+2) {gsl_vector_memcpy (&Uab_col.vector, Uty);}
		else if (b==n_cvt+1) {gsl_vector_memcpy (&Uab_col.vector, Utx);}
		else {
			gsl_vector_const_view UtW_col=gsl_matrix_const_column (UtW, b-1);
			gsl_vector_memcpy (&Uab_col.vector, &UtW_col.vector);
		}
		
		gsl_vector_mul(&Uab_col.vector, Utx);
	}
	
	return;
}



void Calcab (const gsl_matrix *W, const gsl_vector *y, gsl_vector *ab) 
{
	size_t index_ab;
	size_t n_cvt=W->size2;
	
	double d;
	gsl_vector *v_a=gsl_vector_alloc (y->size);
	gsl_vector *v_b=gsl_vector_alloc (y->size);
	
	for (size_t a=1; a<=n_cvt+2; ++a) {
		if (a==n_cvt+1) {continue;}
		
		if (a==n_cvt+2) {gsl_vector_memcpy (v_a, y);}
		else {
			gsl_vector_const_view W_col=gsl_matrix_const_column (W, a-1);
			gsl_vector_memcpy (v_a, &W_col.vector);
		}
		
		for (size_t b=a; b>=1; --b) {		
			if (b==n_cvt+1) {continue;}
			
			index_ab=GetabIndex (a, b, n_cvt);
			
			if (b==n_cvt+2) {gsl_vector_memcpy (v_b, y);}
			else {
				gsl_vector_const_view W_col=gsl_matrix_const_column (W, b-1);
				gsl_vector_memcpy (v_b, &W_col.vector);
			}			
			
			gsl_blas_ddot (v_a, v_b, &d);
			gsl_vector_set(ab, index_ab, d);
		}
	}
	
	gsl_vector_free (v_a);
	gsl_vector_free (v_b);
	return;
}


void Calcab (const gsl_matrix *W, const gsl_vector *y, const gsl_vector *x, gsl_vector *ab) 
{	
	size_t index_ab;
	size_t n_cvt=W->size2;
	
	double d;
	gsl_vector *v_b=gsl_vector_alloc (y->size);
	
	for (size_t b=1; b<=n_cvt+2; ++b) {			
		index_ab=GetabIndex (n_cvt+1, b, n_cvt);
		
		if (b==n_cvt+2) {gsl_vector_memcpy (v_b, y);}
		else if (b==n_cvt+1) {gsl_vector_memcpy (v_b, x);}
		else {
			gsl_vector_const_view W_col=gsl_matrix_const_column (W, b-1);
			gsl_vector_memcpy (v_b, &W_col.vector);
		}
		
		gsl_blas_ddot (x, v_b, &d);
		gsl_vector_set(ab, index_ab, d);
	}
	
	gsl_vector_free (v_b);
	
	return;
}





void LMM::AnalyzeGene (const gsl_matrix *U, const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Utx, const gsl_matrix *W, const gsl_vector *x) 
{
	ifstream infile (file_gene.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading gene expression file:"<<file_gene<<endl; return;}
	
	clock_t time_start=clock();
	
	string line;
	char *ch_ptr;
	
	double lambda_mle=0, lambda_remle=0, beta=0, se=0, p_wald=0, p_lrt=0, p_score=0;
	double logl_H1=0.0, logl_H0=0.0, l_H0;
	int c_phen;
	string rs; //gene id
	double d;
	
	//Calculate basic quantities
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	gsl_vector *y=gsl_vector_alloc (U->size1);
	gsl_vector *Uty=gsl_vector_alloc (U->size2);
	gsl_matrix *Uab=gsl_matrix_alloc (U->size2, n_index);
	gsl_vector *ab=gsl_vector_alloc (n_index);	
		
	//header
	getline(infile, line);
	
	for (size_t t=0; t<ng_total; t++) {
		!safeGetline(infile, line).eof();
		if (t%d_pace==0 || t==ng_total-1) {ProgressBar ("Performing Analysis ", t, ng_total-1);}
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		rs=ch_ptr;
		
		c_phen=0; 
		for (size_t i=0; i<indicator_idv.size(); ++i) {
			ch_ptr=strtok (NULL, " , \t");
			if (indicator_idv[i]==0) {continue;}
			
			d=atof(ch_ptr); 			
			gsl_vector_set(y, c_phen, d);
			
			c_phen++;
		}
		
		time_start=clock();
		gsl_blas_dgemv (CblasTrans, 1.0, U, y, 0.0, Uty);		
		time_UtX+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
	
		//calculate null
		time_start=clock();
		
		gsl_matrix_set_zero (Uab);
		
		CalcUab (UtW, Uty, Uab);
		FUNC_PARAM param0={false, ni_test, n_cvt, eval, Uab, ab, 0};
		
		if (a_mode==2 || a_mode==3 || a_mode==4) {
			CalcLambda('L', param0, l_min, l_max, n_region, l_H0, logl_H0);
		}
		
		//calculate alternative
		CalcUab(UtW, Uty, Utx, Uab);
		FUNC_PARAM param1={false, ni_test, n_cvt, eval, Uab, ab, 0};
		
		//3 is before 1
		if (a_mode==3 || a_mode==4) {
			CalcRLScore (l_H0, param1, beta, se, p_score);
		}
		
		if (a_mode==1 || a_mode==4) {
			CalcLambda ('R', param1, l_min, l_max, n_region, lambda_remle, logl_H1);
			CalcRLWald (lambda_remle, param1, beta, se, p_wald);
		}
		
		if (a_mode==2 || a_mode==4) {
			CalcLambda ('L', param1, l_min, l_max, n_region, lambda_mle, logl_H1);
			p_lrt=gsl_cdf_chisq_Q (2.0*(logl_H1-logl_H0), 1);	
		}
		
		time_opt+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		//store summary data
		SUMSTAT SNPs={beta, se, lambda_remle, lambda_mle, p_wald, p_lrt, p_score};
		sumStat.push_back(SNPs);
    }
	cout<<endl;
	
	gsl_vector_free (y);
	gsl_vector_free (Uty);
	gsl_matrix_free (Uab);
	gsl_vector_free (ab);
	
	infile.close();
	infile.clear();
	
	return;
}





void LMM::AnalyzeBimbam (const gsl_matrix *U, const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Uty, const gsl_matrix *W, const gsl_vector *y) 
{
	igzstream infile (file_geno.c_str(), igzstream::in);
//	ifstream infile (file_geno.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return;}

	clock_t time_start=clock();
	
	string line;
	char *ch_ptr;
	
	double lambda_mle=0, lambda_remle=0, beta=0, se=0, p_wald=0, p_lrt=0, p_score=0;
	double logl_H1=0.0;
	int n_miss, c_phen;
	double geno, x_mean;
	
	//Calculate basic quantities
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;

	gsl_vector *x=gsl_vector_alloc (U->size1);
	gsl_vector *x_miss=gsl_vector_alloc (U->size1);
	gsl_vector *Utx=gsl_vector_alloc (U->size2);
	gsl_matrix *Uab=gsl_matrix_alloc (U->size2, n_index);
	gsl_vector *ab=gsl_vector_alloc (n_index);	
	
	gsl_matrix_set_zero (Uab);
	CalcUab (UtW, Uty, Uab);
//	if (e_mode!=0) {
//		gsl_vector_set_zero (ab);
//		Calcab (W, y, ab);
//	}	
	
	//start reading genotypes and analyze	
	for (size_t t=0; t<indicator_snp.size(); ++t) {
//		if (t>1) {break;}
		!safeGetline(infile, line).eof();
		if (t%d_pace==0 || t==(ns_total-1)) {ProgressBar ("Reading SNPs  ", t, ns_total-1);}
		if (indicator_snp[t]==0) {continue;}
		
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		ch_ptr=strtok (NULL, " , \t");
		ch_ptr=strtok (NULL, " , \t");		
		
		x_mean=0.0; c_phen=0; n_miss=0;
		gsl_vector_set_zero(x_miss);
		for (size_t i=0; i<ni_total; ++i) {
			ch_ptr=strtok (NULL, " , \t");
			if (indicator_idv[i]==0) {continue;}
			
			if (strcmp(ch_ptr, "NA")==0) {gsl_vector_set(x_miss, c_phen, 0.0); n_miss++;}
			else {
				geno=atof(ch_ptr); 				
				
				gsl_vector_set(x, c_phen, geno); 
				gsl_vector_set(x_miss, c_phen, 1.0); 
				x_mean+=geno;
			}
			c_phen++;
		}	
		
		x_mean/=(double)(ni_test-n_miss);
		
		for (size_t i=0; i<ni_test; ++i) {
			if (gsl_vector_get (x_miss, i)==0) {gsl_vector_set(x, i, x_mean);}
			geno=gsl_vector_get(x, i);
			if (x_mean>1) {
				gsl_vector_set(x, i, 2-geno);
			}
		}
		
		
		//calculate statistics
		time_start=clock();
		gsl_blas_dgemv (CblasTrans, 1.0, U, x, 0.0, Utx);		
		time_UtX+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		CalcUab(UtW, Uty, Utx, Uab);
//		if (e_mode!=0) {
//			Calcab (W, y, x, ab);
//		}
		
		time_start=clock();
		FUNC_PARAM param1={false, ni_test, n_cvt, eval, Uab, ab, 0};
		
		//3 is before 1
		if (a_mode==3 || a_mode==4) {
			CalcRLScore (l_mle_null, param1, beta, se, p_score);
		}
		
		if (a_mode==1 || a_mode==4) {
			CalcLambda ('R', param1, l_min, l_max, n_region, lambda_remle, logl_H1);	
			CalcRLWald (lambda_remle, param1, beta, se, p_wald);
		}
		
		if (a_mode==2 || a_mode==4) {
			CalcLambda ('L', param1, l_min, l_max, n_region, lambda_mle, logl_H1);
			p_lrt=gsl_cdf_chisq_Q (2.0*(logl_H1-logl_mle_H0), 1);	
		}			
		
		if (x_mean>1) {beta*=-1;}
		
		time_opt+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		//store summary data
		SUMSTAT SNPs={beta, se, lambda_remle, lambda_mle, p_wald, p_lrt, p_score};
		sumStat.push_back(SNPs);
    }	
	cout<<endl;
	
	gsl_vector_free (x);
	gsl_vector_free (x_miss);
	gsl_vector_free (Utx);
	gsl_matrix_free (Uab);
	gsl_vector_free (ab);
	
	infile.close();
	infile.clear();
	
	return;
}







void LMM::AnalyzePlink (const gsl_matrix *U, const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Uty, const gsl_matrix *W, const gsl_vector *y) 
{
	string file_bed=file_bfile+".bed";
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return;}
	
	clock_t time_start=clock();
	
	char ch[1];
	bitset<8> b;	
	
	double lambda_mle=0, lambda_remle=0, beta=0, se=0, p_wald=0, p_lrt=0, p_score=0;
	double logl_H1=0.0;
	int n_bit, n_miss, ci_total, ci_test;
	double geno, x_mean;
		
	//Calculate basic quantities
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;

	gsl_vector *x=gsl_vector_alloc (U->size1);
	gsl_vector *Utx=gsl_vector_alloc (U->size2);
	gsl_matrix *Uab=gsl_matrix_alloc (U->size2, n_index);	
	gsl_vector *ab=gsl_vector_alloc (n_index);	
	
	gsl_matrix_set_zero (Uab);
	CalcUab (UtW, Uty, Uab);
//	if (e_mode!=0) {
//		gsl_vector_set_zero (ab);
//		Calcab (W, y, ab);
//	}
		
	//calculate n_bit and c, the number of bit for each snp
	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1; }

	//print the first three majic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}
	
	
	for (vector<SNPINFO>::size_type t=0; t<snpInfo.size(); ++t) {
		if (t%d_pace==0 || t==snpInfo.size()-1) {ProgressBar ("Reading SNPs  ", t, snpInfo.size()-1);}
		if (indicator_snp[t]==0) {continue;}
		
		infile.seekg(t*n_bit+3);		//n_bit, and 3 is the number of magic numbers
		
		//read genotypes
		x_mean=0.0;	n_miss=0; ci_total=0; ci_test=0; 
		for (int i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && ci_total==(int)ni_total) {break;}
				if (indicator_idv[ci_total]==0) {ci_total++; continue;}

				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(x, ci_test, 2); x_mean+=2.0; }
					else {gsl_vector_set(x, ci_test, 1); x_mean+=1.0; }
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(x, ci_test, 0); }                                  
					else {gsl_vector_set(x, ci_test, -9); n_miss++; }
				}

				ci_total++;
				ci_test++;
			}
		}
		
		x_mean/=(double)(ni_test-n_miss);
				
		for (size_t i=0; i<ni_test; ++i) {			
			geno=gsl_vector_get(x,i);
			if (geno==-9) {gsl_vector_set(x, i, x_mean); geno=x_mean;}
			if (x_mean>1) {
				gsl_vector_set(x, i, 2-geno);
			}
		}
		
		//calculate statistics
		time_start=clock();
		gsl_blas_dgemv (CblasTrans, 1.0, U, x, 0.0, Utx);
		time_UtX+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		CalcUab(UtW, Uty, Utx, Uab);
//		if (e_mode!=0) {
//			Calcab (W, y, x, ab);
//		}
		
		time_start=clock();
		FUNC_PARAM param1={false, ni_test, n_cvt, eval, Uab, ab, 0};
		
		//3 is before 1, for beta
		if (a_mode==3 || a_mode==4) {
			CalcRLScore (l_mle_null, param1, beta, se, p_score);
		}
		
		if (a_mode==1 || a_mode==4) {
			CalcLambda ('R', param1, l_min, l_max, n_region, lambda_remle, logl_H1);	
			CalcRLWald (lambda_remle, param1, beta, se, p_wald);
		}
		
		if (a_mode==2 || a_mode==4) {
			CalcLambda ('L', param1, l_min, l_max, n_region, lambda_mle, logl_H1);
			p_lrt=gsl_cdf_chisq_Q (2.0*(logl_H1-logl_mle_H0), 1);	
		}		
		
		if (x_mean>1) {beta*=-1;}		
		
		time_opt+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		//store summary data
		SUMSTAT SNPs={beta, se, lambda_remle, lambda_mle, p_wald, p_lrt, p_score};
		sumStat.push_back(SNPs);
    }	
	cout<<endl;
	
	gsl_vector_free (x);
	gsl_vector_free (Utx);
	gsl_matrix_free (Uab);
	gsl_vector_free (ab);
	
	infile.close();
	infile.clear();	
	
	return;
}

// WJA added
#include <assert.h>
void LMM::AnalyzeBGEN (const gsl_matrix *U, const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Uty, const gsl_matrix *W, const gsl_vector *y) 
{
	string file_bgen=file_bgenfile+".bgen";
	ifstream infile (file_bgen.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bgen file:"<<file_bgen<<endl; return;}


	clock_t time_start=clock();
	
	string line;
	
	double lambda_mle=0, lambda_remle=0, beta=0, se=0, p_wald=0, p_lrt=0, p_score=0;
	double logl_H1=0.0;
	int n_miss, c_phen;
	double geno, x_mean;
	
	//Calculate basic quantities
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;

	gsl_vector *x=gsl_vector_alloc (U->size1);
	gsl_vector *x_miss=gsl_vector_alloc (U->size1);
	gsl_vector *Utx=gsl_vector_alloc (U->size2);
	gsl_matrix *Uab=gsl_matrix_alloc (U->size2, n_index);
	gsl_vector *ab=gsl_vector_alloc (n_index);	
	
	gsl_matrix_set_zero (Uab);
	CalcUab (UtW, Uty, Uab);
//	if (e_mode!=0) {
//		gsl_vector_set_zero (ab);
//		Calcab (W, y, ab);
//	}	

	// Read in header
	uint32_t bgen_header_offset;
	uint32_t bgen_header_length;
	uint32_t bgen_nsamples;
	uint32_t bgen_nsnps;
	uint32_t bgen_flags;
	infile.read(reinterpret_cast<char*>(&bgen_header_offset),4);
	infile.read(reinterpret_cast<char*>(&bgen_header_length),4);
	infile.read(reinterpret_cast<char*>(&bgen_nsnps),4);	
	infile.read(reinterpret_cast<char*>(&bgen_nsamples),4);
	infile.ignore(4+bgen_header_length-20);
	infile.read(reinterpret_cast<char*>(&bgen_flags),4);
	bool CompressedSNPBlocks=bgen_flags&0x1;
	bool LongIds=bgen_flags&0x4;
	double bgen_geno_prob_AA, bgen_geno_prob_AB, bgen_geno_prob_BB, bgen_geno_prob_miss;

	uint32_t bgen_N;
	uint16_t bgen_LS;
	uint16_t bgen_LR;
	uint16_t bgen_LC;	
	uint32_t bgen_SNP_pos;
	uint32_t bgen_LA;
	std::string bgen_A_allele;
	uint32_t bgen_LB;
	std::string bgen_B_allele;
	uint32_t bgen_P;
	size_t unzipped_data_size;
	string id;
	string rs;
	string chr;
	std::cout<<"Warning: WJA hard coded SNP missingness threshold of 10%"<<std::endl;


	
	//start reading genotypes and analyze	
	for (size_t t=0; t<indicator_snp.size(); ++t) 
	{
		
//		if (t>1) {break;}
		if (t%d_pace==0 || t==(ns_total-1)) {ProgressBar ("Reading SNPs  ", t, ns_total-1);}
		// read SNP header
		id.clear();
		rs.clear();
		chr.clear();
		bgen_A_allele.clear();
		bgen_B_allele.clear();
		
		infile.read(reinterpret_cast<char*>(&bgen_N),4);
		infile.read(reinterpret_cast<char*>(&bgen_LS),2);

		std::copy_n(std::istreambuf_iterator<char>(infile), static_cast<size_t>(bgen_LS), std::back_inserter(id));
		infile.ignore(1);
		infile.read(reinterpret_cast<char*>(&bgen_LR),2);

		std::copy_n(std::istreambuf_iterator<char>(infile), static_cast<size_t>(bgen_LR), std::back_inserter(rs));
		infile.ignore(1);
		infile.read(reinterpret_cast<char*>(&bgen_LC),2);
		std::copy_n(std::istreambuf_iterator<char>(infile), static_cast<size_t>(bgen_LC), std::back_inserter(chr));
		infile.ignore(1);
		
		infile.read(reinterpret_cast<char*>(&bgen_SNP_pos),4);
		infile.read(reinterpret_cast<char*>(&bgen_LA),4);
		std::copy_n(std::istreambuf_iterator<char>(infile), static_cast<size_t>(bgen_LA), std::back_inserter(bgen_A_allele));
		infile.ignore(1);

		infile.read(reinterpret_cast<char*>(&bgen_LB),4);

		std::copy_n(std::istreambuf_iterator<char>(infile), static_cast<size_t>(bgen_LB), std::back_inserter(bgen_B_allele));
		infile.ignore(1);
		

		uint16_t unzipped_data[3*bgen_N];

		if (indicator_snp[t]==0) {
			if(CompressedSNPBlocks)
				infile.read(reinterpret_cast<char*>(&bgen_P),4);
			else	
				bgen_P=6*bgen_N;
		
			infile.ignore(static_cast<size_t>(bgen_P));

			continue;
		}


		if(CompressedSNPBlocks)
		{

		
			infile.read(reinterpret_cast<char*>(&bgen_P),4);
			uint8_t zipped_data[bgen_P];

			unzipped_data_size=6*bgen_N;
	
			infile.read(reinterpret_cast<char*>(zipped_data),bgen_P);	

			int result=uncompress(reinterpret_cast<Bytef*>(unzipped_data), reinterpret_cast<uLongf*>(&unzipped_data_size), reinterpret_cast<Bytef*>(zipped_data), static_cast<uLong> (bgen_P));
			assert(result == Z_OK);
				
		}
		else
		{

			bgen_P=6*bgen_N;
			infile.read(reinterpret_cast<char*>(unzipped_data),bgen_P);
		}
	
		x_mean=0.0; c_phen=0; n_miss=0;
		gsl_vector_set_zero(x_miss);
		for (size_t i=0; i<ni_total; ++i) {
			if (indicator_idv[i]==0) {continue;}
			
		
				bgen_geno_prob_AA=static_cast<double>(unzipped_data[i*3])/32768.0;
				bgen_geno_prob_AB=static_cast<double>(unzipped_data[i*3+1])/32768.0;
				bgen_geno_prob_BB=static_cast<double>(unzipped_data[i*3+2])/32768.0;
				// WJA
				bgen_geno_prob_miss=1.0-bgen_geno_prob_AA-bgen_geno_prob_AB-bgen_geno_prob_BB;
				if (bgen_geno_prob_miss>0.1) {gsl_vector_set(x_miss, c_phen, 0.0); n_miss++;}
				else {

					bgen_geno_prob_AA=bgen_geno_prob_AA/(1.0-bgen_geno_prob_miss);
					bgen_geno_prob_AB=bgen_geno_prob_AB/(1.0-bgen_geno_prob_miss);
					bgen_geno_prob_BB=bgen_geno_prob_BB/(1.0-bgen_geno_prob_miss);

					geno=2.0*bgen_geno_prob_AA+bgen_geno_prob_AB;		
				
					gsl_vector_set(x, c_phen, geno); 
					gsl_vector_set(x_miss, c_phen, 1.0); 
					x_mean+=geno;
			}
			c_phen++;
		}	

		x_mean/=(double)(ni_test-n_miss);
	
		for (size_t i=0; i<ni_test; ++i) {
			if (gsl_vector_get (x_miss, i)==0) {gsl_vector_set(x, i, x_mean);}
			geno=gsl_vector_get(x, i);
			if (x_mean>1) {
				gsl_vector_set(x, i, 2-geno);
			}
		}
		
		
		//calculate statistics
		time_start=clock();
		gsl_blas_dgemv (CblasTrans, 1.0, U, x, 0.0, Utx);		
		time_UtX+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		CalcUab(UtW, Uty, Utx, Uab);
//		if (e_mode!=0) {
//			Calcab (W, y, x, ab);
//		}
		
		time_start=clock();
		FUNC_PARAM param1={false, ni_test, n_cvt, eval, Uab, ab, 0};
		
		//3 is before 1
		if (a_mode==3 || a_mode==4) {
			CalcRLScore (l_mle_null, param1, beta, se, p_score);
		}
		
		if (a_mode==1 || a_mode==4) {
			CalcLambda ('R', param1, l_min, l_max, n_region, lambda_remle, logl_H1);	
			CalcRLWald (lambda_remle, param1, beta, se, p_wald);
		}
		
		if (a_mode==2 || a_mode==4) {
			CalcLambda ('L', param1, l_min, l_max, n_region, lambda_mle, logl_H1);
			p_lrt=gsl_cdf_chisq_Q (2.0*(logl_H1-logl_mle_H0), 1);	
		}			
		
		if (x_mean>1) {beta*=-1;}
		
		time_opt+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
		//store summary data
		SUMSTAT SNPs={beta, se, lambda_remle, lambda_mle, p_wald, p_lrt, p_score};
		sumStat.push_back(SNPs);
    }	
	cout<<endl;
	
	gsl_vector_free (x);
	gsl_vector_free (x_miss);
	gsl_vector_free (Utx);
	gsl_matrix_free (Uab);
	gsl_vector_free (ab);
	
	infile.close();
	infile.clear();
	
	return;

}





void MatrixCalcLR (const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty, const gsl_vector *K_eval, const double l_min, const double l_max, const size_t n_region, vector<pair<size_t, double> > &pos_loglr) 
{
	double logl_H0, logl_H1, log_lr, lambda0, lambda1;
	
	gsl_vector *w=gsl_vector_alloc (Uty->size);
	gsl_matrix *Utw=gsl_matrix_alloc (Uty->size, 1);	
	gsl_matrix *Uab=gsl_matrix_alloc (Uty->size, 6);
	gsl_vector *ab=gsl_vector_alloc (6);	
	
	gsl_vector_set_zero(ab);
	gsl_vector_set_all (w, 1.0);
	gsl_vector_view Utw_col=gsl_matrix_column (Utw, 0);	
	gsl_blas_dgemv (CblasTrans, 1.0, U, w, 0.0, &Utw_col.vector);		
	
	CalcUab (Utw, Uty, Uab) ;	
	FUNC_PARAM param0={true, Uty->size, 1, K_eval, Uab, ab, 0};	
	
	CalcLambda('L', param0, l_min, l_max, n_region, lambda0, logl_H0);
	
	for (size_t i=0; i<UtX->size2; ++i) {
		gsl_vector_const_view UtX_col=gsl_matrix_const_column (UtX, i);
		CalcUab(Utw, Uty, &UtX_col.vector, Uab);
		FUNC_PARAM param1={false, UtX->size1, 1, K_eval, Uab, ab, 0};
		
		CalcLambda ('L', param1, l_min, l_max, n_region, lambda1, logl_H1);
		log_lr=logl_H1-logl_H0;				
		
		pos_loglr.push_back(make_pair(i,log_lr) );
	}
	
	gsl_vector_free (w);
	gsl_matrix_free (Utw);
	gsl_matrix_free (Uab);
	gsl_vector_free (ab);
	
	return;
}




void CalcLambda (const char func_name, FUNC_PARAM &params, const double l_min, const double l_max, const size_t n_region, double &lambda, double &logf)
{
	if (func_name!='R' && func_name!='L' && func_name!='r' && func_name!='l') {cout<<"func_name only takes 'R' or 'L': 'R' for log-restricted likelihood, 'L' for log-likelihood."<<endl; return;}
	
	vector<pair<double, double> > lambda_lh;
	
	//evaluate first order derivates in different intervals
	double lambda_l, lambda_h, lambda_interval=log(l_max/l_min)/(double)n_region;
	double dev1_l, dev1_h, logf_l, logf_h;
	
	for (size_t i=0; i<n_region; ++i) {
		lambda_l=l_min*exp(lambda_interval*i);
		lambda_h=l_min*exp(lambda_interval*(i+1.0));
		
		if (func_name=='R' || func_name=='r') {
			dev1_l=LogRL_dev1 (lambda_l, &params);
			dev1_h=LogRL_dev1 (lambda_h, &params);
		}
		else {
			dev1_l=LogL_dev1 (lambda_l, &params);
			dev1_h=LogL_dev1 (lambda_h, &params);
		}
		
		if (dev1_l*dev1_h<=0) {
			lambda_lh.push_back(make_pair(lambda_l, lambda_h));
		}
	}
	
	//if derivates do not change signs in any interval
	if (lambda_lh.empty()) {
		if (func_name=='R' || func_name=='r') {
			logf_l=LogRL_f (l_min, &params);
			logf_h=LogRL_f (l_max, &params);
		}
		else {
			logf_l=LogL_f (l_min, &params);
			logf_h=LogL_f (l_max, &params);
		}
		
		if (logf_l>=logf_h) {lambda=l_min; logf=logf_l;} else {lambda=l_max; logf=logf_h;}
	}
	else {
		//if derivates change signs
		int status;
		int iter=0, max_iter=100;
		double l, l_temp;	
		
		gsl_function F;
		gsl_function_fdf FDF;
		
		F.params=&params;
		FDF.params=&params;
		
		if (func_name=='R' || func_name=='r') {
			F.function=&LogRL_dev1;
			FDF.f=&LogRL_dev1;
			FDF.df=&LogRL_dev2;
			FDF.fdf=&LogRL_dev12;
		}
		else {
			F.function=&LogL_dev1;
			FDF.f=&LogL_dev1;
			FDF.df=&LogL_dev2;
			FDF.fdf=&LogL_dev12;
		}
		
		const gsl_root_fsolver_type *T_f;
		gsl_root_fsolver *s_f;
		T_f=gsl_root_fsolver_brent;
		s_f=gsl_root_fsolver_alloc (T_f);
		
		const gsl_root_fdfsolver_type *T_fdf;
		gsl_root_fdfsolver *s_fdf;
		T_fdf=gsl_root_fdfsolver_newton;
		s_fdf=gsl_root_fdfsolver_alloc(T_fdf);	
		
		for (vector<double>::size_type i=0; i<lambda_lh.size(); ++i) {
			lambda_l=lambda_lh[i].first; lambda_h=lambda_lh[i].second;
			
			gsl_root_fsolver_set (s_f, &F, lambda_l, lambda_h);
			
			do {
				iter++;
				status=gsl_root_fsolver_iterate (s_f);
				l=gsl_root_fsolver_root (s_f);
				lambda_l=gsl_root_fsolver_x_lower (s_f);
				lambda_h=gsl_root_fsolver_x_upper (s_f);
				status=gsl_root_test_interval (lambda_l, lambda_h, 0, 1e-1);		
			}
			while (status==GSL_CONTINUE && iter<max_iter); 				
			
			iter=0;
			
			gsl_root_fdfsolver_set (s_fdf, &FDF, l);	
			
			do {
				iter++;
				status=gsl_root_fdfsolver_iterate (s_fdf);
				l_temp=l;
				l=gsl_root_fdfsolver_root (s_fdf);
				status=gsl_root_test_delta (l, l_temp, 0, 1e-5);		
			}
			while (status==GSL_CONTINUE && iter<max_iter && l>l_min && l<l_max); 
			
			l=l_temp;
			if (l<l_min) {l=l_min;}
			if (l>l_max) {l=l_max;}
			if (func_name=='R' || func_name=='r') {logf_l=LogRL_f (l, &params);} else {logf_l=LogL_f (l, &params);}			
			
			if (i==0) {logf=logf_l; lambda=l;}
			else if (logf<logf_l) {logf=logf_l; lambda=l;}
			else {}
		}
		gsl_root_fsolver_free (s_f);	
		gsl_root_fdfsolver_free (s_fdf);		
		
		if (func_name=='R' || func_name=='r') {
			logf_l=LogRL_f (l_min, &params);
			logf_h=LogRL_f (l_max, &params);
		}
		else {
			logf_l=LogL_f (l_min, &params);
			logf_h=LogL_f (l_max, &params);
		}
		
		if (logf_l>logf) {lambda=l_min; logf=logf_l;} 
		if (logf_h>logf) {lambda=l_max; logf=logf_h;}
	}
	
	return;
}





//calculate lambda in the null model
void CalcLambda (const char func_name, const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Uty, const double l_min, const double l_max, const size_t n_region, double &lambda, double &logl_H0)
{
	if (func_name!='R' && func_name!='L' && func_name!='r' && func_name!='l') {cout<<"func_name only takes 'R' or 'L': 'R' for log-restricted likelihood, 'L' for log-likelihood."<<endl; return;}

	size_t n_cvt=UtW->size2, ni_test=UtW->size1;
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	gsl_matrix *Uab=gsl_matrix_alloc (ni_test, n_index);	
	gsl_vector *ab=gsl_vector_alloc (n_index);	
	
	gsl_matrix_set_zero (Uab);
	CalcUab (UtW, Uty, Uab);
//	if (e_mode!=0) {
//		gsl_vector_set_zero (ab);
//		Calcab (W, y, ab);
//	}
		
	FUNC_PARAM param0={true, ni_test, n_cvt, eval, Uab, ab, 0};
	
	CalcLambda(func_name, param0, l_min, l_max, n_region, lambda, logl_H0);
	
	gsl_matrix_free(Uab);	
	gsl_vector_free(ab);	
	
	return;
}
	
	
//obtain REMLE estimate for PVE using lambda_remle
void CalcPve (const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Uty, const double lambda, const double trace_G, double &pve, double &pve_se)
{
	size_t n_cvt=UtW->size2, ni_test=UtW->size1;
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	gsl_matrix *Uab=gsl_matrix_alloc (ni_test, n_index);	
	gsl_vector *ab=gsl_vector_alloc (n_index);	
	
	gsl_matrix_set_zero (Uab);
	CalcUab (UtW, Uty, Uab);
	//	if (e_mode!=0) {
	//		gsl_vector_set_zero (ab);
	//		Calcab (W, y, ab);
	//	}
	
	FUNC_PARAM param0={true, ni_test, n_cvt, eval, Uab, ab, 0};
	
	double se=sqrt(-1.0/LogRL_dev2 (lambda, &param0));
	
	pve=trace_G*lambda/(trace_G*lambda+1.0);
	pve_se=trace_G/((trace_G*lambda+1.0)*(trace_G*lambda+1.0))*se;
	
	gsl_matrix_free (Uab);
	gsl_vector_free (ab);	
	return;
}

//obtain REML estimate for Vg and Ve using lambda_remle
//obtain beta and se(beta) for coefficients
//ab is not used when e_mode==0
void CalcLmmVgVeBeta (const gsl_vector *eval, const gsl_matrix *UtW, const gsl_vector *Uty, const double lambda, double &vg, double &ve, gsl_vector *beta, gsl_vector *se_beta)
{
	size_t n_cvt=UtW->size2, ni_test=UtW->size1;
	size_t n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	
	gsl_matrix *Uab=gsl_matrix_alloc (ni_test, n_index);	
	gsl_vector *ab=gsl_vector_alloc (n_index);	
	gsl_matrix *Pab=gsl_matrix_alloc (n_cvt+2, n_index);
	gsl_vector *Hi_eval=gsl_vector_alloc(eval->size);
	gsl_vector *v_temp=gsl_vector_alloc(eval->size);
	gsl_matrix *HiW=gsl_matrix_alloc(eval->size, UtW->size2);
	gsl_matrix *WHiW=gsl_matrix_alloc(UtW->size2, UtW->size2);
	gsl_vector *WHiy=gsl_vector_alloc(UtW->size2);
	gsl_matrix *Vbeta=gsl_matrix_alloc(UtW->size2, UtW->size2);
	
	gsl_matrix_set_zero (Uab);
	CalcUab (UtW, Uty, Uab);	
	
	gsl_vector_memcpy (v_temp, eval);
	gsl_vector_scale (v_temp, lambda);
	gsl_vector_set_all (Hi_eval, 1.0);
	gsl_vector_add_constant (v_temp, 1.0);
	gsl_vector_div (Hi_eval, v_temp);
	
	//calculate beta
	gsl_matrix_memcpy (HiW, UtW);
	for (size_t i=0; i<UtW->size2; i++) {
		gsl_vector_view HiW_col=gsl_matrix_column(HiW, i);
		gsl_vector_mul(&HiW_col.vector, Hi_eval);
	}
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, HiW, UtW, 0.0, WHiW);
	gsl_blas_dgemv (CblasTrans, 1.0, HiW, Uty, 0.0, WHiy);
	
	int sig;
	gsl_permutation * pmt=gsl_permutation_alloc (UtW->size2);
	LUDecomp (WHiW, pmt, &sig);
	LUSolve (WHiW, pmt, WHiy, beta);
	LUInvert (WHiW, pmt, Vbeta);
		
	//calculate vg and ve
	CalcPab (n_cvt, 0, Hi_eval, Uab, ab, Pab);	
	
	size_t index_yy=GetabIndex (n_cvt+2, n_cvt+2, n_cvt);	
	double P_yy=gsl_matrix_get (Pab, n_cvt, index_yy);	
	
	ve=P_yy/(double)(ni_test-n_cvt);
	vg=ve*lambda;
	
	//with ve, calculate se(beta)
	gsl_matrix_scale(Vbeta, ve);
	
	//obtain se_beta
	for (size_t i=0; i<Vbeta->size1; i++) {
		gsl_vector_set (se_beta, i, sqrt(gsl_matrix_get(Vbeta, i, i) ) );
	}
	
	gsl_matrix_free(Uab);
	gsl_matrix_free(Pab);
	gsl_vector_free(ab);
	gsl_vector_free(Hi_eval);
	gsl_vector_free(v_temp);
	gsl_matrix_free(HiW);
	gsl_matrix_free(WHiW);
	gsl_vector_free(WHiy);
	gsl_matrix_free(Vbeta);
	
	gsl_permutation_free(pmt);
	return;
}

