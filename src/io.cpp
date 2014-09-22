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
#include <string>
#include <iomanip>
#include <bitset>
#include <vector>
#include <map>
#include <set>
#include <cstring>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"

#include "lapack.h"
#include "gzstream.h"
#include "mathfunc.h"

#ifdef FORCE_FLOAT
#include "io_float.h"
#else
#include "io.h"
#endif


using namespace std;



//Print process bar
void ProgressBar (string str, double p, double total)
{
	double progress = (100.0 * p / total); 
	int barsize = (int) (progress / 2.0); 
	char bar[51];
	
	cout<<str;
	for (int i = 0; i <50; i++) {
		if (i<barsize) {bar[i] = '=';}
		else {bar[i]=' ';}
		cout<<bar[i];
	}
	cout<<setprecision(2)<<fixed<<progress<<"%\r"<<flush;
	
	return;
}


//Print process bar (with acceptance ratio)
void ProgressBar (string str, double p, double total, double ratio)
{
	double progress = (100.0 * p / total); 
	int barsize = (int) (progress / 2.0); 
	char bar[51];
	
	cout<<str;
	for (int i = 0; i <50; i++) {
		if (i<barsize) {bar[i] = '=';}
		else {bar[i]=' ';}
		cout<<bar[i];
	}
	cout<<setprecision(2)<<fixed<<progress<<"%    "<<ratio<<"\r"<<flush;
	
	
	return;
}

// in case files are ended with "\r" or "\r\n"
std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

//Read snp file
bool ReadFile_snps (const string &file_snps, set<string> &setSnps)
{
	setSnps.clear();

	ifstream infile (file_snps.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open snps file: "<<file_snps<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	while (getline(infile, line)) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		setSnps.insert(ch_ptr); 
	}
	
	infile.close();
	infile.clear();	
	
	return true;
}


//Read log file
bool ReadFile_log (const string &file_log, double &pheno_mean)
{
	ifstream infile (file_log.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open log file: "<<file_log<<endl; return false;}
	
	string line;
	char *ch_ptr;
	size_t flag=0;
	
	while (getline(infile, line)) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		ch_ptr=strtok (NULL, " , \t");
		
		if (ch_ptr!=NULL && strcmp(ch_ptr, "estimated")==0) {
			ch_ptr=strtok (NULL, " , \t");
			if (ch_ptr!=NULL && strcmp(ch_ptr, "mean")==0) {
				ch_ptr=strtok (NULL, " , \t");
				if (ch_ptr!=NULL && strcmp(ch_ptr, "=")==0) {
					ch_ptr=strtok (NULL, " , \t");
					pheno_mean=atof(ch_ptr);
					flag=1;
				}
			}
		}
		
		if (flag==1) {break;}
	}
	
	infile.close();
	infile.clear();	
	
	return true;
}


//Read bimbam annotation file
bool ReadFile_anno (const string &file_anno, map<string, string> &mapRS2chr, map<string, long int> &mapRS2bp, map<string, double> &mapRS2cM)
{
	mapRS2chr.clear();
	mapRS2bp.clear();
	
	ifstream infile (file_anno.c_str(), ifstream::in);
	if (!infile) {cout<<"error opening annotation file: "<<file_anno<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	string rs;
	long int b_pos;
	string chr;
	double cM;
	
	while (!safeGetline(infile, line).eof()) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		rs=ch_ptr;
		ch_ptr=strtok (NULL, " , \t");
		if (strcmp(ch_ptr, "NA")==0) {b_pos=-9;} else {b_pos=atol(ch_ptr);}
		ch_ptr=strtok (NULL, " , \t");
		if (ch_ptr==NULL || strcmp(ch_ptr, "NA")==0) {chr="-9";} else {chr=ch_ptr;}
		ch_ptr=strtok (NULL, " , \t");
		if (ch_ptr==NULL || strcmp(ch_ptr, "NA")==0) {cM=-9;} else {cM=atof(ch_ptr);}
		
		mapRS2chr[rs]=chr;
		mapRS2bp[rs]=b_pos;
		mapRS2cM[rs]=cM;
	}
	
	infile.close();
	infile.clear();	
	
	return true;
}

//read one column of phenotype
bool ReadFile_column (const string &file_pheno, vector<int> &indicator_idv, vector<double> &pheno, const int &p_column)
{
	indicator_idv.clear();
	pheno.clear();
	
	igzstream infile (file_pheno.c_str(), igzstream::in);
//	ifstream infile (file_pheno.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open phenotype file: "<<file_pheno<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	string id;
	double p;
	while (!safeGetline(infile, line).eof()) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		for (int i=0; i<(p_column-1); ++i) {
			ch_ptr=strtok (NULL, " , \t");	
		}		
		if (strcmp(ch_ptr, "NA")==0) {indicator_idv.push_back(0); pheno.push_back(-9);}		//pheno is different from pimass2
		else {p=atof(ch_ptr); indicator_idv.push_back(1); pheno.push_back(p);}
	}
	
	infile.close();
	infile.clear();	
	
	return true;
}



//Read bimbam phenotype file, p_column=1, 2 ...
bool ReadFile_pheno (const string &file_pheno, vector<vector<int> > &indicator_pheno, vector<vector<double> > &pheno, const vector<size_t> &p_column)
{
	indicator_pheno.clear();
	pheno.clear();
	
	igzstream infile (file_pheno.c_str(), igzstream::in);
//	ifstream infile (file_pheno.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open phenotype file: "<<file_pheno<<endl; return false;}

	string line;
	char *ch_ptr;
  
	string id;
	double p;
	
	vector<double> pheno_row;
	vector<int> ind_pheno_row;
	
	size_t p_max=*max_element(p_column.begin(), p_column.end() );
	map<size_t, size_t> mapP2c;
	for (size_t i=0; i<p_column.size(); i++) {
		mapP2c[p_column[i]]=i;
		pheno_row.push_back(-9);
		ind_pheno_row.push_back(0);
	}	
	
	while (!safeGetline(infile, line).eof()) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		
		size_t i=0;
		while (i<p_max ) {			
			if (mapP2c.count(i+1)!=0) {
				if (strcmp(ch_ptr, "NA")==0) {ind_pheno_row[mapP2c[i+1]]=0; pheno_row[mapP2c[i+1]]=-9;}
				else {p=atof(ch_ptr); ind_pheno_row[mapP2c[i+1]]=1; pheno_row[mapP2c[i+1]]=p;}
			}
			i++;
			ch_ptr=strtok (NULL, " , \t");	
		}
		
		indicator_pheno.push_back(ind_pheno_row);	
		pheno.push_back(pheno_row);			
	}
 
	infile.close();
	infile.clear();	
	
	return true;
}


bool ReadFile_cvt (const string &file_cvt, vector<int> &indicator_cvt, vector<vector<double> > &cvt, size_t &n_cvt)
{
	indicator_cvt.clear();
	
	ifstream infile (file_cvt.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open covariates file: "<<file_cvt<<endl; return false;}
	
	string line;
	char *ch_ptr;
	double d;	
	
	int flag_na=0;	
	
	while (!safeGetline(infile, line).eof()) {
		vector<double> v_d; flag_na=0;
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		while (ch_ptr!=NULL) {
			if (strcmp(ch_ptr, "NA")==0) {flag_na=1; d=-9;}
			else {d=atof(ch_ptr);}
			
			v_d.push_back(d);
			ch_ptr=strtok (NULL, " , \t");	
		}
		if (flag_na==0) {indicator_cvt.push_back(1);} else {indicator_cvt.push_back(0);} 
		cvt.push_back(v_d);
	}
	
	if (indicator_cvt.empty()) {n_cvt=0;}
	else {
		flag_na=0;
		for (vector<int>::size_type i=0; i<indicator_cvt.size(); ++i) {
			if (indicator_cvt[i]==0) {continue;}
			
			if (flag_na==0) {flag_na=1; n_cvt=cvt[i].size();}
			if (flag_na!=0 && n_cvt!=cvt[i].size()) {cout<<"error! number of covariates in row "<<i<<" do not match other rows."<<endl; return false;}
		}
	}
	
	infile.close();
	infile.clear();	
	
	return true;
}



//Read .bim file
bool ReadFile_bim (const string &file_bim, vector<SNPINFO> &snpInfo)
{
	snpInfo.clear();
	
	ifstream infile (file_bim.c_str(), ifstream::in);
	if (!infile) {cout<<"error opening .bim file: "<<file_bim<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	string rs;
	long int b_pos;
	string chr;
	double cM;
	string major;
	string minor;
	
	while (getline(infile, line)) {
		ch_ptr=strtok ((char *)line.c_str(), " \t");
		chr=ch_ptr;
		ch_ptr=strtok (NULL, " \t");
		rs=ch_ptr;
		ch_ptr=strtok (NULL, " \t");
		cM=atof(ch_ptr);
		ch_ptr=strtok (NULL, " \t");
		b_pos=atol(ch_ptr);
		ch_ptr=strtok (NULL, " \t");
		minor=ch_ptr;
		ch_ptr=strtok (NULL, " \t");
		major=ch_ptr;
		
		SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, -9, -9, -9};
		snpInfo.push_back(sInfo);
	}
	
	infile.close();
	infile.clear();	
	return true;
}


//Read .fam file
bool ReadFile_fam (const string &file_fam, vector<vector<int> > &indicator_pheno, vector<vector<double> > &pheno, map<string, int> &mapID2num, const vector<size_t> &p_column)
{
	indicator_pheno.clear();
	pheno.clear();
	mapID2num.clear();	
	
	igzstream infile (file_fam.c_str(), igzstream::in);
	//ifstream infile (file_fam.c_str(), ifstream::in);
	if (!infile) {cout<<"error opening .fam file: "<<file_fam<<endl; return false;}

	string line;
	char *ch_ptr;

	string id;
	int c=0;
	double p;

	vector<double> pheno_row;
	vector<int> ind_pheno_row;
	
	size_t p_max=*max_element(p_column.begin(), p_column.end() );
	map<size_t, size_t> mapP2c;
	for (size_t i=0; i<p_column.size(); i++) {
		mapP2c[p_column[i]]=i;
		pheno_row.push_back(-9);
		ind_pheno_row.push_back(0);
	}	
	
	while (!safeGetline(infile, line).eof()) {
		ch_ptr=strtok ((char *)line.c_str(), " \t");
		ch_ptr=strtok (NULL, " \t");
		id=ch_ptr;
		ch_ptr=strtok (NULL, " \t");
		ch_ptr=strtok (NULL, " \t");
		ch_ptr=strtok (NULL, " \t");
		ch_ptr=strtok (NULL, " \t");
		
		size_t i=0;
		while (i<p_max ) {
			if (mapP2c.count(i+1)!=0 ) {
				if (strcmp(ch_ptr, "NA")==0) {
					ind_pheno_row[mapP2c[i+1]]=0; pheno_row[mapP2c[i+1]]=-9;
				} else {
					p=atof(ch_ptr);
					
					if (p==-9) {ind_pheno_row[mapP2c[i+1]]=0; pheno_row[mapP2c[i+1]]=-9;}
					else {ind_pheno_row[mapP2c[i+1]]=1; pheno_row[mapP2c[i+1]]=p;}
				}
			}
			i++;
			ch_ptr=strtok (NULL, " , \t");	
		}
		
		indicator_pheno.push_back(ind_pheno_row);
		pheno.push_back(pheno_row);				
		
		mapID2num[id]=c; c++;
	}
 
	infile.close();
	infile.clear();	
	return true;
}






//Read bimbam mean genotype file, the first time, to obtain #SNPs for analysis (ns_test) and total #SNP (ns_total)
bool ReadFile_geno (const string &file_geno, const set<string> &setSnps, const gsl_matrix *W, vector<int> &indicator_idv, vector<int> &indicator_snp, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, map<string, string> &mapRS2chr, map<string, long int> &mapRS2bp, map<string, double> &mapRS2cM, vector<SNPINFO> &snpInfo, size_t &ns_test)
{
	indicator_snp.clear();
	snpInfo.clear();
	
	igzstream infile (file_geno.c_str(), igzstream::in);
//	ifstream infile (file_geno.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return false;}

	gsl_vector *genotype=gsl_vector_alloc (W->size1);
	gsl_vector *genotype_miss=gsl_vector_alloc (W->size1);
	gsl_matrix *WtW=gsl_matrix_alloc (W->size2, W->size2);
	gsl_matrix *WtWi=gsl_matrix_alloc (W->size2, W->size2);
	gsl_vector *Wtx=gsl_vector_alloc (W->size2);
	gsl_vector *WtWiWtx=gsl_vector_alloc (W->size2);
	gsl_permutation * pmt=gsl_permutation_alloc (W->size2);
	
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
	int sig;	
	LUDecomp (WtW, pmt, &sig);
	LUInvert (WtW, pmt, WtWi);
	
	double v_x, v_w;
	int c_idv=0;
	
	string line;
	char *ch_ptr;
		
	string rs;
	long int b_pos;
	string chr;
	string major;
	string minor;
	double cM;
  
	double maf, geno, geno_old;
	size_t n_miss;
	size_t n_0, n_1, n_2;
	int flag_poly;
	
	int ni_total=indicator_idv.size();
	int ni_test=0;
	for (int i=0; i<ni_total; ++i) {
		ni_test+=indicator_idv[i];
	}
	ns_test=0;
	
	while (!safeGetline(infile, line).eof()) {		
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		rs=ch_ptr;
		ch_ptr=strtok (NULL, " , \t");
		minor=ch_ptr;
		ch_ptr=strtok (NULL, " , \t");
		major=ch_ptr;
		
		if (setSnps.size()!=0 && setSnps.count(rs)==0) {
			SNPINFO sInfo={"-9", rs, -9, -9, minor, major, -9, -9, -9};
			snpInfo.push_back(sInfo);
			indicator_snp.push_back(0);
			continue;
		}
				
		if (mapRS2bp.count(rs)==0) {chr="-9"; b_pos=-9;cM=-9;}
		else {b_pos=mapRS2bp[rs]; chr=mapRS2chr[rs]; cM=mapRS2cM[rs];}		
				
		maf=0; n_miss=0; flag_poly=0; geno_old=-9;
		n_0=0; n_1=0; n_2=0;
		c_idv=0; gsl_vector_set_zero (genotype_miss);
		for (int i=0; i<ni_total; ++i) {
			ch_ptr=strtok (NULL, " , \t");
			if (indicator_idv[i]==0) {continue;}		

			if (strcmp(ch_ptr, "NA")==0) {gsl_vector_set (genotype_miss, c_idv, 1); n_miss++; c_idv++; continue;}
			
			geno=atof(ch_ptr);
			if (geno>=0 && geno<=0.5) {n_0++;}
			if (geno>0.5 && geno<1.5) {n_1++;}
			if (geno>=1.5 && geno<=2.0) {n_2++;}
			
			gsl_vector_set (genotype, c_idv, geno); 
			
//			if (geno<0) {n_miss++; continue;}
			
			if (flag_poly==0) {geno_old=geno; flag_poly=2;}
			if (flag_poly==2 && geno!=geno_old) {flag_poly=1;}
			
			maf+=geno;
			
			c_idv++;
		}
		maf/=2.0*(double)(ni_test-n_miss);	
		
		SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, n_miss, (double)n_miss/(double)ni_test, maf};
		snpInfo.push_back(sInfo);
		
		if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}
		
		if ( (maf<maf_level || maf> (1.0-maf_level)) && maf_level!=-1 ) {indicator_snp.push_back(0); continue;}
		
		if (flag_poly!=1) {indicator_snp.push_back(0); continue;}
		
		if (hwe_level!=0) {
			if (CalcHWE(n_0, n_2, n_1)<hwe_level) {indicator_snp.push_back(0); continue;}
		}
		
		//filter SNP if it is correlated with W
		for (size_t i=0; i<genotype->size; ++i) {			
			if (gsl_vector_get (genotype_miss, i)==1) {geno=maf*2.0; gsl_vector_set (genotype, i, geno);}		
		}
		
		gsl_blas_dgemv (CblasTrans, 1.0, W, genotype, 0.0, Wtx);
		gsl_blas_dgemv (CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
		gsl_blas_ddot (genotype, genotype, &v_x);
		gsl_blas_ddot (Wtx, WtWiWtx, &v_w);
		
		if (v_w/v_x >= r2_level) {indicator_snp.push_back(0); continue;}
		
		indicator_snp.push_back(1); 
		ns_test++;
	}
	
	gsl_vector_free (genotype);
	gsl_vector_free (genotype_miss);
	gsl_matrix_free (WtW);
	gsl_matrix_free (WtWi);
	gsl_vector_free (Wtx);
	gsl_vector_free (WtWiWtx);
	gsl_permutation_free (pmt);
	
	infile.close();
	infile.clear();	
	
	return true;
}





      
//Read bed file, the first time
bool ReadFile_bed (const string &file_bed, const set<string> &setSnps, const gsl_matrix *W, vector<int> &indicator_idv, vector<int> &indicator_snp, vector<SNPINFO> &snpInfo, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, size_t &ns_test)
{
	indicator_snp.clear();
	size_t ns_total=snpInfo.size();
	
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}

	gsl_vector *genotype=gsl_vector_alloc (W->size1);
	gsl_vector *genotype_miss=gsl_vector_alloc (W->size1);
	gsl_matrix *WtW=gsl_matrix_alloc (W->size2, W->size2);
	gsl_matrix *WtWi=gsl_matrix_alloc (W->size2, W->size2);
	gsl_vector *Wtx=gsl_vector_alloc (W->size2);
	gsl_vector *WtWiWtx=gsl_vector_alloc (W->size2);
	gsl_permutation * pmt=gsl_permutation_alloc (W->size2);
	
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, W, W, 0.0, WtW);
	int sig;	
	LUDecomp (WtW, pmt, &sig);
	LUInvert (WtW, pmt, WtWi);
	
	double v_x, v_w, geno;
	size_t c_idv=0;
	
	char ch[1];
	bitset<8> b;
  	
	size_t ni_total=indicator_idv.size();
	size_t ni_test=0;
	for (size_t i=0; i<ni_total; ++i) {
		ni_test+=indicator_idv[i];
	}
	ns_test=0;
	
	//calculate n_bit and c, the number of bit for each snp
	size_t n_bit;
	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1;}

	//ignore the first three majic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}
	
	double maf;
	size_t n_miss;
	size_t n_0, n_1, n_2, c;	
	
	//start reading snps and doing association test
	for (size_t t=0; t<ns_total; ++t) {
		infile.seekg(t*n_bit+3);		//n_bit, and 3 is the number of magic numbers
		
		if (setSnps.size()!=0 && setSnps.count(snpInfo[t].rs_number)==0) {
			snpInfo[t].n_miss=-9;
			snpInfo[t].missingness=-9;
			snpInfo[t].maf=-9;
			indicator_snp.push_back(0);
			continue;
		}

		//read genotypes
		c=0; maf=0.0; n_miss=0; n_0=0; n_1=0; n_2=0;
		c_idv=0; gsl_vector_set_zero (genotype_miss);
		for (size_t i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && c==ni_total) {break;}
				if (indicator_idv[c]==0) {c++; continue;}
				c++;
				
				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(genotype, c_idv, 2.0); maf+=2.0; n_2++;}
					else {gsl_vector_set(genotype, c_idv, 1.0); maf+=1.0; n_1++;}
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(genotype, c_idv, 0.0); maf+=0.0; n_0++;}                                  
					else {gsl_vector_set(genotype_miss, c_idv, 1); n_miss++; }
				}
				c_idv++;
			}
		}
		maf/=2.0*(double)(ni_test-n_miss);
		
		snpInfo[t].n_miss=n_miss;
		snpInfo[t].missingness=(double)n_miss/(double)ni_test;
		snpInfo[t].maf=maf;
		
		if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}
		
		if ( (maf<maf_level || maf> (1.0-maf_level)) && maf_level!=-1 ) {indicator_snp.push_back(0); continue;}
		
		if ( (n_0+n_1)==0 || (n_1+n_2)==0 || (n_2+n_0)==0) {indicator_snp.push_back(0); continue;}
		
		if (hwe_level!=1) {
			if (CalcHWE(n_0, n_2, n_1)<hwe_level) {indicator_snp.push_back(0); continue;}
		}
			
		
		//filter SNP if it is correlated with W
		for (size_t i=0; i<genotype->size; ++i) {			
			if (gsl_vector_get (genotype_miss, i)==1) {geno=maf*2.0; gsl_vector_set (genotype, i, geno);}		
		}
		
		gsl_blas_dgemv (CblasTrans, 1.0, W, genotype, 0.0, Wtx);
		gsl_blas_dgemv (CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
		gsl_blas_ddot (genotype, genotype, &v_x);
		gsl_blas_ddot (Wtx, WtWiWtx, &v_w);
		
		if (v_w/v_x > r2_level) {indicator_snp.push_back(0); continue;}
		
		indicator_snp.push_back(1); 
		ns_test++;
	}
	
	gsl_vector_free (genotype);
	gsl_vector_free (genotype_miss);
	gsl_matrix_free (WtW);
	gsl_matrix_free (WtWi);
	gsl_vector_free (Wtx);
	gsl_vector_free (WtWiWtx);
	gsl_permutation_free (pmt);
		  
	infile.close();
	infile.clear();	
	
	return true;
}



void ReadFile_kin (const string &file_kin, vector<int> &indicator_idv, map<string, int> &mapID2num, const size_t k_mode, bool &error, gsl_matrix *G) 
{
	igzstream infile (file_kin.c_str(), igzstream::in);
//	ifstream infile (file_kin.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open kinship file: "<<file_kin<<endl; error=true; return;}
	
	size_t ni_total=indicator_idv.size();
	
	gsl_matrix_set_zero (G);
	
	string line;
	char *ch_ptr;	
	double d;
	
	if (k_mode==1) {
		size_t i_test=0, i_total=0, j_test=0, j_total=0;
		while (getline(infile, line)) {
			if (i_total==ni_total) {cout<<"error! number of rows in the kinship file is larger than the number of phentypes."<<endl; error=true;}			
			
			if (indicator_idv[i_total]==0) {i_total++; continue;}
			
			j_total=0; j_test=0;
			ch_ptr=strtok ((char *)line.c_str(), " , \t");
			while (ch_ptr!=NULL) {
				if (j_total==ni_total) {cout<<"error! number of columns in the kinship file is larger than the number of phentypes for row = "<<i_total<<endl; error=true;}
				
				d=atof(ch_ptr);
				if (indicator_idv[j_total]==1) {gsl_matrix_set (G, i_test, j_test, d); j_test++;}				
				j_total++;
				
				ch_ptr=strtok (NULL, " , \t");
			}
			if (j_total!=ni_total) {cout<<"error! number of columns in the kinship file do not match the number of phentypes for row = "<<i_total<<endl; error=true;}
			i_total++; i_test++;			
		}
		if (i_total!=ni_total) {cout<<"error! number of rows in the kinship file do not match the number of phentypes."<<endl; error=true;}
	}	
	else {  
		map<size_t, size_t> mapID2ID;
		size_t c=0;
		for (size_t i=0; i<indicator_idv.size(); i++) {
			if (indicator_idv[i]==1) {mapID2ID[i]=c; c++;}
		}
		
		string id1, id2;
		double Cov_d;
		size_t n_id1, n_id2;
		
		while (getline(infile, line)) {
			ch_ptr=strtok ((char *)line.c_str(), " , \t");
			id1=ch_ptr;
			ch_ptr=strtok (NULL, " , \t");
			id2=ch_ptr;
			ch_ptr=strtok (NULL, " , \t");
			d=atof(ch_ptr);
			if (mapID2num.count(id1)==0 || mapID2num.count(id2)==0) {continue;}
			if (indicator_idv[mapID2num[id1]]==0 || indicator_idv[mapID2num[id2]]==0) {continue;}
			
			n_id1=mapID2ID[mapID2num[id1]];
			n_id2=mapID2ID[mapID2num[id2]];
			
			Cov_d=gsl_matrix_get(G, n_id1, n_id2);
			if (Cov_d!=0 && Cov_d!=d) {cout<<"error! redundant and unequal terms in the kinship file, for id1 = "<<id1<<" and id2 = "<<id2<<endl;}
			else {
				gsl_matrix_set(G, n_id1, n_id2, d);
				gsl_matrix_set(G, n_id2, n_id1, d);
			}
		}
	}
	
	infile.close();
	infile.clear();	
	
	return;
}


void ReadFile_mk (const string &file_mk, vector<int> &indicator_idv, map<string, int> &mapID2num, const size_t k_mode, bool &error, gsl_matrix *G) 
{
	igzstream infile (file_mk.c_str(), igzstream::in);
	if (!infile) {cout<<"error! fail to open file: "<<file_mk<<endl; error=true; return;}

	string file_kin, line;

	size_t i=0;
	while (getline(infile, line)) {
	  file_kin=line.c_str();
	  gsl_matrix_view G_sub=gsl_matrix_submatrix(G, 0, i*G->size1, G->size1, G->size1);
	  ReadFile_kin (file_kin, indicator_idv, mapID2num, k_mode, error, &G_sub.matrix);
	  i++;
	}

	infile.close();
	infile.clear();	
	return;
}


void ReadFile_eigenU (const string &file_ku, bool &error, gsl_matrix *U) 
{
	igzstream infile (file_ku.c_str(), igzstream::in);
//	ifstream infile (file_ku.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open the U file: "<<file_ku<<endl; error=true; return;}
	
	size_t n_row=U->size1, n_col=U->size2, i_row=0, i_col=0;
	
	gsl_matrix_set_zero (U);
	
	string line;
	char *ch_ptr;	
	double d;
	
	while (getline(infile, line)) {
		if (i_row==n_row) {cout<<"error! number of rows in the U file is larger than expected."<<endl; error=true;}			
				
		i_col=0;
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		while (ch_ptr!=NULL) {
			if (i_col==n_col) {cout<<"error! number of columns in the U file is larger than expected, for row = "<<i_row<<endl; error=true;}
			
			d=atof(ch_ptr);
			gsl_matrix_set (U, i_row, i_col, d);			
			i_col++;
			
			ch_ptr=strtok (NULL, " , \t");
		}
		
		i_row++;
	}
		
	infile.close();
	infile.clear();	
	
	return;
}




void ReadFile_eigenD (const string &file_kd, bool &error, gsl_vector *eval) 
{
	igzstream infile (file_kd.c_str(), igzstream::in);
//	ifstream infile (file_kd.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open the D file: "<<file_kd<<endl; error=true; return;}
	
	size_t n_row=eval->size, i_row=0;
	
	gsl_vector_set_zero (eval);
	
	string line;
	char *ch_ptr;	
	double d;
	
	while (getline(infile, line)) {
		if (i_row==n_row) {cout<<"error! number of rows in the D file is larger than expected."<<endl; error=true;}			
		
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		d=atof(ch_ptr);
		
		ch_ptr=strtok (NULL, " , \t");
		if (ch_ptr!=NULL) {cout<<"error! number of columns in the D file is larger than expected, for row = "<<i_row<<endl; error=true;}
		
		gsl_vector_set (eval, i_row, d);
		
		i_row++;
	}
	
	infile.close();
	infile.clear();	
	
	return;
}



//read bimbam mean genotype file and calculate kinship matrix
bool BimbamKin (const string &file_geno, vector<int> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin) 
{
	igzstream infile (file_geno.c_str(), igzstream::in);
	//ifstream infile (file_geno.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	size_t n_miss;
	double d, geno_mean, geno_var;
	
	size_t ni_total=matrix_kin->size1;
	gsl_vector *geno=gsl_vector_alloc (ni_total);
	gsl_vector *geno_miss=gsl_vector_alloc (ni_total);

	size_t ns_test=0;
	for (size_t t=0; t<indicator_snp.size(); ++t) {
		!safeGetline(infile, line).eof();
		if (t%display_pace==0 || t==(indicator_snp.size()-1)) {ProgressBar ("Reading SNPs  ", t, indicator_snp.size()-1);}
		if (indicator_snp[t]==0) {continue;}
		
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		ch_ptr=strtok (NULL, " , \t");
		ch_ptr=strtok (NULL, " , \t");
		
		geno_mean=0.0; n_miss=0; geno_var=0.0;
		gsl_vector_set_all(geno_miss, 0);
		for (size_t i=0; i<ni_total; ++i) {
			ch_ptr=strtok (NULL, " , \t");
			if (strcmp(ch_ptr, "NA")==0) {gsl_vector_set(geno_miss, i, 0); n_miss++;}
			else {
				d=atof(ch_ptr);
				gsl_vector_set (geno, i, d);
				gsl_vector_set (geno_miss, i, 1);
				geno_mean+=d;
				geno_var+=d*d;
			}
		}
		
		geno_mean/=(double)(ni_total-n_miss);
		geno_var+=geno_mean*geno_mean*(double)n_miss;
		geno_var/=(double)ni_total;
		geno_var-=geno_mean*geno_mean;
//		geno_var=geno_mean*(1-geno_mean*0.5);
		
		for (size_t i=0; i<ni_total; ++i) {
			if (gsl_vector_get (geno_miss, i)==0) {gsl_vector_set(geno, i, geno_mean);}
		}		
		
		gsl_vector_add_constant (geno, -1.0*geno_mean);
		
		if (geno_var!=0) {
			if (k_mode==1) {gsl_blas_dsyr (CblasUpper, 1.0, geno, matrix_kin);}
			else if (k_mode==2) {gsl_blas_dsyr (CblasUpper, 1.0/geno_var, geno, matrix_kin);}
			else {cout<<"Unknown kinship mode."<<endl;}
		}
		
		ns_test++;
    }	
	cout<<endl;
	
	gsl_matrix_scale (matrix_kin, 1.0/(double)ns_test);
	
	for (size_t i=0; i<ni_total; ++i) {
		for (size_t j=0; j<i; ++j) {
			d=gsl_matrix_get (matrix_kin, j, i);
			gsl_matrix_set (matrix_kin, i, j, d);
		}
	}
	
	gsl_vector_free (geno);
	gsl_vector_free (geno_miss);
	
	infile.close();
	infile.clear();	
	
	return true;
}







bool PlinkKin (const string &file_bed, vector<int> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin) 
{
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}
		
	char ch[1];
	bitset<8> b;
	
	size_t n_miss, ci_total;
	double d, geno_mean, geno_var;
	
	size_t ni_total=matrix_kin->size1;
	gsl_vector *geno=gsl_vector_alloc (ni_total);

	size_t ns_test=0;
	int n_bit;
	
	//calculate n_bit and c, the number of bit for each snp
	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1; }

	//print the first three majic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}	
	
	for (size_t t=0; t<indicator_snp.size(); ++t) {
		if (t%display_pace==0 || t==(indicator_snp.size()-1)) {ProgressBar ("Reading SNPs  ", t, indicator_snp.size()-1);}
		if (indicator_snp[t]==0) {continue;}
		
		infile.seekg(t*n_bit+3);		//n_bit, and 3 is the number of magic numbers
		
		//read genotypes
		geno_mean=0.0;	n_miss=0; ci_total=0; geno_var=0.0;
		for (int i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && ci_total==ni_total) {break;}

				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(geno, ci_total, 2.0); geno_mean+=2.0; geno_var+=4.0; }
					else {gsl_vector_set(geno, ci_total, 1.0); geno_mean+=1.0; geno_var+=1.0;}
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(geno, ci_total, 0.0); }      
					else {gsl_vector_set(geno, ci_total, -9.0); n_miss++; }
				}

				ci_total++;
			}
		}
				
		geno_mean/=(double)(ni_total-n_miss);
		geno_var+=geno_mean*geno_mean*(double)n_miss;
		geno_var/=(double)ni_total;
		geno_var-=geno_mean*geno_mean;
//		geno_var=geno_mean*(1-geno_mean*0.5);
		
		for (size_t i=0; i<ni_total; ++i) {
			d=gsl_vector_get(geno,i);
			if (d==-9.0) {gsl_vector_set(geno, i, geno_mean);}
		}		
		
		gsl_vector_add_constant (geno, -1.0*geno_mean);
		
		if (geno_var!=0) {
			if (k_mode==1) {gsl_blas_dsyr (CblasUpper, 1.0, geno, matrix_kin);}
			else if (k_mode==2) {gsl_blas_dsyr (CblasUpper, 1.0/geno_var, geno, matrix_kin);}
			else {cout<<"Unknown kinship mode."<<endl;}
		}
		
		ns_test++;
    }	
	cout<<endl;
	
	gsl_matrix_scale (matrix_kin, 1.0/(double)ns_test);
	
	for (size_t i=0; i<ni_total; ++i) {
		for (size_t j=0; j<i; ++j) {
			d=gsl_matrix_get (matrix_kin, j, i);
			gsl_matrix_set (matrix_kin, i, j, d);
		}
	}
	
	gsl_vector_free (geno);
	
	infile.close();
	infile.clear();	
	
	return true;
}





//Read bimbam mean genotype file, the second time, recode "mean" genotype and calculate K
bool ReadFile_geno (const string &file_geno, vector<int> &indicator_idv, vector<int> &indicator_snp, gsl_matrix *UtX, gsl_matrix *K, const bool calc_K)
{
	igzstream infile (file_geno.c_str(), igzstream::in);
//	ifstream infile (file_geno.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	if (calc_K==true) {gsl_matrix_set_zero (K);}
	
	gsl_vector *genotype=gsl_vector_alloc (UtX->size1);
	gsl_vector *genotype_miss=gsl_vector_alloc (UtX->size1);
	double geno, geno_mean;
	size_t n_miss;
	
	int ni_total=(int)indicator_idv.size();
	int ns_total=(int)indicator_snp.size();
	int ni_test=UtX->size1;
	int ns_test=UtX->size2;
	
	int c_idv=0, c_snp=0;
	
	for (int i=0; i<ns_total; ++i) {
		!safeGetline(infile, line).eof();
		if (indicator_snp[i]==0) {continue;}	
		
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		ch_ptr=strtok (NULL, " , \t");
		ch_ptr=strtok (NULL, " , \t");
		
		c_idv=0; geno_mean=0; n_miss=0;
		gsl_vector_set_zero (genotype_miss);
		for (int j=0; j<ni_total; ++j) {
			ch_ptr=strtok (NULL, " , \t");
			if (indicator_idv[j]==0) {continue;}			
			
			if (strcmp(ch_ptr, "NA")==0) {gsl_vector_set (genotype_miss, c_idv, 1); n_miss++;}
			else {			
				geno=atof(ch_ptr);
				gsl_vector_set (genotype, c_idv, geno); 
				geno_mean+=geno;
			}
			c_idv++;
		}
		
		geno_mean/=(double)(ni_test-n_miss);
		
		for (size_t i=0; i<genotype->size; ++i) {			
			if (gsl_vector_get (genotype_miss, i)==1) {geno=0;}
			else {geno=gsl_vector_get (genotype, i); geno-=geno_mean;}
			
			gsl_vector_set (genotype, i, geno);
			gsl_matrix_set (UtX, i, c_snp, geno);
		}
		
		if (calc_K==true) {gsl_blas_dsyr (CblasUpper, 1.0, genotype, K);}
		
		c_snp++;
	}	
	
	if (calc_K==true) {
		gsl_matrix_scale (K, 1.0/(double)ns_test);
		
		for (size_t i=0; i<genotype->size; ++i) {
			for (size_t j=0; j<i; ++j) {
				geno=gsl_matrix_get (K, j, i);
				gsl_matrix_set (K, i, j, geno);
			}
		}
	}
	
	gsl_vector_free (genotype);
	gsl_vector_free (genotype_miss);
	
	infile.clear();
	infile.close();
	
	return true;
}







//Read bimbam mean genotype file, the second time, recode "mean" genotype and calculate K
bool ReadFile_bed (const string &file_bed, vector<int> &indicator_idv, vector<int> &indicator_snp, gsl_matrix *UtX, gsl_matrix *K, const bool calc_K)
{
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}
	
	char ch[1];
	bitset<8> b;
	
	int ni_total=(int)indicator_idv.size();
	int ns_total=(int)indicator_snp.size();
	int ni_test=UtX->size1;
	int ns_test=UtX->size2;
	int n_bit;
	
	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1;}
	
	//print the first three majic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}
	
	if (calc_K==true) {gsl_matrix_set_zero (K);}
	
	gsl_vector *genotype=gsl_vector_alloc (UtX->size1);	
	
	double geno, geno_mean;
	size_t n_miss;	
	int c_idv=0, c_snp=0, c=0;
	
	//start reading snps and doing association test
	for (int t=0; t<ns_total; ++t) {
		if (indicator_snp[t]==0) {continue;}	
		infile.seekg(t*n_bit+3);		//n_bit, and 3 is the number of magic numbers
		
		//read genotypes
		c_idv=0; geno_mean=0.0; n_miss=0; c=0;
		for (int i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && c==ni_total) {break;}				
				if (indicator_idv[c]==0) {c++; continue;}
				c++;
				
				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(genotype, c_idv, 2.0); geno_mean+=2.0;}
					else {gsl_vector_set(genotype, c_idv, 1.0); geno_mean+=1.0;}
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(genotype, c_idv, 0.0); geno_mean+=0.0;}                               
					else {gsl_vector_set(genotype, c_idv, -9.0); n_miss++;}
				}
				c_idv++;
			}
		}
		
		geno_mean/=(double)(ni_test-n_miss);
		
		for (size_t i=0; i<genotype->size; ++i) {		
			geno=gsl_vector_get (genotype, i);
			if (geno==-9) {geno=0;}
			else {geno-=geno_mean;}
			
			gsl_vector_set (genotype, i, geno);
			gsl_matrix_set (UtX, i, c_snp, geno);
		}
		
		if (calc_K==true) {gsl_blas_dsyr (CblasUpper, 1.0, genotype, K);}
		
		c_snp++;
	}	
	
	if (calc_K==true) {
		gsl_matrix_scale (K, 1.0/(double)ns_test);
		
		for (size_t i=0; i<genotype->size; ++i) {
			for (size_t j=0; j<i; ++j) {
				geno=gsl_matrix_get (K, j, i);
				gsl_matrix_set (K, i, j, geno);
			}
		}
	}
	
	gsl_vector_free (genotype);		  
	infile.clear();
	infile.close();
	
	return true;
}





bool ReadFile_est (const string &file_est, const vector<size_t> &est_column, map<string, double> &mapRS2est)
{
	mapRS2est.clear();
	
	ifstream infile (file_est.c_str(), ifstream::in);
	if (!infile) {cout<<"error opening estimated parameter file: "<<file_est<<endl; return false;}
	
	string line;
	char *ch_ptr;
	
	string rs;
	double alpha, beta, gamma, d;
	
	//header
	getline(infile, line);
	
	size_t n=*max_element(est_column.begin(), est_column.end());
	
	while (getline(infile, line)) {
		ch_ptr=strtok ((char *)line.c_str(), " \t");		
		
		alpha=0.0; beta=0.0; gamma=1.0;
		for (size_t i=0; i<n+1; ++i) {
			if (i==est_column[0]-1) {rs=ch_ptr;}
			if (i==est_column[1]-1) {alpha=atof(ch_ptr);}
			if (i==est_column[2]-1) {beta=atof(ch_ptr);}
			if (i==est_column[3]-1) {gamma=atof(ch_ptr);}
			if (i<n) {ch_ptr=strtok (NULL, " \t");}
		}
		
		d=alpha+beta*gamma;
		
		if (mapRS2est.count(rs)==0) {
			mapRS2est[rs]=d;
		}
		else {
			cout<<"the same SNP occurs more than once in estimated parameter file: "<<rs<<endl; return false;
		}
	}
	
	infile.clear();
	infile.close();
	return true;
}



bool CountFileLines (const string &file_input, size_t &n_lines)
{
	igzstream infile (file_input.c_str(), igzstream::in);
	//ifstream infile (file_input.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open file: "<<file_input<<endl; return false;}

	n_lines=count(istreambuf_iterator<char>(infile), istreambuf_iterator<char>(), '\n');
	infile.seekg (0, ios::beg);
	
	return true;
}



//Read gene expression file
bool ReadFile_gene (const string &file_gene, vector<double> &vec_read, vector<SNPINFO> &snpInfo, size_t &ng_total)
{
	vec_read.clear();
	ng_total=0;
	
	ifstream infile (file_gene.c_str(), ifstream::in);
	if (!infile) {cout<<"error! fail to open gene expression file: "<<file_gene<<endl; return false;}
	
	string line;
	char *ch_ptr;
	string rs;
	
	size_t n_idv=0, t=0;
	
	//header
	getline(infile, line);
	
	while (getline(infile, line)) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		rs=ch_ptr;
		
		ch_ptr=strtok (NULL, " , \t");	
		
		t=0;
		while (ch_ptr!=NULL) {
			if (ng_total==0) {
				vec_read.push_back(0);
				t++;
				n_idv++;
			} else {
				vec_read[t]+=atof(ch_ptr);		
				t++;
			}
			
			ch_ptr=strtok (NULL, " , \t");	
		}
		
		if (t!=n_idv) {cout<<"error! number of columns doesn't match in row: "<<ng_total<<endl; return false;}
		
		SNPINFO sInfo={"-9", rs, -9, -9, "-9", "-9", -9, -9, -9};
		snpInfo.push_back(sInfo);
		
		ng_total++;
	}
	
	infile.close();
	infile.clear();	
	
	return true;
}


