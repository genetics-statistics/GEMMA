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
#include "eigenlib.h"

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

		SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, 0, -9, -9, 0, 0, 0};
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
	//eigenlib_dgemm("T", "N", 1.0, W, W, 0.0, WtW);
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
	size_t file_pos;

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

	file_pos=0;
	while (!safeGetline(infile, line).eof()) {
		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		rs=ch_ptr;
		ch_ptr=strtok (NULL, " , \t");
		minor=ch_ptr;
		ch_ptr=strtok (NULL, " , \t");
		major=ch_ptr;

		if (setSnps.size()!=0 && setSnps.count(rs)==0) {
		  SNPINFO sInfo={"-9", rs, -9, -9, minor, major, 0, -9, -9, 0, 0, file_pos};
		  snpInfo.push_back(sInfo);
		  indicator_snp.push_back(0);

		  file_pos++;
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

		SNPINFO sInfo={chr, rs, cM, b_pos, minor, major, n_miss, (double)n_miss/(double)ni_test, maf, ni_test-n_miss, 0, file_pos};
		snpInfo.push_back(sInfo);
		file_pos++;

		if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}

		if ( (maf<maf_level || maf> (1.0-maf_level)) && maf_level!=-1 ) {indicator_snp.push_back(0); continue;}

		if (flag_poly!=1) {indicator_snp.push_back(0); continue;}

		if (hwe_level!=0 && maf_level!=-1) {
			if (CalcHWE(n_0, n_2, n_1)<hwe_level) {indicator_snp.push_back(0); continue;}
		}

		//filter SNP if it is correlated with W
		//unless W has only one column, of 1s
		for (size_t i=0; i<genotype->size; ++i) {
			if (gsl_vector_get (genotype_miss, i)==1) {geno=maf*2.0; gsl_vector_set (genotype, i, geno);}
		}

		gsl_blas_dgemv (CblasTrans, 1.0, W, genotype, 0.0, Wtx);
		gsl_blas_dgemv (CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
		gsl_blas_ddot (genotype, genotype, &v_x);
		gsl_blas_ddot (Wtx, WtWiWtx, &v_w);

		if (W->size2!=1 && v_w/v_x >= r2_level) {indicator_snp.push_back(0); continue;}

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
			snpInfo[t].file_position=t;
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
		snpInfo[t].n_idv=ni_test-n_miss;
		snpInfo[t].n_nb=0;
		snpInfo[t].file_position=t;

		if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}

		if ( (maf<maf_level || maf> (1.0-maf_level)) && maf_level!=-1 ) {indicator_snp.push_back(0); continue;}

		if ( (n_0+n_1)==0 || (n_1+n_2)==0 || (n_2+n_0)==0) {indicator_snp.push_back(0); continue;}

		if (hwe_level!=1 && maf_level!=-1) {
			if (CalcHWE(n_0, n_2, n_1)<hwe_level) {indicator_snp.push_back(0); continue;}
		}


		//filter SNP if it is correlated with W
		//unless W has only one column, of 1s
		for (size_t i=0; i<genotype->size; ++i) {
			if (gsl_vector_get (genotype_miss, i)==1) {geno=maf*2.0; gsl_vector_set (genotype, i, geno);}
		}

		gsl_blas_dgemv (CblasTrans, 1.0, W, genotype, 0.0, Wtx);
		gsl_blas_dgemv (CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
		gsl_blas_ddot (genotype, genotype, &v_x);
		gsl_blas_ddot (Wtx, WtWiWtx, &v_w);

		if (W->size2!=1 && v_w/v_x > r2_level) {indicator_snp.push_back(0); continue;}

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





//read the genotype for one SNP; remember to read empty lines
//geno stores original genotypes without centering
//missing values are replaced by mean
bool Bimbam_ReadOneSNP (const size_t inc, const vector<int> &indicator_idv, igzstream &infile, gsl_vector *geno, double &geno_mean)
{
  size_t ni_total=indicator_idv.size();

  //  if (infile.eof()) {infile.clear();}
  //  infile.seekg(pos);

  string line;
  char *ch_ptr;
  bool flag=false;

  for (size_t i=0; i<inc; i++) {
    !safeGetline(infile, line).eof();
  }

  if (!safeGetline(infile, line).eof()) {
    ch_ptr=strtok ((char *)line.c_str(), " , \t");
    ch_ptr=strtok (NULL, " , \t");
    ch_ptr=strtok (NULL, " , \t");

    geno_mean=0.0;
    double d;
    size_t c_idv=0;
    vector<size_t> geno_miss;

    for (size_t i=0; i<ni_total; ++i) {
      ch_ptr=strtok (NULL, " , \t");
      if (indicator_idv[i]==0) {continue;}

      if (strcmp(ch_ptr, "NA")==0) {
	geno_miss.push_back(c_idv);
      } else {
	d=atof(ch_ptr);
	gsl_vector_set (geno, c_idv, d);
	geno_mean+=d;
      }
      c_idv++;
    }

    geno_mean/=(double)(c_idv-geno_miss.size() );

    for (size_t i=0; i<geno_miss.size(); ++i) {
      gsl_vector_set(geno, geno_miss[i], geno_mean);
    }
    flag=true;
  }

  return flag;
}


//for plink, store SNPs as double too
void Plink_ReadOneSNP (const int pos, const vector<int> &indicator_idv, ifstream &infile, gsl_vector *geno, double &geno_mean)
{
  size_t ni_total=indicator_idv.size(), n_bit;
  if (ni_total%4==0) {n_bit=ni_total/4;}
  else {n_bit=ni_total/4+1;}
  infile.seekg(pos*n_bit+3); //n_bit, and 3 is the number of magic numbers

  //read genotypes
  char ch[1];
  bitset<8> b;

  geno_mean=0.0;
  size_t c=0, c_idv=0;
  vector<size_t> geno_miss;

  for (size_t i=0; i<n_bit; ++i) {
    infile.read(ch,1);
    b=ch[0];
    for (size_t j=0; j<4; ++j) { //minor allele homozygous: 2.0; major: 0.0;
      if ((i==(n_bit-1)) && c==ni_total) {break;}
      if (indicator_idv[c]==0) {c++; continue;}
      c++;

      if (b[2*j]==0) {
	if (b[2*j+1]==0) {
	  gsl_vector_set (geno, c_idv, 2);
	  geno_mean+=2.0;
	} else {
	  gsl_vector_set (geno, c_idv, 1);
	  geno_mean+=1.0;
	}
      } else {
	if (b[2*j+1]==1) {
	  gsl_vector_set (geno, c_idv, 0);
	  geno_mean+=0.0;
	} else {
	  geno_miss.push_back(c_idv);
	}
      }

      c_idv++;
    }
  }

  geno_mean/=(double)(c_idv-geno_miss.size());

  for (size_t i=0; i<geno_miss.size(); ++i) {
    gsl_vector_set(geno, geno_miss[i], geno_mean);
  }

  return;
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
		  if (k_mode==1) {
		    gsl_blas_dsyr (CblasUpper, 1.0, geno, matrix_kin);
		    //eigenlib_dsyr (1.0, geno, matrix_kin);
		  } else if (k_mode==2) {
		    gsl_blas_dsyr (CblasUpper, 1.0/geno_var, geno, matrix_kin);
		    //eigenlib_dsyr (1.0/geno_var, geno, matrix_kin);
		  } else {
		    cout<<"Unknown kinship mode."<<endl;
		  }
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



//compact version of the above function, using uchar instead of gsl_matrix
bool ReadFile_geno (const string &file_geno, vector<int> &indicator_idv, vector<int> &indicator_snp, vector<vector<unsigned char> > &Xt, gsl_matrix *K, const bool calc_K, const size_t ni_test, const size_t ns_test)
{
	igzstream infile (file_geno.c_str(), igzstream::in);
    //	ifstream infile (file_geno.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return false;}

	Xt.clear();
	vector<unsigned char> Xt_row;
	for (size_t i=0; i<ni_test; i++) {
	  Xt_row.push_back(0);
	}

	string line;
	char *ch_ptr;

	if (calc_K==true) {gsl_matrix_set_zero (K);}

	gsl_vector *genotype=gsl_vector_alloc (ni_test);
	gsl_vector *genotype_miss=gsl_vector_alloc (ni_test);
	double geno, geno_mean;
	size_t n_miss;

	size_t ni_total= indicator_idv.size();
	size_t ns_total= indicator_snp.size();

	size_t c_idv=0, c_snp=0;

	for (size_t i=0; i<ns_total; ++i) {
		!safeGetline(infile, line).eof();
		if (indicator_snp[i]==0) {continue;}

		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		ch_ptr=strtok (NULL, " , \t");
		ch_ptr=strtok (NULL, " , \t");

		c_idv=0; geno_mean=0; n_miss=0;
		gsl_vector_set_zero (genotype_miss);
		for (uint j=0; j<ni_total; ++j) {
			ch_ptr=strtok (NULL, " , \t");
			if (indicator_idv[j]==0) {continue;}

			if (strcmp(ch_ptr, "NA")==0) {gsl_vector_set (genotype_miss, c_idv, 1); n_miss++;} else {
				geno=atof(ch_ptr);
				gsl_vector_set (genotype, c_idv, geno);
				geno_mean+=geno;
			}
			c_idv++;
		}

		geno_mean/=(double)(ni_test-n_miss);

		for (size_t j=0; j<genotype->size; ++j) {
			if (gsl_vector_get (genotype_miss, j)==1) {
			  geno=geno_mean;
			} else {
			  geno=gsl_vector_get (genotype, j);
			}

			Xt_row[j]=Double02ToUchar(geno);
			gsl_vector_set (genotype, j, (geno-geno_mean));
		}
		Xt.push_back(Xt_row);

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

	size_t ni_total=indicator_idv.size();
	size_t ns_total=indicator_snp.size();
	size_t ni_test=UtX->size1;
	size_t ns_test=UtX->size2;
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
	size_t c_idv=0, c_snp=0, c=0;

	//start reading snps and doing association test
	for (size_t t=0; t<ns_total; ++t) {
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




//compact version of the above function, using uchar instead of gsl_matrix
bool ReadFile_bed (const string &file_bed, vector<int> &indicator_idv, vector<int> &indicator_snp, vector<vector<unsigned char> > &Xt, gsl_matrix *K, const bool calc_K, const size_t ni_test, const size_t ns_test)
{
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}

	Xt.clear();
	vector<unsigned char> Xt_row;
	for (size_t i=0; i<ni_test; i++) {
	  Xt_row.push_back(0);
	}

	char ch[1];
	bitset<8> b;

	size_t ni_total=indicator_idv.size();
	size_t ns_total=indicator_snp.size();
	int n_bit;

	if (ni_total%4==0) {n_bit=ni_total/4;}
	else {n_bit=ni_total/4+1;}

	//print the first three majic numbers
	for (int i=0; i<3; ++i) {
		infile.read(ch,1);
		b=ch[0];
	}

	if (calc_K==true) {gsl_matrix_set_zero (K);}

	gsl_vector *genotype=gsl_vector_alloc (ni_test);

	double geno, geno_mean;
	size_t n_miss;
	size_t c_idv=0, c_snp=0, c=0;

	//start reading snps and doing association test
	for (size_t t=0; t<ns_total; ++t) {
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
			if (geno==-9) {geno=geno_mean;}

			Xt_row[i]=Double02ToUchar(geno);

			geno-=geno_mean;

			gsl_vector_set (genotype, i, geno);
		}
		Xt.push_back(Xt_row);

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

	igzstream infile (file_gene.c_str(), igzstream::in);
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

		SNPINFO sInfo={"-9", rs, -9, -9, "-9", "-9", 0, -9, -9, 0, 0, 0};
		snpInfo.push_back(sInfo);

		ng_total++;
	}

	infile.close();
	infile.clear();

	return true;
}







// WJA Added
//Read Oxford sample file
bool ReadFile_sample(const string &file_sample, vector<vector<int> > &indicator_pheno, vector<vector<double> > &pheno, const vector<size_t> &p_column, vector<int> &indicator_cvt, vector<vector<double> > &cvt, size_t &n_cvt)
{
	indicator_pheno.clear();
	pheno.clear();
	indicator_cvt.clear();

	igzstream infile (file_sample.c_str(), igzstream::in);

	if (!infile) {cout<<"error! fail to open sample file: "<<file_sample<<endl; return false;}

	string line;
	char *ch_ptr;


	string id;
	double p,d;

	vector<double> pheno_row;
	vector<int> ind_pheno_row;
	int flag_na=0;

	size_t num_cols=0;
	size_t num_p_in_file=0;
	size_t num_cvt_in_file=0;

//	size_t p_max=*max_element(p_column.begin(), p_column.end());

	map<size_t, size_t> mapP2c;
	for (size_t i=0; i<p_column.size(); i++) {
		mapP2c[p_column[i]]=i;
		pheno_row.push_back(-9);
		ind_pheno_row.push_back(0);
	}

	// read header line1
	if(!safeGetline(infile, line).eof()) {
		ch_ptr=strtok((char *)line.c_str(), " \t");
		if(strcmp(ch_ptr, "ID_1")!=0) {return false;}
		ch_ptr=strtok(NULL, " \t");
		if(strcmp(ch_ptr, "ID_2")!=0) {return false;}
		ch_ptr=strtok(NULL, " \t");
		if(strcmp(ch_ptr, "missing")!=0) {return false;}
		while (ch_ptr!=NULL) {
			num_cols++;
			ch_ptr=strtok (NULL, " \t");

		}
		num_cols--;
	}

	vector<map<uint32_t, size_t> > cvt_factor_levels;

	char col_type[num_cols];
	// read header line2
	if(!safeGetline(infile, line).eof()) {
		ch_ptr=strtok ((char *)line.c_str(), " \t");
		if(strcmp(ch_ptr, "0")!=0) {return false;}
		ch_ptr=strtok(NULL, " \t");
		if(strcmp(ch_ptr, "0")!=0) {return false;}
		ch_ptr=strtok(NULL, " \t");
		if(strcmp(ch_ptr, "0")!=0) {return false;}
		size_t it=0;
		ch_ptr=strtok (NULL, " \t");
		if(ch_ptr!=NULL)
			while(ch_ptr!=NULL){
				col_type[it++]=ch_ptr[0];
				if(ch_ptr[0]=='D') {cvt_factor_levels.push_back(map<uint32_t, size_t>());num_cvt_in_file++;}
				if(ch_ptr[0]=='C') {num_cvt_in_file++;}
				if((ch_ptr[0]=='P')||(ch_ptr[0]=='B')) {num_p_in_file++;}
				ch_ptr=strtok(NULL, " \t");
			}

	}

	while (!safeGetline(infile, line).eof()) {

		ch_ptr=strtok ((char *)line.c_str(), " \t");

		for(int it=0;it<3;it++){ch_ptr=strtok(NULL, " \t");}


		size_t i=0;
		size_t p_i=0;
		size_t fac_cvt_i=0;

		while (i<num_cols) {

			if((col_type[i]=='P')||(col_type[i]=='B'))
			{
				if (mapP2c.count(p_i+1)!=0) {
					if (strcmp(ch_ptr, "NA")==0) {ind_pheno_row[mapP2c[p_i+1]]=0; pheno_row[mapP2c[p_i+1]]=-9;}
					else {p=atof(ch_ptr); ind_pheno_row[mapP2c[p_i+1]]=1; pheno_row[mapP2c[p_i+1]]=p;}
				}
				p_i++;
			}
			if(col_type[i]=='D')
			{
				// NOTE THIS DOES NOT CHECK TO BE SURE LEVEL IS INTEGRAL i.e for atoi error
				if (strcmp(ch_ptr, "NA")!=0) {uint32_t level=atoi(ch_ptr); if(cvt_factor_levels[fac_cvt_i].count(level) == 0) {cvt_factor_levels[fac_cvt_i][level]=cvt_factor_levels[fac_cvt_i].size();}}
				fac_cvt_i++;
			}

			ch_ptr=strtok (NULL, " \t");
			i++;
		}


		indicator_pheno.push_back(ind_pheno_row);
		pheno.push_back(pheno_row);

	}
	// close and reopen the file
 	infile.close();
 	infile.clear();

	if(num_cvt_in_file>0)
	{
		igzstream infile2 (file_sample.c_str(), igzstream::in);

		if (!infile2) {cout<<"error! fail to open sample file: "<<file_sample<<endl; return false;}
		// skip header
		safeGetline(infile2, line);
		safeGetline(infile2, line);

		// pull in the covariates now we now the number of factor levels
		while (!safeGetline(infile2, line).eof()) {

			vector<double> v_d; flag_na=0;
			ch_ptr=strtok ((char *)line.c_str(), " \t");

			for(int it=0;it<3;it++){ch_ptr=strtok(NULL, " \t");}


			size_t i=0;
			size_t fac_cvt_i=0;
			size_t num_fac_levels;
			while (i<num_cols) {

				if(col_type[i]=='C')
				{
					if (strcmp(ch_ptr, "NA")==0) {flag_na=1; d=-9;}
					else {d=atof(ch_ptr);}

					v_d.push_back(d);
				}


				if(col_type[i]=='D')
				{
					// NOTE THIS DOES NOT CHECK TO BE SURE LEVEL IS INTEGRAL i.e for atoi error
					num_fac_levels=cvt_factor_levels[fac_cvt_i].size();
					if(num_fac_levels>1)
					{
						if (strcmp(ch_ptr, "NA")==0) {flag_na=1; for(size_t it=0;it<num_fac_levels-1; it++) {v_d.push_back(-9);}}
						else {uint32_t level=atoi(ch_ptr); for(size_t it=0;it<num_fac_levels-1;it++) {cvt_factor_levels[fac_cvt_i][level]==it+1 ? v_d.push_back(1.0) : v_d.push_back(0.0); }}
					}
					fac_cvt_i++;
				}

				ch_ptr=strtok (NULL, " \t");
				i++;
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

		infile2.close();
		infile2.clear();
	}
 	return true;
}



// WJA Added
//Read bgen file, the first time
#include <cstdint>
#include <assert.h>
bool ReadFile_bgen(const string &file_bgen, const set<string> &setSnps, const gsl_matrix *W, vector<int> &indicator_idv, vector<int> &indicator_snp, vector<SNPINFO> &snpInfo, const double &maf_level, const double &miss_level, const double &hwe_level, const double &r2_level, size_t &ns_test)
{

	indicator_snp.clear();

	ifstream infile (file_bgen.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bgen file:"<<file_bgen<<endl; return false;}

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

	// read in header
	uint32_t bgen_snp_block_offset;
	uint32_t bgen_header_length;
	uint32_t bgen_nsamples;
	uint32_t bgen_nsnps;
	uint32_t bgen_flags;
	infile.read(reinterpret_cast<char*>(&bgen_snp_block_offset),4);
	infile.read(reinterpret_cast<char*>(&bgen_header_length),4);
	bgen_snp_block_offset-=4;
	infile.read(reinterpret_cast<char*>(&bgen_nsnps),4);
	bgen_snp_block_offset-=4;
	infile.read(reinterpret_cast<char*>(&bgen_nsamples),4);
	bgen_snp_block_offset-=4;
	infile.ignore(4+bgen_header_length-20);
	bgen_snp_block_offset-=4+bgen_header_length-20;
	infile.read(reinterpret_cast<char*>(&bgen_flags),4);
	bgen_snp_block_offset-=4;
	bool CompressedSNPBlocks=bgen_flags&0x1;
	bool LongIds=bgen_flags&0x4;

	if(!LongIds) {return false;}

	infile.ignore(bgen_snp_block_offset);

	ns_test=0;

	size_t ns_total=static_cast<size_t>(bgen_nsnps);

	snpInfo.clear();
	string rs;
	long int b_pos;
	string chr;
//	double cM;
	string major;
	string minor;
	string id;

	double v_x, v_w;
	int c_idv=0;


	double maf, geno, geno_old;
	size_t n_miss;
	size_t n_0, n_1, n_2;
	int flag_poly;

	double bgen_geno_prob_AA, bgen_geno_prob_AB, bgen_geno_prob_BB, bgen_geno_prob_non_miss;


	size_t ni_total=indicator_idv.size();   // total number of samples in phenotype file
	size_t ni_test=0;   // number of samples to use in test

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

	for (size_t i=0; i<ni_total; ++i) {

		ni_test+=indicator_idv[i];
	}



//	ns_total=1;
	for (size_t t=0; t<ns_total; ++t) {

		id.clear();
		rs.clear();
		chr.clear();
		bgen_A_allele.clear();
		bgen_B_allele.clear();

		infile.read(reinterpret_cast<char*>(&bgen_N),4);
		infile.read(reinterpret_cast<char*>(&bgen_LS),2);

		id.resize(bgen_LS);
		infile.read(&id[0], bgen_LS);

		infile.read(reinterpret_cast<char*>(&bgen_LR),2);
		rs.resize(bgen_LR);
		infile.read(&rs[0], bgen_LR);

		infile.read(reinterpret_cast<char*>(&bgen_LC),2);
		chr.resize(bgen_LC);
		infile.read(&chr[0], bgen_LC);

		infile.read(reinterpret_cast<char*>(&bgen_SNP_pos),4);

		infile.read(reinterpret_cast<char*>(&bgen_LA),4);
		bgen_A_allele.resize(bgen_LA);
		infile.read(&bgen_A_allele[0], bgen_LA);


		infile.read(reinterpret_cast<char*>(&bgen_LB),4);
		bgen_B_allele.resize(bgen_LB);
		infile.read(&bgen_B_allele[0], bgen_LB);


		// should we switch according to MAF?
		minor=bgen_B_allele;
		major=bgen_A_allele;
		b_pos=static_cast<long int>(bgen_SNP_pos);

		uint16_t unzipped_data[3*bgen_N];

		if (setSnps.size()!=0 && setSnps.count(rs)==0) {
			SNPINFO sInfo={"-9", rs, -9, -9, minor, major, -9, -9, -9};
			snpInfo.push_back(sInfo);
			indicator_snp.push_back(0);
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


		maf=0; n_miss=0; flag_poly=0; geno_old=-9;
		n_0=0; n_1=0; n_2=0;
		c_idv=0;
		gsl_vector_set_zero (genotype_miss);
		for (size_t i=0; i<bgen_N; ++i) {
			// CHECK this set correctly!
			if (indicator_idv[i]==0) {continue;}


			bgen_geno_prob_AA=static_cast<double>(unzipped_data[i*3])/32768.0;
			bgen_geno_prob_AB=static_cast<double>(unzipped_data[i*3+1])/32768.0;
			bgen_geno_prob_BB=static_cast<double>(unzipped_data[i*3+2])/32768.0;
			bgen_geno_prob_non_miss=bgen_geno_prob_AA+bgen_geno_prob_AB+bgen_geno_prob_BB;

			//CHECK 0.1 OK
			if (bgen_geno_prob_non_miss<0.9) {gsl_vector_set (genotype_miss, c_idv, 1); n_miss++; c_idv++; continue;}


			bgen_geno_prob_AA/=bgen_geno_prob_non_miss;
			bgen_geno_prob_AB/=bgen_geno_prob_non_miss;
			bgen_geno_prob_BB/=bgen_geno_prob_non_miss;

			geno=2.0*bgen_geno_prob_BB+bgen_geno_prob_AB;
			if (geno>=0 && geno<=0.5) {n_0++;}
			if (geno>0.5 && geno<1.5) {n_1++;}
			if (geno>=1.5 && geno<=2.0) {n_2++;}

			gsl_vector_set (genotype, c_idv, geno);

			// CHECK WHAT THIS DOES
			if (flag_poly==0) {geno_old=geno; flag_poly=2;}
			if (flag_poly==2 && geno!=geno_old) {flag_poly=1;}

			maf+=geno;

			c_idv++;
		}

		maf/=2.0*static_cast<double>(ni_test-n_miss);

		SNPINFO sInfo={chr, rs, -9, b_pos, minor, major, n_miss, (double)n_miss/(double)ni_test, maf};
		snpInfo.push_back(sInfo);

		if ( (double)n_miss/(double)ni_test > miss_level) {indicator_snp.push_back(0); continue;}

		if ( (maf<maf_level || maf> (1.0-maf_level)) && maf_level!=-1 ) {indicator_snp.push_back(0); continue;}

		if (flag_poly!=1) {indicator_snp.push_back(0); continue;}

		if (hwe_level!=0 && maf_level!=-1) {
			if (CalcHWE(n_0, n_2, n_1)<hwe_level) {indicator_snp.push_back(0); continue;}
		}

		//filter SNP if it is correlated with W
		//unless W has only one column, of 1s
		for (size_t i=0; i<genotype->size; ++i) {
			if (gsl_vector_get (genotype_miss, i)==1) {geno=maf*2.0; gsl_vector_set (genotype, i, geno);}
		}

		gsl_blas_dgemv (CblasTrans, 1.0, W, genotype, 0.0, Wtx);
		gsl_blas_dgemv (CblasNoTrans, 1.0, WtWi, Wtx, 0.0, WtWiWtx);
		gsl_blas_ddot (genotype, genotype, &v_x);
		gsl_blas_ddot (Wtx, WtWiWtx, &v_w);

		if (W->size2!=1 && v_w/v_x >= r2_level) {indicator_snp.push_back(0); continue;}

		indicator_snp.push_back(1);
		ns_test++;

	}




	return true;

}


//read oxford genotype file and calculate kinship matrix
bool bgenKin (const string &file_oxford, vector<int> &indicator_snp, const int k_mode, const int display_pace, gsl_matrix *matrix_kin)
{
	string file_bgen=file_oxford;
	ifstream infile (file_bgen.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bgen file:"<<file_bgen<<endl; return false;}


	// read in header
	uint32_t bgen_snp_block_offset;
	uint32_t bgen_header_length;
	uint32_t bgen_nsamples;
	uint32_t bgen_nsnps;
	uint32_t bgen_flags;
	infile.read(reinterpret_cast<char*>(&bgen_snp_block_offset),4);
	infile.read(reinterpret_cast<char*>(&bgen_header_length),4);
	bgen_snp_block_offset-=4;
	infile.read(reinterpret_cast<char*>(&bgen_nsnps),4);
	bgen_snp_block_offset-=4;
	infile.read(reinterpret_cast<char*>(&bgen_nsamples),4);
	bgen_snp_block_offset-=4;
	infile.ignore(4+bgen_header_length-20);
	bgen_snp_block_offset-=4+bgen_header_length-20;
	infile.read(reinterpret_cast<char*>(&bgen_flags),4);
	bgen_snp_block_offset-=4;
	bool CompressedSNPBlocks=bgen_flags&0x1;
//	bool LongIds=bgen_flags&0x4;

	infile.ignore(bgen_snp_block_offset);

	double bgen_geno_prob_AA, bgen_geno_prob_AB, bgen_geno_prob_BB, bgen_geno_prob_non_miss;

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
	double genotype;


	size_t n_miss;
	double d, geno_mean, geno_var;

	size_t ni_total=matrix_kin->size1;
	gsl_vector *geno=gsl_vector_alloc (ni_total);
	gsl_vector *geno_miss=gsl_vector_alloc (ni_total);

	size_t ns_test=0;
	for (size_t t=0; t<indicator_snp.size(); ++t) {

		if (t%display_pace==0 || t==(indicator_snp.size()-1)) {ProgressBar ("Reading SNPs  ", t, indicator_snp.size()-1);}

		id.clear();
		rs.clear();
		chr.clear();
		bgen_A_allele.clear();
		bgen_B_allele.clear();

		infile.read(reinterpret_cast<char*>(&bgen_N),4);
		infile.read(reinterpret_cast<char*>(&bgen_LS),2);

		id.resize(bgen_LS);
		infile.read(&id[0], bgen_LS);

		infile.read(reinterpret_cast<char*>(&bgen_LR),2);
		rs.resize(bgen_LR);
		infile.read(&rs[0], bgen_LR);

		infile.read(reinterpret_cast<char*>(&bgen_LC),2);
		chr.resize(bgen_LC);
		infile.read(&chr[0], bgen_LC);

		infile.read(reinterpret_cast<char*>(&bgen_SNP_pos),4);

		infile.read(reinterpret_cast<char*>(&bgen_LA),4);
		bgen_A_allele.resize(bgen_LA);
		infile.read(&bgen_A_allele[0], bgen_LA);


		infile.read(reinterpret_cast<char*>(&bgen_LB),4);
		bgen_B_allele.resize(bgen_LB);
		infile.read(&bgen_B_allele[0], bgen_LB);




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



		geno_mean=0.0; n_miss=0; geno_var=0.0;
		gsl_vector_set_all(geno_miss, 0);

		for (size_t i=0; i<bgen_N; ++i) {


				bgen_geno_prob_AA=static_cast<double>(unzipped_data[i*3])/32768.0;
				bgen_geno_prob_AB=static_cast<double>(unzipped_data[i*3+1])/32768.0;
				bgen_geno_prob_BB=static_cast<double>(unzipped_data[i*3+2])/32768.0;
				// WJA
				bgen_geno_prob_non_miss=bgen_geno_prob_AA+bgen_geno_prob_AB+bgen_geno_prob_BB;
				if (bgen_geno_prob_non_miss<0.9) {gsl_vector_set(geno_miss, i, 0.0); n_miss++;}
				else {

					bgen_geno_prob_AA/=bgen_geno_prob_non_miss;
					bgen_geno_prob_AB/=bgen_geno_prob_non_miss;
					bgen_geno_prob_BB/=bgen_geno_prob_non_miss;

					genotype=2.0*bgen_geno_prob_BB+bgen_geno_prob_AB;

					gsl_vector_set(geno, i, genotype);
					gsl_vector_set(geno_miss, i, 1.0);
					geno_mean+=genotype;
					geno_var+=genotype*genotype;
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























//read header to determine which column contains which item
bool ReadHeader (const string &line, HEADER &header)
{
  string rs_ptr[]={"rs","RS","snp","SNP","snps","SNPS","snpid","SNPID","rsid","RSID"};
  set<string> rs_set(rs_ptr, rs_ptr+10);
  string chr_ptr[]={"chr","CHR"};
  set<string> chr_set(chr_ptr, chr_ptr+2);
  string pos_ptr[]={"ps","PS","pos","POS","base_position","BASE_POSITION", "bp", "BP"};
  set<string> pos_set(pos_ptr, pos_ptr+8);
  string cm_ptr[]={"cm","CM"};
  set<string> cm_set(cm_ptr, cm_ptr+2);
  string a1_ptr[]={"a1","A1","allele1","ALLELE1"};
  set<string> a1_set(a1_ptr, a1_ptr+4);
  string a0_ptr[]={"a0","A0","allele0","ALLELE0"};
  set<string> a0_set(a0_ptr, a0_ptr+4);

  string z_ptr[]={"z","Z","z_score","Z_SCORE","zscore","ZSCORE"};
  set<string> z_set(z_ptr, z_ptr+6);
  string beta_ptr[]={"beta","BETA","b","B"};
  set<string> beta_set(beta_ptr, beta_ptr+4);
  string sebeta_ptr[]={"se_beta","SE_BETA","se","SE"};
  set<string> sebeta_set(sebeta_ptr, sebeta_ptr+4);
  string chisq_ptr[]={"chisq","CHISQ","chisquare","CHISQUARE"};
  set<string> chisq_set(chisq_ptr, chisq_ptr+4);
  string p_ptr[]={"p","P","pvalue","PVALUE","p-value","P-VALUE"};
  set<string> p_set(p_ptr, p_ptr+6);

  string n_ptr[]={"n","N","ntotal","NTOTAL","n_total","N_TOTAL"};
  set<string> n_set(n_ptr, n_ptr+6);
  string nmis_ptr[]={"nmis","NMIS","n_mis","N_MIS","n_miss","N_MISS"};
  set<string> nmis_set(nmis_ptr, nmis_ptr+6);
  string nobs_ptr[]={"nobs","NOBS","n_obs","N_OBS"};
  set<string> nobs_set(nobs_ptr, nobs_ptr+4);

  string af_ptr[]={"af","AF","maf","MAF","f","F","allele_freq","ALLELE_FREQ","allele_frequency","ALLELE_FREQUENCY"};
  set<string> af_set(af_ptr, af_ptr+10);
  string var_ptr[]={"var","VAR"};
  set<string> var_set(var_ptr, var_ptr+2);

  string ws_ptr[]={"window_size","WINDOW_SIZE","ws","WS"};
  set<string> ws_set(ws_ptr, ws_ptr+4);
  string cor_ptr[]={"cor","COR","r","R"};
  set<string> cor_set(cor_ptr, cor_ptr+4);

  header.rs_col=0; header.chr_col=0; header.pos_col=0; header.a1_col=0; header.a0_col=0; header.z_col=0; header.beta_col=0; header.sebeta_col=0; header.chisq_col=0; header.p_col=0; header.n_col=0; header.nmis_col=0; header.nobs_col=0; header.af_col=0; header.var_col=0; header.ws_col=0; header.cor_col=0; header.coln=0;

  char *ch_ptr;
  string type;
  size_t n_error=0;

  ch_ptr=strtok ((char *)line.c_str(), " , \t");
  while (ch_ptr!=NULL) {
    type=ch_ptr;
    if (rs_set.count(type)!=0) {
      if (header.rs_col==0) {header.rs_col=header.coln+1;} else {cout<<"error! more than two rs columns in the file."<<endl; n_error++;}
    } else if (chr_set.count(type)!=0) {
      if (header.chr_col==0) {header.chr_col=header.coln+1;} else {cout<<"error! more than two chr columns in the file."<<endl; n_error++;}
    } else if (pos_set.count(type)!=0) {
      if (header.pos_col==0) {header.pos_col=header.coln+1;} else {cout<<"error! more than two pos columns in the file."<<endl; n_error++;}
    } else if (cm_set.count(type)!=0) {
      if (header.cm_col==0) {header.cm_col=header.coln+1;} else {cout<<"error! more than two cm columns in the file."<<endl; n_error++;}
    } else if (a1_set.count(type)!=0) {
      if (header.a1_col==0) {header.a1_col=header.coln+1;} else {cout<<"error! more than two allele1 columns in the file."<<endl; n_error++;}
    } else if (a0_set.count(type)!=0) {
      if (header.a0_col==0) {header.a0_col=header.coln+1;} else {cout<<"error! more than two allele0 columns in the file."<<endl; n_error++;}
    } else if (z_set.count(type)!=0) {
      if (header.z_col==0) {header.z_col=header.coln+1;} else {cout<<"error! more than two z columns in the file."<<endl; n_error++;}
    } else if (beta_set.count(type)!=0) {
      if (header.beta_col==0) {header.beta_col=header.coln+1;} else {cout<<"error! more than two beta columns in the file."<<endl; n_error++;}
    } else if (sebeta_set.count(type)!=0) {
      if (header.sebeta_col==0) {header.sebeta_col=header.coln+1;} else {cout<<"error! more than two se_beta columns in the file."<<endl; n_error++;}
    } else if (chisq_set.count(type)!=0) {
      if (header.chisq_col==0) {header.chisq_col=header.coln+1;} else {cout<<"error! more than two z columns in the file."<<endl; n_error++;}
    } else if (p_set.count(type)!=0) {
      if (header.p_col==0) {header.p_col=header.coln+1;} else {cout<<"error! more than two p columns in the file."<<endl; n_error++;}
    } else if (n_set.count(type)!=0) {
      if (header.n_col==0) {header.n_col=header.coln+1;} else {cout<<"error! more than two n_total columns in the file."<<endl; n_error++;}
    } else if (nmis_set.count(type)!=0) {
      if (header.nmis_col==0) {header.nmis_col=header.coln+1;} else {cout<<"error! more than two n_mis columns in the file."<<endl; n_error++;}
    } else if (nobs_set.count(type)!=0) {
      if (header.nobs_col==0) {header.nobs_col=header.coln+1;} else {cout<<"error! more than two n_obs columns in the file."<<endl; n_error++;}
    } else if (ws_set.count(type)!=0) {
      if (header.ws_col==0) {header.ws_col=header.coln+1;} else {cout<<"error! more than two window_size columns in the file."<<endl; n_error++;}
    } else if (af_set.count(type)!=0) {
      if (header.af_col==0) {header.af_col=header.coln+1;} else {cout<<"error! more than two af columns in the file."<<endl; n_error++;}
    } else if (cor_set.count(type)!=0) {
      if (header.cor_col==0) {header.cor_col=header.coln+1;} else {cout<<"error! more than two cor columns in the file."<<endl; n_error++;}
    } else {}

    ch_ptr=strtok (NULL, " , \t");
    header.coln++;
  }

  if (header.cor_col!=0 && header.cor_col!=header.coln) {cout<<"error! the cor column should be the last column."<<endl; n_error++;}

  if (header.rs_col==0) {
    if (header.chr_col!=0 && header.pos_col!=0) {
      cout<<"missing an rs column. rs id will be replaced by chr:pos"<<endl;
    } else {
      cout<<"error! missing an rs column."<<endl; n_error++;
    }
  }

  if (n_error==0) {return true;} else {return false;}
}




//read category file, record mapRS2in
//the category file does not contain a null category
//so if a snp has 0 for all categories, then it is not included in the analysis
bool ReadFile_cat (const string &file_cat, map<string, size_t> &mapRS2cat, size_t &n_vc)
{
  mapRS2cat.clear();

  igzstream infile (file_cat.c_str(), igzstream::in);
  if (!infile) {cout<<"error! fail to open category file: "<<file_cat<<endl; return false;}

  string line;
  char *ch_ptr;

  string rs, chr, a1, a0, pos, cm;
  size_t i_cat;// ns_vc=0;

  //read header
  HEADER header;
  !safeGetline(infile, line).eof();
  ReadHeader (line, header);

  //use the header to count the number of categories
  n_vc=header.coln;
  if (header.rs_col!=0) {n_vc--;}
  if (header.chr_col!=0) {n_vc--;}
  if (header.pos_col!=0) {n_vc--;}
  if (header.cm_col!=0) {n_vc--;}
  if (header.a1_col!=0) {n_vc--;}
  if (header.a0_col!=0) {n_vc--;}

  //read the following lines to record mapRS2cat
  while (!safeGetline(infile, line).eof()) {
    ch_ptr=strtok ((char *)line.c_str(), " , \t");

    i_cat=0;
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
      } else if (atoi(ch_ptr)==1 || atoi(ch_ptr)==0) {
	if (i_cat==0) {
	  if (header.rs_col==0) {
	    rs=chr+":"+pos;
	  }
	}

	if (atoi(ch_ptr)==1 && mapRS2cat.count(rs)==0) {mapRS2cat[rs]=i_cat;}
	i_cat++;
      } else {}

      ch_ptr=strtok (NULL, " , \t");
    }

    //if (mapRS2cat.count(rs)==0) {mapRS2cat[rs]=n_vc+1; ns_vc++;}
  }

  //if (ns_vc>0) {n_vc++;}

  infile.clear();
  infile.close();

  return true;
}




//read bimbam mean genotype file and calculate kinship matrix; this time, the kinship matrix is not centered, and can contain multiple K matrix
bool BimbamKin (const string &file_geno, vector<int> &indicator_idv, vector<int> &indicator_snp, const int k_mode, const int display_pace, const map<string, size_t> &mapRS2cat, map<string, double> &mapRS2var, vector<SNPINFO> &snpInfo, gsl_matrix *matrix_kin)
{
	igzstream infile (file_geno.c_str(), igzstream::in);
	//ifstream infile (file_geno.c_str(), ifstream::in);
	if (!infile) {cout<<"error reading genotype file:"<<file_geno<<endl; return false;}

	string line;
	char *ch_ptr;

	size_t n_miss;
	double d, geno_mean, geno_var;

	size_t ni_test=matrix_kin->size1;
	gsl_vector *geno=gsl_vector_alloc (ni_test);
	gsl_vector *geno_miss=gsl_vector_alloc (ni_test);

	size_t n_vc=matrix_kin->size2/ni_test, i_vc;
	string rs;
	vector<size_t> ns_vec;
	for (size_t i=0; i<n_vc; i++) {
	  ns_vec.push_back(0);
	}

	size_t ns_test=0;
	for (size_t t=0; t<indicator_snp.size(); ++t) {
		!safeGetline(infile, line).eof();
		if (t%display_pace==0 || t==(indicator_snp.size()-1)) {ProgressBar ("Reading SNPs  ", t, indicator_snp.size()-1);}
		if (indicator_snp[t]==0) {continue;}

		ch_ptr=strtok ((char *)line.c_str(), " , \t");
		ch_ptr=strtok (NULL, " , \t");
		ch_ptr=strtok (NULL, " , \t");

		rs=snpInfo[t].rs_number;//this line is new

		geno_mean=0.0; n_miss=0; geno_var=0.0;
		gsl_vector_set_all(geno_miss, 0);

		size_t j=0;
		for (size_t i=0; i<indicator_idv.size(); ++i) {
		  if (indicator_idv[i]==0) {continue;}
			ch_ptr=strtok (NULL, " , \t");
			if (strcmp(ch_ptr, "NA")==0) {gsl_vector_set(geno_miss, i, 0); n_miss++;}
			else {
				d=atof(ch_ptr);
				gsl_vector_set (geno, j, d);
				gsl_vector_set (geno_miss, j, 1);
				geno_mean+=d;
				geno_var+=d*d;
			}
			j++;
		}

		geno_mean/=(double)(ni_test-n_miss);
		geno_var+=geno_mean*geno_mean*(double)n_miss;
		geno_var/=(double)ni_test;
		geno_var-=geno_mean*geno_mean;
//		geno_var=geno_mean*(1-geno_mean*0.5);

		for (size_t i=0; i<ni_test; ++i) {
			if (gsl_vector_get (geno_miss, i)==0) {gsl_vector_set(geno, i, geno_mean);}
		}

		//this line is new; removed
		//gsl_vector_add_constant (geno, -1.0*geno_mean);

		if (geno_var!=0) {
		  mapRS2var[rs]=geno_var;

		  if (k_mode==1) {
		    if (n_vc==1 || mapRS2cat.size()==0 ) {
		      gsl_blas_dsyr (CblasUpper, 1.0, geno, matrix_kin);
                ns_vec[0]++;
		    } else if (mapRS2cat.count(rs)!=0) {
		      i_vc=mapRS2cat.at(rs);
		      ns_vec[i_vc]++;
		      gsl_matrix_view kin_sub=gsl_matrix_submatrix(matrix_kin, 0, ni_test*i_vc, ni_test, ni_test);
		      gsl_blas_dsyr (CblasUpper, 1.0, geno, &kin_sub.matrix);
		    }

		    //eigenlib_dsyr (1.0, geno, matrix_kin);
		  } else if (k_mode==2) {
		    if (n_vc==1 || mapRS2cat.size()==0 ) {
		      gsl_blas_dsyr (CblasUpper, 1.0/geno_var, geno, matrix_kin);
                ns_vec[0]++;
		    } else if (mapRS2cat.count(rs)!=0) {
		      i_vc=mapRS2cat.at(rs);
		      ns_vec[i_vc]++;
		      gsl_matrix_view kin_sub=gsl_matrix_submatrix(matrix_kin, 0, ni_test*i_vc, ni_test, ni_test);
		      gsl_blas_dsyr (CblasUpper, 1.0/geno_var, geno, &kin_sub.matrix);
		    }
		  } else {
		    cout<<"Unknown kinship mode."<<endl;
		  }
		}
		ns_test++;
    }
	cout<<endl;

	for (size_t t=0; t<n_vc; t++) {
	  if (ns_vec[t]!=0) {gsl_matrix_scale (matrix_kin, 1.0/(double)ns_vec[t]);}

	  for (size_t i=0; i<ni_test; ++i) {
	    for (size_t j=0; j<i; ++j) {
	      d=gsl_matrix_get (matrix_kin, j, i+ni_test*t);
	      gsl_matrix_set (matrix_kin, i, j+ni_test*t, d);
	    }
	  }
	}

	gsl_vector_free (geno);
	gsl_vector_free (geno_miss);

	infile.close();
	infile.clear();

	return true;
}







bool PlinkKin (const string &file_bed, vector<int> &indicator_idv, vector<int> &indicator_snp, const int k_mode, const int display_pace, const map<string, size_t> &mapRS2cat, map<string, double> &mapRS2var, vector<SNPINFO> &snpInfo, gsl_matrix *matrix_kin)
{
	ifstream infile (file_bed.c_str(), ios::binary);
	if (!infile) {cout<<"error reading bed file:"<<file_bed<<endl; return false;}

	char ch[1];
	bitset<8> b;

	size_t n_miss, ci_total, ci_test;
	double d, geno_mean, geno_var;

	size_t ni_test=matrix_kin->size1;
	size_t ni_total=indicator_idv.size();
	gsl_vector *geno=gsl_vector_alloc (ni_test);

	size_t ns_test=0;
	int n_bit;

	size_t n_vc=matrix_kin->size2/ni_test, i_vc;
	string rs;
	vector<size_t> ns_vec;
	for (size_t i=0; i<n_vc; i++) {
	  ns_vec.push_back(0);
	}

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

		rs=snpInfo[t].rs_number;//this line is new

		//read genotypes
		geno_mean=0.0;	n_miss=0; ci_total=0; geno_var=0.0; ci_test=0;
		for (int i=0; i<n_bit; ++i) {
			infile.read(ch,1);
			b=ch[0];
			for (size_t j=0; j<4; ++j) {                //minor allele homozygous: 2.0; major: 0.0;
				if ((i==(n_bit-1)) && ci_total==ni_total) {break;}
				if (indicator_idv[ci_total]==0) {ci_total++; continue;}

				if (b[2*j]==0) {
					if (b[2*j+1]==0) {gsl_vector_set(geno, ci_test, 2.0); geno_mean+=2.0; geno_var+=4.0; }
					else {gsl_vector_set(geno, ci_test, 1.0); geno_mean+=1.0; geno_var+=1.0;}
				}
				else {
					if (b[2*j+1]==1) {gsl_vector_set(geno, ci_test, 0.0); }
					else {gsl_vector_set(geno, ci_test, -9.0); n_miss++; }
				}

				ci_test++;
				ci_total++;
			}
		}


		geno_mean/=(double)(ni_test-n_miss);
		geno_var+=geno_mean*geno_mean*(double)n_miss;
		geno_var/=(double)ni_test;
		geno_var-=geno_mean*geno_mean;
//		geno_var=geno_mean*(1-geno_mean*0.5);

		for (size_t i=0; i<ni_test; ++i) {
			d=gsl_vector_get(geno,i);
			if (d==-9.0) {gsl_vector_set(geno, i, geno_mean);}
		}

		//this line is new; removed
		//gsl_vector_add_constant (geno, -1.0*geno_mean);

		if (geno_var!=0) {
		  mapRS2var[rs]=geno_var;
			if (k_mode==1) {
			  if (n_vc==1 || mapRS2cat.size()==0 ) {
			    gsl_blas_dsyr (CblasUpper, 1.0, geno, matrix_kin);
			    ns_vec[0]++;
			  } else if (mapRS2cat.count(rs)!=0) {
			    i_vc=mapRS2cat.at(rs);
			    ns_vec[i_vc]++;
			    gsl_matrix_view kin_sub=gsl_matrix_submatrix(matrix_kin, 0, ni_test*i_vc, ni_test, ni_test);
			    gsl_blas_dsyr (CblasUpper, 1.0, geno, &kin_sub.matrix);
			  }
			} else if (k_mode==2) {
			  if (n_vc==1 || mapRS2cat.size()==0 ) {
			    gsl_blas_dsyr (CblasUpper, 1.0/geno_var, geno, matrix_kin);
			    ns_vec[0]++;
			  } else if (mapRS2cat.count(rs)!=0) {
			    i_vc=mapRS2cat.at(rs);
			    ns_vec[i_vc]++;
			    gsl_matrix_view kin_sub=gsl_matrix_submatrix(matrix_kin, 0, ni_test*i_vc, ni_test, ni_test);
			    gsl_blas_dsyr (CblasUpper, 1.0/geno_var, geno, &kin_sub.matrix);
			  }
			} else {
			  cout<<"Unknown kinship mode."<<endl;
			}
		}

		ns_test++;
    }
	cout<<endl;

	for (size_t t=0; t<n_vc; t++) {
	  if (ns_vec[t]!=0) {gsl_matrix_scale (matrix_kin, 1.0/(double)ns_vec[t]);}

	  for (size_t i=0; i<ni_test; ++i) {
	    for (size_t j=0; j<i; ++j) {
	      d=gsl_matrix_get (matrix_kin, j, i+ni_test*t);
	      gsl_matrix_set (matrix_kin, i, j+ni_test*t, d);
	      //cout<<d<<" ";
	    }
	    //cout<<endl;
	  }
	}

	d=0;
	for (size_t i=0; i<ni_test; ++i) {
	  for (size_t j=0; j<ni_test; ++j) {
	    d+=gsl_matrix_get (matrix_kin, i, j)*gsl_matrix_get (matrix_kin, i, j);
	  }
	}
	d/=(double)ni_test*(double)ni_test;
	//cout<<"trace = "<<scientific<<d-1/(double)ni_test<<endl;



	gsl_vector_free (geno);

	infile.close();
	infile.clear();

	return true;
}



//read var file, store mapRS2var
bool ReadFile_var (const string &file_var, map<string, double> &mapRS2var)
{
  mapRS2var.clear();

  igzstream infile (file_var.c_str(), igzstream::in);
  if (!infile) {cout<<"error! fail to open var file: "<<file_var<<endl; return false;}

  char *ch_ptr;
  string line, rs;
  double var;

  while (!safeGetline(infile, line).eof()) {
    ch_ptr=strtok ((char *)line.c_str(), " , \t");
    rs=ch_ptr;
    ch_ptr=strtok (NULL, " , \t");
    var=atof(ch_ptr);
    mapRS2var[rs]=var;
  }

  return true;
}


//read beta file, use the mapRS2var to select snps (and to provide var if maf/var is not provided in the beta file), calculate q
void ReadFile_beta (const string &file_beta, const int k_mode, const map<string, size_t> &mapRS2cat, const map<string, double> &mapRS2var, gsl_vector *q, gsl_vector *s, size_t &ni_total, size_t &ns_total, size_t &ns_test)
{
  gsl_vector_set_zero(q);
  ni_total=0; ns_total=0; ns_test=0;

  igzstream infile (file_beta.c_str(), igzstream::in);
  if (!infile) {cout<<"error! fail to open beta file: "<<file_beta<<endl; return;}

  string line;
  char *ch_ptr;
  string type;

  string rs, chr, a1, a0, pos, cm;
  double z=0, beta=0, se_beta=0, chisq=0, pvalue=0, zsquare=0, af=0, var_x=0;
  size_t n_total=0, n_mis=0, n_obs=0;

  vector<double> vec_q, vec_s;
  for (size_t i=0; i<q->size; i++) {
    vec_q.push_back(0.0);
    vec_s.push_back(0.0);
  }

  //read header
  HEADER header;
  !safeGetline(infile, line).eof();
  ReadHeader (line, header);

  if (header.n_col==0 ) {
    if (header.nobs_col==0 && header.nmis_col==0) {
      cout<<"error! missing sample size in the beta file."<<endl;
    } else {
      cout<<"total sample size will be replaced by obs/mis sample size."<<endl;
    }
  }

  if (header.z_col==0 && (header.beta_col==0 || header.sebeta_col==0) && header.chisq_col==0 && header.p_col==0) {
    cout<<"error! missing z scores in the beta file."<<endl;
  }

  if (header.af_col==0 && header.var_col==0 && mapRS2var.size()==0) {
    cout<<"error! missing allele frequency in the beta file."<<endl;
  }

  while (!safeGetline(infile, line).eof()) {
    ch_ptr=strtok ((char *)line.c_str(), " , \t");

    z=0; beta=0; se_beta=0; chisq=0; pvalue=0;
    n_total=0; n_mis=0; n_obs=0; af=0; var_x=0;
    for (size_t i=0; i<header.coln; i++) {
      if (header.rs_col!=0 && header.rs_col==i+1) {rs=ch_ptr;}
      if (header.chr_col!=0 && header.chr_col==i+1) {chr=ch_ptr;}
      if (header.pos_col!=0 && header.pos_col==i+1) {pos=ch_ptr;}
      if (header.cm_col!=0 && header.cm_col==i+1) {cm=ch_ptr;}
      if (header.a1_col!=0 && header.a1_col==i+1) {a1=ch_ptr;}
      if (header.a0_col!=0 && header.a0_col==i+1) {a0=ch_ptr;}

      if (header.z_col!=0 && header.z_col==i+1) {z=atof(ch_ptr);}
      if (header.beta_col!=0 && header.beta_col==i+1) {beta=atof(ch_ptr);}
      if (header.sebeta_col!=0 && header.sebeta_col==i+1) {se_beta=atof(ch_ptr);}
      if (header.chisq_col!=0 && header.chisq_col==i+1) {chisq=atof(ch_ptr);}
      if (header.p_col!=0 && header.p_col==i+1) {pvalue=atof(ch_ptr);}

      if (header.n_col!=0 && header.n_col==i+1) {n_total=atoi(ch_ptr);}
      if (header.nmis_col!=0 && header.nmis_col==i+1) {n_mis=atoi(ch_ptr);}
      if (header.nobs_col!=0 && header.nobs_col==i+1) {n_obs=atoi(ch_ptr);}

      if (header.af_col!=0 && header.af_col==i+1) {af=atof(ch_ptr);}
      if (header.var_col!=0 && header.var_col==i+1) {var_x=atof(ch_ptr);}

      ch_ptr=strtok (NULL, " , \t");
    }

    if (header.rs_col==0) {
      rs=chr+":"+pos;
    }

    if (header.n_col==0) {
      n_total=n_mis+n_obs;
    }

    //both z values and beta/se_beta have directions, while chisq/pvalue do not
    if (header.z_col!=0) {
      zsquare=z*z;
    } else if (header.beta_col!=0 && header.sebeta_col!=0) {
      z=beta/se_beta;
      zsquare=z*z;
    } else if (header.chisq_col!=0) {
      zsquare=chisq;
    } else if (header.p_col!=0) {
      zsquare=gsl_cdf_chisq_Qinv (pvalue, 1);
    } else {zsquare=0;}

    //if the snp is also present in cor file, then do calculations
    if (mapRS2var.count(rs)!=0 && (mapRS2cat.size()==0 || mapRS2cat.count(rs)!=0) ) {
      //obtain var_x
      if (k_mode==1) {
	if (header.var_col==0) {
	  if (header.af_col!=0) {
	    var_x=2.0*af*(1.0-af);
	  } else {
	    var_x=mapRS2var.at(rs);
	  }
	}
      } else {
	var_x=1.0;
      }

      //compute q
      if (mapRS2cat.size()!=0) {
	vec_q[mapRS2cat.at(rs) ]+=(zsquare-1.0)*var_x/(double)n_total;
	vec_s[mapRS2cat.at(rs) ]+=var_x;
      } else {
	vec_q[0]+=(zsquare-1.0)*var_x/(double)n_total;
	vec_s[0]+=var_x;
      }

      ni_total=max(ni_total, n_total);
      ns_test++;
    }

    ns_total++;
  }

  //save q
  for (size_t i=0; i<q->size; i++) {
    if (vec_s[i]!=0) {
      gsl_vector_set(q, i, vec_q[i]/vec_s[i]);
    }
    gsl_vector_set(s, i, vec_s[i]);
  }

  infile.clear();
  infile.close();

  return;
}




//read S file: S and Svar
void ReadFile_s (const string &file_s, gsl_matrix *S, gsl_matrix *Svar)
{
  igzstream infile (file_s.c_str(), igzstream::in);
  if (!infile) {cout<<"error! fail to open s file: "<<file_s<<endl; return;}

  string line;
  char *ch_ptr;
  double d;

  for (size_t i=0; i<S->size1; i++) {
    !safeGetline(infile, line).eof();
    ch_ptr=strtok ((char *)line.c_str(), " , \t");
    for (size_t j=0; j<S->size2; j++) {
      d=gsl_matrix_get(S, i, j)+atof(ch_ptr);
      gsl_matrix_set(S, i, j, d);
      ch_ptr=strtok (NULL, " , \t");
    }
  }

  for (size_t i=0; i<Svar->size1; i++) {
    !safeGetline(infile, line).eof();
    ch_ptr=strtok ((char *)line.c_str(), " , \t");
    for (size_t j=0; j<Svar->size2; j++) {
      d=gsl_matrix_get(Svar, i, j)+atof(ch_ptr);
      gsl_matrix_set(Svar, i, j, d);
      ch_ptr=strtok (NULL, " , \t");
    }
  }

  infile.clear();
  infile.close();

  return;
}




void ReadFile_ms (const string &file_ms, gsl_matrix *S, gsl_matrix *Svar)
{
  gsl_matrix_set_zero(S);
  gsl_matrix_set_zero(Svar);

  string file_name;

  igzstream infile (file_ms.c_str(), igzstream::in);
  if (!infile) {cout<<"error! fail to open ms file: "<<file_ms<<endl; return;}

  while (!safeGetline(infile, file_name).eof()) {
    ReadFile_s(file_name, S, Svar);
  }

  infile.clear();
  infile.close();

  return;
}




//read V file: V (i.e. Q)
void ReadFile_v (const string &file_v, gsl_matrix *V)
{
  igzstream infile (file_v.c_str(), igzstream::in);
  if (!infile) {cout<<"error! fail to open v file: "<<file_v<<endl; return;}

  string line;
  char *ch_ptr;
  double d;

  for (size_t i=0; i<V->size1; i++) {
    !safeGetline(infile, line).eof();
    ch_ptr=strtok ((char *)line.c_str(), " , \t");
    for (size_t j=0; j<V->size2; j++) {
      d=gsl_matrix_get(V, i, j)+atof(ch_ptr);
      gsl_matrix_set(V, i, j, d);
      ch_ptr=strtok (NULL, " , \t");
    }
  }

  infile.clear();
  infile.close();

  return;
}


void ReadFile_mv (const string &file_mv, gsl_matrix *V)
{
  gsl_matrix_set_zero(V);

  string file_name;

  igzstream infile (file_mv.c_str(), igzstream::in);
  if (!infile) {cout<<"error! fail to open ms file: "<<file_mv<<endl; return;}

  while (!safeGetline(infile, file_name).eof()) {
    ReadFile_v(file_name, V);
  }

  infile.clear();
  infile.close();

  return;
}


//read q file: q, s and ni_test
void ReadFile_q (const string &file_s, gsl_vector *q_vec, gsl_vector *s_vec, double &df)
{
  igzstream infile (file_s.c_str(), igzstream::in);
  if (!infile) {cout<<"error! fail to open s file: "<<file_s<<endl; return;}

  string line;
  char *ch_ptr;
  double d;

  for (size_t i=0; i<q_vec->size; i++) {
    !safeGetline(infile, line).eof();
    ch_ptr=strtok ((char *)line.c_str(), " , \t");
    d=gsl_vector_get(q_vec, i)+atof(ch_ptr);
    gsl_vector_set(q_vec, i, d);
  }

  for (size_t i=0; i<s_vec->size; i++) {
    !safeGetline(infile, line).eof();
    ch_ptr=strtok ((char *)line.c_str(), " , \t");
    d=gsl_vector_get(s_vec, i)+atof(ch_ptr);
    gsl_vector_set(s_vec, i, d);
  }

  !safeGetline(infile, line).eof();
  ch_ptr=strtok ((char *)line.c_str(), " , \t");
  df=atof(ch_ptr);

  infile.clear();
  infile.close();

  return;
}



void ReadFile_mq (const string &file_mq, gsl_vector *q_vec, gsl_vector *s_vec, double &df)
{
  gsl_vector_set_zero(q_vec);
  gsl_vector_set_zero(s_vec);

  string file_name;

  igzstream infile (file_mq.c_str(), igzstream::in);
  if (!infile) {cout<<"error! fail to open mq file: "<<file_mq<<endl; return;}

  while (!safeGetline(infile, file_name).eof()) {
    ReadFile_q(file_name, q_vec, s_vec, df);
  }

  infile.clear();
  infile.close();

  return;
}
