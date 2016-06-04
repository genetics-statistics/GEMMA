#include "utils.h"
#include <iostream>
#include <cstring>
using namespace std;

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


//read header to determine which column contains which item
bool ReadHeader (const string &line, HEADER &header)
{
  string rs_ptr[]={"rs","RS","snp","SNP","snps","SNPS","snpid","SNPID","rsid","RSID","MarkerName"};
  set<string> rs_set(rs_ptr, rs_ptr+11);
  string chr_ptr[]={"chr","CHR"};
  set<string> chr_set(chr_ptr, chr_ptr+2);
  string pos_ptr[]={"ps","PS","pos","POS","base_position","BASE_POSITION", "bp", "BP"};
  set<string> pos_set(pos_ptr, pos_ptr+8);
  string cm_ptr[]={"cm","CM"};
  set<string> cm_set(cm_ptr, cm_ptr+2);
  string a1_ptr[]={"a1","A1","allele1","ALLELE1","Allele1","INC_ALLELE"};
  set<string> a1_set(a1_ptr, a1_ptr+5);
  string a0_ptr[]={"a0","A0","allele0","ALLELE0","Allele0","a2","A2","allele2","ALLELE2","Allele2","DEC_ALLELE"};
  set<string> a0_set(a0_ptr, a0_ptr+10);

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
  string ncase_ptr[]={"ncase","NCASE","n_case","N_CASE"};
  set<string> ncase_set(ncase_ptr, ncase_ptr+4);
  string ncontrol_ptr[]={"ncontrol","NCONTROL","n_control","N_CONTROL"};
  set<string> ncontrol_set(ncontrol_ptr, ncontrol_ptr+4);

  string af_ptr[]={"af","AF","maf","MAF","f","F","allele_freq","ALLELE_FREQ","allele_frequency","ALLELE_FREQUENCY","Freq.Allele1.HapMapCEU","FreqAllele1HapMapCEU", "Freq1.Hapmap"};
  set<string> af_set(af_ptr, af_ptr+13);
  string var_ptr[]={"var","VAR"};
  set<string> var_set(var_ptr, var_ptr+2);

  string ws_ptr[]={"window_size","WINDOW_SIZE","ws","WS"};
  set<string> ws_set(ws_ptr, ws_ptr+4);
  string cor_ptr[]={"cor","COR","r","R"};
  set<string> cor_set(cor_ptr, cor_ptr+4);

  header.rs_col=0; header.chr_col=0; header.pos_col=0;  header.cm_col=0; header.a1_col=0; header.a0_col=0; header.z_col=0; header.beta_col=0; header.sebeta_col=0; header.chisq_col=0; header.p_col=0; header.n_col=0; header.nmis_col=0; header.nobs_col=0; header.ncase_col=0; header.ncontrol_col=0; header.af_col=0; header.var_col=0; header.ws_col=0; header.cor_col=0; header.coln=0;

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
    } else if (ncase_set.count(type)!=0) {
      if (header.ncase_col==0) {header.ncase_col=header.coln+1;} else {cout<<"error! more than two n_case columns in the file."<<endl; n_error++;}
    } else if (ncontrol_set.count(type)!=0) {
      if (header.ncontrol_col==0) {header.ncontrol_col=header.coln+1;} else {cout<<"error! more than two n_control columns in the file."<<endl; n_error++;}
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
