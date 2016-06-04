#include "utils.h"
#include <iostream>

//map a number 1-(n_cvt+2) to an index between 0 and [(n_c+2)^2+(n_c+2)]/2-1
size_t GetabIndex (const size_t a, const size_t b, const size_t n_cvt) {
	if (a>n_cvt+2 || b>n_cvt+2 || a<=0 || b<=0) {std::cout<<"error in GetabIndex."<<std::endl; return 0;}
	size_t index;
	size_t l, h;
	if (b>a) {l=a; h=b;} else {l=b; h=a;}

	size_t n=n_cvt+2;
	index=(2*n-l+2)*(l-1)/2+h-l;

	return index;
}
