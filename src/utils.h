#ifndef __UTILS_H__
#define __UTILS_H__

#include <cstdlib>

//map a number 1-(n_cvt+2) to an index between 0 and [(n_c+2)^2+(n_c+2)]/2-1
size_t GetabIndex (const size_t a, const size_t b, const size_t n_cvt);

#endif // __UTILS_H__
