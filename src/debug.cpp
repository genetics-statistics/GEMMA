
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"

#include "debug.h"
#include "mathfunc.h"

// Helper function called by macro validate_K(K, check)
void do_validate_K(const gsl_matrix *K, bool do_check, const char *__file, int __line) {
  if (do_check) {
    debug_msg("Validating K");
    if (!checkMatrixEigen(K)) warning_at_msg(__file,__line,"K has small or negative eigenvalues!");
    if (!isMatrixIllConditioned(K)) warning_at_msg(__file,__line,"K is ill conditioned!")    if (!isMatrixSymmetric(K)) fail_at_msg(__file,__line,"K is not symmetric!" );
    if (!isMatrixPositiveDefinite(K)) fail_at_msg(__file,__line,"K is not positive definite!");
;
  }
}
