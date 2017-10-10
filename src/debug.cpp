
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
void do_validate_K(const gsl_matrix *K, bool do_check, bool strict, const char *__file, int __line) {
  if (do_check) {
    // debug_msg("Validating K");
    auto eigenvalues = getEigenValues(K);
    const uint count_small = count_small_values(eigenvalues,EIGEN_MINVALUE);
    if (count_small>1) {
      std::string msg = "K has ";
      msg += std::to_string(count_small);
      msg += " eigenvalues close to zero";
      warning_at_msg(__file,__line,msg);
    }
    if (!isMatrixIllConditioned(eigenvalues))
      warning_at_msg(__file,__line,"K is ill conditioned!");
    if (!isMatrixSymmetric(K))
      fail_at_msg(strict,__file,__line,"K is not symmetric!" );
    bool negative_values = has_negative_values_but_one(eigenvalues);
    if (negative_values) {
      warning_at_msg(__file,__line,"K has more than one negative eigenvalues!");
    }
    if (count_small>0 && negative_values && !isMatrixPositiveDefinite(K))
      fail_at_msg(strict,__file,__line,"K is not positive definite!");
    gsl_vector_free(eigenvalues);
  }
}
