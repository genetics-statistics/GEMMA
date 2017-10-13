
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

static bool debug_mode     = false;
static bool debug_no_check = false;
static bool debug_strict   = false;

void debug_set_debug_mode(bool setting) { debug_mode = setting; }
void debug_set_no_check_mode(bool setting) {debug_no_check = setting; }
void debug_set_strict_mode(bool setting) { debug_strict = setting; }

bool is_debug_mode() { return debug_mode; };
bool is_no_check_mode() { return debug_no_check; };
bool is_strict_mode() { return debug_strict; };

/*
  Helper function to make sure gsl allocations do their job because
  gsl_matrix_alloc does not initiatize values (behaviour that changed
  in GSL2) we introduced a 'strict mode' by initializing the buffer
  with NaNs. This happens in STRICT mode without NO-CHECKS
  (i.e. -strict option).
*/
gsl_matrix *gsl_matrix_safe_alloc(size_t rows,size_t cols) {
  gsl_matrix *m = gsl_matrix_alloc(rows,cols);
  enforce_msg(m,"Not enough memory"); // just to be sure when there is no error handler set
  if (debug_strict && !debug_no_check) {
    gsl_matrix_set_all(m, nan(""));
  }
  return m;
}

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
    const bool negative_values = has_negative_values_but_one(eigenvalues);
    if (negative_values) {
      warning_at_msg(__file,__line,"K has more than one negative eigenvalues!");
    }
    if (count_small>1 && negative_values && !isMatrixPositiveDefinite(K))
      fail_at_msg(strict,__file,__line,"K is not positive definite!");
    gsl_vector_free(eigenvalues);
  }
}
