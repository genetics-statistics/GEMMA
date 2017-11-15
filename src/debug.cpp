
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
static bool debug_check    = true;  // check data/algorithms
static bool debug_strict   = false; // fail on error
static bool debug_quiet    = false;
static uint debug_issue    = 0;     // track github issues
static bool debug_legacy   = false; // legacy mode

void debug_set_debug_mode(bool setting) { debug_mode = setting; }
void debug_set_no_check_mode(bool setting) {debug_check = !setting; }
void debug_set_strict_mode(bool setting) { debug_strict = setting; }
void debug_set_quiet_mode(bool setting) { debug_quiet = setting; }
void debug_set_issue(uint issue) { debug_issue = issue; }
void debug_set_legacy_mode(bool setting) { debug_legacy = setting; }

bool is_debug_mode() { return debug_mode; };
bool is_no_check_mode() { return !debug_check; };
bool is_check_mode() { return debug_check; };
bool is_strict_mode() { return debug_strict; };
bool is_quiet_mode() { return debug_quiet; };
bool is_issue(uint issue) { return issue == debug_issue; };
bool is_legacy_mode() { return debug_legacy; };

/*
  Helper function to make sure gsl allocations do their job because
  gsl_matrix_alloc does not initiatize values (behaviour that changed
  in GSL2) we introduced a 'strict mode' by initializing the buffer
  with NaNs. This happens when NO-CHECKS is not set (default) and with
  DEBUG (i.e. -debug option).
*/
gsl_matrix *gsl_matrix_safe_alloc(size_t rows,size_t cols) {
  gsl_matrix *m = gsl_matrix_alloc(rows,cols);
  enforce_msg(m,"Not enough memory"); // just to be sure when there is no error handler set
  if (is_check_mode() && is_debug_mode()) {
    gsl_matrix_set_all(m, nan(""));
  }
  return m;
}

int gsl_matrix_safe_memcpy (gsl_matrix *dest, const gsl_matrix *src) {
  enforce(dest->size1 == src->size1);
  enforce(dest->size2 == src->size2);
  return gsl_matrix_memcpy(dest,src);
}

void do_gsl_matrix_safe_free (gsl_matrix *m, const char *__pretty_function, const char *__file, int __line) {
  enforce(m);
  if (is_check_mode() && is_debug_mode()) {
    bool has_NaN = has_nan(m);
    bool has_Inf = has_inf(m);
    if (has_NaN || has_Inf) {
      std::string msg = "Matrix (size ";
      msg += std::to_string(m->size1);
      msg += "x";
      msg += std::to_string(m->size2);
      msg += ")";
      if (has_Inf)
        warnfail_at_msg(is_strict_mode(),__pretty_function,__file,__line,(msg+" contains Infinite on free!").c_str());
      if (has_NaN)
        warnfail_at_msg(is_strict_mode(),__pretty_function,__file,__line,(msg+" contains NaN on free!").c_str());
    }
  }
  return gsl_matrix_free(m);
}

int gsl_vector_safe_memcpy (gsl_vector *dest, const gsl_vector *src) {
  enforce(dest->size == src->size);
  return gsl_vector_memcpy(dest,src);
}

void do_gsl_vector_safe_free (gsl_vector *v, const char *__pretty_function, const char *__file, int __line) {
  enforce(v);
  if (is_check_mode() && is_debug_mode()) {
    bool has_NaN = has_nan(v);
    bool has_Inf = has_inf(v);
    if (has_NaN || has_Inf) {
      std::string msg = "Vector (size ";
      msg += std::to_string(v->size);
      msg += ")";
      if (has_Inf)
        warnfail_at_msg(is_strict_mode(),__pretty_function,__file,__line,(msg+" contains Infinite on free!").c_str());
      if (has_NaN)
        warnfail_at_msg(is_strict_mode(),__pretty_function,__file,__line,(msg+" contains NaN on free!").c_str());
    }
  }
  return gsl_vector_free(v);
}

/*
  Helper function to make sure gsl allocations do their job because
  gsl_vector_alloc does not initiatize values (behaviour that changed
  in GSL2) we introduced a 'strict mode' by initializing the buffer
  with NaNs. This happens when NO-CHECKS is not set and with DEBUG
  (i.e. -debug option).
*/
gsl_vector *gsl_vector_safe_alloc(size_t n) {
  gsl_vector *v = gsl_vector_alloc(n);
  enforce_msg(v,"Not enough memory"); // just to be sure when there is no error handler set
  if (is_check_mode() && is_debug_mode()) {
    gsl_vector_set_all(v, nan(""));
  }
  return v;
}

char *do_strtok_safe(char *tokenize, const char *delimiters, const char *__pretty_function, const char *__file, int __line) {
  auto token = strtok(tokenize,delimiters);
  if (token == NULL && (is_debug_mode() || is_strict_mode()))
    fail_at_msg(__file,__line,string("strtok failed in ") + __pretty_function);
  return token;
}

// Helper function called by macro validate_K(K, check). K is validated
// unless -no-check option is used.
void do_validate_K(const gsl_matrix *K, const char *__pretty_function, const char *__file, int __line) {
  if (is_check_mode()) {
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
      warnfail_at_msg(is_strict_mode(),__pretty_function,__file,__line,"K is not symmetric!" );
    const bool negative_values = has_negative_values_but_one(eigenvalues);
    if (negative_values) {
      warning_at_msg(__file,__line,"K has more than one negative eigenvalues!");
    }
    if (count_small>1 && negative_values && !isMatrixPositiveDefinite(K))
      warnfail_at_msg(is_strict_mode(),__pretty_function,__file,__line,"K is not positive definite!");
    gsl_vector_free(eigenvalues);
  }
}
