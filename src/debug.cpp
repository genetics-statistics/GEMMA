/*
    Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright © 2011-2017, Xiang Zhou
    Copyright © 2017, Peter Carbonetto
    Copyright © 2017, Pjotr Prins

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <vector>
#include <fenv.h>

#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"

#include "debug.h"
#include "mathfunc.h"

static bool debug_mode      = false;
static bool debug_data_mode = false;
static bool debug_check     = false;  // check data/algorithms
static bool debug_fpe_check = true;  // check floating point errors (intel hardware)
static bool debug_strict    = false; // fail on error, more rigorous checks
static bool debug_quiet     = false;
static uint debug_issue     = 0;     // track github issues
static bool debug_legacy    = false; // legacy mode

void debug_set_debug_mode(bool setting) { debug_mode = setting; }
void debug_set_debug_data_mode(bool setting) { debug_data_mode = setting; }
void debug_set_check_mode(bool setting) {debug_check = setting; }
void debug_set_no_check_mode(bool setting) {debug_check = !setting; }
void debug_set_no_fpe_check_mode(bool setting) {debug_fpe_check = !setting; }
void debug_set_strict_mode(bool setting) { debug_strict = setting; }
void debug_set_quiet_mode(bool setting) { debug_quiet = setting; }
void debug_set_issue(uint issue) { debug_issue = issue; }
void debug_set_legacy_mode(bool setting) { debug_legacy = setting; }

bool is_debug_mode() { return debug_mode; };
bool is_debug_data_mode() { return debug_data_mode; };
bool is_no_check_mode() { return !debug_check; };
bool is_check_mode() { return debug_check; };
bool is_fpe_check_mode() { return debug_fpe_check; };
bool is_strict_mode() { return debug_strict; };
bool is_quiet_mode() { return debug_quiet; };
bool is_issue(uint issue) { return issue == debug_issue; };
bool is_legacy_mode() { return debug_legacy; };

#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>

void sighandler(int signum)
{
  cout << R"(
FATAL ERROR: GEMMA caused a floating point error which suggests machine boundaries were reached.

You can disable floating point tests with the -no-check switch (use at your own risk!)
)" << endl;
  signal(signum, SIG_DFL);
  kill(getpid(), signum); // should force a core dump
}

/*
   Force the floating point processor to throw an exception when the result of
   a double/float computation is overflow, underflow, NaN or inf. In principle
   this is an Intel hardware feature that does not slow down computations.
*/

#if defined(__APPLE__) && defined(__MACH__)

// Public domain polyfill for feenableexcept on OS X
// http://www-personal.umich.edu/~williams/archive/computation/fe-handling-example.c

inline int feenableexcept(unsigned int excepts)
{
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
    // previous masks
    unsigned int old_excepts;

    if (fegetenv(&fenv)) {
        return -1;
    }
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    // unmask
    fenv.__control &= ~new_excepts;
    fenv.__mxcsr   &= ~(new_excepts << 7);

    return fesetenv(&fenv) ? -1 : old_excepts;
}

inline int fedisableexcept(unsigned int excepts)
{
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
    // all previous masks
    unsigned int old_excepts;

    if (fegetenv(&fenv)) {
        return -1;
    }
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    // mask
    fenv.__control |= new_excepts;
    fenv.__mxcsr   |= new_excepts << 7;

    return fesetenv(&fenv) ? -1 : old_excepts;
}

#endif

void enable_segfpe() {
  if (!is_fpe_check_mode() || is_legacy_mode()) return;
  #ifdef __GNUC__
    #if defined(__x86_64__)
  // debug_msg("enable segfpe hardware floating point error detection");
       signal(SIGFPE, sighandler);
       feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
    #endif
  #endif
}

void disable_segfpe() {
  if (!is_fpe_check_mode() || is_legacy_mode()) return;
  #ifdef __GNUC__
    #if defined(__x86_64__)
  // debug_msg("disable segfpe");
      fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
    #endif
  #endif
}

void write(const char *s, const char *msg) {
  if (!is_debug_data_mode()) return;
  cout << s << ": " << msg << endl;
}

void write(const double d, const char *msg) {
  if (!is_debug_data_mode()) return;
  cout << std::setprecision(6) << d << ": " << msg << endl;
}

void write(const gsl_vector *v, const char *msg) {
  if (!is_debug_data_mode()) return;
  if (msg) cout << "// " << msg << endl;
  cout << "// vector size: " << v->size << endl;
  cout << "double " << msg << "[] = {";
  for (size_t i=0; i < v->size; i++) {
    cout << gsl_vector_get(v,i) << ",";
  }
  cout << "}" << endl;
}

void write(const gsl_matrix *m, const char *msg) {
  if (!is_debug_data_mode()) return;
  if (msg) cout << "// " << msg << endl;
  // Matrices are stored in row-major order, meaning that each row of
  // elements forms a contiguous block in memory. This is the standard
  // “C-language ordering” of two-dimensional arrays. The number of
  // rows is size1.
  auto rows = m->size1; // see https://www.gnu.org/software/gsl/manual/html_node/Matrices.html#Matrices
  auto cols = m->size2;
  auto tda = m->tda;

  cout << "// matrix size: " << cols << " cols, " << rows << " rows," << tda << " tda" << endl;
  cout << "double " << msg << "[] = {";
  for (size_t row=0; row < rows; row++) {
    for (size_t col=0; col < cols; col++) {
      // cout << "(" << i << "," << j << ")";
      cout << gsl_matrix_safe_get(m,row,col);
      cout << ",";
    }
    cout << "// row " << row << endl;
  }
  cout << "}" << endl;
}

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

void do_gsl_matrix_safe_free (gsl_matrix *m, const char *__pretty_function, const char *__file, int __line, bool warn_only) {
  enforce(m);
  if (is_strict_mode() && is_check_mode() && is_debug_mode()) {
    bool has_NaN = has_nan(m);
    bool has_Inf = has_inf(m);
    if (has_NaN || has_Inf) {
      write(m);
      std::string msg = "Matrix (size ";
      msg += std::to_string(m->size1);
      msg += "x";
      msg += std::to_string(m->size2);
      msg += ")";
      if (warn_only) {
        if (has_Inf)
          warning_at_msg(__file,__line,(msg+" contains Infinite on free!").c_str());
        if (has_NaN)
          warning_at_msg(__file,__line,(msg+" contains NaN on free!").c_str());
      }
      else {
        if (has_Inf)
          warnfail_at_msg(is_strict_mode(),__pretty_function,__file,__line,(msg+" contains Infinite on free!").c_str());
        if (has_NaN)
          warnfail_at_msg(is_strict_mode(),__pretty_function,__file,__line,(msg+" contains NaN on free!").c_str());
      }
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
  if (is_strict_mode() && is_check_mode() && is_debug_mode()) {
    bool has_NaN = has_nan(v);
    bool has_Inf = has_inf(v);
    if (has_NaN || has_Inf) {
      write(v);
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

double do_gsl_matrix_safe_get(const gsl_matrix * m, const size_t row, const size_t col,
                              const char *__pretty_function, const char *__file, int __line) {
  enforce(m);
  if (!is_legacy_mode() && (is_debug_mode() || is_check_mode() || is_strict_mode())) {
    auto rows = m->size1; // see above write function
    auto cols = m->size2;
    if (col >= cols || row >= rows) {
      std::string msg = "Matrix out of bounds (" + std::to_string(rows) + "," + std::to_string(cols) + ") ";
      msg += std::to_string(row);
      msg += "r,";
      msg += std::to_string(col);
      fail_at_msg(__file,__line,msg.c_str());
    }
  }
  return gsl_matrix_get(m,row,col);
}


char *do_strtok_safe(char *tokenize, const char *delimiters, const char *__pretty_function, const char *__file, int __line,
                     const char *infile) {
  auto token = strtok(tokenize,delimiters);
  if (token == NULL) {
    if (infile)
      fail_at_msg(__file,__line,string("Parsing input file '") + infile + "' failed in function " + __pretty_function);
    else
      fail_at_msg(__file,__line,string("Parsing input file failed in function ") + __pretty_function);
  }
  return token;
}

// Helper function called by macro validate_K(K, check). K is validated
// unless -no-check option is used.
void do_validate_K(const gsl_matrix *K, const char *__pretty_function, const char *__file, int __line) {
  if (is_check_mode()) {
    // debug_msg("Validating K");
    auto eigenvalues = getEigenValues(K);
    const uint count_small = count_abs_small_values(eigenvalues,EIGEN_MINVALUE);
    if (count_small>1) {
      std::string msg = "K has ";
      msg += std::to_string(count_small);
      msg += " eigenvalues close to zero";
      warning_at_msg(__file,__line,msg);
    }
    if (isMatrixIllConditioned(eigenvalues))
      warning_at_msg(__file,__line,"K is ill conditioned!");
    if (!isMatrixSymmetric(K))
      warnfail_at_msg(is_strict_mode(),__pretty_function,__file,__line,"K is not symmetric!" );
    const bool negative_eigen_values = has_negative_values_but_one(eigenvalues);
    if (negative_eigen_values) {
      warning_at_msg(__file,__line,"K has more than one negative eigenvalues!");
    }
    if (count_small>1 && negative_eigen_values && !isMatrixPositiveDefinite(K))
      warnfail_at_msg(is_strict_mode(),__pretty_function,__file,__line,"K is not positive definite!");
    gsl_vector_free(eigenvalues);
  }
}
