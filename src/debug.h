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

#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <assert.h>
#include <iostream>
#include <csignal>

#include "gsl/gsl_matrix.h"

void gemma_gsl_error_handler (const char * reason,
                              const char * file,
                              int line, int gsl_errno);

void debug_set_debug_mode(bool setting);
void debug_set_debug_data_mode(bool setting);
void debug_set_check_mode(bool setting);
void debug_set_no_check_mode(bool setting);
void debug_set_no_fpe_check_mode(bool setting);
void debug_set_strict_mode(bool setting);
void debug_set_quiet_mode(bool setting);
void debug_set_issue(uint issue);
void debug_set_legacy_mode(bool setting);

bool is_debug_mode();
bool is_debug_data_mode();
bool is_no_check_mode();
bool is_check_mode();
bool is_fpe_check_mode();
bool is_strict_mode();
bool is_quiet_mode();
bool is_issue(uint issue);
bool is_legacy_mode();

void enable_segfpe();
void disable_segfpe();

#define check_int_mult_overflow(m,n) \
  { auto x = m * n;                                      \
    enforce_msg(x / m == n, "multiply integer overflow"); }

void write(const double d, const char *msg = "");
void write(const char *s, const char *msg = "");
void write(const gsl_vector *v, const char *msg = "");
void write(const gsl_matrix *m, const char *msg = "");

gsl_matrix *gsl_matrix_safe_alloc(size_t rows,size_t cols);
int gsl_matrix_safe_memcpy (gsl_matrix *dest, const gsl_matrix *src);
void gsl_matrix_safe_free (gsl_matrix *v);
void gsl_matrix_warn_free (gsl_matrix *v);
void do_gsl_matrix_safe_free (gsl_matrix *m, const char *__pretty_function, const char *__file, int __line, bool warn_only);

double do_gsl_matrix_safe_get (const gsl_matrix * m, const size_t i, const size_t j, const char *__pretty_function, const char *__file, int __line);
#define gsl_matrix_safe_get(m,i,j) do_gsl_matrix_safe_get(m, i, j,__SHOW_FUNC,__FILE__,__LINE__);

gsl_vector *gsl_vector_safe_alloc(size_t n);
int gsl_vector_safe_memcpy (gsl_vector *dest, const gsl_vector *src);
void gsl_vector_safe_free (gsl_vector *v);
void do_gsl_vector_safe_free (gsl_vector *v, const char *__pretty_function, const char *__file, int __line);

char *do_strtok_safe(char *tokenize, const char *delimiters, const char *__pretty_function, const char *__file, int __line, const char *infile = NULL);
#define strtok_safe2(string,delimiters,infile) do_strtok_safe(string,delimiters,__SHOW_FUNC,__FILE__,__LINE__,infile)
#define strtok_safe(string,delimiters) do_strtok_safe(string,delimiters,__SHOW_FUNC,__FILE__,__LINE__)

// Validation routines
void do_validate_K(const gsl_matrix *K, const char*__pretty_func, const char *__file, int __line);

#define ROUND(f) round(f * 10000.)/10000
#define validate_K(K) do_validate_K(K,__SHOW_FUNC,__FILE__,__LINE__)

#define warning_at_msg(__file,__line,msg) cerr << "**** WARNING: " << msg << " in " << __file << " at line " << __line << endl;

inline void warnfail_at_msg(bool strict, const char *__function, const char *__file, int __line, const char *msg) {
  if (strict)
    std::cerr << "**** STRICT FAIL: ";
  else
    std::cerr << "**** WARNING: ";
  std::cerr << msg << " in " << __file << " at line " << __line << " in " << __function << std::endl;
  if (strict)
    std::raise(SIGINT); // keep stack trace for gdb
}

inline void fail_at_msg(const char *__file, int __line, std::string msg) {
  std::cerr << "**** FAILED: " << msg << " in " << __file << " at line " << __line << std::endl;
  std::raise(SIGINT); // keep stack trace for gdb
}

# ifndef __ASSERT_VOID_CAST
# define __ASSERT_VOID_CAST (void)
# endif

inline void fail_msg(const char *msg) {
  std::cerr << "**** FAILED: " << msg << std::endl;
  std::raise(SIGINT); // keep stack trace for gdb
}

inline void fail_msg(std::string msg) {
  std::cerr << "**** FAILED: " << msg << std::endl;
  std::raise(SIGINT); // keep stack trace for gdb
}

#define info_msg(msg) if (!is_quiet_mode()) cerr << "**** INFO: " << msg << "." << endl;
#define msg(msg) info_msg(msg);

#if defined NDEBUG
  #define __SHOW_FUNC __func__

  #define warning_msg(msg) cerr << "**** WARNING: " << msg << endl;
  #define debug_msg(msg)
  #define assert_issue(is_issue, expr)

#else // DEBUG

  #define __SHOW_FUNC __func__

  #define warning_msg(msg) cerr << "**** WARNING: " << msg << " in " << __FILE__ << " at line " << __LINE__ << " in " << __func__ << endl;
  #define debug_msg(msg) (is_debug_mode() && cerr << "**** DEBUG: " << msg << " in " << __FILE__ << " at line " << __LINE__ << " in " << __func__ << endl);
  #define assert_issue(is_issue, expr) \
    ((is_issue) ? enforce_msg(expr,"FAIL: ISSUE assert") : __ASSERT_VOID_CAST(0))

#endif // NDEBUG

// enforce works like assert but also when NDEBUG is set (i.e., it
// always works). enforce_msg prints message instead of expr

/* This prints an "Assertion failed" message and aborts.  */
inline void __enforce_fail(const char *__assertion, const char *__file,
                    unsigned int __line,
                    const char *__function)
{
  std::cout << "ERROR: Enforce failed for " << __assertion << " in " << __file << " at line " << __line << " in " << __function << std::endl;
  std::raise(SIGINT); // keep stack trace for gdb
  // exit(1);
}

#define enforce(expr)                                                          \
  ((expr)                                                                      \
       ? __ASSERT_VOID_CAST(0)                                                 \
       : __enforce_fail(__STRING(expr), __FILE__, __LINE__, __SHOW_FUNC))

#define enforce_msg(expr, msg)                                                 \
  ((expr) ? __ASSERT_VOID_CAST(0)                                              \
          : __enforce_fail(msg, __FILE__, __LINE__, __SHOW_FUNC))

#define enforce_str(expr, msg)                                                 \
  ((expr)                                                                      \
       ? __ASSERT_VOID_CAST(0)                                                 \
       : __enforce_fail((msg).c_str(), __FILE__, __LINE__, __SHOW_FUNC))

#define enforce_is_int(s) \
  enforce_str(std::regex_match(s, std::regex("^[-+]?[0-9]+$")),s + " not an integer")

#define enforce_is_float(s)                                             \
  enforce_str(std::regex_match(s, std::regex("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$")),s + " not a float")

// Helpers to create a unique varname per MACRO
#define COMBINE1(X, Y) X##Y
#define COMBINE(X, Y) COMBINE1(X, Y)

#define enforce_gsl(expr)                                                      \
  auto COMBINE(res, __LINE__) = (expr);                                        \
  (COMBINE(res, __LINE__) == 0                                                 \
       ? __ASSERT_VOID_CAST(0)                                                 \
       : __enforce_fail(gsl_strerror(COMBINE(res, __LINE__)), __FILE__,         \
                        __LINE__, __SHOW_FUNC))

#define enforce_fexists(fn, msg)                                               \
  if (!fn.empty())                                                             \
    enforce_msg(stat(fn.c_str(), &fileInfo) == 0,                              \
                ((std::string(__STRING(fn)) + " " + fn + ": " + msg).c_str()));

#define gsl_matrix_safe_free(m) \
  do_gsl_matrix_safe_free(m,__SHOW_FUNC,__FILE__,__LINE__,false);
#define gsl_matrix_warn_free(m) \
  do_gsl_matrix_safe_free(m,__SHOW_FUNC,__FILE__,__LINE__,true);
#define gsl_vector_safe_free(v) \
  do_gsl_vector_safe_free(v,__SHOW_FUNC,__FILE__,__LINE__);

#endif
