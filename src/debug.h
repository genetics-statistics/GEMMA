#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <assert.h>
#include <iostream>

#include "gsl/gsl_matrix.h"

void gemma_gsl_error_handler (const char * reason,
                              const char * file,
                              int line, int gsl_errno);


// Validation routines
void do_validate_K(const gsl_matrix *K, bool do_check, bool strict, const char *__file, int __line);

#define ROUND(f) round(f * 10000.)/10000
#define validate_K(K,check,strict) do_validate_K(K,check,strict,__FILE__,__LINE__)

#define warning_at_msg(__file,__line,msg) cerr << "**** WARNING: " << msg << " in " << __file << " at line " << __line << endl;

inline void fail_at_msg(bool strict, const char *__file, int __line, const char *msg) {
  if (strict)
    std::cerr << "**** STRICT FAIL: ";
  else
    std::cerr << "**** WARNING: ";
  std::cerr << msg << " in " << __file << " at line " << __line << std::endl;
  if (strict)
    exit(1);
}
#if defined NDEBUG

#define warning_msg(msg) cerr << "**** WARNING: " << msg << endl;
#define debug_msg(msg)
#define assert_issue(is_issue, expr)

#else // DEBUG

#define warning_msg(msg) cerr << "**** WARNING: " << msg << " in " << __FILE__ << " at line " << __LINE__ << " in " << __FUNCTION__ << endl;
#define debug_msg(msg) cerr << "**** DEBUG: " << msg << " in " << __FILE__ << " at line " << __LINE__ << " in " << __FUNCTION__ << endl;
#define assert_issue(is_issue, expr) \
  ((is_issue) ? enforce_msg(expr,"FAIL: ISSUE assert") : __ASSERT_VOID_CAST(0))

#endif

// enforce works like assert but also when NDEBUG is set (i.e., it
// always works). enforce_msg prints message instead of expr

/* This prints an "Assertion failed" message and aborts.  */
inline void __enforce_fail(const char *__assertion, const char *__file,
                    unsigned int __line,
                    const char *__function)
{
  std::cout << "ERROR: Enforce failed for " << __assertion << " in " << __file << " at line " << __line << " in " << __PRETTY_FUNCTION__ << std::endl;
  exit(1);
}

#define __ASSERT_FUNCTION __PRETTY_FUNCTION__

#define enforce(expr)                                                          \
  ((expr)                                                                      \
       ? __ASSERT_VOID_CAST(0)                                                 \
       : __enforce_fail(__STRING(expr), __FILE__, __LINE__, __ASSERT_FUNCTION))

#define enforce_msg(expr, msg)                                                 \
  ((expr) ? __ASSERT_VOID_CAST(0)                                              \
          : __enforce_fail(msg, __FILE__, __LINE__, __ASSERT_FUNCTION))

#define enforce_str(expr, msg)                                                 \
  ((expr)                                                                      \
       ? __ASSERT_VOID_CAST(0)                                                 \
       : __enforce_fail((msg).c_str(), __FILE__, __LINE__, __ASSERT_FUNCTION))

// Helpers to create a unique varname per MACRO
#define COMBINE1(X, Y) X##Y
#define COMBINE(X, Y) COMBINE1(X, Y)

#define enforce_gsl(expr)                                                      \
  auto COMBINE(res, __LINE__) = (expr);                                        \
  (COMBINE(res, __LINE__) == 0                                                 \
       ? __ASSERT_VOID_CAST(0)                                                 \
       : __enforce_fail(gsl_strerror(COMBINE(res, __LINE__)), __FILE__,         \
                        __LINE__, __ASSERT_FUNCTION))

#endif
