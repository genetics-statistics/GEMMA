#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <assert.h>
#include <iostream>

// enforce works like assert but also when NDEBUG is set (i.e., it
// always works). enforce_msg prints message instead of expr

#define ROUND(f) round(f * 10000.)/10000
#if defined NDEBUG
#define debug_msg(msg)
#define assert_issue(is_issue, expr)
#else
#define debug_msg(msg) cout << "**** DEBUG: " << msg << endl;
#define assert_issue(is_issue, expr) \
  ((is_issue) ? enforce_msg(expr,"FAIL: ISSUE assert") : __ASSERT_VOID_CAST(0))

#endif

/* This prints an "Assertion failed" message and aborts.  */
inline void __enforce_fail(const char *__assertion, const char *__file,
                    unsigned int __line,
                    const char *__function)
{
  std::cout << "ERROR: Enforce failed for " << __assertion << " in " << __file << " at line " << __line << " in " << __function << std::endl;
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
