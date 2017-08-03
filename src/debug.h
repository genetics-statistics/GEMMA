#ifndef __DEBUG_H__
#define __DEBUG_H__

// enforce works like assert but also when NDEBUG is set (i.e., it
// always works). enforce_msg prints message instead of expr

#if defined NDEBUG
#define debug_msg(msg)
#else
#define debug_msg(msg) cout << "**** DEBUG: " << msg << endl;
#endif

/* This prints an "Assertion failed" message and aborts.  */
extern "C" void __assert_fail(const char *__assertion, const char *__file,
                              unsigned int __line,
                              const char *__function) __THROW
    __attribute__((__noreturn__));

#define __ASSERT_FUNCTION __PRETTY_FUNCTION__

#define enforce(expr)                                                          \
  ((expr)                                                                      \
       ? __ASSERT_VOID_CAST(0)                                                 \
       : __assert_fail(__STRING(expr), __FILE__, __LINE__, __ASSERT_FUNCTION))

#define enforce_msg(expr, msg)                                                 \
  ((expr) ? __ASSERT_VOID_CAST(0)                                              \
          : __assert_fail(msg, __FILE__, __LINE__, __ASSERT_FUNCTION))

#define enforce_str(expr, msg)                                                 \
  ((expr)                                                                      \
       ? __ASSERT_VOID_CAST(0)                                                 \
       : __assert_fail((msg).c_str(), __FILE__, __LINE__, __ASSERT_FUNCTION))

// Helpers to create a unique varname per MACRO
#define COMBINE1(X, Y) X##Y
#define COMBINE(X, Y) COMBINE1(X, Y)

#define enforce_gsl(expr)                                                      \
  auto COMBINE(res, __LINE__) = (expr);                                        \
  (COMBINE(res, __LINE__) == 0                                                 \
       ? __ASSERT_VOID_CAST(0)                                                 \
       : __assert_fail(gsl_strerror(COMBINE(res, __LINE__)), __FILE__,         \
                       __LINE__, __ASSERT_FUNCTION))

#endif
