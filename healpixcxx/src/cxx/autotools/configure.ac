AC_INIT([Healpix_cxx], [3.0])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
LT_INIT
AC_CONFIG_MACRO_DIR([m4])

AC_DEFUN([AX_CHECK_COMPILE_FLAG],
[AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_[]_AC_LANG_ABBREV[]flags_$4_$1])dnl
AC_CACHE_CHECK([whether _AC_LANG compiler accepts $1], CACHEVAR, [
  ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
  _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $4 $1"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
    [AS_VAR_SET(CACHEVAR,[yes])],
    [AS_VAR_SET(CACHEVAR,[no])])
  _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])
AS_IF([test x"AS_VAR_GET(CACHEVAR)" = xyes],
  [m4_default([$2], :)],
  [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CHECK_COMPILE_FLAGS

AC_PROG_CC_C99
AC_OPENMP
CPPFLAGS="$CPPFLAGS $OPENMP_CFLAGS"
CFLAGS="$CFLAGS $OPENMP_CFLAGS"
AX_CHECK_COMPILE_FLAG([-march=native],[CC="$CC -march=native"])
AX_CHECK_COMPILE_FLAG([-fno-tree-vectorize],[CFLAGS="$CFLAGS -fno-tree-vectorize"])
AX_CHECK_COMPILE_FLAG([-ffast-math],[CFLAGS="$CFLAGS -ffast-math"])
AX_CHECK_COMPILE_FLAG([-fomit-frame-pointer],[CFLAGS="$CFLAGS -fomit-frame-pointer"])

AC_PROG_CXX
AC_LANG_PUSH([C++])
AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"
AX_CHECK_COMPILE_FLAG([-march=native],[CXX="$CXX -march=native"])
AX_CHECK_COMPILE_FLAG([-fno-tree-vectorize],[CXXFLAGS="$CXXFLAGS -fno-tree-vectorize"])
AX_CHECK_COMPILE_FLAG([-ffast-math],[CXXFLAGS="$CXXFLAGS -ffast-math"])
AX_CHECK_COMPILE_FLAG([-fomit-frame-pointer],[CXXFLAGS="$CXXFLAGS -fomit-frame-pointer"])
AC_LANG_POP([C++])

AC_CHECK_LIB([cfitsio],[ffgnrwll],,AC_MSG_ERROR([could not find the cfitsio library]),[-lm])
AC_CHECK_HEADERS([fitsio.h],,AC_MSG_ERROR([could not find the cfitsio header file]),)

AC_PROG_LIBTOOL

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
