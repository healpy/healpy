AC_INIT([healpix_cxx], [3.0])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

dnl
dnl Needed for linking on Windows.
dnl Protect with m4_ifdef because AM_PROG_AR is required in
dnl autoconf >= 1.12 when using -Wall, but the macro is
dnl absent in old versions of autoconf.
dnl
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

LT_INIT
AC_CONFIG_MACRO_DIR([m4])

dnl
dnl By default, install the headers into a subdirectory of
dnl ${prefix}/include to avoid possible header filename collisions.
dnl
includedir="${includedir}/${PACKAGE_NAME}"

dnl
dnl Enable silent build rules if this version of Automake supports them
dnl
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

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

dnl
dnl Introduce --enable-native-optimizations command line argument to turn on
dnl -march=native compiler flag, disabled by default.
dnl
AC_ARG_ENABLE(
    [native-optimizations],
    [AS_HELP_STRING([--enable-native-optimizations], [Enable non-portable optimizations for your own CPU by compiling with -march=native @<:@default=no@:>@])]
)

AC_PROG_CC_C99
AC_OPENMP
CPPFLAGS="$CPPFLAGS $OPENMP_CFLAGS"
CFLAGS="$CFLAGS $OPENMP_CFLAGS"
AS_IF(
    [test "x$enable_native_optimizations" = "xyes"],
    [AX_CHECK_COMPILE_FLAG([-march=native],[CC="$CC -march=native"])]
)
AX_CHECK_COMPILE_FLAG([-fno-tree-vectorize],[CFLAGS="$CFLAGS -fno-tree-vectorize"])
AX_CHECK_COMPILE_FLAG([-ffast-math],[CFLAGS="$CFLAGS -ffast-math"])
AX_CHECK_COMPILE_FLAG([-fomit-frame-pointer],[CFLAGS="$CFLAGS -fomit-frame-pointer"])

AC_PROG_CXX
AC_LANG_PUSH([C++])
AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"
AS_IF(
    [test "x$enable_native_optimizations" = "xyes"],
    [AX_CHECK_COMPILE_FLAG([-march=native],[CXX="$CXX -march=native"])]
)
AX_CHECK_COMPILE_FLAG([-fno-tree-vectorize],[CXXFLAGS="$CXXFLAGS -fno-tree-vectorize"])
AX_CHECK_COMPILE_FLAG([-ffast-math],[CXXFLAGS="$CXXFLAGS -ffast-math"])
AX_CHECK_COMPILE_FLAG([-fomit-frame-pointer],[CXXFLAGS="$CXXFLAGS -fomit-frame-pointer"])
AC_LANG_POP([C++])

dnl
dnl Determine the compiler flags (CFITSIO_CFLAGS)
dnl and linker flags (CFITSIO_LIBS) required to build against libcfitsio.
dnl
dnl Add help messages to describe the variables CFITSIO_CFLAGS and CFITSIO_LIBS,
dnl which can be provided as environment variables or as command-line arguments.
dnl These are output variables of PKG_CHECK_MODULES, but they can be overridden.
dnl
AC_ARG_VAR([CFITSIO_CFLAGS], [compiler flags needed for building with cfitsio, if not found by pkg-config])
AC_ARG_VAR([CFITSIO_LIBS], [linker flags needed for building with cfitsio, if not found by pkg-config])
dnl
dnl Now check for cfitsio using pkg-config. If cfitsio is not found with
dnl pkg-config, CFITSIO_LIBS defaults to "-lcfitsio".
dnl
PKG_CHECK_MODULES([CFITSIO], [cfitsio], , [CFITSIO_LIBS=-lcfitsio])
dnl
dnl Linking against libcfitsio requires libm too.
dnl
CFITSIO_LIBS="$CFITSIO_LIBS -lm"
dnl
dnl With pkg-config < 0.24, PKG_CHECK_MODULES does not AC_SUBST the _CFLAGS
dnl and _LIBS variables.
dnl
AC_SUBST([CFITSIO_CFLAGS])
AC_SUBST([CFITSIO_LIBS])
dnl
dnl Perform a test compile and link to determine if cfitsio is usable. If
dnl pkg-config did not find cfitsio, then the user may still have provided
dnl CFITSIO_CFLAGS and CFITSIO_LIBS manually. If pkg-config did find cfitsio,
dnl it does not hurt to test that cfitsio is usable.
dnl
save_CFLAGS=$CFLAGS
save_LIBS=$LIBS
CFLAGS="$CFLAGS $CFITSIO_CFLAGS"
LIBS="$LIBS $CFITSIO_LIBS"
AC_CHECK_LIB([cfitsio], [ffgnrwll], [:], AC_MSG_ERROR([could not find the cfitsio library]))
AC_CHECK_HEADERS([fitsio.h], [:], AC_MSG_ERROR([could not find the cfitsio header file]))
CFLAGS=$save_CFLAGS
LIBS=$save_LIBS

AC_PROG_LIBTOOL

dnl
dnl Create pkgconfig .pc file.
dnl
AX_CREATE_PKGCONFIG_INFO(,,,,[$OPENMP_CFLAGS $CFITSIO_CFLAGS])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
