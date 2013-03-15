AC_INIT(config/config.auto.in)

AC_CHECK_PROG([uname_found],[uname],[1],[0])
if test $uname_found -eq 0 ; then
    echo "No uname found; setting system type to unknown."
    system="unknown"
else
    system=`uname -s`-`uname -r`
fi
AC_LANG([C])
AC_PROG_CC_C99

AC_TRY_COMPILE([], [@%:@ifndef __INTEL_COMPILER
choke me
@%:@endif], [ICC=[yes]], [ICC=[no]])

if test $ICC = yes; then GCC=no; fi
CCTYPE=unknown
if test $GCC = yes; then CCTYPE=gcc; fi
if test $ICC = yes; then CCTYPE=icc; fi
AC_OPENMP

AC_LANG_PUSH([C++])

AC_TRY_COMPILE([], [@%:@ifndef __INTEL_COMPILER
choke me
@%:@endif], [ICPC=[yes]], [ICPC=[no]])

if test $ICPC = yes; then GXX=no; fi
CXXTYPE=unknown
if test $GXX = yes; then CXXTYPE=g++; fi
if test $ICPC = yes; then CXXTYPE=icpc; fi

AC_OPENMP

AC_LANG_POP([C++])

AC_ARG_WITH(libcfitsio,
  [  --with-libcfitsio=<path>         prefix to the cfitsio installation],
  [CPPFLAGS="-I$withval/include $CPPFLAGS"
   LDFLAGS="-L$withval/lib $LDFLAGS"])
AC_ARG_WITH(libcfitsio-include,
  [  --with-libcfitsio-include=<path> path containing the cfitsio headers],
  [CPPFLAGS="-I$withval $CPPFLAGS"])
AC_ARG_WITH(libcfitsio-lib,
  [  --with-libcfitsio-lib=<path>     path containing the cfitsio library],
  [   LDFLAGS="-L$withval $LDFLAGS"])

SILENT_RULE=".SILENT:"
AC_ARG_ENABLE(noisy-make,
  [  --enable-noisy-make     enable detailed make output],
  [if test "$enableval" = yes; then
     SILENT_RULE=""
   fi])

ENABLE_DEBUG=no
AC_ARG_ENABLE(debug,
  [  --enable-debug          enable generation of debugging symbols],
  [if test "$enableval" = yes; then
     ENABLE_DEBUG=yes
   fi])

case $CCTYPE in
  gcc)
    CCFLAGS="-O3 -ffast-math -fomit-frame-pointer -fno-tree-vectorize"
    GCCVERSION="`$CC -dumpversion 2>&1`"
    echo "Using gcc version $GCCVERSION"
    AC_SUBST(GCCVERSION)
    case $system in
      Darwin-*)
        ;;
      *)
        CCFLAGS="$CCFLAGS -ffunction-sections -fdata-sections"
        ;;
    esac
    changequote(,)
    gcc43=`echo $GCCVERSION | grep -c '4\.[3456789]'`
    changequote([,])
    if test $gcc43 -gt 0; then
      CCFLAGS="$CCFLAGS -march=native"
    fi
    ;;
  icc)
    CCFLAGS="-O3 -ip -Wbrief -Wall -vec-report0 -openmp-report0 -wd383,981,1419,1572"
    ;;
  *)
    CCFLAGS="-O2"
    # Don't do anything now
    ;;
esac

case $CXXTYPE in
  g++)
    CXXCFLAGS="-O2 -ffast-math -fomit-frame-pointer"
    GXXVERSION="`$CXX -dumpversion 2>&1`"
    echo "Using g++ version $GXXVERSION"
    AC_SUBST(GXXVERSION)
    case $system in
      Darwin-*)
        ;;
      *)
        CXXCFLAGS="$CXXCFLAGS -ffunction-sections -fdata-sections"
        LDFLAGS="$LDFLAGS -Wl,--gc-sections"
        ;;
    esac
    changequote(,)
    gxx43=`echo $GXXVERSION | grep -c '4\.[3456789]'`
    changequote([,])
    if test $gxx43 -gt 0; then
      CXXCFLAGS="$CXXCFLAGS -march=native"
    fi
    ;;
  icpc)
    CXXCFLAGS="-O3 -ip -Wbrief -Wall -vec-report0 -openmp-report0 -wd383,444,981,1572"
    ;;
  *)
    CXXCFLAGS="-O2"
    ;;
esac

case $system in
  Darwin-*)
    ARCREATE="libtool -static -o"
    ;;
  *)
    ARCREATE="ar cr"
    ;;
esac

CCFLAGS="$CCFLAGS $OPENMP_CFLAGS"
CXXCFLAGS="$CXXCFLAGS $OPENMP_CXXFLAGS"

if test $ENABLE_DEBUG = yes; then
  CCFLAGS="$CCFLAGS -g"
  CXXCFLAGS="$CXXCFLAGS -g"
fi

CCFLAGS_NO_C="$CCFLAGS $CPPFLAGS"
CXXCFLAGS_NO_C="$CXXCFLAGS $CPPFLAGS"

LDCCFLAGS="$LDFLAGS $CCFLAGS"
LDCXXFLAGS="$LDFLAGS $CXXCFLAGS"

AC_CHECK_LIB([cfitsio],[ffgnrwll],,AC_MSG_ERROR([could not find the cfitsio library]),[-lm])
AC_CHECK_HEADERS([fitsio.h],,AC_MSG_ERROR([could not find the cfitsio header file]),)

AC_SUBST(SILENT_RULE)
AC_SUBST(CC)
AC_SUBST(CCFLAGS_NO_C)
AC_SUBST(CXX)
AC_SUBST(CXXCFLAGS_NO_C)
AC_SUBST(LDCCFLAGS)
AC_SUBST(LDCXXFLAGS)
AC_SUBST(LIBS)
AC_SUBST(ARCREATE)

AC_OUTPUT(config/config.auto)
