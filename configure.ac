AC_PREREQ(2.57)
AC_INIT([rphast], 1.5, [mjhubisz@cornell.edu], [], [http://compgen.cshl.edu/phast])
#AC_CONFIG_SRCDIR([RPHAST/src/rph_init.c])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   AC_MSG_ERROR([Could not determine R_HOME.])
fi

AC_MSG_CHECKING(for CFLAGS)
#AC_PROG_CC
CC=`"${R_HOME}/bin/R" CMD config CC`
#RCFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
RCFLAGS="-g"
#CFLAGS="${RCFLAGS} ${CFLAGS}"
CFLAGS=" -g"
AC_MSG_RESULT(${CFLAGS})
AC_PROG_CC
AC_MSG_RESULT(${CFLAGS})

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([stddef.h stdlib.h string.h strings.h sys/time.h unistd.h], ,
[AC_MSG_ERROR(
[Couldn't find all necessary system headers])])

AC_CHECK_HEADER(pcre.h, , [AC_MSG_ERROR([Couldn't find pcre.h... download/install at www.pcre.org])])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([atexit bzero floor gethostname gettimeofday memset pow re_comp regcomp rint sqrt strcasecmp strdup strstr strtol])

AC_SEARCH_LIBS([pcre_compile], [pcre], [],
 [AC_MSG_RESULT(no)
  AC_MSG_ERROR([ Need to install pcre ])])

# use vecLib if available
AC_MSG_CHECKING(for vecLib)
OLD_CFLAGS=${CFLAGS}
OLD_LIBS=${LIBS}
CFLAGS="${CFLAGS} -DVECLIB"
LIBS="${LIBS} -framework vecLib"
AC_LINK_IFELSE( [AC_LANG_PROGRAM()], 
		[AC_MSG_RESULT(yes)
                 LAPACK_FLAG="-DVECLIB"
                 LAPACK_LIBS="-framework vecLib"],
		[ AC_MSG_RESULT(no)
                  CFLAGS=${OLD_CFLAGS}])
CFLAGS=${OLD_CFLAGS}
#echo "VECLIB=${VECLIB} CFLAGS=${CFLAGS}"

# if not veclib, try to use R's lapack
if test -z "$LAPACK_FLAG"; then
  LAPACK_LIBS=`"${R_HOME}/bin/R" CMD config LAPACK_LIBS`
  BLAS_LIBS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
  FLIBS=`"${R_HOME}/bin/R" CMD config FLIBS`
  SHLIB=`"${R_HOME}/bin/R" CMD config SHLIB_LDFLAGS`
  FPIC=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
  LIBS="${FPIC} ${SHLIB} ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS}"

  AC_MSG_CHECKING(for R's Lapack)
  AC_LINK_IFELSE([AC_LANG_PROGRAM()],
                 [AC_MSG_RESULT(yes)
                  LAPACK_LIBS="\${LAPACK_LIBS} \${BLAS_LIBS} \${FLIBS}"
		  LAPACK_FLAG="-DR_LAPACK"
                 ],
                 [AC_MSG_RESULT(no)
                  AC_MSG_ERROR(Error linking to R's lapack library.  Check R installation.)])
fi

LIBS="${OLD_LIBS}"

OLD_CFLAGS=${CFLAGS}
CFLAGS="${CFLAGS} -fno-strict-aliasing"
AC_MSG_CHECKING(for -fno-strict-aliasing)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM()], AC_MSG_RESULT(yes), [AC_MSG_RESULT(no)
		                           CFLAGS=${OLD_CFLAGS}])

AC_C_INLINE

AC_SUBST(LAPACK_FLAG)
AC_SUBST(LAPACK_LIBS)
AC_SUBST(CFLAGS)

AC_OUTPUT(src/Makevars)

cat <<EOF

$PACKAGE-$VERSION configuration summary:

    C preprocessor .................. : ${CPP}
    C compiler ...................... : ${CC}
    C compiler flags ................ : ${CFLAGS}
    Extra libraries ................. : ${LIBS}
    Lapack flag ..................... : ${LAPACK_FLAG}
    Lapack libs ..................... : ${LAPACK_LIBS}
    BLAS libs ....................... : ${BLAS_LIBS}
    FLIBS ........................... : ${FLIBS}

EOF

dnl end configure.ac

