AC_INIT([rphast], 0.0.1, [mjhubisz@cornell.edu], [], [http://compgen.bscb.cornell.edu/phast])
AC_CONFIG_SRCDIR([RPHAST/src/rph_init.c])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   AC_MSG_ERROR([Could not determine R_HOME.])
fi

AC_PROG_CC
CC=`"${R_HOME}/bin/R" CMD config CC`
RCFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CFLAGS="${RCFLAGS} ${CFLAGS}"

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([stddef.h stdlib.h string.h strings.h sys/time.h unistd.h])


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


# use -pthread if available otherwise -lpthread
AC_MSG_CHECKING(for pthread)

OLD_LIBS=${LIBS}
LIBS=""
# first check if we need -pthread at all
AC_LINK_IFELSE([
#include <pthread.h>
int main() {
  pthread_t x;
  pthread_join(x, NULL);
  return 0;
}
], AC_MSG_RESULT(yes), 

# now see if we can use -pthread, which is preferred over -lpthread?
[ LIBS=-pthread; AC_LINK_IFELSE([
#include <pthread.h>
int main() {
  pthread_t x;
  pthread_join(x, NULL);
  return 0;
}
], AC_MSG_RESULT(${LIBS}), [ LIBS=-lpthread; AC_LINK_IFELSE([
#include <pthread.h>
int main() {
  pthread_t x;
  pthread_join(x, NULL);
  return 0;
]}, AC_MSG_RESULT(${LIBS}), [AC_MSG_RESULT(no)
                            AC_MSG_ERROR(pthread libs required)])])])

PTHREAD_LIBS=${LIBS}
LIBS=${OLD_LIBS}

# use vecLib if available
AC_MSG_CHECKING(for vecLib)
OLD_CFLAGS=${CFLAGS}
OLD_LIBS=${LIBS}
CFLAGS="${CFLAGS} -DVECLIB"
LIBS="${LIBS} -framework vecLib"
AC_LINK_IFELSE( AC_LANG_PROGRAM(), 
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
  AC_LINK_IFELSE(AC_LANG_PROGRAM(),
                 [AC_MSG_RESULT(yes)
                  LAPACK_LIBS="\${LAPACK_LIBS} \${BLAS_LIBS} \${FLIBS}"
		  LAPACK_FLAG="-DR_LAPACK"
                 ],
                 [AC_MSG_RESULT(no)
                  AC_MSG_ERROR(Error linking to R's lapack library.  Check R installation.)])
fi


LIBS="${OLD_LIBS}"

PHAST_VERSION="\"`cat src/version`\""

AC_C_INLINE

AC_SUBST(PHAST_VERSION)
AC_SUBST(INLINE_FLAG)
AC_SUBST(PTHREAD_LIBS)
AC_SUBST(LAPACK_FLAG)
AC_SUBST(LAPACK_LIBS)

AC_OUTPUT(src/Makevars)