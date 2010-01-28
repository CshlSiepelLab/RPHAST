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
OLD_LIBS=${LIBS}
LIBS="-pthread"
AC_TRY_LINK(, [int i; i++;], PTHREAD_LIBS=${LIBS}, oops="failed")
if test -n "${oops}" ; then
  LIBS="-lpthread"
  AC_TRY_LINK(, [int i; i++;], PTHREAD_LIBS=${LIBS}, oops="again")
  if test "${oops}" = "again" ; then
   echo "pthread not found"
  fi
fi
LIBS="${OLD_LIBS}"


# use -fgnu89-inline if necessary
OLDCFLAGS="${CFLAGS}"
CFLAGS="${CFLAGS} -fgnu89-inline"
AC_TRY_COMPILE(, [int i; i++;], INLINE_FLAG="-fgnu89-inline")
CFLAGS="${OLD_CFLAGS}"
PHAST_VERSION="\"`cat src/version`\""


AC_C_INLINE

AC_SUBST(PHAST_VERSION)
AC_SUBST(INLINE_FLAG)
AC_SUBST(PTHREAD_LIBS)

AC_OUTPUT(src/Makevars)
