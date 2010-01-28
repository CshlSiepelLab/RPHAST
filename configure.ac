AC_INIT([rphast], 0.0.1, [mjhubisz@cornell.edu])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   AC_MSG_ERROR([Could not determine R_HOME.])
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
RCFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CFLAGS="${RCFLAGS} ${CFLAGS}"

# use -fgnu89-inline if necessary
OLDCFLAGS="${CFLAGS}"
CFLAGS="${CFLAGS} -fgnu89-inline"
AC_TRY_LINK(, [int i; i++;], INLINE_FLAG="-fgnu89-inline")
CFLAGS="${OLD_CFLAGS}"
PHAST_VERSION="\"`cat src/version`\""

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
echo "LIBS=${LIBS}"
echo "PTHREAD_LIBS=${PTHREAD_LIBS}"

AC_C_INLINE

AC_SUBST(PHAST_VERSION)
AC_SUBST(INLINE_FLAG)
AC_SUBST(PTHREAD_LIBS)

AC_OUTPUT(src/Makevars)
