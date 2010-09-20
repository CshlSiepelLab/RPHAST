AC_PREREQ(2.57)
AC_INIT([rphast], 1.1, [mjhubisz@cornell.edu], [], [http://compgen.bscb.cornell.edu/phast])
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

OLD_CFLAGS=${CFLAGS}
CFLAGS="${CFLAGS} -fno-strict-aliasing"
AC_MSG_CHECKING(for -fno-strict-aliasing)
AC_COMPILE_IFELSE(AC_LANG_PROGRAM(), AC_MSG_RESULT(yes), [AC_MSG_RESULT(no)
		                           CFLAGS=${OLD_CFLAGS}])

AC_C_INLINE

##################################################
## The following comes from pcre's configure.ac ##
## (edited for RPHAST purposes)                 ##
##################################################

m4_define(pcre_major, [8])
m4_define(pcre_minor, [02])
m4_define(pcre_prerelease, [])
m4_define(pcre_date, [2010-03-19])

# Libtool shared library interface versions (current:revision:age)
m4_define(libpcre_version, [0:1:0])
m4_define(libpcreposix_version, [0:0:0])
m4_define(libpcrecpp_version, [0:0:0])

AC_CONFIG_HEADERS(src/pcre/config.h)

# This was added at the suggestion of libtoolize (03-Jan-10)
AC_CONFIG_MACRO_DIR([src/pcre/m4])

# Check for a 64-bit integer type
AC_TYPE_INT64_T

# Handle --enable-rebuild-chartables
AC_ARG_ENABLE(rebuild-chartables,
              AS_HELP_STRING([--enable-rebuild-chartables],
                             [rebuild character tables in current locale]),
              , enable_rebuild_chartables=no)

# Handle --enable-utf8 (disabled by default)
AC_ARG_ENABLE(utf8,
              AS_HELP_STRING([--enable-utf8],
                             [enable UTF-8 support (incompatible with --enable-ebcdic)]),
              , enable_utf8=unset)

# Handle --enable-unicode-properties
AC_ARG_ENABLE(unicode-properties,
              AS_HELP_STRING([--enable-unicode-properties],
                             [enable Unicode properties support (implies --enable-utf8)]),
              , enable_unicode_properties=no)

# Separate newline options
ac_pcre_newline=lf
AC_ARG_ENABLE(newline-is-cr,
              AS_HELP_STRING([--enable-newline-is-cr],
                             [use CR as newline character]),
              ac_pcre_newline=cr)
AC_ARG_ENABLE(newline-is-lf,
              AS_HELP_STRING([--enable-newline-is-lf],
                             [use LF as newline character (default)]),
              ac_pcre_newline=lf)
AC_ARG_ENABLE(newline-is-crlf,
              AS_HELP_STRING([--enable-newline-is-crlf],
                             [use CRLF as newline sequence]),
              ac_pcre_newline=crlf)
AC_ARG_ENABLE(newline-is-anycrlf,
              AS_HELP_STRING([--enable-newline-is-anycrlf],
                             [use CR, LF, or CRLF as newline sequence]),
              ac_pcre_newline=anycrlf)
AC_ARG_ENABLE(newline-is-any,
              AS_HELP_STRING([--enable-newline-is-any],
                             [use any valid Unicode newline sequence]),
              ac_pcre_newline=any)
enable_newline="$ac_pcre_newline"

# Handle --enable-bsr-anycrlf
AC_ARG_ENABLE(bsr-anycrlf,
              AS_HELP_STRING([--enable-bsr-anycrlf],
                             [\R matches only CR, LF, CRLF by default]),
              , enable_bsr_anycrlf=no)

# Handle --enable-ebcdic
AC_ARG_ENABLE(ebcdic,
              AS_HELP_STRING([--enable-ebcdic],
                             [assume EBCDIC coding rather than ASCII; incompatible with --enable-utf8; use only in (uncommon) EBCDIC environments; it implies --enable-rebuild-chartables]),
              , enable_ebcdic=no)

# Handle --disable-stack-for-recursion
AC_ARG_ENABLE(stack-for-recursion,
              AS_HELP_STRING([--disable-stack-for-recursion],
                             [don't use stack recursion when matching]),
              , enable_stack_for_recursion=yes)

# Handle --enable-pcregrep-libz
AC_ARG_ENABLE(pcregrep-libz,
              AS_HELP_STRING([--enable-pcregrep-libz],
                             [link pcregrep with libz to handle .gz files]),
              , enable_pcregrep_libz=no)

# Handle --enable-pcregrep-libbz2
AC_ARG_ENABLE(pcregrep-libbz2,
              AS_HELP_STRING([--enable-pcregrep-libbz2],
                             [link pcregrep with libbz2 to handle .bz2 files]),
              , enable_pcregrep_libbz2=no)

# Handle --enable-pcretest-libreadline
AC_ARG_ENABLE(pcretest-libreadline,
              AS_HELP_STRING([--enable-pcretest-libreadline],
                             [link pcretest with libreadline]),
              , enable_pcretest_libreadline=no)

# Handle --with-posix-malloc-threshold=NBYTES
AC_ARG_WITH(posix-malloc-threshold,
            AS_HELP_STRING([--with-posix-malloc-threshold=NBYTES],
                           [threshold for POSIX malloc usage (default=10)]),
            , with_posix_malloc_threshold=10)

# Handle --with-link-size=N
AC_ARG_WITH(link-size,
            AS_HELP_STRING([--with-link-size=N],
                           [internal link size (2, 3, or 4 allowed; default=2)]),
            , with_link_size=2)

# Handle --with-match-limit=N
AC_ARG_WITH(match-limit,
            AS_HELP_STRING([--with-match-limit=N],
                           [default limit on internal looping (default=10000000)]),
            , with_match_limit=10000000)

# Handle --with-match-limit_recursion=N
#
# Note: In config.h, the default is to define MATCH_LIMIT_RECURSION
# symbolically as MATCH_LIMIT, which in turn is defined to be some numeric
# value (e.g. 10000000). MATCH_LIMIT_RECURSION can otherwise be set to some
# different numeric value (or even the same numeric value as MATCH_LIMIT,
# though no longer defined in terms of the latter).
#
AC_ARG_WITH(match-limit-recursion,
            AS_HELP_STRING([--with-match-limit-recursion=N],
                           [default limit on internal recursion (default=MATCH_LIMIT)]),
            , with_match_limit_recursion=MATCH_LIMIT)

# Make sure that if enable_unicode_properties was set, that UTF-8 support
# is enabled.
#
if test "x$enable_unicode_properties" = "xyes"
then
  if test "x$enable_utf8" = "xno"
  then
    AC_MSG_ERROR([support for Unicode properties requires UTF-8 support])
  fi
  enable_utf8=yes
fi

if test "x$enable_utf8" = "xunset"
then
  enable_utf8=no
fi

# Make sure that if enable_ebcdic is set, rebuild_chartables is also enabled.
# Also check that UTF-8 support is not requested, because PCRE cannot handle
# EBCDIC and UTF-8 in the same build. To do so it would need to use different
# character constants depending on the mode.
#
if test "x$enable_ebcdic" = "xyes"
then
  enable_rebuild_chartables=yes
  if test "x$enable_utf8" = "xyes"
  then
    AC_MSG_ERROR([support for EBCDIC and UTF-8 cannot be enabled at the same time])
  fi
fi

# Convert the newline identifier into the appropriate integer value.
case "$enable_newline" in
  lf)      ac_pcre_newline_value=10   ;;
  cr)      ac_pcre_newline_value=13   ;;
  crlf)    ac_pcre_newline_value=3338 ;;
  anycrlf) ac_pcre_newline_value=-2   ;;
  any)     ac_pcre_newline_value=-1   ;;
  *)
  AC_MSG_ERROR([invalid argument \"$enable_newline\" to --enable-newline option])
  ;;
esac

# Check argument to --with-link-size
case "$with_link_size" in
  2|3|4) ;;
  *)
  AC_MSG_ERROR([invalid argument \"$with_link_size\" to --with-link-size option])
  ;;
esac

AH_TOP([
/* On Unix-like systems config.h.in is converted by "configure" into config.h.
Some other environments also support the use of "configure". PCRE is written in
Standard C, but there are a few non-standard things it can cope with, allowing
it to run on SunOS4 and other "close to standard" systems.

If you are going to build PCRE "by hand" on a system without "configure" you
should copy the distributed config.h.generic to config.h, and then set up the
macro definitions the way you need them. You must then add -DHAVE_CONFIG_H to
all of your compile commands, so that config.h is included at the start of
every source.

Alternatively, you can avoid editing by using -D on the compiler command line
to set the macro values. In this case, you do not have to set -DHAVE_CONFIG_H.

PCRE uses memmove() if HAVE_MEMMOVE is set to 1; otherwise it uses bcopy() if
HAVE_BCOPY is set to 1. If your system has neither bcopy() nor memmove(), set
them both to 0; an emulation function will be used. */])


# The files below are C++ header files.
pcre_have_type_traits="0"
pcre_have_bits_type_traits="0"

# Using AC_SUBST eliminates the need to include config.h in a public .h file
AC_SUBST(pcre_have_type_traits)
AC_SUBST(pcre_have_bits_type_traits)

# Conditional compilation
AM_CONDITIONAL(WITH_REBUILD_CHARTABLES, test "x$enable_rebuild_chartables" = "xyes")

# Checks for typedefs, structures, and compiler characteristics.

AC_C_CONST

# Checks for library functions.

AC_CHECK_FUNCS(bcopy memmove strerror)

# Check for the availability of libz (aka zlib)

AC_CHECK_HEADERS([zlib.h], [HAVE_ZLIB_H=1])
AC_CHECK_LIB([z], [gzopen], [HAVE_LIBZ=1])

# Check for the availability of libbz2. Originally we just used AC_CHECK_LIB,
# as for libz. However, this had the following problem, diagnosed and fixed by
# a user:
#
#   - libbz2 uses the Pascal calling convention (WINAPI) for the functions
#     under Win32.
#   - The standard autoconf AC_CHECK_LIB fails to include "bzlib.h",
#     therefore missing the function definition.
#   - The compiler thus generates a "C" signature for the test function.
#   - The linker fails to find the "C" function.
#   - PCRE fails to configure if asked to do so against libbz2.
#
# Solution:
#
#   - Replace the AC_CHECK_LIB test with a custom test.

AC_CHECK_HEADERS([bzlib.h], [HAVE_BZLIB_H=1])
# Original test
# AC_CHECK_LIB([bz2], [BZ2_bzopen], [HAVE_LIBBZ2=1])
#
# Custom test follows

AC_MSG_CHECKING([for libbz2])
OLD_LIBS="$LIBS"
LIBS="$LIBS -lbz2"
AC_LINK_IFELSE( AC_LANG_PROGRAM([[
#ifdef HAVE_BZLIB_H
#include <bzlib.h>
#endif]],
[[return (int)BZ2_bzopen("conftest", "rb");]]),
[AC_MSG_RESULT([yes]);HAVE_LIBBZ2=1; break;],
AC_MSG_RESULT([no]))
LIBS="$OLD_LIBS"

# Check for the availabiity of libreadline

AC_CHECK_HEADERS([readline/readline.h], [HAVE_READLINE_H=1])
AC_CHECK_HEADERS([readline/history.h], [HAVE_HISTORY_H=1])
AC_CHECK_LIB([readline], [readline], [HAVE_LIB_READLINE=1])

# This facilitates -ansi builds under Linux
dnl AC_DEFINE([_GNU_SOURCE], [], [Enable GNU extensions in glibc])

PCRE_STATIC_CFLAG=""
if test "x$enable_shared" = "xno" ; then
  AC_DEFINE([PCRE_STATIC], [1], [
    Define if linking statically (TODO: make nice with Libtool)])
  PCRE_STATIC_CFLAG="-DPCRE_STATIC"
fi
AC_SUBST(PCRE_STATIC_CFLAG)

# Here is where pcre specific defines are handled

if test "$enable_utf8" = "yes"; then
  AC_DEFINE([SUPPORT_UTF8], [], [
    Define to enable support for the UTF-8 Unicode encoding. This will
    work even in an EBCDIC environment, but it is incompatible with
    the EBCDIC macro. That is, PCRE can support *either* EBCDIC code
    *or* ASCII/UTF-8, but not both at once.])
fi

if test "$enable_unicode_properties" = "yes"; then
  AC_DEFINE([SUPPORT_UCP], [], [
    Define to enable support for Unicode properties])
fi

if test "$enable_stack_for_recursion" = "no"; then
  AC_DEFINE([NO_RECURSE], [], [
    PCRE uses recursive function calls to handle backtracking while
    matching. This can sometimes be a problem on systems that have
    stacks of limited size. Define NO_RECURSE to get a version that
    doesn't use recursion in the match() function; instead it creates
    its own stack by steam using pcre_recurse_malloc() to obtain memory
    from the heap. For more detail, see the comments and other stuff
    just above the match() function. On systems that support it,
    "configure" can be used to set this in the Makefile
    (use --disable-stack-for-recursion).])
fi

if test "$enable_pcregrep_libz" = "yes"; then
  AC_DEFINE([SUPPORT_LIBZ], [], [
    Define to allow pcregrep to be linked with libz, so that it is
    able to handle .gz files.])
fi

if test "$enable_pcregrep_libbz2" = "yes"; then
  AC_DEFINE([SUPPORT_LIBBZ2], [], [
    Define to allow pcregrep to be linked with libbz2, so that it is
    able to handle .bz2 files.])
fi

if test "$enable_pcretest_libreadline" = "yes"; then
  AC_DEFINE([SUPPORT_LIBREADLINE], [], [
    Define to allow pcretest to be linked with libreadline.])
fi

AC_DEFINE_UNQUOTED([NEWLINE], [$ac_pcre_newline_value], [
  The value of NEWLINE determines the newline character sequence. On
  systems that support it, "configure" can be used to override the
  default, which is 10. The possible values are 10 (LF), 13 (CR),
  3338 (CRLF), -1 (ANY), or -2 (ANYCRLF).])

if test "$enable_bsr_anycrlf" = "yes"; then
  AC_DEFINE([BSR_ANYCRLF], [], [
    By default, the \R escape sequence matches any Unicode line ending
    character or sequence of characters. If BSR_ANYCRLF is defined, this is
    changed so that backslash-R matches only CR, LF, or CRLF. The build-
    time default can be overridden by the user of PCRE at runtime. On
    systems that support it, "configure" can be used to override the
    default.])
fi

AC_DEFINE_UNQUOTED([LINK_SIZE], [$with_link_size], [
  The value of LINK_SIZE determines the number of bytes used to store
  links as offsets within the compiled regex. The default is 2, which
  allows for compiled patterns up to 64K long. This covers the vast
  majority of cases. However, PCRE can also be compiled to use 3 or 4
  bytes instead. This allows for longer patterns in extreme cases. On
  systems that support it, "configure" can be used to override this default.])

AC_DEFINE_UNQUOTED([POSIX_MALLOC_THRESHOLD], [$with_posix_malloc_threshold], [
  When calling PCRE via the POSIX interface, additional working storage
  is required for holding the pointers to capturing substrings because
  PCRE requires three integers per substring, whereas the POSIX
  interface provides only two. If the number of expected substrings is
  small, the wrapper function uses space on the stack, because this is
  faster than using malloc() for each call. The threshold above which
  the stack is no longer used is defined by POSIX_MALLOC_THRESHOLD. On
  systems that support it, "configure" can be used to override this
  default.])

AC_DEFINE_UNQUOTED([MATCH_LIMIT], [$with_match_limit], [
  The value of MATCH_LIMIT determines the default number of times the
  internal match() function can be called during a single execution of
  pcre_exec(). There is a runtime interface for setting a different
  limit. The limit exists in order to catch runaway regular
  expressions that take for ever to determine that they do not match.
  The default is set very large so that it does not accidentally catch
  legitimate cases. On systems that support it, "configure" can be
  used to override this default default.])

AC_DEFINE_UNQUOTED([MATCH_LIMIT_RECURSION], [$with_match_limit_recursion], [
  The above limit applies to all calls of match(), whether or not they
  increase the recursion depth. In some environments it is desirable
  to limit the depth of recursive calls of match() more strictly, in
  order to restrict the maximum amount of stack (or heap, if
  NO_RECURSE is defined) that is used. The value of
  MATCH_LIMIT_RECURSION applies only to recursive calls of match(). To
  have any useful effect, it must be less than the value of
  MATCH_LIMIT. The default is to use the same value as MATCH_LIMIT.
  There is a runtime method for setting a different limit. On systems
  that support it, "configure" can be used to override the default.])

AC_DEFINE([MAX_NAME_SIZE], [32], [
  This limit is parameterized just in case anybody ever wants to
  change it. Care must be taken if it is increased, because it guards
  against integer overflow caused by enormously large patterns.])

AC_DEFINE([MAX_NAME_COUNT], [10000], [
  This limit is parameterized just in case anybody ever wants to
  change it. Care must be taken if it is increased, because it guards
  against integer overflow caused by enormously large patterns.])

AH_VERBATIM([PCRE_EXP_DEFN], [
/* If you are compiling for a system other than a Unix-like system or
   Win32, and it needs some magic to be inserted before the definition
   of a function that is exported by the library, define this macro to
   contain the relevant magic. If you do not define this macro, it
   defaults to "extern" for a C compiler and "extern C" for a C++
   compiler on non-Win32 systems. This macro apears at the start of
   every exported function that is part of the external API. It does
   not appear on functions that are "external" in the C sense, but
   which are internal to the library. */
#undef PCRE_EXP_DEFN])

if test "$enable_ebcdic" = "yes"; then
  AC_DEFINE_UNQUOTED([EBCDIC], [], [
    If you are compiling for a system that uses EBCDIC instead of ASCII
    character codes, define this macro as 1. On systems that can use
    "configure", this can be done via --enable-ebcdic. PCRE will then
    assume that all input strings are in EBCDIC. If you do not define
    this macro, PCRE will assume input strings are ASCII or UTF-8 Unicode.
    It is not possible to build a version of PCRE that supports both
    EBCDIC and UTF-8.])
fi

# Platform specific issues
NO_UNDEFINED=
EXPORT_ALL_SYMBOLS=
case $host_os in
  cygwin* | mingw* )
    if test X"$enable_shared" = Xyes; then
      NO_UNDEFINED="-no-undefined"
      EXPORT_ALL_SYMBOLS="-Wl,--export-all-symbols"
    fi
    ;;
esac

# The extra LDFLAGS for each particular library
# (Note: The libpcre*_version bits are m4 variables, assigned above)

EXTRA_LIBPCRE_LDFLAGS="$EXTRA_LIBPCRE_LDFLAGS \
                       $NO_UNDEFINED -version-info libpcre_version"

EXTRA_LIBPCREPOSIX_LDFLAGS="$EXTRA_LIBPCREPOSIX_LDFLAGS \
                            $NO_UNDEFINED -version-info libpcreposix_version"

EXTRA_LIBPCRECPP_LDFLAGS="$EXTRA_LIBPCRECPP_LDFLAGS \
                          $NO_UNDEFINED -version-info libpcrecpp_version \
                          $EXPORT_ALL_SYMBOLS"

AC_SUBST(EXTRA_LIBPCRE_LDFLAGS)
AC_SUBST(EXTRA_LIBPCREPOSIX_LDFLAGS)
AC_SUBST(EXTRA_LIBPCRECPP_LDFLAGS)

# When we run 'make distcheck', use these arguments.
DISTCHECK_CONFIGURE_FLAGS="--enable-cpp --enable-unicode-properties"
AC_SUBST(DISTCHECK_CONFIGURE_FLAGS)

# Check that, if --enable-pcregrep-libz or --enable-pcregrep-libbz2 is
# specified, the relevant library is available.

if test "$enable_pcregrep_libz" = "yes"; then
  if test "$HAVE_ZLIB_H" != "1"; then
    echo "** Cannot --enable-pcregrep-libz because zlib.h was not found"
    exit 1
  fi
  if test "$HAVE_LIBZ" != "1"; then
    echo "** Cannot --enable-pcregrep-libz because libz was not found"
    exit 1
  fi
  LIBZ="-lz"
fi
AC_SUBST(LIBZ)

if test "$enable_pcregrep_libbz2" = "yes"; then
  if test "$HAVE_BZLIB_H" != "1"; then
    echo "** Cannot --enable-pcregrep-libbz2 because bzlib.h was not found"
    exit 1
  fi
  if test "$HAVE_LIBBZ2" != "1"; then
    echo "** Cannot --enable-pcregrep-libbz2 because libbz2 was not found"
    exit 1
  fi
  LIBBZ2="-lbz2"
fi
AC_SUBST(LIBBZ2)

# Similarly for --enable-pcretest-readline

if test "$enable_pcretest_libreadline" = "yes"; then
  if test "$HAVE_READLINE_H" != "1"; then
    echo "** Cannot --enable-pcretest-readline because readline/readline.h was not found."
    exit 1
  fi
  if test "$HAVE_HISTORY_H" != "1"; then
    echo "** Cannot --enable-pcretest-readline because readline/history.h was not found."
    exit 1
  fi
  LIBREADLINE="-lreadline"
fi
AC_SUBST(LIBREADLINE)

AC_CONFIG_FILES(src/pcre/pcre.h)

AC_OUTPUT()

dnl end configure.ac

##################################################
## End pcre configure.ac section                ##
##################################################


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

