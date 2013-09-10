#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.60])
AC_INIT([P2NFFT], [0.1], [sebastian.banert@s2008.tu-chemnitz.de])

AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([*      Configuring solver P2NFFT                               *])
AC_MSG_NOTICE([****************************************************************])

AC_CONFIG_SRCDIR([init.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([1.7 -Wall foreign])

AC_LANG([C])

# Checks for the MPI C compiler
AX_PROG_CC_MPI(,,AC_MSG_FAILURE([The solver P2NNFT requires an MPI C compiler.]))

AC_PROG_RANLIB
AM_PROG_CC_C_O
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# Checks for header files.
AC_CHECK_HEADERS([gsl/gsl_sf_bessel.h gsl/gsl_sf_gamma.h],[],[AC_MSG_ERROR([Required header files for GNU Scientific Library not found.])])

# Set up FCS library arguments.
AX_FCS_CONFIG_SOLVER_ARGS

# Set up FCS library.
AX_FCS_CONFIG_SOLVER

# Checks for structures.

# Checks for compiler characteristics.
AC_PROG_CC_STDC

# Checks for library functions.
p2nfft_ldadd=
AC_SEARCH_LIBS([cos],[m],[test "${ac_cv_search_cos}" != "none required" && p2nfft_ldadd="${ac_cv_search_cos} ${p2nfft_ldadd}"],[AC_MSG_ERROR([Required library for math functions not found.])])
AC_SEARCH_LIBS([cblas_dgemm],[gslcblas],[test "${ac_cv_search_cblas_dgemm}" != "none required" && p2nfft_ldadd="${ac_cv_search_cblas_dgemm} ${p2nfft_ldadd}"],[AC_MSG_ERROR([Required library for GNU Scientific Library not found.])])
AC_SEARCH_LIBS([gsl_blas_dgemm],[gsl],[test "${ac_cv_search_gsl_blas_dgemm}" != "none required" && p2nfft_ldadd="${ac_cv_search_gsl_blas_dgemm} ${p2nfft_ldadd}"],[AC_MSG_ERROR([Required library for GNU Scientific Library not found.])])

AX_LIB_FFTW3

if test "x$ax_lib_fftw3_mpi" = "xno"; then
  AC_MSG_ERROR([You don't seem to have the FFTW3-MPI library installed. You can ]
      [download it from http://www.fftw.org. If you have installed FFTW3-MPI, ]
      [make sure that this configure script can find it. See ./configure ]
      [--help for more information.])
fi

# threads are not yet implemented
#if  test "x$enable_threads" = "xyes" -a "x$ax_lib_fftw3_threads" = "xno"; then
#  AC_MSG_ERROR([You don't seem to have the threaded FFTW 3 library installed.])
#fi

AX_LIB_PFFT

AX_LIB_PNFFT

LIBS="$LIBS $fftw3_mpi_LIBS $fftw3_LIBS"
CPPFLAGS="$CPPFLAGS $fftw3_mpi_CPPFLAGS $fftw3_CPPFLAGS"
LDFLAGS="$LDFLAGS $fftw3_mpi_LDFLAGS $fftw3_LDFLAGS"

# Checks for system services

# we save the cache here because the tests below fail more commonly

# Add P2NFFT solver stuff to main FCS package.
AX_FCS_PACKAGE_RESET
AX_FCS_PACKAGE_ADD([p2nfft_LIBS],[-lfcs_p2nfft])
AX_FCS_PACKAGE_ADD([p2nfft_LIBS_A],[libfcs_p2nfft.a])
AX_FCS_PACKAGE_ADD([p2nfft_LDADD],[${p2nfft_ldadd}])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
