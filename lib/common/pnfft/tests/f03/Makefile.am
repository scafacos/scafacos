# Get Fortran compile rules that include preprocessing.
include $(top_srcdir)/build-aux/fortran-rules.am

# Directory of pnfft.f (which is build first)
AM_FCCPPFLAGS = -I$(top_builddir)/api

# Libraries to add to all programs that are built.
LDADD = $(top_builddir)/lib@PNFFT_PREFIX@pnfft@PREC_SUFFIX@.la $(pfft_LIBS) $(fftw3_mpi_LIBS) $(fftw3_LIBS)

# noinst_LTLIBRARIES = libtests.la

# search for Fortran source files to build objects
AM_DEFAULT_SOURCE_EXT = .f90

check_PROGRAMS = 

# These programs are built by 'make check' and may be tested afterwards.
#if ENABLE_COMMON_PNFFT
# do not compile checks
#else
# check_PROGRAMS += \
# 	pnfft_test
#endif
