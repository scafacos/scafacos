# Get Fortran compile rules that include preprocessing.
include $(top_srcdir)/build-aux/fortran-rules.am

# Directory of ipnfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel 

# Directory of pnfft.h
AM_CPPFLAGS += -I$(top_srcdir)/api

# Directory of pnfft.f (which is build first)
AM_CPPFLAGS += -I$(top_builddir)/api

# Same flags for the Fortran preprocessor.
AM_FCCPPFLAGS = $(AM_CPPFLAGS)

# noinst_LTLIBRARIES = libtests.la

# These programs are built by 'make check' and may be tested afterwards.
check_PROGRAMS = simple_test \
		 pnfft_check \
		 pnfft_check_vs_pfft \
		 pnfft_check_grad \
		 pnfft_check_2d \
		 pnfft_test \
		 pnfft_test_adv

# check_PROGRAMS = fortran_ghost_test fortran_simple_test sample_ghost minimal_test simple_test simple_test_4d

# fortran_ghost_test_SOURCES = fortran_ghost_test.F90
# fortran_simple_test_SOURCES = fortran_simple_test.F90
# sample_ghost_SOURCES = sample_ghost.c
# minimal_test_SOURCES = minimal_test.c
simple_test_SOURCES = simple_test.c
# simple_test_4d_SOURCES = simple_test_4d.c

pnfft_check_SOURCES = pnfft_check.c
pnfft_check_vs_pfft_SOURCES = pnfft_check_vs_pfft.c
pnfft_check_grad_SOURCES = pnfft_check_grad.c
pnfft_check_2d_SOURCES = pnfft_check_2d.c
pnfft_test_SOURCES = pnfft_test.c
pnfft_test_adv_SOURCES = pnfft_test_adv.c

# Libraries to add to all programs that are built.
LDADD = $(top_builddir)/lib@PNFFT_PREFIX@pnfft@PREC_SUFFIX@.la $(pfft_LIBS) $(fftw3_mpi_LIBS) $(fftw3_LIBS)

# Group local sources into convenience lib.
# libtests_la_SOURCES = \
#       fortran_ghost_test.F90 \
# 	fortran_simple_test.F90 \
# 	minimal_test.c \
# 	sample_ghost.c \
# 	simple_test.c \
# 	simple_test_4d.c


EXTRA_DIST = \
	plot_m_vs_error.py

