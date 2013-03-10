# Subdirectories
SUBDIRS = 

# Enable Fortran subdirectory only on request.
if ENABLE_FORTRAN
SUBDIRS += f03
endif

# Directory of ipnfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel 

# Directory of pnfft.h
AM_CPPFLAGS += -I$(top_srcdir)/api

# Libraries to add to all programs that are built.
LDADD = $(top_builddir)/lib@PNFFT_PREFIX@pnfft@PREC_SUFFIX@.la $(pfft_LIBS) $(fftw3_mpi_LIBS) $(fftw3_LIBS)

# These programs are built by 'make check' and may be tested afterwards.
check_PROGRAMS = \
	simple_test

#if ENABLE_COMMON_PNFFT
# do not compile checks
#else
check_PROGRAMS += \
	pnfft_check \
	pnfft_check_vs_pfft \
	pnfft_check_grad \
	pnfft_check_2d \
	pnfft_test \
	pnfft_test_adv
#endif

EXTRA_DIST = \
	plot_m_vs_error.py

