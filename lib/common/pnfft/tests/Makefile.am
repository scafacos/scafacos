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
	simple_test \
	simple_test_c2r

#if ENABLE_COMMON_PNFFT
# do not compile checks
#else
check_PROGRAMS += \
	check_trafo check_trafo_transposed \
	check_trafo_grad check_trafo_transposed_grad \
	check_trafo_2d check_trafo_transposed_2d \
	check_adj check_adj_transposed \
	check_vs_pfft \
	pnfft_test \
	pnfft_test_adv \
	check_trafo_vs_ndft check_trafo_grad_vs_ndft check_adj_vs_ndft \
	check_trafo_vs_ndft_transposed_2d \
	simple_test_c2r_c2c_compare
#endif

EXTRA_DIST = \
	plot_m_vs_error.py

