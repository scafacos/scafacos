# Remove files created during configure.
DISTCLEANFILES = fcs-package.info

ACLOCAL_AMFLAGS = -I m4

if ENABLE_LIBRARY_INSTALL
lib_LIBRARIES = libfcs_p2nfft.a
else
noinst_LIBRARIES = libfcs_p2nfft.a
endif

# Directory of pnfft.h
AM_CFLAGS = -I$(top_srcdir)/../../lib/common/pnfft/api

# Directory of pfft.h
AM_CFLAGS += -I$(top_srcdir)/../../lib/common/pfft/api

# Directory of fcs.h
AM_CFLAGS += -I$(top_srcdir)/../../src

# Directory of fcs_config.h
AM_CFLAGS += -I$(top_builddir)/../..

#LIBS_P2NFFTTEST = ../direct/libdirect.a

# Libraries to add to all programs that are built.

libfcs_p2nfft_a_CPPFLAGS = -I$(top_srcdir)/../../src -I$(top_srcdir)/../../lib -I$(top_srcdir)/../../lib/common/fcs-common
libfcs_p2nfft_a_SOURCES = \
	types.h \
	utils.h \
	constants.h \
	p2nfft.h \
	init.c init.h \
	tune.c tune.h \
	run.c run.h \
	parameters.c parameters.h \
	nearfield.c nearfield.h \
	kernels.c kernels.h \
	regularization.c regularization.h \
	cg_cos_coeff.c cg_cos_coeff.h \
	cg_cos_err.c cg_cos_err.h

#libfcs_p2nfft_a_LIBADD = $(pnfft_LIBS) $(pfft_LIBS) $(fftw3_mpi_LIBS) $(fftw3_LIBS) ../common/near/libnear.a ../direct/libdirect.a

#noinst_PROGRAMS = tests/test_random
#tests_test_random_SOURCES = tests/test_random.c
#tests_test_random_LDADD = libfcs_p2nfft.a \
#                          ../../src/libfcs.a \
#                          ../common/pnfft/.libs/libpnfft.a \
#                          ../common/pfft/.libs/libpfft.a \
#                          ../direct/libdirect.a \
#                          ../common/near/libnear.a

