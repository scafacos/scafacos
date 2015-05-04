# Get linker flags from scafacos.mk generated during configure.
# * The 'indirection' of the @...@ substitution allows to inject an include
#   command into the Makefile without being expanded by Automake.
# * The variable definition ensures that Automake adds this line before the
#   build rules.
# * Using the command as a variable name 'activates' the deferred expansion of
#   the SCAFACOS_MK_LDADD and SCAFACOS_MK_DEPS variables from scafacos.mk.
DUMMY_SCAFACOS_MK_@EVAL_INCLUDE_PACKAGE_SCAFACOS_MK@ = DUMMY

# Get Fortran compile rules that include preprocessing.
include $(top_srcdir)/build-aux/fortran-rules.am

# Include paths for preprocessor and Fortran compiler.
#AM_CPPFLAGS = -I$(top_srcdir)/src/fortran -I$(top_builddir)/src/fortran
AM_FCFLAGS = $(FC_MODINC)$(top_srcdir)/src/fortran $(FC_MODINC)$(top_builddir)/src/fortran

# Libraries.
LDADD = $(SCAFACOS_MK_LDADD)

# Programs to build.
check_PROGRAMS = interface_test
#check_PROGRAMS = interface_test
interface_test_SOURCES = interface_test.f90
interface_test_DEPENDENCIES = $(SCAFACOS_MK_DEPS)

#test_pp3mg_SOURCES = test_pp3mg.f90
#test_pepc_SOURCES = test_pepc.f90

#test_pp3mg_LDADD = $(LIBS_pp3mg)
#test_pepc_LDADD =  $(LIBS_pepc)

# Scripts (distributed) needed for tests.
dist_check_SCRIPTS =
if ENABLE_DIRECT
dist_check_SCRIPTS += start_direct.sh
endif
if ENABLE_EWALD
dist_check_SCRIPTS += start_ewald.sh
endif
if ENABLE_FMM
dist_check_SCRIPTS += start_fmm.sh
endif
if ENABLE_MMM1D
dist_check_SCRIPTS += start_mmm1d.sh
endif
#don't know working parameters, so test is disabled
#if ENABLE_MMM2D
#dist_check_SCRIPTS += start_mmm2d.sh
#endif
if ENABLE_PEPC
dist_check_SCRIPTS += start_pepc.sh
endif
if ENABLE_P2NFFT
dist_check_SCRIPTS += start_p2nfft.sh
endif
if ENABLE_P3M
dist_check_SCRIPTS += start_p3m.sh
endif
if ENABLE_VMG
dist_check_SCRIPTS += start_vmg.sh
endif

# Tests to run.
TESTS = $(dist_check_SCRIPTS)

# Clean up after ourselves.
CLEANFILES = \
	load_*.dat \
	timing_*.dat \
	parts_all.dat \
	domains.dat \
	run.out \
	F_test.* \
	*.file
