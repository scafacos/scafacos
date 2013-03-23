
# Get linker flags from scafacos-fcs.mk. (The assignment and evaluation is
# required to force Automake to include 'scafacos-fcs.mk' early enough.)
DUMMY:=$(eval include $(top_builddir)/package/scafacos-fcs.mk)

# Preprocessor flags.
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/lib -I$(top_srcdir)/lib/common/fcs-common

# Libraries.
TEST_C_LDADD = $(SCAFACOS_MK_LDADD)

# Programs to build.
check_PROGRAMS = interface_test
interface_test_DEPENDENCIES = $(SCAFACOS_MK_DEPS)

# Programs to build.
check_PROGRAMS += multiple_test
multiple_test_DEPENDENCIES = $(SCAFACOS_MK_DEPS)

if ENABLE_DIRECT
check_PROGRAMS += test_direct
test_direct_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
endif

if ENABLE_MMM1D
check_PROGRAMS += test_mmm1d
test_mmm1d_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
endif

if ENABLE_MMM2D
check_PROGRAMS += test_mmm2d
test_mmm2d_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
endif

if ENABLE_P3M
check_PROGRAMS += test_p3m
test_p3m_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
endif

if ENABLE_P2NFFT
check_PROGRAMS += test_p2nfft
test_p2nfft_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
check_PROGRAMS += test_p2nfft_compute_near
test_p2nfft_compute_near_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
if ENABLE_DIRECT
check_PROGRAMS += test_p2nfft_solver_only
test_p2nfft_solver_only_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
check_PROGRAMS += test_p2nfft_nonperiodic
test_p2nfft_nonperiodic_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
check_PROGRAMS += test_p2nfft_nonperiodic_minimal
test_p2nfft_nonperiodic_minimal_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
endif
endif

if ENABLE_VMG
check_PROGRAMS += test_vmg_periodic test_vmg_open
test_vmg_periodic_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
test_vmg_periodic_LDADD = $(TEST_C_LDADD)
test_vmg_open_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
test_vmg_open_LDADD = $(TEST_C_LDADD)
endif

if ENABLE_PEPC
check_PROGRAMS += test_pepc
test_pepc_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
endif

#if ENABLE_FMM
#check_PROGRAMS += test_fmm
#test_fmm_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
#endif

# Quick fix for test_p2nfft_solver_only
if ENABLE_P2NFFT
if ENABLE_DIRECT
AM_CPPFLAGS += -I$(top_srcdir)/lib/p2nfft
endif
endif

# Scripts (distributed) needed for tests.
dist_check_SCRIPTS =
if ENABLE_DIRECT
dist_check_SCRIPTS += start_direct.sh
endif
if ENABLE_FMM
dist_check_SCRIPTS += start_fmm.sh
endif
if ENABLE_MMM1D
dist_check_SCRIPTS += start_mmm1d.sh
endif
if ENABLE_P2NFFT
dist_check_SCRIPTS += start_p2nfft.sh
if ENABLE_DIRECT
dist_check_SCRIPTS += start_p2nfft_nonperiodic.sh
endif
endif
if ENABLE_P3M
dist_check_SCRIPTS += start_p3m.sh
endif
if ENABLE_PP3MG_PMG
dist_check_SCRIPTS += start_pp3mg.sh
endif
if ENABLE_VMG
dist_check_SCRIPTS += start_vmg.sh
endif
if ENABLE_PEPC
dist_check_SCRIPTS += start_pepc.sh
endif

#if ENABLE_FMM
#dist_check_SCRIPTS += start_fmm.sh
#endif

check_PROGRAMS += test_comm
test_comm_DEPENDENCIES = $(SCAFACOS_MK_DEPS)
dist_check_SCRIPTS += start_comm.sh

LDADD = $(TEST_C_LDADD)

# Tests to run with 'make check'.
TESTS = $(dist_check_SCRIPTS)

# Clean up after ourselves.
CLEANFILES = \
	load_*.dat \
	timing_*.dat \
	parts_all.dat \
	domains.dat \
	run.out \
	*.file \
	C_test.xyz

clean-local:
	-rm -rf diag
