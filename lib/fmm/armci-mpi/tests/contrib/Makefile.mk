#
# Copyright (C) 2010. See COPYRIGHT in top-level directory.
#

check_PROGRAMS += \
                  tests/contrib/armci-perf         \
                  tests/contrib/armci-test         \
                  # end

TESTS          += \
                  tests/contrib/armci-perf         \
                  tests/contrib/armci-test         \
                  # end

tests_contrib_armci_perf_LDADD = libarmci.la -lm
tests_contrib_armci_test_LDADD = libarmci.la -lm

include tests/contrib/cg/Makefile.mk
include tests/contrib/lu/Makefile.mk
include tests/contrib/transp1D/Makefile.mk
include tests/contrib/non-blocking/Makefile.mk
