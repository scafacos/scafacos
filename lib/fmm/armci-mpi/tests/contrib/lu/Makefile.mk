#
# Copyright (C) 2010. See COPYRIGHT in top-level directory.
#

check_PROGRAMS += tests/contrib/lu/lu tests/contrib/lu/lu-block tests/contrib/lu/lu-b-bc

TESTS          += tests/contrib/lu/lu-block tests/contrib/lu/lu-b-bc

tests_contrib_lu_lu_SOURCES = tests/contrib/lu/lu.c tests/contrib/lu/lu_timing.c
tests_contrib_lu_lu_LDADD = libarmci.la -lm

tests_contrib_lu_lu_block_SOURCES = tests/contrib/lu/lu-block.c tests/contrib/lu/lu_timing.c
tests_contrib_lu_lu_block_LDADD = libarmci.la -lm

tests_contrib_lu_lu_b_bc_SOURCES = tests/contrib/lu/lu-b-bc.c tests/contrib/lu/lu_timing.c
tests_contrib_lu_lu_b_bc_LDADD = libarmci.la -lm
