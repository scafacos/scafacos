#
# Copyright (C) 2010. See COPYRIGHT in top-level directory.
#

check_PROGRAMS += tests/contrib/non-blocking/overlap tests/contrib/non-blocking/simple

TESTS          += tests/contrib/non-blocking/simple

tests_contrib_non_blocking_overlap_SOURCES = tests/contrib/non-blocking/overlap.c
tests_contrib_non_blocking_overlap_LDADD = libarmci.la -lm

tests_contrib_non_blocking_simple_SOURCES = tests/contrib/non-blocking/simple.c
tests_contrib_non_blocking_simple_LDADD = libarmci.la -lm
