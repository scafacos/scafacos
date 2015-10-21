#
# Copyright (C) 2010. See COPYRIGHT in top-level directory.
#

check_PROGRAMS += tests/contrib/transp1D/transp1D-c

TESTS          += tests/contrib/transp1D/transp1D-c

tests_contrib_transp1D_transp1D_c_SOURCES = tests/contrib/transp1D/transp1D-c.c
tests_contrib_transp1D_transp1D_c_LDADD = libarmci.la
