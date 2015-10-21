#
# Copyright (C) 2010. See COPYRIGHT in top-level directory.
#

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/src

lib_LTLIBRARIES = libarmci.la

# Needed to connect with the GA build system
noinst_LTLIBRARIES = libarmcii.la

libarmci_la_SOURCES = src/buffer.c        \
                      src/debug.c         \
                      src/groups.c        \
                      src/internals.c     \
                      src/malloc.c        \
                      src/gmr.c           \
                      src/message.c       \
                      src/message_gop.c   \
                      src/mutex.c         \
                      src/mutex_hdl_queue.c \
                      src/onesided.c      \
                      src/onesided_nb.c   \
                      src/rmw.c           \
                      src/strided.c       \
                      src/topology.c      \
                      src/util.c          \
                      src/value_ops.c     \
                      src/vector.c        \
                      src/init_finalize.c \
                      src/conflict_tree.c \
                      src/parmci.c

libarmci_la_LDFLAGS = -version-info $(libarmci_abi_version)

libarmcii_la_SOURCES = $(libarmci_la_SOURCES)
libarmcii_la_LDFLAGS = $(libarmci_abi_version)

include_HEADERS = src/armci.h src/message.h src/armcix.h

bin_PROGRAMS =
check_PROGRAMS = 
TESTS = 
XFAIL_TESTS = 

MPIEXEC = mpiexec -n 2
TESTS_ENVIRONMENT = $(MPIEXEC)

include benchmarks/Makefile.mk
include tests/Makefile.mk

.PHONY: checkprogs
checkprogs: $(check_PROGRAMS)
