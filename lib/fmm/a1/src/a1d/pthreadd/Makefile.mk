#
# Copyright (C) 2010 by Argonne National Laboratory.
#     See COPYRIGHT in top-level directory.
#

liba1_la_SOURCES += $(top_srcdir)/src/a1d/pthreadd/pthreadd_initialize.c \
	$(top_srcdir)/src/a1d/pthreadd/pthreadd_finalize.c \
	$(top_srcdir)/src/a1d/pthreadd/pthreadd_collectives.c \
	$(top_srcdir)/src/a1d/pthreadd/pthreadd_misc.c 

noinst_HEADERS += $(top_srcdir)/src/a1d/pthreadd/pthreaddimpl.h
