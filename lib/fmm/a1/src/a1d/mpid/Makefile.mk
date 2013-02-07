#
# Copyright (C) 2010 by Argonne National Laboratory.
#     See COPYRIGHT in top-level directory.
#

liba1_la_SOURCES += $(top_srcdir)/src/a1d/mpid/mpid_param.c \
    $(top_srcdir)/src/a1d/mpid/mpid_initialize.c \
	$(top_srcdir)/src/a1d/mpid/mpid_finalize.c \
	$(top_srcdir)/src/a1d/mpid/mpid_malloc.c \
	$(top_srcdir)/src/a1d/mpid/mpid_free.c \
	$(top_srcdir)/src/a1d/mpid/mpid_flush.c \
	$(top_srcdir)/src/a1d/mpid/mpid_flush_all.c \
	$(top_srcdir)/src/a1d/mpid/mpid_put.c \
	$(top_srcdir)/src/a1d/mpid/mpid_get.c \
	$(top_srcdir)/src/a1d/mpid/mpid_collectives.c \
	$(top_srcdir)/src/a1d/mpid/mpid_misc.c 

noinst_HEADERS += $(top_srcdir)/src/a1d/mpid/mpidimpl.h
