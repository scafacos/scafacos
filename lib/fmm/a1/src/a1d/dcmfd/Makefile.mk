#
# Copyright (C) 2010 by Argonne National Laboratory.
#     See COPYRIGHT in top-level directory.
#

liba1_la_SOURCES += $(top_srcdir)/src/a1d/dcmfd/dcmfd_param.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_initialize.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_finalize.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_malloc.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_free.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_flush.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_flush_group.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_put.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_puts.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_putv.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_get.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_gets.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_getv.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_putacc.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_putaccs.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_putaccv.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_putmodv.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_rmw.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_collectives.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_misc.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_util.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_requestpool.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_handlepool.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_bufferpool.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_counter.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_mutex.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_cht.c \
	$(top_srcdir)/src/a1d/dcmfd/dcmfd_wait.c

noinst_HEADERS += $(top_srcdir)/src/a1d/dcmfd/dcmfdimpl.h
