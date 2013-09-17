#
# Copyright (C) 2010 by Argonne National Laboratory.
#     See COPYRIGHT in top-level directory.
#

liba1_la_SOURCES += src/a1d/dcmfd/dcmfd_param.c \
	src/a1d/dcmfd/dcmfd_initialize.c \
	src/a1d/dcmfd/dcmfd_finalize.c \
	src/a1d/dcmfd/dcmfd_malloc.c \
	src/a1d/dcmfd/dcmfd_free.c \
	src/a1d/dcmfd/dcmfd_flush.c \
	src/a1d/dcmfd/dcmfd_flush_group.c \
	src/a1d/dcmfd/dcmfd_put.c \
	src/a1d/dcmfd/dcmfd_puts.c \
	src/a1d/dcmfd/dcmfd_putv.c \
	src/a1d/dcmfd/dcmfd_get.c \
	src/a1d/dcmfd/dcmfd_gets.c \
	src/a1d/dcmfd/dcmfd_getv.c \
	src/a1d/dcmfd/dcmfd_putacc.c \
	src/a1d/dcmfd/dcmfd_putaccs.c \
	src/a1d/dcmfd/dcmfd_putaccv.c \
	src/a1d/dcmfd/dcmfd_putmodv.c \
	src/a1d/dcmfd/dcmfd_rmw.c \
	src/a1d/dcmfd/dcmfd_collectives.c \
	src/a1d/dcmfd/dcmfd_misc.c \
	src/a1d/dcmfd/dcmfd_util.c \
	src/a1d/dcmfd/dcmfd_requestpool.c \
	src/a1d/dcmfd/dcmfd_handlepool.c \
	src/a1d/dcmfd/dcmfd_bufferpool.c \
	src/a1d/dcmfd/dcmfd_counter.c \
	src/a1d/dcmfd/dcmfd_mutex.c \
	src/a1d/dcmfd/dcmfd_cht.c \
	src/a1d/dcmfd/dcmfd_wait.c

noinst_HEADERS += src/a1d/dcmfd/dcmfdimpl.h
