#
# Copyright (C) 2010 by Argonne National Laboratory.
#     See COPYRIGHT in top-level directory.
#

liba1_la_SOURCES += src/a1d/mpid/mpid_param.c \
	src/a1d/mpid/mpid_initialize.c \
	src/a1d/mpid/mpid_finalize.c \
	src/a1d/mpid/mpid_malloc.c \
	src/a1d/mpid/mpid_free.c \
	src/a1d/mpid/mpid_flush.c \
	src/a1d/mpid/mpid_flush_all.c \
	src/a1d/mpid/mpid_put.c \
	src/a1d/mpid/mpid_get.c \
	src/a1d/mpid/mpid_collectives.c \
	src/a1d/mpid/mpid_misc.c 

noinst_HEADERS += src/a1d/mpid/mpidimpl.h
