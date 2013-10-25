#
# Copyright (C) 2010 by Argonne National Laboratory.
#     See COPYRIGHT in top-level directory.
#

liba1_la_SOURCES += src/a1d/pthreadd/pthreadd_initialize.c \
	src/a1d/pthreadd/pthreadd_finalize.c \
	src/a1d/pthreadd/pthreadd_collectives.c \
	src/a1d/pthreadd/pthreadd_misc.c 

noinst_HEADERS += src/a1d/pthreadd/pthreaddimpl.h
