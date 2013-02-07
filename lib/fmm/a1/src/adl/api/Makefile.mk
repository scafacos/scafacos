#
# Copyright (C) 2010 by Argonne National Laboratory.
#     See COPYRIGHT in top-level directory.
#

liba1_la_SOURCES += $(top_srcdir)/src/adl/api/adl_initialize.c \
                    $(top_srcdir)/src/adl/api/adl_finalize.c \
                    $(top_srcdir)/src/adl/api/adl_malloc.c \
                    $(top_srcdir)/src/adl/api/adl_free.c \
                    $(top_srcdir)/src/adl/api/adl_flush.c \
                    $(top_srcdir)/src/adl/api/adl_flush_group.c \
                    $(top_srcdir)/src/adl/api/adl_put.c \
                    $(top_srcdir)/src/adl/api/adl_puts.c \
                    $(top_srcdir)/src/adl/api/adl_putv.c \
                    $(top_srcdir)/src/adl/api/adl_get.c \
                    $(top_srcdir)/src/adl/api/adl_gets.c \
                    $(top_srcdir)/src/adl/api/adl_getv.c \
                    $(top_srcdir)/src/adl/api/adl_putacc.c \
                    $(top_srcdir)/src/adl/api/adl_putaccs.c \
                    $(top_srcdir)/src/adl/api/adl_putaccv.c \
                    $(top_srcdir)/src/adl/api/adl_putmodv.c \
                    $(top_srcdir)/src/adl/api/adl_rmw.c \
                    $(top_srcdir)/src/adl/api/adl_wait.c \
                    $(top_srcdir)/src/adl/api/adl_collectives.c \
                    $(top_srcdir)/src/adl/api/adl_misc.c \
                    $(top_srcdir)/src/adl/api/adl_counter.c \
                    $(top_srcdir)/src/adl/api/adl_mutex.c \
                    $(top_srcdir)/src/adl/api/adl_handle.c
