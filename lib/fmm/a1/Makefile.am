#
# Copyright (C) 2010 by Argonne National Laboratory.
#     See COPYRIGHT in top-level directory.
#

# Remove files created during configure.
DISTCLEANFILES = fcs-package.info

EXTRA_DIST = maint/version

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/src/include -I$(top_builddir)/src/include

lib_LTLIBRARIES = liba1.la
liba1_la_SOURCES =
liba1_la_LDFLAGS = -version-info $(liba1_abi_version)

noinst_HEADERS = $(top_builddir)/src/include/a1.h $(top_builddir)/src/include/a1d.h $(top_builddir)/src/include/a1u.h $(top_builddir)/src/include/armci.h

include src/adl/util/Makefile.mk
include src/adl/api/Makefile.mk
include src/adl/interfaces/armci/Makefile.mk
#TODO# include src/adl/interfaces/global-arrays/Makefile.mk
include src/adl/interfaces/armci/Makefile.mk

if a1_device_dcmfd
include src/a1d/dcmfd/Makefile.mk
endif

#if a1_device_pthreadd
#include src/a1d/pthreadd/Makefile.mk
#endif

#if a1_device_mpid
#include src/a1d/mpid/Makefile.mk
#endif

#if a1_device_geminid
#include src/a1d/geminid/Makefile.mk
#endif

#if a1_device_lapid
#include src/a1d/lapid/Makefile.mk
#endif

#if a1_device_ibd
#include src/a1d/ibd/Makefile.mk
#endif

docs:
	if test "@DOXYGEN@" != "" ; then @DOXYGEN@ a1-doxygen.cfg ; fi
