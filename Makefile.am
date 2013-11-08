###################################################
# Makefile for the ScaFaCos fcs interface         #
#                                                 #
###################################################

# Traverse these subdirectories, the current one first (for config.h).
SUBDIRS = . lib src package test

if ENABLE_DOC
SUBDIRS += doc
endif

DIST_SUBDIRS = $(SUBDIRS)

# Let aclocal find macro files in the m4 directory.
ACLOCAL_AMFLAGS = -I m4

# Generated public headers that get installed.
nodist_include_HEADERS = fcs_config.h fcs_fconfig.h

# These files are not automatically distributed.
EXTRA_DIST = COPYING.GPL COPYING.LGPL
EXTRA_DIST += build-aux/fortran-depcomp

EXTRA_DIST += examples/Makefile examples/fmm_test.F90

# Remove files created during configure.
DISTCLEANFILES = fcs-package.info

#################################################################
# Documentation
#################################################################
if ENABLE_DOC
.PHONY: FORCE doc manual doc-upload
doc manual doc-upload: FORCE
	cd doc; $(MAKE) --print-directory $@
FORCE:
endif
