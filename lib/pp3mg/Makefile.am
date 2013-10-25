# Remove files created during configure.
DISTCLEANFILES = fcs-package.info

# Subdirectories to descend into.
SUBDIRS = pmg pp3mg

sublibs = pmg/libpmg.la pp3mg/libpp3mg.la

# Let aclocal find the add-on macros.
ACLOCAL_AMFLAGS = -I m4

if ENABLE_LIBRARY_INSTALL
lib_LTLIBRARIES = libfcs_pp3mg.la
else
noinst_LTLIBRARIES = libfcs_pp3mg.la
endif

libfcs_pp3mg_la_SOURCES =

libfcs_pp3mg_la_LIBADD = $(sublibs)
