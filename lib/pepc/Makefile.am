# Let aclocal find the add-on macros.
ACLOCAL_AMFLAGS = -I m4

# Distribute script for Fortran module dependency generation.
EXTRA_DIST = build-aux/fortran-depcomp

# Remove files created during configure.
DISTCLEANFILES = fcs-package.info

# Subdirectories to descend into.
SUBDIRS = src sl_pepc openpa

sublibs = src/liblpepc.la sl_pepc/libsl_pepc.la openpa/src/libopa.la

if ENABLE_LIBRARY_INSTALL
lib_LTLIBRARIES = libfcs_pepc.la
else
noinst_LTLIBRARIES = libfcs_pepc.la
endif

libfcs_pepc_la_SOURCES =

libfcs_pepc_la_LIBADD = $(sublibs)
