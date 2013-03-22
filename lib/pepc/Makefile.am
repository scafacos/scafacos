# Let aclocal find the add-on macros.
ACLOCAL_AMFLAGS = -I m4

# Import central add_lib_objects support.
include $(top_srcdir)/build-aux/lib_a_support.am

# Distribute script for Fortran module dependency generation.
EXTRA_DIST = build-aux/fortran-depcomp

# Remove files created during configure.
DISTCLEANFILES = fcs-package.info

# Subdirectories to descend into.
SUBDIRS = src sl_pepc openpa

sublibs = src/liblpepc.a sl_pepc/libsl_pepc.a openpa/src/libopa.a

if ENABLE_LIBRARY_INSTALL
lib_LIBRARIES = libfcs_pepc.a
else
noinst_LIBRARIES = libfcs_pepc.a
endif

libfcs_pepc_a_SOURCES =

libfcs_pepc.a: $(sublibs)
	$(mkinstalldirs) .libs
	rm -f $@
	for l in $^ ; do l="./$${l}" lib="$${l##*/}" ; dir="$${l%$${lib}}" ; $(add_dir_lib_objects) ; done
	$(RANLIB) $@
