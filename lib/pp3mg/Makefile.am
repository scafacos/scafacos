# Remove files created during configure.
DISTCLEANFILES = fcs-package.info

# Subdirectories to descend into.
SUBDIRS = pmg pp3mg

sublibs = pmg/libpmg.a pp3mg/libpp3mg.a

# Let aclocal find the add-on macros.
ACLOCAL_AMFLAGS = -I m4

# Import central add_lib_objects support.
include $(top_srcdir)/build-aux/lib_a_support.am

if ENABLE_LIBRARY_INSTALL
lib_LIBRARIES = libfcs_pp3mg.a
else
noinst_LIBRARIES = libfcs_pp3mg.a
endif

libfcs_pp3mg_a_SOURCES =

libfcs_pp3mg.a: $(sublibs)
	$(mkinstalldirs) .libs
	rm -f $@
	for l in $^ ; do l="./$${l}" lib="$${l##*/}" ; dir="$${l%$${lib}}" ; $(add_lib_objects) ; done
	$(RANLIB) $@
