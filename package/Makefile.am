# Let aclocal find the add-on macros.
ACLOCAL_AMFLAGS = -I m4

# Import central add_lib_objects and add_subdir_lib_objects support.
include $(top_srcdir)/build-aux/lib_a_support.am

# Install pkg-config file in correct location.
pkgconfigdir = $(libdir)/pkg-config
pkgconfig_DATA = scafacos-fcs.pc

if ENABLE_SINGLE_LIB

fcslibs = $(SCAFACOS_MK_LIBS_A)

lib_LIBRARIES = libfcs.a
libfcs_a_SOURCES =

libfcs.a: $(fcslibs)
	$(mkinstalldirs) .libs
	rm -f $@
	for l in $^ ; do dir=.; lib=$$l; $(add_dir_lib_objects) ; done
	$(RANLIB) $@

all-local: libfcs.a

clean-local: clean-lib
clean-lib:
	-rm -f libfcs.a
	-rm -rf .libs

endif
