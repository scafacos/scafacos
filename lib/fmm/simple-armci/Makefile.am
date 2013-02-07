
# Import central add_lib_objects and add_subdir_lib_objects support.
include $(top_srcdir)/build-aux/lib_a_support.am

lib_LIBARIES =
noinst_LIBRARIES = libsimple_armci.a

noinst_LIBRARIES += libparmci.a

SUBDIRS = generic
libparmci_a_CPPFLAGS = -I$(srcdir)/generic
simple_armci_libs = libparmci.a generic/liba1g.a

if ENABLE_SIMPLE_ARMCI_DCMFD
SUBDIRS += dcmfd
libparmci_a_CPPFLAGS += -DDCMFD
simple_armci_libs += dcmfd/liba1d.a
endif

if ENABLE_SIMPLE_ARMCI_PAMID
SUBDIRS += pamid
libparmci_a_CPPFLAGS += -DPAMID
simple_armci_libs += pamid/liba1d.a
endif

libparmci_a_SOURCES = \
  parmci.c parmci.h

libsimple_armci_a_SOURCES =

libsimple_armci.a: $(simple_armci_libs)
	$(mkinstalldirs) .libs
	rm -f $@
	dir=.; lib=libparmci.a; $(add_subdir_lib_objects)
	dir=generic; lib=liba1g.a; $(add_subdir_lib_objects)
if ENABLE_SIMPLE_ARMCI_DCMFD
	dir=dcmfd; lib=liba1d.a; $(add_subdir_lib_objects)
endif
if ENABLE_SIMPLE_ARMCI_PAMID
	dir=pamid; lib=liba1d.a; $(add_subdir_lib_objects)
endif
	$(RANLIB) $@

all-local: libsimple_armci.a

clean-local: clean-lib
clean-lib:
	-rm -f libsimple_armci.a
	-rm -rf .libs
