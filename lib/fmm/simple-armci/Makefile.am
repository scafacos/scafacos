
lib_LTLIBARIES =
noinst_LTLIBRARIES = libsimple_armci.la

noinst_LTLIBRARIES += libparmci.la

SUBDIRS = generic
libparmci_la_CPPFLAGS = -I$(srcdir)/generic
simple_armci_libs = libparmci.la generic/liba1g.la

if ENABLE_SIMPLE_ARMCI_DCMFD
SUBDIRS += dcmfd
libparmci_la_CPPFLAGS += -DDCMFD
simple_armci_libs += dcmfd/liba1d.la
endif

if ENABLE_SIMPLE_ARMCI_PAMID
SUBDIRS += pamid
libparmci_la_CPPFLAGS += -DPAMID
simple_armci_libs += pamid/liba1d.la
endif

libparmci_la_SOURCES = \
  parmci.c parmci.h

libsimple_armci_la_SOURCES =

libsimple_armci_la_LIBADD = $(simple_armci_libs)
