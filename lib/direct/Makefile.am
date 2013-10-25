if ENABLE_SINGLE_LIB
noinst_LTLIBRARIES = libfcs_direct.la
else
lib_LTLIBRARIES = libfcs_direct.la
endif

libfcs_direct_la_CPPFLAGS = -I$(top_srcdir)/lib
libfcs_direct_la_SOURCES = \
  directc.c directc.h \
  z_tools.h
