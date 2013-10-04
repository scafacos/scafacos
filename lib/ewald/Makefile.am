if ENABLE_SINGLE_LIB
noinst_LTLIBRARIES = libfcs_ewald.la
else
lib_LTLIBRARIES = libfcs_ewald.la
endif

libfcs_ewald_la_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/lib -I$(top_srcdir)/lib/common/fcs-common
libfcs_ewald_la_SOURCES = ewald.c
