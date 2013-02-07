if ENABLE_SINGLE_LIB
noinst_LIBRARIES = libfcs_ewald.a
else
lib_LIBRARIES = libfcs_ewald.a
endif

libfcs_ewald_a_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/lib -I$(top_srcdir)/lib/common/fcs-common
libfcs_ewald_a_SOURCES = ewald.c
