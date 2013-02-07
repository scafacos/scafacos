if ENABLE_SINGLE_LIB
noinst_LIBRARIES = libfcs_direct.a
else
lib_LIBRARIES = libfcs_direct.a
endif

libfcs_direct_a_CPPFLAGS = -I$(top_srcdir)/lib
libfcs_direct_a_SOURCES = \
  directc.c directc.h \
  z_tools.h
