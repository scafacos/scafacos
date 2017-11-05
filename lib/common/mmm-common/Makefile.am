if ENABLE_SINGLE_LIB
noinst_LTLIBRARIES = libfcs_mmm.la
else
lib_LTLIBRARIES = libfcs_mmm.la
endif

libfcs_mmm_la_CPPFLAGS = -I$(top_srcdir)/lib/common/fcs-common

libfcs_mmm_la_SOURCES = \
    mmm-common.c mmm-common.h \
    specfunc.c specfunc.h
