if ENABLE_SINGLE_LIB
noinst_LTLIBRARIES = libfcs_mmm1d.la
else
lib_LTLIBRARIES = libfcs_mmm1d.la
endif

libfcs_mmm1d_la_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/lib -I$(top_srcdir)/lib/common/fcs-common
libfcs_mmm1d_la_SOURCES = \
	init.c init.h \
	parameters.c parameters.h \
	tune.c tune.h \
	run.c run.h \
	mmm1d.h \
	types.h
