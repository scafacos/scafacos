if ENABLE_SINGLE_LIB
noinst_LIBRARIES = libfcs_mmm2d.a
else
lib_LIBRARIES = libfcs_mmm2d.a
endif

libfcs_mmm2d_a_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/lib -I$(top_srcdir)/lib/common/fcs-common
libfcs_mmm2d_a_SOURCES = \
	init.c init.h \
	parameters.c parameters.h \
	tune.c tune.h \
	communication.c communication.h \
	run.c run.h \
	types.h \
	mmm2d.h \
	../common/mmm-common/mmm-common.c \
	../common/mmm-common/mmm-common.h \
	../common/mmm-common/specfunc.c \
	../common/mmm-common/specfunc.h

