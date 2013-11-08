if ENABLE_SINGLE_LIB
noinst_LTLIBRARIES = libfcs_p3m.la
else
lib_LTLIBRARIES = libfcs_p3m.la
endif

libfcs_p3m_la_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/lib -I$(top_srcdir)/lib/common/fcs-common @fftw3_CPPFLAGS@
libfcs_p3m_la_SOURCES = \
	utils.h utils.c \
	communication.h communication.c \
	caf.c caf.h \
	fft.c fft.h \
	influence_function.c influence_function.h \
	init.c init.h \
	parameters.c parameters.h \
	prepare.c prepare.h \
	timing.c timing.h \
	tune.c tune.h \
	tune_broadcast.c tune_broadcast.h \
	run.c run.h \
	p3m.h \
	types.h \
	error_estimate.c \
	error_estimate.h
