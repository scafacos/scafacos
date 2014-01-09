if ENABLE_SINGLE_LIB
noinst_LTLIBRARIES = libfcs_p3m.la
else
lib_LTLIBRARIES = libfcs_p3m.la
endif

libfcs_p3m_la_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/lib -I$(top_srcdir)/lib/common/fcs-common $(fftw3_CPPFLAGS)
libfcs_p3m_la_SOURCES = \
	utils.hpp \
	communication.hpp communication.cpp \
	caf.cpp caf.hpp \
	fft.cpp fft.hpp \
	influence_function.cpp influence_function.hpp \
	init.cpp init.hpp \
	parameters.cpp parameters.hpp \
	prepare.cpp prepare.hpp \
	timing.cpp timing.hpp \
	tune.cpp tune.hpp \
	tune_broadcast.cpp tune_broadcast.hpp \
	run.cpp run.hpp \
	p3m.hpp \
	types.hpp \
	error_estimate.cpp error_estimate.hpp
