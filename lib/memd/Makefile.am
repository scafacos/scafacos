if ENABLE_SINGLE_LIB
noinst_LTLIBRARIES = libfcs_memd.la
else
lib_LTLIBRARIES = libfcs_memd.la
endif

libfcs_memd_la_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/lib -I$(top_srcdir)/lib/common/fcs-common
libfcs_memd_la_SOURCES = \
	helper_functions.h helper_functions.c \
	communication.h communication.c \
	getter_setter.h getter_setter.c \
	init.h init.c \
	interpolate.h interpolate.c \
	initial_solution.h initial_solution.c \
	forces.h forces.c \
	run.h run.c \
	observables.h observables.c \
	data_types.h \
	memd.h
