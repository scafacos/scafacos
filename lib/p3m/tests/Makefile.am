DUMMY:=$(eval include $(top_builddir)/package/scafacos-fcs.mk)

check_PROGRAMS = compute_error
compute_error_SOURCES = compute_error.cpp
compute_error_LDADD = $(SCAFACOS_MK_LDADD) $(top_builddir)/src/libfcs.la
compute_error_CPPFLAGS = -I $(top_srcdir)/src -I $(srcdir)/.. $(fftw3_CPPFLAGS)
