include $(top_srcdir)/src/common-rules.am
include $(srcdir)/filelist.am

AM_FCFLAGS += $(FC_FIXEDFORM)
noinst_LTLIBRARIES = libfmmgradt.la

libfmmgradt_la_SOURCES = fmmgradtmain.f fmmgradtp00.f
if ENABLE_FMM_UNROLLED
libfmmgradt_la_SOURCES += $(filelist)
endif

# Header dependencies.
*.$(OBJEXT): $(top_srcdir)/src/fmm.h $(top_builddir)/fconfig.h
