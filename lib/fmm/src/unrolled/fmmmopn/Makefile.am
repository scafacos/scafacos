include $(top_srcdir)/src/common-rules.am
include $(srcdir)/filelist.am

AM_FCFLAGS += $(FC_FIXEDFORM)
noinst_LTLIBRARIES = libfmmmopn.la

libfmmmopn_la_SOURCES = fmmmopnmain.f fmmmopnp00.f
if ENABLE_FMM_UNROLLED
libfmmmopn_la_SOURCES += $(filelist)
endif

# Header dependencies.
*.$(OBJEXT): $(top_srcdir)/src/fmm.h $(top_builddir)/fconfig.h
