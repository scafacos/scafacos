include $(top_srcdir)/src/common-rules.am
include $(srcdir)/filelist.am

AM_FCFLAGS += $(FC_FIXEDFORM)
noinst_LTLIBRARIES = libfmmoopn.la

libfmmoopn_la_SOURCES = fmmoopnmain.f fmmoopnp00.f
if ENABLE_FMM_UNROLLED
libfmmoopn_la_SOURCES += $(filelist)
endif

# Header dependencies.
*.$(OBJEXT): $(top_srcdir)/src/fmm.h $(top_builddir)/fconfig.h
