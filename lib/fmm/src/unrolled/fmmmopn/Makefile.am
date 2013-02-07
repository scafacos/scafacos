include $(top_srcdir)/src/common-rules.am
include $(srcdir)/filelist.am

AM_FCFLAGS += $(FC_FIXEDFORM)
noinst_LIBRARIES = libfmmmopn.a

libfmmmopn_a_SOURCES = fmmmopnmain.f fmmmopnp00.f
if ENABLE_FMM_UNROLLED
libfmmmopn_a_SOURCES += $(filelist)
endif

# Header dependencies.
*.$(OBJEXT): $(top_srcdir)/src/fmm.h $(top_builddir)/fconfig.h
