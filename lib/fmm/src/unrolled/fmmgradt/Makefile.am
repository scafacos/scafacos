include $(top_srcdir)/src/common-rules.am
include $(srcdir)/filelist.am

AM_FCFLAGS += $(FC_FIXEDFORM)
noinst_LIBRARIES = libfmmgradt.a

libfmmgradt_a_SOURCES = fmmgradtmain.f fmmgradtp00.f
if ENABLE_FMM_UNROLLED
libfmmgradt_a_SOURCES += $(filelist)
endif

# Header dependencies.
*.$(OBJEXT): $(top_srcdir)/src/fmm.h $(top_builddir)/fconfig.h
