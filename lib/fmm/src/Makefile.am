# Build the current directory first, so the needed modules exist.
SUBDIRS = \
	. \
	errorcontrol/0d \
	errorcontrol/1d \
	errorcontrol/2d \
	errorcontrol/3d
SUBDIRS += \
	unrolled/fmmoopn \
	unrolled/fmmmopn \
	unrolled/fmmgradt

# Get Fortran compile rules that include preprocessing.
include $(top_srcdir)/build-aux/fortran-rules.am

LTOBJEXT = lo

# Create object files in subdirs, so the inference
# rules from fortran-rules.am still apply.
AUTOMAKE_OPTIONS = subdir-objects

FCCPPFLAGS = $(FC_MODINC).
default_FCFLAGS = $(FC_INTEGER_SIZE_8) $(FC_NO_RANGE_CHECK)
AM_FCFLAGS = $(default_FCFLAGS)

noinst_LTLIBRARIES = libfmm.la

addrarithm_files = addrarithm/diffcpointers.c \
addrarithm/cptrsize.c \
addrarithm/ptroffset.c

fortran_sources = \
fmmkinds.f \
fmmkinds.h \
mp_constants.f90 \
memalloc.f90 \
fmmwigner.f \
memalloc.f90 \
mp_wrapper.f90 \
mp_wrapper.h \
mp_info.f \
fmmalloc.f \
pvlist.f90 \
errorcontrol/fmmgetgzyx.f \
neighbors/getcjpa.f90 \
neighbors/caladdr2.f90 \
neighbors/getdist.f90 \
neighbors/getdistms.f90 \
neighbors/getneighbors.f90 \
fmm_fcs_binding.f90 \
fmminit.f \
fmm.f \
fmm_tune.f \
cfmm_tune.f \
fmm_cbindings.f90 \
fmm_cbindings.h \
fmm.h

libfmm_la_SOURCES = $(fortran_sources) \
dummy_malloc.c

libfmm_la_SOURCES += $(addrarithm_files)

if ENABLE_IBM_F_INTRINSICS
libfmm_la_SOURCES += pass2trfrqdcach.f90
pass2trfrqdcach.$(LTOBJEXT) : FCFLAGS:=$(FCFLAGS) -qarch=450d
else
libfmm_la_SOURCES += pass2trfrqdcach.legacy.f90
endif

#include passes/pass2/Makefile.include
#libfmm_la_SOURCES += $(pass2_sources)

# These files are to be treated as fixed form source.
fmm.$(LTOBJEXT) \
fmminit.$(LTOBJEXT) \
fmmalloc.$(LTOBJEXT) \
fmmkinds.$(LTOBJEXT) \
fmmwigner.$(LTOBJEXT) \
fmm_tune.$(LTOBJEXT) \
cfmm_tune.$(LTOBJEXT) \
mp_info.$(LTOBJEXT) \
errorcontrol/fmmgetgzyx.$(LTOBJEXT) : \
   AM_FCFLAGS = $(default_FCFLAGS) $(FC_FIXEDFORM)

# These files are large and should be built with optimization disabled.
errorcontrol/fmmgetgzyx.$(LTOBJEXT) \
neighbors/getdist.$(LTOBJEXT) neighbors/getdistms.$(LTOBJEXT) : \
   FCFLAGS += -O0

# The MPI/ARMCI constants may not be compiled with changed default integer size.
mp_constants.$(LTOBJEXT) : \
   AM_FCFLAGS =

# Header dependencies.
fmmkinds.$(LTOBJEXT): fmm.h ../fconfig.h
fmm.$(LTOBJEXT): fmm.h ../fconfig.h
fmmalloc.$(LTOBJEXT): fmm.h ../fconfig.h
mp_info.$(LTOBJEXT): fmm.h ../fconfig.h
pvlist.$(LTOBJEXT): fmm.h ../fconfig.h
fmm.$(LTOBJEXT): fmm.h ../fconfig.h
mp_wrapper.$(LTOBJEXT): fmm.h ../fconfig.h
if !ENABLE_IBM_F_INTRINSICS
pass2trfrqdcach.legacy.$(LTOBJEXT): fmm.h ../fconfig.h
endif

# Use generated module dependency rules.
$(eval -include $(builddir)/module.deps)

$(builddir)/module.deps: $(fortran_sources) fmm.h ../fconfig.h
	SRCDIR="$(srcdir)" FCCPP="$(FCCPP) -I$(builddir) -I$(srcdir) -I$(top_builddir) $(FCCPPFLAGS)" CPP="$(CPP)" OBJEXT="$(LTOBJEXT)" MODEXT="$(FC_MODEXT)" \
	$(SHELL) $(top_srcdir)/build-aux/fortran-depcomp $^ > $@

DISTCLEANFILES = module.deps
