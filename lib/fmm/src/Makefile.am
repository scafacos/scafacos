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

# Create object files in subdirs, so the inference
# rules from fortran-rules.am still apply.
AUTOMAKE_OPTIONS = subdir-objects

FCCPPFLAGS = $(FC_MODINC).
default_FCFLAGS = $(FC_INTEGER_SIZE_8) $(FC_NO_RANGE_CHECK)
AM_FCFLAGS = $(default_FCFLAGS)

noinst_LIBRARIES = libfmm.a

addrarithm_files = addrarithm/diffcpointers.c \
addrarithm/cptrsize.c \
addrarithm/ptroffset.c

libfmm_a_SOURCES = \
fmmkinds.f \
fmmkinds.h \
mp_constants.f90 \
memalloc.f90 \
fmmwigner.f \
dummy_malloc.c \
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

libfmm_a_SOURCES += $(addrarithm_files)

if ENABLE_IBM_F_INTRINSICS
libfmm_a_SOURCES += pass2trfrqdcach.f90
pass2trfrqdcach.$(OBJEXT) : FCFLAGS:=$(FCFLAGS) -qarch=450d
else
libfmm_a_SOURCES += pass2trfrqdcach.legacy.f90
endif

#include passes/pass2/Makefile.include
#libfmm_a_SOURCES += $(pass2_sources)

# These files are to be treated as fixed form source.
fmm.$(OBJEXT) \
fmminit.$(OBJEXT) \
fmmalloc.$(OBJEXT) \
fmmkinds.$(OBJEXT) \
fmmwigner.$(OBJEXT) \
fmm_tune.$(OBJEXT) \
cfmm_tune.$(OBJEXT) \
mp_info.$(OBJEXT) \
errorcontrol/fmmgetgzyx.$(OBJEXT) : \
   AM_FCFLAGS = $(default_FCFLAGS) $(FC_FIXEDFORM)

# These files are large and should be built with optimization disabled.
errorcontrol/fmmgetgzyx.$(OBJEXT) \
neighbors/getdist.$(OBJEXT) neighbors/getdistms.$(OBJEXT) : \
   FCFLAGS += -O0

# The MPI/ARMCI constants may not be compiled with changed default integer size.
mp_constants.$(OBJEXT) : \
   AM_FCFLAGS =

# Header dependencies.
fmmkinds.$(OBJEXT): fmm.h ../fconfig.h
fmm.$(OBJEXT): fmm.h ../fconfig.h
fmmalloc.$(OBJEXT): fmm.h ../fconfig.h
mp_info.$(OBJEXT): fmm.h ../fconfig.h
pvlist.$(OBJEXT): fmm.h ../fconfig.h
fmm.$(OBJEXT): fmm.h ../fconfig.h
mp_wrapper.$(OBJEXT): fmm.h ../fconfig.h
if !ENABLE_IBM_F_INTRINSICS
pass2trfrqdcach.legacy.$(OBJEXT): fmm.h ../fconfig.h
endif

# Use generated module dependency rules.
$(eval -include $(builddir)/module.deps)

$(builddir)/module.deps: $(libfmm_a_SOURCES) fmm.h ../fconfig.h
	SRCDIR="$(srcdir)" CPP="cpp -I$(builddir) -I$(srcdir) -I$(top_builddir) $(FCCPPFLAGS)" OBJEXT="$(OBJEXT)" MODEXT="$(FC_MODEXT)" \
	$(SHELL) $(top_srcdir)/build-aux/fortran-depcomp $^ > $@

DISTCLEANFILES = module.deps
