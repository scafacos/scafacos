# Let aclocal find the add-on macros.
ACLOCAL_AMFLAGS = -I m4

# Distribute script for Fortran module dependency generation.
EXTRA_DIST = build-aux/fortran-depcomp

# Remove files created during configure.
DISTCLEANFILES = fcs-package.info

SUBDIRS =
if ENABLE_FMM_ARMCI
SUBDIRS += armci
endif
if ENABLE_FMM_A1
SUBDIRS += a1
endif
if ENABLE_FMM_SIMPLE_ARMCI
SUBDIRS += simple-armci
endif
if ENABLE_FMM_ARMCIMPI
SUBDIRS += armci-mpi
endif
SUBDIRS += sl_fmm src

DIST_SUBDIRS = $(SUBDIRS)

fmmlibs = \
	src/libfmm.la \
	src/errorcontrol/0d/lib0d.la \
	src/errorcontrol/3d/lib3d.la \
	src/errorcontrol/1d/lib1d.la \
	src/errorcontrol/2d/lib2d.la \
	src/unrolled/fmmmopn/libfmmmopn.la \
	src/unrolled/fmmgradt/libfmmgradt.la \
	src/unrolled/fmmoopn/libfmmoopn.la \
	sl_fmm/libsl_fmm.la

if ENABLE_FMM_ARMCI
fmmlibs += \
	armci/libarmci.la
endif

if ENABLE_FMM_A1
fmmlibs += \
	a1/liba1.la
endif

if ENABLE_FMM_SIMPLE_ARMCI
fmmlibs += \
	simple-armci/libsimple_armci.a
endif

if ENABLE_FMM_ARMCIMPI
fmmlibs += \
	armci-mpi/libarmci.la
endif


if ENABLE_LIBRARY_INSTALL
lib_LTLIBRARIES = libfcs_fmm.la
else
noinst_LTLIBRARIES = libfcs_fmm.la
endif

libfcs_fmm_la_SOURCES =

libfcs_fmm_la_LIBADD = $(fmmlibs)
