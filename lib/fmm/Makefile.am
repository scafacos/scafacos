# Let aclocal find the add-on macros.
ACLOCAL_AMFLAGS = -I m4

# Import central add_lib_objects and add_subdir_lib_objects support.
include $(top_srcdir)/build-aux/lib_a_support.am

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
SUBDIRS += sl_fmm src

DIST_SUBDIRS = $(SUBDIRS)

fmmlibs = \
	src/libfmm.a \
	src/errorcontrol/0d/lib0d.a \
	src/errorcontrol/3d/lib3d.a \
	src/errorcontrol/1d/lib1d.a \
	src/errorcontrol/2d/lib2d.a \
	src/unrolled/fmmmopn/libfmmmopn.a \
	src/unrolled/fmmgradt/libfmmgradt.a \
	src/unrolled/fmmoopn/libfmmoopn.a \
	sl_fmm/libsl_fmm.a

if ENABLE_FMM_ARMCI
fmmlibs += \
	armci/.libs/libarmci.a
endif

if ENABLE_FMM_A1
fmmlibs += \
	a1/.libs/liba1.a
endif

if ENABLE_FMM_SIMPLE_ARMCI
fmmlibs += \
	simple-armci/libsimple_armci.a
endif


if ENABLE_LIBRARY_INSTALL
lib_LIBRARIES = libfcs_fmm.a
else
noinst_LIBRARIES = libfcs_fmm.a
endif

libfcs_fmm_a_SOURCES =

libfcs_fmm.a: $(fmmlibs)
	$(mkinstalldirs) .libs
	rm -f $@
	dir=src; lib=libfmm.a; $(add_subdir_lib_objects)
	dir=src/errorcontrol/0d; lib=lib0d.a; $(add_lib_objects)
	dir=src/errorcontrol/1d; lib=lib1d.a; $(add_lib_objects)
	dir=src/errorcontrol/2d; lib=lib2d.a; $(add_lib_objects)
	dir=src/errorcontrol/3d; lib=lib3d.a; $(add_lib_objects)
	dir=src/unrolled/fmmgradt; lib=libfmmgradt.a; $(add_lib_objects)
	dir=src/unrolled/fmmmopn; lib=libfmmmopn.a; $(add_lib_objects)
	dir=src/unrolled/fmmoopn; lib=libfmmoopn.a; $(add_lib_objects)
	dir=sl_fmm; lib=libsl_fmm.a; $(add_lib_objects)
if ENABLE_FMM_ARMCI
	dir=armci/.libs; lib=../.libs/libarmci.a; $(add_subdir_lib_objects)
endif
if ENABLE_FMM_A1
	dir=a1; $(add_dir_objects)
endif
if ENABLE_FMM_SIMPLE_ARMCI
	dir=simple-armci; lib=libsimple_armci.a; $(add_subdir_lib_objects)
endif
	$(RANLIB) $@

all-local: libfcs_fmm.a

clean-local: clean-lib
clean-lib:
	-rm -f libfcs_fmm.a
	-rm -rf .libs
