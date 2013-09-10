if ENABLE_FORTRAN_INTERFACE
SUBDIRS = fortran
endif

if ENABLE_SINGLE_LIB
noinst_LIBRARIES = libfcs.a
else
lib_LIBRARIES = libfcs.a
endif

# Headers that will be installed.
include_HEADERS = \
  fcs.h \
  FCSDefinitions.h \
  FCSResult_p.h \
  FCSInterface_p.h

libfcs_a_CPPFLAGS = -I$(top_srcdir)/lib/common/fcs-common -I$(top_srcdir)/lib

libfcs_a_SOURCES = \
  FCSInterface.h \
  FCSInterface.c \
  FCSResult.h \
  FCSResult.c

if ENABLE_DIRECT
libfcs_a_SOURCES += \
  fcs_direct.h \
  fcs_direct.c
include_HEADERS += \
  fcs_direct_p.h
endif

if ENABLE_EWALD
libfcs_a_SOURCES += \
  fcs_ewald.h
include_HEADERS += \
  fcs_ewald_p.h
endif

if ENABLE_FMM
libfcs_a_SOURCES += \
  fcs_fmm.h \
  fcs_fmm.c
include_HEADERS += \
  fcs_fmm_p.h
endif

if ENABLE_MEMD
libfcs_a_SOURCES += \
  fcs_memd.h \
  fcs_memd.c
include_HEADERS += \
  fcs_memd_p.h
endif

if ENABLE_MMM1D
libfcs_a_SOURCES += \
  fcs_mmm1d.h \
  fcs_mmm1d.c
include_HEADERS += \
  fcs_mmm1d_p.h
endif

if ENABLE_MMM2D
libfcs_a_SOURCES += \
  fcs_mmm2d.h \
  fcs_mmm2d.c
include_HEADERS += \
  fcs_mmm2d_p.h
endif

if ENABLE_P2NFFT
libfcs_a_SOURCES += \
  fcs_p2nfft.h \
  fcs_p2nfft.c
include_HEADERS += \
  fcs_p2nfft_p.h \
  fcs_p2nfft_wrappers.h \
  fcs_p2nfft_shortcuts.h
endif

if ENABLE_P3M
libfcs_a_SOURCES += \
  fcs_p3m.h \
  fcs_p3m.c
include_HEADERS += \
  fcs_p3m_p.h
endif

if ENABLE_PEPC
libfcs_a_SOURCES += \
  fcs_pepc.h \
  fcs_pepc.c
include_HEADERS += \
  fcs_pepc_p.h
endif

if ENABLE_PP3MG
libfcs_a_SOURCES += \
  fcs_pp3mg.h \
  fcs_pp3mg.c
include_HEADERS += \
  fcs_pp3mg_p.h
endif

if ENABLE_VMG
libfcs_a_SOURCES += \
  fcs_vmg.h \
  fcs_vmg.c
include_HEADERS += \
  fcs_vmg_p.h
endif

if ENABLE_WOLF
libfcs_a_SOURCES += \
  fcs_wolf.h \
  fcs_wolf.c
include_HEADERS += \
  fcs_wolf_p.h
endif
