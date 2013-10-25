# Directories
SUBDIRS = 
if ENABLE_COMMON_NEAR
SUBDIRS += common/near
endif
if ENABLE_COMMON_GRIDSORT
SUBDIRS += common/gridsort
endif
if ENABLE_COMMON_RESORT
SUBDIRS += common/resort
endif
SUBDIRS += common/fcs-common
if ENABLE_COMMON_FFTW
SUBDIRS += common/fftw-3.3
endif
if ENABLE_COMMON_PFFT
SUBDIRS += common/pfft
endif
if ENABLE_COMMON_PNFFT
SUBDIRS += common/pnfft
endif
if ENABLE_DIRECT
SUBDIRS += direct
endif
if ENABLE_EWALD
SUBDIRS += ewald
endif
if ENABLE_FMM
SUBDIRS += fmm
endif
if ENABLE_MEMD
SUBDIRS += memd
endif
if ENABLE_MMM1D
SUBDIRS += mmm1d
endif
if ENABLE_MMM2D
SUBDIRS += mmm2d
endif
if ENABLE_P2NFFT
SUBDIRS += p2nfft
endif
if ENABLE_P3M
SUBDIRS += p3m
endif
if ENABLE_PEPC
SUBDIRS += pepc
endif
if ENABLE_PP3MG
SUBDIRS += pp3mg
endif
if ENABLE_VMG
SUBDIRS += vmg
endif
if ENABLE_WOLF
SUBDIRS += wolf
endif

# For disabled solvers, disregard the whole tree below.
# (This means that the person creating a distribution needs
# to enable all solvers.)
DIST_SUBDIRS = $(SUBDIRS)
