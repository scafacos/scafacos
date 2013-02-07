# include the list of codelets

include $(top_srcdir)/rdft/simd/codlist.mk

ALL_CODELETS = $(SIMD_CODELETS)
BUILT_SOURCES= $(SIMD_CODELETS) $(CODLIST)
DISTCLEANFILES= $(CODLIST)
EXTRA_DIST = $(BUILT_SOURCES) genus.c
INCLUDE_SIMD_HEADER="\#include SIMD_HEADER"
XRENAME=XSIMD
SOLVTAB_NAME = XSIMD(solvtab_rdft)

# include special rules for regenerating codelets.
include $(top_srcdir)/support/Makefile.codelets

