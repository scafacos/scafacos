# Directory of ipnfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel 

# Headers that are installed.
if ENABLE_HEADER_INSTALL
include_HEADERS = pnfft.h
nodist_include_HEADERS = pnfft.f pnfft.f03
else
noinst_HEADERS = pnfft.h
endif
	
noinst_LTLIBRARIES = libapi.la

EXTRA_DIST = f03-api.sh f03-wrap.sh genf03-api.pl genf03-wrap.pl
CLEANFILES = pnfft.f03 pnfft.f03.in pnfftl.f03 pnfftl.f03.in f03-wrap.c

# Group local sources into convenience lib.
libapi_la_SOURCES = \
        api-basic.c \
        api-adv.c \
        api-guru.c \
	pnfft.f.in \
	fortran-api.c \
	fortran-mangling.h \
	fortran-prototypes.h \
	fortran-wrappers.h \
	$(builddir)/f03-wrap.c

# pnfft.f03.in: pnfft.h f03-api.sh genf03-api.pl
# 	sh $(srcdir)/f03-api.sh d f l > $@
# 
# pnfft.f03: pnfft.f03.in
# 	sed 's/C_MPI_FINT/@C_MPI_FINT@/' $(srcdir)/pnfft.f03.in > $@
# 
# f03-wrap.c: pnfft.h f03-wrap.sh genf03-wrap.pl
# 	sh $(srcdir)/f03-wrap.sh > $@


# provide separate header for long precision, since some compilers do not support it
pnfft.f03: pnfft.f03.in
	sed 's/C_MPI_FINT/@C_MPI_FINT@/' pnfft.f03.in > $@

pnfftl.f03: pnfftl.f03.in
	sed 's/C_MPI_FINT/@C_MPI_FINT@/' pnfftl.f03.in > $@


if TRUE
# enable only in maintainer mode?
pnfft.f03.in: pnfft.h f03-api.sh genf03-api.pl
	(cd $(srcdir); sh f03-api.sh d f > $(abs_builddir)/$@)

pnfftl.f03.in: pnfft.h f03-api.sh genf03-api.pl
	(cd $(srcdir); sh f03-api.sh l > $(abs_builddir)/$@)

f03-wrap.c: pnfft.h f03-wrap.sh genf03-wrap.pl
	(cd $(srcdir); sh f03-wrap.sh >  $(abs_builddir)/$@)
endif
