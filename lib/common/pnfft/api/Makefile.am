# Directory of ipnfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel 

# Headers that are installed.
if ENABLE_HEADER_INSTALL
include_HEADERS = pnfft.h
nodist_include_HEADERS = pnfft.f
else
noinst_HEADERS = pnfft.h
endif
	
noinst_LTLIBRARIES = libapi.la

# Group local sources into convenience lib.
libapi_la_SOURCES = \
        api-basic.c \
        api-adv.c \
        api-guru.c \
	pnfft.f.in \
	fortran-api.c \
	fortran-mangling.h \
	fortran-prototypes.h \
	fortran-wrappers.h

