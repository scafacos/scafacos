# Directory of ipnfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel 

# Headers that are installed.
if ENABLE_HEADER_INSTALL
include_HEADERS = pnfft.h
# provide separate header for long precision, since some compilers do not support it
nodist_include_HEADERS = pnfft.f03 pnfftl.f03
else
noinst_HEADERS = pnfft.h
endif
	
noinst_LTLIBRARIES = libapi.la

EXTRA_DIST = f03-api.sh f03-wrap.sh genf03-api.pl genf03-wrap.pl build-f03-interface.sh

# Group local sources into convenience lib.
libapi_la_SOURCES = \
        api-basic.c \
        api-adv.c \
        api-guru.c \
	pnfft.f03.in \
	pnfftl.f03.in \
	f03-wrap.c

