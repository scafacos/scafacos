AM_CPPFLAGS = -I$(top_srcdir)/kernel -I$(top_srcdir)/dft	\
-I$(top_srcdir)/rdft -I$(top_srcdir)/api
AM_CFLAGS = $(STACK_ALIGN_CFLAGS)

if OPENMP
FFTWOMPLIB = libfcs_fftw3@PREC_SUFFIX@_omp.la
else
FFTWOMPLIB = 
endif

if THREADS
if COMBINED_THREADS
noinst_LTLIBRARIES = libfcs_fftw3@PREC_SUFFIX@_threads.la
else
lib_LTLIBRARIES = libfcs_fftw3@PREC_SUFFIX@_threads.la $(FFTWOMPLIB)
endif
else
lib_LTLIBRARIES = $(FFTWOMPLIB)
endif

# pkgincludedir = $(includedir)/fftw3@PREC_SUFFIX@
# pkginclude_HEADERS = threads.h

libfcs_fftw3@PREC_SUFFIX@_threads_la_SOURCES = api.c conf.c threads.c	\
threads.h dft-vrank-geq1.c ct.c rdft-vrank-geq1.c hc2hc.c		\
vrank-geq1-rdft2.c f77api.c f77funcs.h
libfcs_fftw3@PREC_SUFFIX@_threads_la_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
libfcs_fftw3@PREC_SUFFIX@_threads_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@
if !COMBINED_THREADS
libfcs_fftw3@PREC_SUFFIX@_threads_la_LIBADD = ../libfcs_fftw3@PREC_SUFFIX@.la
endif

libfcs_fftw3@PREC_SUFFIX@_omp_la_SOURCES = api.c conf.c openmp.c	\
threads.h dft-vrank-geq1.c ct.c rdft-vrank-geq1.c hc2hc.c	\
vrank-geq1-rdft2.c f77api.c f77funcs.h
libfcs_fftw3@PREC_SUFFIX@_omp_la_CFLAGS = $(AM_CFLAGS) $(OPENMP_CFLAGS)
libfcs_fftw3@PREC_SUFFIX@_omp_la_LDFLAGS = -version-info @SHARED_VERSION_INFO@
if !COMBINED_THREADS
libfcs_fftw3@PREC_SUFFIX@_omp_la_LIBADD = ../libfcs_fftw3@PREC_SUFFIX@.la
endif
