# Directory of ipnfft.h
AM_CPPFLAGS = -I$(top_srcdir)/kernel 

# Directory of pnfft.h
AM_CPPFLAGS += -I$(top_srcdir)/api 

noinst_LTLIBRARIES = libkernel.la

# Group local sources into convenience lib.
libkernel_la_SOURCES = \
	ndft-parallel.c \
	assign.c \
	matrix_D.c \
	matrix_D.h \
	bessel_i0.c \
	bessel_i0.h \
	bessel_i1.c \
	bessel_i1.h \
	bspline.c \
	bspline.h \
	sinc.c \
	sinc.h \
	malloc.c \
	timer.c \
	check.c \
	ipnfft.h
