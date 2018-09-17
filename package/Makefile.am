# Let aclocal find the add-on macros.
ACLOCAL_AMFLAGS = -I m4

# Install pkg-config file in correct location.
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = scafacos.pc

# Install cmake package file.
datadir = $(datarootdir)/cmake/scafacos
data_DATA = scafacosConfig.cmake

if ENABLE_SINGLE_LIB

fcslibs = $(SCAFACOS_MK_LIBS_A)

lib_LTLIBRARIES = libfcs.la

libfcs_la_SOURCES =

libfcs_la_LIBADD = $(fcslibs)

endif
