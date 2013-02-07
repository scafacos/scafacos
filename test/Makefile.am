###################################################
# Makefile for the ScaFaCos fcs test cases        #
#                                                 #
###################################################

# Build current directory first, for the inp_data symlink.

SUBDIRS = . c

if HAVE_IOMANIP
SUBDIRS += generic
endif

if ENABLE_FORTRAN_INTERFACE
SUBDIRS += fortran
endif

# Add a symlink to the in-tree inp_data in a VPATH build.
all-local: make-data-symlink defs
make-data-symlink:
	-@if test "$(srcdir)" != . && test ! -r inp_data; then \
	  ln -s "$(srcdir)/inp_data" .; \
	fi

distclean-local: clean-data-symlink
clean-data-symlink:
	-@if test "$(srcdir)" != .; then \
	  rm -f inp_data; \
	fi

.PHONY: make-data-symlink clean-data-symlink


# Distribute all input data, but not any SVN files.
dist-hook:
	cp -pR $(srcdir)/inp_data $(distdir)/inp_data
	rm -rf `find $(distdir)/inp_data -name .svn`
